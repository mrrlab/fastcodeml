
#include "OptSQP.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>
#include <iomanip>


// ----------------------------------------------------------------------
//	Class members definition: OptSQP
// ----------------------------------------------------------------------
double OptSQP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	mH1Optimization = (mN == mNumTimes+5);
	
	allocateMemory();
	
#ifdef SCALE_OPT_VARIABLES
	// set the scaling
	int i;
	// get a better scaling for the branch lengths variables.
	// it is often near ~0.25, multiply it by 4 so it is "more" around 1
	// in the new space representation
	mUpperBound.assign(mNumTimes, 200.0);
	
	
	i = mNumTimes + 1; 		// v1
	mUpperBound[i] = 20.0;
	
	i = mNumTimes + 2; 		// w0
	mUpperBound[i] = 50.0;
	
	
	// shrink the w2 variable between 0 and 1 so it is about the same scale as 
	// the other variables in the new space representation
	i = mNumTimes+4; 		// w2
	if (mH1Optimization)
	{
		mLowerBound[i] = 0.0;
		mUpperBound[i] = 1.0;
	}
	
#endif // SCALE_OPT_VARIABLES
	
	double maxl = 1e7;
	SQPminimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------
void OptSQP::allocateMemory(void)
{
	mSizeVect = mN*sizeof(double);
	
	mXEvaluator.resize(mN);
	mSpace.resize(2*mN*mN + mN*7);

	mWorkSpaceVect = &mSpace[0];
	mWorkSpaceMat = mWorkSpaceVect + mN;
	
	mGradient = mWorkSpaceMat + mN*mN;
	mP = mGradient + mN;
	mHessian = mP + mN;
	
	mSk = mHessian + mN*mN;
	mYk = mSk + mN;
	
	mXPrev = mYk + mN;
	mGradPrev = mXPrev + mN;
	
	// set active set counters to zero
	mActiveSet.resize(mN, 0);
}


#ifdef SCALE_OPT_VARIABLES
// ----------------------------------------------------------------------
void OptSQP::scaleVariables(double *aX)
{
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		
		double x = slb + (aX[i] - lb) * (sub-slb)/(ub-lb);
		
		x = min2(x, sub);
		x = max2(x, slb);
		
		aX[i] = x;
	}
}


// ----------------------------------------------------------------------
void OptSQP::unscaleVariables(double *aX)
{
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		double x;
		
		x = lb + (aX[i] - slb) * (ub-lb)/(sub-slb);
		
		x = min2(x, ub);
		x = max2(x, lb);
		
		aX[i] = x;
	}
}
#endif // SCALE_OPT_VARIABLES


// ----------------------------------------------------------------------
void OptSQP::SQPminimizer(double *aF, double *aX)
{
#ifdef SCALE_OPT_VARIABLES
	scaleVariables(aX);
#endif // SCALE_OPT_VARIABLES
	
	//double df = 0.0;
	*aF = evaluateFunction(aX, mTrace);
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "Initial point for SQP:";
		mModel->printVar(mXEvaluator, *aF);
	}
	
	// compute current gradient
	computeGradient(aX, *aF, mGradient);
	
	// initialize the hessian matrix
	hessianInitialization();
	
	// bounds for the QP subproblem
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	mQPsolver.reset(new BOXCQP(mN, &localLowerBound[0], &localUpperBound[0]));	
	
		
	// ----------------------------------------- main loop
	bool convergenceReached = false;
	for (mStep = 0; !convergenceReached; ++mStep)
	{
		// update local bounds
		memcpy(&localLowerBound[0], &mLowerBound[0], mSizeVect);
		memcpy(&localUpperBound[0], &mUpperBound[0], mSizeVect);
		daxpy_(&mN, &minus_one, aX, &I1, &localLowerBound[0], &I1);
		daxpy_(&mN, &minus_one, aX, &I1, &localUpperBound[0], &I1);
		
		// save current parameters
		double f_prev = *aF;
		memcpy(mGradPrev, mGradient, mSizeVect);
		memcpy(mXPrev, aX, mSizeVect);
		
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Quadratic program solving..." << std::endl;
		
		// solve quadratic program to get the search direction		
		bool QPsolutionOnBorder;
		bool QP_converged = mQPsolver->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder, mWorkSpaceVect);
		
		if (!QP_converged) // take the projected gradient direction
		{
			hessianInitialization();
			memcpy(mP, mGradient, mSizeVect);
			dscal_(&mN, &minus_one, mP, &I1);
			#pragma omp parallel for
			for (int i(0); i<mN; ++i)
			{
				mP[i] = max2(mP[i], localLowerBound[i]);
				mP[i] = min2(mP[i], localUpperBound[i]);
			}
		}
	
#if 1
		// try to extend the limits to the boundaries (-> global line search)
		double alpha = 1e16;
		for (int i(0); i<mN; ++i)
		{
			double l = localLowerBound[i];
			double u = localUpperBound[i];
			double p = mP[i];
			if (p < -1e-8)
				alpha = min2(alpha, l/p);
			else if (p > 1e-8)
				alpha = min2(alpha, u/p);
		}
#else		
		double alpha = 1.0;
#endif		 
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "<g,p> = " << ddot_(&mN, mGradient, &I1, mP, &I1) << std::endl;
			std::cout << "Line Search with a_max = " << alpha << "..." << std::endl;
		}

		
		
		// line search
		lineSearch(&alpha, aX, aF);
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			if (aX[i] < mLowerBound[i])
				aX[i] = mLowerBound[i];
			if (aX[i] > mUpperBound[i])
				aX[i] = mUpperBound[i];
		}
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "Step length found:" << alpha << std::endl;
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *aF);
		}
		
		// check convergence
		double df = f_prev - *aF;
		//double diff_x_norm = dnrm2_(&mN, mSk, &I1);
		convergenceReached =  fabs(df) < mAbsoluteError //&& diff_x_norm < mAbsoluteError)
							|| mStep >= mMaxIterations;
		
#if 0
		// gradient free components
		memcpy(mWorkSpaceVect, mGradient, mSizeVect);
		#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			if (  (fabs(aX[i]-mLowerBound[i])<1e-4 && mGradient[i] > 0.0)
				||(fabs(aX[i]-mUpperBound[i])<1e-4 && mGradient[i] < 0.0))
			{
				mWorkSpaceVect[i] = 0.0;
			}
		}
		double free_gradient_norm = dnrm2_(&mN, mWorkSpaceVect, &I1);
		std::cout << "|free gradient| = " << std::scientific << std::setprecision(12) << free_gradient_norm << std::endl;
		if (convergenceReached && free_gradient_norm > 1e-1)
		{
			std::cout << "Not converged, reinitialize the hessian." << std::endl;
			hessianInitialization();
			convergenceReached = false;
		}
#endif
		
		if (!convergenceReached)
		{
			// update the system
			
			computeGradient(aX, *aF, mGradient);
				
			memcpy(mSk, aX, mSizeVect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
			
			memcpy(mYk, mGradient, mSizeVect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "BFGS update..." << std::endl;
		
			// update the B matrix
			BFGSupdate();
		
			// update the active set
			double tolerance_active_set = 1e-4;
			activeSetUpdate(aX, tolerance_active_set);
		}
	}
#ifdef SCALE_OPT_VARIABLES
	unscaleVariables(aX);
#endif // SCALE_OPT_VARIABLES
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunction(const double *aX, bool aTrace)
{
	memcpy(&mXEvaluator[0], aX, mSizeVect);
#ifdef SCALE_OPT_VARIABLES
	unscaleVariables(&mXEvaluator[0]);
#endif // SCALE_OPT_VARIABLES
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if (mStopIfBigger && f >= mThreshold) throw FastCodeMLEarlyStopLRT();
	
	return -f;
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunctionForLineSearch(const double* aX, double aAlpha)
{
	memcpy(mWorkSpaceVect, aX, mSizeVect);
	daxpy_(&mN, &aAlpha, mP, &I1, mWorkSpaceVect, &I1);
	return evaluateFunction(mWorkSpaceVect, mTrace);
}


// ----------------------------------------------------------------------
void OptSQP::computeGradient(const double *aX, double aF0, double *aGrad)
{
	volatile double eh;
	double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], aX, mSizeVect);
	int i;
	double *delta = mWorkSpaceVect;
	
#ifdef SCALE_OPT_VARIABLES
	double *x = mWorkSpaceMat;
	memcpy(x, aX, mSizeVect);
	unscaleVariables(&mXEvaluator[0]);
	unscaleVariables(x);
#else
	const double *x = aX;
#endif // SCALE_OPT_VARIABLES
	
	// branch lengths
	for (i=0; i<mNumTimes; ++i)
	{
		eh = sqrt_eps * ( 1.0 + x[i] );
		if( x[i] + eh > mUpperBoundUnscaled[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		delta[i] = mXEvaluator[i] - x[i];
	}

	for (i=0; i<mNumTimes; ++i)
	{
		if (mActiveSet[i] == 0)
		{
			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-aF0)/delta[i];
		}
		// otherwise we don't change it
	}
	
	// other variables
	memcpy(&mXEvaluator[0], x, mSizeVect);
	for(; i<mN; ++i)
	{
		if (mActiveSet[i] == 0)
		{
			eh = sqrt_eps * ( 1.0 + fabs(x[i]) );
			if ( x[i] + eh > mUpperBoundUnscaled[i] )
				eh = -eh;
			mXEvaluator[i] += eh;
			eh = mXEvaluator[i] - x[i];

			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);

			aGrad[i] = (f-aF0)/eh;
			mXEvaluator[i] = x[i];
		}
	}
#ifdef SCALE_OPT_VARIABLES
	#pragma omp parallel for
	for (int j(0); j<mN; ++j)
	{
		double slb = mLowerBound[j];
		double sub = mUpperBound[j];
		double lb = mLowerBoundUnscaled[j];
		double ub = mUpperBoundUnscaled[j];
		
		aGrad[j] *= (ub-lb)/(sub-slb);
	}
#endif // SCALE_OPT_VARIABLES
}


// ----------------------------------------------------------------------
void OptSQP::hessianInitialization(void)
{
	// initialize hessian matrix to identity
	int n_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&n_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
	
#ifdef NON_IDENTITY_HESSIAN
	int i_;
	for (i_ = 0; i_<mNumTimes; ++i_)
	{
		mHessian[i_*diag_stride]  = 3.0;	// Htt
	}
	--i_;
	++i_; mHessian[i_*diag_stride] = 1.0; // Hv0v0
	++i_; mHessian[i_*diag_stride] = 2.0; // Hv1v1
	++i_; mHessian[i_*diag_stride] = 1.0; // Hww
	++i_; mHessian[i_*diag_stride] = 1.5; // Hkk
#endif // NON_IDENTITY_HESSIAN
	
	
#ifdef SCALE_OPT_VARIABLES
	// change the space of the hessian approximation representation
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		
		double scale = (ub-lb)/(sub-slb);
		scale = scale*scale;
		
		mHessian[i*diag_stride] *= scale;
	}		
#endif // SCALE_OPT_VARIABLES
}


// ----------------------------------------------------------------------
void OptSQP::BFGSupdate(void)
{
	// local variables
	double ys, sBs;
	double *Bs, *BssB, *yy;
	char trans = 'N';
	
	int n_sq = mN*mN;
	
	// compute vector B*mSk
	Bs = mWorkSpaceVect;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	sBs = ddot_(&mN, mSk, &I1, Bs,  &I1);
	ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
	
	bool reset_hessian = false;
#if 1	
	// Powell-SQP update:
	// change y so the matrix is positive definite
	double sigma = 0.2; // empirical value found by Powell
	double rho = ys / sBs;
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "ys = " << ys << std::endl;
	if (rho < sigma)
	{
		if (rho < 0.0)
		{
			reset_hessian = true;
		}
		else
		{
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "BFGS update leading to a non positive-definite matrix, performing Powell SQP:" << std::endl;
		
			double powell_factor = (1.0-sigma) / (1.0 - rho);
			dscal_(&mN, &powell_factor, mYk, &I1);
			powell_factor = 1.0 - powell_factor;
			daxpy_(&mN, &powell_factor, Bs, &I1, mYk, &I1);
		
			ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "ys = " << ys << " after Powell modification." << std::endl;
		}
	}
#endif

	if (reset_hessian)
	{
		hessianInitialization();
	}
	else
	{
		// compute Matrix B*mSk * mSk^T*B
		BssB = mWorkSpaceMat;
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			double prefactor = - Bs[i] / sBs;
			dcopy_(&mN, &Bs[0], &I1, &BssB[i*mN], &I1);
			dscal_(&mN, &prefactor, &BssB[i*mN], &I1);
		}
	
		// add the BssB / sBs contribution
		daxpy_(&n_sq, &D1, BssB, &I1, mHessian, &I1);
	
		// compute matrix y**T * y
		yy = mWorkSpaceMat;
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			double prefactor = mYk[i] / ys;
			dcopy_(&mN, &mYk[0], &I1, &yy[i*mN], &I1);
			dscal_(&mN, &prefactor, &yy[i*mN], &I1);
		}
	
		// add the yy / ys contribution
		daxpy_(&n_sq, &D1, yy, &I1, mHessian, &I1);
	
#if 1
		// make the diagonal more important in order to avoid non positive definite matrix, 
		// due to roundoff errors; this also speeds up the computation
		int diag_stride = mN+1;
		double factor = 1.0 + ((mN>30) ? 0.1:0.0); //1e-1*(1.0-1.0/static_cast<double>(mN));
		double inv_factor = 1.0/factor;
		dscal_(&n_sq, &inv_factor, mHessian, &I1);
		dscal_(&mN, &factor, mHessian, &diag_stride);
#endif
#if 0
		double *H = mWorkSpaceMat;
		memcpy(H, mHessian, mN*mSizeVect);
	
		double vl = -1e16;
		double vu = 0.1;
	
		double accuracy = 1e-8;
		int number_eigen_values;
		double *eigen_values = mWorkSpaceVect;
	
		std::vector<double> eigen_vectors(mN*mN);
		std::vector<int>	isuppz(2*mN);
		std::vector<double> work(1);
		std::vector<int> iwork(1);
		int lwork = -1;
		int liwork = -1;
		int info;
		int i1__not_used__, i2__not_used__;
		std::cout << "Workspace_querry:" << std::endl;
		dsyevr_("V", "V", "U"
			,&mN, H ,&mN
			,&vl, &vu
			,&i1__not_used__, &i2__not_used__
			,&accuracy, &number_eigen_values, eigen_values, &eigen_vectors[0], &mN
		    ,&isuppz[0], &work[0], &lwork, &iwork[0], &liwork, &info);
		lwork = static_cast<int>(work[0]);
		liwork = iwork[0];
		work.resize(lwork);
		iwork.resize(liwork);
	
		std::cout << "EigenValue solver:" << std::endl;
		dsyevr_("V", "V", "U"
			,&mN, H ,&mN
			,&vl, &vu
			,&i1__not_used__, &i2__not_used__
			,&accuracy, &number_eigen_values, eigen_values, &eigen_vectors[0], &mN
		    ,&isuppz[0], &work[0], &lwork, &iwork[0], &liwork, &info);
		
		std::cout << "Negative Eigen Values: " << std::endl;
		for (int i(0); i<number_eigen_values; ++i)
		{
			std::cout << std::setprecision(14) << eigen_values[i] << " ";
		}
		std::cout << std::endl;
    
		//double diag_to_add = (eigen_values[0] > 1e-2) ? 0.1 : 0.1 + eigen_values[0];
		//daxpy_(&mN, &diag_to_add, &D1, &I0, mHessian, &I1);
#endif
#if 0 // measure the condition number of the BFGS hessian approximation (experimental purpose)
		double *H = mWorkSpaceMat;
		memcpy(H, mHessian, mN*mSizeVect);
	
		double accuracy = 1e-8;
		int number_eigen_values;
		double *eigen_values = mWorkSpaceVect;
	
		std::vector<double> work(1);
		std::vector<int> iwork(1);
		int lwork = -1;
		int liwork = -1;
		int info;
		
		std::cout << "Workspace_querry:" << std::endl;
		
		dsyevd_("N", "U", &mN, H, &mN, eigen_values
               ,&work[0], &lwork, &iwork[0], &liwork, &info);
		
		lwork = static_cast<int>(work[0]);
		liwork = iwork[0];
		work.resize(lwork);
		iwork.resize(liwork);
	
		std::cout << "EigenValue solver:" << std::endl;
		dsyevd_("N", "U", &mN, H, &mN, eigen_values
               ,&work[0], &lwork, &iwork[0], &liwork, &info);
		
		std::cout << "Condition number = " << eigen_values[mN-1] / eigen_values[0] << std::endl;
#endif
	}
}


// ----------------------------------------------------------------------
void OptSQP::activeSetUpdate(const double *aX, const double aTolerance)
{
	// number of iterations where we skip the gradient computation (component i)
	const int max_count_lower = static_cast<const int>(1.3*log (static_cast<double>(mN)/10.)) + (mN>30 ? 1:0);
	const int max_count_upper = (mN > 30 ? 1 : 0);
	 
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		if (mActiveSet[i] > 0)
		{
			// reduce counters
			--mActiveSet[i];
		}
		else
		{
#ifdef SCALE_OPT_VARIABLES
			const double active_set_tol = aTolerance * (mUpperBound[i]-mLowerBound[i])/(mUpperBoundUnscaled[i]-mLowerBoundUnscaled[i]);
			const double y_tolerance = (mUpperBoundUnscaled[i]-mLowerBoundUnscaled[i])/(mUpperBound[i]-mLowerBound[i]) *1e-3*static_cast<double>(mN)/8.0; //aTolerance;
#else
			const double active_set_tol = aTolerance;
			const double y_tolerance = 1e-3*static_cast<double>(mN)/8.0; //aTolerance;
#endif // SCALE_OPT_VARIABLES
			if (fabs(mYk[i]) < y_tolerance)
			{			
				// update active set so we can reduce the gradient computation				
				if (aX[i] <= mLowerBound[i] + active_set_tol && mGradient[i] >= 0.0)
				{
					mActiveSet[i] = max_count_lower;
					if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
					std::cout << "\tVariable " << i << " in the (lower) active set.\n";
				}
				else if (aX[i] >= mUpperBound[i] - active_set_tol && mGradient[i] <= 0.0)
				{
					mActiveSet[i] = max_count_upper;
					if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
					std::cout << "\tVariable " << i << " in the (upper) active set.\n";
				}
			}
		}
	}
}

#ifndef STRONG_WOLFE_LINE_SEARCH
// ----------------------------------------------------------------------
void OptSQP::lineSearch(double *aAlpha, double *aX, double *aF)
{
	// constants for the Wolfe condition
	double c1 (2e-1);
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	double phi_0 = *aF, phi, phi_prev;
	double a_prev = 0.;
	double phi_a_prime;
	double a = *aAlpha;
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi_prev = phi_0;
	phi = evaluateFunctionForLineSearch(aX, a);
	
	double sigma, sigma_bas;
	int max_iter_back, max_iter_up;
	
	// we take a step dependant on the problem size:
	// if the problem is large, we are able to spend more time 
	// to find a better solution.
	// otherwise, we consider that the solution is sufficiently 
	// good and continue. 
	// The time of line search should be small compared to the 
	// gradient computation
	
	max_iter_back = max_iter_up = static_cast<int> (ceil( 3.*log(mN+10.) ));
	sigma_bas 	= pow(1e-5, 1./static_cast<double>(max_iter_back));
	
	
	// begin by a backtrace
	int iter = 0;
	while (phi > phi_0 + phi_0_prime*a*c1 && iter < max_iter_back)
	{
		++iter;
		a_prev = a;
		phi_prev = phi;
		sigma = sigma_bas * (0.9 + 0.2*randFrom0to1());
		a *= sigma;
		phi = evaluateFunctionForLineSearch(aX, a);
	}
	
	// compute the derivative
	double eh = sqrt(DBL_EPSILON);
	if ( a+eh >= 1.0 ) {eh = -eh;}
	phi_a_prime = (phi - evaluateFunctionForLineSearch(aX, a+eh))/eh;
	
	iter = 0;
	if (phi_a_prime < 0.0 && a != *aAlpha)
	{
		double a0 = a_prev;
		while (phi < phi_prev && iter < max_iter_back)
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			sigma = sigma_bas * (0.85 + 0.3*randFrom0to1());
			a = a + sigma*(a0-a);
			phi = evaluateFunctionForLineSearch(aX, a);
		}
		if (phi_prev < phi)
		{
			a = a_prev;
			phi = evaluateFunctionForLineSearch(aX, a);
		}
	}
	else
	{
		sigma_bas = 0.7;
		while (phi < phi_prev && iter < max_iter_up)
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			sigma = sigma_bas * (0.7 + 0.6*randFrom0to1());
			a *= sigma;
			phi = evaluateFunctionForLineSearch(aX, a);
		}
		if (phi_prev < phi)
		{
			a = a_prev;
			phi = evaluateFunctionForLineSearch(aX, a);
		}
	}
	
	
	*aF = phi;
	*aAlpha = a;
	daxpy_(&mN, aAlpha, mP, &I1, aX, &I1);
}

#else
// ----------------------------------------------------------------------
void OptSQP::lineSearch(double *aAlpha, double *aX, double *aF)
{
	// constants for the Wolfe condition
	// note that there must be 0 < c1 < c2 < 1
	const double c1 (2e-1), c2 (0.3);
	const double amax = *aAlpha;
	const double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	const double phi_0 = *aF;
	double phi, phi_prev;
	double a_prev = 0.0;
	
	double a = randFrom0to1();
		
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi = phi_prev = phi_0;
	
	int iter = 0;
	
	const int max_iter = static_cast<int> (ceil( 3.*log(mN+10.) ));
	
	while (iter < max_iter)
	{
		++iter;
		phi_prev = phi;
		phi = evaluateFunctionForLineSearch(aX, a);
		
		if (mVerbose >= VERBOSE_MORE_DEBUG)
		std::cout << "DEBUG LINE SEARCH: phi = " << phi << " for a = " << a << std::endl; 
		
		if (phi > phi_0 + c1*a*phi_0_prime || ((phi > phi_prev) && (iter > 1)) )
		{
			// solution between a_prev and a
			a = zoom(a_prev, a, aX, phi_0, phi_0_prime, phi_prev, c1, c2);
			break;
		}
		
		// compute the derivative at point a
		double eh = sqrt(DBL_EPSILON);
		if ( a+eh >= 1.0 ) {eh = -eh;}
		double phi_a_prime = (evaluateFunctionForLineSearch(aX, a+eh) - phi)/eh;
		
		if (fabs(phi_a_prime) <= -c2*phi_0_prime)
		{
			// wolfe conditions are satisfied, stop searching
			break;
		}
		if (phi_a_prime >= 0.0)
		{
			// solution between a and a_prev
			a = zoom(a, a_prev, aX, phi_0, phi_0_prime, phi, c1, c2);
			break;
		}
		double sigma = randFrom0to1();
		a_prev = a;
		a = amax + sigma*(a-amax);
	}
	
	
	*aF = evaluateFunctionForLineSearch(aX, a);
	*aAlpha = a;
	daxpy_(&mN, aAlpha, mP, &I1, aX, &I1);
}


// ----------------------------------------------------------------------
double OptSQP::zoom(double aAlo, double aAhi, const double *aX, const double& aPhi0, const double& aPhi0Prime, const double& aAphiLo, const double& c1, const double& c2)
{
	double a, phi_a_prime;
	double philo = aAphiLo;
	a = 0.5*(aAlo+aAhi);
	const double tolerance = 0.4 / static_cast<double>(mN);
	double phi = aPhi0;
	
	while (fabs(aAhi-aAlo) > tolerance)
	{
		double tmp = 0.5;
		a = tmp*aAlo + (1.-tmp)*aAhi;
		phi = evaluateFunctionForLineSearch(aX, a);
		if (mVerbose >= VERBOSE_MORE_DEBUG)
			std::cout << "DEBUG ZOOM: phi = " << phi << " for a = " << a << " alo: " << aAlo << " ahi: " << aAhi << " philo: " << philo  << std::endl;
		 
		if (phi > aPhi0 + a*c1*aPhi0Prime || phi > philo)
		{
			aAhi = a;
		}	
		else
		{
			double eh = sqrt(DBL_EPSILON);
			if ( a+eh >= 1.0 ) {eh = -eh;}
			phi_a_prime = (evaluateFunctionForLineSearch(aX, a+eh) - phi)/eh;
			
			if (fabs(phi_a_prime) <= -c2*aPhi0Prime)
				return a;
				
			if (phi_a_prime*(aAhi-aAlo) >= 0.0)
				aAhi = aAlo;
				
			aAlo = a;
			philo = phi;
		}
	}
	
	// make sure to have a small enough step (only if near 0)
	if (a <= tolerance && phi > aPhi0)
	{
		int step_decrease = 0;
		const int max_step_decrease = 10;
		while (step_decrease++ < max_step_decrease)
		{
			a *= 0.2;
			phi = evaluateFunctionForLineSearch(aX, a);
			if (phi < aPhi0) step_decrease = max_step_decrease+1;
		}
	}
	return a;
}
#endif // STRONG_WOLFE_LINE_SEARCH


