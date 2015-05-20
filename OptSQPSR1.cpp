
#include "OptSQPSR1.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>
#include <assert.h>


// ----------------------------------------------------------------------
//	Class members definition: OptSQP
// ----------------------------------------------------------------------
double OptSQPSR1::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	allocateMemory();
	
#ifdef SCALE_OPT_VARIABLES_SR1
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
	if(mN > i)
	{
		mLowerBound[i] = 0.0;
		mUpperBound[i] = 1.0;
	}
#endif // SCALE_OPT_VARIABLES_SR1
	
	double maxl = 1e7;
	SQPminimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------
void OptSQPSR1::allocateMemory(void)
{
	mSizeVect = mN*sizeof(double);
	
	mXEvaluator.resize(mN);
	mSpace.resize(2*mN*mN + mN*8);
	
	
	mWorkSpaceVect = &mSpace[0];
	mWorkSpaceMat = mWorkSpaceVect + mN;
	
	mGradient = mWorkSpaceMat + mN*mN;
	mP = mGradient + mN;
	mD = mP + mN;
	mHessian = mD + mN;
	
	mSk = mHessian + mN*mN;
	mYk = mSk + mN;
	
	mXPrev = mYk + mN;
	mGradPrev = mXPrev + mN;
	
	// set active set counters to zero
	mActiveSet.resize(mN, 0);
}


#ifdef SCALE_OPT_VARIABLES_SR1
// ----------------------------------------------------------------------
void OptSQPSR1::scaleVariables(double *x)
{
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		
		double x_ = slb + (x[i] - lb) * (sub-slb)/(ub-lb);
		
		if (x_ < slb)
			x_ = slb;
		if (x_ > sub)
			x_ = sub;
		x[i] = x_;
	}
}


// ----------------------------------------------------------------------
void OptSQPSR1::unscaleVariables(double *x)
{
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		double x_;
		
		// assume slb = sb = 0
		if (i<mNumTimes)
			x_ = x[i] * (ub-lb)/(sub-slb);
		else
			x_ = lb + (x[i] - slb) * ub/sub;
		
		if (x_ < lb)
			x_ = lb;
		if (x_ > ub)
			x_ = ub;
		x[i] = x_;
	}
}
#endif // SCALE_OPT_VARIABLES_SR1


// ----------------------------------------------------------------------
void OptSQPSR1::SQPminimizer(double *f, double *x)
{
#ifdef SCALE_OPT_VARIABLES_SR1
	scaleVariables(x);
#endif // SCALE_OPT_VARIABLES_SR1
	
	*f = evaluateFunction(x, mTrace);
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "Solution after Bootstrap:";
		mModel->printVar(mXEvaluator, *f);
	}
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	// bounds for the QP subproblem (l-x0 <= p <= u-x0)
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	mQPsolver.reset(new BOXCQP(mN, &localLowerBound[0], &localUpperBound[0]));	
	
	// initialize hessian matrix to identity
	int N_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&N_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
	
	
#ifdef SCALE_OPT_VARIABLES_SR1
	// change the space of the hessian approximation representation
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		
		double scale = (ub-lb)/(sub-slb);
		scale = scale*scale;
		
		mHessian[i*(mN+1)] *= scale;
	}		
#endif // SCALE_OPT_VARIABLES_SR1

	
	// main loop
	bool convergenceReached = false;
	for(mStep = 0; !convergenceReached; ++mStep)
	{
		// update local bounds
		memcpy(&localLowerBound[0], &mLowerBound[0], mSizeVect);
		memcpy(&localUpperBound[0], &mUpperBound[0], mSizeVect);
		daxpy_(&mN, &minus_one, x, &I1, &localLowerBound[0], &I1);
		daxpy_(&mN, &minus_one, x, &I1, &localUpperBound[0], &I1);
		
		
		// save current parameters
		double f_prev = *f;
		memcpy(mGradPrev, mGradient, mSizeVect);
		memcpy(mXPrev, x, mSizeVect);
		
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Quadratic program solving..." << std::endl;
		
		// solve quadratic program to get the search direction		
		solveUndefinedQP(&localLowerBound[0], &localUpperBound[0]);
		
		
		// line search
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		double alpha = 1.0;
		lineSearch(&alpha, x, f);
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			if (x[i] < mLowerBound[i])
				x[i] = mLowerBound[i];
			if (x[i] > mUpperBound[i])
				x[i] = mUpperBound[i];
		}
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "Step length found:" << alpha << std::endl;
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *f);
		}		
		
		// check convergence
		convergenceReached =   fabs(f_prev - *f) < mAbsoluteError
							|| mStep >= mMaxIterations;
		
		if (!convergenceReached)
		{
			// update the system
			
			computeGradient(x, *f, mGradient);
				
			memcpy(mSk, x, mSizeVect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
			memcpy(mYk, mGradient, mSizeVect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
		
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "SR1 update..." << std::endl;
			
		
			// update the B matrix
			SR1update();

		
			std::cout << "Hessian diagonal at step " << mStep << ":\n";
			for(int i=0; i<mN; ++i)
			{
				std::cout << mHessian[i*(mN+1)] << " ";
			}
			std::cout << std::endl;
		
			// update the active set
			//const int max_count_lower = (mN > 30 ? static_cast<const int>(log(static_cast<double>(mN))) : 1);
			const int max_count_lower = static_cast<const int>(1.3*log (static_cast<double>(mN)/10.)) + (mN>30 ? 1:0);
			const int max_count_upper = (mN > 30 ? 1 : 0);
	 
			#pragma omp parallel for
			for(int i=0; i<mN; ++i)
			{
				if (mActiveSet[i] > 0)
				{
					// reduce counters for active sets
					--mActiveSet[i];
				}
				else
				{
					const double active_set_tol = 1e-4 * (mUpperBound[i]-mLowerBound[i])/(mUpperBoundUnscaled[i]-mLowerBoundUnscaled[i]);
					// update active set so we can reduce the gradient computation				
					if (x[i] <= mLowerBound[i] + active_set_tol && mGradient[i] >= 0.0)
					{
						mActiveSet[i] = max_count_lower;
						if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
						std::cout << "Variable " << i << " in the (lower) active set.\n";
					}
					if (x[i] >= mUpperBound[i] - active_set_tol && mGradient[i] <= 0.0)
					{
						mActiveSet[i] = max_count_upper;
						if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
						std::cout << "Variable " << i << " in the (upper) active set.\n";
					}
				}
			}
		}
	}
#ifdef SCALE_OPT_VARIABLES_SR1
	unscaleVariables(x);
#endif // SCALE_OPT_VARIABLES_SR1
}


// ----------------------------------------------------------------------
double OptSQPSR1::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, mSizeVect);
#ifdef SCALE_OPT_VARIABLES_SR1
	unscaleVariables(&mXEvaluator[0]);
#endif // SCALE_OPT_VARIABLES_SR1
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if (mStopIfBigger && f >= mThreshold) throw FastCodeMLEarlyStopLRT();
	
	return -f;
}


// ----------------------------------------------------------------------
double OptSQPSR1::evaluateFunctionForLineSearch(const double* x, double alpha)
{
#if 0
	memcpy(mWorkSpaceVect, x, mSizeVect);
	daxpy_(&mN, &alpha, mP, &I1, mWorkSpaceVect, &I1);
#else
	memcpy(mWorkSpaceVect, x, mSizeVect);
	
	double a_2 = alpha*0.5;
	double a_sq_2 = square(alpha)*0.5;
	daxpy_(&mN, &a_2, mP, &I1, mWorkSpaceVect, &I1);
	daxpy_(&mN, &a_sq_2, mD, &I1, mWorkSpaceVect, &I1);
#endif
	return evaluateFunction(mWorkSpaceVect, mTrace);
}


// ----------------------------------------------------------------------
void OptSQPSR1::computeGradient(const double *x, double f0, double *aGrad)
{
	volatile double eh;
	double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], x, mSizeVect);
	size_t i;
	double *delta = mWorkSpaceVect;
	
	
#ifdef SCALE_OPT_VARIABLES_SR1
	double *x_ = mWorkSpaceMat;
	memcpy(x_, x, mSizeVect);
	unscaleVariables(&mXEvaluator[0]);
	unscaleVariables(x_);
#else
	const double *x_ = x;
#endif // SCALE_OPT_VARIABLES_SR1
	
	// branch lengths
	for(i=0; i<static_cast<size_t>(mNumTimes); ++i)
	{
		eh = sqrt_eps * ( 1.0 + x_[i] );
		if( x_[i] + eh > mUpperBoundUnscaled[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		delta[i] = mXEvaluator[i] - x_[i];
	}

	for(i=0; i<static_cast<size_t>(mNumTimes); ++i)
	{
		if(mActiveSet[i] == 0)
		{
			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-f0)/delta[i];
		}
		// otherwise we don't change it
	}
	
	// other variables
	memcpy(&mXEvaluator[0], x_, mSizeVect);
	for(; i<static_cast<size_t>(mN); ++i)
	{
		if (mActiveSet[i] == 0)
		{
			eh = sqrt_eps * ( 1.0 + fabs(x_[i]) );
			if ( x_[i] + eh > mUpperBoundUnscaled[i] )
				eh = -eh;
			mXEvaluator[i] += eh;
			eh = mXEvaluator[i] - x_[i];

			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-f0)/eh;
			mXEvaluator[i] = x_[i];
		}
	}
#ifdef SCALE_OPT_VARIABLES_SR1
	#pragma omp parallel for
	for (size_t j(0); j<mN; ++j)
	{
		double slb = mLowerBound[j];
		double sub = mUpperBound[j];
		double lb = mLowerBoundUnscaled[j];
		double ub = mUpperBoundUnscaled[j];
		
		aGrad[j] *= (ub-lb)/(sub-slb);
	}
#endif // SCALE_OPT_VARIABLES_SR1
}


// ----------------------------------------------------------------------
void OptSQPSR1::SR1update(void)
{
	// local variables
	double *v, *Bs;
	double vs;
	double *vvT;
	char trans = 'N';
	
	int mN_sq = mN*mN;
	
	// compute vector B*mSk
	Bs = mWorkSpaceVect;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	// compute vector v = -y + Bs
	v = Bs;
	daxpy_(&mN, &minus_one, mYk, &I1, v, &I1);
	
	vs = -ddot_(&mN, v, &I1, mSk, &I1);
	
	// only update if vs is big enough
	if (fabs(vs) > 1e-8)
	{
		double inverse_vs = 1.0/vs;
	
		// compute Matrix v.v^T / vs
		vvT = mWorkSpaceMat;
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			double prefactor = v[i] * inverse_vs;
			dcopy_(&mN, v, &I1, &vvT[i*mN], &I1);
			dscal_(&mN, &prefactor, &vvT[i*mN], &I1);
		}
	
		// add the v.v^T / vs contribution if it is not too big compared to B_prev
		
		double Frob_prev = dnrm2_(&mN_sq, mHessian, &I1);
		double Frob_modif = dnrm2_(&mN_sq, vvT, &I1);
		
		if (Frob_modif < 1e8 * (Frob_prev+1.0))
			daxpy_(&mN_sq, &D1, vvT, &I1, mHessian, &I1);
	}
}


#ifndef STRONG_WOLFE_LINE_SEARCH_SR1
// ----------------------------------------------------------------------
void OptSQPSR1::lineSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	double c1 (2e-4);
	double phi_0_prime = 0.5*ddot_(&mN, mP, &I1, mGradient, &I1);
	double phi_0 = *f, phi, phi_prev;
	double a_prev = 0.0;
	double phi_a_prime;
	double a = *aalpha;
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi_prev = phi_0;
	phi = evaluateFunctionForLineSearch(x, a);
	
	double sigma, sigma_bas;
	int maxIterBack, maxIterUp;
	
	// we take a step dependant on the problem size:
	// if the problem is large, we are able to spend more time 
	// to find a better solution.
	// otherwise, we consider that the solution is sufficiently 
	// good and continue. 
	// The time of line search should be small compared to the 
	// gradient computation
	
	maxIterBack = maxIterUp = static_cast<int> (ceil( 3.*log(mN+10.) ));
	sigma_bas 	= pow(1e-3, 1./static_cast<double>(maxIterBack));
	
	
	// begin by a backtrace
	size_t iter = 0;
	while(phi > phi_0 + phi_0_prime*a*c1 && iter < static_cast<size_t>(maxIterBack))
	{
		++iter;
		a_prev = a;
		phi_prev = phi;
		//sigma = 0.3+0.3*randFrom0to1();
		sigma = sigma_bas * (0.9 + 0.2*randFrom0to1());
		a *= sigma;
		phi = evaluateFunctionForLineSearch(x, a);
		std::cout << "\tDEBUG: a = " << a << ", phi = " << phi << std::endl; 
	}
	
	// compute the derivative
	double eh = sqrt(DBL_EPSILON);
	if ( a+eh >= 1.0 ) {eh = -eh;}
	phi_a_prime = (phi - evaluateFunctionForLineSearch(x, a+eh))/eh;
	
	iter = 0;
	if (phi_a_prime < 0.0 && fabs(a - *aalpha) > 1e-8)
	{
		double a0 = a_prev;
		while(phi < phi_prev && iter < static_cast<size_t>(maxIterBack))
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			//sigma = 0.3+0.4*randFrom0to1();
			sigma = sigma_bas * (0.85 + 0.3*randFrom0to1());
			a = a + sigma*(a0-a);
			phi = evaluateFunctionForLineSearch(x, a);
		}
		if (phi_prev < phi)
		{
			a = a_prev;
			phi = evaluateFunctionForLineSearch(x, a);
		}
	}
	else
	{
		sigma_bas = 0.7;
		while(phi < phi_prev && iter < static_cast<size_t>(maxIterUp))
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			//sigma = 0.5+0.5*sqrt(randFrom0to1());
			sigma = sigma_bas * (0.7 + 0.6*randFrom0to1());
			a *= sigma;
			phi = evaluateFunctionForLineSearch(x, a);
		}
		if (phi_prev < phi)
		{
			a = a_prev;
			phi = evaluateFunctionForLineSearch(x, a);
		}
	}
	
	
	*f = phi;
	*aalpha = a;
	daxpy_(&mN, aalpha, mP, &I1, x, &I1);
	a = a*a*0.5;
	daxpy_(&mN, &a, mD, &I1, x, &I1);
}

#else
// ----------------------------------------------------------------------
void OptSQPSR1::lineSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	// note that there must be 0 < c1 < c2 < 1
	const double c1 (2e-1), c2 (0.3);
	const double amax = *aalpha;
	const double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	const double phi_0 = *f;
	double phi, phi_prev;
	double a_prev = 0.0;
	double phi_a_prime;
	
	double a = randFrom0to1();
		
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi = phi_prev = phi_0;
	
	double sigma, sigma_bas;
	int iter = 0;
	
	const int max_iter = static_cast<int> (ceil( 3.*log(mN+10) ));
	sigma_bas = pow(1e-2, 1./static_cast<double>(max_iter));
	
	while(iter < max_iter)
	{
		++iter;
		phi_prev = phi;
		phi = evaluateFunctionForLineSearch(x, a);
		
		if (mVerbose >= VERBOSE_MORE_DEBUG)
		std::cout << "DEBUG LINE SEARCH: phi = " << phi << " for a = " << a << std::endl; 
		
		if (phi > phi_0 + c1*a*phi_0_prime || ((phi > phi_prev) && (iter > 1)) )
		{
			// solution between a_prev and a
			a = zoom(a_prev, a, x, phi_0, phi_0_prime, phi_prev, c1, c2);
			break;
		}
		
		// compute the derivative at point a
		double eh = sqrt(DBL_EPSILON);
		if ( a+eh >= 1.0 ) {eh = -eh;}
		phi_a_prime = (evaluateFunctionForLineSearch(x, a+eh) - phi)/eh;
		
		if (fabs(phi_a_prime) <= -c2*phi_0_prime)
		{
			// wolfe conditions are satisfied, stop searching
			break;
		}
		if (phi_a_prime >= 0.0)
		{
			// solution between a and a_prev
			a = zoom(a, a_prev, x, phi_0, phi_0_prime, phi, c1, c2);
			break;
		}
		//sigma = sigma_bas * (0.9 + 0.2*randFrom0to1());
		sigma = randFrom0to1();
		//sigma = (sigma>0.5) ? square(sigma) : sqrt(sigma);
		a_prev = a;
		a = amax + sigma*(a-amax);
	}
	
	
	*f = evaluateFunctionForLineSearch(x,a);
	*aalpha = a;
	daxpy_(&mN, aalpha, mP, &I1, x, &I1);
}


// ----------------------------------------------------------------------
double OptSQPSR1::zoom(double alo, double ahi, double *x, const double& phi_0, const double& phi_0_prime, const double& aphi_lo, const double& c1, const double& c2)
{
	double a, phi, phi_a_prime;
	double philo = aphi_lo;
	a = 0.5*(alo+ahi);
	while( fabs(ahi-alo) > 0.01 )
	{
		double tmp = 0.5;
		a = tmp*alo + (1.-tmp)*ahi;
		phi = evaluateFunctionForLineSearch(x, a);
		if (mVerbose >= VERBOSE_MORE_DEBUG)
		std::cout << "DEBUG ZOOM: phi = " << phi << " for a = " << a << " alo: " << alo << " ahi: " << ahi << " philo: " << philo  << std::endl;
		 
		if (phi > phi_0 + a*c1*phi_0_prime || phi > philo)
		{
			ahi = a;
		}	
		else
		{
			double eh = sqrt(DBL_EPSILON);
			if ( a+eh >= 1.0 ) {eh = -eh;}
			phi_a_prime = (evaluateFunctionForLineSearch(x, a+eh) - phi)/eh;
			
			if (fabs(phi_a_prime) <= -c2*phi_0_prime)
				return a;
				
			if (phi_a_prime*(ahi-alo) >= 0.0)
				ahi = alo;
				
			alo = a;
			philo = phi;
		}
	}
	
	return a;
}
#endif // STRONG_WOLFE_LINE_SEARCH_SR1


// ----------------------------------------------------------------------
void OptSQPSR1::solveUndefinedQP(const double *localLowerBound, const double *localUpperBound)
{
	int mN_sq = mN*mN;
	// create a working copy of the hessian
	double *eigenVectors = mWorkSpaceMat;
	memcpy(eigenVectors, mHessian, mN*mSizeVect);
	
	// find the eigenvalues and eigenvectors of B
	const char jobz = 'V';
	const char uplo = 'U';
	std::vector<double> eigenValues_(mN);
	double *eigenValues = &eigenValues_[0];
	
	const int lwork = 1 + 6*mN + 2*mN*mN;
	std::vector<double> work(lwork);
	
	const int liwork = 3+5*mN;
	std::vector<int> iwork(liwork);
	
	int info;
	
	dsyevd_(&jobz
		   ,&uplo
		   ,&mN
		   ,eigenVectors
		   ,&mN
		   ,eigenValues
		   ,&work[0]
		   ,&lwork
		   ,&iwork[0]
		   ,&liwork
		   ,&info);
	
	assert(info == 0);
	
	// find lower index M such that the M first eigenvalues are not strictly positive
	// M=-1 means no eigenvalue is not strictly positive
	int M = -1;
	while (M < mN && eigenValues[M+1] <= 0.0) {++M;}
	int number_positive_eigenvalues = mN - M - 1;
	
	std::cout << "\t\t" << M << " negative eigen values, " << number_positive_eigenvalues << " strictly positive." << std::endl;
	
	// -- solve the QP in convex part
	
	if (M < mN-1)
	{
		// form the modified hessian matrix
		double *lambda_p_S = &work[0];
		double *convex_hessian = lambda_p_S + mN_sq;
		memcpy(lambda_p_S, eigenVectors, mN*mSizeVect);
		
		double *gradient = convex_hessian + mN_sq;
		memcpy(gradient, mGradient, mSizeVect);
		//dcopy_(&mN, &D0, &I0, gradient, &I1);
		//dcopy_(&number_positive_eigenvalues, &mGradient[M+1], &I1, &gradient[M+1], &I1);

		//#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			double prefactor = (i <= M) ? 0.0001 : eigenValues[i];
			dscal_(&mN, &prefactor, &lambda_p_S[i*mN], &I1);
			
			char trans = 'T';
			dgemv_(&trans, &mN, &mN, &D1, eigenVectors, &mN, &lambda_p_S[i*mN], &I1, &D0, &convex_hessian[i*mN], &I1);
		}
		
		// use BOXQP to solve the sub problem
		bool solution_on_border;
		mQPsolver->solveQP(convex_hessian, gradient, &mN, mP, &solution_on_border);
	}
	else
	{
		dcopy_(&mN, &D0, &I0, mP, &I1);
	}
	
	
	dcopy_(&mN, &D0, &I0, mD, &I1);
	// -- weighted sum of the the negative curvatures directions
	if (M != -1)
	{
		double *negative_curv_direction = &work[0];
		dcopy_(&mN, &D0, &I0, negative_curv_direction, &I1);
		for (int i=0; i<=M; ++i)
		{
			// - |lambdai| * <g,Si>
			const double proportion = - ddot_(&mN, mGradient, &I1, &eigenVectors[i*mN], &I1);
			std::cout << "\tProportion for eigen value " << i << "(" << eigenValues[i] << ") : " << proportion << std::endl;
			daxpy_(&mN, &proportion, &eigenVectors[i*mN], &I1, negative_curv_direction, &I1);
		}
	
		// find largest a such that l <= a*mD <= u
		double a = 1e16;
		for (int i=0; i<mN; ++i)
		{
			double di = negative_curv_direction[i];
			double maxa;
			if (di < 0.0)
			{
				maxa = localLowerBound[i] / di;
			}
			else
			{
				maxa = localUpperBound[i] / di;
			}
			a = min2(a, maxa);
		}
		std::cout << "\tNegative space: alpha found = " << a << std::endl;
		// update mD
		daxpy_(&mN, &a, negative_curv_direction, &I1, mD, &I1);
	}
}


