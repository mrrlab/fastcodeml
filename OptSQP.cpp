
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
	// allocate workspace
	mN = static_cast<int>(aVars.size());
	allocateMemory();
	// compute maximum log-likelihood
	double minimum_opposite_likelihood = 1e7;
	SQPminimizer(&minimum_opposite_likelihood, &aVars[0]);
	return -minimum_opposite_likelihood;
}


// ----------------------------------------------------------------------
void OptSQP::allocateMemory(void)
{
	mSizeVect = mN*sizeof(double);
	
	// allocate the double precision workspace
	mXEvaluator.resize(mN);
	mSpace.resize(2*mN*mN + mN*7);

	// assign working arrays to allocated space
	mWorkSpaceVect	= &mSpace[0];
	mWorkSpaceMat	= mWorkSpaceVect	+ mN;
	mGradient		= mWorkSpaceMat		+ mN*mN;
	mP				= mGradient			+ mN;
	mHessian		= mP				+ mN;
	mSk				= mHessian			+ mN*mN;
	mYk				= mSk				+ mN;
	mXPrev			= mYk				+ mN;
	mGradPrev		= mXPrev			+ mN;
	
	// set active set counters to zero
	mActiveSet.resize(mN, 0);
}


// ----------------------------------------------------------------------
void OptSQP::SQPminimizer(double *aF, double *aX)
{
	*aF = evaluateFunction(aX, mTrace);
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "Initial point for SQP:";
		mModel->printVar(mXEvaluator, *aF);
	}
	
	// compute gradient at x0
	computeGradient(aX, *aF, mGradient);
	
	// initialize the hessian matrix to identity
	hessianInitialization();
	
	// allocate local bounds for the QP subproblem
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
		
		// solve quadratic program to get the search direction	
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Solving quadratic program" << std::endl;
		bool QPsolutionOnBorder;
		bool QP_converged = mQPsolver->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder, NULL);
		
		// take the projected gradient direction in case of unsuccesfull QP solution
		// unsuccesful QP probably due to bad hessian matrix (e.g. not positive definite) -> reinitialize it
		if (!QP_converged) 
		{
			hessianInitialization();
			memcpy(mP, mGradient, mSizeVect);
			dscal_(&mN, &minus_one, mP, &I1);
			#pragma omp parallel for
			for (int i=0; i<mN; ++i)
			{
				mP[i] = max2(mP[i], localLowerBound[i]);
				mP[i] = min2(mP[i], localUpperBound[i]);
			}
		}
	
		// extend the upper limit of line search to the boundaries
		double alpha = 1e16;
		for (int i=0; i<mN; ++i)
		{
			double l = localLowerBound[i];
			double u = localUpperBound[i];
			double p = mP[i];
			if (p < -1e-8)
				alpha = min2(alpha, l/p);
			else if (p > 1e-8)
				alpha = min2(alpha, u/p);
		}
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "Line Search with a_max = " << alpha << "..." << std::endl;
		}
		
		// perform line search along mP
		lineSearch(&alpha, aX, aF);
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			double x_ = aX[i];
			x_ = min2(x_, mUpperBound[i]);
			x_ = max2(x_, mLowerBound[i]);
			aX[i] = x_;
		}
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "Step length found:" << alpha << std::endl;
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *aF);
		}
		
		// check convergence
		const double df = f_prev - *aF;
		
		if (mHimmelblauTermination)
		{
			const double diff_x_norm = fabs(mSk[idamax_(&mN, mSk, &I1)]);
			convergenceReached = (fabs(df) < mAbsoluteError && diff_x_norm < mAbsoluteError);
			convergenceReached = convergenceReached || mStep >= mMaxIterations;
		}
		else
		{
			convergenceReached =  fabs(df) < mAbsoluteError	|| mStep >= mMaxIterations;
		}
		
		// update variables before next iteration
		if (!convergenceReached)
		{
			// update gradient			
			computeGradient(aX, *aF, mGradient);
			
			// compute s = x - x_prev
			memcpy(mSk, aX, mSizeVect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
			
			// compute y = g - g_prev
			memcpy(mYk, mGradient, mSizeVect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
			
			// update the B matrix
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "BFGS update..." << std::endl;
			BFGSupdate();
		
			// update the active set
			const double tolerance_active_set = 1e-4;
			activeSetUpdate(aX, tolerance_active_set);
		}
	}
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunction(const double *aX, bool aTrace)
{
	memcpy(&mXEvaluator[0], aX, mSizeVect);
	const double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
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
	const double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], aX, mSizeVect);
	int i;
	double *delta = mWorkSpaceVect;
	
	const double *x = aX;
	
	// branch lengths
	for (i=0; i<mNumTimes; ++i)
	{
		eh = sqrt_eps * ( 1.0 + x[i] );
		if( x[i] + eh > mUpperBound[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		delta[i] = mXEvaluator[i] - x[i];
	}
	for (i=0; i<mNumTimes; ++i)
	{
		// compute this component if we know it can change
		if (mActiveSet[i] == 0)
		{
			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-aF0)/delta[i];
		}
	}
	
	// other parameters
	memcpy(&mXEvaluator[0], x, mSizeVect);
	for(; i<mN; ++i)
	{
		if (mActiveSet[i] == 0)
		{
			eh = sqrt_eps * ( 1.0 + fabs(x[i]) );
			if ( x[i] + eh > mUpperBound[i] )
				eh = -eh;
			mXEvaluator[i] += eh;
			eh = mXEvaluator[i] - x[i];

			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-aF0)/eh;
			mXEvaluator[i] = x[i];
		}
	}
}


// ----------------------------------------------------------------------
void OptSQP::hessianInitialization(void)
{
	// initialize hessian matrix to identity
	const int n_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&n_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
}


// ----------------------------------------------------------------------
void OptSQP::BFGSupdate(void)
{
	const char not_transposed = 'N';
	const int n_sq = mN*mN;
	
	// compute vector B*mSk
	double *Bs = mWorkSpaceVect;
	dgemv_(&not_transposed, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	const double sBs = ddot_(&mN, mSk, &I1, Bs,  &I1);
	double ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
	
	bool reset_hessian = false;

	// Powell-SQP update:
	// change y so the matrix is positive definite
	const double sigma = 0.2; // empirical value found by Powell
	const double rho = ys / sBs;
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "ys = " << std::scientific << ys << std::fixed << std::endl;
	}
	if (rho < sigma)
	{
		if (rho < 0.0)
		{
			reset_hessian = true;
		}
		else
		{
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "BFGS update leading to a non positive-definite matrix, performing Powell-SQP update:" << std::endl;
		
			double powell_factor = (1.0-sigma) / (1.0 - rho);
			dscal_(&mN, &powell_factor, mYk, &I1);
			powell_factor = 1.0 - powell_factor;
			daxpy_(&mN, &powell_factor, Bs, &I1, mYk, &I1);
		
			ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			{
				std::cout << "ys = " << std::scientific  << " after Powell modification." << std::endl;
				std::cout << std::fixed;
			}
		}
	}

	if (reset_hessian)
	{
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "\tReset the hessian matrix approximation." << std::endl;
		hessianInitialization();
	}
	else
	{
		// compute Matrix B*mSk * mSk^T*B
		double *BssB = mWorkSpaceMat;
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			const double prefactor = - Bs[i] / sBs;
			dcopy_(&mN, &Bs[0], &I1, &BssB[i*mN], &I1);
			dscal_(&mN, &prefactor, &BssB[i*mN], &I1);
		}
		
		// add the BssB / sBs contribution
		daxpy_(&n_sq, &D1, BssB, &I1, mHessian, &I1);
	
		// compute matrix y**T * y
		double *yy = mWorkSpaceMat;
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			const double prefactor = mYk[i] / ys;
			dcopy_(&mN, &mYk[0], &I1, &yy[i*mN], &I1);
			dscal_(&mN, &prefactor, &yy[i*mN], &I1);
		}
	
		// add the yy / ys contribution
		daxpy_(&n_sq, &D1, yy, &I1, mHessian, &I1);

		// measure the condition number of the BFGS hessian approximation and modify the hessian if badly conditioned
		
		// estimate the condition number
		double *hessian = mWorkSpaceMat;
		memcpy(hessian, mHessian, mN*mSizeVect);
		
		// hessian 1-norm
		double norm_hessian = 0.0;
		for (int i(0); i<mN; ++i)
		{
			const double sum_current_column = dasum_(&mN, hessian+i*mN, &I1);
			norm_hessian = max2(norm_hessian, sum_current_column); 
		}
		
		// LU decomposition of the hessian (needed for condition number estimate)
		const char uplo = 'U';
		int info;
		dpotrf_(&uplo, &mN, hessian, &mN, &info);
		if (info != 0)
			std::cerr << "\tError during the Choleski decomposition of hessian. INFO = " << info << std::endl;
		
		// compute condition number
		std::vector<double> work(3*mN);
		std::vector<int> iwork(mN);
		double reciprocal_condition_number;
		dpocon_(&uplo, &mN, hessian, &mN, &norm_hessian, &reciprocal_condition_number, &work[0], &iwork[0], &info);
		if (info != 0)
			std::cerr << "\tError during the condition number estimation of hessian. INFO = " << info << std::endl;
		double condition_number = 1.0 / reciprocal_condition_number;
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Condition number before = " << condition_number << std::endl;

		// make the diagonal more important in order to avoid non positive definite matrix, 
		// due to roundoff errors; 
		// this also speeds up the computation as the condition number is reduced by a factor of 10^3 in some cases!
		if (condition_number > 1e3)
		{
			// values obtained from experiments
			const double off_diagonal_scaling = min2(1./1.1, 15.35 / (log(condition_number) + 7.67));
			double *diagonal_backup = mWorkSpaceVect;
			const int diagonal_stride = mN+1;
			dcopy_(&mN, mHessian, &diagonal_stride, diagonal_backup, &I1);
			dscal_(&n_sq, &off_diagonal_scaling, mHessian, &I1);
			dcopy_(&mN, diagonal_backup, &I1, mHessian, &diagonal_stride);
		}
	}
}


// ----------------------------------------------------------------------
void OptSQP::activeSetUpdate(double *aX, const double aTolerance)
{
	// number of iterations where we skip the gradient computation (component i)
	const int max_count_lower = static_cast<const int>(1.3*log (static_cast<double>(mN)/10.)) + 1;
	const int max_count_upper = (mN > 30 ? 1 : 0);
	 
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		if (mActiveSet[i] > 0)
		{
			// reduce counters
			--mActiveSet[i];
			mSk[i] = 0.0;
		}
		else
		{
			const double active_set_tol = aTolerance;
			const double y_tolerance = min2(1e-3*static_cast<double>(mN)/8.0, 1e-1);
			
			if (fabs(mYk[i]) < y_tolerance)
			{
				// update active set so we can reduce the gradient computation				
				if (aX[i] <= mLowerBound[i] + active_set_tol && mGradient[i] >= 0.0)
				{
					// put the variable in the active set
					mActiveSet[i] = max_count_lower;
					
					// make it equal to the boundary so mSk[i] is null at next iteration
					// this should avoid badly conditioned matrix updates
					aX[i] = mLowerBound[i];
					
					if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
					std::cout << "\tVariable " << i << " in the (lower) active set.\n";
				}
				else if (aX[i] >= mUpperBound[i] - active_set_tol && mGradient[i] <= 0.0)
				{
					// put the variable in the active set
					mActiveSet[i] = max_count_upper;
					
					// make it equal to the boundary so mSk[i] is null at next iteration
					// this should avoid badly conditioned matrix updates
					aX[i] = mUpperBound[i];
					
					if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
					std::cout << "\tVariable " << i << " in the (upper) active set.\n";
				}
				else if (fabs(mSk[i]) < active_set_tol && fabs(mGradient[i]) < active_set_tol)
				{
					// put the variable in the active set
					mActiveSet[i] = max_count_upper;
					if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
					std::cout << "\tVariable " << i << " not moving: put it in the active set.\n";
				}
			}
		}
	}
}

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



