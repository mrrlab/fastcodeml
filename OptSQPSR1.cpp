
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
	mSpace.resize(3*mN*mN + mN*8);
	
	
	mWorkSpaceVect = &mSpace[0];
	mWorkSpaceMat = mWorkSpaceVect + mN;
	
	mGradient = mWorkSpaceMat + mN*mN;
	mP = mGradient + mN;
	mHessian = mP + mN;
	mInverseHessian = mHessian + mN*mN;
	
	mSk = mInverseHessian + mN*mN;
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
	
	dcopy_(&N_sq, &D0, &I0, mInverseHessian, &I1);
	dcopy_(&mN, &D1, &I0, mInverseHessian, &diag_stride);
	
	
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
		mHessianInverse[i*(mN+1)] /= scale;
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
			std::cout << "Computing new search direction..." << std::endl;
		
		// solve quadratic program to get the search direction		
		computeSearchDirection(x);
		
		
		// line search
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		// extend the limits to the boundaries
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
			
		
			// update the Hessian informations
			SR1update();
		
#if 0
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
#endif
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
	memcpy(mWorkSpaceVect, x, mSizeVect);
	daxpy_(&mN, &alpha, mP, &I1, mWorkSpaceVect, &I1);

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
	const int mN_sq = mN*mN;
	const double eps1 = 1e-6;
	const double eps2 = 1e-6;
	char trans = 'N';
	
	// compute vector Bs
	double *Bs = mWorkSpaceVect;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	// compute vector v = y - Bs
	double *v = Bs;
	daxpy_(&mN, &minus_one, mYk, &I1, v, &I1);
	dscal_(&mN, &minus_one, v, &I1);
	
	const double vs = ddot_(&mN, v, &I1, mSk, &I1);
	const double threshold_abs_vs = eps1 * dnrm2_(&mN, mSk, &I1) * dnrm2_(&mN, v, &I1);
	
	// compute vector Hy
	double *Hy = mP; // unused space
	dgemv_(&trans, &mN, &mN, &D1, mInverseHessian, &mN, mYk, &I1, &D0, Hy, &I1); 
	
	// compute vector w = s - Hy
	double *w = Hy;
	daxpy_(&mN, &minus_one, mSk, &I1, w, &I1);
	dscal_(&mN, &minus_one, v, &I1);
	
	const double wy = ddot_(&mN, w, &I1, mYk, &I1);
	const double threshold_abs_wy = eps2 * dnrm2_(&mN, mYk, &I1) * dnrm2_(&mN, w, &I1);
	
	// only update if well defined
	if (fabs(vs) > threshold_abs_vs && fabs(wy) > threshold_abs_wy)
	{
		// --- Hessian update
		double inverse_vs = 1.0/vs;
	
		// compute Matrix v.v^T / vs
		double *vvT = mWorkSpaceMat;
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			double prefactor = v[i] * inverse_vs;
			dcopy_(&mN, v, &I1, &vvT[i*mN], &I1);
			dscal_(&mN, &prefactor, &vvT[i*mN], &I1);
		}
		daxpy_(&mN_sq, &D1, vvT, &I1, mHessian, &I1);
		
		// --- Inverse Hessian update
		double inverse_wy = 1.0/wy;
	
		// compute Matrix w.w^T / wy
		double *wwT = mWorkSpaceMat;
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			double prefactor = w[i] * inverse_wy;
			dcopy_(&mN, w, &I1, &wwT[i*mN], &I1);
			dscal_(&mN, &prefactor, &wwT[i*mN], &I1);
		}
		daxpy_(&mN_sq, &D1, wwT, &I1, mInverseHessian, &I1);
	}
	else
	{
		std::cout << "\tSkipping SR1 update" << std::cout;
	}
}


#ifndef STRONG_WOLFE_LINE_SEARCH_SR1
// ----------------------------------------------------------------------
void OptSQPSR1::lineSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	double c1 (2e-4);
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
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
// TODO: use mQPsolver

void OptSQPSR1::computeSearchDirection(const double *aX)
{
	int mN_sq = mN*mN;
	const double tolerance_min_eigenvalue = -1e-3;
	const double tolerance_QP = 1e-2;
	bool use_projected_gradient (false);
	
	// compute the projected gradient
	double *projected_gradient = mWorkSpaceVect;
	memcpy(projected_gradient, mGradient, mSizeVect);
	dscal_(&mN, &minus_one, projected_gradient, &I1);
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		double p = projected_gradient[i];
		double l = mLowerBound[i] - aX[i];
		double u = mUpperBound[i] - aX[i];
		p = min2(max2(p, l), u);
		projected_gradient[i] = p;
	}
	
	// compute minimum eigenvalue and corresponding direction
	char V('V'), I('I'), U('U');		
	const int il = 1;
	const int iu = 1;
	const double abs_tol_eigenvalues = 1e-8;
	int num_eigenvalues_found;
	double *min_eigenvalue_direction = mWorkSpaceMat;
	double *eigenvalues = min_eigenvalue_direction + mN;
	
	std::vector<double> work(1);
	std::vector<int> iwork(1);
	int lwork, liwork, info;
	
	// workspace querry
	lwork = -1; liwork = -1;
	dsyevr_(&V, &I, &U, &mN, mHessian, &mN, NULL, NULL, 
       	&il, &iu, &abs_tol_eigenvalues, &num_eigenvalues_found,
       	eigenvalues, min_eigenvalue_direction, &mN, NULL,
		&work[0], &lwork, &iwork[0], &liwork, &info);
	
	lwork = static_cast<int> (work[0]);
	liwork = iwork[0];
	work.resize(lwork);
	iwork.resize(liwork);
	
	// compute eigenvalue/eigenvector
	dsyevr_(&V, &I, &U, &mN, mHessian, &mN, NULL, NULL, 
       	&il, &iu, &abs_tol_eigenvalues, &num_eigenvalues_found,
       	eigenvalues, min_eigenvalue_direction, &mN, NULL,
		&work[0], &lwork, &iwork[0], &liwork, &info);
	
	
	double lambda_min = eigenvalues[0];
	std::cout << "\tLambda_min = " << lambda_min << std::endl;
	
	if (lambda_min > tolerance_min_eigenvalue)
	{
		double *hessian_matrix;
		if (lambda_min < tolerance_QP)
		{
			const double diagonal = tolerance_QP - lambda_min + 1.0;
			const int diagonal_stride = mN+1;
			hessian_matrix = mWorkSpaceMat;
			memcpy(hessian_matrix, mHessian, mN_sq*sizeof(double));
			daxpy_(&mN, &D1, &diagonal, &I0, hessian_matrix, &diagonal_stride);
		}
		else
		{
			hessian_matrix = mHessian;
		}
		// consider the matrix to be positive definite
		// solve Quadratic Program
		bool QP_solution_on_border;
		bool QP_converged = mQPsolver->solveQP(hessian_matrix, mGradient, &mN, mP, &QP_solution_on_border, NULL);
		if (!QP_converged)
		{
			use_projected_gradient = true;
		}
		else
		{
			std::cout << "\tQP direction chosen." << std::endl;
		}
	}
	else
	{
		// choose the negative curvature direction
		memcpy(mP, min_eigenvalue_direction, mSizeVect);
		// scale it
		double scale = max2(dnrm2_(&mN, mSk, &I1), 1e-1);		
		if (ddot_(&mN, mP, &I1, mGradient, &I1) > 0.0)
		{
			scale = -scale;
		}
		dscal_(&mN, &scale, mP, &I1);
		
		// take its projection
		#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			double p = mP[i];
			double l = mLowerBound[i] - aX[i];
			double u = mUpperBound[i] - aX[i];
			p = min2(max2(p, l), u);
			mP[i] = p;
		}
		
		scale = dnrm2_(&mN, mP, &I1);
		
		double NC_g = ddot_(&mN, mP, &I1, mGradient, &I1);
		double grad_norm = dnrm2_(&mN, mGradient, &I1);
		std::cout << "\t angle: " << NC_g/grad_norm/scale << std::endl;
		if (NC_g/grad_norm/scale > -1e-1) 
		{
			use_projected_gradient = true;
		}
		else
		{
			std::cout << "\tProjected negative curvature direction chosen." << std::endl;
		}
	}
	
	if (use_projected_gradient)
	{
		memcpy(mP, projected_gradient, mSizeVect);
		std::cout << "\tProjected gradient chosen." << std::endl;
	}
}


#if 0 // backup: alg. http://research.sabanciuniv.edu/15897/1/sr1nc.pdf
void OptSQPSR1::computeSearchDirection(const double *aX)
{
	int mN_sq = mN*mN;
	
	const double tau = 2.0;
	const double epsilon = 1e-6;
	
	// compute the active gradient (opposite direction)
	double *active_gradient = mWorkSpaceVect;
	memcpy(active_gradient, mGradient, mSizeVect);
	
	const double tol = 1e-4;
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		double ag = -active_gradient[i];
		double x = aX[i];
		double l = mLowerBound[i];
		double u = mUpperBound[i];
		if (x-l < tol)
			ag = max2(ag, 0.0);
		else if(u-x < tol)
			ag = min2(ag, 0.0);
					
		active_gradient[i] = ag;
	}
	
	// compute quasi Newton direction
	double *newton_direction = mWorkSpaceMat;
	char trans = 'N';
	dgemv_(&trans, &mN, &mN, &minus_one, mInverseHessian, &mN, mGradient, &I1, &D0, newton_direction, &I1);
	
	// keep only the active part
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		double d = newton_direction[i];
		double x = aX[i];
		double l = mLowerBound[i];
		double u = mUpperBound[i];
		
		if (x-l < tol)
			d = max2(d, 0.0);
		else if(u-x < tol)
			d = min2(d, 0.0);
			
		newton_direction[i] = d;
	}
	
	// compute negative curvature direction in case of not p.d. matrix (suspected)
	double y_s = ddot_(&mN, mSk, &I1, mYk, &I1);
	double g_nd = ddot_(&mN, newton_direction, &I1, mGradient, &I1);
	double *negative_curv_direction = newton_direction+mN;
	
	if (g_nd > 0.0 || y_s < 0.0)
	{
		// compute negative eigenvalue and negative curvature direction
		char V('V'), I('I'), U('U');		
		const int il = 1;
		const int iu = 1;
		const double abs_tol_eigenvalues = 1e-8;
		int num_eigenvalues_found;
		double *eigenvalues = negative_curv_direction+mN;
		
		std::vector<double> work(1);
		std::vector<int> iwork(1);
		int lwork, liwork, info;
		
		// workspace querry
		lwork = -1; liwork = -1;
		dsyevr_(&V, &I, &U, &mN, mHessian, &mN, NULL, NULL, 
        	&il, &iu, &abs_tol_eigenvalues, &num_eigenvalues_found,
        	eigenvalues, negative_curv_direction, &mN, NULL,
			&work[0], &lwork, &iwork[0], &liwork, &info);
		
		lwork = static_cast<int> (work[0]);
		liwork = iwork[0];
		
		work.resize(lwork);
		iwork.resize(liwork);
		
		// compute eigenvalue/eigenvector
		dsyevr_(&V, &I, &U, &mN, mHessian, &mN, NULL, NULL, 
        	&il, &iu, &abs_tol_eigenvalues, &num_eigenvalues_found,
        	eigenvalues, negative_curv_direction, &mN, NULL,
			&work[0], &lwork, &iwork[0], &liwork, &info);
		
		
		// keep only the active part
		#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			double d = negative_curv_direction[i];
			double x = aX[i];
			double l = mLowerBound[i];
			double u = mUpperBound[i];
			
			if (x-l < tol)
				d = max2(d, 0.0);
			else if(u-x < tol)
				d = min2(d, 0.0);
				
			negative_curv_direction[i] = d;
		}
		
		// select the good sign
		if (ddot_(&mN, mGradient, &I1, negative_curv_direction, &I1) > 0.0)
		{
			dscal_(&mN, &minus_one, negative_curv_direction, &I1);
		}
		
		// normalize it
		double scale = 1.0 / dnrm2_(&mN, negative_curv_direction, &I1);
		dscal_(&mN, &scale, negative_curv_direction, &I1);
		
		std::cout << "\tNegative curvature direction computed. Corresponding eigenvalue: " << eigenvalues[0] << std::endl;
	}
	else
	{
		dcopy_(&mN, &D0, &I0, negative_curv_direction, &I1);
	}
	// choose direction
	double *Bd = negative_curv_direction+mN;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, negative_curv_direction, &I1, &D0, Bd, &I1); 
	
	double dBd = ddot_(&mN, Bd, &I1, negative_curv_direction, &I1);
	double d_g = ddot_(&mN, mGradient, &I1, negative_curv_direction, &I1);
	double norm_nd = dnrm2_(&mN, newton_direction, &I1);
	
	std::cout << "\tg_nd = " << g_nd << std::endl;
	
	if (g_nd < tau * norm_nd * (d_g + 0.5*dBd))
	{
		memcpy(mP, newton_direction, mSizeVect);
		std::cout << "\tNewton direction chosen" << std::endl;
	}
	else
	{
		const double norm_ag = dnrm2_(&mN, active_gradient, &I1);
		if (fabs(d_g) < epsilon*norm_ag)
		{
			memcpy(mP, active_gradient, mSizeVect);
			std::cout << "\tGradient direction chosen" << std::endl;
		}
		else
		{
			memcpy(mP, negative_curv_direction, mSizeVect);
			std::cout << "\tNegative direction chosen" << std::endl;
		}
	}
}
#endif

