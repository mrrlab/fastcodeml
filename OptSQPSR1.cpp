
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
	mHessian = mP + mN;
	
	mSk =  mHessian + mN*mN;
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
			std::cout << "Computing new search direction..." << std::endl;
		
		// solve approximately the quadratic program to get the search direction		
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
			x[i] = min2(mUpperBound[i], max2(mLowerBound[i], x[i]));
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
			
		
			// update the Hessian information
			SR1update();
		
			// update the active set
			// TODO
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
	
	// only update if well defined
	if (fabs(vs) > threshold_abs_vs)
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

void OptSQPSR1::computeSearchDirection(const double *aX)
{

	char trans = 'N';
	
	// --- compute the generalized Cauchy point (GCP)
	
	double *x = mWorkSpaceVect;
	double *g = mWorkSpaceMat;
	double *d = g + mN;
	double *Bd = d + mN;
	double *b = Bd + mN;
	
	std::vector<int> active_set_cauchy(mN, 0);
	std::vector<int> J; J.reserve(mN);
	
	// x vector
	memcpy(x, aX, mSizeVect);
	// g = gradient - Bx
	memcpy(g, mGradient, mSizeVect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, x, &I1, &D1, g, &I1);
	// projected gradient in the limit of a small gradient
	const double tol_active_set = 1e-4;
	memcpy(d, mGradient, mSizeVect);
	dscal_(&mN, &minus_one, d, &I1);
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		double d_ = d[i];
		const double x_ = x[i];
		const double l_ = mLowerBound[i];
		const double u_ = mUpperBound[i];
		d_ = (x_-l_ < tol_active_set) ? max2(d_, 0.0) : d_;
		d_ = (u_-x_ < tol_active_set) ? min2(d_, 0.0) : d_;
		d[i] = d_;
	}
	// Bd
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, d, &I1, &D0, Bd, &I1);
	double f_prime = ddot_(&mN, mGradient, &I1, d, &I1);
	double f_double_prime = ddot_(&mN, d, &I1, Bd, &I1);
	
	if (f_prime < 0.0)
	{
		bool GCP_found = false;
		while(!GCP_found)
		{
			// find the next break point
			double delta_t = 1e16;
			const double tol_d = 1e-6;
			for (int i(0); i<mN; ++i)
			{
				const double l_ = mLowerBound[i]-x[i];
				const double u_ = mUpperBound[i]-x[i];
				const double d_ = d[i];
				if (d_ < -tol_d)
					delta_t = min2(delta_t, l_/d_);			
				else if (d_ > tol_d)
					delta_t = min2(delta_t, u_/d_);
			}
			// find the set J
			for (int i(0); i<J.size(); ++i){active_set_cauchy[J[i]] = 1;}
			J.clear();
			for (int i(0); i<mN; ++i)
			{
				const double x_next = x[i] + delta_t*d[i];
				if (   fabs(x_next - mLowerBound[i]) < tol_active_set
					|| fabs(x_next - mUpperBound[i]) < tol_active_set )
				{
					if (active_set_cauchy[i] != 1)
					{
						J.push_back(i);
					}
				}
			}
			// verify if the GCP is in this interval
			double neg_ratio_f = -f_prime/f_double_prime;
			if (f_double_prime > 0.0 && 0.0 < neg_ratio_f && neg_ratio_f < delta_t)
			{
				daxpy_(&mN, &neg_ratio_f, d, &I1, x, &I1);
				GCP_found = true;
			}
			else
			{
				// update line derivatives
				dcopy_(&mN, &D0, &I0, b, &I1);
				for (int i(0); i<J.size(); ++i)
				{
					int line = J[i];
					daxpy_(&mN, &d[line], &mHessian[line*mN], &I1, b, &I1);
				}
				daxpy_(&mN, &delta_t, d, &I1, x, &I1);
				f_prime += delta_t*f_double_prime;
				f_prime -= ddot_(&mN, b, &I1, x, &I1);
			
				f_double_prime -= 2.0 * ddot_(&mN, b, &I1, d, &I1);
				for (int i(0); i<J.size(); ++i)
				{
					int line = J[i];
					f_prime -= d[line]*g[line];
					f_double_prime += d[line]*b[line];
					d[line] = 0.0;
				}
				if (f_prime >= 0.0)
					GCP_found = true;
			}
		}
	}
	for (int i(0); i<J.size(); ++i){active_set_cauchy[J[i]] = 1;}
	J.clear();
	for (int i(0); i<mN; ++i){if (active_set_cauchy[i] == 1) J.push_back(i);} // J is now the active set
	
	std::cout << "\tCauchy point calculated, refining the solution with a conjugate gradient method..." << std::endl;
	
	// --- refine the solution with a conjugate gradient algorithm
	
	double *dx = mWorkSpaceMat;
	double *r = dx + mN;
	double *p = r + mN;
	double *y = p + mN;
	
	// compute residual
	memcpy(dx, x, mSizeVect);
	daxpy_(&mN, &minus_one, aX, &I1, dx, &I1);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, dx, &I1, &D1, r, &I1);
	daxpy_(&mN, &minus_one, mGradient, &I1, r, &I1);
	for (int i(0); i<J.size(); ++i) {r[J[i]] = 0.0;}
	
	// set p = 0
	dcopy_(&mN, &D0, &I0, p, &I1);
	
	const double eta_sq = 1e-2;
	
	double rho_1 = 1.0;
	double rho_2 = ddot_(&mN, r, &I1, r, &I1);
	
	bool CG_converged = false;
	while (!CG_converged)
	{
		std::cout << "\t CG: rho_2 = " << rho_2 << std::endl;
		if (rho_2 < eta_sq)
		{
			break;
		}
		
		double beta = rho_2 / rho_1;
		dscal_(&mN, &beta, p, &I1);
		daxpy_(&mN, &D1, r, &I1, p, &I1);
	
		// compute y = Bp (restricted to the active set
		dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, p, &I1, &D1, y, &I1);
		for (int i(0); i<J.size(); ++i) {y[J[i]] = 0.0;}
	
		// compute alpha1, largest positive real s.t. l < x+alpha1*p < u
		double alpha_1 = 1e16;
		const double tol_p = 1e-4;
		for (int i(0); i<mN; ++i)
		{
			if (active_set_cauchy[i] == 1)
				continue;
			
			double l_ = mLowerBound[i] - x[i];
			double u_ = mUpperBound[i] - x[i];
			double p_ = p[i];
			if (p_ < -tol_p)
				alpha_1 = min2(alpha_1, l_/p_);			
			else if (p_ > tol_p)
				alpha_1 = min2(alpha_1, u_/p_);
		}
		
		const double py = ddot_(&mN, p, &I1, y, &I1);
		if (py < 0.0)
		{
			daxpy_(&mN, &alpha_1, p, &I1, x, &I1);
			CG_converged = true;
		}
		else
		{
			double alpha_2 = rho_2 / py;
			if (alpha_2 > alpha_1)
			{
				daxpy_(&mN, &alpha_1, p, &I1, x, &I1);
				CG_converged = true;
			}
			else
			{
				daxpy_(&mN, &alpha_2, p, &I1, x, &I1);
				alpha_2 = -alpha_2;
				daxpy_(&mN, &alpha_2, y, &I1, r, &I1);
				
				rho_1 = rho_2;
				rho_2 = ddot_(&mN, r, &I1, r, &I1);
			}
		}
	}
	
	memcpy(mP, x, mSizeVect);
	daxpy_(&mN, &minus_one, aX, &I1, x, &I1);
}

