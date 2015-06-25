
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
		computeSearchDirection(x, &localLowerBound[0], &localUpperBound[0]);
		
		
		// line search
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		
		// try to extend the limits to the boundaries (-> global line search)
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
	
	// compute vector v = y - Bs
	double *v = mWorkSpaceVect;
	memcpy(v, mYk, mSizeVect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, mSk, &I1, &D1, v, &I1); 	
	
	const double vs = ddot_(&mN, v, &I1, mSk, &I1);
	const double threshold_abs_vs = eps1 * dnrm2_(&mN, mSk, &I1) * dnrm2_(&mN, v, &I1);
	
	// update if satisfies conditions
	double *vvT = mWorkSpaceMat;
	bool skip_update = false;
	if (fabs(vs) > threshold_abs_vs)
	{
		// --- Hessian update
		const double inverse_vs = 1.0/vs;
	
		// compute Matrix v.v^T / vs
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			double prefactor = v[i] * inverse_vs;
			//daxpy_(&mN, &prefactor, v, &I1, &mHessian[i*mN], &I1);
			dcopy_(&mN, v, &I1, &vvT[i*mN], &I1);
			dscal_(&mN, &prefactor, &vvT[i*mN], &I1);
		}
		const double norm_update = dnrm2_(&mN_sq, vvT, &I1);
		skip_update = (norm_update > 1e8);
	}
	else
	{
		skip_update = true;
	}
		
	if (!skip_update)
	{
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
void OptSQPSR1::lineSearch(double *aAlpha, double *aX, double *aF)
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
	// TODO:
	// Consider the case when the second Wolfe condition is not satisfied
	
	*aF = evaluateFunctionForLineSearch(aX, a);
	*aAlpha = a;
	daxpy_(&mN, aAlpha, mP, &I1, aX, &I1);
}


// ----------------------------------------------------------------------
double OptSQPSR1::zoom(double aAlo, double aAhi, const double *aX, const double& aPhi0, const double& aPhi0Prime, const double& aAphiLo, const double& c1, const double& c2)
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
#endif // STRONG_WOLFE_LINE_SEARCH_SR1


// ----------------------------------------------------------------------

void OptSQPSR1::computeSearchDirection(const double *aX, const double *aLocalLowerBound, const double *aLocalUpperBound)
{
	const int maximum_iterations_cg = mN*2;
	const char trans = 'N';
	
	const double epsilon = 0.1;
	
	double Delta_sq;
	const double Delta_sq_min = 0.5;
	const double theta = 1e-8;
	if (mStep == 0)
	{
		Delta_sq = 5.0;
	}
	else
	{
		Delta_sq = max2(Delta_sq_min, 10.0*ddot_(&mN, mSk, &I1, mSk, &I1));
	} 
	
	double *p = mWorkSpaceVect;
	double *r = mWorkSpaceMat;
	double *w = r + mN;
	double *next_iterate = w+mN;
	double *projected_gradient = next_iterate+mN;
	
	// compute the projected gradient direction
	memcpy(projected_gradient, mGradient, mSizeVect);
	dscal_(&mN, &minus_one, projected_gradient, &I1);
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		double pg = projected_gradient[i];
		pg = min2(pg, aLocalUpperBound[i]);
		pg = max2(pg, aLocalLowerBound[i]);
		projected_gradient[i] = pg;
	}
	
	const double threshold_residual = square(epsilon)*ddot_(&mN, projected_gradient, &I1, projected_gradient, &I1);
	
	// start with a null starting point
	dcopy_(&mN, &D0, &I0, mP, &I1);
	
	// residual
	memcpy(r, projected_gradient, mSizeVect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, mP, &I1, &D1, r, &I1);
	
	double rho_prev = 1.0;
	double rho = ddot_(&mN, r, &I1, r, &I1);
	
	bool cg_converged = false;
	int step_cg = 0;
	while (step_cg++ < maximum_iterations_cg && !cg_converged)
	{
		std::cout << "\trho = " << rho << std::endl;
		// test stopping criteria
		if (rho < threshold_residual)
		{
			cg_converged = true;
			std::cout << "\tCG convergence reached" << std::endl;
		}
		else
		{
			// compute conjugate gradient direction
			if (step_cg == 1)
			{
				memcpy(p, r, mSizeVect);
				//dscal_(&mN, &minus_one, p, &I1);
			} 
			else
			{
				const double beta = rho/rho_prev;
				dscal_(&mN, &beta, p, &I1);
				//daxpy_(&mN, &minus_one, r, &I1, p, &I1);
				daxpy_(&mN, &D1, r, &I1, p, &I1);
			}
			// compute max allowed step
			const double pp = ddot_(&mN, p, &I1, p, &I1);
			const double dp = ddot_(&mN, mP, &I1, p, &I1);
			const double dd = ddot_(&mN, mP, &I1, mP, &I1);
			const double xi = sqrt(square(dp) - pp*dd + pp*Delta_sq);
			double alpha = (xi-dp) / pp;
			const double tol_p = 1e-12;
			for (int i(0); i<mN; ++i)
			{
				const double l = aLocalLowerBound[i] - mP[i];
				const double u = aLocalUpperBound[i] - mP[i];
				const double p_ = p[i];
				if (p_ < -tol_p)
					alpha = min2(alpha, l/p_);
				else if (p_ > tol_p)
					alpha = min2(alpha, u/p_);
			}
			const double alpha_max = alpha;
			// compute curvature informations
			dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, p, &I1, &D0, w, &I1);
			const double gamma = ddot_(&mN, p, &I1, w, &I1);
			// consider the bounds and negative curvatures
			if (gamma > 0)
			{
				alpha = min2(alpha, rho/gamma);
			}
			else if (step_cg > 1)
			{
				cg_converged = true;
				std::cout << "\tCG convergence reached: negative curvature found" << std::endl;
			}
			if (!cg_converged)
			{
				// compute new iterate
				memcpy(next_iterate, mP, mSizeVect);
				daxpy_(&mN, &alpha, p, &I1, next_iterate, &I1);
				const double gd = ddot_(&mN, mGradient, &I1, next_iterate, &I1);
				const double threshold_gd = -theta*dnrm2_(&mN, mGradient, &I1) * dnrm2_(&mN, next_iterate, &I1);
				if (gd > threshold_gd)
				{
					cg_converged = true;
					std::cout << "\tCG convergence reached: next iterate not descent enough" << std::endl;
				}
				else if (alpha == alpha_max)
				{
					memcpy(mP, next_iterate, mSizeVect);
					cg_converged = true;
					std::cout << "\tCG convergence reached: reached border" << std::endl;
				}
				else
				{
					std::cout << "\talpha = " << alpha << ", gamma = " << gamma << std::endl;
					memcpy(mP, next_iterate, mSizeVect);
					alpha = -alpha;
					daxpy_(&mN, &alpha, w, &I1, r, &I1);
					rho_prev = rho;
					rho = ddot_(&mN, r, &I1, r, &I1);
				}
			}
		}
	}
}

