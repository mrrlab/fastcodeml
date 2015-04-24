
#include "OptTrustRegion.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>
#include <iomanip>

// ----------------------------------------------------------------------
//	Class members definition: OptTrustRegion
// ----------------------------------------------------------------------
double OptTrustRegion::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	alocateMemory();
	
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	// set the scaling
	int i;
	// get a better scaling for the branch lengths variables.
	// it is often near ~0.25, multiply it by 4 so it is "more" around 1
	// in the new space representation
	mUpperBound.assign(mNumTimes, 500.0);
	
	
	i = mNumTimes + 1; 		// v1
	mUpperBound[i] = 20.0;
	
	i = mNumTimes + 2; 		// w0
	mUpperBound[i] = 20.0;
	
	
	// shrink the w2 variable between 0 and 1 so it is about the same scale as 
	// the other variables in the new space representation
	i = mNumTimes+4; 		// w2
	if(mN > i)
	{
		mLowerBound[i] = 0.0;
		mUpperBound[i] = 1.0;
	}
#endif // SCALE_OPT_TRUST_REGION_VARIABLES
	
	double maxl = 1e7;
	SQPminimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------
void OptTrustRegion::alocateMemory(void)
{
	size_vect = mN*sizeof(double);
	
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


#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
// ----------------------------------------------------------------------
void OptTrustRegion::scaleVariables(double *x)
{
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
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
void OptTrustRegion::unscaleVariables(double *x)
{
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
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
#endif // SCALE_OPT_TRUST_REGION_VARIABLES


// ----------------------------------------------------------------------
void OptTrustRegion::SQPminimizer(double *f, double *x)
{
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	scaleVariables(x);
#endif // SCALE_OPT_TRUST_REGION_VARIABLES
	
	double f_prev;
	*f = evaluateFunction(x, mTrace);
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	// bounds for the QP subproblem
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	mQPsolver.reset(new BOXCQP(mN, &localLowerBound[0], &localUpperBound[0]));	
	
	// initialize hessian matrix to identity
	int N_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&N_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
	
	
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	// change the space of the hessian approximation representation
	//#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		
		double scale = (ub-lb)/(sub-slb);
		scale = scale*scale;		
		mHessian[i*diag_stride] *= scale;
	}		
#endif // SCALE_OPT_TRUST_REGION_VARIABLES

	// trust region algorithm parameters
	const double threshold_acceptance_ratio = 0.2;
	const double max_trust_region_radius = 1.2;
	double trust_region_radius = 0.1;
	double improvement_ratio = 1.0;
	
	// main loop
	bool convergenceReached = false;
	for(mStep = 0; !convergenceReached; ++mStep)
	{
	
		// save current parameters
		f_prev = *f;
		memcpy(mGradPrev, mGradient, size_vect);
		memcpy(mXPrev, x, size_vect);
		
		
		// compute the next step using a trust region strategy	
		int rejected_trust_region_radius = 0;
		const int max_rejected_trust_region_radius = 20;
		bool found_improvement = false;
		while (rejected_trust_region_radius < max_rejected_trust_region_radius
			   && not found_improvement)
		{
			++rejected_trust_region_radius;
			
			// --- compute the step in trust region
			// update local bounds
			memcpy(&localLowerBound[0], &mLowerBound[0], size_vect);
			memcpy(&localUpperBound[0], &mUpperBound[0], size_vect);
			daxpy_(&mN, &minus_one, x, &I1, &localLowerBound[0], &I1);
			daxpy_(&mN, &minus_one, x, &I1, &localUpperBound[0], &I1);
			#pragma omp parallel for
			for(size_t i(0); i<mN; ++i)
			{
				double& l = localLowerBound[i];
				double& u = localUpperBound[i];
				l = max2(l, -trust_region_radius);
				u = min2(u,  trust_region_radius);
			}
		
			// solve quadratic program in l-infinity norm
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "conjugate gradient Steihaug algorithm..." << std::endl;
			CG_Steihaug(&localLowerBound[0], &localUpperBound[0], mP);
		
			// candidate solution
			double *x_candidate = mWorkSpaceVect;
			memcpy(x_candidate, x, size_vect);
			daxpy_(&mN, &D1, mP, &I1, x_candidate, &I1);
			double f_candidate = evaluateFunction(x_candidate, mTrace);
		
			// --- choose the trust region radius
			improvement_ratio = computeRatio(f_prev, mP, f_candidate);
			std::cout << "improvement_ratio : " << std::setprecision(4) << std::scientific << improvement_ratio << std::endl;
			
			// compute the l-infinity norm of the previous step
			double length_prev_step = 0.0;
			for(size_t i(0); i<mN; ++i)
			{
				if (fabs(mP[i]) > length_prev_step)
					length_prev_step = fabs(mP[i]);
			}
			std::cout << "trust_region_radius : " << std::setprecision(4) << std::scientific << trust_region_radius << ", |mP| : " << std::setprecision(4) << std::scientific << length_prev_step << std::endl;
			
			if (improvement_ratio < 0.25)
			{
				// reduce trust region radius
				trust_region_radius = 0.25*length_prev_step;
			}
			else if (improvement_ratio > 0.75 && fabs(trust_region_radius - length_prev_step) <= 1e-8)
			{
				// increase trust region radius
				trust_region_radius = min2(2.0*trust_region_radius, max_trust_region_radius);
				std::cout << trust_region_radius << " <= " << max_trust_region_radius << std::endl;
			}
			// otherwise we don't change the trust region radius
		
			// update the solution
			if (improvement_ratio > threshold_acceptance_ratio)
			{
				memcpy(x, x_candidate, size_vect);
				*f = f_candidate;
				found_improvement = true;
			}
		}			
		
#if 0		
		// line search
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line search:";
		
		
		double alpha = 1e8;
		// find biggest possible step size
		double low, high, pi, candidate;
			
		for(size_t i(0); i<mN; ++i)
		{
			low = mLowerBound[i] - x[i];
			high = mUpperBound[i] - x[i];
			pi = mP[i];
				
			if(fabs(pi) > 1e-8)
			{
				candidate = pi > 0 ? high/pi : low/pi;
				alpha = candidate < alpha ? candidate : alpha;
			}
		}
		
		
		lineSearch(&alpha, x, f);
#endif	
				
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *f);
		}	
		
		
		// check convergence
		convergenceReached =   fabs(f_prev - *f) < mAbsoluteError
							|| mStep >= mMaxIterations;
		
		if (not convergenceReached)
		{
			// update the system
			computeGradient(x, *f, mGradient);
				
			memcpy(mSk, x, size_vect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
			memcpy(mYk, mGradient, size_vect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
		
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "hessian update..." << std::endl;
			
			// update the B matrix
			hessianUpdate();
		
			// update the active set
			const int max_count_lower = static_cast<const int>(1.3*log (static_cast<double>(mN)/10.)) + (mN>30 ? 1:0);
			const int max_count_upper = (mN > 30 ? 1 : 0);
	 
			#pragma omp parallel for
			for(size_t i(0); i<mN; ++i)
			{
				if (mActiveSet[i] > 1)
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
						mActiveSet[i] = 1+max_count_lower;
						if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
						std::cout << "Variable " << i << " in the (lower) active set.\n";
					}
					else if (x[i] >= mUpperBound[i] - active_set_tol && mGradient[i] <= 0.0)
					{
						mActiveSet[i] = 1+max_count_upper;
						if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
						std::cout << "Variable " << i << " in the (upper) active set.\n";
					}
					else
					{
						mActiveSet[i] = 0;
					}
				}
			}
		}
	}
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	unscaleVariables(x);
#endif // SCALE_OPT_TRUST_REGION_VARIABLES
}


// ----------------------------------------------------------------------
double OptTrustRegion::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, size_vect);
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	unscaleVariables(&mXEvaluator[0]);
#endif // SCALE_OPT_TRUST_REGION_VARIABLES
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if (mStopIfBigger && f >= mThreshold) throw FastCodeMLEarlyStopLRT();
	
	return -f;
}


// ----------------------------------------------------------------------
double OptTrustRegion::evaluateFunctionForLineSearch(const double* x, double alpha)
{
	memcpy(mWorkSpaceVect, x, size_vect);
	daxpy_(&mN, &alpha, mP, &I1, mWorkSpaceVect, &I1);
	return evaluateFunction(mWorkSpaceVect, mTrace);
}


// ----------------------------------------------------------------------
double OptTrustRegion::computeRatio(double f0, double *dx, double f)
{	
	double *dx_proj = mWorkSpaceVect;
	memcpy(dx_proj, dx, size_vect);
	
	#pragma omp parallel for
	for (size_t i(0); i<mN; ++i)
	{
		if (mActiveSet[i] > 0)
		{
			dx_proj[i] = 0.0;
		}
	}
		
	// compute vector B*dx
	double *Bdx = mWorkSpaceMat;
	char trans = 'N';
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, dx_proj, &I1, &D0, Bdx, &I1); 
	
	
	#pragma omp parallel for
	for (size_t i(0); i<mN; ++i)
	{
		if (mActiveSet[i] > 0)
		{
			Bdx[i] = 0.0;
		}
	}
		
	double dm = ddot_(&mN, dx, &I1, mGradient, &I1);
	dm += ddot_(&mN, Bdx, &I1, dx,  &I1);
	if (fabs(dm) < 1e-16) {dm = 1e-16;}
	return (f-f0)/dm;
}


// ----------------------------------------------------------------------
void OptTrustRegion::computeGradient(const double *x, double f0, double *aGrad)
{
	volatile double eh;
	double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], x, size_vect);
	size_t i;
	double *delta = mWorkSpaceVect;
	
	
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	double *x_ = mWorkSpaceMat;
	memcpy(x_, x, size_vect);
	unscaleVariables(&mXEvaluator[0]);
	unscaleVariables(x_);
#else
	const double *x_ = x;
#endif // SCALE_OPT_TRUST_REGION_VARIABLES
	
	// branch lengths
	for(i=0; i<mNumTimes; ++i)
	{
		eh = sqrt_eps * ( 1.0 + x_[i] );
		if( x_[i] + eh > mUpperBoundUnscaled[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		delta[i] = mXEvaluator[i] - x_[i];
	}

	for(i=0; i<mNumTimes; ++i)
	{
		if(mActiveSet[i] < 2)
		{
			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-f0)/delta[i];
		}
		// otherwise we don't change it
	}
	
	// other variables
	memcpy(&mXEvaluator[0], x_, size_vect);
	for(; i<mN; ++i)
	{
		if (mActiveSet[i] < 2)
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
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
	#pragma omp parallel for
	for (size_t j(0); j<mN; ++j)
	{
		double slb = mLowerBound[j];
		double sub = mUpperBound[j];
		double lb = mLowerBoundUnscaled[j];
		double ub = mUpperBoundUnscaled[j];
		
		aGrad[j] *= (ub-lb)/(sub-slb);
	}
#endif // SCALE_OPT_TRUST_REGION_VARIABLES
}


// ----------------------------------------------------------------------
void OptTrustRegion::hessianUpdate(void)
{
	// local variables
	double *v, *Bs;
	double vs, inverse_vs;
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
	inverse_vs = 1.0/vs;
	
	// compute Matrix v.v^T / vs
	vvT = mWorkSpaceMat;
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		double prefactor = v[i] * inverse_vs;
		dcopy_(&mN, v, &I1, &vvT[i*mN], &I1);
		dscal_(&mN, &prefactor, &vvT[i*mN], &I1);
	}
	
	// add the v.v^T / vs contribution
	daxpy_(&mN_sq, &D1, vvT, &I1, mHessian, &I1);
}


/// find_t_CG_S
/// find t > 0 such that |p+td| = Delta
///
/// @param[in]		localLowerBound The lower bound for p+td
/// @param[in]		localUpperBound	The Upper bound for p+td
/// @param[in]		p				vector p
/// @param[in]		d				vector d
/// @param[in,out]	t				in: high value. out: the solution
///
void find_t_CG_S(const int N, const std::vector<int>& ActiveSet, const double *localLowerBound, const double *localUpperBound, const double *p, const double *d, double &t)
{
	for (size_t i(0); i<N; ++i)
	{
		if (ActiveSet[i] < 2)
		{
			double low  = localLowerBound[i] - p[i];
			double high = localUpperBound[i] - p[i];
			double di = d[i];
			if (fabs(di) < 1e-6)
			{
				continue;
			}
			if (di > 0.0)
			{
				if (t > high/di)
					t = high/di;
				else if (t < low/di)
					t = low/di;
			}
			else
			{
				if (t < high/di)
					t = high/di;
				else if (t > low/di)
					t = low/di;
			}
		}
	}
}


// ----------------------------------------------------------------------
void OptTrustRegion::CG_Steihaug(const double *localLowerBound, const double *localUpperBound, double *search_direction)
{
	const double r_tol = 1e-8;
	char trans = 'N';
	
	double *p = search_direction;	// search direction
	double *r = mWorkSpaceMat;		// residual
	double *d = r + mN;				// improvement direction 
	double *proj_grad = d+mN;
	
	// set p0 = 0
	dcopy_(&mN, &D0, &I0, p, &I1);
	// set r0 = gradient
	memcpy(proj_grad, mGradient, size_vect);
	#pragma omp parallel for
	for (size_t i(0); i<mN; ++i)
	{
		if (mActiveSet[i] > 0)
		{
			proj_grad[i] = 0.0;
		}
	}
	memcpy(r, proj_grad, size_vect);
	// set d0 = -gradient
	memcpy(d, mGradient, size_vect);
	dscal_(&mN, &minus_one, d, &I1);
	
	
	double r0_norm = dnrm2_(&mN, r, &I1);
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "\tr0_norm = " << std::setprecision(4) << std::scientific << r0_norm << std::endl;
	
	if (r0_norm < r_tol)
	{
		return;
	}
	bool search_direction_found = false;
	while (not search_direction_found)
	{
		// compute dT B d
		double *Bd = mWorkSpaceVect;
		dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, d, &I1, &D0, Bd, &I1);
		#pragma omp parallel for
		for (size_t i(0); i<mN; ++i)
		{
			if (mActiveSet[i] < 1)
			{
				Bd[i] = 0.0;
			}
		}
		double dBd = ddot_(&mN, d, &I1, Bd, &I1);
		
		if (dBd < 0.0)
		{
			// find t such that it minimizes m(p_)
			// for p_ = p+td
			// There is only two possibilities: t<0 and t>0, we compare both
			double t_positive = 1e16;
			find_t_CG_S(mN, mActiveSet, localLowerBound, localUpperBound, p, d, t_positive);
			double M_positive = t_positive * ( t_positive*dBd
											 + ddot_(&mN, d, &I1, proj_grad, &I1)
											 + ddot_(&mN, d, &I1, Bd, &I1));			
			
			dscal_(&mN, &minus_one, d, &I1);
			double t_negative = 1e16;
			find_t_CG_S(mN, mActiveSet, localLowerBound, localUpperBound, p, d, t_negative);
			
			double M_negative = t_negative * ( t_negative*dBd
											 - ddot_(&mN, d, &I1, proj_grad, &I1)
											 - ddot_(&mN, d, &I1, Bd, &I1));
			
			if (M_negative < M_positive)
			{
				daxpy_(&mN, &t_negative, d, &I1, p, &I1);
			}
			else
			{
				t_positive = -t_positive; // d has been set to -d
				daxpy_(&mN, &t_positive, d, &I1, p, &I1);
			}			
			return;
		}
		
		// compute a new direction
		double rr = ddot_(&mN, r, &I1, r, &I1);
		double alpha = rr / dBd;
		double *p_next = mWorkSpaceVect;
		memcpy(p_next, p, size_vect);
		daxpy_(&mN, &alpha, d, &I1, p_next, &I1);
		
		// verify if it satisfies the constraints
		bool constraints_satisfied = true;
		for (size_t i(0); i<mN; ++i)
		{
			double low  = localLowerBound[i];
			double high = localUpperBound[i];
			if (p_next[i] < low || p_next[i] > high)
			{
				constraints_satisfied = false;
				break;
			}
		}
		if (not constraints_satisfied)
		{
			// find t such that |p+td| satisfies the bounds
			double t = alpha;
			find_t_CG_S(mN, mActiveSet, localLowerBound, localUpperBound, p, d, t);
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "\tconstraints not satisfied, t found : " << std::setprecision(4) << std::scientific << t << std::endl;
			// update the solution
			daxpy_(&mN, &t, d, &I1, p, &I1);
			return;
		}
		// otherwise keep the new update
		memcpy(p, p_next, size_vect);
		// update r
		daxpy_(&mN, &alpha, Bd, &I1, r, &I1);
		double r_norm = dnrm2_(&mN, r, &I1);
		if (r_norm < r_tol*r0_norm)
		{
			return;
		}
		// update d
		double beta = ddot_(&mN, r, &I1, r, &I1)/rr;
		dscal_(&mN, &beta, d, &I1);
		daxpy_(&mN, &D1, r, &I1, d, &I1);
	}
}



// ----------------------------------------------------------------------
void OptTrustRegion::lineSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	double c1 (2e-1);
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	double phi_0 = *f, phi, phi_prev;
	double a_prev = *aalpha;
	double phi_a_prime;
	double a = *aalpha;
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi_prev = phi_0;
	a_prev = 0.;
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
	
	maxIterBack = maxIterUp = static_cast<int> (ceil( 3.*log(mN+10) ));
	sigma_bas 	= pow(1e-3, 1./static_cast<double>(maxIterBack));
	
	
	// begin by a backtrace
	size_t iter = 0;
	while(phi > phi_0 + phi_0_prime*a*c1 && iter < maxIterBack)
	{
		++iter;
		a_prev = a;
		phi_prev = phi;
		//sigma = 0.3+0.3*randFrom0to1();
		sigma = sigma_bas * (0.9 + 0.2*randFrom0to1());
		a *= sigma;
		phi = evaluateFunctionForLineSearch(x, a);
	}
	
	// compute the derivative
	double eh = sqrt(DBL_EPSILON);
	if ( a+eh >= 1.0 ) {eh = -eh;}
	phi_a_prime = (phi - evaluateFunctionForLineSearch(x, a+eh))/eh;
	
	iter = 0;
	if (phi_a_prime < 0.0 && a != *aalpha)
	{
		double a0 = a_prev;
		while(phi < phi_prev && iter < maxIterBack)
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
		while(phi < phi_prev && iter < maxIterUp)
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

