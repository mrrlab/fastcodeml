
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
	mUpperBound.assign(mNumTimes, 100.0);
	
	
	i = mNumTimes + 1; 		// v1
	mUpperBound[i] = 10.0;
	
	i = mNumTimes + 2; 		// w0
	mUpperBound[i] = 10.0;
	
	
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
	mSpace.resize(2*mN*mN + mN*8);
	
	
	mWorkSpaceVect = &mSpace[0];
	mWorkSpaceMat = mWorkSpaceVect + mN;
	
	mGradient = mWorkSpaceMat + mN*mN;
	mProjectedGradient = mGradient + mN;
	mP = mProjectedGradient + mN;
	mHessian = mP + mN;
	
	mSk = mHessian + mN*mN;
	mYk = mSk + mN;
	
	mXPrev = mYk + mN;
	mGradPrev = mXPrev + mN;
	
	// set active set counters to zero
	mActiveSet.resize(mN, 0);
	mLocalActiveSet.resize(mN, false);
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
	f_prev = *f;
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	// bounds for the QP subproblem
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	
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
	const double gamma1 = 0.3;
	const double gamma2 = 1.8;
	double trust_region_radius = 0.5;
	const double max_trust_region_radius = 5.0;
	double improvement_ratio = 1.0;
	std::vector<double> x_candidate_(mN);
	
	const double activeSetTolerance = 1e-5;
	
	// main loop
	bool convergenceReached = false;
	for(mStep = 0; mStep < mMaxIterations; ++mStep)
	{
	
		// ---- check convergence
		
		// compute the projected gradient
		memcpy(mProjectedGradient, mGradient, size_vect);
		#pragma omp parallel for 
		for (size_t i(0); i<mN; ++i)
		{
			if (mActiveSet[i] != 0)
			{
				mProjectedGradient[i] = 0.0;
			}
		}
		double projected_gradient_norm = dnrm2_(&mN, mProjectedGradient, &I1);
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Convergence check: |g_projected| = " << projected_gradient_norm << std::endl;
			
		
		if (projected_gradient_norm < mAbsoluteError)
			break;
				
		
		// ---- Compute next step
		
		// update local bounds
		dcopy_(&mN, &mLowerBound[0], &I1, &localLowerBound[0], &I1);
		daxpy_(&mN, &minus_one, x, &I1, &localLowerBound[0], &I1);
		dcopy_(&mN, &mUpperBound[0], &I1, &localUpperBound[0], &I1);
		daxpy_(&mN, &minus_one, x, &I1, &localUpperBound[0], &I1);
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			double& l = localLowerBound[i];
			double& u = localUpperBound[i];
			l = max2(l, -trust_region_radius);
			u = min2(u,  trust_region_radius);
		}
		
		// compute the general Cauchy point
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Compute the General Cauchy Point..." << std::endl;
		generalCauchyPoint(&localLowerBound[0], &localUpperBound[0], mP);
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "General cauchy Point computed. " << std::endl;
				
		// compute the local active set for the GCP
		#pragma omp parallel for
		for (size_t i(0); i<mN; ++i)
		{
			double x_GCP_i = x[i] + mP[i];
			if (   fabs(x_GCP_i - localLowerBound[i]) < activeSetTolerance
				|| fabs(x_GCP_i - localUpperBound[i]) < activeSetTolerance )
			{
				mLocalActiveSet[i] = true;
			}
			else
			{
				mLocalActiveSet[i] = false;
			}
		}
		
		// refine the new iterate using a modified conjugate gradient algorithm
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Refine the General Cauchy Point using modified conjugate gradient..." << std::endl;
		modifiedConjugateGradient(&localLowerBound[0], &localUpperBound[0]);
		
		
		// candidate solution
		double *x_candidate = &x_candidate_[0];
		memcpy(x_candidate, x, size_vect);
		daxpy_(&mN, &D1, mP, &I1, x_candidate, &I1);
		#pragma omp parallel for
		for (size_t i(0); i<mN; ++i)
		{
			if (x_candidate[i] < mLowerBound[i])
				x_candidate[i] = mLowerBound[i];
			else if (x_candidate[i] > mUpperBound[i])
				x_candidate[i] = mUpperBound[i];
		}
		double f_candidate = evaluateFunction(x_candidate, mTrace);
		
		// compute the infinity-norm of mP
		double mP_infinity_norm = 0.0;
		int highest_direction = 0;
		for (size_t i(0); i<mN; ++i)
		{
			double abs_mPi = fabs(mP[i]);
			if (abs_mPi > mP_infinity_norm)
			{
				mP_infinity_norm = abs_mPi;
				highest_direction = i;
			}
		}
		
		std::cout << "The limiting variable is the variable " << highest_direction << "/" << mN << std::endl;
		
		// compute the improvement ratio
		double improvement_ratio = computeRatio(*f, mP, f_candidate);
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Improvement ratio of the candidate function : " << improvement_ratio << std::endl;
		
		// change the trust region radius accordingly
		if (improvement_ratio <= threshold_acceptance_ratio)
		{
			trust_region_radius *= gamma1;
		}
		else if (improvement_ratio > 0.95 && fabs(mP_infinity_norm-trust_region_radius) < 1e-3*trust_region_radius)
		{
			trust_region_radius = min2(gamma2*trust_region_radius, max_trust_region_radius);
		}
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Trust region radius : " << trust_region_radius << ", |mP| = " << mP_infinity_norm << std::endl;
		
		// update the solution if needed
		if (improvement_ratio > threshold_acceptance_ratio)
		{
			// save previous solution
			f_prev = *f;
			memcpy(mGradPrev, mGradient, size_vect);
			memcpy(mXPrev, x, size_vect);
			

			// update solution
			memcpy(x, x_candidate, size_vect);
			*f = evaluateFunction(x, mTrace);
			computeGradient(x, *f, mGradient);
			
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			{
				std::cout << "New Solution:";
				mModel->printVar(mXEvaluator, *f);
			}	
				
			// update the hessian matrix
			memcpy(mSk, x, size_vect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
			memcpy(mYk, mGradient, size_vect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "hessian update..." << std::endl;
			hessianUpdate();
			
		}
		
		// update the active set
		#pragma omp parallel for
		for (size_t i(0); i<mN; ++i)
		{
			double x_i = x[i];
			if (   fabs(x_i - mLowerBound[i]) < activeSetTolerance
				|| fabs(x_i - mUpperBound[i]) < activeSetTolerance )
			{
				mActiveSet[i] = 1;
			}
			else
			{
				mActiveSet[i] = 0;
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
	dm += 0.5*ddot_(&mN, Bdx, &I1, dx,  &I1);
	double df = f-f0;
	if (fabs(dm) < 1e-16) {dm = 1e-16;}
	
	return (df < 0.0)? fabs(df/dm) : -fabs(df/dm);
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
		if( x_[i] + 1e8*eh > mUpperBoundUnscaled[i] )
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
#ifdef TRUST_REGION_SR1_MATRIX_UPDATE
	// local variables
	double *v, *Bs;
	double vs, inverse_vs;
	double *vvT;
	char trans = 'N';
	
	int mN_sq = mN*mN;
	
	// compute vector v = y-B*mSk
	v = mWorkSpaceVect;
	memcpy(v, mYk, size_vect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, mSk, &I1, &D1, v, &I1); 	
	
	
	vs = ddot_(&mN, v, &I1, mSk, &I1);
	
	std::cout << "v.s = " << vs << std::endl;
	inverse_vs = 1.0/vs;
	
	// compute Matrix v.v^T / vs
	vvT = mWorkSpaceMat;
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		double prefactor = v[i] * inverse_vs;
		//std::cout << "i=" << i << ", prefactor=" << prefactor << ", v[i]="<< v[i] << ", y[i] = " << mYk[i] << ", s[i]=" << mSk[i] << std::endl;
		dcopy_(&mN, v, &I1, &vvT[i*mN], &I1);
		dscal_(&mN, &prefactor, &vvT[i*mN], &I1);
	}
	
	
	// add the v.v^T / vs contribution only if the modification is not too large
	// Frobenius norm:
	double modification_norm = dnrm2_(&mN_sq, vvT, &I1)/static_cast<double>(mN);
	double previous_hessian_norm = dnrm2_(&mN_sq, mHessian, &I1)/static_cast<double>(mN);
	
	std::cout << "\tModif norm: " << std::scientific << std::setprecision(12) << modification_norm << std::endl;
	std::cout << "\tprev norm: " << std::scientific << std::setprecision(12) << previous_hessian_norm << std::endl;
	
	const double epsilon_ = 1e-8;
	const double gamma_	= 1e8;
	
	if (  (fabs(vs) < epsilon_) || (modification_norm > gamma_*(previous_hessian_norm+1.0))  )
		std::cout << "------------ Matrix update skipped! ---------------" << std::endl;
	else
		daxpy_(&mN_sq, &D1, vvT, &I1, mHessian, &I1);
	
#endif // TRUST_REGION_SR1_MATRIX_UPDATE
}	


// ----------------------------------------------------------------------
/// find_t_GCP
/// find higher t > 0 such that (l <= p+td <= u)
///
/// @param[in]		N				Size of the vectors
/// @param[in]		localLowerBound The lower bound for p+td
/// @param[in]		localUpperBound	The Upper bound for p+td
/// @param[in]		p				vector p
/// @param[in]		d				vector d
/// @param[in,out]	t				in: high value. out: the solution
///
void find_t_GCP(const int N, const double *localLowerBound, const double *localUpperBound, const double *p, const double *d, double &t)
{
	for (size_t i(0); i<N; ++i)
	{
		double di = d[i];
		if (fabs(di) > 1e-16)
		{
			double low  = localLowerBound[i] - p[i];
			double high = localUpperBound[i] - p[i];
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
/// updateJSet
/// Find the variables that "block" the progression of the general cauchy point calculation, i.e. 
/// variables x=p+td : x==l or x==u (componentwise)
///
/// @param[in]		N				Size of the vectors
/// @param[in]		localLowerBound The lower bound for p+td
/// @param[in]		localUpperBound	The Upper bound for p+td
/// @param[in]		p				vector p
/// @param[in]		d				vector d
/// @param[in]		t				parameter t
/// @param[out]		aSetJ			The set J containing the "blocking variables"
///
void updateJSet(const int N, const double *localLowerBound, const double *localUpperBound, const double *p, const double *d, const double t, std::vector<int>& aSetJ, int step_GCP)
{
	const double tolerance = 1e-8;
	#pragma omp parallel for
	for (size_t i(0); i<N; ++i)
	{
		double xi = p[i] + t*d[i];
		if ( aSetJ[i] == 0 && ((fabs(xi-localLowerBound[i]) < tolerance) || (fabs(xi-localUpperBound[i]) < tolerance))  )
		{
			aSetJ[i] = step_GCP;
				//std::cout << "\t\tVariable " << i << " in the set J (step " << step_GCP << ")." << std::endl;
		}
	}
}


// ----------------------------------------------------------------------
void OptTrustRegion::generalCauchyPoint(const double *localLowerBound, const double *localUpperBound, double *cauchy_point_from_x)
{
	char trans = 'N';
	
	// Initialization
	double *d = mWorkSpaceVect;
	double *g = mWorkSpaceMat;
	double *Bd = g+mN;
	double *active_part_d = Bd+mN;
	double f_prime, f_double_prime;
	
	// set J
	// contains the variables indices which block the progression of the alg., i.e. which reach the local bounds.
	// setJ[i] is 0 if the variable is not blocking, n if it blocks at step n
	std::vector<int> setJ(mN, 0);	
	int step_GCP = 1;
	
	dcopy_(&mN, &D0, &I0, cauchy_point_from_x, &I1);
	
	memcpy(d, mProjectedGradient, size_vect);
	//memcpy(d, mGradient, size_vect);
	//double scale_d = -sqrt(DBL_EPSILON)/dnrm2_(&mN, mGradient, &I1);
	dscal_(&mN, &minus_one, d, &I1);	
	memcpy(g, mGradient, size_vect);
	
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, d, &I1, &D0, Bd, &I1);
	f_prime = ddot_(&mN, g, &I1, d, &I1);
	f_double_prime = ddot_(&mN, Bd, &I1, d, &I1);
	
	std::cout << "\tDEBUG : f_prime = " << f_prime << ", f_double_prime = " << f_double_prime << std::endl;
	
	
	if (f_prime < 0)
	{
		bool GCP_found = false;
		while (not GCP_found)
		{
			// find the next breakpoint
			double delta_t = 1e16;
			find_t_GCP(mN, localLowerBound, localUpperBound, cauchy_point_from_x, d, delta_t);
			
			//if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				//std::cout << "\tdelta_t found :" << delta_t << std::endl;
				
			// update J set
			//if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				//std::cout << "\tUpdating set J:" << std::endl;
			updateJSet(mN, localLowerBound, localUpperBound, cauchy_point_from_x, d, delta_t, setJ, step_GCP);
			
			// test if the GCP has been found
			const double ratio_f = -f_prime/f_double_prime;
			if ( (f_double_prime > 0.0) && (0.0 < ratio_f && ratio_f < delta_t) )
			{
				daxpy_(&mN, &ratio_f, d, &I1, cauchy_point_from_x, &I1);
				GCP_found = true;
			}
			else
			{
				// update line derivatives
				memcpy(active_part_d, d, size_vect);
				#pragma omp parallel for
				for (size_t i(0); i<mN; ++i)
				{
					if (setJ[i] != step_GCP)
					{
						active_part_d[i] = 0.0;
					}
				}
				dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, active_part_d, &I1, &D0, Bd, &I1);
				daxpy_(&mN, &delta_t, d, &I1, cauchy_point_from_x, &I1);
				f_prime += delta_t*f_double_prime - ddot_(&mN, Bd, &I1, cauchy_point_from_x, &I1) - ddot_(&mN, active_part_d, &I1, g, &I1);
				f_double_prime += ddot_(&mN, Bd, &I1, active_part_d, &I1) - 2.0*ddot_(&mN, Bd, &I1, d, &I1);
				
				//#pragma omp parallel for
				//std::cout << "d = ";
				for (size_t i(0); i<mN; ++i)
				{
					if (setJ[i] == step_GCP)
					{
						d[i] = 0.0;
					}
					//std::cout << d[i] << " ";
				}
				//std::cout << std::endl;
				
				if (f_prime >= 0.0 || step_GCP >= 2*mN)
					GCP_found = true;
				
				//std::cout << "\tf_prime : " << std::setprecision(15) << f_prime << std::endl; 
				++ step_GCP;
			}
		}
	}
}


// ----------------------------------------------------------------------
void OptTrustRegion::modifiedConjugateGradient(const double *localLowerBound, const double *localUpperBound)
{
	char trans = 'N';

	// initialization
	double *r = mWorkSpaceVect;
	double *p = mWorkSpaceMat;
	double *y = p+mN;
	
	double rho_1, rho_2;
	
	double eta_sq = dnrm2_(&mN, mProjectedGradient, &I1);
	eta_sq *= min2(0.1, sqrt(eta_sq));
	eta_sq *= eta_sq;
	
	// set r <- -B * cauchy_from_x - g
#if 0
	memcpy(r, mGradient, size_vect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, mP, &I1, &minus_one, r, &I1);
#else
	memcpy(y, mP, size_vect);
	#pragma omp parallel for
	for (size_t i(0); i<mN; ++i)
	{
		if (mLocalActiveSet[i])
		{
			y[i] = 0.0;
		}
	}
	memcpy(r, mProjectedGradient, size_vect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, y, &I1, &minus_one, r, &I1);
#endif
	dcopy_(&mN, &D0, &I0, p, &I1);
	
	// restrict the computaion to the free variables, i.e. not in the local active set
	#pragma omp parallel for
	for (size_t i(0); i<mN; ++i)
	{
		if (mLocalActiveSet[i])
		{
			r[i] = 0.0;
		}
	}
	
	rho_1 = 1.0;
	rho_2 = ddot_(&mN, r, &I1, r, &I1);
	
	bool CG_terminated = (rho_2 < eta_sq);
	while (not CG_terminated)
	{
		std::cout << "\t\t rho_2 = " << rho_2 << "/" << eta_sq << std::endl;
		double beta = rho_2 / rho_1;
		// set p <- r + beta * p
		dscal_(&mN, &beta, p, &I1);
		daxpy_(&mN, &D1, r, &I1, p, &I1);
		// set y <- Bp 
		dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, p, &I1, &D0, y, &I1);
		// compute largest a1 s.t. l<=x+a1*p<=u for the FREE variables
		double a1 = 1e16;
		const double tol__ = 1e-8;
		for (size_t i(0); i<mN; ++i)
		{
			if (not mLocalActiveSet[i])
			{
				double pi = p[i];
				double maxa;
				if (pi < -tol__)
				{
					maxa = (localLowerBound[i]-mP[i]) / pi;
				}
				else if (pi > tol__)
				{
					maxa = (localUpperBound[i]-mP[i]) / pi;
				}
				
				if (a1 > maxa)
				{
					a1 = maxa;
				}
			}
		}
		
		double py = ddot_(&mN, p, &I1, y, &I1);
		if (py <= 0.0)
		{
			daxpy_(&mN, &a1, p, &I1, mP, &I1);
			CG_terminated = true;
		}
		else
		{
			double a2 = rho_2/py;
			if (a2 > a1)
			{
				daxpy_(&mN, &a1, p, &I1, mP, &I1);
				CG_terminated = true;
			}
			else
			{
				daxpy_(&mN, &a2, p, &I1, mP, &I1);
				a2 = -a2;
				daxpy_(&mN, &a2, y, &I1, r, &I1);
				rho_1 = rho_2;
				rho_2 = ddot_(&mN, r, &I1, r, &I1);
				CG_terminated = (rho_2 < eta_sq);
			}
		}		
	}
}


