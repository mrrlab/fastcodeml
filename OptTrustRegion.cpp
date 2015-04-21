
#include "OptTrustRegion.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>


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
		
		std::cout << slb << " " << " " << sub << " " << lb << " " << ub << std::endl;
		mHessian[i*diag_stride] *= scale;
		std::cout << i << " : " << mHessian[i*diag_stride] << std::endl;
	}		
#endif // SCALE_OPT_TRUST_REGION_VARIABLES

	// trust region algorithm parameters
	double eta1, _eta1_, _eta2_;
	double gamma1, gamma2, gamma3;
	eta1 = 0.25;
	_eta1_ = 0.25;
	_eta2_ = 0.75;
	gamma1 = 2.5e-1;
	gamma2 = 7.5e-1;
	gamma3 = 2.0e+0;
	double trust_region_radius = 0.05;
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
		bool trust_region_step_found = false;
		int trust_region_tries = 0;
		int trust_region_max_tries = 10;
		while(!trust_region_step_found)
		{
			++ trust_region_tries;
			// --- choose the trust region radius
			if (mStep > 0)
			{
				double retrospective_ratio = improvement_ratio; // computeRatio(*f, mP, f_prev);
				std::cout << "retrospective_ratio : " << retrospective_ratio << std::endl;
				if (retrospective_ratio < _eta1_)
				{
					// take trust radius between gamma1 and gamma2 * radius (randomly)
					trust_region_radius = (gamma1 + randFrom0to1()*(gamma2-gamma1))*trust_region_radius;
				}
				else if (retrospective_ratio > _eta2_)
				{
					// take trust radius between 1 and gamma3 * radius (randomly)
					trust_region_radius = (1.0 + randFrom0to1()*(gamma3-1.0))*trust_region_radius;
				}
				else
				{
					// take trust radius between gamma2 and 1 * radius (randomly)
					trust_region_radius = (gamma2 + randFrom0to1()*(1.0-gamma2))*trust_region_radius;
				}
			}
	
	
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
			bool QPsolutionOnBorder;
			double *x_candidate = mWorkSpaceVect;
			
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Quadratic program solving..." << std::endl;
			mQPsolver->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder);
		
			// candidate solution
			memcpy(x_candidate, x, size_vect);
			daxpy_(&mN, &D1, mP, &I1, x_candidate, &I1);
			double f_candidate = evaluateFunction(x_candidate, mTrace);
		
			// --- update the solution accordingly to the expected improvement
			improvement_ratio = computeRatio(*f, mP, f_candidate);
			std::cout << "improvement_ratio : " << improvement_ratio << std::endl;
			std::cout << "trust_region_radius : " << trust_region_radius << std::endl;
			
		
			if (improvement_ratio >= eta1)
			{
				memcpy(x, x_candidate, size_vect);
				*f = f_candidate;
				trust_region_step_found = true;
			}
			else
			{
				if (trust_region_tries > trust_region_max_tries)
				{
					memcpy(x, x_candidate, size_vect);
					*f = f_candidate;
					break;
				}
			}
		}
		
				
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
				std::cout << "BFGS update..." << std::endl;
			
			// update the B matrix
			BFGSupdate();
		
			// update the active set
			const int max_count_lower = static_cast<const int>(1.3*log (static_cast<double>(mN)/10.)) + (mN>30 ? 1:0);
			const int max_count_upper = (mN > 30 ? 1 : 0);
	 
			#pragma omp parallel for
			for(size_t i(0); i<mN; ++i)
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
	double dm = ddot_(&mN, dx, &I1, mGradient, &I1);
	
	// compute vector B*dx
	double *Bdx = mWorkSpaceMat;
	char trans = 'N';
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, dx, &I1, &D0, Bdx, &I1); 
		
	
	dm += ddot_(&mN, Bdx, &I1, dx,  &I1);
	if (fabs(dm) < 1e-12) {dm = 1e-12;}
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
		if(mActiveSet[i] == 0)
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
void OptTrustRegion::BFGSupdate(void)
{
	// local variables
	double ys, sBs, inverse_sBs, inverse_ys, theta, sigma, theta_tmp;
	double *Bs, *BssB, *yy;
	char trans = 'N';
	
	int mN_sq = mN*mN;
	
	// compute vector B*mSk
	Bs = mWorkSpaceVect;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	sBs = ddot_(&mN, mSk, &I1, Bs,  &I1);
	ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
	

	// Powell-SQP update:
	// change y so the matrix is positive definite
	sigma = 0.2; // empirical value
	
	if (ys < sigma * sBs)
	{
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "BFGS update leading to a non positive-definite matrix, performing Powell SQP:" << std::endl;
		
		theta = sBs - ys;
		
		if (fabs(theta) > 1e-8)
		{
			do
			{
				theta_tmp = (1.0 - sigma) * sBs / theta;
				sigma = 0.9*sigma;
			} while(theta_tmp >= 1.0);
			
			theta = theta_tmp;
			theta_tmp = 1.0 - theta;
			
			dscal_(&mN, &theta, mYk, &I1);
			daxpy_(&mN, &theta_tmp, Bs, &I1, mYk, &I1);
			ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
		}
	}
	
	// compute Matrix B*mSk * mSk^T*B
	BssB = mWorkSpaceMat;
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		double prefactor = - Bs[i] / sBs;
		dcopy_(&mN, &Bs[0], &I1, &BssB[i*mN], &I1);
		dscal_(&mN, &prefactor, &BssB[i*mN], &I1);
	}
	
	// add the BssB / sBs contribution
	daxpy_(&mN_sq, &D1, BssB, &I1, mHessian, &I1);
	
	// compute matrix y**T * y
	yy = mWorkSpaceMat;
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		double prefactor = mYk[i] / ys;
		dcopy_(&mN, &mYk[0], &I1, &yy[i*mN], &I1);
		dscal_(&mN, &prefactor, &yy[i*mN], &I1);
	}
	
	// add the yy / ys contribution
	daxpy_(&mN_sq, &D1, yy, &I1, mHessian, &I1);
	
	// make the diagonal more important in order to avoid non positive definite matrix, 
	// due to roundoff errors
	int diag_stride = mN+1;
	double factor = 1.1;
	double inv_factor = 1.0/factor;
	dscal_(&mN_sq, &inv_factor, mHessian, &I1);
	dscal_(&mN, &factor, mHessian, &diag_stride);
}


