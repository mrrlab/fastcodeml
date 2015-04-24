
#include "OptSQP.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>


// ----------------------------------------------------------------------
//	Class members definition: OptSQP
// ----------------------------------------------------------------------
double OptSQP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	alocateMemory();
	
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
	if(mN > i)
	{
		mLowerBound[i] = 0.0;
		mUpperBound[i] = 1.0;
	}
#endif
	
	double maxl = 1e7;
	SQPminimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------
void OptSQP::alocateMemory(void)
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


#ifdef SCALE_OPT_VARIABLES
// ----------------------------------------------------------------------
void OptSQP::scaleVariables(double *x)
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
void OptSQP::unscaleVariables(double *x)
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
#endif // SCALE_OPT_VARIABLES


// ----------------------------------------------------------------------
void OptSQP::SQPminimizer(double *f, double *x)
{
#ifdef SCALE_OPT_VARIABLES
	scaleVariables(x);
#endif // SCALE_OPT_VARIABLES
	
	double f_prev;
	*f = evaluateFunction(x, mTrace);
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "Solution after Bootstrap:";
		mModel->printVar(mXEvaluator, *f);
	}
	
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
	
	
#ifdef SCALE_OPT_VARIABLES
	// change the space of the hessian approximation representation
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		double slb = mLowerBound[i];
		double sub = mUpperBound[i];
		double lb = mLowerBoundUnscaled[i];
		double ub = mUpperBoundUnscaled[i];
		
		double scale = (ub-lb)/(sub-slb);
		scale = scale*scale;
		
		mHessian[i*(mN+1)] *= scale;
	}		
#endif // SCALE_OPT_VARIABLES

		
	// main loop
	bool convergenceReached = false;
	for(mStep = 0; !convergenceReached; ++mStep)
	{
		// update local bounds
		memcpy(&localLowerBound[0], &mLowerBound[0], size_vect);
		memcpy(&localUpperBound[0], &mUpperBound[0], size_vect);
		daxpy_(&mN, &minus_one, x, &I1, &localLowerBound[0], &I1);
		daxpy_(&mN, &minus_one, x, &I1, &localUpperBound[0], &I1);
		
		// save current parameters
		f_prev = *f;
		memcpy(mGradPrev, mGradient, size_vect);
		memcpy(mXPrev, x, size_vect);
		
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Quadratic program solving..." << std::endl;
		
		// solve quadratic program to get the search direction		
		bool QPsolutionOnBorder;
		mQPsolver->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder);
		
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		
		double alpha = 1.0;
#if 0
		// attempt to get a longer step if possible. 
		if(not QPsolutionOnBorder)
		{
			// find biggest possible step size
			alpha = 1e8;
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
		}
#endif
		
		// line search
		lineSearch(&alpha, x, f);
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			if (x[i] < mLowerBound[i])
				x[i] = mLowerBound[i];
			if (x[i] > mUpperBound[i])
				x[i] = mUpperBound[i];
		}
		
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Step length found:" << alpha << std::endl;
		
				
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *f);
		}
		
#if 0
		// try to select a better solution around the new point
		const size_t max_rand_tries = static_cast<const size_t>(sqrt(mN));
		
		double new_f(*f);
		double move_distance = alpha*dnrm2_(&mN, mP, &I1);
		for(size_t rand_try(0); rand_try<max_rand_tries; ++rand_try)
		{
			// build a random point
			for(size_t i(0); i<mN; ++i)
			{
				double my_rand_x_i = x[i];
				
				if (mActiveSet[i] == 0)
				{
					double tmp = 0.5-randFrom0to1();
					tmp *= square(tmp);
					my_rand_x_i += move_distance * (0.3*tmp);
				
					if (my_rand_x_i < mLowerBound[i])
						my_rand_x_i = mLowerBound[i];
					if (my_rand_x_i > mUpperBound[i])
						my_rand_x_i = mUpperBound[i];
				}
				
				mXEvaluator[i] = my_rand_x_i;
			}
			new_f = -mModel->computeLikelihood(mXEvaluator, mTrace);

			std::cout << "RANDOM TRY: f=" << *f << ", attempt : " << new_f << std::endl;
			double df = new_f - *f;
			//std::cout << "df" << df << ", prob barrier : " << exp(-df/Temperature) << std::endl;
			if (df <= 0.0 || randFrom0to1() <= exp(-df/Temperature))
			{
				memcpy(x, &mXEvaluator[0], size_vect);
				*f = new_f;
			}			
			Temperature *= 0.9;
		}
		*f = evaluateFunction(x, mTrace);
#endif		
		
		
		
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
			//const int max_count_lower = (mN > 30 ? static_cast<const int>(log(static_cast<double>(mN))) : 1);
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
#ifdef SCALE_OPT_VARIABLES
	unscaleVariables(x);
#endif // SCALE_OPT_VARIABLES
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, size_vect);
#ifdef SCALE_OPT_VARIABLES
	unscaleVariables(&mXEvaluator[0]);
#endif // SCALE_OPT_VARIABLES
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if (mStopIfBigger && f >= mThreshold) throw FastCodeMLEarlyStopLRT();
	
	return -f;
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunctionForLineSearch(const double* x, double alpha)
{
	memcpy(mWorkSpaceVect, x, size_vect);
	daxpy_(&mN, &alpha, mP, &I1, mWorkSpaceVect, &I1);
	return evaluateFunction(mWorkSpaceVect, mTrace);
}


// ----------------------------------------------------------------------
void OptSQP::computeGradient(const double *x, double f0, double *aGrad)
{
	volatile double eh;
	double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], x, size_vect);
	size_t i;
	double *delta = mWorkSpaceVect;
	
	
#ifdef SCALE_OPT_VARIABLES
	double *x_ = mWorkSpaceMat;
	memcpy(x_, x, size_vect);
	unscaleVariables(&mXEvaluator[0]);
	unscaleVariables(x_);
#else
	const double *x_ = x;
#endif // SCALE_OPT_VARIABLES
	
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
#ifdef SCALE_OPT_VARIABLES
	#pragma omp parallel for
	for (size_t j(0); j<mN; ++j)
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
void OptSQP::BFGSupdate(void)
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
	
#ifdef SELECTIVE_SIZING_STRATEGY
	double ss = ddot_(&mN, mSk, &I1, mSk, &I1);
	double YS_SS = ys/ss;
	double SBS_SS = sBs/ss;
	
	double thetaSS, gammaSS;
	
	double eps1 = 1e-3;
	double eps2 = 7e-1;
	double tau1 = 0.5;
	double tau2 = 1e6;
	
	bool apply_sizing_factor = false;
	
	// compute the sizing factor
	if (mStep == 0)
	{
		gammaSS = ys/sBs;
		gammaSS = max2(gammaSS, eps2);
		apply_sizing_factor = true;
	}
	else
	{
		thetaSS = min2(tau1, tau2*ss);
		gammaSS  = ( (1.0-thetaSS)*mYS_SS_prev  + (thetaSS)*YS_SS );
		gammaSS /= ( (1.0-thetaSS)*mSBS_SS_prev + (thetaSS)*SBS_SS);
		if (gammaSS < 1.0-eps1)
		{
			gammaSS = max2(eps2, gammaSS);
			apply_sizing_factor = true;
		}
	}
	
	// apply the sizing factor if needed
	if (apply_sizing_factor)
	{
		dscal_(&mN_sq, &gammaSS, mHessian, &I1);
	
		// update the BFGS variables
		sBs *= gammaSS;
		dscal_(&mN, &gammaSS, Bs, &I1);
	}
	
	// save previous variables for next update
	mYS_SS_prev = YS_SS;
	mSBS_SS_prev = SBS_SS;	
#endif // SELECTIVE_SIZING_STRATEGY	
	
	
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


#ifndef STRONG_WOLFE_LINE_SEARCH
// ----------------------------------------------------------------------
void OptSQP::lineSearch(double *aalpha, double *x, double *f)
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

#else
// ----------------------------------------------------------------------
void OptSQP::lineSearch(double *aalpha, double *x, double *f)
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
double OptSQP::zoom(double alo, double ahi, double *x, const double& phi_0, const double& phi_0_prime, const double& aphi_lo, const double& c1, const double& c2)
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
#endif // STRONG_WOLFE_LINE_SEARCH

