
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
	mSpace.resize(2*mN*mN + mN*9);
	
	
	mWorkSpaceVect = &mSpace[0];
	mWorkSpaceMat = mWorkSpaceVect + mN;
	
	mGradient = mWorkSpaceMat + mN*mN;
	mProjectedGradient = mGradient + mN;
	mP = mProjectedGradient + mN;
	mHessian = mP + mN;
	
	mSk =  mHessian + mN*mN;
	mYk = mSk + mN;
	
	mXPrev = mYk + mN;
	mGradPrev = mXPrev + mN;
	
	// set active set counters to zero
	mActiveSet.resize(mN, 0);
	mFixedVariables.resize(mN, 0);
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
	*f = evaluateFunction(x, mTrace);
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "Solution after Bootstrap:";
		mModel->printVar(mXEvaluator, *f);
	}
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	// compute the projected gradient
	memcpy(mProjectedGradient, mGradient, mSizeVect);
	dscal_(&mN, &minus_one, mProjectedGradient, &I1);
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		double p = mProjectedGradient[i];
		p = max2(p, mLowerBound[i]-x[i]);
		p = min2(p, mUpperBound[i]-x[i]);
		mProjectedGradient[i] = p;
	}
	
	// bounds for the QP subproblem (l-x0 <= p <= u-x0)
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	
	// initialize hessian matrix to identity
	int N_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&N_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
	
	// ------------------------------------------- main loop
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
		
		// find the fixed variables
		updateFixedVariables(x);
		
		// decide if we use an SPG iteration or a CG iteration
		double *gi = mWorkSpaceVect;
		memcpy(gi, mProjectedGradient, mSizeVect);
		#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			if (mFixedVariables[i] == 1)
			{
				gi[i] = 0.0;
			}
		}
		const double norm_gi = dnrm2_(&mN, gi, &I1);
		const double norm_pg = dnrm2_(&mN, mProjectedGradient, &I1);
		const double mu = 0.8;
		
		if (norm_gi < mu*norm_pg)
		{
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "SPG iteration..." << std::endl;
			// perform a spg iteration
			spectralProjectedGradientIteration(x, f);
		}
		else
		{
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Approximate quadratic program solving ..." << std::endl;
				
			// solve approximately the quadratic program to get the search direction		
			computeSearchDirection(x, &localLowerBound[0], &localUpperBound[0]);
			// line search
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "Line Search..." << std::endl;
		
			lineSearch(x, f);
		}
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "New solution:";
			mModel->printVar(mXEvaluator, *f);
		}
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			x[i] = min2(mUpperBound[i], max2(mLowerBound[i], x[i]));
		}
		
		// check convergence
#if 1
		const int index_amax_pg = idamax_(&mN, mProjectedGradient, &I1);
		const double norm_projected_gradient = fabs(mProjectedGradient[index_amax_pg]);
		//std::cout << "|projected gradient| = " << norm_projected_gradient << std::endl;
		
		convergenceReached =   (fabs(f_prev - *f) < mAbsoluteError) && (norm_projected_gradient < mAbsoluteError)
							|| mStep >= mMaxIterations;
#else		
		convergenceReached =   fabs(f_prev - *f) < mAbsoluteError
							|| mStep >= mMaxIterations;
#endif
		
		if (!convergenceReached)
		{
			// update the system
			
			computeGradient(x, *f, mGradient);
			
			memcpy(mProjectedGradient, mGradient, mSizeVect);
			dscal_(&mN, &minus_one, mProjectedGradient, &I1);
			#pragma omp parallel for
			for (int i(0); i<mN; ++i)
			{
				double p = mProjectedGradient[i];
				p = max2(p, localLowerBound[i]);
				p = min2(p, localUpperBound[i]);
				mProjectedGradient[i] = p;
			}
				
			memcpy(mSk, x, mSizeVect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
			memcpy(mYk, mGradient, mSizeVect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
		
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "SR1 update..." << std::endl;
			
		
			// update the Hessian information
			SR1update();
		
			// update the active set
			activeSetUpdate(x, 1e-4);
		}
	}
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
		const double norm_hessian = dnrm2_(&mN_sq, mHessian, &I1);
		skip_update = (norm_update > 1e5*norm_hessian);
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

// ----------------------------------------------------------------------

void OptSQPSR1::activeSetUpdate(double *aX, const double aTolerance)
{
	// number of iterations where we skip the gradient computation (component i)
	const int max_count_lower = static_cast<const int>(1.3*log (static_cast<double>(mN)/10.)) + 1;//(mN>30 ? 1:0);
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
void OptSQPSR1::updateFixedVariables(const double *aX)
{
	const double tolerance_fixed_variable = 1e-6;
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		const double x = aX[i];
		if (  (x-mLowerBound[i] < tolerance_fixed_variable)
		   || (mUpperBound[i]-x < tolerance_fixed_variable))
		{
			mFixedVariables[i] = 1;
		}
		else
		{
			mFixedVariables[i] = 0;
		}
	}
}


// ----------------------------------------------------------------------
void OptSQPSR1::spectralProjectedGradientIteration(double *aX, double *aF)
{
	// spectral gradient coefficient
	double spg_coeff;
	const double spg_coeff_min = 1e-4;
	const double spg_coeff_max = 1e3;
	
	const double ys = ddot_(&mN, mSk, &I1, mYk, &I1);
	if ( (mStep == 0) || (ys < 0.0) )
	{
		spg_coeff = 1.0;
	}
	else
	{
		spg_coeff = ys / ddot_(&mN, mSk, &I1, mSk, &I1);
		spg_coeff = min2(spg_coeff_max, max2(spg_coeff_min, spg_coeff));
	}
	
	// search direction
	spg_coeff = -spg_coeff;
	memcpy(mP, mGradient, mSizeVect);
	dscal_(&mN, &spg_coeff, mP, &I1);
	#pragma omp parallel for
	for (int i(0); i<mN; ++i)
	{
		const double l = mLowerBound[i]-aX[i];
		const double u = mUpperBound[i]-aX[i];
		double d = mP[i];
		d = min2(u, max2(l, d));
		mP[i] = d;
	}
	
	// backtracking line search
	const double phi0 = *aF;
	const double phi0_prime = ddot_(&mN, mGradient, &I1, mP, &I1);
	double alpha = 1.0;
	backtrackingLineSearch(aX, aF, alpha, phi0, phi0_prime);
}



// ----------------------------------------------------------------------
void OptSQPSR1::computeSearchDirection(const double *aX, const double *aLocalLowerBound, const double *aLocalUpperBound)
{
	char trans = 'N';
	
	// set of fixed variables to avoid long loops of O(N)
	std::vector<int> J; J.reserve(mN);
	for (int i(0); i<mN; ++i)
	{
		if (mFixedVariables[i] == 1)
		{
			J.push_back(i);
		}
	}
	
	// "trust region" radius
	double Delta_sq;
	const double Delta_sq_min = sqrt(mN)*5e-1;
	if (mStep == 0)
	{
		Delta_sq = 5.0;
	}
	else
	{
		Delta_sq = max2(Delta_sq_min, 10.0*ddot_(&mN, mSk, &I1, mSk, &I1));
	}
	
	const int maximum_iterations_cg = mN+10;
	
	const double epsilon = 0.1;
	const double theta = 1e-8;
		
	double *p = mWorkSpaceVect;
	double *r = mWorkSpaceMat;
	double *w = r + mN;
	double *next_iterate = w+mN;
	double *gi = next_iterate+mN;	// restricted projected gradient
	
	memcpy(gi, mProjectedGradient, mSizeVect);
	for (int i(0); i<J.size(); ++i) {gi[J[i]] = 0.0;}
	
	const double threshold_residual = square(epsilon)*ddot_(&mN, gi, &I1, gi, &I1);
	
	// start at point x
	dcopy_(&mN, &D0, &I0, mP, &I1);
	
	// residual
	memcpy(r, gi, mSizeVect);
	dgemv_(&trans, &mN, &mN, &minus_one, mHessian, &mN, mP, &I1, &D1, r, &I1);
	for (int i(0); i<J.size(); ++i)	{r[J[i]] = 0.0;}
	
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
				for (int i(0); i<J.size(); ++i)	{p[J[i]] = 0.0;}
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
				if (mFixedVariables[i] == 0)
				{
					const double l = aLocalLowerBound[i] - mP[i];
					const double u = aLocalUpperBound[i] - mP[i];
					const double p_ = p[i];
					if (p_ < -tol_p)
						alpha = min2(alpha, l/p_);
					else if (p_ > tol_p)
						alpha = min2(alpha, u/p_);
				}
			}
			const double alpha_max = alpha;
			// compute curvature informations
			dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, p, &I1, &D0, w, &I1);
			for (int i(0); i<J.size(); ++i)	{w[J[i]] = 0.0;}
			const double gamma = ddot_(&mN, p, &I1, w, &I1);
			// consider the bounds and negative curvatures
			if (gamma > 0)
			{
				alpha = min2(alpha, rho/gamma);
			}
			else if (step_cg > 1)
			{
				memcpy(mP, mProjectedGradient, mSizeVect);
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


// ----------------------------------------------------------------------
void OptSQPSR1::lineSearch(double *aX, double *aF)
{
	// compute data required for line search
	const double phi0 = *aF;
	const double phi0_prime = ddot_(&mN, mGradient, &I1, mP, &I1);
	
	std::cout << "\tphi0_prime = " << phi0_prime << std::endl;
	
	double alpha_max = 1e16;
	const double tol_division = 1e-12;
	for (int i(0); i<mN; ++i)
	{
		if (mFixedVariables[i] != 1)
		{
			const double l = mLowerBound[i] - aX[i];
			const double u = mUpperBound[i] - aX[i];
			double p = mP[i];
			if (p < -tol_division)
				alpha_max = min2(alpha_max, l/p);
			else if (p > tol_division)
				alpha_max = min2(alpha_max, u/p);
		}
	}
	double alpha = min2(1.0, alpha_max);
	double phi; 
	
	if (alpha_max > 1.0) 	// x+mP is inside the domain
	{
		phi = evaluateFunctionForLineSearch(aX, 1.0);
		if (phi <= (phi0 + mGamma*phi0_prime))
		{
			// compute the line derivative at point alpha
			double eh = sqrt(DBL_EPSILON);
			double phi_a_prime = (evaluateFunctionForLineSearch(aX, 1.0+eh) - phi)/eh;
			std::cout << "\tphia_prime = " << phi_a_prime << std::endl;
			if (phi_a_prime >= mBeta*phi0_prime)
			{
				alpha = 1.0;
				daxpy_(&mN, &alpha, mP, &I1, aX, &I1);
				*aF = phi;
			}
			else
			{
				std::cout << "\tExtrapolation..." << std::endl;
				// extrapolation
				extrapolatingLineSearch(aX, aF, alpha, alpha_max);
			}
		}
		else
		{	
			std::cout << "\tBacktrace..." << std::endl;
			// backtracking
			backtrackingLineSearch(aX, aF, alpha, phi0, phi0_prime);
		}
	}
	else				// x+mP is outside the domain
	{
		phi = evaluateFunctionForLineSearch(aX, alpha_max);
		std::cout << "\tphi_max = " << phi << std::endl;
		if (phi < phi0)
		{
			std::cout << "\tExtrapolation..." << std::endl;
			// extrapolation
			extrapolatingLineSearch(aX, aF, alpha, alpha_max);
		}
		else
		{
			std::cout << "\tBacktrace..." << std::endl;
			// backtracking
			backtrackingLineSearch(aX, aF, alpha, phi0, phi0_prime);
		}
	}
}


// ----------------------------------------------------------------------
void OptSQPSR1::backtrackingLineSearch(double *aX, double *aF, const double aAlpha, const double aPhi0, const double aPhi0_prime)
{
	double alpha = aAlpha;
	const double s1 = 0.3, s2 = 0.9;
	
	double phi = evaluateFunctionForLineSearch(aX, alpha);
	int counter = 0;
	const int max_iter_backtrace = 15;
	while ( (counter++ < max_iter_backtrace) && (phi > aPhi0 + alpha * mGamma * aPhi0_prime) )
	{
		std::cout << "\t\talpha = " << alpha << std::endl;
		alpha *= (s1 + randFrom0to1()*(s2-s1));
		phi = evaluateFunctionForLineSearch(aX, alpha);
	}
	*aF = phi;
	daxpy_(&mN, &alpha, mP, &I1, aX, &I1);
}


// ----------------------------------------------------------------------
void OptSQPSR1::extrapolatingLineSearch(double *aX, double *aF, const double aAlpha, const double aAlphaMax)
{
	const double N = 2.0;
	double alpha = aAlpha;
	double alpha_trial;
	
	bool converged = false;
	while (!converged)
	{
		if (alpha < aAlphaMax && N*alpha > aAlphaMax)
			alpha_trial = aAlphaMax;
		else
			alpha_trial = N*alpha;
	
		double *p_a = mWorkSpaceMat;
		double *p_a_trial = p_a + mN;
	
		memcpy(p_a, aX, mSizeVect);
		memcpy(p_a_trial, aX, mSizeVect);
		daxpy_(&mN, &alpha, mP, &I1, p_a, &I1);
		daxpy_(&mN, &alpha_trial, mP, &I1, p_a_trial, &I1);
		#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			const double l = mLowerBound[i];
			const double u = mUpperBound[i];
			double p_a_ = p_a[i];
			double p_a_trial_ = p_a_trial[i];
			p_a_ = max2(min2(p_a_, u), l);
			p_a_trial_ = max2(min2(p_a_trial_, u), l);
			p_a[i] = p_a_;
			p_a_trial[i] = p_a_trial_;
		}
		double *diff = p_a_trial+mN;
		memcpy(diff, p_a_trial, mSizeVect);
		daxpy_(&mN, &minus_one, p_a, &I1, diff, &I1);
		int index_max = idamax_(&mN, diff, &I1);
		const double inf_norm_diff = fabs(diff[index_max]);
	
		if (alpha >= aAlphaMax && inf_norm_diff < mAbsoluteError)
		{
			memcpy(aX, p_a, mSizeVect);
			*aF = evaluateFunction(aX, mTrace); 
			converged = true;
			std::cout << "\t\tExtrapolating: Converged after reaching border" << std::endl;
		}
		else
		{
			double phi_a = evaluateFunction(p_a, mTrace); 
			double phi_a_trial = evaluateFunction(p_a_trial, mTrace);
			if (phi_a_trial >= phi_a)
			{
				memcpy(aX, p_a, mSizeVect);
				*aF = phi_a;
				converged = true;
				std::cout << "\t\tExtrapolating: Converged after increasing step" << std::endl;
			}
			else
			{
				alpha = alpha_trial;
			}
		}
	}
}





