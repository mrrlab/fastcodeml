
#include "OptSQP.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>


// ----------------------------------------------------------------------
//	Class members definition: BOXCQP
// ----------------------------------------------------------------------


void BOXCQP::solveQP(const double *B, const double *d, double *x)
{
	std::vector<int> IPIV(mN);
	int INFO;
	
	// initialize parameters
	memcpy(mLHS, B, mN*sizeof(double));
	
	// null lagrangian multiplicators
	dcopy_(&mN, &D0, &I0, mLambda, &I1);
	dcopy_(&mN, &D0, &I0, mMu, &I1);

	// solution of the unconstrained problem
	memcpy(x, d, mN*sizeof(double));
	dscal_(&mN, &minus_one, x, &I1);
	memcpy(mLHS, B, mN*mN*sizeof(double));
	dgesv(&mN, &I1, mLHS, &mN, &IPIV[0], x, &mN, &INFO);
	
	if(INFO != 0)
		std::cout << "Error: couldn't solve the initial linear system in BOXCQP. INFO: " << INFO << std::endl;
	
	
	// verify if the solution is valide
	bool convergenceReached(true);
	#pragma omp parallel for reduction(&&: convergenceReached)
	for(size_t i(0); i<mN; ++i)
	{
		convergenceReached = convergenceReached && (ma[i] <= x[i] && x[i] <= mb[i] );
	}
	
	// main loop
	for(size_t step(0); !convergenceReached; ++step)
	{
		// --- update the sets
		updateSets(x);
		
		// --- update parameters
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			switch(mSets[i])
			{
				case LSET:
					x[i] 		= ma[i];
					mMu[i] 		= 0.0;
					break;
				
				case USET:
					x[i] 		= mb[i];
					mLambda[i] 	= 0.0;
					break;
				
				case SSET:
					mMu[i] 		= 0.0;
					mLambda[i] 	= 0.0;
					break;
			};
		}
		
		// setup the known vectors; set their values to zero where we don't know it. 
		memcpy(mLambda_known, 	mLambda, 	mN*sizeof(double));
		memcpy(mMu_known, 		mMu, 		mN*sizeof(double));
		memcpy(mx_known, 		x, 			mN*sizeof(double));
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			switch(mSets[i])
			{
				case LSET:
					mLambda_known[i] = 0.0;
					break;
				
				case USET:
					mMu_known[i] = 0.0;
					break;
				
				case SSET:
					mx_known[i] = 0.0;
					break;
			};
		}
		
		// --- solve new linear system
		
		// setup the left hand side matrix
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			switch(mSets[i])
			{
				case LSET:
					for(size_t j(0); j<mN; ++j) {dcopy_(&mN, &D0, &I0, &mLHS[i*mN], &I1);}
					mLHS[i*(mN+1)] = -1.0;
					break;
				case USET:	
					for(size_t j(0); j<mN; ++j) {dcopy_(&mN, &D0, &I0, &mLHS[i*mN], &I1);}
					mLHS[i*(mN+1)] = 1.0;
					break;
				case SSET:
					for(size_t j(0); j<mN; ++j) {dcopy_(&mN, &B[i*mN], &I1, &mLHS[i*mN], &I1);}
					break;
			};
		}
		
		// setup the right hand side vector
		memcpy(mRHS, d, mN*sizeof(double));
		char trans = 'N';
		dgemv_(&trans, &mN, &mN, &D1, B, &mN, mx_known, &I1, &D1, mRHS, &I1);
		daxpy_(&mN, &D1, mMu_known, &I1, mRHS, &I1);
		daxpy_(&mN, &minus_one, mLambda_known, &I1, mRHS, &I1);
		dscal_(&mN, &minus_one, mRHS, &I1);
		
		// solve linear system
		dgesv(&mN, &I1, mLHS, &mN, &IPIV[0], mRHS, &mN, &INFO);
		if(INFO != 0)
		std::cout << "Error: couldn't solve the linear system in BOXCQP. INFO: " << INFO << std::endl;
		
		// update solutions
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			switch(mSets[i])
			{
				case LSET:
					mLambda[i] = mRHS[i];
					break;
				
				case USET:
					mMu[i] = mRHS[i];
					break;
				
				case SSET:
					x[i] = mRHS[i];
					break;
			};
		}
		
		// --- verify validity of the solution
		convergenceReached = true;
		#pragma omp parallel for reduction(&&: convergenceReached)
		for(size_t i(0); i<mN; ++i)
		{
			// local variable
			bool iVariableCorrect;
			switch(mSets[i])
			{
				case LSET:
					iVariableCorrect = mLambda[i] >= 0.0;
					break;
				case USET:
					iVariableCorrect = mMu[i] >= 0.0;
					break;
				case SSET:
					iVariableCorrect = (x[i] >= ma[i]) && (x[i] <= mb[i]);
					break;
			};
			convergenceReached = convergenceReached && iVariableCorrect;
		}	
	}
}


// ----------------------------------------------------------------------
void BOXCQP::alocateMemory(void)
{
	mSets.resize(mN);
	mSpace.resize(mN*mN + 6*mN);
	
	mLambda = &mSpace[0];
	mMu 	= mLambda + mN;
	mx_known 		= mMu + mN;
	mMu_known 		= mx_known + mN;
	mLambda_known 	= mMu_known + mN;
	mRHS	= mLambda_known + mN;
	mLHS	= mRHS + mN;
}


// ----------------------------------------------------------------------
void BOXCQP::updateSets(double *ax)
{
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		if(  (ax[i] < ma[i])   ||   (ax[i] == ma[i] && mLambda[i] >= 0.0) )
			mSets[i] = LSET;
		else if( (ax[i] > mb[i])   ||   (ax[i] == mb[i] && mMu[i] >= 0.0) )
			mSets[i] = USET;
		else
			mSets[i] = SSET;
	}
}
// ----------------------------------------------------------------------






// ----------------------------------------------------------------------
//	Class members definition: OptSQP
// ----------------------------------------------------------------------
double OptSQP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	alocateMemory();
	
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
}


// ----------------------------------------------------------------------
void OptSQP::SQPminimizer(double *f, double *x)
{
	double f_prev;
	*f = evaluateFunction(x, mTrace);
	
	double alpha = 1.0;
	double scale_s;
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	mQPsolver.reset(new BOXCQP(mN, &localLowerBound[0], &localUpperBound[0]));	
	
	// initialize hessian matrix to identity
	int N_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&N_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
	
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
		
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Quadratic program solving..." << std::endl;
		
		// solve quadratic program to get the search direction		
		mQPsolver->solveQP(mHessian, mGradient, mP);
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		
		alpha = 1.0;
		
		// line search
		lineSearch(&alpha, x, f);
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			if(x[i] < mLowerBound[i])
				x[i] = mLowerBound[i];
			if(x[i] > mUpperBound[i])
				x[i] = mUpperBound[i];
		}
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Step length found:" << alpha << std::endl;
		
		memcpy(&mXEvaluator[0], x, size_vect);
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *f);
		}
		
		// update the system
		computeGradient(x, *f, mGradient);
		
		memcpy(mSk, x, size_vect);
		daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
		memcpy(mYk, mGradient, size_vect);
		daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "BFGS update..." << std::endl;
			
		BFGSupdate();
		
		// check convergence
		convergenceReached = fabs(f_prev - *f) < mAbsoluteError;
	}
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, size_vect);
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if(mStopIfBigger && f >= mThreshold) throw nlopt::forced_stop();
	
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
	double f;
	memcpy(&mXEvaluator[0], x, size_vect);
	
	for(size_t i(0); i<mN; ++i)
	{
		eh = sqrt(DBL_EPSILON) * ( 1.0 + fabs(x[i]) );
		//eh = 1e-7 * ( 1.0 + fabs(x[i]) );
		if( x[i] + eh > mUpperBound[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		eh = mXEvaluator[i] - x[i];
		f = -mModel->computeLikelihood(mXEvaluator, false);
		aGrad[i] = (f-f0)/eh;
		mXEvaluator[i] = x[i];
	}
}


// ----------------------------------------------------------------------
void OptSQP::BFGSupdate(void)
{
	// local variables
	double ys, sBs, inverse_sBs, inverse_ys, theta, sigma, theta_tmp;
	double *Bs, *BssB, *yy;
	char trans = 'N';
	
	int mN_sq = mN*mN;
	
#if 0	
	// apply scaling against roundoff errors
	double scale_s = double(mN) / dnrm2_(&mN, mSk, &I1);
	dscal_(&mN, &scale_s, mSk, &I1);
	dscal_(&mN, &scale_s, mYk, &I1);
#endif
	
	// compute vector B*mSk
	Bs = mWorkSpaceVect;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	sBs = ddot_(&mN, mSk, &I1, Bs,  &I1);
	ys  = ddot_(&mN, mSk, &I1, mYk, &I1);
	
	
	
	// Powell-SQP update:
	// change y so the matrix is positive definite
	sigma = 0.2; // empirical value
	
	if(ys < sigma * sBs)
	{
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "BFGS update leading to a non positive-definite matrix, performing Powell SQP:" << std::endl;
		
		theta = sBs - ys;
		
		if(fabs(theta) > 1e-8)
		{
			do
			{
				theta_tmp = (1.0 - sigma) * sBs / theta;
				sigma = 0.9*sigma;
			}while(theta_tmp >= 1.0);
			
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
	double c1 (1e-1);
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	double phi_0 = *f, phi, phi_prev;
	double a_prev = *aalpha;
	double phi_a_prime;
	double a = *aalpha;
	
	if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi_prev = phi_0;
	a_prev = 0.;
	phi = evaluateFunctionForLineSearch(x, a);
	
	double sigma = 0.3;
	
	// begin by a backtrace
	size_t iter = 0;
	while(phi > phi_0 + phi_0_prime*a*c1 && iter < 10)
	{
		++iter;
		a_prev = a;
		phi_prev = phi;
		sigma = 0.3+0.3*randFrom0to1();
		a *= sigma;
		phi = evaluateFunctionForLineSearch(x, a);
	}
	
	// compute the derivative
	double eh = sqrt(DBL_EPSILON);
	if( a+eh >= 1.0 ) {eh = -eh;}
	phi_a_prime = (phi - evaluateFunctionForLineSearch(x, a+eh))/eh;
	
	iter = 0;
	if(phi_a_prime < 0.0 && a != *aalpha)
	{
		double a0 = a_prev;
		while(phi < phi_prev && iter < 10)
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			sigma = 0.3+0.4*randFrom0to1();
			a = a + sigma*(a0-a);
			phi = evaluateFunctionForLineSearch(x, a);
		}
		if(phi_prev < phi)
		{
			phi = phi_prev;
			a = a_prev;
		}
	}
	else
	{
		while(phi < phi_prev && iter < 10)
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			sigma = 0.5+0.5*sqrt(randFrom0to1());
			a *= sigma;
			phi = evaluateFunctionForLineSearch(x, a);
		}
		if(phi_prev < phi)
		{
			phi = phi_prev;
			a = a_prev;
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
	// strong Wolfe conditions constants
	const double c1(1e-2);
	const double c2(0.9);
	
	mNumCallZoom = 0;
	
	// local variables
	double a_prev = 0.;
	double a = *aalpha;
	double amax = 1.0;
	double phi_0 = *f;
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1); 
	double phi_a_prime;
	
	double phi = evaluateFunctionForLineSearch(x, a);
	double phi_prev = phi_0;
	
	if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime : " << phi_0_prime << std::endl;
	
	
	// first perform a backtrace to  avoid too large values
	double back_multiplier = 0.3;
	int back_iter(0);
	for(back_iter = 0; back_iter < 10 && phi > phi_0 + c1*a*phi_0_prime; ++back_iter)
	{
		a *= back_multiplier;
		phi = evaluateFunctionForLineSearch(x, a);
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "LINE_SEARCH_DEBUG: backtracking: " << a << " " << phi << std::endl;
	}
	
	if(back_iter==0)
		a *= back_multiplier;
	amax = a/back_multiplier;
	a_prev = a*back_multiplier;
	
	bool WolfeSatisfied(false);
	int iter = 0;
	while(!WolfeSatisfied && iter<10)
	{
		if( (phi > phi_0 + c1*a*phi_0_prime) || (phi >= phi_prev && iter > 0) )
		{
			*aalpha = zoom(a_prev, a, x, phi_0, phi_0_prime, c1, c2);
			WolfeSatisfied = true;
		}
		else
		{
			// evaluate the derivative phi_a_prime
			double eh = sqrt(DBL_EPSILON);
			if( a+eh >= 1.0 ) {eh = -eh;}
			phi_a_prime = (phi - evaluateFunctionForLineSearch(x, a+eh))/eh;
		
			if(fabs(phi_a_prime) <= -c2*phi_0_prime)
			{
				*aalpha = a;
				WolfeSatisfied = true;
			}
			else if( phi_a_prime >= 0.0 )
			{
				*aalpha = a = zoom(a, a_prev, x, phi_0, phi_0_prime, c1, c2);
				WolfeSatisfied = true;
			}
			else
			{
				double propa = 0.5;
				a = propa*a + (1.-propa)*amax;
				++iter;
			}
		}
		phi = evaluateFunctionForLineSearch(x, a);
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "LINE_SEARCH_INFO: strong wolfe line search: " << a << " " << phi << std::endl;
	}
	
	
	*aalpha = a;
	daxpy_(&mN, aalpha, mP, &I1, x, &I1);	
	*f = phi;
}


// ----------------------------------------------------------------------
double OptSQP::zoom(double low, double high, const double *x, 
					double const& phi_0, double const& phi_0_prime, 
					double const& c1, double const& c2)
{
	// bisection
	double proplow = 0.5;
	double a = proplow*low + (1.0-proplow)*high;
	
	++ mNumCallZoom;
	if(mNumCallZoom > 5)
		return a;
	
	double phi = evaluateFunctionForLineSearch(x, a);
	double phihigh = evaluateFunctionForLineSearch(x, high);
	
	if( phi > phi_0 + c1*a*phi_0_prime || phi >= phihigh )
		return zoom(low, a, x, phi_0, phi_0_prime, c1, c2);
	else
	{
		// evaluate the derivative phi_a_prime to verify the second Wolfe condition
		double eh = sqrt(DBL_EPSILON);
		if( a+eh >= 1.0 ) {eh = -eh;}
		double phi_a_prime = (phi - evaluateFunctionForLineSearch(x, a+eh))/eh;
		
		if( fabs(phi_a_prime) <= -c2*phi_0_prime)
			return a;
		if( phi_a_prime*(high-low) >= 0.0 )
			return zoom(a, low, x, phi_0, phi_0_prime, c1, c2);
		return zoom(a, high, x, phi_0, phi_0_prime, c1, c2);
	}
	
}

#endif // STRONG_WOLFE_LINE_SEARCH

