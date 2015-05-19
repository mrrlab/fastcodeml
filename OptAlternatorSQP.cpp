#include "OptAlternatorSQP.h"

// ----------------------------------------------------------------------
//	Class members definition: OptSQP
// ----------------------------------------------------------------------
double OptAlternatorSQP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	allocateMemory();
	
	double maxl = 1e7;
	AlternatorSQPminimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------
void OptAlternatorSQP::allocateMemory(void)
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
void OptAlternatorSQP::AlternatorSQPminimizer(double *f, double *x)
{
	*f = evaluateFunction(x, mTrace);
	
	int roundsFullSpace		= 3;
	int roundsBranchSpace	= 3;
	
	double alpha = 1.0;
	
	int branchSpaceCounter 	= 0;
	int fullSpaceCounter	= 0;
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	std::vector<double> localLowerBound(mN);
	std::vector<double> localUpperBound(mN);
	mQPsolverFull.reset(new BOXCQP(mN, &localLowerBound[0], &localUpperBound[0]));
	mQPsolverBL.reset(new BOXCQP(mNumTimes, &localLowerBound[0], &localUpperBound[0]));
	
	
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
		double f_prev = *f;
		memcpy(mGradPrev, mGradient, size_vect);
		memcpy(mXPrev, x, size_vect);
		
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Quadratic program solving..." << std::endl;
		
		// solve quadratic program to get the search direction	
		bool QPsolutionOnBorder;	
		switch(mSearchSpace)
		{
			case SPACE_FULL:
				mQPsolverFull->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder);
				
				++fullSpaceCounter;
				if(fullSpaceCounter > roundsFullSpace)
				{
					mSearchSpace = SPACE_BRANCHES_ONLY;
					fullSpaceCounter = 0;
				}
				break;
				
			case SPACE_BRANCHES_ONLY:
				mQPsolverBL->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder);
				// make the search direction only along the branch lengths
				{
				int mNExtra = mN-mNumTimes;
				dcopy_(&mNExtra, &D0, &I0, &mP[mNumTimes], &I1);
				}
				++branchSpaceCounter;
				if(branchSpaceCounter > roundsBranchSpace)
				{
					mSearchSpace = SPACE_FULL;
					branchSpaceCounter = 0;
				}
				break;
			default:
				throw "Wrong mSearchSpace";
		}
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		
		alpha = 1.0;
		
		// line search
		lineSearch(&alpha, x, f);
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
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
				
		memcpy(mSk, x, size_vect);
		daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
		computeGradient(x, *f, mGradient);
		
		memcpy(mYk, mGradient, size_vect);
		daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "BFGS update..." << std::endl;
			
		BFGSupdate();
		
		// check convergence
		convergenceReached =   fabs(f_prev - *f) < mAbsoluteError
							|| mStep >= mMaxIterations;
	}
}


// ----------------------------------------------------------------------
double OptAlternatorSQP::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, size_vect);
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if(mStopIfBigger && f >= mThreshold) throw nlopt::forced_stop();
	
	return -f;
}


// ----------------------------------------------------------------------
double OptAlternatorSQP::evaluateFunctionForLineSearch(const double* x, double alpha)
{
	memcpy(mWorkSpaceVect, x, size_vect);
	daxpy_(&mN, &alpha, mP, &I1, mWorkSpaceVect, &I1);
	return evaluateFunction(mWorkSpaceVect, mTrace);
}


// ----------------------------------------------------------------------
void OptAlternatorSQP::computeGradient(const double *x, double f0, double *aGrad)
{
	volatile double eh;
	double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], x, size_vect);
	size_t i;
	double *delta = &mWorkSpaceVect[0];
	
	// branch lengths
	for(i=0; i<mNumTimes; ++i)
	{
		eh = sqrt_eps * ( 1.0 + x[i] );
		if( x[i] + eh > mUpperBound[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		delta[i] = mXEvaluator[i] - x[i];
	}
	
	for(i=0; i<mNumTimes; ++i)
	{
		f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
		aGrad[i] = (f-f0)/delta[i];
	}
	
	// other variables
	memcpy(&mXEvaluator[0], x, size_vect);
	for(; i<mN; ++i)
	{
		eh = sqrt_eps * ( 1.0 + fabs(x[i]) );
		if( x[i] + eh > mUpperBound[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		eh = mXEvaluator[i] - x[i];
		f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
		aGrad[i] = (f-f0)/eh;
		mXEvaluator[i] = x[i];
	}
}


// ----------------------------------------------------------------------
void OptAlternatorSQP::BFGSupdate(void)
{
	// local variables
	double ys, sBs, theta, sigma, theta_tmp;
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
	for(int i=0; i<mN; ++i)
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
	for(int i=0; i<mN; ++i)
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


// ----------------------------------------------------------------------
void OptAlternatorSQP::lineSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	double c1 (1e-1);
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	double phi_0 = *f, phi, phi_prev;
	double a_prev = 0.0;
	double phi_a_prime;
	double a = *aalpha;
	
	if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi_prev = phi_0;
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
			a = a_prev;
			phi = evaluateFunctionForLineSearch(x, a);
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
			a = a_prev;
			phi = evaluateFunctionForLineSearch(x, a);
		}
	}
	
	
	*f = phi;
	*aalpha = a;
	daxpy_(&mN, aalpha, mP, &I1, x, &I1);
}


// ----------------------------------------------------------------------
#if 0
void OptAlternatorSQP::getSpaceProperties(int& idFirstVar, int& idLastVar) const
{
	switch(mSearchSpace)
	{
		case SPACE_FULL:
			idFirstVar	= 0;
			idLastVar	= mN;
			break;
			
		case SPACE_BRANCHES_ONLY:
			idFirstVar	= 0;
			idLastVar	= mNumTimes;
			break;
		
		case SPACE_EXTRA_ONLY:
			idFirstVar	= mNumTimes;
			idLastVar	= mN;
			break;		
	}
}
#endif
