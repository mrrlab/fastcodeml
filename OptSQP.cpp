
#include "OptSQP.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>


// ----------------------------------------------------------------------
//	Class members definition: BOXCQP
// ----------------------------------------------------------------------


void BOXCQP::solveQP(const double *B, const double *d, const int *LDA, double *x, bool *aSolutionOnBorder)
{
	std::vector<int> IPIV(mN);
	int INFO;
	char trans = 'N';
	
	// initialize parameters
	
	// null lagrangian multiplicators
	dcopy_(&mN, &D0, &I0, mLambda, &I1);
	dcopy_(&mN, &D0, &I0, mMu, &I1);

	// solution of the unconstrained problem
	memcpy(x, d, mN*sizeof(double));
	dscal_(&mN, &minus_one, x, &I1);
	
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
		memcpy(&mLHS[i*mN], &B[i**LDA], mN*sizeof(double));
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
	
	// If the solution is already valid, it means it is within the bounds 
	*aSolutionOnBorder = convergenceReached;
	
	// main loop
	for(size_t step(0); !convergenceReached; ++step)
	{
		// --- update the sets
		updateSets(x);
		
#ifdef USE_SUBMATRIX_QP
		
		std::vector<int>::iterator it, jt;
		int i, j, k, Nsub;
		
		// --- update parameters
		
		for(it=mListLset.begin(); it!=mListLset.end(); ++it)
		{
			i = *it;
			x[i]	= ma[i];
			mMu[i]	= 0.0;
		}
		for(it=mListUset.begin(); it!=mListUset.end(); ++it)
		{
			i = *it;
			x[i]		= mb[i];
			mLambda[i]	= 0.0;
		}
		for(it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			mMu[i] 		= 0.0;
			mLambda[i] 	= 0.0;
		}
		
		// --- solve new linear system to get first x
		
		// setup the left hand side matrix and right hand side vector
		Nsub = 0;
		k = 0;
		for(it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			// LHS
			for(jt=mListSset.begin(); jt!=mListSset.end(); ++jt)
			{
				j = *jt;
				mLHS[k] = B[i**LDA+j];
				++k;
			}
			// RHS
			// local variable
			double rhs_tmp( - d[i] );
			for(jt=mListLset.begin(); jt!=mListLset.end(); ++jt)
			{
				j = *jt;
				rhs_tmp -= B[i**LDA + j] * ma[j];
			}
			for(jt=mListUset.begin(); jt!=mListUset.end(); ++jt)
			{
				j = *jt;
				rhs_tmp -= B[i**LDA + j] * mb[j];
			}
			mRHS[Nsub] = rhs_tmp;
			++Nsub;
		}
		
		// solve linear subsystem
		dgesv(&Nsub, &I1, mLHS, &Nsub, &IPIV[0], mRHS, &Nsub, &INFO);
		if(INFO != 0)
		std::cout << "Error: couldn't solve the linear system in BOXCQP. INFO: " << INFO << std::endl;
		
		// update parameters
		// x
		k = 0;
		for(it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			x[i] = mRHS[k];
			++k;
		}
		// lambda
		for(it=mListLset.begin(); it!=mListLset.end(); ++it)
		{
			i = *it;
			mLambda[i] = ddot_(&mN, &B[i**LDA], &I1, x, &I1) + d[i];
		}
		// mu
		for(it=mListUset.begin(); it!=mListUset.end(); ++it)
		{
			i = *it;
			mMu[i] = -ddot_(&mN, &B[i**LDA], &I1, x, &I1) - d[i];
		}
		
#else // NDEF USE_SUBMATRIX_QP
		
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
					dcopy_(&mN, &D0, &I0, &mLHS[i*mN], &I1);
					mLHS[i*(mN+1)] = -1.0;
					break;
				case USET:	
					dcopy_(&mN, &D0, &I0, &mLHS[i*mN], &I1);
					mLHS[i*(mN+1)] = 1.0;
					break;
				case SSET:
					dcopy_(&mN, &B[i**LDA], &I1, &mLHS[i*mN], &I1);
					break;
			};
		}
		
		// setup the right hand side vector
		memcpy(mRHS, d, mN*sizeof(double));
		dgemv_(&trans, &mN, &mN, &D1, B, LDA, mx_known, &I1, &D1, mRHS, &I1);
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
#endif // USE_SUBMATRIX_QP
		
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
#ifdef USE_SUBMATRIX_QP
	mSpace.resize(mN*mN + 3*mN);
	mListLset.reserve(mN);
	mListUset.reserve(mN);
	mListSset.reserve(mN);
#else
	mSpace.resize(mN*mN + 6*mN);
#endif
	
	mLambda = &mSpace[0];
	mMu 	= mLambda + mN;
#ifndef USE_SUBMATRIX_QP
	mx_known 		= mMu + mN;
	mMu_known 		= mx_known + mN;
	mLambda_known 	= mMu_known + mN;
	mRHS	= mLambda_known + mN;
#else
	mRHS	= mMu + mN;
#endif
	mLHS	= mRHS + mN;
}


// ----------------------------------------------------------------------
void BOXCQP::updateSets(double *ax)
{
#ifdef USE_SUBMATRIX_QP
	mListLset.clear();
	mListUset.clear();
	mListSset.clear();
	
	for(size_t i(0); i<mN; ++i)
	{
		if(  (ax[i] < ma[i])   ||   (ax[i] == ma[i] && mLambda[i] >= 0.0) )
		{
			mSets[i] = LSET;
			mListLset.push_back(i);
		}
		else if( (ax[i] > mb[i])   ||   (ax[i] == mb[i] && mMu[i] >= 0.0) )
		{
			mSets[i] = USET;
			mListUset.push_back(i);
		}
		else
		{
			mSets[i] = SSET;
			mListSset.push_back(i);
		}
	}
#else
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
#endif // USE_SUBMATRIX_QP
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
		bool QPsolutionOnBorder;
		mQPsolver->solveQP(mHessian, mGradient, &mN, mP, &QPsolutionOnBorder);
		
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Line Search..." << std::endl;
		
		
		alpha = 1.0;
#if 0
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
			if(x[i] < mLowerBound[i])
				x[i] = mLowerBound[i];
			if(x[i] > mUpperBound[i])
				x[i] = mUpperBound[i];
		}
		
		if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Step length found:" << alpha << std::endl;
		
		//*f = evaluateFunction(x, mTrace);
		
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
		convergenceReached =   fabs(f_prev - *f) < mAbsoluteError
							|| mStep >= mMaxIterations;
	}						
}


// ----------------------------------------------------------------------
double OptSQP::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, size_vect);
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if(mStopIfBigger && f >= mThreshold) throw FastCodeMLEarlyStopLRT();
	
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
	double c1 (2e-1);
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
	if( a+eh >= 1.0 ) {eh = -eh;}
	phi_a_prime = (phi - evaluateFunctionForLineSearch(x, a+eh))/eh;
	
	iter = 0;
	if(phi_a_prime < 0.0 && a != *aalpha)
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
		if(phi_prev < phi)
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

#else
// ----------------------------------------------------------------------
void OptSQP::lineSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	double c1 (2e-1), c2 (0.1);
	double amax = *aalpha;
	double phi_0_prime = ddot_(&mN, mP, &I1, mGradient, &I1);
	double phi_0 = *f, phi, phi_prev;
	double a_prev = 0.0;//*aalpha;
	double phi_a_prime;
	
	double a = randFrom0to1();
		
	if(mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi = phi_prev = phi_0;
	
	double sigma, sigma_bas;
	int maxIter, iter=0;
	
	maxIter = static_cast<int> (ceil( 3.*log(mN+10) ));
	sigma_bas 	= pow(1e-2, 1./static_cast<double>(maxIter));
	
	while(iter < maxIter)
	{
		++iter;
		phi_prev = phi;
		phi = evaluateFunctionForLineSearch(x, a);
		std::cout << "DEBUG LINE SEARCH: phi = " << phi << " for a = " << a << std::endl; 
		
		if(phi > phi_0 + c1*a*phi_0_prime || ((phi > phi_prev) && (iter > 1)) )
		{
			a = zoom(a_prev, a, x, phi_0, phi_0_prime, phi_prev, c1, c2);
			break;
		}
		
		// compute the derivative at point a
		double eh = sqrt(DBL_EPSILON);
		if( a+eh >= 1.0 ) {eh = -eh;}
		phi_a_prime = (evaluateFunctionForLineSearch(x, a+eh) - phi)/eh;
		
		if(fabs(phi_a_prime) <= -c2*phi_0_prime)
		{
			// wolfe conditions are satisfied, stop searching
			break;
		}
		if(phi_a_prime >= 0.0)
		{
			std::cout << "case 3, phi_a_prime: " << phi_a_prime << std::endl;
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
double OptSQP::zoom(double alo, double ahi, double *x, const double& phi_0, const double& phi_0_prime, const double& phi_lo, const double& c1, const double& c2)
{
	double a, phi, phi_a_prime;
	double philo = phi_lo;
	a = 0.5*(alo+ahi);
	while( fabs(ahi-alo) > 0.01 )
	{
		double tmp = 0.5;
		a = tmp*alo + (1.-tmp)*ahi;
		phi = evaluateFunctionForLineSearch(x, a);
		std::cout << "DEBUG ZOOM: phi = " << phi << " for a = " << a << " alo: " << alo << " ahi: " << ahi << " philo: " << philo  << std::endl; 
		if(phi > phi_0 + a*c1*phi_0_prime || phi > philo)
		{
			ahi = a;
		}	
		else
		{
			double eh = sqrt(DBL_EPSILON);
			if( a+eh >= 1.0 ) {eh = -eh;}
			phi_a_prime = (evaluateFunctionForLineSearch(x, a+eh) - phi)/eh;
			
			if(fabs(phi_a_prime) <= -c2*phi_0_prime)
				return a;
				
			if(phi_a_prime*(ahi-alo) >= 0.0)
				ahi = alo;
				
			alo = a;
			philo = phi;
		}
	}
	
	return a;
}
#endif // STRONG_WOLFE_LINE_SEARCH

