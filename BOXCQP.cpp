
#include "blas.h"
#include "lapack.h"
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include "BOXCQP.h"

// ----------------------------------------------------------------------
//	Class members definition: BOXCQP
// ----------------------------------------------------------------------

bool BOXCQP::solveQP(const double *aB, const double *aD, const int *aLDB, double *aX, bool *aSolutionOnBorder, double *aUnconstrainedDirection)
{
	int info;
	#ifdef USE_SCALING_COND
	const char fact = 'E';
	const char uplo = 'U';
	char equed;
	double rcond, ferr, berr;
	#else
	std::vector<int> IPIV(mN);
	#endif
	
	// initialize parameters
	
	// null lagrangian multiplicators
	dcopy_(&mN, &D0, &I0, mLambda, &I1);
	dcopy_(&mN, &D0, &I0, mMu, &I1);

	// solution of the unconstrained problem
	memcpy(aX, aD, mN*sizeof(double));
	dscal_(&mN, &minus_one, aX, &I1);
	
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
		memcpy(&mLHS[i*mN], &aB[i**aLDB], mN*sizeof(double));
	
	#ifdef USE_SCALING_COND
	dposvx_(&fact, &uplo, &mN, &I1, mLHS, &mN 
		   ,mSubmatrixFact, &mN, &equed
		   ,mDiagScaling, aX, &mN
		   ,mSubSolution, &mN, &rcond, &ferr, &berr
		   ,&mDWorkSpace[0], &mIWorkspace[0], &info);
	memcpy(aX, mSubSolution, mN*sizeof(double));
	if (info != 0 && info != (mN+1))
		std::cout << "Error: couldn't solve the linear system in BOXCQP. info: " << info << std::endl;
	else
		std::cout << "\tequilibrated: " << equed << ", cond = " << 1.0/rcond << ", ferr = " << std::scientific << ferr << ", berr = " << berr << std::fixed << std::endl;
	#else
	dgesv_(&mN, &I1, mLHS, &mN, &IPIV[0], aX, &mN, &info);
	#endif
	
	if (info != 0)
		std::cout << "Error: couldn't solve the initial linear system in BOXCQP. info: " << info << std::endl;
	
	if (aUnconstrainedDirection != NULL)
	{
		memcpy(aUnconstrainedDirection, aX, mN*sizeof(double));
	}
	
	// verify if the solution has converged
	bool convergence_reached(true);
	#pragma omp parallel for reduction(&&: convergence_reached)
	for (int i=0; i<mN; ++i)
	{
		convergence_reached = convergence_reached && (mLowerBounds[i] <= aX[i] && aX[i] <= mUpperBounds[i] );
	}
	
	// If the solution is already valid, it means it is within the bounds 
	*aSolutionOnBorder = convergence_reached;
	
	// main loop
	size_t max_step = (mN < 100) ? 100+mN*mN : mN*10;
	for (size_t step(0); !convergence_reached && step < max_step; ++step)
	{
		// --- update the sets
		updateSets(aX);
		
#ifdef USE_SUBMATRIX_QP
		
		std::vector<int>::iterator it, jt;
		int i, j, k, Nsub;
		
		// --- update parameters
		
		for (it=mListLset.begin(); it!=mListLset.end(); ++it)
		{
			i = *it;
			aX[i]	= mLowerBounds[i];
			mMu[i]	= 0.0;
		}
		for (it=mListUset.begin(); it!=mListUset.end(); ++it)
		{
			i = *it;
			aX[i]		= mUpperBounds[i];
			mLambda[i]	= 0.0;
		}
		for (it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			mMu[i] 		= 0.0;
			mLambda[i] 	= 0.0;
		}
		
		// --- solve new linear system to get first x
		
		// setup the left hand side matrix and right hand side vector
		Nsub = 0;
		k = 0;
		for (it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			// LHS (left-hand-side submatrix, positive definite)
			for (jt=mListSset.begin(); jt!=mListSset.end(); ++jt)
			{
				j = *jt;
				mLHS[k] = aB[i**aLDB+j];
				++k;
			}
			// RHS (right-hand-side subvector)
			// local variable
			double rhs_tmp( - aD[i] );
			for (jt=mListLset.begin(); jt!=mListLset.end(); ++jt)
			{
				j = *jt;
				rhs_tmp -= aB[i**aLDB + j] * mLowerBounds[j];
			}
			for (jt=mListUset.begin(); jt!=mListUset.end(); ++jt)
			{
				j = *jt;
				rhs_tmp -= aB[i**aLDB + j] * mUpperBounds[j];
			}
			mRHS[Nsub] = rhs_tmp;
			++Nsub;
		}
		
		// solve linear subsystem
		#ifdef USE_SCALING_COND
		dposvx_(&fact, &uplo, &Nsub, &I1, mLHS, &Nsub 
			   ,mSubmatrixFact, &Nsub, &equed
			   ,mDiagScaling, mRHS, &Nsub
			   ,mSubSolution, &Nsub, &rcond, &ferr, &berr
			   ,&mDWorkSpace[0], &mIWorkspace[0], &info);
		memcpy(mRHS, mSubSolution, Nsub*sizeof(double));
		if (info != 0 && info != (Nsub+1))
			std::cout << "Error: couldn't solve the linear system in BOXCQP. info: " << info << std::endl;
		#else
		dgesv_(&Nsub, &I1, mLHS, &Nsub, &IPIV[0], mRHS, &Nsub, &info);
		if (info != 0)
			std::cout << "Error: couldn't solve the linear system in BOXCQP. info: " << info << std::endl;
		#endif
		
		// update parameters
		// x
		k = 0;
		for (it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			aX[i] = mRHS[k];
			++k;
		}
		// lambda
		for (it=mListLset.begin(); it!=mListLset.end(); ++it)
		{
			i = *it;
			mLambda[i] = ddot_(&mN, &aB[i**aLDB], &I1, aX, &I1) + aD[i];
		}
		// mu
		for (it=mListUset.begin(); it!=mListUset.end(); ++it)
		{
			i = *it;
			mMu[i] = -ddot_(&mN, &aB[i**aLDB], &I1, aX, &I1) - aD[i];
		}
		
#else // NDEF USE_SUBMATRIX_QP
		
		// --- update parameters
		
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			switch (mSets[i])
			{
				case LSET:
					aX[i] 		= mLowerBounds[i];
					mMu[i] 		= 0.0;
					break;
				
				case USET:
					aX[i] 		= mUpperBounds[i];
					mLambda[i] 	= 0.0;
					break;
				
				case SSET:
					mMu[i] 		= 0.0;
					mLambda[i] 	= 0.0;
					break;
			};
		}
		
		// setup the known vectors; set their values to zero where we don't know it. 
		memcpy(mLambdaKnown, 	mLambda, 	mN*sizeof(double));
		memcpy(mMuKnown, 		mMu, 		mN*sizeof(double));
		memcpy(mXKnown, 		aX, 		mN*sizeof(double));
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			switch (mSets[i])
			{
				case LSET:
					mLambdaKnown[i] = 0.0;
					break;
				
				case USET:
					mMuKnown[i] = 0.0;
					break;
				
				case SSET:
					mXKnown[i] = 0.0;
					break;
			};
		}

		
		// --- solve new linear system
				
		// setup the left hand side matrix
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			switch (mSets[i])
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
					dcopy_(&mN, &aB[i**aLDB], &I1, &mLHS[i*mN], &I1);
					break;
			};
		}
		
		// setup the right hand side vector
		memcpy(mRHS, aD, mN*sizeof(double));
		dgemv_("N", &mN, &mN, &D1, aB, aLDB, mXKnown, &I1, &D1, mRHS, &I1);
		daxpy_(&mN, &D1, mMuKnown, &I1, mRHS, &I1);
		daxpy_(&mN, &minus_one, mLambdaKnown, &I1, mRHS, &I1);
		dscal_(&mN, &minus_one, mRHS, &I1);
		
		// solve linear system
		dgesv_(&mN, &I1, mLHS, &mN, &IPIV[0], mRHS, &mN, &info);
		if (info != 0)
		std::cout << "Error: couldn't solve the linear system in BOXCQP. info: " << info << std::endl;
		
		// update solutions
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
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
					aX[i] = mRHS[i];
					break;
			};
		}
#endif // USE_SUBMATRIX_QP
		
		// --- verify validity of the solution
		convergence_reached = true;
		#pragma omp parallel for reduction(&&: convergence_reached)
		for (int i=0; i<mN; ++i)
		{
			// local variable
			bool current_variable_correct = false;
			switch (mSets[i])
			{
				case LSET:
					current_variable_correct = mLambda[i] >= 0.0;
					break;
				case USET:
					current_variable_correct = mMu[i] >= 0.0;
					break;
				case SSET:
					current_variable_correct = (aX[i] >= mLowerBounds[i]) && (aX[i] <= mUpperBounds[i]);
					break;
			};
			convergence_reached = convergence_reached && current_variable_correct;
		}	
	}
	
	if (!convergence_reached)
	{
		std::cout << "\t\tBOXCQP: convergence not reached." << std::endl;
	}
	else
	{
		#pragma omp parallel for
		for (int i=0; i<mN; ++i)
		{
			double x = aX[i];
			double l = mLowerBounds[i];
			double u = mUpperBounds[i];
			x = (x < l) ? l:x;
			x = (x > u) ? u:x;
			aX[i] = x;
		}
	}
	
	return convergence_reached;
}


// ----------------------------------------------------------------------
void BOXCQP::allocateMemory(void)
{
	mSets.resize(mN);
#ifdef USE_SUBMATRIX_QP
	mSpace.resize(2*mN*mN + 5*mN);
	mListLset.reserve(mN);
	mListUset.reserve(mN);
	mListSset.reserve(mN);
	
	mDWorkSpace.resize(3*mN);
	mIWorkspace.resize(mN);
#else
	mSpace.resize(mN*mN + 6*mN);
#endif
	
	mLambda = &mSpace[0];
	mMu 	= mLambda + mN;
#ifndef USE_SUBMATRIX_QP
	mXKnown 		= mMu + mN;
	mMuKnown 		= mXKnown + mN;
	mLambdaKnown 	= mMuKnown + mN;
	mRHS			= mLambdaKnown + mN;
#else
	mSubmatrixFact	= mMu + mN;
	mDiagScaling	= mSubmatrixFact + mN*mN;
	mSubSolution	= mDiagScaling + mN;
	mRHS			= mSubSolution + mN;
#endif
	mLHS	= mRHS + mN;
}


// ----------------------------------------------------------------------
void BOXCQP::updateSets(double *aX)
{
#ifdef USE_SUBMATRIX_QP
	mListLset.clear();
	mListUset.clear();
	mListSset.clear();
	
	for(int i=0; i<mN; ++i)
	{
		if ( (aX[i] < mLowerBounds[i])   ||   ( fabs(aX[i] - mLowerBounds[i]) < 1e-8 && mLambda[i] >= 0.0) )
		{
			mSets[i] = LSET;
			mListLset.push_back(i);
		}
		else if ( (aX[i] > mUpperBounds[i])   ||   ( fabs(aX[i] - mUpperBounds[i]) < 1e-8  && mMu[i] >= 0.0) )
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
	for(int i=0; i<mN; ++i)
	{
		if ( (aX[i] < mLowerBounds[i])   ||   ( fabs(aX[i] - mLowerBounds[i]) < 1e-8 && mLambda[i] >= 0.0) )
			mSets[i] = LSET;
		else if ( (aX[i] > mUpperBounds[i])   ||   ( fabs(aX[i] - mUpperBounds[i]) < 1e-8 && mMu[i] >= 0.0) )
			mSets[i] = USET;
		else
			mSets[i] = SSET;
	}
#endif // USE_SUBMATRIX_QP
}
// ----------------------------------------------------------------------


