
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
	std::vector<int> IPIV(mN);
	int INFO;
	//char trans = 'N';
	
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
	
	dgesv(&mN, &I1, mLHS, &mN, &IPIV[0], aX, &mN, &INFO);
	
	if (INFO != 0)
		std::cout << "Error: couldn't solve the initial linear system in BOXCQP. INFO: " << INFO << std::endl;
	
	if (aUnconstrainedDirection != NULL)
	{
		memcpy(aUnconstrainedDirection, aX, mN*sizeof(double));
	}
	
	// verify if the solution is valide
	bool convergenceReached(true);
	#pragma omp parallel for reduction(&&: convergenceReached)
	for(int i=0; i<mN; ++i)
	{
		convergenceReached = convergenceReached && (mLowerBounds[i] <= aX[i] && aX[i] <= mUpperBounds[i] );
	}
	
	// If the solution is already valid, it means it is within the bounds 
	*aSolutionOnBorder = convergenceReached;
	
	// main loop
	size_t max_step = 100+mN*mN;
	for(size_t step(0); !convergenceReached && step < max_step; ++step)
	{
		// --- update the sets
		updateSets(aX);
		
#ifdef USE_SUBMATRIX_QP
		
		std::vector<int>::iterator it, jt;
		int i, j, k, Nsub;
		
		// --- update parameters
		
		for(it=mListLset.begin(); it!=mListLset.end(); ++it)
		{
			i = *it;
			aX[i]	= mLowerBounds[i];
			mMu[i]	= 0.0;
		}
		for(it=mListUset.begin(); it!=mListUset.end(); ++it)
		{
			i = *it;
			aX[i]		= mUpperBounds[i];
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
				mLHS[k] = aB[i**aLDB+j];
				++k;
			}
			// RHS
			// local variable
			double rhs_tmp( - aD[i] );
			for(jt=mListLset.begin(); jt!=mListLset.end(); ++jt)
			{
				j = *jt;
				rhs_tmp -= aB[i**aLDB + j] * mLowerBounds[j];
			}
			for(jt=mListUset.begin(); jt!=mListUset.end(); ++jt)
			{
				j = *jt;
				rhs_tmp -= aB[i**aLDB + j] * mUpperBounds[j];
			}
			mRHS[Nsub] = rhs_tmp;
			++Nsub;
		}
		
		// solve linear subsystem
		
		dgesv(&Nsub, &I1, mLHS, &Nsub, &IPIV[0], mRHS, &Nsub, &INFO);
		if (INFO != 0)
		std::cout << "Error: couldn't solve the linear system in BOXCQP. INFO: " << INFO << std::endl;
		
		// update parameters
		// x
		k = 0;
		for(it=mListSset.begin(); it!=mListSset.end(); ++it)
		{
			i = *it;
			aX[i] = mRHS[k];
			++k;
		}
		// lambda
		for(it=mListLset.begin(); it!=mListLset.end(); ++it)
		{
			i = *it;
			mLambda[i] = ddot_(&mN, &aB[i**aLDB], &I1, aX, &I1) + aD[i];
		}
		// mu
		for(it=mListUset.begin(); it!=mListUset.end(); ++it)
		{
			i = *it;
			mMu[i] = -ddot_(&mN, &aB[i**aLDB], &I1, aX, &I1) - aD[i];
		}
		
#else // NDEF USE_SUBMATRIX_QP
		
		// --- update parameters
		
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			switch(mSets[i])
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
		memcpy(mXKnown, 		aX, 			mN*sizeof(double));
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			switch(mSets[i])
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
		for(int i=0; i<mN; ++i)
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
		dgesv(&mN, &I1, mLHS, &mN, &IPIV[0], mRHS, &mN, &INFO);
		if (INFO != 0)
		std::cout << "Error: couldn't solve the linear system in BOXCQP. INFO: " << INFO << std::endl;
		
		// update solutions
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
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
		convergenceReached = true;
		#pragma omp parallel for reduction(&&: convergenceReached)
		for(int i=0; i<mN; ++i)
		{
			// local variable
			bool iVariableCorrect = false;
			switch(mSets[i])
			{
				case LSET:
					iVariableCorrect = mLambda[i] >= 0.0;
					break;
				case USET:
					iVariableCorrect = mMu[i] >= 0.0;
					break;
				case SSET:
					iVariableCorrect = (aX[i] >= mLowerBounds[i]) && (aX[i] <= mUpperBounds[i]);
					break;
			};
			convergenceReached = convergenceReached && iVariableCorrect;
		}	
	}
	
	if (!convergenceReached)
	{
		std::cout << "\t\tBOXCQP: convergence not reached." << std::endl;
	}
	else
	{
		#pragma omp parallel for
		for (int i(0); i<mN; ++i)
		{
			double x = aX[i];
			double l = mLowerBounds[i];
			double u = mUpperBounds[i];
			x = (x < l) ? l:x;
			x = (x > u) ? u:x;
			aX[i] = x;
		}
	}
	
	return convergenceReached;
}


// ----------------------------------------------------------------------
void BOXCQP::allocateMemory(void)
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
	mXKnown 		= mMu + mN;
	mMuKnown 		= mXKnown + mN;
	mLambdaKnown 	= mMuKnown + mN;
	mRHS			= mLambdaKnown + mN;
#else
	mRHS	= mMu + mN;
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


