
#include "BOXCQP.h"

#include "blas.h"
#include "lapack.h"
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

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
	
	if (INFO != 0)
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
		if (INFO != 0)
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
		if (INFO != 0)
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
		if ( (ax[i] < ma[i])   ||   (ax[i] == ma[i] && mLambda[i] >= 0.0) )
		{
			mSets[i] = LSET;
			mListLset.push_back(i);
		}
		else if ( (ax[i] > mb[i])   ||   (ax[i] == mb[i] && mMu[i] >= 0.0) )
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
		if ( (ax[i] < ma[i])   ||   (ax[i] == ma[i] && mLambda[i] >= 0.0) )
			mSets[i] = LSET;
		else if ( (ax[i] > mb[i])   ||   (ax[i] == mb[i] && mMu[i] >= 0.0) )
			mSets[i] = USET;
		else
			mSets[i] = SSET;
	}
#endif // USE_SUBMATRIX_QP
}
// ----------------------------------------------------------------------


