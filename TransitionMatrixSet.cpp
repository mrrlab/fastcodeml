
#include "TransitionMatrixSet.h"

void TransitionMatrixSet::computeMatrixSet(const TransitionMatrix& aQw0,
											 const TransitionMatrix& aQ1,
											 double aSbg,
											 double aSfg,
											 unsigned int aFgBranch,
											 const std::vector<double>& aParams)
{
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aFgBranch, aQw0, aQ1, aSbg, aSfg, aParams) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(static)
#endif
	for(int branch=0; branch < (int)mNumMatrices; ++branch)
	{
		const double t = (branch == (int)aFgBranch) ? aParams[branch]/aSfg : aParams[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);

		mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
	}
	mMatrices[aFgBranch+mNumMatrices*2] = &mMatrixSpace[(mNumMatrices+aFgBranch)*MATRIX_SLOT];
}


void TransitionMatrixSet::computeMatrixSet(const TransitionMatrix& aQw0,
											 const TransitionMatrix& aQ1,
											 const TransitionMatrix& aQw2,
											 double aSbg,
											 double aSfg,
											 unsigned int aFgBranch,
											 const std::vector<double>& aParams)
{
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aFgBranch, aQw0, aQ1, aSbg, aSfg, aParams) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(static)
#endif
	for(int branch=0; branch < (int)mNumMatrices; ++branch)
	{
		const double t = (branch == (int)aFgBranch) ? aParams[branch]/aSfg : aParams[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);

		mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*3] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
	}

	aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*mNumMatrices+aFgBranch)*MATRIX_SLOT], aParams[aFgBranch]/aSfg);
	mMatrices[aFgBranch+mNumMatrices*2] = &mMatrixSpace[(2*mNumMatrices+aFgBranch)*MATRIX_SLOT];
	mMatrices[aFgBranch+mNumMatrices*3] = &mMatrixSpace[(2*mNumMatrices+aFgBranch)*MATRIX_SLOT];
}

