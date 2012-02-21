
#include "ProbabilityMatrixSet.h"

void ProbabilityMatrixSet::initForH0(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
		mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
	}
	mMatrices[mFgBranch+mNumMatrices*2] = &mMatrixSpace[(mNumMatrices+mFgBranch)*MATRIX_SLOT];
}


void ProbabilityMatrixSet::initForH1(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
		mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];

		// These will be overwritten outside the loop, so skip them
		if(branch != mFgBranch)
		{
			mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
			mMatrices[branch+mNumMatrices*3] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		}
	}
	mMatrices[mFgBranch+mNumMatrices*2] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
	mMatrices[mFgBranch+mNumMatrices*3] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
}


void ProbabilityMatrixSet::computeMatrixSetH0(const TransitionMatrix& aQw0,
											 const  TransitionMatrix& aQ1,
											 double aSbg,
											 double aSfg,
											 const  std::vector<double>& aParams)
{
	int num_matrices = mNumMatrices;
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQw0, aQ1, aSbg, aSfg, aParams, num_matrices) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(static)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
		const double t = (branch == mFgBranch) ? aParams[branch]/aSfg : aParams[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);

		//mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
	}
	//mMatrices[mFgBranch+mNumMatrices*2] = &mMatrixSpace[(mNumMatrices+mFgBranch)*MATRIX_SLOT];
}


void ProbabilityMatrixSet::computeMatrixSetH1(const TransitionMatrix& aQw0,
											 const  TransitionMatrix& aQ1,
											 const  TransitionMatrix& aQw2,
											 double aSbg,
											 double aSfg,
											 const  std::vector<double>& aParams)
{
	int num_matrices = mNumMatrices;
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQw0, aQ1, aSbg, aSfg, aParams, num_matrices) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(static)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
		const double t = (branch == mFgBranch) ? aParams[branch]/aSfg : aParams[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);

		//mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*3] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
	}

	aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], aParams[mFgBranch]/aSfg);
	//mMatrices[mFgBranch+mNumMatrices*2] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
	//mMatrices[mFgBranch+mNumMatrices*3] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
}

