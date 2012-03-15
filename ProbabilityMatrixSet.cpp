
#include "ProbabilityMatrixSet.h"

void ProbabilityMatrixSet::initForH0(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
		//mMatrices[branch+num_matrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		//mMatrices[branch+num_matrices*1] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];
		//mMatrices[branch+num_matrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];

		mMatrices[branch+num_matrices*0] = mMatrices[branch+num_matrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+num_matrices*1] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];
	}
	mMatrices[mFgBranch+num_matrices*2] = &mMatrixSpace[(num_matrices+mFgBranch)*MATRIX_SLOT];
}


void ProbabilityMatrixSet::initForH1(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
		mMatrices[branch+num_matrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+num_matrices*1] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];

		// These will be overwritten outside the loop, so skip them
		if(branch != mFgBranch)
		{
			mMatrices[branch+num_matrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
			mMatrices[branch+num_matrices*3] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];
		}
		else
		{
			mMatrices[mFgBranch+num_matrices*2] = mMatrices[mFgBranch+num_matrices*3] = &mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT];
		}
	}
	//mMatrices[mFgBranch+num_matrices*2] = &mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT];
	//mMatrices[mFgBranch+num_matrices*3] = &mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT];
}

#define PREV_TIME

void ProbabilityMatrixSet::computeMatrixSetH0(const TransitionMatrix& aQw0,
											 const  TransitionMatrix& aQ1,
											 bool   aAnyMatrixChanged,
											 double aSbg,
											 double aSfg,
											 const  std::vector<double>& aParams)
{
	int num_matrices = mNumMatrices;
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQw0, aQ1, aSbg, aSfg, aParams, num_matrices, aAnyMatrixChanged) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(static)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
		const double t = (branch == mFgBranch) ? aParams[branch]/aSfg : aParams[branch]/aSbg;

#ifdef PREV_TIME
		if(aAnyMatrixChanged || isDifferent(t, mPrevTime[branch]))
		{
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
#if 0
		if(aChangedQ1)
		{
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
		else if(aChangedQw0)
		{
			if(isDifferent(t, mPrevTime[branch])) aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
		else if(isDifferent(t, mPrevTime[branch]))
		{
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
#endif
#else
		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
#endif
		//mMatrices[branch+mNumMatrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*1] = &mMatrixSpace[(mNumMatrices+branch)*MATRIX_SLOT];
		//mMatrices[branch+mNumMatrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
	}
	//mMatrices[mFgBranch+mNumMatrices*2] = &mMatrixSpace[(mNumMatrices+mFgBranch)*MATRIX_SLOT];
}


void ProbabilityMatrixSet::computeMatrixSetH1(const TransitionMatrix& aQw0,
											 const  TransitionMatrix& aQ1,
											 const  TransitionMatrix& aQw2,
											 bool   aChangedQw2,
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

#ifdef PREV_TIME
	const double t = aParams[mFgBranch]/aSfg;
	if(aChangedQw2 || isDifferent(t, mPrevTime[mFgBranch]))
	{
		aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], t);
		mPrevTime[mFgBranch] = t;
	}
#else
	aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], aParams[mFgBranch]/aSfg);
#endif
	//mMatrices[mFgBranch+mNumMatrices*2] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
	//mMatrices[mFgBranch+mNumMatrices*3] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
}

