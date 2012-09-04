
#include "ProbabilityMatrixSet.h"
#ifdef SAVE_OCTAVE
#include <cstdio>
extern void SaveToOctave(const double *CVariable, char *OctaveVariable, FILE *FilePointer, int Rows, int Columns);
#endif

void ProbabilityMatrixSet::initForH0(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
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
}

//#define PREV_TIME

void ProbabilityMatrixSet::computeMatrixSetH0(const TransitionMatrix& aQw0,
											 const  TransitionMatrix& aQ1,
											 bool   aAnyMatrixChanged,
											 double aSbg,
											 double aSfg,
											 const  std::vector<double>& aParams)
{
	const int num_matrices = mNumMatrices;
	const double* params = &aParams[0];

#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQw0, aQ1, aSbg, aSfg, params, num_matrices, aAnyMatrixChanged) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(runtime)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
		const double t = (branch == mFgBranch) ? params[branch]/aSfg : params[branch]/aSbg;

#ifdef PREV_TIME
		if(aAnyMatrixChanged || isDifferent(t, mPrevTime[branch]))
		{
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
#ifdef SAVE_OCTAVE
FILE *fp = fopen("m.oct", "a");
SaveToOctave(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], "Yw0", fp, 61, 61);
SaveToOctave(&mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], "Y1", fp, 61, 61);
SaveToOctave(&aParams[branch], "T", fp, 1, 1);
SaveToOctave(&t, "TS", fp, 1, 1);
fclose(fp);
#endif
#if 0
		if(aChangedQ1)
		{
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
		else if(aChangedQw0)
		{
			if(isDifferent(t, mPrevTime[branch])) aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
		else if(isDifferent(t, mPrevTime[branch]))
		{
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			aQ1.computeFullTransitionMatrix(&mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
			mPrevTime[branch] = t;
		}
#endif
#else
		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
#endif
	}
}

void ProbabilityMatrixSet::computePartialMatrixSetH0(const TransitionMatrix& aQw0,
						    const TransitionMatrix& aQ1,
							double aSbg,
							double aSfg,
						    const std::vector<double>& aParams,
							unsigned int aBranch)
{
	// Save the value
	memcpy(mSave1, &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
	memcpy(mSave2, &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));

	// Compute the two matrices at the new t value
	const double* params = &aParams[0];

	const double t = (static_cast<int>(aBranch) == mFgBranch) ? params[aBranch]/aSfg : params[aBranch]/aSbg;

	aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], t);
	aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], t);
}

void ProbabilityMatrixSet::restoreSavedMatrixH0(unsigned int aBranch)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSave1, N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSave2, N*N*sizeof(double));
}

void ProbabilityMatrixSet::computeMatrixSetH1(const  TransitionMatrix& aQw0,
											  const  TransitionMatrix& aQ1,
											  const  TransitionMatrix& aQw2,
											  bool   aChangedQw2,
											  double aSbg,
											  double aSfg,
											  const  std::vector<double>& aParams)
{
	// To speedup access to variables
	const int num_matrices = mNumMatrices;
	const double* params = &aParams[0];

#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQw0, aQ1, aSbg, aSfg, params, num_matrices) schedule(static)
#else
	#pragma omp parallel for default(shared) schedule(runtime)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
		const double t = (branch == mFgBranch) ? params[branch]/aSfg : params[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
	}

#ifdef PREV_TIME
	const double t = aParams[mFgBranch]/aSfg;
	if(aChangedQw2 || isDifferent(t, mPrevTime[mFgBranch]))
	{
		aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT], t);
		mPrevTime[mFgBranch] = t;
	}
#else
	aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT], params[mFgBranch]/aSfg);
#endif
}

