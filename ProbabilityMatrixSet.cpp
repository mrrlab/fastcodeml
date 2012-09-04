
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
											 double aSbg,
											 double aSfg,
											 const  std::vector<double>& aParams)
{
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
}

void ProbabilityMatrixSet::computePartialMatrixSetH0(const TransitionMatrix& aQw0,
													 const TransitionMatrix& aQ1,
													 double aSbg,
													 double aSfg,
													 const std::vector<double>& aParams,
													 size_t aBranch)
{
	// Compute the two matrices at the new t value
	const double t = (static_cast<int>(aBranch) == mFgBranch) ? aParams[aBranch]/aSfg : aParams[aBranch]/aSbg;

#ifdef _MSC_VER
	#pragma omp parallel sections default(none) shared(aBranch, aQw0, aQ1)
#else
	#pragma omp parallel sections default(shared)
#endif
	{
		#pragma omp section
		{
			// Save the value and compute the value for the new branch length
			memcpy(mSaveQw0, &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], t);
		}
		#pragma omp section
		{
			// Save the value and compute the value for the new branch length
			memcpy(mSaveQ1, &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
			aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], t);
		}
	}
}

void ProbabilityMatrixSet::restoreSavedMatrixH0(size_t aBranch)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQw0, N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQ1,  N*N*sizeof(double));
}


void ProbabilityMatrixSet::computeMatrixSetH1(const  TransitionMatrix& aQw0,
											  const  TransitionMatrix& aQ1,
											  const  TransitionMatrix& aQw2,
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

	aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT], params[mFgBranch]/aSfg);
}

void ProbabilityMatrixSet::computePartialMatrixSetH1(const  TransitionMatrix& aQw0,
													 const  TransitionMatrix& aQ1,
													 const  TransitionMatrix& aQw2,
													 double aSbg,
													 double aSfg,
													 const std::vector<double>& aParams,
													 size_t aBranch)
{
	// Compute the two matrices at the new t value
	const double t = (static_cast<int>(aBranch) == mFgBranch) ? aParams[aBranch]/aSfg : aParams[aBranch]/aSbg;

#ifdef _MSC_VER
	#pragma omp parallel sections default(none) shared(aBranch, aQw0, aQ1, aQw2)
#else
	#pragma omp parallel sections default(shared)
#endif
	{
		#pragma omp section
		{
			// Save the value and compute the value for the new branch length
			memcpy(mSaveQw0, &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
			aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], t);
		}
		#pragma omp section
		{
			// Save the value and compute the value for the new branch length
			memcpy(mSaveQ1,  &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
			aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], t);
		}
		#pragma omp section
		{
			if(static_cast<int>(aBranch) == mFgBranch)
			{
				memcpy(mSaveQw2, &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], N*N*sizeof(double));
				aQw2.computeFullTransitionMatrix(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], t);
			}
		}
	}
}


void ProbabilityMatrixSet::restoreSavedMatrixH1(size_t aBranch)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQw0, N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQ1,  N*N*sizeof(double));
	if(static_cast<int>(aBranch) == mFgBranch)
	{
		memcpy(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], mSaveQw2, N*N*sizeof(double));
	}
}

