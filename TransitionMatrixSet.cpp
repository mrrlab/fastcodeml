
#include "TransitionMatrixSet.h"


void TransitionMatrixSet::computeMatrixSetH0(const TransitionMatrix& aQw0,
											 const TransitionMatrix& aQ1,
											 double aSbg,
											 double aSfg,
											 unsigned int aFgBranch,
											 const std::vector<double>& aParams)
{
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aFgBranch, aQw0, aQ1, aSbg, aSfg, aParams)
#else
	#pragma omp parallel for default(shared)
#endif
	for(int branch=0; branch < (int)mNumMatrices; ++branch)
	{
		if(branch == (int)aFgBranch)
		{
			aQw0.computeFullTransitionMatrix(mMatrixSpace+0*mNumMatrices*N*N+branch*N*N, aParams[branch]/aSfg);
			aQ1.computeFullTransitionMatrix(mMatrixSpace+1*mNumMatrices*N*N+branch*N*N,  aParams[branch]/aSfg);
		}
		else
		{
			aQw0.computeFullTransitionMatrix(mMatrixSpace+0*mNumMatrices*N*N+branch*N*N, aParams[branch]/aSbg);
			aQ1.computeFullTransitionMatrix(mMatrixSpace+1*mNumMatrices*N*N+branch*N*N,  aParams[branch]/aSbg);
		}
	}

	//unsigned int k;
	//for(k=0; k < mNumMatrices*2; ++k) mMatrices[k] = mMatrixSpace+k*N*N;
	//for(k=0; k < mNumMatrices; ++k)   mMatrices[k+mNumMatrices*2] = mMatrixSpace+k*N*N;
	
	for(unsigned int k=0; k < mNumMatrices; ++k)
	{
		mMatrices[k+mNumMatrices*0] = mMatrixSpace+k*N*N;
		mMatrices[k+mNumMatrices*1] = mMatrixSpace+(mNumMatrices+k)*N*N;
		mMatrices[k+mNumMatrices*2] = mMatrixSpace+k*N*N;
	}
	mMatrices[aFgBranch+mNumMatrices*2] = mMatrixSpace+(mNumMatrices+aFgBranch)*N*N;
}


void TransitionMatrixSet::computeMatrixSetH1(const TransitionMatrix& aQw0,
											 const TransitionMatrix& aQ1,
											 const TransitionMatrix& aQw2,
											 double aSbg,
											 double aSfg,
											 unsigned int aFgBranch,
											 const std::vector<double>& aParams)
{
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aFgBranch, aQw0, aQ1, aQw2, aSbg, aSfg, aParams)
#else
	#pragma omp parallel for default(shared)
#endif
	for(int branch=0; branch < (int)mNumMatrices; ++branch)
	{
		if(branch == (int)aFgBranch)
		{
			aQw0.computeFullTransitionMatrix(mMatrixSpace+0*mNumMatrices*N*N+branch*N*N, aParams[branch]/aSfg);
			aQ1.computeFullTransitionMatrix(mMatrixSpace+1*mNumMatrices*N*N+branch*N*N,  aParams[branch]/aSfg);
		}
		else
		{
			aQw0.computeFullTransitionMatrix(mMatrixSpace+0*mNumMatrices*N*N+branch*N*N, aParams[branch]/aSbg);
			aQ1.computeFullTransitionMatrix(mMatrixSpace+1*mNumMatrices*N*N+branch*N*N,  aParams[branch]/aSbg);
		}
	}

	for(unsigned int k=0; k < mNumMatrices; ++k)
	{
		mMatrices[k+mNumMatrices*0] = mMatrixSpace+k*N*N;
		mMatrices[k+mNumMatrices*1] = mMatrixSpace+(mNumMatrices+k)*N*N;
		mMatrices[k+mNumMatrices*2] = mMatrixSpace+k*N*N;
		mMatrices[k+mNumMatrices*3] = mMatrixSpace+(mNumMatrices+k)*N*N;
	}

	aQw2.computeFullTransitionMatrix(mMatrixSpace+2*mNumMatrices*N*N+aFgBranch*N*N, aParams[aFgBranch]/aSfg);
	mMatrices[aFgBranch+mNumMatrices*2] = mMatrixSpace+(2*mNumMatrices+aFgBranch)*N*N;
	mMatrices[aFgBranch+mNumMatrices*3] = mMatrixSpace+(2*mNumMatrices+aFgBranch)*N*N;
}

