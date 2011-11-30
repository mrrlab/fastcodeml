
#ifndef TRANSITIONMATRIXSET_H
#define TRANSITIONMATRIXSET_H
#include <iostream>

#include <vector>
#include "MatrixSize.h"
#include "TransitionMatrix.h"
#include "AlignedMalloc.h"
#include "CodonFrequencies.h"
#include "MathSupport.h"

#ifdef USE_LAPACK
#include "blas.h"
#endif
#ifdef USE_MKL_VML
#include <mkl_vml_functions.h>
#endif

/// Routines that manage the set of transition matrices for all branches of a tree.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-04-05 (initial version)
///     @version 1.0
///
///
class TransitionMatrixSet
{
public:

	/// Create matrix set
	///
	/// @param[in] aNumMatrices The number of matrices to be managed (is the number of branches of the tree)
	/// @param[in] aNumSets How many sets to allocate (one set is composed by the bg and fg matrices for one of the tree traversals)
	///
	TransitionMatrixSet(unsigned int aNumMatrices, unsigned int aNumSets) : mNumMatrices(aNumMatrices), mNumSets(aNumSets)
	{
		mMatrixSpace  = (double *)alignedMalloc(sizeof(double)*aNumSets*aNumMatrices*MATRIX_SLOT, CACHE_LINE_ALIGN);
		mMatrices     = (double**)alignedMalloc(sizeof(double*)*aNumSets*aNumMatrices, CACHE_LINE_ALIGN);
		CodonFrequencies* cf = CodonFrequencies::getInstance();
		mInvCodonFreq = cf->getInvCodonFrequencies();
	}

	/// Destructor.
	///
	~TransitionMatrixSet()
	{
		alignedFree(mMatrixSpace);
		alignedFree(mMatrices);
	}

	/// Compute the three sets of matrices for the H0 hypothesis
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w1
	///
	///	@param[in] aQw0 The mQw0 transition matrix
	///	@param[in] aQ1 The mQ1 transition matrix
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (k, w0, 02, p0+p1, p0/(p0+p1), w2)
	///
	void computeMatrixSetH0(const TransitionMatrix& aQw0,
						    const TransitionMatrix& aQ1,
							double aSbg,
							double aSfg,
						    unsigned int aFgBranch,
						    const std::vector<double>& aParams);

	/// Compute the four sets of matrices for the H1 hypothesis
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w2
	/// - set 3: w1, w2
	///
	///	@param[in] aQw0 The mQw0 transition matrix
	///	@param[in] aQ1 The mQ1 transition matrix
	///	@param[in] aQw2 The mQw2 transition matrix
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (k, w0, 02, p0+p1, p0/(p0+p1), w2)
	///
	void computeMatrixSetH1(const TransitionMatrix& aQw0,
						    const TransitionMatrix& aQ1,
						    const TransitionMatrix& aQw2,
							double aSbg,
							double aSfg,
						    unsigned int aFgBranch,
						    const std::vector<double>& aParams);
#ifndef NEW_LIKELIHOOD
	///	Multiply the aGin vector by the precomputed exp(Q*t) matrix
	///
	/// @param[in] aSetIdx Which set to use (starts from zero)
	/// @param[in] aBranch Which branch
	/// @param[in] aGin The input vector to be multiplied by the matrix exponential
	/// @param[out] aGout The resulting vector
	///
	void doTransition(unsigned int aSetIdx, unsigned int aBranch, const double* aGin, double* aGout) const
	{
#ifdef USE_LAPACK
#ifdef USE_DSYRK
		dsymv_("U", &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);

		elementWiseMult(aGout, mInvCodonFreq);

#elif defined(USE_DGEMM)
		dgemv_("N", &N, &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);
#else
		dgemv_("T", &N, &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);
#endif
#else
		for(int r=0; r < N; ++r)
		{
			double x = 0;
			for(int c=0; c < N; ++c) x += mMatrices[aSetIdx*mNumMatrices+aBranch][r*N+c]*aGin[c];
			aGout[r] = x;
		}
#endif
	}
#else

	///	Multiply the aMin fat-vector by the precomputed exp(Q*t) matrix
	///
	/// @param[in] aSetIdx Which set to use (starts from zero)
	/// @param[in] aBranch Which branch
	/// @param[in] aNumSites Number of sites composing the fat-vector
	/// @param[in] aMin The input fat-vector to be multiplied by the matrix exponential
	/// @param[out] aMout The resulting fat-vector
	///
	void doTransition(unsigned int aSetIdx, unsigned int aBranch, int aNumSites, const double* aMin, double* aMout) const
	{
#ifdef USE_LAPACK
#ifdef USE_DSYRK
	
	dsymm_("L", "U", &N, &aNumSites, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aMin, &VECTOR_SLOT, &D0, aMout, &VECTOR_SLOT);

#if 0
	for(int c=0; c < aNumSites; ++c)
	{
		for(int r=0; r < N; ++r)
		{
			aMout[c*VECTOR_SLOT+r] *= mInvCodonFreq[r];
		}
	}
#endif

#ifdef USE_MKL_VML
	for(int c=0; c < aNumSites; ++c)
	{
		vdMul(N, &aMout[c*VECTOR_SLOT], mInvCodonFreq, &aMout[c*VECTOR_SLOT]);
	}
#else
	for(int r=0; r < N; ++r)
	{
		dscal_(&aNumSites, &mInvCodonFreq[r], aMout+r, &VECTOR_SLOT);
	}
#endif


#elif defined(USE_DGEMM)
		dgemm_( "N",
				"N",
				&N,
				&aNumSites,
				&N,
				&D1,
				mMatrices[aSetIdx*mNumMatrices+aBranch],
				&N,
				aMin,
				&VECTOR_SLOT,
				&D0,
				aMout,
				&VECTOR_SLOT);
#else
		dgemm_( "T",
				"N",
				&N,
				&aNumSites,
				&N,
				&D1,
				mMatrices[aSetIdx*mNumMatrices+aBranch],
				&N,
				aMin,
				&VECTOR_SLOT,
				&D0,
				aMout,
				&VECTOR_SLOT);
#endif
#else
		for(int r=0; r < N; ++r)
		{
			for(int c=0; c < aNumSites; ++c)
			{
				double x = 0;
				for(int k=0; k < N; ++k) x += mMatrices[aSetIdx*mNumMatrices+aBranch][r*N+k]*aMin[c*VECTOR_SLOT+k]; // aMin is transposed
				aMout[c*VECTOR_SLOT+r] = x; // also aMout is transposed
			}
		}
#endif
	}
#endif

	/// Return the number of sets contained.
	///
	/// @return The number of sets contained
	///
	unsigned int size(void) const {return mNumSets;}

private:
	unsigned int	mNumMatrices;		///< Number of matrices in each set
	unsigned int	mNumSets;			///< Number of sets
	double*			mMatrixSpace;		///< Starts of the matrix storage area
	double**		mMatrices;			///< Access to the matrix set
	const double*	mInvCodonFreq;		///< Inverse of the codon frequencies
};


#endif

