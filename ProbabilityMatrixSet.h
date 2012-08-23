
#ifndef PROBABILITYMATRIXSET_H
#define PROBABILITYMATRIXSET_H

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
#ifdef SAVE_OCTAVE
extern void SaveToOctave(const double *CVariable, char *OctaveVariable, FILE *FilePointer, int Rows, int Columns); //CHHS
static int x = 1;
#endif

// Uncomment to move the element wise multiplication up one level
#define BUNDLE_ELEMENT_WISE_MULT

//extern int num_matvect;
/// Set of probability matrices for all branches of a tree.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-04-05 (initial version)
///     @version 1.0
///
class ProbabilityMatrixSet
{
public:

	/// Create matrix set
	///
	/// @param[in] aNumMatrices The number of matrices to be managed (is the number of branches of the tree)
	/// @param[in] aNumSets How many sets to allocate (one set is composed by the bg and fg matrices for one of the tree traversals)
	///
	ProbabilityMatrixSet(size_t aNumMatrices, unsigned int aNumSets) : mNumMatrices(static_cast<int>(aNumMatrices)), mNumSets(aNumSets), mFgBranch(0)
	{
		mMatrixSpace  = static_cast<double*>(alignedMalloc(sizeof(double)*aNumSets*aNumMatrices*MATRIX_SLOT, CACHE_LINE_ALIGN));
		mMatrices     = static_cast<double**>(alignedMalloc(sizeof(double*)*aNumSets*aNumMatrices, CACHE_LINE_ALIGN));
		mInvCodonFreq = CodonFrequencies::getInstance()->getInvCodonFrequencies();
		mPrevTime     = new double[aNumMatrices];
	}

	/// Destructor.
	///
	~ProbabilityMatrixSet()
	{
		alignedFree(mMatrixSpace);
		alignedFree(mMatrices);
		delete [] mPrevTime;
	}

	/// Return the number of sets contained in this ProbabilityMatrixSet
	///
	/// @return The number of sets
	///
	unsigned int size(void) const {return mNumSets;}

	/// Initialize the set for a given foreground branch number for H0
	///
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	///
	void initForH0(unsigned int aFgBranch);

	/// Initialize the set for a given foreground branch number for H1
	///
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	///
	void initForH1(unsigned int aFgBranch);

	/// Compute the three sets of matrices for the H0 hypothesis.
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w1
	///
	///	@param[in] aQw0 The mQw0 transition matrix
	///	@param[in] aQ1 The mQ1 transition matrix
	///	@param[in] aAnyMatrixChanged Set to true if any matrix changed
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (k, w0, 02, p0+p1, p0/(p0+p1), w2)
	///
	void computeMatrixSetH0(const TransitionMatrix& aQw0,
						    const TransitionMatrix& aQ1,
							bool  aAnyMatrixChanged,
							double aSbg,
							double aSfg,
						    const std::vector<double>& aParams);

	/// Compute the four sets of matrices for the H1 hypothesis.
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w2
	/// - set 3: w1, w2
	///
	///	@param[in] aQw0 The mQw0 transition matrix
	///	@param[in] aQ1 The mQ1 transition matrix
	///	@param[in] aQw2 The mQw2 transition matrix
	///	@param[in] aChangedQw2 Set to true if Qw2 matrix changed
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (k, w0, 02, p0+p1, p0/(p0+p1), w2)
	///
	void computeMatrixSetH1(const  TransitionMatrix& aQw0,
						    const  TransitionMatrix& aQ1,
						    const  TransitionMatrix& aQw2,
						    bool   aChangedQw2,
							double aSbg,
							double aSfg,
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
		//num_matvect++;
#ifdef USE_LAPACK
		dsymv_("U", &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);

#if !defined(BUNDLE_ELEMENT_WISE_MULT)
	elementWiseMult(aGout, mInvCodonFreq);
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

	///	Multiply the aGin vector by the precomputed exp(Q*t) matrix
	///
	/// @param[in] aSetIdx Which set to use (starts from zero)
	/// @param[in] aBranch Which branch
	/// @param[in] aCodon The leaf codon id (will supplant aGin)
	/// @param[out] aGout The resulting vector
	///
	void doTransitionAtLeaf(unsigned int aSetIdx, unsigned int aBranch, int aCodon, double* aGout) const
	{

#ifdef SAVE_OCTAVE
		if(x) {
			FILE* fp = fopen("sasa.oct", "w");
			SaveToOctave(mMatrices[aSetIdx*mNumMatrices+aBranch], "M", fp, N, N); //CHHS
			fclose(fp);
			x = 0;
			int i;
			int k = 5;
			for(i=0; i < k; ++i) printf("%10.7lf\n", mMatrices[aSetIdx*mNumMatrices+aBranch][k*N+i]);
			for(i=k; i < N; ++i) printf("%10.7lf\n", mMatrices[aSetIdx*mNumMatrices+aBranch][i*N+k]);
		}
#endif

#ifdef USE_LAPACK
		// Simply copy the symmetric matrix column
		// instead of: dsymv_("U", &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);
		int i = 0;
		for(; i < aCodon; ++i) aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][aCodon*N+i];
		for(; i < N; ++i)      aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][i*N+aCodon];

#if !defined(BUNDLE_ELEMENT_WISE_MULT)
		elementWiseMult(aGout, mInvCodonFreq);
#endif
#else
		for(int r=0; r < N; ++r)
		{
			aGout[r] = mMatrices[aSetIdx*mNumMatrices+aBranch][r*N+aCodon];
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
	
	dsymm_("L", "U", &N, &aNumSites, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aMin, &VECTOR_SLOT, &D0, aMout, &VECTOR_SLOT);

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

private:
	int				mNumMatrices;		///< Number of matrices in each set (should be int)
	unsigned int	mNumSets;			///< Number of sets
	int				mFgBranch;			///< Foreground branch number (should be int)
	double*			mMatrixSpace;		///< Starts of the matrix storage area
	double**		mMatrices;			///< Access to the matrix set
	const double*	mInvCodonFreq;		///< Inverse of the codon frequencies
	double*			mPrevTime;			///< Previous computation times (to check if they are sufficiently different to require a recomputation)
};

#endif

