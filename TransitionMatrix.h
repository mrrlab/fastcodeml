
#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

#include <cstring>
#include <cmath>
#include "MatrixSize.h"

/// If the time is in absolute value less than this, consider it zero
///
static const double NEAR_ZERO_TIME = 1e-100;

#ifdef USE_LAPACK
#include "blas.h"
#endif


/// The transition matrix plus its eigen decomposition.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-02-23 (initial version)
///     @version 1.0
///
class TransitionMatrix
{
public:
	/// Constructor
	///
	/// @param[in] aDim The matrix size (should be <= N)
	///
	TransitionMatrix(int aDim=N)
	{
		mDim = aDim;
		memset(mQ, 0, N*N*sizeof(double));
	}

	/// Store the precomputed codon frequency array, its square root and indication of the non null values.
	///
	/// @param[in] aCodonFreq The codon frequency array
	/// @param[in] aNumGoodFreq Number of not vanishing codon frequencies
	/// @param[in] aSqrtCodonFreq Square root of the frequency of the corresponding codon
	/// @param[in] aGoodFreq Flag to mark the corresponding frequncy "not small"
	///
	void setCodonFrequencies(const double* aCodonFreq, unsigned int aNumGoodFreq, const double* aSqrtCodonFreq, const bool* aGoodFreq)
	{
		mCodonFreq      = aCodonFreq;
		mNumGoodFreq	= aNumGoodFreq;
		mSqrtCodonFreq	= aSqrtCodonFreq;
		mGoodFreq		= aGoodFreq;
	}

	/// Fill the Q matrix and return the matrix scale value.
	///
	/// @param[in] aOmega The omega value.
	/// @param[in] aK The k value.
	/// @param[in] aCodonFreq The codon frequency array
	///
	/// @return The Q matrix scale value.
	///
	double fillQ(double aOmega, double aK);

	/// Fill the Q matrix and return the matrix scale value. Optimized routine to be used for omega == 1
	///
	/// @param[in] aK The k value.
	/// @param[in] aCodonFreq The codon frequency array
	///
	/// @return The Q matrix scale value.
	///
	double fillQ(double aK);

	/// Compute the eigendecomposition of the Q matrix.
	/// The used codon frequencies should be already loaded using setCodonFrequencies()
	/// The results are stored internally
	///
	void eigenQREV(void);

#ifdef CHECK_ALGO
	/// Check the eigen decomposition.
	///
	/// @param[in] aFull If true also partially print the matrices
	///
	void checkEigen(bool aFull=false) const;

	/// Print the Q matrix
	///
	/// @param[in] aMaxRow Max number of rows to print
	/// @param[in] aMaxCol Max number of columns to print. If missing or equal zero, then takes the same value as aMaxRow
	///
	void print(unsigned int aMaxRow=6, unsigned int aMaxCol=0) const;

	/// Print all the decomposed matrices (U, V, D)
	///
	/// @param[in] aMaxRow Max number of rows to print
	/// @param[in] aMaxCol Max number of columns to print. If missing or equal zero, then takes the same value as aMaxRow
	///
	void printDecomposed(unsigned int aMaxRow=6, unsigned int aMaxCol=0) const;
#endif
	/// Store in an external matrix the result of exp(Q*t)
	///
	/// @param[out] aOut The matrix where the result should be stored (size: N*N) under USE_LAPACK it is stored transposed
	/// @param[in] aT The time to use in the computation (it is always > 0)
	///
	inline void computeFullTransitionMatrix(double* aOut, double aT) const
	{
#if defined(USE_LAPACK) && defined(USE_DGEMM)

		double tmp[N*N];
		memcpy(tmp, mV, sizeof(double)*N*N);
		//CHHS Compute D*V by multiplying row 1 by e^(first root), row 2 by e^(second root) etc.
		int i, j;
		for(i=j=0; i < N; ++i, j+=N)
		{
			//CHHS DSCAL(N,DA,DX,INCX); y = alpha * y
			double expt = exp(aT*mD[i]);
			dscal_(&N, &expt, tmp+j, &I1); //CHHS Scale single row    
		}
		dgemm_("T", "T", &N, &N, &N, &D1, mU, &N, tmp, &N, &D0, aOut, &N);

#else

		// The first iteration of the loop (k == 0) is split out to initialize aOut
		double *p = aOut;
		double expt = exp(aT * mD[0]);

		for(int i=0; i < N; ++i)
		{
			double uexpt = mU[i*N] * expt;

			for(int j=0; j < N; ++j)
			{
				*p++ = uexpt * mV[j];
			}
		}

		// The subsequent iterations are computed normally
		for(int k = 1; k < N; ++k)
		{
			p = aOut;
			expt = exp(aT * mD[k]);

			for(int i=0; i < N; ++i)
			{
				double uexpt = mU[i*N + k] * expt;

				for(int j=0; j < N; ++j)
				{
					*p++ += uexpt * mV[k*N + j];
				}
			}
		}
#endif
	}

private:
	/// Compute the eigendecomposition
	///
	/// @param[in,out] aU Matrix to be decomposed on input, Eigenvectors on output
	/// @param[in] aDim The matrix dimension
	/// @param[out] aR The eigenvalues
	/// @param[out] aWork A working area used only for non lapack version
	///
	void inline eigenRealSymm(double* aU, int aDim, double* aR, double* aWork);


private:
	int    mDim;		///< The matrix size (should be <= N)
	double mQ[N*N];		///< The Q matrix
	double mU[N*N];		///< The left adjusted eigenvectors matrix
	double mV[N*N];		///< The right adjusted eigenvectors matrix
	double mD[N];		///< The matrix eigenvalues

	const double*	mCodonFreq;		///< Experimental codon frequencies
	int				mNumGoodFreq;	///< Number of codons whose frequency is not zero
	const double*	mSqrtCodonFreq;	///< Square Root of experimental codon frequencies
	const bool*		mGoodFreq;		///< True if the corresponding codon frequency is not small
};

#endif


