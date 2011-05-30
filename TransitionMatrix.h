
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

	/// Fill the Q matrix and return the matrix scale value.
	///
	/// @param[in] aOmega The omega value.
	/// @param[in] aK The k value.
	/// @param[in] aCodonFreq The codon frequency array
	///
	/// @return The Q matrix scale value.
	///
	double fillQ(double aOmega, double aK, const double* aCodonFreq);

	/// Fill the Q matrix and return the matrix scale value. Optimized routine to be used for omega == 1
	///
	/// @param[in] aK The k value.
	/// @param[in] aCodonFreq The codon frequency array
	///
	/// @return The Q matrix scale value.
	///
	double fillQ(double aK, const double* aCodonFreq);

	/// Compute the eigendecomposition of the Q matrix.
	/// The results are stored internally
	///
	/// @param[in] aNumGoodFreq Number of not vanishing codon frequencies
	/// @param[in] aSqrtCodonFreq Square root of the frequency of the corresponding codon
	/// @param[in] aGoodFreq Flag to mark the corresponding frequncy "not small"
	///
	void eigenQREV(int aNumGoodFreq, const double* aSqrtCodonFreq, const bool* aGoodFreq);

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

	/// Store in an external matrix the result of exp(Q*t)
	///
	/// @param[out] aOut The matrix where the result should be stored (size: N*N) under USE_LAPACK it is stored transposed
	/// @param[in] aT The time to use in the computation (it is always > 0)
	/// @param[in,out] aWorkarea Temporary array N*N that should be zeroed before first usage
	///
#if defined(USE_LAPACK) && defined(USE_DGEMM)
	inline void computeFullTransitionMatrix(double* aOut, double aT, double* aWorkarea) const
	{
		for(int i=0; i < N; ++i) aWorkarea[i*(N+1)] = exp(aT * mD[i]);
		double tmp[N*N];
		dgemm_("N", "T", &N, &N, &N, &D1, aWorkarea, &N,  mV, &N, &D0,  tmp, &N);
		dgemm_("T", "N", &N, &N, &N, &D1, mU,        &N, tmp, &N, &D0, aOut, &N);
#else
	inline void computeFullTransitionMatrix(double* aOut, double aT, double* /*aWorkarea*/) const
	{
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
};

#endif


