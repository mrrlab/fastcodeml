
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
	/// @param[out] aOut The matrix where the result should be stored (size: N*N)
	/// @param[in] aT The time to use in the computation (it is always > 0)
	///
	inline void computeFullTransitionMatrix(double* aOut, double aT) const
	{
#if 0
		memset(aOut, 0, N*N*sizeof(double));

		for(int k = 0; k < N; ++k)
		{
			double *p = aOut;
			double expt = exp(aT * mD[k]);

			for(int i=0; i < N; ++i)
			{
				double uexpt = mU[i*N + k] * expt;

				for(int j=0; j < N; ++j)
				{
					*p++ += uexpt * mV[k*N + j];
				}
			}
		}
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

#if 0
	/// Compute matrix exponential from the decomposed matrix and the given time.
	///
	/// @param[in] aT The time value
	/// @param[in] aGin The input vector to be multiplied by the matrix exponential
	/// @param[out] aGout The resulting vector
	///
	inline void expMat(double aT, const double* aGin, double* aGout) const
	{
		if(aT > NEAR_ZERO_TIME || aT < -NEAR_ZERO_TIME)
		{
			// M = U * diag(exp(t*D)) * V
			// gout = M * gin
			// gout = U * diag(exp(t*D)) * (V * gin)
#ifdef USE_LAPACK
			double tmp[N];
			dgemv_("T", &N, &N, &D1, mV, &N, aGin, &I1, &D0, tmp, &I1);
			for(int k = 0; k < N; ++k) tmp[k] *= exp(aT * mD[k]);
			dgemv_("T", &N, &N, &D1, mU, &N, tmp, &I1, &D0, aGout, &I1);
#else
#ifdef OLD_METHOD
			double P[N*N];
			memset(P, 0, N*N*sizeof(double));

			unsigned int k, i, j;
			double expt, uexpt, *pP;

			for(k = 0; k < N; ++k)
			{
				for(i=0, pP=P, expt=exp(aT * mD[k]); i < N; ++i)
				{
					for(j=0, uexpt=mU[i*N + k] * expt; j < N; ++j)
					{
						*pP++ += uexpt * mV[k*N + j];
					}
				}
			}

			for(i=0; i < N; ++i)
			{
				double x = 0;
				for(j=0; j < N; ++j) x += P[i*N+j]*aGin[j];
				aGout[i] = x;
			}
#else
			double tmp[N];
			unsigned int i, j;

			for(i=0; i < N; ++i)
			{
				double x = 0;
				for(j=0; j < N; ++j) x += mV[i*N+j]*aGin[j];
				tmp[i] = x * exp(aT * mD[i]);
			}
			for(i=0; i < N; ++i)
			{
				double x = 0;
				for(j=0; j < N; ++j) x += mU[i*N+j]*tmp[j];
				aGout[i] = x;
			}
#endif
#endif
		}
		else
		{
			// Time near zero, no transform of the vector
			memcpy(aGout, aGin, N*sizeof(double));
		}
	}
#endif

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
	double mD[N];		///< The diagonalized matrix
};

#endif


