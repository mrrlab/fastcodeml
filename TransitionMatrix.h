
#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

#include <cstring>
#include <cmath>
#include <vector>
#include "MatrixSize.h"
#ifdef USE_MKL_VML
#include <mkl_vml_functions.h>
#endif

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
	TransitionMatrix() : mDim(N), mGoodFreq(N, true)
	{
		memset(mQ, 0, N*N*sizeof(double));
		
		// Just to be sure
		mCodonFreq = 0;
		mNumGoodFreq = 0;
		mSqrtCodonFreq = 0;
		//mGoodFreq = 0;
	}

	/// Store the precomputed codon frequency array, its square root and indication of the non null values.
	///
	/// @param[in] aCodonFreq The codon frequency array
	/// @param[in] aNumGoodFreq Number of not vanishing codon frequencies
	/// @param[in] aSqrtCodonFreq Square root of the frequency of the corresponding codon
	/// @param[in] aGoodFreq Flag to mark the corresponding frequncy "not small"
	///
	//inline void setCodonFrequencies(const double* aCodonFreq, unsigned int aNumGoodFreq, const double* aSqrtCodonFreq, const bool* aGoodFreq)
	inline void setCodonFrequencies(const double* aCodonFreq, unsigned int aNumGoodFreq, const double* aSqrtCodonFreq, const std::vector<bool>& aGoodFreq)

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
	/// Depending on the definition of USE_DSYRK use the old or the new method
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
#if defined(USE_LAPACK) && defined(USE_DGEMM) && !defined(USE_DSYRK)

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

#elif defined(USE_LAPACK) && defined(USE_DSYRK)

		double tmp[N*N];
		double expt[N];

		aT /= 2.;
#ifndef USE_MKL_VML
		for(int c=0; c < N; ++c)
		{
			expt[c] = exp(aT*mD[c]); // So it is exp(D*T/2)
		}
#else
		vdLinearFrac(N, mD, mD, aT, 0.0, 0.0, 1.0, tmp);
		vdExp(N, tmp, expt);
#endif

#if 0
		memcpy(tmp, mV, N*N*sizeof(double));
#endif
		for(int c=0; c < N; ++c)
		{
#if 1
			for(int r=0; r < N; ++r)
			{
				tmp[r*N+c] = expt[c]*mV[r*N+c];
			}
#else
			dscal_(&N, &expt[c], tmp+c, &N);
#endif
		}
		dsyrk_("U", "T", &N, &N, &D1, tmp, &N, &D0, aOut, &N);

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
	/// Order suggested by icc to improve locality
	/// 'mV, mCodonFreq, mQ, mDim, mSqrtCodonFreq, mD, mNumGoodFreq, mU, mGoodFreq'
	double			mV[N*N];		///< The right adjusted eigenvectors matrix (with the new method instead contains pi^1/2*R where R are the autovectors)
	const double*	mCodonFreq;		///< Experimental codon frequencies
	double			mQ[N*N];		///< The Q matrix
	int				mDim;			///< The matrix size (should be <= N)
	const double*	mSqrtCodonFreq;	///< Square Root of experimental codon frequencies
	double			mD[N];			///< The matrix eigenvalues
	//const bool*		mGoodFreq;	///< True if the corresponding codon frequency is not small
	int				mNumGoodFreq;	///< Number of codons whose frequency is not zero
	double			mU[N*N];		///< The left adjusted eigenvectors matrix
	std::vector<bool>
					mGoodFreq;		///< True if the corresponding codon frequency is not small
};

#endif


