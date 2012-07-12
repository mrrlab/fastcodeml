
#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

#include <iostream>

#include <cstring>
#include <cmath>
#include <vector>
#include <bitset>
#include "MatrixSize.h"
#include "CompilerHints.h"
#include "CodonFrequencies.h"

#ifdef USE_MKL_VML
#include <mkl_vml_functions.h>
#endif

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
	TransitionMatrix()
	{
		// Initialize Q matrix to all zeroes (so only non-zero values are written)
#ifdef USE_S_MATRIX
		memset(mS, 0, N*N*sizeof(double));
#else
		memset(mQ, 0, N*N*sizeof(double));
#endif		
		// Initialize the codons' frequencies
		CodonFrequencies* cf = CodonFrequencies::getInstance();
		mCodonFreq = cf->getCodonFrequencies();
		mNumGoodFreq = cf->getNumGoodCodons();
		mSqrtCodonFreq = cf->getSqrtCodonFrequencies();
		cf->cloneGoodCodonIndicators(mGoodFreq);

#ifndef USE_LAPACK
		// Fill the identity matrix (to be used when time is zero)
		memset(mIdentity, 0, N*N*sizeof(double));
		for(int i=0; i < N; ++i) mIdentity[i*(1+N)] = 1.0;
#endif
	}

	/// Fill the Q (or the S) matrix and return the matrix scale value.
	///
	/// @param[in] aOmega The omega value.
	/// @param[in] aK The k value.
	///
	/// @return The Q matrix scale value.
	///
	double fillMatrix(double aOmega, double aK);

	/// Fill the Q (or the S) matrix and return the matrix scale value. Optimized routine to be used for omega == 1
	///
	/// @param[in] aK The k value.
	///
	/// @return The Q matrix scale value.
	///
	double fillMatrix(double aK);

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

	/// Check the eigen decomposition on the reduced matrix.
	///
	/// @param[in] aDim Dimension of the reduced matrix
	/// @param[in] aPrev The reduced matrix before the eigendecomposition
	/// @param[in] aFull If true also partially print the matrices
	///
	void checkReducedEigen(int aDim, const double* aPrev, bool aFull=false) const;

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
	void computeFullTransitionMatrix(double* RESTRICT aOut, double aT) const
	{
std::cerr << "(*)P" << std::endl;
#ifdef USE_LAPACK

		double ALIGN64 tmp[N*N64]; // The rows are padded to 64 to increase performance
		double ALIGN64 expt[N];

		double tm = aT / 2.;
#ifndef USE_MKL_VML
		// Manual unrolling gives the best results here.
		// So it is exp(D*T/2). Remember, the eigenvalues are stored in reverse order
		for(int c=0; c < N-1; )
		{
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
		}
		expt[N-1] = exp(tm*mD[0]);

		for(int r=0; r < N; ++r)
		{
			for(int c=0; c < N; ++c)
			{
				tmp[r*N64+c] = expt[c]*mV[r*N+c];
				//tmp[r*N+c] = expt[c]*mV[r*N+c];
			}
		}
#else
		// Manual unrolling gives the best results here.
		// Remember, the eigenvalues are stored in reverse order
        for(int j=0; j < N-1; )
        {
            tmp[j] = tm*mD[N-1-j]; ++j; 
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
        }
		tmp[60] = tm*mD[0];

		vdExp(N, tmp, expt);
		for(int r=0; r < N; ++r)
		{
			vdMul(N, expt, &mV[r*N], &tmp[r*N64]);
			//vdMul(N, expt, &mV[r*N], &tmp[r*N]);
		}
#endif

		dsyrk_("U", "T", &N, &N, &D1, tmp, &N64, &D0, aOut, &N);
		//dsyrk_("U", "T", &N, &N, &D1, tmp, &N, &D0, aOut, &N);

#else
		// if time is zero or almost zero, the transition matrix become an identity matrix
		if(aT < 1e-100)
		{
			memcpy(aOut, mIdentity, N*N*sizeof(double));
			return;
		}

		// The first iteration of the loop (k == 0) is split out to initialize aOut
		double *p = aOut;
		double expt = exp(aT * mD[N-1]); // Remember, the eigenvalues are stored in reverse order
		for(int i=0; i < N; ++i)
		{
			const double uexpt = mU[i*N] * expt;

			for(int j=0; j < N; ++j)
			{
				*p++ = uexpt * mV[j];
			}
		}

		// The subsequent iterations are computed normally
		for(int k = 1; k < N; ++k)
		{
			p = aOut;
			expt = exp(aT * mD[N-1-k]); // Remember, the eigenvalues are stored in reverse order

			for(int i=0; i < N; ++i)
			{
				const double uexpt = mU[i*N + k] * expt;

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
	void eigenRealSymm(double* RESTRICT aU, int aDim, double* RESTRICT aR, double* RESTRICT aWork);


#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4324) // Padding added due to alignment request
#endif

private:
	// Order suggested by icc to improve locality
	// 'mV, mU, mSqrtCodonFreq, mNumGoodFreq, mQ, mD, mCodonFreq, mGoodFreq'
	double ALIGN64	mV[N*N];		///< The right adjusted eigenvectors matrix (with the new method instead contains pi^1/2*R where R are the autovectors)
	double ALIGN64	mU[N*N];		///< The left adjusted eigenvectors matrix
#ifndef USE_LAPACK
	double ALIGN64	mIdentity[N*N];	///< Pre-filled identify matix
#endif
	const double*	mSqrtCodonFreq;	///< Square root of experimental codon frequencies
	int				mNumGoodFreq;	///< Number of codons whose frequency is not zero
#ifndef USE_S_MATRIX
	double ALIGN64	mQ[N*N];		///< The Q matrix
#else
	double ALIGN64	mS[N*N];		///< The S matrix (remember Q = S*Pi)
#endif
	double ALIGN64	mD[N];			///< The matrix eigenvalues stored in reverse order
	const double*	mCodonFreq;		///< Experimental codon frequencies
	std::bitset<N>	mGoodFreq;		///< True if the corresponding codon frequency is not small
	
#ifdef _MSC_VER
    #pragma warning(pop)
#endif
};

#endif


