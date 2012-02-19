
#ifndef CODONFREQUENCIES_H
#define CODONFREQUENCIES_H

#include <vector>
#include <bitset>
#include <cmath>
#include <cfloat>
#include "MatrixSize.h"
#include "Types.h"

/// If codon probability is greater than this value, the codon is marked as "good codon".
///
static const double GOOD_CODON_THRESHOLD = 1e-100;

/// Compute and distribute codon empirical frequencies.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-10-24 (initial version)
///     @version 1.0
///
class CodonFrequencies
{
public:
	/// Return a pointer to the singleton instance
	///
	/// @return The pointer to the instance
	///
	static CodonFrequencies* getInstance(void);

	/// Codon empirical frequencies models
	///
	enum CodonFrequencyModel
	{
		CODON_FREQ_MODEL_UNIF=0,	///< All codon probabilities are equal
		CODON_FREQ_MODEL_F3X4		///< F3x4 model
	};

	/// Compute the codon frequencies from the codon count
	///
	/// @param[in] aCodonCount Count of condons of a certain type
	/// @param[in] aModel Codon frequency model to use
	///
	void setCodonFrequencies(const std::vector<unsigned int>& aCodonCount, CodonFrequencyModel aModel)
	{
		// Compute mCodonFrequencies based on the selected model
		switch(aModel)
		{
		default:
		case CODON_FREQ_MODEL_UNIF:
			mCodonFrequencies.assign(N, 1./(double)N);
			break;

		case CODON_FREQ_MODEL_F3X4:
			setCodonFrequenciesF3x4(aCodonCount);
			break;
		}

		// Support values needed for the eigensolver
		mNumGoodCodons = 0;
		for(int k=0; k < N; ++k)
		{
			mCodonFreqSqrt[k] = sqrt(mCodonFrequencies[k]);

			// Count the number of valid codons
			if(mCodonFrequencies[k] > GOOD_CODON_THRESHOLD)
			{
				mGoodCodon.set(k);
				++mNumGoodCodons;
				mCodonFreqInv[k] = 1./mCodonFrequencies[k];
			}
			else
			{
				mGoodCodon.reset(k);
				mCodonFreqInv[k] = DBL_MAX;
			}
		}
	}

	/// Return a pointer to the codon frequencies array
	///
	/// @return Pointer to the codon frequencies array 
	///
	const double* getCodonFrequencies(void) const {return &mCodonFrequencies[0];}

	/// Return a pointer to the codon square root of frequencies array
	///
	/// @return Pointer to the codon square root of frequencies array 
	///
	const double* getSqrtCodonFrequencies(void) const {return &mCodonFreqSqrt[0];}

	/// Return a pointer to the codon inverse of frequencies array
	///
	/// @return Pointer to the codon inverse of frequencies array 
	///
	const double* getInvCodonFrequencies(void) const {return &mCodonFreqInv[0];}

	/// Return the number of non-zero codon frequencies
	///
	/// @return Number of non-zero codon frequencies
	///
	unsigned int getNumGoodCodons(void) const {return mNumGoodCodons;}

	/// Clone the array of codon not null indicators
	///
	/// @param[out] aGoodCodon Bit set indicator array
	///
	void cloneGoodCodonIndicators(std::bitset<N>& aGoodCodon) const {aGoodCodon = mGoodCodon;}

private:
	/// Set codon frequencies according to the F3x4 model
	///
	/// @param[in] aCodonCount The count of each codon occurrencies
	///
	void setCodonFrequenciesF3x4(const std::vector<unsigned int>& aCodonCount);

	/// Convert the codon number in the 1 to 64 range to 1 to 61
	///
	/// @param[in] aId64 Codon id (0 to 63)
	///
	/// @return The corresponding id in the range 0 to 60
	///
	int codon64to61(unsigned int aId64) const;

private:
	static CodonFrequencies*	mInstance;					///< Pointer to the singleton instance
	SSEAlignedDoubleVector		mCodonFrequencies;			///< Experimental codon frequencies
	SSEAlignedDoubleVector		mCodonFreqSqrt;				///< Square root of experimental codon frequencies
	SSEAlignedDoubleVector		mCodonFreqInv;				///< Inverse of experimental codon frequencies (must be aligned for SSE instructions)
	std::bitset<N>				mGoodCodon;					///< True if the corresponding codon frequency is not small
	unsigned int				mNumGoodCodons;				///< Number of codons whose frequency is not zero

protected:
	/// Protected constructor
	///
	CodonFrequencies() : mCodonFrequencies(N, 1./static_cast<double>(N)), mCodonFreqSqrt(N, 1./sqrt(static_cast<double>(N))), mCodonFreqInv(N, static_cast<double>(N)), mNumGoodCodons(N)
	{
		mGoodCodon.set();
	}
};


#endif

