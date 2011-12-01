
#include <cstring>
#include "CodonFrequencies.h"


CodonFrequencies* CodonFrequencies::mInstance = 0;

CodonFrequencies* CodonFrequencies::getInstance(void)
{
	if(!mInstance) mInstance = new CodonFrequencies();
	return mInstance;
}


int CodonFrequencies::codon64to61(unsigned int aId64) const
{
	if(aId64 > 63 || aId64 == 10 || aId64 == 11 || aId64 == 14) return -1;

	if(aId64 > 14) return aId64-3;
	if(aId64 > 11) return aId64-2;
	return aId64;
}


void CodonFrequencies::setCodonFrequenciesF3x4(const std::vector<unsigned int>& aCodonCount)
{
	int k, j;

#ifdef CHECK_ALGO
	// Print the table of codon counts
	for(k=0; k < 64; ++k)
	{
		int id = codon64to61(k);
		if(id < 0) std::cerr << std::setw(5) << 0;
		else       std::cerr << std::setw(5) << aCodonCount[id];
		if(k % 4 == 3) std::cerr << std::endl;
	}
#endif
	// Compute the 3x4 table
	double fb3x4sg[12];

	memset(fb3x4sg, 0, 12*sizeof(double));

    for(k = 0; k < 64; k++)
    {
		int kk = codon64to61(k);
		if(kk < 0) continue;

        fb3x4sg[0 * 4 + k / 16]      += aCodonCount[kk];
        fb3x4sg[1 * 4 + (k / 4) % 4] += aCodonCount[kk];
        fb3x4sg[2 * 4 + k % 4]       += aCodonCount[kk];
    }

    for(j = 0; j < 3; j++)
    {
        double t = 0;
		for(k=0; k < 4; ++k) t += fb3x4sg[j*4+k];
		for(k=0; k < 4; ++k) fb3x4sg[j*4+k] /= t;
    }

#ifdef CHECK_ALGO
	for(k=0; k < 12; ++k)
	{
		std::cerr << std::fixed << std::setprecision(6) << fb3x4sg[k];
		if(k % 4 == 3) std::cerr << std::endl;
	}
	std::cerr << std::endl;
#endif

	// Compute codon frequency from the 3x4 table
	for(k=0; k < 64; ++k)
	{
		int kk = codon64to61(k);
		if(kk < 0) continue;

		mCodonFrequencies[kk] = fb3x4sg[k / 16] * fb3x4sg[4 + (k / 4) % 4] * fb3x4sg[8 + k % 4];
	}
	double t = 0;
	for(k=0; k < N; ++k) t += mCodonFrequencies[k];
	for(k=0; k < N; ++k) mCodonFrequencies[k] /= t;

#ifdef CHECK_ALGO
	for(k=0; k < N; ++k)
	{
		std::cerr << std::fixed << std::setprecision(10) << mCodonFrequencies[k] << ' ';
		if(k % 4 == 3) std::cerr << std::endl;
	}
	std::cerr << std::endl;
#endif
}

