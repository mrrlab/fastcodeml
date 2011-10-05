
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "BayesTest.h"

BayesTest::BayesTest(size_t aNumSites) : mNumSites(aNumSites)
{
	mSiteClassProb.resize(4*mNumSites);
}


void BayesTest::computeBEB(void)
{
	// Test
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		mSiteClassProb[0*mNumSites+site] = rand()/(double)RAND_MAX;
		mSiteClassProb[1*mNumSites+site] = rand()/(double)RAND_MAX;
		mSiteClassProb[2*mNumSites+site] = rand()/(double)RAND_MAX;
		mSiteClassProb[3*mNumSites+site] = rand()/(double)RAND_MAX;

		double tot = mSiteClassProb[0*mNumSites+site] +
					 mSiteClassProb[1*mNumSites+site] +
					 mSiteClassProb[2*mNumSites+site] +
					 mSiteClassProb[3*mNumSites+site];

		mSiteClassProb[0*mNumSites+site] /= tot;
		mSiteClassProb[1*mNumSites+site] /= tot;
		mSiteClassProb[2*mNumSites+site] /= tot;
		mSiteClassProb[3*mNumSites+site] /= tot;
	}
}

void BayesTest::computeNEB(void)
{
}

void BayesTest::printPositiveSelSites(unsigned int aFgBranch)
{
	bool print_title = true;

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if is a type 2 site with prob > 95%
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > 0.95)
		{
			// Put a title the firts time
			if(print_title)
			{
				std::cerr << "Printing positive sel sites for branch " << aFgBranch << std::endl;
				print_title = false;
			}
			
			// Print site number and probability
			std::cerr << std::setw(6) << site << " " << std::fixed << std::setprecision(5) << prob << std::endl;
		}
	}
}



