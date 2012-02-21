
#include <iostream>
#include <iomanip>
#include "BayesTest.h"


// The minimum value for class 2 sites probability to be a positive selection site.
static const double MIN_PROB = 0.95;

BayesTest::BayesTest(size_t aNumSites) : mNumSites(aNumSites)
{
	mSiteClassProb.resize(4*mNumSites);
}


void BayesTest::computeBEB(void)
{
	// Test
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		mSiteClassProb[0*mNumSites+site] = randFrom0to1();
		mSiteClassProb[1*mNumSites+site] = randFrom0to1();
		mSiteClassProb[2*mNumSites+site] = randFrom0to1();
		mSiteClassProb[3*mNumSites+site] = randFrom0to1();

		double tot = mSiteClassProb[0*mNumSites+site] +
					 mSiteClassProb[1*mNumSites+site] +
					 mSiteClassProb[2*mNumSites+site] +
					 mSiteClassProb[3*mNumSites+site];

		mSiteClassProb[0*mNumSites+site] /= tot;
		mSiteClassProb[1*mNumSites+site] /= tot;
		mSiteClassProb[2*mNumSites+site] /= tot;
		mSiteClassProb[3*mNumSites+site] /= tot;

		// Just to be sure at least some site appears
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > 0.9 && prob < 0.95)
		{
			mSiteClassProb[2*mNumSites+site] += 0.025;
			mSiteClassProb[3*mNumSites+site] += 0.025;
			mSiteClassProb[0*mNumSites+site] -= 0.025;
			mSiteClassProb[1*mNumSites+site] -= 0.025;
		}
	}
}

void BayesTest::printPositiveSelSites(size_t aFgBranch) const
{
	bool print_title = true;

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if is a type 2 site with prob > 95%
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > MIN_PROB)
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


void BayesTest::extractPositiveSelSites(std::vector<unsigned int>& aPositiveSelSites, std::vector<double>& aPositiveSelSitesProb) const
{
	// Prepare the output vectors
	aPositiveSelSites.clear();
	aPositiveSelSitesProb.clear();

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if it is a type 2 site with prob > 95%
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > MIN_PROB)
		{
			aPositiveSelSites.push_back(site);
			aPositiveSelSitesProb.push_back(prob);
		}
	}
}
