
#include <iostream>
#include <iomanip>
#include "BayesTest.h"
#include "VerbosityLevels.h"


/// The minimum value for class 2 sites probability to be a positive selection site.
///
static const double MIN_PROB = 0.95;

BayesTest::BayesTest(size_t aNumSites, unsigned int aVerbose) : mSiteClassProb(4*aNumSites), mNumSites(aNumSites), mVerbose(aVerbose), mPriors(aNumSites*BEB_NUM_CAT)
{
	// Searching range for w0
	const static double w0b0 = 0.;
	const static double w0b1 = 1.;

	// Searching range for w2
	const static double w2b0 = 1.;
	const static double w2b1 = 11.;

	// Initialize the w0 and w2 values to be tested
	for(int i=0; i < BEB_N1D; i++)
	{
		mPara[0][i] = mPara[1][i] = -1;						// p0 & p1
		mPara[2][i] = w0b0 + (i+0.5)*(w0b1-w0b0)/BEB_N1D;	// w0
		mPara[3][i] = w2b0 + (i+0.5)*(w2b1-w2b0)/BEB_N1D;	// w2
	}
}


double BayesTest::getGridParams(BranchSiteModelAltHyp& aModel)
{
	// Omega for foreground and background branches
	double omega_fg;
	double omega_bg;

	BranchSiteModelAltHyp::CodonClass codon_class;

	// Calculating f(x_h|w) 
	// Order of site classes for iw or f(x_h|w):
	//                     fore   back     #sets
	//   site class 0:      w0     w0        10
	//   site class 1:      w1=1   w1=1       1
	//   site class 2a:     w0     w2       100
	//   site class 2b:     w1=1   w2        10

	for(int iw=0; iw < BEB_NUM_CAT; ++iw)
	{
		if(iw < BEB_N1D)										// class 0: w0 w0
		{
			omega_fg = omega_bg = mPara[2][iw];
			codon_class = BranchSiteModelAltHyp::CODON_CLASS_0;
		}
		else if(iw == BEB_N1D)									// class 1: w1 w1
		{
			omega_fg = omega_bg = 1.;
			codon_class = BranchSiteModelAltHyp::CODON_CLASS_1;
		}
		else if(iw < (BEB_N1D+1+BEB_N1D*BEB_N1D))				// class 2a: w0 w2
		{                                  
            omega_fg = mPara[2][(iw-BEB_N1D-1) / BEB_N1D];
            omega_bg = mPara[3][(iw-BEB_N1D-1) % BEB_N1D];
			codon_class = BranchSiteModelAltHyp::CODON_CLASS_2a;
		}
		else													// class 2b: w1 w2
		{                                                       
            omega_fg = 1.;
            omega_bg = mPara[3][iw-BEB_N1D-1-BEB_N1D*BEB_N1D];
			codon_class = BranchSiteModelAltHyp::CODON_CLASS_2b;
		}

		// Do It
		aModel.computeLikelihoodForBEB(codon_class, omega_fg, omega_bg, &mPriors[iw*mNumSites]);

		std::cerr << std::setw(3) << iw << ' ';
		std::cerr << std::fixed << std::setw(8) << std::setprecision(4) << omega_fg << ' ';
		std::cerr << std::fixed << std::setw(8) << std::setprecision(4) << omega_bg << " = " << std::setprecision(16);
		std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+0] << ' ';        // Print only the first three sites
		std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+1] << ' ';
		std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+2] << std::endl;
	}

	const std::vector<double>& site_mult = aModel.getSiteMultiplicity();
	double scale = 0.;
	for(int site=0; site < mNumSites; ++site)
	{
		int k;
		double fh = mPriors[site];
		for(k=1; k < BEB_NUM_CAT; ++k)
		{
			if(mPriors[k*mNumSites+site] > fh) fh = mPriors[k*mNumSites+site];
		}
		for(k=0; k < BEB_NUM_CAT; ++k)
		{
			mPriors[k*mNumSites+site] = exp(mPriors[k*mNumSites+site]-fh);
		}
		scale += fh*site_mult[site];
	}

	return scale;
}

void BayesTest::computeBEB(BranchSiteModelAltHyp& aModel)
{
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << std::endl << "Computing BEB" << std::endl;

	// Enable it as soon as the dependency list problem as been solved.
	//	double s1 = getGridParams(aModel);

	////////////////////////////////////////////////////////////////////////////////////////
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
			std::cerr << std::setw(6) << site << " " << std::fixed << std::setprecision(6) << prob << std::endl;
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
