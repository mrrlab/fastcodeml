
#include <iostream>
#include <iomanip>
#include "BayesTest.h"
#include "VerbosityLevels.h"
#include "TreeAndSetsDependencies.h"


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
	for(unsigned int i=0; i < BEB_N1D; i++)
	{
		mPara[0][i] = mPara[1][i] = -1;						// p0 & p1
		mPara[2][i] = w0b0 + (i+0.5)*(w0b1-w0b0)/BEB_N1D;	// w0
		mPara[3][i] = w2b0 + (i+0.5)*(w2b1-w2b0)/BEB_N1D;	// w2
	}
}


double BayesTest::getGridParams(BranchSiteModelAltHyp& aModel, size_t aFgBranch)
{
	// Omega for foreground and background branches and the corresponding codon class
	double omega_fg;
	double omega_bg;

	// Get the optimized kappa
	std::vector<double> vars;
	aModel.getVariables(vars);
	size_t kappa_idx = vars.size()-2;
	const double kappa = vars[kappa_idx];

	// Compute the corresponding Q matrices for foreground and background branches
	TransitionMatrix q_fg, q_bg;
	double scale_q_fg, scale_q_bg;

	// Create the dependency list for forest likelihood computation
	Forest& forest = aModel.getForest();
	TreeAndSetsDependencies dep(forest, mVerbose);
	dep.computeDependencies(1, true);

	// Initialize the probability list
	ProbabilityMatrixSetBEB beb_set(forest.getNumBranches());
	beb_set.initializeSet(forest.adjustFgBranchIdx(aFgBranch));

	// Calculating f(x_h|w) 
	// Order of site classes for iw or f(x_h|w):
	//                     fore   back     #sets
	//   site class 0:      w0     w0        10
	//   site class 1:      w1=1   w1=1       1
	//   site class 2a:     w0     w2       100
	//   site class 2b:     w1=1   w2        10

	for(unsigned int iw=0; iw < BEB_NUM_CAT; ++iw)
	{
		if(iw < BEB_N1D)										// class 0: w0 w0
		{
			omega_fg = omega_bg = mPara[2][iw];
		}
		else if(iw == BEB_N1D)									// class 1: w1 w1
		{
			omega_fg = omega_bg = 1.;
		}
		else if(iw < (BEB_N1D+1+BEB_N1D*BEB_N1D))				// class 2a: w0 w2
		{                                  
            omega_fg = mPara[2][(iw-BEB_N1D-1) / BEB_N1D];
            omega_bg = mPara[3][(iw-BEB_N1D-1) % BEB_N1D];
		}
		else													// class 2b: w1 w2
		{                                                       
            omega_fg = 1.;
            omega_bg = mPara[3][iw-BEB_N1D-1-BEB_N1D*BEB_N1D];
		}

		// Fill the matrices and compute their eigen decomposition.
#ifdef _MSC_VER
		#pragma omp parallel sections default(none) shared(omega_fg, omega_bg, kappa, scale_q_fg, scale_q_bg, q_fg, q_bg)
#else
		#pragma omp parallel sections default(shared)
#endif
		{
			#pragma omp section
			{
				scale_q_fg = q_fg.fillMatrix(omega_fg, kappa);
				q_fg.eigenQREV();
			} 
			#pragma omp section
			{
				scale_q_bg = q_bg.fillMatrix(omega_bg, kappa);
				q_bg.eigenQREV();
			}
		}

		// Fill the matrix set with the new matrices and the times computed before
		beb_set.fillMatrixSet(q_fg, q_bg, scale_q_bg, scale_q_fg, vars);

		// Compute likelihoods
		CacheAlignedDoubleVector likelihoods(forest.getNumSites());
		forest.computeLikelihoods(beb_set, likelihoods, dep.getDependencies());
		memcpy(&mPriors[iw*mNumSites], &likelihoods[0], forest.getNumSites()*sizeof(double));

		std::cerr << std::setw(3) << iw << ' ';
		std::cerr << std::fixed << std::setw(8) << std::setprecision(4) << omega_fg << ' ';
		std::cerr << std::fixed << std::setw(8) << std::setprecision(4) << omega_bg << " = " << std::setprecision(16);
		std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+0] << ' ';        // Print only the first three sites
		std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+1] << ' ';
		std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+2] << std::endl;
	}

	const std::vector<double>& site_mult = aModel.getSiteMultiplicity();
	double scale = 0.;
	for(size_t site=0; site < mNumSites; ++site)
	{
		size_t k;
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

void BayesTest::computeBEB(BranchSiteModelAltHyp& aModel, size_t aFgBranch)
{
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << std::endl << "Computing BEB for " << aFgBranch << std::endl;

	// Enable it as soon as the dependency list problem as been solved.
	double scale1 = getGridParams(aModel, aFgBranch);
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << std::endl << "Scale is " << scale1 << std::endl;

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
