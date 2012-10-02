
#ifndef BAYESTEST_H
#define BAYESTEST_H

#include <vector>
#include <cstdlib>
#include "BranchSiteModel.h"

/// Few constants
const static unsigned int BEB_N1D = 10;	///< Number of categories for w0 and w2
const static unsigned int BEB_DIMS = 4;	///< Number of codon classes (0, 1, 2a, 2b)
const static unsigned int BEB_NUM_CAT = BEB_N1D + 1 + BEB_N1D*BEB_N1D + BEB_N1D; ///< Total number of categories for w0 and w2 (it is com.ncatG in codeml.c)

/// Tests to find the sites under positive selection.
///
///  @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///  @date 2010-12-22 (initial version)
///  @version 1.0
///
class BayesTest
{
public:
	/// Constructor.
	///
	/// @param[in] aNumSites Number of sites
	/// @param[in] aVerbose The verbosity level
	///
	explicit BayesTest(size_t aNumSites, unsigned int aVerbose=0);

	/// Destructor.
	///
	~BayesTest() {}
	
	/// Bayes Empirical Bayes (BEB) test.
	///
	/// @todo Missing computeBEB routine. The values that are output are simulated.
	///
	void computeBEB(BranchSiteModelAltHyp& aModel);
	
	/// Print the sites under positive selection.
	///
	/// @param[in] aFgBranch Identifier of the branch marked as foreground branch
	///
	void printPositiveSelSites(size_t aFgBranch) const;

	/// Extract the sites under positive selection and the corresponding probabilities.
	///
	/// @param[out] aPositiveSelSites Vector of sites under positive selection
	/// @param[out] aPositiveSelSitesProb Corresponding probabilities
	///
	void extractPositiveSelSites(std::vector<unsigned int>& aPositiveSelSites, std::vector<double>& aPositiveSelSitesProb) const;

private:
	/// Generate a double random number between 0 and 1
	///
	/// @return The random number
	///
	static inline double randFrom0to1(void) {return static_cast<double>(rand())/static_cast<double>(RAND_MAX);}

	/// This sets up the grid (mPara[][]) according to the priors.  
	/// It calculates the probability of data at each site given w: f(f_h|w).  
	/// This is calculated using the branch model (NSsites = 0 model = 2), with 
	/// BayesEB=2 used to force the use of the correct scale factors in GetPMatBranch().
	///
	///@verbatim
	/// Order of site classes for iw or f(x_h|w):
	///                     fore   back     num.sets
	/// Branchsite A (121 sets)
	///   site class 0:      w0     w0        10
	///   site class 1:      w1=1   w1=1       1
	///   site class 2a:     w0     w2       100
	///   site class 2b:     w1=1   w2        10
	///@endverbatim
	///
	double getGridParams(BranchSiteModelAltHyp& aModel);

private:
	std::vector<double> mSiteClassProb;				///< Probability of a site to pertain to a given class (one row per class (4 classes), one column per site).
	size_t				mNumSites;					///< Number of sites.
	double				mPara[BEB_DIMS][BEB_N1D];	///< Parameters for w0, w1, w2 prior computation
	unsigned int		mVerbose;					///< If greather than zero prints more info
	std::vector<double>	mPriors;					///< Computed priors (each points to a list, one for each site)
};

#endif

