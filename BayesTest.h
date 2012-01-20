
#ifndef BAYESTEST_H
#define BAYESTEST_H

#include <vector>

/// Tests to find the sites under positive selection.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
///
class BayesTest
{
public:
	/// Constructor.
	///
	/// @param[in] aNumSites Number of sites
	///
	BayesTest(size_t aNumSites);

	/// Destructor.
	///
	~BayesTest() {}
	
	/// Bayes Empirical Bayes (BEB) test.
	///
	/// @todo Missing computeBEB routine. The values that are output are simulated.
	///
	void computeBEB(void);
	
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
	std::vector<double> mSiteClassProb;		///< Probability of a site to pertain to a given class (one row per class (4 classes), one column per site).
	size_t				mNumSites;			///< Number of sites.
};

#endif

