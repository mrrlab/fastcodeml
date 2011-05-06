
#ifndef BAYESTEST_H
#define BAYESTEST_H

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
	BayesTest(unsigned int aNumSites);

	/// Destructor.
	///
	~BayesTest();
	
	/// Bayes Empirical Bayes (BEB) test.
	///
	/// @todo Missing computeBEB routine. The values that are output are simulated.
	///
	void computeBEB(void);

	/// Naive Empirical Bayes (NEB) test.
	///
	/// @todo Missing computeNEB routine
	///
	void computeNEB(void);
	
	/// Print the sites under positive selection.
	///
	/// @param[in] aFgBranch Identifier of the branch marked as foreground branch
	///
	void printPositiveSelSites(unsigned int aFgBranch);


private:
	double* mSiteClassProb;		///< Probability of a site to pertain to a given class (one row per class (4 classes), one column per site).
	unsigned int mNumSites;		///< Number of sites.
};

#endif

