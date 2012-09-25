
#ifndef WRITERESULTS_H
#define WRITERESULTS_H

#include <map>
#include <vector>
#include <utility>

/// Write all results to a given file.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-24 (initial version)
///     @version 1.0
///
///
class WriteResults
{
public:
	/// Constructor
	///
	/// @param[in] aFilename The filename to which the output should go
	///
	explicit WriteResults(const char* aFilename) : mFilename(aFilename) {}

	/// Output the results to the file given or to the stdout
	///
	/// @param[in] aOutputToStdout If true the output is not written to the file but redirect to the screen
	///
	void outputResults(bool aOutputToStdout=false);

	/// Save the likelihood for later printing.
	///
	/// @param[in] aFgBranch The foreground branch to which the log-likelihood refers
	/// @param[in] aLnL The log-likelihood
	/// @param[in] aHypothesis The hypothesis (0 for H0 and 1 for H1) for which the log-likelihood has been computed
	///
	void saveLnL(size_t aFgBranch, double aLnL, unsigned int aHypothesis);

	/// Save the positive selection sites and corresponding probabilities for later printing.
	///
	/// @param[in] aFgBranch The foreground branch to which the sites refers
	/// @param[in] aPositiveSelSites The site in which the positive selection appears as computed by the BEB test
	/// @param[in] aPositiveSelSitesProb The corresponding probability as computed by the BEB test
	///
	void savePositiveSelSites(size_t aFgBranch, const std::vector<unsigned int>& aPositiveSelSites, const std::vector<double>& aPositiveSelSitesProb);


private:
	const char*					mFilename;			///< The file to which the output should be redirect. If null, no printing appear
	std::map<size_t, double>	mLnL[2];			///< The log-likelihood for the given fg branch and for the two hypothesis
	std::map<size_t, std::pair<std::vector<unsigned int>, std::vector<double> > >
								mPositiveSelSites;	///< The sites under positive selection and the corresponding probabilities for a given fg branch
};

#endif

