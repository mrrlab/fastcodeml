
#ifndef WRITERESULTS_H
#define WRITERESULTS_H

#include <map>
#include <set>
#include <vector>
#include <utility>

/// Write all results to a given file.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-24 (initial version)
///     @version 1.1
///
///
class WriteResults {
public:
  /// Constructor
  ///
  /// @param[in] aFilename The filename to which the output should go. If null
  /// no printing will happen
  ///
  explicit WriteResults(const char *aFilename) : mFilename(aFilename) {}

  /// Destructor.
  ///
  ~WriteResults() {
    mLnL[0].clear();
    mLnL[1].clear();
    mPositiveSelSites.clear();
  }

  /// Output the results to the file given
  /// @param[in] aIbSet Set of internal branches
  /// @param[in] aBranchAll Whether fg branches are just internal or internals+leaves
  ///
  void outputResults(std::set<int> aIbSet, bool aBranchAll);

  /// Save the likelihood for later printing.
  ///
  /// @param[in] aFgBranch The foreground branch to which the log-likelihood
  /// refers
  /// @param[in] aLnL The log-likelihood
  /// @param[in] aHypothesis The hypothesis (0 for H0 and 1 for H1) for which
  /// the log-likelihood has been computed
  ///
  void saveLnL(size_t aFgBranch, double aLnL, unsigned int aHypothesis);


  /// Save the positive selection sites and corresponding probabilities for
  /// later printing. 
  ///
  /// @param[in] aFgBranch The foreground branch to which the sites refers
  /// @param[in] aPositiveSelSites The site in which the positive selection
  /// appears as computed by the BEB test
  /// @param[in] aPositiveSelSitesProb The corresponding probability as computed
  /// by the BEB test
  ///
  
  void savePositiveSelSites(size_t aFgBranch,
                            const std::vector<unsigned int> &aPositiveSelSites,
                            const std::vector<double> &aPositiveSelSitesProb);

  /// Check if anything will be output at the end.
  ///
  /// @return True if the output to file is enabled.
  ///
  
  bool isWriteResultsEnabled(void) const { return mFilename != NULL; }

  /// Returns the correct index order to have the sites printed in the correct
  /// order.
  ///
  /// @param[in] aSites The vector of site numbers
  ///
  /// @return The list of indices that gives the sites in the correct order.
  ///
  const std::vector<size_t> &
  orderSites(const std::vector<unsigned int> &aSites) const;

  /// Save the parameters string for later printing.
  ///
  /// @param[in] aFgBranch The foreground branch to which the log-likelihood
  /// refers
  /// @param[in] aParamStr The parameters string
  /// @param[in] aHypothesis The hypothesis (0 for H0 and 1 for H1) for which
  /// the log-likelihood has been computed
  ///
  void saveParameters(size_t aFgBranch, std::string &aParamStr,
                      unsigned int aHypothesis);

private:
  const char *mFilename; ///< The file to which the results should be written.
  /// If null, no printing appear
  std::map<size_t, double> mLnL[2]; ///< The log-likelihood for the given fg
  /// branch and for the two hypothesis
  std::map<size_t, std::pair<std::vector<unsigned int>, std::vector<double> > >
      mPositiveSelSites; ///< The sites under positive selection and the
  /// corresponding probabilities for a given fg branch
  mutable std::vector<size_t>
      mSiteOrder; ///< The new site+prob order computed by orderSites routine

  std::map<size_t, std::string> mParamStr[2]; ///< The parameters string for the
  /// given fg branch and for the two
  /// hypothesis
};



/// Write all results to a given file for the case of fg branches marked in nwk file.

class WriteResultsMfg {
public:
  /// Constructor
  ///
  /// @param[in] aFilename The filename to which the output should go. If null
  /// no printing will happen
  ///
  explicit WriteResultsMfg(const char *aFilename) : mFilename(aFilename) {}

  /// Destructor.
  ///
  ~WriteResultsMfg() {
    mLnL[0].clear();
    mLnL[1].clear();
    mPositiveSelSites.clear();
  }

  /// Output the results to the file given
  /// @param[in] aFgSet The set of marked fg branches
  ///
  void outputResults(std::set<int> aFgSet);

  /// Save the likelihood for later printing.
  ///
  /// @param[in] aFgBranchSet The marked foreground branches to which the log-likelihood
  /// refers
  /// @param[in] aLnL The log-likelihood
  /// @param[in] aHypothesis The hypothesis (0 for H0 and 1 for H1) for which
  /// the log-likelihood has been computed
  ///
  void saveLnL(std::set<int> aFgBranchSet, double aLnL, unsigned int aHypothesis);


  /// Save the positive selection sites and corresponding probabilities for
  /// later printing. 
  ///
  /// @param[in] aFgBranchSet The marked foreground branches to which the sites refers
  /// @param[in] aPositiveSelSites The site in which the positive selection
  /// appears as computed by the BEB test
  /// @param[in] aPositiveSelSitesProb The corresponding probability as computed
  /// by the BEB test
  ///
  
  void savePositiveSelSites(std::set<int> aFgBranchSet,
                            const std::vector<unsigned int> &aPositiveSelSites,
                            const std::vector<double> &aPositiveSelSitesProb);

  /// Check if anything will be output at the end.
  ///
  /// @return True if the output to file is enabled.
  ///
  
  bool isWriteResultsEnabled(void) const { return mFilename != NULL; }

  /// Returns the correct index order to have the sites printed in the correct
  /// order.
  ///
  /// @param[in] aSites The vector of site numbers
  ///
  /// @return The list of indices that gives the sites in the correct order.
  ///
  const std::vector<size_t> &
  orderSites(const std::vector<unsigned int> &aSites) const;

  /// Save the parameters string for later printing.
  ///
  /// @param[in] aFgBranchSet The marked foreground branches to which the log-likelihood
  /// refers
  /// @param[in] aParamStr The parameters string
  /// @param[in] aHypothesis The hypothesis (0 for H0 and 1 for H1) for which
  /// the log-likelihood has been computed
  ///
  void saveParameters(std::set<int> aFgBranchSet, std::string &aParamStr,
                      unsigned int aHypothesis);

private:
  const char *mFilename; ///< The file to which the results should be written.
  /// If null, no printing appear
  std::map<size_t, double> mLnL[2]; ///< The log-likelihood for the given marekd fg
  /// branches and for the two hypothesis
  std::map<size_t, std::pair<std::vector<unsigned int>, std::vector<double> > >
      mPositiveSelSites; ///< The sites under positive selection and the
  /// corresponding probabilities for a given marked fg branches
  mutable std::vector<size_t>
      mSiteOrder; ///< The new site+prob order computed by orderSites routine

  std::map<size_t, std::string> mParamStr[2]; ///< The parameters string for the
  /// given marked fg branches and for the two
  /// hypothesis
};


#endif
