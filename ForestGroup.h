#ifndef FORESTGROUP_H
#define FORESTGROUP_H

#include "Forest.h"
#include "CmdLine.h"
#include "Newick.h"
#include "Phylip.h"
#include "Exceptions.h"
#include "VerbosityLevels.h"
#include "BranchSiteModel.h"
#include "BayesTest.h"
#include "CodonFrequencies.h"

/// Defines a group of forests and functions associated with them.
///
/// This corresponds to a set of input file <Tree, alignment> file pairs.
/// This class is intended to be the top level class to manage anything related
/// to groups of forests or actions on them. As such, it manages the single threaded
/// solution procedure.
///
///
///     @author Kareem Ali - UNIL
///     @date 2014-06-18 (initial version)
///     @version 1.1.1
///
///

class ForestGroup
{
public:
    /// Constructor - default initializes the data structures.
    ForestGroup() : mForests(), mCodonFrequencies(), mNullHypModels(),
        mAltHypModels(), mBayesTests(), mCurrentForestIndex(0) {}

    /// Destructor
    ~ForestGroup();

    /// Initializer - initializes the data structures
    ///
    /// @param[in]  aCmdLine The command line object
    /// @exception  FastCodeMLFatal For no forests to analyze (no tree/alignment files)
    void initForests(const CmdLine &aCmdLine);

    /// Getter (non const) - returns non-const pointer to owned internal
    /// structures - do not delete. Forests are ordered as in CmdLine.mTreeFiles.
    ///
    /// @return Vector of pointers to the forests. Do not delete these.
    std::vector<Forest*> getForests();

    /// Get the pointer to the forest at index ii (below functions keyed on the same index ii).
    ///
    /// @return Pointer to the forest at index ii
    const Forest* getForestPtr(size_t ii) const {return mForests[ii]; }

    /// Get the codon frequencies object for the forest at index ii
    ///
    /// @return Pointer to the codon frequencies object at index ii
    CodonFrequencies* getCodonFrequenciesPtr(size_t ii) { return mCodonFrequencies[ii]; }

    /// Get a pointer to the null hypothesis test for forest at index ii
    ///
    /// @return Pointer to the branch site null hypothesis test object at index ii
    BranchSiteModelNullHyp* getNullHypothesisTest(size_t ii) { return mNullHypModels[ii]; }

    /// Get a pointer to the alternative hypothesis test for forest at index ii
    ///
    /// @return Pointer to the branch site alternative hypothesis test object at index ii
    BranchSiteModelAltHyp* getAltHypothesisTest(size_t ii) { return mAltHypModels[ii]; }

    /// Get a pointer to the bayes hypothesis test for forest at index ii
    ///
    /// @return Pointer to the bayes hypothesis test object at index ii
    BayesTest* getBayesTest(size_t ii) { return mBayesTests[ii]; }

    /// Get the size of the forest group
    ///
    /// @return The size of the forest group (the number of species trees / alignment pairs).
    size_t size() const {return mForests.size(); }

    /// Solve a forest
    ///
    ///	@param[in,out] aForest           The forest to solve
    /// @param[in]     aForestIndex      The index of the forest in the forest group.
    /// @param[in]     aCmdLine          The command line object
    ///
    /// @return        The results as a string from the WriteResults object
    std::string solveForest(
        Forest &aForest, size_t aForestIndex, const CmdLine &aCmdLine);

    /// Output the total number of internal branches in the group
    ///
    /// @param[in] The command line object
    ///
    /// @return    The total number of internal branches
    size_t getTotalNumInternalBranches(const CmdLine &aCmdLine) const;

private:
    std::vector<Forest*>                 mForests;            ///< Defining member - a forest group is a vector of pointers to forests. These will be ordered as in the aCmdLine mTreeFiles/mGeneFiles objects after init.
    std::vector<CodonFrequencies*>       mCodonFrequencies;   /// The codon frequencies objects, associated with the forests. Forest at index 1 corresponds to codon frequencies at index 1 and so on.

    std::vector<BranchSiteModelNullHyp*> mNullHypModels;      /// Null hypothesis branch site models
    std::vector<BranchSiteModelAltHyp*>  mAltHypModels;       /// Alternative hypothesis site models
    std::vector<BayesTest*>              mBayesTests;         /// The bayes test objects

    size_t                               mCurrentForestIndex; ///< Counter for current forest number (incremented in solveForest)

    /// Adds a forest to the forest group, for given tree file and gene file.
    ///
    /// @param[in]  cmd       The command line object
    /// @param[in]  treeFile  The tree file to add
    /// @param[in]  geneFile  The gene file associated with the tree, above
    void addForest(const CmdLine &cmd, const char *treeFile, const char *geneFile);

    /// Assignment operator is blocked since dependencies do not implement it.
    ForestGroup& operator=(const ForestGroup &other);

    /// Blocking copy constructor for the same reason as assignment, above.
    ForestGroup(const ForestGroup &other);
};

#endif


