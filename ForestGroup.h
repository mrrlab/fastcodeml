#ifndef FORESTGROUP_H
#define FORESTGROUP_H

#include "Forest.h"
#include "CmdLine.h"
#include "Newick.h"
#include "Phylip.h"
#include "Exceptions.h"
#include "VerbosityLevels.h"

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
    ForestGroup() : mForests(), mCurrentForestIndex(0) {}

    /// Destructor
    ~ForestGroup();

    /// Initializer - initializes the data structures
    /// @param[in]  aCmdLine The command line object
    /// @exception  FastCodeMLFatal For no forests to analyze (no tree/alignment files)
    void initForests(const CmdLine &aCmdLine);

    /// Getter (non const) - returns non-const pointer to owned internal
    /// structures - do not delete. Forests are ordered as in CmdLine.mTreeFiles.
    ///
    /// @return Vector of pointers to the forests. Do not delete these.
    std::vector<Forest*> getForests();

    /// Solve a forest
    ///
    ///	@param[in,out] aForest   The forest to solve
    /// @param[in]     aCmdLine  The command line object
    /// @param[in]     aTreeFile The name of the tree file we are solving
    /// @param[in]     aGeneFile The name of the alignment file we are solving with
    /// @return        The results as a string from the WriteResults object
    std::string solveForest(Forest &aForest, const CmdLine &aCmdLine,
        const std::string &aTreeFile, const std::string &aGeneFile);

private:
    std::vector<Forest*> mForests;            ///< Defining member - a forest group is a vector of pointers to forests. These will be ordered as in the aCmdLine mTreeFiles/mGeneFiles objects after init.

    size_t               mCurrentForestIndex; ///< Counter for current forest number.

    /// Adds a forest to the forest group, for given tree file and gene file.
    /// @param[in]  cmd       The command line object
    /// @param[in]  treeFile  The tree file to add
    /// @param[in]  geneFile  The gene file associated with the tree, above
    void addForest(
        const CmdLine &cmd,
        const char *treeFile,
        const char *geneFile);

    /// Assignment operator is blocked since dependencies do not implement it.
    ForestGroup& operator=(const ForestGroup &other);

    /// Blocking copy constructor for the same reason as assignment, above.
    ForestGroup(const ForestGroup &other);
};

#endif


