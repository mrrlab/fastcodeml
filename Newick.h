
#ifndef NEWICK_H
#define NEWICK_H

#include <string>
#include "TreeNode.h"
#include "PhyloTree.h"

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
//#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wlong-long"
#endif

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_parse_tree.hpp>

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
//#pragma GCC diagnostic pop
#endif

/// %Newick format file specific functionalities.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-02-15 (initial version)
///     @version 1.1
///
class Newick : public PhyloTree {
public:
  /// Constructor.
  ///
  /// @param[in] aVerboseLevel Set the verbosity level
  /// - 0: No messages
  /// - 3: Print the tree structure
  ///
  explicit Newick(unsigned int aVerboseLevel = 0) : PhyloTree(aVerboseLevel) {}

  /// Load a phylo tree definition from a Newick formatted file.
  ///
  /// @param[in] aFilename The filename
  ///
  /// @exception FastCodeMLFatal For errors like cannot open the file or cannot
  /// parse the tree
  ///
  virtual void readFile(const char *aFilename);

  /// Print the phylogenetic tree completed with all the info loaded in Newick
  /// format.
  ///
  /// @param[in] aOut Output stream
  /// @param[in] aNode The node from which to start. If null starts from the
  /// root.
  ///
  virtual void printTreeUnformatted(std::ostream &aOut,
                                    TreeNode *aNode = NULL) const;

  /// Print the phylogenetic tree completed with all the info loaded in the same
  /// format as read in and annotated with the branch numbers.
  ///
  /// @param[in] aOut Output stream
  /// @param[in] aNode The node from which to start. If null starts from the
  /// root.
  /// @param[in] aInternalBranch Internal branch identifier to annotate the
  /// current branch.
  /// @param[in] whether leaves should be labeled or not.
  /// @param[in] whether branch numbers should be shown or not.
  ///
  /// @return The new branch id
  ///
  virtual int printTreeAnnotated(std::ostream &aOut, TreeNode *aNode = NULL,
                                 int aInternalBranch = 0,
                                 bool wLeaves = false, bool bNumber=true) const;

  /// Print the phylogenetic tree completed with all the info loaded in the same
  /// format as read in and annotated with the branch numbers.
  ///
  /// @param[in] aOut Output stream
  /// @param[in] aNode The node from which to start. If null starts from the
  /// root.
  /// @param[in] aInternalBranch Internal branch identifier to annotate the
  /// current branch.
  /// @param[in] whether leaves should be labeled or not.
  /// @param[in] array of branch lengths after estimation to be labeled in the
  /// tree.
  /// @param[in] whether branch numbers should be shown or not.
  ///
  /// @return The new branch id
  ///
  virtual int
  printTreeAnnotatedWithEstLens(std::ostream &aOut, TreeNode *aNode = NULL,
                                int aInternalBranch = 0, bool wLeaves = false,
                                std::vector<double> *mVar = NULL, bool bNumber=true) const;

private:
  /// Load a phylo tree definition from a Newick formatted string.
  ///
  /// @param[in] aTreeAsString The string to be loaded
  ///
  /// @exception FastCodeMLFatal For parsing errors
  ///
  void loadTreeFromString(const std::string &aTreeAsString);

  /// Access to the Boost::Spirit parse tree
  ///
  typedef boost::spirit::classic::tree_match<char const *>::tree_iterator
      ParseTreeIteratorType;

  /// Parse the parse tree and build the phylo tree structure.
  ///
  /// @param[in] aTreeIterator Tree iterator
  /// @param[in,out] aNode Node from which the construction should continue
  ///
  /// @exception FastCodeMLFatal For parsing errors
  ///
  void evaluateTreeNode(ParseTreeIteratorType const &aTreeIterator,
                        TreeNode *aNode);

  /// Print an indented form of the parse tree.
  ///
  /// @param[in] aTreeIterator Tree iterator
  /// @param[in] aIndent Indent level (each level increases by aIndentIncrement)
  /// @param[in] aIndentIncrement Indent level increment at each level
  ///
  void printTree(ParseTreeIteratorType const &aTreeIterator,
                 unsigned int aIndent = 0, unsigned int aIndentIncrement = 2);

private:
  std::string mParsedPortion; ///< In case of error part of the string
  /// successfully parsed else it is empty.
};

#endif
