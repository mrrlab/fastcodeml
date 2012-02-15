
#ifndef NEWICK_H
#define NEWICK_H

#include <string>
#include "TreeNode.h"
#include "PhyloTree.h"

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_parse_tree.hpp>

/// Newick specific functionalities.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-02-15 (initial version)
///     @version 1.0
///
class Newick : public PhyloTree
{
public:
	/// Constructor.
	///
	/// @param[in] aVerboseLevel Set the verbosity level
	/// - 0: No messages
	/// - 3: Print the tree structure
	///
	explicit Newick(unsigned int aVerboseLevel=0) : PhyloTree(aVerboseLevel) {}

	/// Load a phylo tree definition from a Newick formatted file.
	///
	/// @param[in] aFilename The filename
	///
	/// @exception FastCodeMLFatalNoMsg For errors like cannot open the file
	///
	virtual void loadTreeFile(const char *aFilename);

	/// Load a phylo tree definition from a Newick formatted string.
	///
	/// @param[in] aTreeAsString The string to be loaded
	///
	/// @exception FastCodeMLFatalNoMsg For errors like cannot open the file
	///
	virtual void loadTreeFromString(const std::string& aTreeAsString);

	/// Print the phylogenetic tree completed with all the info loaded in Newick format.
	///
	/// @param[in] aOut Output stream
	/// @param[in] aNode The node from which to start. If null starts from the root.
	///
	virtual void printTreeUnformatted(std::ostream& aOut, TreeNode *aNode=0) const;


private:

	/// Access to the parse tree
	///
	typedef boost::spirit::classic::tree_match<char const*> ParseTreeMatchType;
	typedef ParseTreeMatchType::tree_iterator ParseTreeIteratorType;

	/// Parse the parse tree and build the phylo tree structure.
	///
	/// @param[in] aTreeIterator Tree iterator
	/// @param[in,out] aNode Node from which the construction should continue
	///
	/// @exception FastCodeMLFatalNoMsg For errors like cannot open the file
	///
	void evaluateTreeNode(ParseTreeIteratorType const& aTreeIterator, TreeNode *aNode);

	/// Print an indented form of the parse tree.
	///
	/// @param[in] aTreeIterator Tree iterator
	/// @param[in] aIndent Indent level (each level increases by two spaces)
	///
	void printTree(ParseTreeIteratorType const& aTreeIterator, int aIndent);


private:
	std::string	mParsedPortion;		///< In case of error part of the string successfully parsed else it is empty.
};

#endif


