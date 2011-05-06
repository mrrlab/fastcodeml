
#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#include <string>
#include <vector>
#include "TreeNode.h"
#include "ForestNode.h"
#include "TransitionMatrix.h"

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_parse_tree.hpp>

typedef boost::spirit::classic::tree_match<char const*> ParseTreeMatchType;
typedef ParseTreeMatchType::tree_iterator ParseTreeIteratorType;

/// Phylogenetic tree public interface.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-08-30 (initial version)
///     @version 1.0
///
class PhyloTree
{
public:
	/// Constructor.
	///
	/// @param[in] aVerboseLevel Set the verbosity level
	/// - 0: No messages
	/// - 3: Print the tree structure
	///
	PhyloTree(int aVerboseLevel=0);

	/// Destructor.
	///
	~PhyloTree();

	/// Load a phylo tree definition from a Newick formatted file. Errors are signalled by exceptions.
	///
	/// @param[in] aFilename The filename
	///
	void loadTree(const char *aFilename);

	/// Load a phylo tree definition from a Newick formatted string. Errors are signalled by exceptions.
	///
	/// @param[in] aTreeAsString The string to be loaded
	///
	void loadTreeFromSting(const std::string& aTreeAsString);

	/// Show the parsing error point. The output is valid only if loadTree() ends in error.
	///
	/// @return The text parsed so far without error
	///
	std::string loadTreeParseError(void) const {return mParsedPortion;}

	/// Print the phylogenetic tree completed with all the info loaded.
	///
	void printFormattedTree(void) const;

	/// Print the phylogenetic tree completed with all the info loaded in Newick format.
	///
	/// @param[in] aNode The node from which to start. If null starts from the root.
	///
	void printNewickTree(TreeNode *aNode=0) const;

	/// Load the list of species in the given array.
	///
	/// @param[out] aSpeciesList Vector of species names as read at the leaves of the tree
	///
	void getSpecies(std::vector<std::string>& aSpeciesList) const;

	/// Number of tree branches.
	///
	/// @return The number of tree branches
	///
	unsigned int getNumBranches(void) const {return mInternalNodes.size()+mLeavesSpecies.size();}

	/// Get the marker in the Newick file (the one starting with '#') for the given internal branch.
	///
	/// @param[in] aInternalBranchIdx Index of the internal branch to access.
	///
	/// @return The string of the marker or an empty string.
	///
	const std::string getMarkerOnNode(unsigned int aInternalBranchIdx) const {return mInternalNodes[aInternalBranchIdx]->getType();}

	/// Clone the tree using ForestNode. Called without aTreeNode starts from the tree root.
	///
	/// @param[out] aForestNode The ForestNode that becomes the root of the cloned tree
	/// @param[in] aTreeId The tree running id.
	/// @param[in] aTreeNode The node from which to start the cloning in the tree. If not present starts from the root
	/// @param[in] aNodeId The node running id. For the root it is UINT_MAX.
	///
	/// @return The node id to the next node
	///
	unsigned int cloneTree(ForestNode* aForestNode, unsigned int aTreeId, const TreeNode* aTreeNode=0, unsigned int aNodeId=0) const;

private:
	/// Parse the parse tree and build the phylo tree structure.
	///
	/// @param[in] aTreeIterator Tree iterator
	/// @param[in,out] aNode Node from which the construction should continue
	///
	void evaluateTreeNode(ParseTreeIteratorType const& aTreeIterator, TreeNode *aNode);

	/// Print an indented form of the parse tree.
	///
	/// @param[in] aTreeIterator Tree iterator
	/// @param[in] aIndent Indent level (each level increases by two spaces)
	///
	void printTree(ParseTreeIteratorType const& aTreeIterator, int aIndent);

	/// Fill the list of Species (the leaves of the tree).
	///
	/// @param[in,out] aNode The tree node (recursively called function starting from the tree root)
	///
	void fillSpecies(TreeNode *aNode);

	/// Fill the list of internal nodes.
	///
	/// @param[in,out] aNode The tree  node (recursively called function starting from the tree root)
	///
	void fillInternalBranches(TreeNode *aNode);

private:
	std::string				mParsedPortion;		///< In case of error part of the string successfully parsed else it is empty.
	TreeNode				mTreeRoot;			///< The root of the phylogenetic tree in memory
	int						mVerboseLevel;		///< The verbosity level
	std::vector<TreeNode *>	mLeavesSpecies;		///< The list of the tree leaves
	std::vector<TreeNode *> mInternalNodes;		///< The list of the tree internal nodes
};

#endif


