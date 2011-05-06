
#ifndef TREENODE_H
#define TREENODE_H

#include <string>
#include <vector>
#include "MatrixSize.h"

/// Internal representation of the Phylogenetic tree
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-08-31 (initial version)
///     @version 1.0
///
class TreeNode
{
public:
	/// Constructor.
	///
	TreeNode() : mBranchLength(0.), mParent(0) {}

	/// Destructor.
	///
	~TreeNode() {}

	/// Add the node label.
	///
	/// @param[in] aLabel The node label
	///
	void addLabel(const std::string& aLabel) {mNodeName = aLabel;}

	/// Add the node type (the string after the '#' symbol in label).
	///
	/// @param[in] aMark The node type string
	///
	void addType(const std::string& aMark) {mNodeMark = aMark;}

	/// Add the branch length.
	///
	/// @param[in] aBranchLength The branch length
	///
	void addLen(double aBranchLength) {mBranchLength = aBranchLength;}

	/// Add a child to the node.
	///
	/// @param[in] aNode The child node
	///
	void addChild(TreeNode* aNode) {mChildrenList.push_back(aNode);}

	/// Add the node parent.
	///
	/// @param[in] aNode The node parent
	///
	void addParent(TreeNode* aNode) {mParent = aNode;}

	/// Print the node info indented by the amount requested.
	///
	/// @param[in] aIndent The number of spaces to indent the printing (each level increases by 3)
	///
	void printFormatted(int aIndent) const;

	/// Print the node in the format of Newick tree.
	///
	void printNode(void) const;

	/// Clear the node.
	///
	void clearNode(void);

	/// Get the node label.
	///
	/// @return The node label
	///
	const std::string getLabel(void) const {return mNodeName;}

	/// Get the node type from file.
	///
	/// @return The node marker
	///
	const std::string getType(void) const {return mNodeMark;}

	/// Check if the node is a leaf.
	///
	/// @return True if the node has no descendants
	///	
	bool isLeaf(void) const {return mChildrenList.empty();}

	/// Return one of the node children.
	///
	/// @param[in] aIdx Index of the child to be returned
	/// @return Pointer to the child node or NULL if the index is out of range
	///
	TreeNode* getChild(unsigned int aIdx) const
	{
		if(aIdx >= mChildrenList.size()) return 0;
		return mChildrenList[aIdx];
	}

	/// Get the branch length.
	///
	/// @return The branch length
	///
	double getLen(void) const { return mBranchLength; }


private:
	std::vector<TreeNode *>	mChildrenList;	///< List of the node children
	std::string				mNodeName;		///< Node label
	std::string				mNodeMark;		///< Node type or empty if the '#' part is not present
	double					mBranchLength;	///< Length of the branch leading to this node as read from the file (not valid for the root)
	TreeNode*				mParent;		///< Pointer to the node parent (null for the root)
};

#endif

