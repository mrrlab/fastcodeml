
#include <iostream>
#include <cstring>
#include "TreeNode.h"


void TreeNode::printFormatted(int aIndent) const
{
	// Indent the current entry
	for(int k=0; k < aIndent; ++k) std::cerr << ' ';

	// Print the node info
	if(mNodeMark.empty())
		std::cerr << '<' << mNodeName << "> " << mBranchLength;
	else
		std::cerr << '<' << mNodeName << ">#" << mNodeMark << ' ' << mBranchLength;

	// End line
	std::cerr << std::endl;

	// Print the children (indented by 3 spaces)
	std::vector<TreeNode *>::const_iterator in=mChildrenList.begin();
	const std::vector<TreeNode *>::const_iterator end=mChildrenList.end();
	for(; in != end; ++in) (*in)->printFormatted(aIndent+3);
}


void TreeNode::printNode(void) const
{
	// Print the node info
	std::cerr << mNodeName;
	if(!mNodeMark.empty()) std::cerr << '#' << mNodeMark;
}


void TreeNode::clearNode(void)
{
	// Recursively remove the nodes
	if(!mChildrenList.empty())
	{
		std::vector<TreeNode *>::iterator in = mChildrenList.begin();
		const std::vector<TreeNode *>::iterator end = mChildrenList.end();
		for(; in != end; ++in) (*in)->clearNode();
		for(in=mChildrenList.begin(); in != end; ++in) delete (*in);
	}

	// Special case for the root node
	if(!mParent)
	{
		mChildrenList.clear();
	}
}

