#include <iostream>
#include <fstream>
#include "PhyloTree.h"
#include "Exceptions.h"

PhyloTree::~PhyloTree()
{
	mTreeRoot.clearNode();
	mLeavesSpecies.clear();
	mInternalNodes.clear();
}


#ifdef CHECK_ALGO
void PhyloTree::printFormattedTree(void) const
{
	mTreeRoot.printFormatted(0);
}
#endif


void PhyloTree::fillSpecies(TreeNode *aNode)
{
	if(aNode->isLeaf())
	{
		mLeavesSpecies.push_back(aNode);
	}
	else
	{
		TreeNode *m;
		for(int idx=0; (m = aNode->getChild(idx)) != 0; ++idx) fillSpecies(m);
	}
}

const std::vector<std::string>& PhyloTree::getSpecies(void) const
{
	static std::vector<std::string> species_list;
	species_list.clear();

	std::vector<TreeNode *>::const_iterator is = mLeavesSpecies.begin();
	const std::vector<TreeNode *>::const_iterator end = mLeavesSpecies.end();
	for(; is != end; ++is)
	{
		std::string label = (*is)->getLabel();
		species_list.push_back(label);
	}

	return species_list;
}

void PhyloTree::fillInternalBranches(TreeNode *aNode)
{
	if(!aNode->isLeaf() && aNode != &mTreeRoot)
	{
		mInternalNodes.push_back(aNode);
	}

	TreeNode *m;
	for(int idx=0; (m = aNode->getChild(idx)) != 0; ++idx) fillInternalBranches(m);
}


size_t PhyloTree::getMarkedInternalBranch(void) const
{
	const size_t nin = mInternalNodes.size();
	size_t marked_branch = 0;
	for(; marked_branch < nin; ++marked_branch)
	{
		if(!mInternalNodes[marked_branch]->getType().empty()) break;
	}

	if(marked_branch >= nin) return UINT_MAX;
	return marked_branch;
}



unsigned int PhyloTree::cloneTree(ForestNode* aForestNode, unsigned int aTreeId, size_t aNumSites, CacheAlignedDoubleVector& aProbVectors, const TreeNode* aTreeNode, unsigned int aNodeId) const
{
	unsigned int id;

	// Start with the tree root
	if(aTreeNode == 0)
	{
		aTreeNode = &mTreeRoot;
		aNodeId   = UINT_MAX;
		aForestNode->mParent = 0;
		id = 0;
	}
	else
	{
		id = aNodeId+1;
	}

	// Set the root values
	aForestNode->mBranchId = aNodeId;
	aForestNode->mOwnTree = aTreeId;

	// Set the internal branch identifier
	size_t int_id;
	const size_t nn = mInternalNodes.size();
	for(int_id=0; int_id < nn; ++int_id) if(aTreeNode == mInternalNodes[int_id]) break;
	aForestNode->mInternalNodeId = (int_id < nn) ? static_cast<unsigned int>(int_id) : UINT_MAX;

#ifndef NEW_LIKELIHOOD
	// Set the pointers.        The sequence is: Branch -> Set -> Site -> 1:N
	// Set the pointers (best). The sequence is: Branch -> Site -> Set -> 1:N
	for(int i=0; i < Nt; ++i)
	{
		//aForestNode->mProb[i] = &aProbVectors[VECTOR_SLOT*(aNumSites*Nt*id+aNumSites*i+aTreeId)]; 
		aForestNode->mProb[i] = &aProbVectors[VECTOR_SLOT*(aNumSites*Nt*id+aTreeId*Nt+i)]; 
	}
#endif

	// Recurse
	TreeNode *m;
	for(int idx=0; (m = aTreeNode->getChild(idx)) != 0; ++idx)
	{
		ForestNode* rn = new ForestNode;
		aForestNode->mChildrenList.push_back(rn);
		aForestNode->mChildrenCount++;
#ifndef NEW_LIKELIHOOD
		aForestNode->mOtherTreeProb.push_back(0);
#endif
		rn->mParent = aForestNode;
		id = cloneTree(rn, aTreeId, aNumSites, aProbVectors, m, id);
	}

	return id;
}

unsigned int PhyloTree::collectGlobalTreeData(std::vector<std::string>& aNodeNames, std::vector<double>& aBranchLengths, size_t* aMarkedIntBranch, const TreeNode* aTreeNode, unsigned int aNodeId) const
{
	unsigned int id;

	// Start with the tree root
	if(aTreeNode == 0)
	{
		aTreeNode = &mTreeRoot;
		aNodeId = UINT_MAX;
		id = 0;
		*aMarkedIntBranch = getMarkedInternalBranch();
	}
	else
	{
		id = aNodeId+1;
	}

	// Save the node values
	aNodeNames.push_back(aTreeNode->getLabel());
	aBranchLengths.push_back(aTreeNode->getLen());

	// Recurse
	TreeNode *m;
	for(int idx=0; (m = aTreeNode->getChild(idx)) != 0; ++idx)
	{
		id = collectGlobalTreeData(aNodeNames, aBranchLengths, aMarkedIntBranch, m, id);
	}

	return id;
}


