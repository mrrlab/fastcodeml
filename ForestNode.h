#ifndef FORESTNODE_H
#define FORESTNODE_H

#include <vector>
#include <iostream>
#include <climits>
#include <cstring>
#include "MatrixSize.h"


/// One node of the forest generated by Forest.
/// The node name is stored in Forest.mNodeNames and the
/// length of the branch leading to this node as read from the file is in Forest.mBranchLengths
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-02-23 (initial version)
///     @version 1.0
///
///
struct ForestNode
{
	unsigned int				mInternalNodeId;			///< Internal node identifier to mark a branch as foreground. UINT_MAX means not an internal node
	unsigned int				mNodeId;					///< An unique index to access the branch length array (starts from zero at the first non-root node)
	unsigned int				mOwnTree;					///< Per tree identifier
	ForestNode*					mParent;					///< Pointer to the node parent (null for the root)
#ifndef NEW_LIKELIHOOD
	//double						mProb0[N*Nt];				///< Codons probability array (called g in the pseudocode) (can be computed by concurrent tree traversals)
	double*						mProb[Nt];
#endif
	std::vector<ForestNode *>	mChildrenList;				///< List of the node children
	std::vector<int>			mSubtreeCodonsSignature;	///< List of codon idx for the subtree rooted at this node (after reduction it is emptied)
#ifndef NEW_LIKELIHOOD
	std::vector<double *>		mOtherTreeProb;				///< Pointers to other tree precomputed mProb, zero if not used, or local array if used from other tree
#endif
	/// Constructor
	///
	ForestNode()
	{
		mParent = 0;
#ifndef NEW_LIKELIHOOD
		memset(mProb, 0, Nt*sizeof(double*));
#endif
	}

	/// Destructor
	///
	~ForestNode()
	{
		// Delete children if in the same tree. Delete partial Prob arrays if not pointer to other tree partial Prob array
		std::vector<ForestNode*>::iterator irn;
		unsigned int i;
		for(irn=mChildrenList.begin(), i=0; irn != mChildrenList.end(); ++irn, ++i)
		{
			if(mOwnTree == (*irn)->mOwnTree)
			{
				delete (*irn);

#ifndef NEW_LIKELIHOOD
				delete [] mOtherTreeProb[i];
#endif
			}
		}

		// Clean all arrays
		mChildrenList.clear();
		mSubtreeCodonsSignature.clear();
#ifndef NEW_LIKELIHOOD
		mOtherTreeProb.clear();
#endif
	}

	/// Copy constructor and assignment
	///
	/// @param[in] aNode Node that has to be assigned to the current node
	///
	ForestNode(const ForestNode& aNode)
	{
		mChildrenList  = aNode.mChildrenList;
		mParent        = aNode.mParent;
		mSubtreeCodonsSignature = aNode.mSubtreeCodonsSignature;
#ifndef NEW_LIKELIHOOD
		memcpy(mProb, aNode.mProb, Nt*sizeof(double*));
#endif
		mInternalNodeId = aNode.mInternalNodeId;
		mNodeId         = aNode.mNodeId;
		mOwnTree        = aNode.mOwnTree;
#ifndef NEW_LIKELIHOOD
		mOtherTreeProb  = aNode.mOtherTreeProb;
#endif
	}

	/// Assignment operator
	///
	/// @param[in] aNode Node that has to be assigned to the current node
	///
	/// @return The node itself
	///
	ForestNode& operator=(const ForestNode& aNode)
	{
		// Make sure not same object
		if(this != &aNode)
		{
			mChildrenList  = aNode.mChildrenList;
			mParent        = aNode.mParent;
			mSubtreeCodonsSignature = aNode.mSubtreeCodonsSignature;
#ifndef NEW_LIKELIHOOD
			memcpy(mProb, aNode.mProb, Nt*sizeof(double*));
#endif
			mInternalNodeId = aNode.mInternalNodeId;
			mNodeId         = aNode.mNodeId;
			mOwnTree        = aNode.mOwnTree;
#ifndef NEW_LIKELIHOOD
			mOtherTreeProb  = aNode.mOtherTreeProb;
#endif
		}

		// Return ref for multiple assignment
		return *this;
	}

	/// Print from this node down
	///
	/// @param[in] aNodeNames The list of node names
	/// @param[in] aOut Output stream
	/// @param[in] aIndent Initial number of indent spaces
	/// @param[in] aIncrement The indent amount is incremented by this value at each level
	///
	void print(const std::vector<std::string>& aNodeNames, std::ostream& aOut=std::cerr, unsigned int aIndent=0, unsigned int aIncrement=3) const
	{
		unsigned int i;

		// Indent
		for(i=0; i < aIndent; ++i) aOut << ' ';

		// Print the name
		aOut << '<' << ((mNodeId  != UINT_MAX) ? aNodeNames[mNodeId+1] : aNodeNames[0]) << "> ";
	
		// Print the ID
		if(mInternalNodeId != UINT_MAX) aOut << '(' << mInternalNodeId << '|' << mNodeId << ") ";
		else                            aOut << '('                    << '|' << mNodeId << ") ";

		// Print the indexes of the codons accumulated till this node
		std::vector<int>::const_iterator ig;
		for(ig=mSubtreeCodonsSignature.begin(); ig != mSubtreeCodonsSignature.end(); ++ig) aOut << *ig << ' ';
		aOut << std::endl;

		// Print the subtree
		std::vector<ForestNode*>::const_iterator irn;
		for(irn=mChildrenList.begin(), i=0; irn != mChildrenList.end(); ++irn, ++i)
		{
			// If the subtree is on the same tree, then print it, otherwise print only the subtree root node name.
			if(mOwnTree == (*irn)->mOwnTree)
			{
				(*irn)->print(aNodeNames, aOut, aIndent+aIncrement, aIncrement);
			}
			else
			{
				for(i=0; i < aIndent+aIncrement; ++i) aOut << ' ';
				i = (*irn)->mNodeId;
				if(i == UINT_MAX) i = 0;
				aOut << '[' << aNodeNames[i] << ']' << std::endl;
			}
		}
	}

	/// Create a list of pointers to leaves.
	///
	/// @param[out] aLeafsList Pointers to leaves are pushed to this vector
	///
	void pushLeaf(std::vector<ForestNode*>& aLeafsList)
	{
		if(mChildrenList.empty())
		{
			aLeafsList.push_back(this);
		}
		else
		{
			std::vector<ForestNode*>::const_iterator irn;
			for(irn=mChildrenList.begin(); irn != mChildrenList.end(); ++irn) (*irn)->pushLeaf(aLeafsList);
		}
	}

	/// Fills the mSubtreeCodonsSignature list with the ordered union of its children's lists.
	///
	void gatherCodons(void)
	{
		std::vector<ForestNode*>::const_iterator irn;
		for(irn=mChildrenList.begin(); irn != mChildrenList.end(); ++irn)
		{
			(*irn)->gatherCodons();
			mSubtreeCodonsSignature.insert(mSubtreeCodonsSignature.end(), (*irn)->mSubtreeCodonsSignature.begin(), (*irn)->mSubtreeCodonsSignature.end());
		}
	}

	/// Count the total branches in the forest
	///
	/// @param[in] aAggressiveStrategy If true use the aggressive simplification strategy
	///
	unsigned int countBranches(bool aAggressiveStrategy=false) const
	{
		unsigned int cnt = 0;
		unsigned int i;

		// Visit the subtrees
		std::vector<ForestNode*>::const_iterator irn;
		for(irn=mChildrenList.begin(), i=0; irn != mChildrenList.end(); ++irn, ++i)
		{
			// If the subtree is on the same tree, then print it, otherwise print only the subtree root node name.
			if(mOwnTree == (*irn)->mOwnTree)
			{
				cnt += (*irn)->countBranches(aAggressiveStrategy)+1;
			}
			else if(!aAggressiveStrategy)
			{
				cnt += 1;
			}
		}
	
		return cnt;
	}

	/// Tests if the given child is in the same tree or not
	///
	/// @param[in] aChildIdx Index of the child in the list of children
	///
	/// @return True if it is in the same tree
	///
	inline bool isSameTree(unsigned int aChildIdx)
	{
		return (mChildrenList[aChildIdx]->mOwnTree == mOwnTree);
	}
};

#endif

