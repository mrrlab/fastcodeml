#include <iostream>
#include <fstream>
#include <cstring>
#include <climits>
#include "PhyloTree.h"
#include "NewickGrammar.h"
#include "MatrixSize.h"
#include "MathSupport.h"
#include "Exceptions.h"

void PhyloTree::printTree(ParseTreeIteratorType const& aTreeIterator, int aIindent)
{
	int k;

	for(k=0; k < aIindent; ++k) std::cout << " ";

    if(aTreeIterator->value.id() == NewickGrammar::treeID)           std::cout << "TREE";
    else if(aTreeIterator->value.id() == NewickGrammar::nodelistID)  std::cout << "NODELST";
    else if(aTreeIterator->value.id() == NewickGrammar::subtreeID)   std::cout << "SUBTREE";
    else if(aTreeIterator->value.id() == NewickGrammar::fulllabelID) std::cout << "FLAB";
    else if(aTreeIterator->value.id() == NewickGrammar::labelID)     std::cout << "LABEL";
    else if(aTreeIterator->value.id() == NewickGrammar::typeID)      std::cout << "TYPE";
    else if(aTreeIterator->value.id() == NewickGrammar::branchlenID) std::cout << "BRANCHLEN";
    else if(aTreeIterator->value.id() == NewickGrammar::cblenID)     std::cout << "CBLEN";
    else															 std::cout << "????";

	if(aTreeIterator->value.begin() != aTreeIterator->value.end())
	{
		std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
		std::cout << " " << label_name << std::endl;
	}
	else
	{
		std::cout << std::endl;
	}

	for(ParseTreeIteratorType ic=aTreeIterator->children.begin(); ic != aTreeIterator->children.end(); ++ic) printTree(ic, aIindent+2);
}


void PhyloTree::evaluateTreeNode(ParseTreeIteratorType const& aTreeIterator, TreeNode *aNode)
{
    unsigned int k;

    if(aTreeIterator->value.id() == NewickGrammar::treeID)
    {
        //std::cout << "tree (" << std::endl;
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
        //std::cout << ")" << std::endl;
    }
    else if(aTreeIterator->value.id() == NewickGrammar::nodelistID)
    {
        //std::cout << "node_list {" << std::endl;
		//TreeNode* x = new TreeNode;
		//n->addChild(x);
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
        //std::cout << "}" << std::endl;
    }
    else if(aTreeIterator->value.id() == NewickGrammar::subtreeID)
    {
        //std::cout << "subtree [" << std::endl;
		TreeNode* x = new TreeNode;
		aNode->addChild(x);
		x->addParent(aNode);
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, x);
        //std::cout << "]" << std::endl;
    }
    else if(aTreeIterator->value.id() == NewickGrammar::fulllabelID)
    {
        //std::cout << "full_label " << std::endl;
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    }
    else if(aTreeIterator->value.id() == NewickGrammar::labelID)
    {
        if(aTreeIterator->children.size() == 0)
		{
			std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
            //std::cout << "label " << label_name << std::endl;
			aNode->addLabel(label_name);
        }
		else
		{
			for(k=0; k < aTreeIterator->children.size(); ++k)
            	evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    	}
	}
    else if(aTreeIterator->value.id() == NewickGrammar::typeID)
    {
        if(aTreeIterator->children.size() == 0)
		{
			std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
            //std::cout << "type " << label_name << std::endl;
			aNode->addType(label_name);
        }
		else
		{
			for(k=0; k < aTreeIterator->children.size(); ++k)
            	evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    	}
	}
    else if(aTreeIterator->value.id() == NewickGrammar::branchlenID)
    {
        if(aTreeIterator->children.size() == 0)
		{
			std::string blen(aTreeIterator->value.begin(), aTreeIterator->value.end());
            //std::cout << "branch_length " << atof(blen.c_str())<< std::endl;
			aNode->addLen(atof(blen.c_str()));
        }
		else
		{
			for(k=0; k < aTreeIterator->children.size(); ++k)
            	evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
		}
    }
    else if(aTreeIterator->value.id() == NewickGrammar::cblenID)
    {
        //std::cout << "colon_plus_len " << std::endl;
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    }
    else
    {
        std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
        std::cerr << "*** " << label_name << std::endl;
		throw FastCodeMLFatalNoMsg();
    }
}

PhyloTree::PhyloTree(int aVerboseLevel) : mVerboseLevel(aVerboseLevel)
{
}

PhyloTree::~PhyloTree()
{
	// Remove the tree
	mTreeRoot.clearNode();
}


void PhyloTree::loadTree(const char *aFilename)
{
    std::ifstream in(aFilename);
	if(!in)
	{
		std::cerr << "Cannot open " << aFilename << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	// Cleaning step. Should not be necessary after fixing the grammar...
	// Find the tree definition start
    std::string str;
    size_t p1 = std::string::npos;
    while(getline(in, str))
	{
        p1 = str.find_first_not_of(" \t\r");
        if(p1 != std::string::npos && str[p1] == '(') break;
	}
	in.close();
	if(p1 == std::string::npos)
	{
		std::cerr << "File " << aFilename << " is empty" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	// Find the tree definition end
	size_t p2 = str.rfind(";");
    if(p2 == std::string::npos)
	{
		std::cerr << "File " << aFilename << " is empty" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	// Pass the tree part to the parser
	loadTreeFromSting(str.substr(p1, p2-p1+1));
}


void PhyloTree::loadTreeFromSting(const std::string& aTreeAsString)
{
	// Clear the tree and initialize the parser
	mTreeRoot.clearNode();
    NewickGrammar tree;

	// Parse the tree
	tree_parse_info<> info = pt_parse(aTreeAsString.c_str(), tree);

	if(info.full)
	{
		mParsedPortion.clear(); // To signal no error
		evaluateTreeNode(info.trees.begin(), &mTreeRoot);

		// Fill the leaves
		mLeavesSpecies.clear();
		fillSpecies(&mTreeRoot);

		// Fill internal branches
		mInternalNodes.clear();
		fillInternalBranches(&mTreeRoot);

		// Dump the tree
		if(mVerboseLevel >= 3)
		{
			std::cerr << "Tree as read in PhyloTree" << std::endl;
			printTree(info.trees.begin(), 0);
			mTreeRoot.printFormatted(0);
			std::vector<TreeNode *>::const_iterator isp;
			for(isp=mLeavesSpecies.begin(); isp != mLeavesSpecies.end(); ++isp) std::cout << (*isp)->getLabel() << std::endl;
		}
	}
	else
	{
		mParsedPortion.assign(aTreeAsString.c_str(), info.stop);

		std::cerr << "Parsing failed" << std::endl;
		std::cerr << "----------------------" << std::endl;
		std::cerr << mParsedPortion << std::endl;
		std::cerr << "----------------------" << std::endl << std::endl;
		throw FastCodeMLFatalNoMsg();
	}
}


void PhyloTree::printFormattedTree(void) const
{
	mTreeRoot.printFormatted(0);
}


void PhyloTree::printNewickTree(TreeNode *aNode) const
{
	TreeNode *m;
	unsigned int idx;

	// Special case for the root
	if(!aNode)
	{
		std::cout << '(';
		for(idx=0; (m = mTreeRoot.getChild(idx)) != 0; ++idx)
		{
			if(idx > 0) std::cout << ',';
			printNewickTree(m);
		}
		std::cout << ')';
		mTreeRoot.printNode();
		std::cout << ";" << std::endl;
	}
	else if(aNode->isLeaf())
	{
		aNode->printNode();
	}
	else
	{
		std::cout << '(';
		for(idx=0; (m = aNode->getChild(idx)) != 0; ++idx)
		{
			if(idx > 0) std::cout << ',';
			printNewickTree(m);
		}
		std::cout << ')';
		aNode->printNode();
	}
}


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


void PhyloTree::getSpecies(std::vector<std::string>& aSpeciesList) const
{
	std::vector<TreeNode *>::const_iterator is;
	for(is=mLeavesSpecies.begin(); is != mLeavesSpecies.end(); ++is)
	{
		std::string label = (*is)->getLabel();
		aSpeciesList.push_back(label);
	}
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


unsigned int PhyloTree::cloneTree(ForestNode* aForestNode, unsigned int aTreeId, const TreeNode* aTreeNode, unsigned int aNodeId) const
{
	unsigned int id;

	// Start with the tree root
	if(aTreeNode == 0)
	{
		aTreeNode = &mTreeRoot;
		aNodeId   = UINT_MAX;
		id = 0;
	}
	else
	{
		id = aNodeId+1;
	}

	// Set the root values
	aForestNode->mNodeName = aTreeNode->getLabel();
	aForestNode->mBranchLength = aTreeNode->getLen();
	aForestNode->mNodeId = aNodeId;
	aForestNode->mOwnTree = aTreeId;

	// Set the internal branch identifier
	unsigned int int_id;
	unsigned int nn = mInternalNodes.size();
	for(int_id=0; int_id < nn; ++int_id) if(aTreeNode == mInternalNodes[int_id]) break;
	aForestNode->mInternalNodeId = (int_id < nn) ? int_id : UINT_MAX;

	// Recurse
	TreeNode *m;
	for(int idx=0; (m = aTreeNode->getChild(idx)) != 0; ++idx)
	{
		ForestNode* rn = new ForestNode;
		aForestNode->mChildrenList.push_back(rn);
		aForestNode->mChildSameTree.push_back(true);
		aForestNode->mOtherTreeProb.push_back(0);
		rn->mParent = aForestNode;
		id = cloneTree(rn, aTreeId, m, id);
	}

	return id;
}

