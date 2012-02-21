#include <iostream>
#include <fstream>
#include "Newick.h"
#include "NewickGrammar.h"
#include "Exceptions.h"
// Access to the Boost::Spirit parse tree
//
//		typedef boost::spirit::classic::tree_match<char const*> ParseTreeMatchType;
//		typedef ParseTreeMatchType::tree_iterator ParseTreeIteratorType;
// has been changed into:
//		typedef boost::spirit::classic::tree_match<char const*>::tree_iterator ParseTreeIteratorType;

void Newick::printTree(ParseTreeIteratorType const& aTreeIterator, unsigned int aIndent, unsigned int aIndentIncrement)
{
	// Indent the level
	for(unsigned int k=0; k < aIndent; ++k) std::cout << ' ';

	// Print the node
    if(aTreeIterator->value.id() == NewickGrammar::treeID)           std::cout << "TREE";
    else if(aTreeIterator->value.id() == NewickGrammar::nodelistID)  std::cout << "NODELST";
    else if(aTreeIterator->value.id() == NewickGrammar::subtreeID)   std::cout << "SUBTREE";
    else if(aTreeIterator->value.id() == NewickGrammar::fulllabelID) std::cout << "FLAB";
    else if(aTreeIterator->value.id() == NewickGrammar::labelID)     std::cout << "LABEL";
    else if(aTreeIterator->value.id() == NewickGrammar::typeID)      std::cout << "TYPE";
    else if(aTreeIterator->value.id() == NewickGrammar::branchlenID) std::cout << "BRANCHLEN";
    else if(aTreeIterator->value.id() == NewickGrammar::cblenID)     std::cout << "CBLEN";
    else															 std::cout << "????";

	// Print the corresponding label if any
	if(aTreeIterator->value.begin() != aTreeIterator->value.end())
	{
		std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
		std::cout << ' ' << label_name << std::endl;
	}
	else
	{
		std::cout << std::endl;
	}

	// Recurse on the children
	for(ParseTreeIteratorType ic=aTreeIterator->children.begin(); ic != aTreeIterator->children.end(); ++ic) printTree(ic, aIndent+aIndentIncrement, aIndentIncrement);
}


void Newick::evaluateTreeNode(ParseTreeIteratorType const& aTreeIterator, TreeNode *aNode)
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


void Newick::readFile(const char *aFilename)
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
	size_t p2 = str.rfind(';');
    if(p2 == std::string::npos)
	{
		std::cerr << "File " << aFilename << " is empty" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	// Pass the tree part to the parser
	loadTreeFromString(str.substr(p1, p2-p1+1));
}


void Newick::loadTreeFromString(const std::string& aTreeAsString)
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
			printTree(info.trees.begin());
			mTreeRoot.printFormatted(0);
			std::vector<TreeNode *>::const_iterator isp=mLeavesSpecies.begin();
			std::vector<TreeNode *>::const_iterator end=mLeavesSpecies.end();
			for(; isp != end; ++isp) std::cout << (*isp)->getLabel() << std::endl;
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


void Newick::printTreeUnformatted(std::ostream& aOut, TreeNode *aNode) const
{
	TreeNode *m;
	unsigned int idx;

	// Special case for the root
	if(!aNode)
	{
		aOut << '(';
		for(idx=0; (m = mTreeRoot.getChild(idx)) != 0; ++idx)
		{
			if(idx > 0) aOut << ',';
			printTreeUnformatted(aOut, m);
		}
		aOut << ')';
		mTreeRoot.printNode();
		aOut << ';' << std::endl;
	}
	else if(aNode->isLeaf())
	{
		aNode->printNode();
	}
	else
	{
		aOut << '(';
		for(idx=0; (m = aNode->getChild(idx)) != 0; ++idx)
		{
			if(idx > 0) aOut << ',';
			printTreeUnformatted(aOut, m);
		}
		aOut << ')';
		aNode->printNode();
	}
}

