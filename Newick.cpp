
#ifdef _MSC_VER
    #pragma warning(disable: 4512) // warning C4512: 'boost::spirit::classic::optional<S>' : assignment operator could not be generated
    #pragma warning(disable: 4503) // warning C4503: '...' : decorated name length exceeded, name was truncated
#endif

#include <iostream>
#include <fstream>
#include "Newick.h"
#include "NewickGrammar.h"
#include "Exceptions.h"
#include "VerbosityLevels.h"

// Access to the Boost::Spirit parse tree
//
//		typedef boost::spirit::classic::tree_match<char const*> ParseTreeMatchType;
//		typedef ParseTreeMatchType::tree_iterator ParseTreeIteratorType;
// has been changed into:
//		typedef boost::spirit::classic::tree_match<char const*>::tree_iterator ParseTreeIteratorType;

void Newick::printTree(ParseTreeIteratorType const& aTreeIterator, unsigned int aIndent, unsigned int aIndentIncrement)
{
	// Indent the level
	for(unsigned int k=0; k < aIndent; ++k) std::cerr << ' ';

	// Print the node (after getting the node id)
	boost::spirit::classic::parser_id id = aTreeIterator->value.id();

         if(id == NewickGrammar::treeID)      std::cerr << "TREE";
    else if(id == NewickGrammar::nodelistID)  std::cerr << "NODELST";
    else if(id == NewickGrammar::subtreeID)   std::cerr << "SUBTREE";
    else if(id == NewickGrammar::fulllabelID) std::cerr << "FLAB";
    else if(id == NewickGrammar::labelID)     std::cerr << "LABEL";
    else if(id == NewickGrammar::typeID)      std::cerr << "TYPE";
    else if(id == NewickGrammar::branchlenID) std::cerr << "BRANCHLEN";
    else if(id == NewickGrammar::cblenID)     std::cerr << "CBLEN";
    else									  std::cerr << "????";

	// Print the corresponding label if any
	if(aTreeIterator->value.begin() != aTreeIterator->value.end())
	{
		std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
		std::cerr << ' ' << label_name << std::endl;
	}
	else
	{
		std::cerr << std::endl;
	}

	// Recurse on the children
	for(ParseTreeIteratorType ic=aTreeIterator->children.begin(); ic != aTreeIterator->children.end(); ++ic) printTree(ic, aIndent+aIndentIncrement, aIndentIncrement);
}


void Newick::evaluateTreeNode(ParseTreeIteratorType const& aTreeIterator, TreeNode *aNode)
{
    unsigned int k;

	// Get the node id
	boost::spirit::classic::parser_id id = aTreeIterator->value.id();

    //if(aTreeIterator->value.id() == NewickGrammar::treeID)
    if(id == NewickGrammar::treeID)
    {
        //std::cerr << "tree (" << std::endl;
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
        //std::cerr << ")" << std::endl;
    }
    else if(id == NewickGrammar::nodelistID)
    {
        //std::cerr << "node_list {" << std::endl;
		//TreeNode* x = new TreeNode;
		//n->addChild(x);
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
        //std::cerr << "}" << std::endl;
    }
    else if(id == NewickGrammar::subtreeID)
    {
        //std::cerr << "subtree [" << std::endl;
		TreeNode* x = new TreeNode;
		aNode->addChild(x);
		x->addParent(aNode);
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, x);
        //std::cerr << "]" << std::endl;
    }
    else if(id == NewickGrammar::fulllabelID)
    {
        //std::cerr << "full_label " << std::endl;
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    }
    else if(id == NewickGrammar::labelID)
    {
        if(aTreeIterator->children.size() == 0)
		{
			std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
            //std::cerr << "label " << label_name << std::endl;
			aNode->addLabel(label_name);
        }
		else
		{
			for(k=0; k < aTreeIterator->children.size(); ++k)
            	evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    	}
	}
    else if(id == NewickGrammar::typeID)
    {
        if(aTreeIterator->children.size() == 0)
		{
			std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
            //std::cerr << "type " << label_name << std::endl;
			aNode->addType(label_name);
        }
		else
		{
			for(k=0; k < aTreeIterator->children.size(); ++k)
            	evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    	}
	}
    else if(id == NewickGrammar::branchlenID)
    {
        if(aTreeIterator->children.size() == 0)
		{
			std::string blen(aTreeIterator->value.begin(), aTreeIterator->value.end());
            //std::cerr << "branch_length " << atof(blen.c_str())<< std::endl;
			aNode->addLen(atof(blen.c_str()));
        }
		else
		{
			for(k=0; k < aTreeIterator->children.size(); ++k)
            	evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
		}
    }
    else if(id == NewickGrammar::cblenID)
    {
        //std::cerr << "colon_plus_len " << std::endl;
        for(k=0; k < aTreeIterator->children.size(); ++k)
            evaluateTreeNode(aTreeIterator->children.begin()+k, aNode);
    }
    else
    {
        std::string label_name(aTreeIterator->value.begin(), aTreeIterator->value.end());
		std::ostringstream o;
        o << "*** " << label_name << std::endl;
		throw FastCodeMLFatal(o);
    }
}


void Newick::readFile(const char *aFilename)
{
	// Open the tree file
    std::ifstream in(aFilename);
	if(!in)
	{
		std::ostringstream o;
		o << "Cannot open " << aFilename << std::endl;
		throw FastCodeMLFatal(o);
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
		std::ostringstream o;
		o << "File " << aFilename << " is empty" << std::endl;
		throw FastCodeMLFatal(o);
	}

	// Find the tree definition end
	size_t p2 = str.rfind(';');
    if(p2 == std::string::npos)
	{
		std::ostringstream o;
		o << "File " << aFilename << " is empty" << std::endl;
		throw FastCodeMLFatal(o);
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
		if(mVerboseLevel >= VERBOSE_MORE_DEBUG)
		{
			std::cerr << "Tree as read in PhyloTree" << std::endl;
			printTree(info.trees.begin());
			mTreeRoot.printFormatted();
			std::vector<TreeNode *>::const_iterator isp(mLeavesSpecies.begin());
			const std::vector<TreeNode *>::const_iterator end(mLeavesSpecies.end());
			for(; isp != end; ++isp) std::cerr << (*isp)->getLabel() << std::endl;
		}
	}
	else
	{
		mParsedPortion.assign(aTreeAsString.c_str(), info.stop);

		std::ostringstream o;
		o << "Parsing failed\n"
		  << "----------------------\n"
		  << mParsedPortion
		  << "\n----------------------" << std::endl << std::endl;
		throw FastCodeMLFatal(o);
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
		for(idx=0; (m = mTreeRoot.getChild(idx)) != NULL; ++idx)
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
		for(idx=0; (m = aNode->getChild(idx)) != NULL; ++idx)
		{
			if(idx > 0) aOut << ',';
			printTreeUnformatted(aOut, m);
		}
		aOut << ')';
		aNode->printNode();
	}
}

