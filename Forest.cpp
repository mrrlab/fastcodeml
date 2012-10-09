
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <boost/dynamic_bitset.hpp>
#include "Forest.h"
#include "ForestNode.h"
#include "Exceptions.h"
#include "MathSupport.h"
#include "MatrixSize.h"
#include "CompilerHints.h"
#include "VerbosityLevels.h"

// Initialize the mask table so the index corresponds to the bit position
const unsigned char ForestNode::mMaskTable[MAX_NUM_CHILDREN] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};


void Forest::loadTreeAndGenes(const PhyloTree& aTree, const Genes& aGenes, CodonFrequencies::CodonFrequencyModel aCodonFrequencyModel)
{
	// Collect global data that refers to the tree and that should not be duplicated on each tree of the forest
	aTree.collectGlobalTreeData(mNodeNames, mBranchLengths, &mMarkedInternalBranch);

	// Number of branches of one tree
	mNumBranches = aTree.getNumBranches();

	// Count the number of unique sites and get the multiplicity of each of them
	mNumSites = aGenes.getNumSites();
	const std::vector<unsigned int>& mult = aGenes.getSiteMultiplicity();

	// Initialize the count of codon types
	std::vector<unsigned int> codon_count(N, 0);

	// Initialize the array of all probability vectors
	mProbs.assign(mNumSites*(mNumBranches+1)*Nt*VECTOR_SLOT, 0.0);
#ifdef NEW_LIKELIHOOD
	mProbsOut.assign(mNumSites*(mNumBranches+1)*Nt*VECTOR_SLOT, 0.0);
#endif

	// Count of tree's leaves
	size_t num_leaves = 0;

	// Allocate the list of pointers to leaves
	std::vector<ForestNode*> leaves;

	// Clone tree inside the forest
	mRoots.resize(mNumSites);
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Create a copy of the tree
		aTree.cloneTree(&mRoots[site], site, mNumSites, mProbs);

		// Create a list of pointers to leaves
		leaves.clear();
		mRoots[site].pushLeaf(leaves);
		num_leaves = leaves.size();

		// Add codon code to leaves
		std::vector<ForestNode*>::const_iterator il=leaves.begin();
		const std::vector<ForestNode*>::const_iterator end=leaves.end();
		for(; il != end; ++il)
		{
			// Node id (adjusted so root is 0)
			unsigned int node = (*il)->mBranchId+1;

			// Get the codon index and add it to the node signature
			int codon = aGenes.getCodonIdx(mNodeNames[node], site);
			if(codon < 0)
			{
				std::ostringstream o;
				o << "Invalid codon at site " << site << " node: " << node << " " << mNodeNames[node] << std::endl;
				throw FastCodeMLFatal(o);
			}
			(*il)->mPreprocessingSupport->mSubtreeCodonsSignature.push_back(codon);

			// Record the codon number on the leaf node (to use a simpler transition for leaves)
			(*il)->mLeafCodon = static_cast<short>(codon);

			// Set leaves probability vector (Nt copies)
#ifdef NEW_LIKELIHOOD
			for(int set=0; set < Nt; ++set) mProbs[VECTOR_SLOT*(node*(Nt*mNumSites)+set*mNumSites+site)+codon] = 1.0
#ifdef USE_GLOBAL_SCALING	
				*GLOBAL_SCALING_FACTOR
#endif
				; // The rest already zeroed by assign()
#else
			for(int set=0; set < Nt; ++set) (*il)->mProb[set][codon] = 1.0
#ifdef USE_GLOBAL_SCALING	
				*GLOBAL_SCALING_FACTOR
#endif
				; // The rest already zeroed by assign()
#endif
			// Count codons
			codon_count[codon] += mult[site];
		}

		// Combine the subtrees signatures going up to the root
		mRoots[site].gatherCodons();
	}

	// Set the number of internal branches
	mNumInternalBranches = mNumBranches - num_leaves;

	// Set the site multeplicity
	mSiteMultiplicity.assign(mult.begin(), mult.end());

	// Set the codon frequencies and related values needed for the eigensolver
	CodonFrequencies* cf = CodonFrequencies::getInstance();
	cf->setCodonFrequencies(codon_count, aCodonFrequencyModel);
	mCodonFreq     = cf->getCodonFrequencies();
	mInvCodonFreq  = cf->getInvCodonFrequencies();
	mInv2CodonFreq = cf->getCodonFreqInv2();

	// Set the mapping from internal branch number to branch number (the last tree has no pruned subtrees)
	std::map<unsigned int, unsigned int> map_internal_to_branchID;
	mapInternalToBranchIdWalker(&mRoots[mNumSites-1], map_internal_to_branchID);

	// Transform the map into a table (better for performance)
	mTableInternalToBranchID.resize(map_internal_to_branchID.size());
	std::map<unsigned int, unsigned int>::const_iterator im=map_internal_to_branchID.begin();
	const std::map<unsigned int, unsigned int>::const_iterator end=map_internal_to_branchID.end();
	for(; im != end; ++im) mTableInternalToBranchID[im->first] = im->second;

#ifdef NEW_LIKELIHOOD
	postLoad();
#endif
}

#ifdef NEW_LIKELIHOOD
void Forest::postLoad(void)
{
    // Prepare the list of node id's by level
    std::vector<ForestNode*> next_level;
    std::vector<ForestNode*> curr_level;
    std::vector<ForestNode*> level_nodes;

    // First level is the root (but it is not added because no processing is done on it)
    //level_nodes.push_back(&mRoots[mNumSites-1]);
    mNodesByLevel.clear();
    //mNodesByLevel.push_back(level_nodes);
    curr_level.push_back(&mRoots[mNumSites-1]);

    // Continue with all levels till reaching the leaves
    for(;; curr_level = next_level)
    {
        // Empty temporary arrays
        next_level.clear();
        level_nodes.clear();

        // Put in a list all the children of the current level nodes
        std::vector<ForestNode*>::const_iterator il=curr_level.begin();
        const std::vector<ForestNode*>::const_iterator end=curr_level.end();
        for(; il != end; ++il)
        {
            if(!(*il)->mChildrenList.empty()) next_level.insert(next_level.end(), (*il)->mChildrenList.begin(), (*il)->mChildrenList.end());
        }

        // No children, the last level was the leaves level
        if(next_level.empty()) break;

        // Add the list of node pointers of this level
        for(il=next_level.begin(); il != next_level.end(); ++il)
        {
            level_nodes.push_back(*il);
        }
        mNodesByLevel.push_back(level_nodes);
    }

#if 0
	// Show the tree before balancing
	std::cerr << std::endl;
	std::vector< std::vector<ForestNode*> >::const_reverse_iterator rinbl;
	unsigned int level = 0;
	for(rinbl=mNodesByLevel.rbegin(); rinbl != mNodesByLevel.rend(); ++rinbl,++level)
	{
		std::cerr << "Level " << level << ": ";
		std::vector<ForestNode*>::const_iterator ifn;
		for(ifn=rinbl->begin(); ifn != rinbl->end(); ++ifn)
		{
			std::cerr << (*ifn)->mBranchId << ((*ifn)->mChildrenList.empty() ? "* " : "  ");
		}
		std::cerr << std::endl;
	}
#endif

	// Try to balance the tree (ie move leaves to fill underfull levels)
	for(bool found=true; found;)
	{
		// Find the level with the maximum number of leaves
		size_t max_len   = 0;
		unsigned int max_level = 0;
		unsigned int max_leaf  = 0;
		unsigned int level = 0;
		std::vector< std::vector<ForestNode*> >::iterator inbl=mNodesByLevel.begin();
		const std::vector< std::vector<ForestNode*> >::iterator end=mNodesByLevel.end();
		for(level=0; inbl != end; ++inbl,++level)
		{
			unsigned int num_leaves = 0;
			unsigned int leaf = 0, i=0;
			std::vector<ForestNode*>::const_iterator ifn=inbl->begin();
			const std::vector<ForestNode*>::const_iterator end=inbl->end();
			for(; ifn != end; ++ifn,++i)
			{
				if((*ifn)->mChildrenList.empty()) {++num_leaves; leaf = i;}
			}
			if(num_leaves == 0) continue;

			size_t len = inbl->size();
			if(len > max_len) {max_len = len; max_level = level; max_leaf = leaf;}

		}

		// Find the first level that can inglobate the leave from level 'max_level' and index 'max_leaf'
		found = false;
		for(inbl=mNodesByLevel.begin()+max_level+1,level=max_level+1; inbl != mNodesByLevel.end(); ++inbl,++level)
		{
			const size_t len = inbl->size();
			if(len < max_len-1)
			{
				mNodesByLevel[level].push_back(mNodesByLevel[max_level][max_leaf]);
				mNodesByLevel[max_level].erase(mNodesByLevel[max_level].begin()+max_leaf);
				found = true;
				break;
			}
		}
	}


#if 0
	// Show the tree after balancing
	std::cerr << std::endl;
	for(rinbl=mNodesByLevel.rbegin(),level=0; rinbl != mNodesByLevel.rend(); ++rinbl,++level)
	{
		std::cerr << "Level " << level << ": ";
		std::vector<ForestNode*>::const_iterator ifn;
		for(ifn=rinbl->begin(); ifn != rinbl->end(); ++ifn)
		{
			std::cerr << (*ifn)->mBranchId << ((*ifn)->mChildrenList.empty() ? "* " : "  ");
		}
		std::cerr << std::endl;
	}
#endif

	// Record the dependencies between branches
	mFatVectorTransform.setBranchDependencies(mNodesByLevel);
}
#endif

#if 0
void Forest::reduceSubtrees(void)
{
	// Setup dependency vectors
	std::vector<unsigned int> empty_vector;
	mTreeDependencies.resize(mNumSites, empty_vector);
	mTreeRevDependencies.resize(mNumSites, empty_vector);

	// Try to merge equal subtrees
	// Trees at the beginning of the forest point to trees ahead
	// (this way a delete does not choke with pointers pointing to freed memory)
	int ns = static_cast<int>(mNumSites);
	for(int i=ns-1; i > 0; --i)
	{
		for(int j=i-1; j >= 0; --j)
		{
			reduceSubtreesWalker(&mRoots[i], &mRoots[j]);
		}
	}
}
#endif

bool Forest::reduceSubtrees(unsigned int aBlocks)
{
	// Setup dependency vectors
	std::vector<unsigned int> empty_vector;
	mTreeDependencies.resize(mNumSites, empty_vector);
	mTreeRevDependencies.resize(mNumSites, empty_vector);

	// Make integer the number of sites and the number of blocks (otherwise the countdown does not work)
	const int ns = static_cast<int>(mNumSites);
	const int num_blocks = static_cast<int>(aBlocks);

	if(num_blocks < 1 || num_blocks >= ns/2)
	{
		// No reduction at all
		return false;
	}
	if(num_blocks == 1)
	{
		// Try to merge equal subtrees (this is the usual global check and reduction)
		// Trees at the beginning of the forest point to trees ahead
		// (this way a delete does not choke with pointers pointing to freed memory)
		for(int i=ns-1; i > 0; --i)
		{
			for(int j=i-1; j >= 0; --j)
			{
				reduceSubtreesWalker(&mRoots[i], &mRoots[j]);
			}
		}
	}
	else
	{
		const int bsize = ns / num_blocks;
		int start = -1;
		int end   =  0;

		for(int block=0; block < num_blocks; ++block)
		{
			end = start+1;
			start += bsize;
			if((ns-1-start) < bsize) start = ns-1;

			for(int i=start; i > end; --i)
			{
				for(int j=i-1; j >= end; --j)
				{
					reduceSubtreesWalker(&mRoots[i], &mRoots[j]);
				}
			}
		}
	}

	return true;
}


void Forest::reduceSubtreesWalker(ForestNode* aNode, ForestNode* aNodeDependent)
{
	unsigned int i;
	const unsigned int nc = aNode->mChildrenCount;
	for(i=0; i < nc; ++i)
	{
		// If one of the two has been already reduced, do nothing
		if(!(aNode->isSameTree(i)) || !(aNodeDependent->isSameTree(i))) continue;

		// Check if same subtree
		if(aNode->mChildrenList[i]->mPreprocessingSupport->mSubtreeCodonsSignature == aNodeDependent->mChildrenList[i]->mPreprocessingSupport->mSubtreeCodonsSignature)
		{
			delete aNodeDependent->mChildrenList[i];
			aNodeDependent->mChildrenList[i] = aNode->mChildrenList[i];
			aNodeDependent->markNotSameTree(i);

			// Record dependencies
			mTreeDependencies[aNodeDependent->mOwnTree].push_back(aNode->mOwnTree);		// [tj] can be done after: t1 t2 t3
			mTreeRevDependencies[aNode->mOwnTree].push_back(aNodeDependent->mOwnTree);	// [tj] should be ready before: t1 t2 t3
		}
	}

	// Recurse
	for(i=0; i < nc; ++i)
	{
		// If one of the two has been already reduced, do nothing
		if(!(aNode->isSameTree(i)) || !(aNodeDependent->isSameTree(i))) continue;

		reduceSubtreesWalker(aNode->mChildrenList[i], aNodeDependent->mChildrenList[i]);
	}
}


void Forest::cleanReductionWorkingData(ForestNode* aNode)
{
	if(!aNode)
	{
		// Invoke on all the trees in the forest
		for(size_t i=0; i < mNumSites; ++i) cleanReductionWorkingData(&mRoots[i]);
	}
	else
	{
		// Clean myself
		delete aNode->mPreprocessingSupport;
		aNode->mPreprocessingSupport = 0;

		// Clean the children
		const unsigned int nc = aNode->mChildrenCount;
		for(unsigned int i = 0; i < nc; ++i)
		{
			if(aNode->isSameTree(i)) cleanReductionWorkingData(aNode->mChildrenList[i]);
		}
	}
}

 std::ostream& operator<< (std::ostream& aOut, const Forest& aForest)
{
	size_t i;

	aOut << std::endl;
	aOut << "Num branches:       " << std::setw(7) << aForest.mNumBranches << std::endl;
	aOut << "Unique sites:       " << std::setw(7) << aForest.mNumSites << std::endl;
	aOut << "Total branches:     " << std::setw(7) << aForest.mNumBranches*aForest.mNumSites << std::endl;

	// Count total branches on the reduced forest
	unsigned int cnt = 0;
	unsigned int cntAggressive = 0;
	for(i=0; i < aForest.mNumSites; ++i)
	{
		const ForestNode& n = aForest.mRoots[i];
		cnt += n.countBranches();
		cntAggressive += n.countBranches(true);
	}
	aOut << "Reduced branches:   " << std::fixed << std::setw(7) << cnt << std::setw(8) << std::setprecision(2) << static_cast<double>(cnt*100.)/static_cast<double>(aForest.mNumBranches*aForest.mNumSites) << '%' << std::endl;
	aOut << "Aggressive reduct.: " << std::fixed << std::setw(7) << cntAggressive << std::setw(8) << std::setprecision(2) << static_cast<double>(cntAggressive*100.)/static_cast<double>(aForest.mNumBranches*aForest.mNumSites) << '%' << std::endl;
	aOut << std::endl;

	// Print forest
	if(aForest.mVerbose >= VERBOSE_DSTRUCT_DUMP)
	{
		for(i=0; i < aForest.mNumSites; ++i)
		{
			aOut << "=== Site " << i << " ===" << std::endl;
			aForest.mRoots[i].print(aForest.getNodeNames(), aOut);
			aOut << std::endl;
		}
	}

	return aOut;
}

void Forest::getEffortPerSite(std::vector<unsigned int>& aEfforts, unsigned int aCostAtLeaf, unsigned int aCostIntern, unsigned int aCostPtr) const
{
	// Initialize effort array
	aEfforts.clear();
	aEfforts.reserve(mNumSites);

	// Get the total cost per site
	for(size_t i=0; i < mNumSites; ++i)
	{
		unsigned int total_cost = mRoots[i].getCost(aCostAtLeaf, aCostIntern, aCostPtr);
		aEfforts.push_back(total_cost);
	}
}

void Forest::setTimesFromLengths(std::vector<double>& aTimes, const ForestNode* aNode) const
{
	if(!aNode) aNode = &mRoots[mNumSites-1];
	else
	{
		const unsigned int id = aNode->mBranchId;
		aTimes[id] = mBranchLengths[id+1];
	}

	std::vector<ForestNode *>::const_iterator ifn=aNode->mChildrenList.begin();
	for(; ifn != aNode->mChildrenList.end(); ++ifn)
	{
		setTimesFromLengths(aTimes, *ifn);
	}
}


void Forest::setLengthsFromTimes(const std::vector<double>& aTimes, ForestNode* aNode)
{
	std::vector<ForestNode *>::const_iterator ifnp;

	// Get all forest connections
	if(!aNode)
	{
		std::vector<ForestNode>::iterator ifn=mRoots.begin();
		for(; ifn != mRoots.end(); ++ifn)
		{
			for(ifnp=ifn->mChildrenList.begin(); ifnp != ifn->mChildrenList.end(); ++ifnp)
			{
				setLengthsFromTimes(aTimes, *ifnp);
			}
		}
	}
	else
	{
		const unsigned int idx = aNode->mBranchId+1;
		mBranchLengths[idx] = aTimes[aNode->mBranchId];

		for(ifnp=aNode->mChildrenList.begin(); ifnp != aNode->mChildrenList.end(); ++ifnp)
		{
			setLengthsFromTimes(aTimes, *ifnp);
		}
	}
}


#ifdef NEW_LIKELIHOOD
// Compute likelihood with the new "Long Vector" approach
//
void Forest::computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int /*aHyp*/)
{
	// Initialize variables
    const unsigned int num_sets = aSet.size();
    //aLikelihoods.assign(num_sets*mNumSites, 1.0);
	//aLikelihoods.resize(num_sets*mNumSites);

	// For each level of the tree (except the root)
	unsigned int level=0;
    std::vector< std::vector<ForestNode*> >::reverse_iterator inbl;
    for(inbl=mNodesByLevel.rbegin(); inbl != mNodesByLevel.rend(); ++inbl,++level)
    {
		const int num_sites = static_cast<int>(inbl->size());
        const int len       = num_sites*num_sets;

#ifdef _MSC_VER
        #pragma omp parallel for default(none) shared(aSet, len, inbl, num_sets, num_sites, level) schedule(guided)
#else
        #pragma omp parallel for default(shared) schedule(guided)
#endif
        for(int i=0; i < len; ++i)
        {
            // Compute probability vector along this branch (for the given set) (reordered to give a 2% speedup)
			const unsigned int set_idx  = i / num_sites;
			const unsigned int site_idx = i - set_idx * num_sites; // Was: unsigned int set_idx = i % num_sets;
            const unsigned int branch   = ((*inbl)[site_idx])->mBranchId;
			const size_t       start    = VECTOR_SLOT*(mNumSites*Nt*(branch+1)+mNumSites*set_idx+mFatVectorTransform.getLowerIndex(branch));

			// For each branch, except the root, compute the transition
            aSet.doTransition(set_idx,
							  branch,
							  static_cast<int>(mFatVectorTransform.getCount(branch)),
							  &mProbs[start],
							  &mProbsOut[start]);
        }

		// Combine the results to have the input for the next round
		mFatVectorTransform.postCompact(mProbsOut, mProbs, level, num_sets);
    }

	// Compute the final likelyhood
	const int num_sites = static_cast<int>(mNumSites);
    const int len       = num_sites*num_sets;

#ifdef _MSC_VER
    #pragma omp parallel for default(none) shared(len, num_sites, aLikelihoods) schedule(guided)
#else
    #pragma omp parallel for default(shared) schedule(guided)
#endif
    for(int i=0; i < len; ++i)
    {
		const unsigned int set_idx = i / num_sites;
		const unsigned int site    = i - set_idx * num_sites; // Was: unsigned int site_idx = i % num_sites;
		const size_t       start   = VECTOR_SLOT*(set_idx*mNumSites+site);

		// Take the result from branch 0 (the root)
        aLikelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, &mProbs[start]);
    }
}
#endif

void Forest::mapInternalToBranchIdWalker(const ForestNode* aNode, std::map<unsigned int, unsigned int>& aMapInternalToBranchID)
{
	const unsigned int nc = aNode->mChildrenCount;
	for(unsigned int i=0; i < nc; ++i)
	{
		ForestNode *m = aNode->mChildrenList[i];

		if(m->mInternalNodeId != UINT_MAX) aMapInternalToBranchID[m->mInternalNodeId] = m->mBranchId;

		mapInternalToBranchIdWalker(m, aMapInternalToBranchID);
	}
}


#ifndef NEW_LIKELIHOOD
void Forest::addAggressiveReduction(ForestNode* aNode)
{
	if(aNode)
	{
		const unsigned int nc = aNode->mChildrenCount;
		for(unsigned int i=0; i < nc; ++i)
		{
			ForestNode *m = aNode->mChildrenList[i];

			if(aNode->isSameTree(i))
			{
				addAggressiveReduction(m);
			}
			else
			{
				ForestNode *other = m->mParent;

				// Add the array on the other side
				if(!other->mOtherTreeProb[i]) other->mOtherTreeProb[i] = static_cast<double*>(alignedMalloc(VECTOR_SLOT*Nt*sizeof(double), CACHE_LINE_ALIGN));

				// Add the pointer here
				aNode->mOtherTreeProb[i] = other->mOtherTreeProb[i];
			}
		}
	}
	else
	{
		for(size_t i=0; i < mNumSites; ++i)
		{
			addAggressiveReduction(&mRoots[i]);
		}
	}
}
#endif


#ifdef NON_RECURSIVE_VISIT
	
void Forest::prepareNonRecursiveVisit(void)
{
	// Clean the list for non-recursive visit to the trees. Clear also the list of respective parents
	mVisitTree.clear();
	mVisitTreeParents.clear();

	// Visit each site tree to collect threading pointers
	const unsigned int ns = mNumSites;
	for(unsigned int i=0; i < ns; ++i)
	{
		std::vector<ForestNode*> visit_list;
		std::vector<ForestNode*> parent_list;

		prepareNonRecursiveVisitWalker(&mRoots[i], 0, i, visit_list, parent_list);

		mVisitTree.push_back(visit_list);
		mVisitTreeParents.push_back(parent_list);
	}
}

void Forest::prepareNonRecursiveVisitWalker(ForestNode* aNode, ForestNode* aParentNode, unsigned int aSite, std::vector<ForestNode*>& aVisitList, std::vector<ForestNode*>& aParentList)
{
	// Collect the number of children of the current node
	const unsigned int nc = aNode->mChildrenCount;

	// Check if it is a leaf or a node on another tree
	if(nc != 0 && aNode->mOwnTree == aSite)
	{
		// Internal nodes
		bool first = true;
		for(unsigned int i=0; i < nc; ++i)
		{
			ForestNode* n = aNode->mChildrenList[i];

			// Mark the first child
			n->mFirstChild = first;
			first = false;

			// Save the child position in the parent node
			n->mChildIdx = i;

			// Visit the subtree starting here
			prepareNonRecursiveVisitWalker(n, aNode, aSite, aVisitList, aParentList);
		}
	}

	// If it is not the root node
	// Store the nodes in the visit order except the root that should not be visited
	// Store also the respective parent node
	if(aParentNode)
	{	
		aVisitList.push_back(aNode);
		aParentList.push_back(aParentNode);
	}
}



void Forest::computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int aHyp)
{
	// To speedup the inner OpenMP parallel loop, this is precomputed
	// so in the call to computeLikelihoodsWalkerTC &mRoots[site] becomes tmp_roots+site
	ForestNode* tmp_roots = &mRoots[0];

	ListDependencies::iterator ivs=mDependenciesClassesAndTrees[aHyp].begin();
	const ListDependencies::iterator end=mDependenciesClassesAndTrees[aHyp].end();
	for(; ivs != end; ++ivs)
	{
		// Things that do not change in the parallel loop
		const int len = static_cast<int>(ivs->size());
		const unsigned int* tmp_ivs = &(*ivs)[0];

#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(aSet, len, tmp_ivs, tmp_roots, aLikelihoods) schedule(static)
#else
		#pragma omp parallel for default(shared) schedule(static)
#endif
		for(int i=0; i < len; ++i)
		{
			// Compute likelihood array at the root of one tree (the access order is the fastest)
			const unsigned int tmp     = tmp_ivs[i];
			const unsigned int site    = getSiteNum(tmp);
			const unsigned int set_idx = getSetNum(tmp);

			computeLikelihoodsWalkerNR(aSet, set_idx, site);

			aLikelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, tmp_roots[site].mProb[set_idx]);
		}
	}
}

void Forest::computeLikelihoodsWalkerNR(const ProbabilityMatrixSet& aSet, unsigned int aSetIdx, unsigned int aSiteIdx)
{
	unsigned int nv = mVisitTree[aSiteIdx].size();
	for(unsigned int j=0; j < nv; ++j)
	{
		ForestNode* n = mVisitTree[aSiteIdx][j];

		const unsigned int branch_id	= n->mBranchId;
		double* node_prob				= n->mProb[aSetIdx];
		double* other_tree_prob			= n->mParent->mOtherTreeProb[n->mChildIdx];

		if(n->mOwnTree == aSiteIdx)
		{
			double* res_prob = n->mParent->mProb[aSetIdx];

			if(n->mFirstChild)
			{
				aSet.doTransition(aSetIdx, branch_id, node_prob, res_prob);
				if(other_tree_prob) memcpy(other_tree_prob+VECTOR_SLOT*aSetIdx, res_prob, N*sizeof(double));
			}
			else
			{
				double ALIGN64 temp[N];
				double* x = other_tree_prob ? other_tree_prob+VECTOR_SLOT*aSetIdx : temp;
				aSet.doTransition(aSetIdx, branch_id, node_prob, x);
				elementWiseMult(res_prob, x);
			}
		}
		else
		{
			double* res_prob = mVisitTreeParents[aSiteIdx][j]->mProb[aSetIdx];

			if(n->mFirstChild)
			{
				memcpy(res_prob, other_tree_prob+VECTOR_SLOT*aSetIdx, N*sizeof(double));
			}
			else
			{
				elementWiseMult(res_prob, other_tree_prob+VECTOR_SLOT*aSetIdx);
			}
		}
	}
}
#endif

#ifdef NEW_LIKELIHOOD
void Forest::prepareNewReduction(ForestNode* aNode)
{
	if(aNode)
	{
		const unsigned int nc = aNode->mChildrenCount;
		for(unsigned int i=0; i < nc; ++i)
		{
			ForestNode* n = aNode->mChildrenList[i];

			if(aNode->isSameTree(i))
			{
				mFatVectorTransform.setNodeExists(n->mBranchId, aNode->mOwnTree);

				prepareNewReduction(n);
			}
			else
			{
				mFatVectorTransform.setNodeReuses(n->mBranchId, aNode->mOwnTree, n->mOwnTree);
			}
		}
	}
	else
	{
		// Initialize the intermediate list used to compute the list of commands
		mFatVectorTransform.initNodeStatus(mNumBranches, mNumSites);

		// Visit each site tree
		size_t ns = mNumSites;
		for(size_t i=0; i < ns; ++i) prepareNewReduction(&mRoots[i]);

		// Print few statistics on the transformation
		//mFatVectorTransform.printCountGoodElements();
		//mFatVectorTransform.printBranchVisitSequence();
		//mFatVectorTransform.printNodeStatus();

		// Compact the matrix (this creates the lists of operations needed)
		mFatVectorTransform.compactMatrix();

		// Print the commands
		//mFatVectorTransform.printCommands();

		// Do the initial move
		//crc(mProbs, mNumSites);
		mFatVectorTransform.preCompactLeaves(mProbs);
		//crc(mProbs, mNumSites);
	}
}


void Forest::prepareNewReductionNoReuse(void)
{
	mFatVectorTransform.initNodeStatusMinimal(mNumBranches, mNumSites);
}
#endif

#ifdef CHECK_ALGO
unsigned int Forest::checkForest(bool aCheckId, const ForestNode* aNode, unsigned int aSite, unsigned int aNodeId) const
{
	if(!aNode)
	{
		std::cerr << "========== Start forest check" << std::endl;
		size_t nsites = mRoots.size();

		// Visit each site tree
		for(size_t i=0; i < nsites; ++i)
		{
			checkForest(aCheckId, &mRoots[i], (unsigned int)i, UINT_MAX);
		}
		std::cerr << "==========   End forest check" << std::endl;
		return UINT_MAX;
	}
	else
	{
		// Check the node
		if(aNode->mOwnTree != aSite) std::cerr << "[site: " << std::setw(5) << aSite << "] mOwnTree mismatch " << aNode->mOwnTree << " should be: " << aSite << std::endl;

		// Check node ID
		if(aCheckId && aNode->mBranchId != aNodeId) std::cerr << "[site: " << std::setw(5) << aSite << "] aNodeId mismatch " << aNode->mBranchId << " should be: " << aNodeId << std::endl;

		unsigned int id = aNodeId+1;

		unsigned int idx;
		std::vector<ForestNode *>::const_iterator icl=aNode->mChildrenList.begin();
		const std::vector<ForestNode *>::const_iterator end=aNode->mChildrenList.end();
		for(idx=0; icl != end; ++icl,++idx)
		{
			if(aNode->isSameTree(idx))
			{
				// Check parent pointer
				if((*icl)->mParent != aNode)
					std::cerr << "[site: " << std::setw(5) << aSite << "] mParent mismatch " << aNode->mParent << " should be: " << aNode << std::endl;

				id = checkForest(aCheckId, *icl, aSite, id);
			}
			else if((*icl)->mOwnTree <= aNode->mOwnTree)
			{
				std::cerr << "[site: " << std::setw(5) << aSite << "] Strange intersite pointer. This: " << aNode->mOwnTree << " children: " << (*icl)->mOwnTree << std::endl;
			}
			else
			{
				//id = checkForest(aCheckId, *icl, aSite, id);
			}
		}
		return id;
	}
}

void chkleaves(CacheAlignedDoubleVector& p, int slot, const char *filename)
{
	std::ofstream out(filename, std::ios_base::trunc | std::ios_base::out);
	if(!out.good()) return;

	int nslots = p.size()/slot;

	for(int i=0; i < nslots; ++i)
	{
		for(int j=0; j < N; ++j)
		{
			if(p[i*slot+j] > 0.1)
			{
				out << std::setw(5) << i << std::setw(3) << j << std::endl;
				break;
			}
		}
	}
	out.close();
}

void crc(const std::vector<double>& v, unsigned int nsites)
{
	union { double value; unsigned char bytes[sizeof(double)]; } data;
	int  c1 = 52845; 
	int  c2 = 22719;
	unsigned int nv = v.size();

	for(unsigned int n=0; n < (nv/(VECTOR_SLOT*nsites*Nt)); ++n)
	{
		std::cerr << std::setw(2) << n << ": ";
		for(unsigned int s=0; s < nsites; ++s)
		{
			long sum = 0;
			int r = 55665;
			for(int j=0; j < N; ++j)
			{
				data.value = v[n*(Nt*nsites*VECTOR_SLOT)+s*(VECTOR_SLOT)+j];

				for(unsigned int k = 0; k < sizeof(double); ++k)
				{
					unsigned char cipher = (data.bytes[k] ^ (r >> 8));
					r = (cipher + r) * c1 + c2;
					sum += cipher;
				}
			}
			if(sum == 0xfb9d) std::cerr << "  -  ";
			else  			  std::cerr << std::hex << sum << ' ';
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
}
#endif


#if !defined(NON_RECURSIVE_VISIT) && !defined(NEW_LIKELIHOOD)

void Forest::computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, const ListDependencies& aDependencies)
{
	// To speedup the inner OpenMP parallel loop, this is precomputed
	// so in the call to computeLikelihoodsWalkerTC &mRoots[site] becomes tmp_roots+site
	const ForestNode* tmp_roots = &mRoots[0];
	double* likelihoods = &aLikelihoods[0];

	ListDependencies::const_iterator ivs=aDependencies.begin();
	const ListDependencies::const_iterator end=aDependencies.end();
	for(; ivs != end; ++ivs)
	{
		// Things that do not change in the parallel loop
		const int len = static_cast<int>(ivs->size());
		const unsigned int* tmp_ivs = &(*ivs)[0];

#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(aSet, len, tmp_ivs, tmp_roots, likelihoods) schedule(static)
#else
		#pragma omp parallel for default(shared)
#endif
		for(int i=0; i < len; ++i)
		{
#ifndef _MSC_VER
			#pragma omp task untied
#endif
			{
				// Compute likelihood array at the root of one tree (the access order is the fastest)
				const unsigned int tmp     = tmp_ivs[i];
				const unsigned int site    = TreeAndSetsDependencies::getSiteNum(tmp);
				const unsigned int set_idx = TreeAndSetsDependencies::getSetNum(tmp);

				const double* g = computeLikelihoodsWalkerTC(tmp_roots+site, aSet, set_idx);

				likelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, g);
			}
		}
	}
}

double* Forest::computeLikelihoodsWalkerTC(const ForestNode* aNode, const ProbabilityMatrixSet& aSet, unsigned int aSetIdx)
{
	double* anode_prob    = aNode->mProb[aSetIdx];
	const unsigned int nc = aNode->mChildrenCount;

	// Shortcut (on the leaves return immediately the probability vector)
	if(nc == 0) return anode_prob;

	bool first = true;
	for(unsigned int idx=0; idx < nc; ++idx)
	{
		// Copy to local var to avoid aliasing
		double* anode_other_tree_prob = aNode->mOtherTreeProb[idx];

		// If the node is in the same tree recurse and eventually save the value, else use the value
		if(aNode->isSameTree(idx))
		{
			const ForestNode *m = aNode->mChildrenList[idx];
			const unsigned int branch_id = m->mBranchId;
			const int leaf_codon = m->mLeafCodon;

			if(leaf_codon >= 0)
			{
				if(first)
				{
					aSet.doTransitionAtLeaf(aSetIdx, branch_id, leaf_codon, anode_prob);
					if(anode_other_tree_prob) memcpy(anode_other_tree_prob+VECTOR_SLOT*aSetIdx, anode_prob, N*sizeof(double));
					first = false;
				}
				else
				{
					double ALIGN64 temp[N];
					double* x = anode_other_tree_prob ? anode_other_tree_prob+VECTOR_SLOT*aSetIdx : temp;
					aSet.doTransitionAtLeaf(aSetIdx, branch_id, leaf_codon, x);
					elementWiseMult(anode_prob, x);
				}
			}
			else
			{
				if(first)
				{
					aSet.doTransition(aSetIdx, branch_id, computeLikelihoodsWalkerTC(m, aSet, aSetIdx), anode_prob);
					if(anode_other_tree_prob) memcpy(anode_other_tree_prob+VECTOR_SLOT*aSetIdx, anode_prob, N*sizeof(double));
					first = false;
				}
				else
				{
					double ALIGN64 temp[N];
					double* x = anode_other_tree_prob ? anode_other_tree_prob+VECTOR_SLOT*aSetIdx : temp;
					aSet.doTransition(aSetIdx, branch_id, computeLikelihoodsWalkerTC(m, aSet, aSetIdx), x);
					elementWiseMult(anode_prob, x);
				}
			}
		}
		else
		{
			if(first)
			{
				memcpy(anode_prob, anode_other_tree_prob+VECTOR_SLOT*aSetIdx, N*sizeof(double));
				first = false;
			}
			else
			{
				elementWiseMult(anode_prob, anode_other_tree_prob+VECTOR_SLOT*aSetIdx);
			}
		}
	}
#ifdef USE_LAPACK
	switch(nc)
	{
	case 1:
		elementWiseMult(anode_prob, mInvCodonFreq);
		break;
	case 2:
		elementWiseMult(anode_prob, mInv2CodonFreq);
		break;
	case 3:
		elementWiseMult(anode_prob, mInv2CodonFreq);
		elementWiseMult(anode_prob, mInvCodonFreq);
		break;
	case 4:
		elementWiseMult(anode_prob, mInv2CodonFreq);
		elementWiseMult(anode_prob, mInv2CodonFreq);
		break;
	default:
		for(unsigned int idx=0; idx < nc; ++idx) elementWiseMult(anode_prob, mInvCodonFreq);
		break;
	}
#endif

#ifdef CHECK_ALGO
	for(unsigned int dx=0; dx < nc; ++dx)
	{
			ForestNode *m = aNode->mChildrenList[dx];
			const unsigned int branch_id = m->mBranchId;
			std::cerr << branch_id << ' ';
	}
	std::cerr << '<' << aNode->mBranchId << ' ' << mNodeNames[aNode->mBranchId+1] << "> ";
	for(int j=0; j<N; ++j) std::cerr << anode_prob[j] << ' ';
	std::cerr << std::endl;
#endif

	return anode_prob;
}

#endif

#ifdef USE_DAG

void Forest::loadForestIntoDAG(unsigned int aMaxCopies, unsigned int aCopyId, const ForestNode* aNode)
{
	if(aNode)
	{
		// For all the node children
		const unsigned int nc = aNode->mChildrenCount;
		for(unsigned int i = 0; i < nc; ++i)
		{
			// Load a dependency (children should be computed before parent)
			mDAG.loadDependency(aCopyId, reinterpret_cast<const void*>(aNode->mChildrenList[i]), reinterpret_cast<const void*>(aNode));

			// If the dependant is on the same tree, continue recursively
			if(aNode->isSameTree(i)) loadForestIntoDAG(aMaxCopies, aCopyId, aNode->mChildrenList[i]);
		}
	}
	else
	{
		// Load all the copies requested
		for(unsigned int cp=0; cp < aMaxCopies; ++cp)
		{
			// Load all the individual trees
			for(int i=static_cast<int>(mNumSites)-1; i >= 0; --i)
			{
				loadForestIntoDAG(aMaxCopies, cp, &mRoots[i]);
			}
		}

		// Finalize the loading
		mDAG.endLoadDependencies();
	}
}
#endif

