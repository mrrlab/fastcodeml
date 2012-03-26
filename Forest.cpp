
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include "Forest.h"
#include "ForestNode.h"
#include "Exceptions.h"
#include "MathSupport.h"
#include "MatrixSize.h"
#include "CompilerHints.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// Initialize the mask table so the index corresponds to the bit position
const unsigned short ForestNode::mMaskTable[MAX_NUM_CHILDREN] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};


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
			(*il)->mPreprocessingSupport->mSubtreeCodonsSignature.push_back(codon);

			// Set leaves probability vector (Nt copies)
#ifdef NEW_LIKELIHOOD
			for(int set=0; set < Nt; ++set) mProbs[VECTOR_SLOT*(node*(Nt*mNumSites)+set*mNumSites+site)+codon] = 1.0; // The rest already zeroed by assign()
#else
			for(int set=0; set < Nt; ++set) (*il)->mProb[set][codon] = 1.0; // The rest already zeroed by assign()
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
	mCodonFreq = cf->getCodonFrequencies();

	// Set the mapping from internal branch number to branch number (the last tree has no pruned subtrees)
	std::map<unsigned int, unsigned int> map_internal_to_branchID;
	mapInternalToBranchIdWalker(&mRoots[mNumSites-1], map_internal_to_branchID);

	// Transform the map into a table (better for performance)
	mTableInternalToBranchID.resize(map_internal_to_branchID.size());
	std::map<unsigned int, unsigned int>::const_iterator im=map_internal_to_branchID.begin();
	const std::map<unsigned int, unsigned int>::const_iterator end=map_internal_to_branchID.end();
	for(; im != end; ++im) mTableInternalToBranchID[im->first] = im->second;

#ifdef NEW_LIKELIHOOD
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
		unsigned int max_len   = 0;
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

			unsigned int len = inbl->size();
			if(len > max_len) {max_len = len; max_level = level; max_leaf = leaf;}

		}

		// Find the first level that can inglobate the leave from level 'max_level' and index 'max_leaf'
		found = false;
		for(inbl=mNodesByLevel.begin()+max_level+1,level=max_level+1; inbl != mNodesByLevel.end(); ++inbl,++level)
		{
			const unsigned int len = inbl->size();
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

#endif
}


void Forest::reduceSubtrees(void)
{
	// Setup dependency vectors
	std::vector<unsigned int> empty_vector;
	mTreeDependencies.resize(mNumSites, empty_vector);
	mTreeRevDependencies.resize(mNumSites, empty_vector);

	// Try to merge equal subtrees
	// Trees at the beginning of the forest point to trees ahead
	// (this way a delete does not choke with pointers pointing to freed memory) 
	for(int i=mNumSites-1; i > 0; --i)
	{
		for(int j=i-1; j >= 0; --j)
		{
			reduceSubtreesWalker(&mRoots[i], &mRoots[j]);
		}
	}
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

void Forest::measureEffort(std::vector<unsigned int>& aEffort)
{
	// Initialize effort array
	aEffort.clear();
	aEffort.reserve(mNumSites);

	// Effort is number of branches plus one to have always non zero values
	for(size_t i=0; i < mNumSites; ++i)
	{
		unsigned int cnt_aggressive = mRoots[i].countBranches(true)+1;
		aEffort.push_back(cnt_aggressive);
	}
}


void Forest::printEffortByGroup(const std::vector<unsigned int>& aEffort, unsigned int aHyp)
{
	std::vector<unsigned int> core_effort;
#ifdef _OPENMP
	unsigned int nthreads = omp_get_max_threads();
#else
	unsigned int nthreads = 1;
#endif

	const unsigned int num_classes = mDependenciesClassesAndTrees[aHyp].size();
	for(unsigned int k=0; k < num_classes; ++k)
	{
		const unsigned int class_num_sites = mDependenciesClassesAndTrees[aHyp][k].size();
		core_effort.assign(nthreads, 0);

		unsigned int sites_per_core_base = class_num_sites/nthreads;
		unsigned int sites_per_core_plus = class_num_sites-sites_per_core_base*nthreads;
		if(mVerbose >= 1) std::cerr << nthreads << std::setw(4) << sites_per_core_base << std::setw(4) << sites_per_core_plus << ' ';
		for(unsigned int j=0; j < class_num_sites; ++j)
		{
			unsigned int site = mDependenciesClassesAndTrees[aHyp][k][j].first;
			
			unsigned int idx = 0;
			if(sites_per_core_plus == 0)
			{
				idx = j/sites_per_core_base;
			}
			else if(j >= (sites_per_core_base+1)*sites_per_core_plus)
			{
				idx = (j-(sites_per_core_base+1)*sites_per_core_plus)/sites_per_core_base+sites_per_core_plus;
			}
			else
			{
				idx = j/(sites_per_core_base+1);
			}
			core_effort[idx] += aEffort[site];
		}
		if(mVerbose >= 1)
		{
			std::cerr << "Trees in class " << std::setw(3) << k << ": " << std::setw(4) << class_num_sites << " |";
			for(unsigned int i=0; i < nthreads; ++i) std::cerr << std::setw(4) << core_effort[i];
			std::cerr << std::endl;
		}
	}
}

unsigned int Forest::totalEffort(const std::vector<unsigned int>& aEffort, unsigned int aHyp)
{
	unsigned int total_effort = 0;

	std::vector<unsigned int> core_effort;
#ifdef _OPENMP
	unsigned int nthreads = omp_get_max_threads();
#else
	unsigned int nthreads = 1;
#endif

	const unsigned int num_classes = mDependenciesClassesAndTrees[aHyp].size();
	for(unsigned int k=0; k < num_classes; ++k)
	{
		const unsigned int class_num_sites = mDependenciesClassesAndTrees[aHyp][k].size();
		core_effort.assign(nthreads, 0);

		unsigned int sites_per_core_base = class_num_sites/nthreads;
		unsigned int sites_per_core_plus = class_num_sites-sites_per_core_base*nthreads;

		for(unsigned int j=0; j < class_num_sites; ++j)
		{
			unsigned int site = mDependenciesClassesAndTrees[aHyp][k][j].first;
			
			unsigned int idx = 0;
			if(sites_per_core_plus == 0)
			{
				idx = j/sites_per_core_base;
			}
			else if(j >= (sites_per_core_base+1)*sites_per_core_plus)
			{
				idx = (j-(sites_per_core_base+1)*sites_per_core_plus)/sites_per_core_base+sites_per_core_plus;
			}
			else
			{
				idx = j/(sites_per_core_base+1);
			}
			core_effort[idx] += aEffort[site];
		}

		unsigned int max_effort = 0;
		for(unsigned int i=0; i < nthreads; ++i) if(core_effort[i] > max_effort) max_effort = core_effort[i];
		total_effort += max_effort;
	}

	return total_effort;
}


void Forest::prepareDependencies(bool aForceSerial)
{
	// Compute the dependencies: for each class list all sites that should be done at this level
	groupByDependency(aForceSerial);
	
	// Compute effort per site
	std::vector<unsigned int> effort;
	measureEffort(effort);

	// Print count and total effort before balancing
	if(!aForceSerial && mVerbose >= 1)
	{
		std::cerr << std::endl << "Before balancing (Classes&Trees)" << std::endl;
		printDependenciesClassesAndTrees();
		std::cerr << std::endl;
		std::cerr << "H0 total: " << totalEffort(effort, 0) << std::endl;
		std::cerr << "H1 total: " << totalEffort(effort, 1) << std::endl;
	}

	// Balance the precomputed dependency lists for both hypothesis
	bool done_ct = balanceDependenciesClassesAndTrees(aForceSerial, 0, false) &&
				   balanceDependenciesClassesAndTrees(aForceSerial, 1, false);

	// Equalize max effort inside each class
	balanceEffort(effort, 0);
	balanceEffort(effort, 1);

	// Print classes after balancing
	if(!aForceSerial && done_ct && mVerbose >= 1)
	{
		std::cerr << std::endl << "After balancing (Classes&Trees)" << std::endl;
		printDependenciesClassesAndTrees();
		std::cerr << std::endl;
		std::cerr << "H0 total: " << totalEffort(effort, 0) << std::endl;
		std::cerr << "H1 total: " << totalEffort(effort, 1) << std::endl;
	}

	// Compute the effort (ie the number of branches+1) for each thread in each class
	//
	/// @todo Finish the intra class balancing using effort values
	//
	if(!aForceSerial && mVerbose >= 1)
	{
		std::cerr << std::endl << "Hypothesis H0" << std::endl;
		printEffortByGroup(effort, 0);

		std::cerr << std::endl << "Hypothesis H1" << std::endl;
		printEffortByGroup(effort, 1);
	}
}

void Forest::groupByDependency(bool aForceSerial)
{
	std::vector< std::vector<unsigned int> > class_dependencies;

	// If no dependencies
	if(aForceSerial)
	{
		std::vector<unsigned int> v(mNumSites);

		for(unsigned int k=0; k < (unsigned int)mNumSites; ++k) v[k] = (unsigned int)mNumSites-k-1; // Remember: prior (could) point to subsequent

		class_dependencies.push_back(v);
	}
	else
	{
		unsigned int i, j;

		// Prepare the search of dependencies
		std::vector<bool> done(mNumSites, false);	// The sites that has dependencies satisfied in the previous level
		std::vector<bool> prev;						// Dependencies till the previous level
		std::vector<unsigned int> v;				// Temporary list of sites

		// Mark trees without dependencies
		// mTreeDependencies[tj] can be done after: t1 t2 t3
		for(i=0; i < (unsigned int)mNumSites; ++i)
		{
			if(mTreeDependencies[i].empty())
			{
				done[i] = true;
				v.push_back(i);
			}
		}

		// Prepare the dependency list
		class_dependencies.push_back(v);
		prev = done;

		// Start to find trees with one, two, ... dependencies
		for(unsigned int numdep=1;; ++numdep)
		{
			v.clear();
			bool all_done = true;
			for(i=0; i < mNumSites; ++i)
			{
				// If tree i has been already processed skip it
				if(prev[i]) continue;
				all_done = false;

				unsigned int nc = mTreeDependencies[i].size();
				bool all = true;
				for(j=0; j < nc; ++j) if(!prev[mTreeDependencies[i][j]]) {all = false; break;}
				if(all)
				{
					v.push_back(i);
					done[i] = true;
				}
			}
			if(all_done) break;
			class_dependencies.push_back(v);
			prev = done;
		}
	}

	// Number of dependencies classes
	unsigned int nc = class_dependencies.size();
	
	// One dependency classes
	std::vector<std::pair<unsigned int, unsigned int> > one_class;

	// Transform the list in two lists one for each hypothesis and multiply the entries by the number of codon classes
	for(unsigned int h=0; h <= 1; ++h)
	{
		// Prepare the dependency list for H0: 3 codon classes; H1: 4 codon classes
		unsigned int num_classes = h ? 4 : 3;

		mDependenciesClassesAndTrees[h].clear();
		for(unsigned int i=0; i < nc; ++i)
		{
			// Prepare the dependency classe
			one_class.clear();

			// Number of trees in the class
			const unsigned int nt = class_dependencies[i].size();
			for(unsigned int cl=0; cl < num_classes; ++cl)
			{
				for(unsigned int j=0; j < nt; ++j)
				{
					one_class.push_back(std::pair<unsigned int, unsigned int>(class_dependencies[i][j], cl));
				}
			}
			mDependenciesClassesAndTrees[h].push_back(one_class);
		}
	}
}


bool Forest::balanceDependenciesClassesAndTrees(bool aForceSerial, int aHyp, bool aGreedy)
{
	// Do nothing if no dependencies or no parallel execution
	if(aForceSerial) return false;

#ifndef _OPENMP
	return false;
#else
	// At each level collect the 'jolly' threads (trees that are not preconditions for trees in classes above)
	// This step makes sense only if run multithread and if there are more than one class
	const unsigned int num_threads = omp_get_max_threads();
	const unsigned int num_classes = mDependenciesClassesAndTrees[aHyp].size();
	if(num_threads < 2 || num_classes < 2) return false;

	// This set will contain the site & class pairs that can be postponed without problem
	std::set<std::pair<unsigned int, unsigned int> > jolly_sites;

	// Check if the current level can be balanced
	for(unsigned int k=0; k < num_classes; ++k)
	{
		// Try to have num sites at this level multiple of number of threads
		const unsigned int num_jolly = jolly_sites.size();
		const unsigned int class_num_sites = mDependenciesClassesAndTrees[aHyp][k].size();
		const unsigned int over = class_num_sites % num_threads;

		// Compute how many sites to add to have a multiple of num threads
		unsigned int needed_add = 0;
		unsigned int resid_jolly = num_jolly;
		if(over)
		{
			unsigned int min_add = num_threads - over;
			if(min_add <= num_jolly)
			{
				needed_add = min_add;
				resid_jolly -= needed_add;
			}
		}
		if(aGreedy)
		{
			needed_add += (resid_jolly/num_threads)*num_threads;
		}

		// Compute how many sites could became a jolly
		unsigned int new_possible_jolly_sites = 0;
		for(unsigned int s=0; s < class_num_sites; ++s)
		{
			// mTreeRevDependencies[tj] should be ready before: t1 t2 t3
			unsigned int site = mDependenciesClassesAndTrees[aHyp][k][s].first;
			if(mTreeRevDependencies[site].empty()) ++new_possible_jolly_sites;
		}

		// Compute how many sites to remove to have a multiple of num threads
		unsigned int needed_remove = 0;
		resid_jolly = new_possible_jolly_sites;
		if(class_num_sites > num_threads && new_possible_jolly_sites >= over)
		{
			needed_remove = over;
			resid_jolly -= needed_remove;
		}
		if(aGreedy && resid_jolly >= 2*num_threads)
		{
			needed_remove += ((resid_jolly - num_threads)/num_threads)*num_threads;
		}

		// There are at least the number needed to add in the jolly list, so add them
		if(needed_add > 0)
		{
			for(unsigned int j=0; j < needed_add; ++j)
			{
				std::set<std::pair<unsigned int, unsigned int> >::iterator it = jolly_sites.begin();
				mDependenciesClassesAndTrees[aHyp][k].push_back(*it);
				jolly_sites.erase(it);
			}
		}
		else if(needed_remove > 0)
		{
			// Save what remains at this level
			std::vector<std::pair<unsigned int, unsigned int> > new_content;

			for(unsigned int s=0; s < class_num_sites; ++s)
			{
				unsigned int site = mDependenciesClassesAndTrees[aHyp][k][s].first;
				if(needed_remove && mTreeRevDependencies[site].empty())
				{
					jolly_sites.insert(mDependenciesClassesAndTrees[aHyp][k][s]);
					--needed_remove;
				}
				else
				{
					new_content.push_back(mDependenciesClassesAndTrees[aHyp][k][s]);
				}
			}
			mDependenciesClassesAndTrees[aHyp][k].swap(new_content);
		}
	}

	// If there are still jolly sites, add them to the last class
	if(!jolly_sites.empty())
	{
		mDependenciesClassesAndTrees[aHyp][num_classes-1].insert(mDependenciesClassesAndTrees[aHyp][num_classes-1].end(), jolly_sites.begin(), jolly_sites.end());
	}

	return true;
#endif
}


void Forest::balanceEffort(const std::vector<unsigned int>& aEffort, unsigned int aHyp)
{
#ifndef _OPENMP
	return;
#else
	// At each level collect the 'jolly' threads (trees that are not preconditions for trees in classes above)
	// This step makes sense only if run multithread and if there are more than one class
	const unsigned int num_threads = omp_get_max_threads();
	const unsigned int num_classes = mDependenciesClassesAndTrees[aHyp].size();
	if(num_threads < 2 || num_classes < 2) return;

	// Balance each level
	for(unsigned int k=1; k < num_classes; ++k)
	{
		const int class_num_sites = static_cast<int>(mDependenciesClassesAndTrees[aHyp][k].size());
		const unsigned int over   = class_num_sites % num_threads;
		const unsigned int blocks = class_num_sites / num_threads;

		// Class 0 and classes without full blocks or with only one level (full or partial) cannot be equalized
		if(blocks == 0) continue;
		if(blocks == 1 && over == 0) continue;
		//unsigned int levels = blocks + (over ? 1 : 0);

		// Create helper list of indices per thread
		std::vector<unsigned int> cnt;
		std::vector< std::vector< unsigned int> > counts;
		for(unsigned int j=0; j < num_threads; ++j) counts.push_back(cnt);

#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(class_num_sites, counts) schedule(static)
#else
		#pragma omp parallel for default(shared) schedule(static)
#endif
		for(int i=0; i < class_num_sites; ++i) counts[omp_get_thread_num()].push_back(i);

		unsigned int min_effort = UINT_MAX;
		unsigned int max_effort = 0;
		for(unsigned int nt=0; nt < num_threads; ++nt)
		{
			unsigned int ns = counts[nt].size();
			unsigned int total_effort_thread = 0;
			for(unsigned int j=0; j < ns; ++j)
			{
				unsigned int idx = counts[nt][j];
				unsigned int site = mDependenciesClassesAndTrees[aHyp][k][idx].first;
				total_effort_thread += aEffort[site];
			}
			if(total_effort_thread < min_effort) min_effort = total_effort_thread;
			if(total_effort_thread > max_effort) max_effort = total_effort_thread;
		}

		// Print the results for checking
		std::cerr << "H" << aHyp << " Balancing " << std::setw(3) << k << " Effort range: " << max_effort-min_effort << std::endl;
		//for(unsigned int nt=0; nt < num_threads; ++nt)
		//{
		//	std::cerr << "Thread" << std::setw(3) << nt << ": ";
		//	for(unsigned int j=0; j < counts[nt].size(); ++j) std::cerr << std::setw(5) << counts[nt][j];
		//	std::cerr << std::endl;
		//}
	}
#endif
}


void Forest::printDependenciesClassesAndTrees(void)
{
	for(unsigned int h=0; h <= 1; ++h)
	{
		std::cerr << std::endl << "Hypothesis H" << h << std::endl;
		const unsigned int num_classes = mDependenciesClassesAndTrees[h].size();
		for(unsigned int numdep=0; numdep < num_classes; ++numdep)
		{
			std::cerr << "Trees in class " << std::setw(3) << numdep << ": " << std::setw(4) << mDependenciesClassesAndTrees[h][numdep].size() << std::endl;
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
		cnt += aForest.mRoots[i].countBranches();
		cntAggressive += aForest.mRoots[i].countBranches(true);
	}
	aOut << "Reduced branches:   " << std::fixed << std::setw(7) << cnt << std::setw(8) << std::setprecision(2) << static_cast<double>(cnt*100.)/static_cast<double>(aForest.mNumBranches*aForest.mNumSites) << '%' << std::endl;
	aOut << "Aggressive reduct.: " << std::fixed << std::setw(7) << cntAggressive << std::setw(8) << std::setprecision(2) << static_cast<double>(cntAggressive*100.)/static_cast<double>(aForest.mNumBranches*aForest.mNumSites) << '%' << std::endl;
	aOut << std::endl;

	// Print forest
	if(aForest.mVerbose >= 2)
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
		const int num_sites = inbl->size();
        const int len       = num_sites*num_sets;

#ifdef _MSC_VER
        #pragma omp parallel for default(none) shared(aSet, len, inbl, num_sets, num_sites, level) schedule(static)
#else
        #pragma omp parallel for default(shared) schedule(static)
#endif
        for(int i=0; i < len; ++i)
        {
            // Compute probability vector along this branch (for the given set) (reordered to give a 2% speedup)
			const unsigned int set_idx  = i / num_sites;
			const unsigned int site_idx = i - set_idx * num_sites; // Was: unsigned int set_idx = i % num_sets;
            const unsigned int branch   = ((*inbl)[site_idx])->mBranchId;
			const unsigned int start    = VECTOR_SLOT*(mNumSites*Nt*(branch+1)+mNumSites*set_idx+mFatVectorTransform.getLowerIndex(branch));

			// For each branch, except the root, compute the transition
            aSet.doTransition(set_idx,
							  branch,
							  mFatVectorTransform.getCount(branch),
							  &mProbs[start],
							  &mProbsOut[start]);
        }

		// Combine the results to have the input for the next round
		mFatVectorTransform.postCompact(mProbsOut, mProbs, level, num_sets);
    }

	// Compute the final likelyhood
	const int num_sites = mNumSites;
    const int len       = num_sites*num_sets;

#ifdef _MSC_VER
    #pragma omp parallel for default(none) shared(len, num_sites, aLikelihoods) schedule(static)
#else
    #pragma omp parallel for default(shared) schedule(static)
#endif
    for(int i=0; i < len; ++i)
    {
		const unsigned int set_idx = i / num_sites;
		const unsigned int site    = i - set_idx * num_sites; // Was: unsigned int site_idx = i % num_sites;
		const unsigned int start   = VECTOR_SLOT*(set_idx*mNumSites+site);

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
	ListDependencies::iterator ivs=mDependenciesClassesAndTrees[aHyp].begin();
	const ListDependencies::iterator end=mDependenciesClassesAndTrees[aHyp].end();
	for(; ivs != end; ++ivs)
	{
		const int len = ivs->size();

#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(aSet, len, ivs, aLikelihoods) schedule(static)
#else
		#pragma omp parallel for default(shared) schedule(static)
#endif
		for(int i=0; i < len; ++i)
		{
			// Compute likelihood array at the root of one tree (the access order is the fastest)
			const unsigned int site    = (*ivs)[i].first;
			const unsigned int set_idx = (*ivs)[i].second;

			computeLikelihoodsWalkerNR(aSet, set_idx, site);

			aLikelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, mRoots[site].mProb[set_idx]);
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
		unsigned int ns = mNumSites;
		for(unsigned int i=0; i < ns; ++i) prepareNewReduction(&mRoots[i]);

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

void Forest::computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int aHyp)
{
	ListDependencies::iterator ivs=mDependenciesClassesAndTrees[aHyp].begin();
	const ListDependencies::iterator end=mDependenciesClassesAndTrees[aHyp].end();
	for(; ivs != end; ++ivs)
	{
		const int len = ivs->size();

#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(aSet, len, ivs, aLikelihoods) schedule(static)
#else
		#pragma omp parallel for default(shared) schedule(static)
#endif
		for(int i=0; i < len; ++i)
		{
			// Compute likelihood array at the root of one tree (the access order is the fastest)
			const unsigned int site    = (*ivs)[i].first;
			const unsigned int set_idx = (*ivs)[i].second;

			double* g = computeLikelihoodsWalkerTC(&mRoots[site], aSet, set_idx);

			aLikelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, g);
		}
	}
}


double* Forest::computeLikelihoodsWalkerTC(ForestNode* aNode, const ProbabilityMatrixSet& aSet, unsigned int aSetIdx)
{
	bool first = true;
	double* anode_prob = aNode->mProb[aSetIdx];

	const unsigned int nc = aNode->mChildrenCount;
	for(unsigned int idx=0; idx < nc; ++idx)
	{
		// Copy to local var to avoid aliasing
		double* anode_other_tree_prob = aNode->mOtherTreeProb[idx];

		// If the node is in the same tree recurse and eventually save the value, else use the value
		if(aNode->isSameTree(idx))
		{
			ForestNode *m = aNode->mChildrenList[idx];
			const unsigned int branch_id = m->mBranchId;

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

	return anode_prob;
}

#endif

