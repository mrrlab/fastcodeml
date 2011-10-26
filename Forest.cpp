
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <algorithm>
#include "Forest.h"
#include "ForestNode.h"
#include "Exceptions.h"
#include "MathSupport.h"
#include "MatrixSize.h"

const unsigned char ForestNode::mMaskTable[MAX_NUM_CHILDREN] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80};

#if 0
void crc(const std::vector<double>& v, unsigned int nsites)
{
	union { double value; unsigned char bytes[sizeof(double)]; } data;
	int  c1 = 52845; 
	int  c2 = 22719;
	unsigned int nv = v.size();

	for(unsigned int n=0; n < (nv/(N*nsites*Nt)); ++n)
	{
		std::cerr << std::setw(2) << n << ": ";
		for(unsigned int s=0; s < nsites; ++s)
		{
			long sum = 0;
			int r = 55665;
			for(int j=0; j < N; ++j)
			{
				data.value = v[n*(Nt*nsites*N)+s*(N)+j];

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


void Forest::loadTreeAndGenes(const PhyloTree& aTree, const Genes& aGenes, bool aIgnoreFreq)
{
	// Check coherence between tree and genes
	checkCoherence(aTree, aGenes);

	// Collect global data that refers to the tree and that should not be duplicated on each tree of the forest
	aTree.collectGlobalTreeData(mNodeNames, mBranchLengths, &mMarkedInternalBranch);

	// Number of branches of one tree
	mNumBranches = aTree.getNumBranches();

	// Count the number of unique sites
	mNumSites = aGenes.getNumSites();
	const unsigned int* mult = aGenes.getSiteMultiplicity();

	// Initialize the count of codon types
	std::vector<unsigned int> codon_count(N, 0);

	// Initialize the array of all probability vectors
	mProbs.assign(mNumSites*(mNumBranches+1)*Nt*VECTOR_SLOT, 0.0);
#ifdef NEW_LIKELIHOOD
	mProbsOut.assign(mNumSites*(mNumBranches+1)*Nt*N, 0.0);
#endif

	// Count of tree's leaves
	size_t num_leaves = 0;

	// Clone tree inside the forest
	mRoots.resize(mNumSites);
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Create a copy of the tree
		aTree.cloneTree(&mRoots[site], site, mNumSites, mProbs);

		// Create a list of pointers to leaves
		std::vector<ForestNode*> leaves;
		leaves.clear();
		mRoots[site].pushLeaf(leaves);
		num_leaves = leaves.size();

		// Add codon code to leaves
		std::vector<ForestNode*>::const_iterator il;
		for(il=leaves.begin(); il != leaves.end(); ++il)
		{
			// Node id (adjusted so root is 0)
			unsigned int node = (*il)->mBranchId+1;

			// Get the codon index and add it to the node signature
			int codon = aGenes.getCodonIdx(mNodeNames[node], site);
			(*il)->mSubtreeCodonsSignature.push_back(codon);

			// Set leaves probability vector (Nt copies)
#ifdef NEW_LIKELIHOOD
			for(int set=0; set < Nt; ++set) mProbs[node*(Nt*mNumSites*N)+set*(mNumSites*N)+site*(N)+codon] = 1.0; // The rest already zeroed by assign()
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
	mSiteMultiplicity.assign(mult, mult+mNumSites);

	// Set the codon frequencies and related values needed for the eigensolver
	CodonFrequencies* cf = CodonFrequencies::getInstance();
	cf->setCodonFrequencies(codon_count, (aIgnoreFreq) ? CodonFrequencies::CODON_FREQ_MODEL_UNIF : CodonFrequencies::CODON_FREQ_MODEL_F3X4);
	mCodonFreq = cf->getCodonFrequencies();

	// Set the mapping from internal branch number to branch number (the last tree has no pruned subtrees)
	std::map<unsigned int, unsigned int> map_internal_to_branchID;
	mapInternalToBranchIdWalker(&mRoots[mNumSites-1], map_internal_to_branchID);

	// Transform the map into a table (better for performance)
	mTableInternalToBranchID.resize(map_internal_to_branchID.size());
	std::map<unsigned int, unsigned int>::const_iterator im;
	for(im=map_internal_to_branchID.begin(); im != map_internal_to_branchID.end(); ++im)
	{
		mTableInternalToBranchID[im->first] = im->second;
	}

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
        std::vector<ForestNode*>::const_iterator il;
        for(il=curr_level.begin(); il != curr_level.end(); ++il)
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
		std::vector< std::vector<ForestNode*> >::iterator inbl;
		unsigned int level = 0;
		for(inbl=mNodesByLevel.begin(),level=0; inbl != mNodesByLevel.end(); ++inbl,++level)
		{
			unsigned int num_leaves = 0;
			unsigned int leaf = 0, i=0;
			std::vector<ForestNode*>::const_iterator ifn;
			for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn,++i)
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
			unsigned int len = inbl->size();
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
	// Trees at the beginning of the forest point to trees ahead
	// (this way a delete does not choke with pointers pointing to freed memory) 
	int i, j;

	// Try to merge equal subtrees
	for(i=mNumSites-1; i > 0; --i)
	{
		for(j=i-1; j >= 0; --j)
		{
			reduceSubtreesWalker(&mRoots[i], &mRoots[j]);
		}
	}
}


void Forest::reduceSubtreesWalker(ForestNode* aNode, ForestNode* aNodeDependent)
{
	size_t i;
	size_t nc = aNode->mChildrenList.size();
	for(i=0; i < nc; ++i)
	{
		// If one of the two has been already reduced, do nothing
		if(!(aNode->isSameTree(i)) || !(aNodeDependent->isSameTree(i))) continue;

		// Check if same subtree
		if(aNode->mChildrenList[i]->mSubtreeCodonsSignature == aNodeDependent->mChildrenList[i]->mSubtreeCodonsSignature)
		{
			delete aNodeDependent->mChildrenList[i];
			aNodeDependent->mChildrenList[i] = aNode->mChildrenList[i];
			aNodeDependent->resetSameTreeFlag(i);
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
		aNode->mSubtreeCodonsSignature.clear();

		// Clean the children
		std::vector<ForestNode *>::iterator icl;
		unsigned int i = 0;
		for(icl=aNode->mChildrenList.begin(); icl != aNode->mChildrenList.end(); ++icl,++i)
		{
			if(aNode->isSameTree(i)) cleanReductionWorkingData(*icl);
		}
	}
}


void Forest::groupByDependency(bool aForceSerial)
{
	size_t i;
	unsigned int k;
	mDependenciesClasses.clear();

	// If no dependencies
	if(aForceSerial)
	{
		std::vector<unsigned int> v(mNumSites);

		for(k=0; k < (unsigned int)mNumSites; ++k) v[k] = (unsigned int)mNumSites-k-1; // Remember: prior (could) point to subsequent

		mDependenciesClasses.push_back(v);
		if(mVerbose >= 1) std::cerr << std::endl << "Trees in class  0: " << std::setw(3) << v.size() << std::endl;

		return;
	}

	// Collect dependencies
	std::vector< std::set<unsigned int> > dependencies(mNumSites);
	for(i=0; i < mNumSites; ++i)
	{
		std::set<unsigned int> dep;
		groupByDependencyWalker(&mRoots[i], dep);
		dependencies[i] = dep;
	}

	// Prepare the search of dependencies
	std::vector<bool> done;			// The sites that has dependencies satisfied in the previous level
	std::vector<bool> prev;			// Dependencies till the previous level
	done.assign(mNumSites, false);

	// Trees without dependencies
	std::vector<unsigned int> v;
	std::vector< std::set<unsigned int> >::iterator is;
	for(is=dependencies.begin(),k=0; is != dependencies.end(); ++is,++k)
	{
		if(is->empty())
		{
			v.push_back(k);
			done[k] = true;
		}
	}
	prev = done;
	mDependenciesClasses.push_back(v);
	if(mVerbose >= 1) std::cerr << std::endl << "Trees in class   0: " << std::setw(4) << v.size() << std::endl;

	// Start to find trees with one, two, ... dependencies
	for(size_t numdep=1;; ++numdep)
	{
		v.clear();
		bool all_done = true;
		std::vector< std::set<unsigned int> >::reverse_iterator ris;
		for(ris=dependencies.rbegin(),k=mNumSites-1; ris != dependencies.rend(); ++ris,--k)
		{
			// If tree k has been already processed skip it
			if(done[k]) continue;

			all_done = false;
			bool all = true;
			std::set<unsigned int>::iterator iv;
			for(iv=ris->begin(); iv != ris->end(); ++iv)
			{
				if(!prev[*iv])
				{
					all = false;
					break;
				}
			}
			if(all)
			{
				v.push_back(k);
				done[k] = true;
			}
		}
		if(all_done) break;
		mDependenciesClasses.push_back(v);
		prev = done;
		if(mVerbose >= 1) std::cerr << "Trees in class " << std::setw(3) << numdep << ": " << std::setw(4) << v.size() << std::endl;
	}
}


void Forest::groupByDependencyWalker(ForestNode* aNode, std::set<unsigned int>& aDependency)
{
	size_t i;
	size_t nc = aNode->mChildrenList.size();
	for(i=0; i < nc; ++i)
	{
		if(aNode->isSameTree(i))
		{
			groupByDependencyWalker(aNode->mChildrenList[i], aDependency);
		}
		else
		{
			unsigned int id = aNode->mChildrenList[i]->mOwnTree;
			aDependency.insert(id);
		}
	}
}


 std::ostream& operator<< (std::ostream& aOut, const Forest& aObj)
{
	size_t i;

	aOut << std::endl;
	aOut << "Num branches:       " << std::setw(7) << aObj.mNumBranches << std::endl;
	aOut << "Unique sites:       " << std::setw(7) << aObj.mNumSites << std::endl;
	aOut << "Total branches:     " << std::setw(7) << aObj.mNumBranches*aObj.mNumSites << std::endl;

	// Count total branches on the reduced forest
	unsigned int cnt = 0;
	unsigned int cntAggressive = 0;
	for(i=0; i < aObj.mNumSites; ++i)
	{
		cnt += aObj.mRoots[i].countBranches();
		cntAggressive += aObj.mRoots[i].countBranches(true);
	}
	aOut << "Reduced branches:   " << std::fixed << std::setw(7) << cnt << std::setw(8) << std::setprecision(2) << (double)(cnt*100.)/(double)(aObj.mNumBranches*aObj.mNumSites) << '%' << std::endl;
	aOut << "Aggressive reduct.: " << std::fixed << std::setw(7) << cntAggressive << std::setw(8) << std::setprecision(2) << (double)(cntAggressive*100.)/(double)(aObj.mNumBranches*aObj.mNumSites) << '%' << std::endl;
	aOut << std::endl;

	// Print forest
	if(aObj.mVerbose >= 2)
	{
		for(i=0; i < aObj.mNumSites; ++i)
		{
			aOut << "=== Site " << i << " ===" << std::endl;
			aObj.mRoots[i].print(aObj.getNodeNames(), aOut);
			aOut << std::endl;
		}
	}

	return aOut;
}

void Forest::exportForest(const char* aFilename, unsigned int aCounter) const
 {
	 std::vector< std::pair<int, int> > node_from;
	 std::vector< std::pair<int, int> > node_to;
	 std::vector< double >				branch_length;

	 // Get all forest connections
	 std::vector<ForestNode>::const_iterator ifn;
	 for(ifn=mRoots.begin(); ifn != mRoots.end(); ++ifn)
	 {
		 exportForestWalker(&(*ifn), mBranchLengths, node_from, node_to, branch_length);
	 }

	// Remove duplicated nodes
	std::set< std::pair<int, int> > vertices;
	std::vector< std::pair<int, int> >::const_iterator ip;
	for(ip = node_from.begin(); ip != node_from.end(); ++ip) vertices.insert(*ip);
	for(ip = node_to.begin(); ip != node_to.end(); ++ip) vertices.insert(*ip);

	// Convert to node indices
	std::map<std::pair<int, int>, int> map;
	int idx;
	std::set< std::pair<int, int> >::const_iterator iv;
	for(iv=vertices.begin(), idx=1; iv != vertices.end(); ++iv, ++idx)
	{
		std::pair<int, int> p(iv->first, iv->second);
		map[p] = idx;
	}

	// Map values to branches (identified by the end node)
	std::map<std::pair<int, int>, double> map_value;
	std::map<std::pair<int, int>, bool> map_same;
	for(size_t i=0; i < node_to.size(); ++i)
	{
		map_value[node_to[i]] = branch_length[i];
	}

	// Check if the filename contains a %03d format
	char temp_filename[1024];
	if(strrchr(aFilename, '%'))
	{
		sprintf(temp_filename, aFilename, aCounter);
		aFilename = temp_filename;
	}
	else if(strrchr(aFilename, '@'))
	{
		char z[1024];
		strncpy(z, aFilename, 1023);
		z[1023] = '\0';
		char *p = strrchr(z, '@');
		*p = '%';
		sprintf(temp_filename, z, aCounter);
		aFilename = temp_filename;
	}

	// Open the file and write the forest
	std::ofstream net(aFilename, std::ios_base::trunc | std::ios_base::out);
	if(!net.good())
	{
		std::cerr << "Cannot create net file <" << aFilename << ">" << std::endl;
	}
	else
	{
		net << "graph [\n";
		net << "comment \"Created by FastCodeML\"\n";
		net << "directed 1\n";
		net << "Version 1\n";

		for(iv=vertices.begin(), idx=1; iv != vertices.end(); ++iv, ++idx)
		{
			net << "node [\n";
			net << "   id " << idx << '\n';
			if(iv->second == 0)
			{
				net << "   label \"Root " << iv->first << "\"\n";
				net << "   type 0\n";
			}
			else
			{
				std::string s = mNodeNames[iv->second];
				if(s.empty()) net << "   label \"" << iv->second << "\"\n";
				else          net << "   label \"" << s << "\"\n";
				net << "   type 1\n";
			}
			net << "]\n";
		}

		for(size_t i=0; i < node_from.size(); ++i)
		{
			std::pair<int, int> pf(node_from[i].first, node_from[i].second);
			std::pair<int, int> pt(node_to[i].first,   node_to[i].second);
			bool same_tree = node_from[i].first == node_to[i].first;

			net << "edge [\n";
			net << "   source " <<  map[pf] << '\n';
			net << "   target " <<  map[pt] << '\n';
			net << "   label \"" << std::fixed << std::setprecision(2) << map_value[pt] << (same_tree ? "" : "+") << "\"\n"; // Trailing plus means not same tree link
			net << "]\n";
		}
		net << "]\n";
		net.close();
	}
}


void Forest::exportForestWalker(const ForestNode* aNode,
								const std::vector<double>& aBranchLengths,
								std::vector< std::pair<int, int> >& aNodeFrom,
								std::vector< std::pair<int, int> >& aNodeTo,
								std::vector< double >& aLength) const
{
	int my_node_id = aNode->mBranchId+1;
	int my_tree_id = aNode->mOwnTree;

	std::vector<ForestNode *>::const_iterator ifn;
	for(ifn=aNode->mChildrenList.begin(); ifn != aNode->mChildrenList.end(); ++ifn)
	{
		int your_node_id = (*ifn)->mBranchId+1;
		int your_tree_id = (*ifn)->mOwnTree;

		std::pair<int, int> p_from(my_tree_id, my_node_id);
		std::pair<int, int> p_to(your_tree_id, your_node_id);

		aNodeFrom.push_back(p_from);
		aNodeTo.push_back(p_to);
		aLength.push_back(mBranchLengths[your_node_id]);

		if(your_tree_id == my_tree_id) exportForestWalker(*ifn, aBranchLengths, aNodeFrom, aNodeTo, aLength);
	}
}

void Forest::checkCoherence(const PhyloTree& aTree, const Genes& aGenes) const
{
	// Get the two species list
	std::vector<std::string> tree_species_list;
	aTree.getSpecies(tree_species_list);
	std::vector<std::string> genes_species_list;
	aGenes.getSpecies(genes_species_list);

	// Should at least have the same number of species
	if(tree_species_list.size() != genes_species_list.size()) throw FastCodeMLFatal("Different number of species in tree and genes");

	// Create correspondence between species names
	std::vector<std::string>::const_iterator is1, is2;
	for(is1=tree_species_list.begin(); is1 != tree_species_list.end(); ++is1)
	{
		bool found = false;
		for(is2=genes_species_list.begin(); is2 != genes_species_list.end(); ++is2)
		{
			if(*is1 == *is2) {found = true; break;}
		}
		if(!found) throw FastCodeMLFatal("Mismatch between species in tree and genes");
	}
}


void Forest::setTimesFromLengths(std::vector<double>& aTimes, const ForestNode* aNode) const
{
	if(!aNode) aNode = &mRoots[mNumSites-1];
	else
	{
		unsigned int idx = aNode->mBranchId+1;
		aTimes[aNode->mBranchId] = mBranchLengths[idx];
	}

	std::vector<ForestNode *>::const_iterator ifn;
	for(ifn=aNode->mChildrenList.begin(); ifn != aNode->mChildrenList.end(); ++ifn)
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
		std::vector<ForestNode>::iterator ifn;
		for(ifn=mRoots.begin(); ifn != mRoots.end(); ++ifn)
		{
			for(ifnp=ifn->mChildrenList.begin(); ifnp != ifn->mChildrenList.end(); ++ifnp)
			{
				setLengthsFromTimes(aTimes, *ifnp);
			}
		}
	}
	else
	{
		unsigned int idx = aNode->mBranchId+1;
		mBranchLengths[idx] = aTimes[aNode->mBranchId];

		for(ifnp=aNode->mChildrenList.begin(); ifnp != aNode->mChildrenList.end(); ++ifnp)
		{
			setLengthsFromTimes(aTimes, *ifnp);
		}
	}
}


#ifndef NEW_LIKELIHOOD
// Compute likelihood with the original approach
//
void Forest::computeLikelihood(const TransitionMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods)
{
	unsigned int num_sets = aSet.size();
	//aLikelihoods.assign(num_sets*mNumSites, 1.0);
	aLikelihoods.resize(num_sets*mNumSites);

	std::vector< std::vector<unsigned int> >::iterator ivs;
	for(ivs=mDependenciesClasses.begin(); ivs != mDependenciesClasses.end(); ++ivs)
	{
		int len = ivs->size()*num_sets;

#ifdef _MSC_VER
		#pragma omp parallel for if(len > 3) default(none) shared(aSet, len, ivs, num_sets, aLikelihoods)
#else
		#pragma omp parallel for default(shared)
#endif
		for(int i=0; i < len; ++i)
		{
			// Compute likelihood array at the root of one tree
			unsigned int site_idx = i / num_sets;
			unsigned int site = (*ivs)[site_idx];
			unsigned int set_idx = i - site_idx * num_sets; // Was: unsigned int set_idx = i % num_sets;

			double* g = computeLikelihoodWalker(&mRoots[site], aSet, set_idx);

			aLikelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, g);
		}
	}
}


double* Forest::computeLikelihoodWalker(ForestNode* aNode, const TransitionMatrixSet& aSet, unsigned int aSetIdx)
{
	//std::vector<ForestNode *>::iterator in;
	bool first = true;
	unsigned int idx;
	//for(in=aNode->mChildrenList.begin(), idx=0; in != aNode->mChildrenList.end(); ++in, ++idx)
	unsigned int nc = aNode->mChildrenList.size();
	for(idx=0; idx < nc; ++idx)
	{
		ForestNode *m = aNode->mChildrenList[idx];

		// If the node is in the same tree recurse, else use the value
		if(aNode->isSameTree(idx))
		{
			if(first)
			{
				//aSet.doTransition(aSetIdx, m->mBranchId, computeLikelihoodWalker(m, aSet, aSetIdx), aNode->mProb0+N*aSetIdx);
				aSet.doTransition(aSetIdx, m->mBranchId, computeLikelihoodWalker(m, aSet, aSetIdx), aNode->mProb[aSetIdx]);
				//if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mOtherTreeProb[idx]+N*aSetIdx, aNode->mProb0+N*aSetIdx, N*sizeof(double));
				if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mOtherTreeProb[idx]+N*aSetIdx, aNode->mProb[aSetIdx], N*sizeof(double));
				first = false;
			}
			else
			{
				double temp[N];
				double* x = aNode->mOtherTreeProb[idx] ? aNode->mOtherTreeProb[idx]+N*aSetIdx : temp;
				aSet.doTransition(aSetIdx, m->mBranchId, computeLikelihoodWalker(m, aSet, aSetIdx), x);
#ifdef USE_MKL_VML
				vdMul(N, aNode->mProb[aSetIdx], x, aNode->mProb[aSetIdx]);
#else
				//for(int i=0; i < N; ++i) aNode->mProb0[i+N*aSetIdx] *= x[i];
				for(int i=0; i < N; ++i) aNode->mProb[aSetIdx][i] *= x[i];
#endif
			}
		}
		else
		{
			if(first)
			{
				//if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mProb0+N*aSetIdx, aNode->mOtherTreeProb[idx]+N*aSetIdx, N*sizeof(double));
				if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mProb[aSetIdx], aNode->mOtherTreeProb[idx]+N*aSetIdx, N*sizeof(double));
				//else aSet.doTransition(aSetIdx, m->mBranchId, m->mProb0+N*aSetIdx, aNode->mProb0+N*aSetIdx);
				else aSet.doTransition(aSetIdx, m->mBranchId, m->mProb[aSetIdx], aNode->mProb[aSetIdx]);
				first = false;
			}
			else
			{
				double temp[N];
				double* x;
				if(aNode->mOtherTreeProb[idx]) 
				{
					x = aNode->mOtherTreeProb[idx]+N*aSetIdx;
				}
				else
				{
					//aSet.doTransition(aSetIdx, m->mBranchId, m->mProb0+N*aSetIdx, temp);
					aSet.doTransition(aSetIdx, m->mBranchId, m->mProb[aSetIdx], temp);
					x = temp;
				}
#ifdef USE_MKL_VML
				vdMul(N, aNode->mProb[aSetIdx], x, aNode->mProb[aSetIdx]);
#else
				//for(int i=0; i < N; ++i) aNode->mProb0[i+N*aSetIdx] *= x[i];
				for(int i=0; i < N; ++i) aNode->mProb[aSetIdx][i] *= x[i];
#endif
			}
		}
	}

	//return aNode->mProb0+N*aSetIdx;
	return aNode->mProb[aSetIdx];
}
#else
// Compute likelihood with the new "Long Vector" approach
//
void Forest::computeLikelihood(const TransitionMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods)
{
	// Initialize variables
    unsigned int num_sets = aSet.size();
    aLikelihoods.assign(num_sets*mNumSites, 1.0);

	// For each level of the tree (except the root)
	unsigned int level=0;
    std::vector< std::vector<ForestNode*> >::reverse_iterator inbl;
    for(inbl=mNodesByLevel.rbegin(); inbl != mNodesByLevel.rend(); ++inbl,++level)
    {
        int len = inbl->size()*num_sets;
#ifdef _MSC_VER
        #pragma omp parallel for default(none) shared(aSet, len, inbl, num_sets, level)
#else
        #pragma omp parallel for default(shared)
#endif
        for(int i=0; i < len; ++i)
        {
            // Compute probability vector along this branch (for the given set)
			unsigned int site_idx = i / num_sets;
			unsigned int set_idx  = i - site_idx * num_sets; // Was: unsigned int set_idx = i % num_sets;
            unsigned int branch   = ((*inbl)[site_idx])->mBranchId;
			unsigned int start    = N*mNumSites*Nt*(branch+1)+N*mNumSites*set_idx+N*mFatVectorTransform.getLowerIndex(branch);

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
    int len = mNumSites*num_sets;
#ifdef _MSC_VER
    #pragma omp parallel for default(none) shared(len, num_sets, aLikelihoods)
#else
    #pragma omp parallel for default(shared)
#endif
    for(int i=0; i < len; ++i)
    {
        unsigned int site    = i / num_sets;
        unsigned int set_idx = i - site * num_sets; // Was: unsigned int set_idx = i % num_sets;
		unsigned int start   = set_idx*mNumSites*N+site*N;

		// Take the result from branch 0 (the root)
        aLikelihoods[set_idx*mNumSites+site] = dot(mCodonFreq, &mProbs[start]);
    }
}
#endif

void Forest::mapInternalToBranchIdWalker(const ForestNode* aNode, std::map<unsigned int, unsigned int>& aMapInternalToBranchID)
{
	std::vector<ForestNode *>::const_iterator in;
	for(in=aNode->mChildrenList.begin(); in != aNode->mChildrenList.end(); ++in)
	{
		ForestNode *m = *in;

		if(m->mInternalNodeId != UINT_MAX) aMapInternalToBranchID[m->mInternalNodeId] = m->mBranchId;

		mapInternalToBranchIdWalker(m, aMapInternalToBranchID);
	}
}


#ifndef NEW_LIKELIHOOD
void Forest::addAggressiveReduction(ForestNode* aNode)
{
	if(aNode)
	{
		unsigned int i;
		std::vector<ForestNode *>::iterator in;
		for(in=aNode->mChildrenList.begin(), i=0; in != aNode->mChildrenList.end(); ++in, ++i)
		{
			ForestNode *m = *in;

			if(aNode->isSameTree(i))
			{
				addAggressiveReduction(m);
			}
			else
			{
				ForestNode *other = m->mParent;

				// Add the array on the other side
				if(!other->mOtherTreeProb[i]) other->mOtherTreeProb[i] = new double[N*Nt];

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

#if 0
void Forest::addAggressiveReductionWalker(ForestNode* aNode)
{
	unsigned int i;
	//std::vector<ForestNode *>::iterator in;
	unsigned int nc = aNode->mChildrenList.size();
	//for(in=aNode->mChildrenList.begin(), i=0; in != aNode->mChildrenList.end(); ++in, ++i)
	for(i=0; i < nc; ++i)
	{
		ForestNode *m = aNode->mChildrenList[i];

		if(aNode->isSameTree(i))
		{
			addAggressiveReductionWalker(m);
		}
		else
		{
			ForestNode *other = m->mParent;

			// Add the array on the other side
			if(!other->mOtherTreeProb[i]) other->mOtherTreeProb[i] = new double[N*Nt];

			// Add the pointer here
			aNode->mOtherTreeProb[i] = other->mOtherTreeProb[i];
		}
	}
}
#endif
#endif

#ifdef NEW_LIKELIHOOD
void Forest::prepareNewReduction(ForestNode* aNode)
{
	if(aNode)
	{
		unsigned int i;
		//std::vector<ForestNode *>::iterator icl;
		unsigned int nc = aNode->mChildrenList.size();
		//for(icl=aNode->mChildrenList.begin(),i=0; icl != aNode->mChildrenList.end(); ++icl,++i)
		for(i=0; i < nc; ++i)
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
		for(size_t i=0; i < mNumSites; ++i) prepareNewReduction(&mRoots[i]);

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
		std::vector<ForestNode *>::const_iterator icl;
		for(icl=aNode->mChildrenList.begin(),idx=0; icl != aNode->mChildrenList.end(); ++icl,++idx)
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
#endif

