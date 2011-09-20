
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

void Forest::loadTreeAndGenes(const PhyloTree& aTree, const Genes& aGenes, bool aIgnoreFreq)
{
	// Check coherence between tree and genes
	checkCoherence(aTree, aGenes);

	// Collect global data that refers to the tree and that should not be duplicated on each tree of the forest
	aTree.collectGlobalTreeData(mNodeNames, mBranchLengths, &mMarkedInternalBranch);

	// Number of branches of one tree
	mNumBranches = aTree.getNumBranches();

	// Count the number of unique sites
	size_t nsites = aGenes.getNumSites();
	const unsigned int* mult = aGenes.getSiteMultiplicity();

	// Initialize the count of codon types
	memset(mCodonCount, 0, N*sizeof(unsigned int));

	// Initialize the array of all probability vectors
	mProbs.assign(nsites*(mNumBranches+1)*Nt*N, 0.0);
	mProbsOut.assign(nsites*(mNumBranches+1)*Nt*N, 0.0);

	// Count of tree's leaves
	size_t num_leaves = 0;

	// Clone tree inside the forest
	mRoots.resize(nsites);
	for(unsigned int j=0; j < nsites; ++j)
	{
		// Create a copy of the tree
		aTree.cloneTree(&mRoots[j], j, nsites, mProbs);

		// Create a list of pointers to leaves
		std::vector<ForestNode*> leaves;
		leaves.clear();
		mRoots[j].pushLeaf(leaves);
		num_leaves = leaves.size();

		// Add codon code to leaves
		std::vector<ForestNode*>::const_iterator il;
		for(il=leaves.begin(); il != leaves.end(); ++il)
		{
			// Node id (adjusted so root is 0)
			unsigned int id = (*il)->mNodeId+1;

			// Get the codon index and add it to the node signature
			int codon = aGenes.getCodonIdx(mNodeNames[id], j);
			(*il)->mSubtreeCodonsSignature.push_back(codon);

			// Set leaves probability vector (Nt copies)
#ifdef NEW_LIKELIHOOD
			for(int k=0; k < Nt; ++k) mProbs[N*nsites*Nt*id+N*nsites*k+N*j+codon] = 1.0; // The rest already zeroed by assign()
#else
			for(int k=0; k < Nt; ++k) (*il)->mProb[k][codon] = 1.0; // The rest already zeroed by assign()
#endif

			// Count codons
			mCodonCount[codon] += mult[j];
		}

		// Combine the subtrees signatures going up to the root
		mRoots[j].gatherCodons();
	}

	// Set the number of internal branches
	mNumInternalBranches = mNumBranches - num_leaves;

	// Set the site multeplicity
	mSiteMultiplicity.resize(nsites);
#ifdef _MSC_VER
        #pragma omp parallel for default(none) shared(mult, nsites)
#else
        #pragma omp parallel for default(shared)
#endif
	for(int i=0; i < (int)nsites; ++i)
	{
		mSiteMultiplicity[i] = (double)mult[i];
	}

	// Set the codon frequencies and related values needed for the eigensolver
	if(aIgnoreFreq)
		setCodonFrequenciesUnif();
	else
		setCodonFrequenciesF3x4();

	// Set the mapping from internal branch number to branch number
	mapInternalToBranchIdWalker(&mRoots[0]);

    // Prepare the list of node id's by level
    std::vector<ForestNode*> next_level;
    std::vector<ForestNode*> curr_level;
    std::vector<ForestNode*> level_nodes;

    // First level is the root (but it is not added because no processing is done on it)
    //level_nodes.push_back(&mRoots[0]);
    mNodesByLevel.clear();
    //mNodesByLevel.push_back(level_nodes);
    curr_level.push_back(&mRoots[0]);

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

	// Push only the branch id's
	mBranchByLevel.clear();
	std::vector< std::vector<ForestNode*> >::reverse_iterator inbl;
	for(inbl=mNodesByLevel.rbegin(); inbl != mNodesByLevel.rend(); ++inbl)
	{
		std::vector<unsigned int> v;
		std::vector<ForestNode*>::iterator ifn;
		for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
		{
			v.push_back((*ifn)->mNodeId);
		}
		mBranchByLevel.push_back(v);
	}
}


void Forest::reduceSubtrees(void)
{
	// Trees at the beginning of the forest point to trees ahead
	// (this way a delete does not choke with pointers pointing to freed memory) 
	int i, j;
	int nsites = (int)mRoots.size();

	// Try to merge equal subtrees
	for(i=nsites-1; i > 0; --i)
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
		if(aNode->mChildrenList[i]->mOwnTree != aNode->mOwnTree) continue;
		if(aNodeDependent->mChildrenList[i]->mOwnTree != aNodeDependent->mOwnTree) continue;

		// Check if same subtree
		if(aNode->mChildrenList[i]->mSubtreeCodonsSignature == aNodeDependent->mChildrenList[i]->mSubtreeCodonsSignature)
		{
			delete aNodeDependent->mChildrenList[i];
			aNodeDependent->mChildrenList[i] = aNode->mChildrenList[i];
		}
	}

	// Recurse
	for(i=0; i < nc; ++i)
	{
		// If one of the two has been already reduced, do nothing
		if(aNode->mChildrenList[i]->mOwnTree != aNode->mOwnTree) continue;
		if(aNodeDependent->mChildrenList[i]->mOwnTree != aNodeDependent->mOwnTree) continue;

		reduceSubtreesWalker(aNode->mChildrenList[i], aNodeDependent->mChildrenList[i]);
	}
}


void Forest::cleanReductionWorkingData(ForestNode* aNode)
{
	if(!aNode)
	{
		// Invoke on all the trees in the forest
		size_t nsites = mRoots.size();
		for(size_t i=0; i < nsites; ++i) cleanReductionWorkingData(&mRoots[i]);
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
			if((*icl)->mOwnTree == aNode->mOwnTree) cleanReductionWorkingData(*icl);
		}
	}
}


void Forest::groupByDependency(bool aForceSerial)
{
	size_t i;
	unsigned int k;
	size_t nsites = mRoots.size();
	mDependenciesClasses.clear();

	// If no dependencies
	if(aForceSerial)
	{
		std::vector<unsigned int> v(nsites);

		for(k=0; k < (unsigned int)nsites; ++k) v[k] = (unsigned int)nsites-k-1; // Remember: prior (could) point to subsequent

		mDependenciesClasses.push_back(v);
		if(mVerbose >= 1) std::cerr << std::endl << "Trees in class  0: " << std::setw(3) << v.size() << std::endl;

		return;
	}

	// Collect dependencies
	std::vector< std::set<unsigned int> > dependencies(nsites);
	for(i=0; i < nsites; ++i)
	{
		std::set<unsigned int> dep;
		groupByDependencyWalker(&mRoots[i], dep);
		dependencies[i] = dep;
	}

	// Prepare the search of dependencies
	std::vector<bool> done;			// The sites that has dependencies satisfied in the previous level
	std::vector<bool> prev;			// Dependencies till the previous level
	done.assign(nsites, false);

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
	if(mVerbose >= 1) std::cerr << std::endl << "Trees in class  0: " << std::setw(3) << v.size() << std::endl;

	// Start to find trees with one, two, ... dependencies
	for(size_t numdep=1;; ++numdep)
	{
		v.clear();
		bool all_done = true;
		std::vector< std::set<unsigned int> >::reverse_iterator ris;
		for(ris=dependencies.rbegin(),k=nsites-1; ris != dependencies.rend(); ++ris,--k)
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
		if(mVerbose >= 1) std::cerr << "Trees in class " << std::setw(2) << numdep << ": " << std::setw(3) << v.size() << std::endl;
	}
}


void Forest::groupByDependencyWalker(ForestNode* aNode, std::set<unsigned int>& aDependency)
{
	size_t i;
	size_t nc = aNode->mChildrenList.size();
	for(i=0; i < nc; ++i)
	{
		if(aNode->mChildrenList[i]->mOwnTree == aNode->mOwnTree)
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
	aOut << "Unique sites:       " << std::setw(7) << aObj.mRoots.size() << std::endl;
	aOut << "Total branches:     " << std::setw(7) << aObj.mNumBranches*aObj.mRoots.size() << std::endl;

	// Count total branches on the reduced forest
	unsigned int cnt = 0;
	unsigned int cntAggressive = 0;
	for(i=0; i < aObj.mRoots.size(); ++i)
	{
		cnt += aObj.mRoots[i].countBranches();
		cntAggressive += aObj.mRoots[i].countBranches(true);
	}
	aOut << "Reduced branches:   " << std::fixed << std::setw(7) << cnt << std::setw(8) << std::setprecision(2) << (double)(cnt*100.)/(double)(aObj.mNumBranches*aObj.mRoots.size()) << '%' << std::endl;
	aOut << "Aggressive reduct.: " << std::fixed << std::setw(7) << cntAggressive << std::setw(8) << std::setprecision(2) << (double)(cntAggressive*100.)/(double)(aObj.mNumBranches*aObj.mRoots.size()) << '%' << std::endl;
	aOut << std::endl;

	// Print forest
	if(aObj.mVerbose >= 2)
	{
		for(i=0; i < aObj.mRoots.size(); ++i)
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
		strcpy(z, aFilename);
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
	int my_node_id = (aNode->mNodeId == UINT_MAX) ? 0 : aNode->mNodeId+1;
	int my_tree_id = aNode->mOwnTree;

	std::vector<ForestNode *>::const_iterator ifn;
	for(ifn=aNode->mChildrenList.begin(); ifn != aNode->mChildrenList.end(); ++ifn)
	{
		int your_node_id = ((*ifn)->mNodeId == UINT_MAX) ? 0 : (*ifn)->mNodeId+1;
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
	if(!aNode) aNode = &mRoots[0];
	else
	{
		unsigned int idx = (aNode->mNodeId == UINT_MAX) ? 0 : aNode->mNodeId+1;
		aTimes[aNode->mNodeId] = mBranchLengths[idx];
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
		unsigned int idx = (aNode->mNodeId == UINT_MAX) ? 0 : aNode->mNodeId+1;
		mBranchLengths[idx] = aTimes[aNode->mNodeId];

		for(ifnp=aNode->mChildrenList.begin(); ifnp != aNode->mChildrenList.end(); ++ifnp)
		{
			setLengthsFromTimes(aTimes, *ifnp);
		}
	}
}


#ifndef NEW_LIKELIHOOD
// Compute likelihood with the original approach
//
void Forest::computeLikelihood(const TransitionMatrixSet& aSet, std::vector<double>& aLikelihoods)
{
	unsigned int num_sets = aSet.size();
	size_t num_sites = mRoots.size();
	aLikelihoods.resize(num_sets*num_sites, 1.0);

	std::vector< std::vector<unsigned int> >::iterator ivs;
	for(ivs=mDependenciesClasses.begin(); ivs != mDependenciesClasses.end(); ++ivs)
	{
		int len = ivs->size()*num_sets;

#ifdef _MSC_VER
		#pragma omp parallel for if(len > 3) default(none) shared(aSet, len, ivs, num_sets, num_sites, aLikelihoods)
#else
		#pragma omp parallel for default(shared)
#endif
		for(int i=0; i < len; ++i)
		{
			// Compute likelihood array at the root of one tree
			unsigned int site = ivs->at(i / num_sets);
			unsigned int set_idx = i % num_sets;
			double* g = computeLikelihoodWalker(&mRoots[site], aSet, set_idx);
			aLikelihoods[set_idx*num_sites+site] = dot(mCodonFrequencies, g);
		}
	}
}


double* Forest::computeLikelihoodWalker(ForestNode* aNode, const TransitionMatrixSet& aSet, unsigned int aSetIdx)
{
	std::vector<ForestNode *>::iterator in;
	bool first = true;
	int idx;
	for(in=aNode->mChildrenList.begin(), idx=0; in != aNode->mChildrenList.end(); ++in, ++idx)
	{
		ForestNode *m = *in;

		// If the node is in the same tree recurse, else use the value
		if((*in)->mOwnTree == aNode->mOwnTree)
		{
			if(first)
			{
				//aSet.doTransition(aSetIdx, m->mNodeId, computeLikelihoodWalker(m, aSet, aSetIdx), aNode->mProb0+N*aSetIdx);
				aSet.doTransition(aSetIdx, m->mNodeId, computeLikelihoodWalker(m, aSet, aSetIdx), aNode->mProb[aSetIdx]);
				//if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mOtherTreeProb[idx]+N*aSetIdx, aNode->mProb0+N*aSetIdx, N*sizeof(double));
				if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mOtherTreeProb[idx]+N*aSetIdx, aNode->mProb[aSetIdx], N*sizeof(double));
				first = false;
			}
			else
			{
				double temp[N];
				double* x = aNode->mOtherTreeProb[idx] ? aNode->mOtherTreeProb[idx]+N*aSetIdx : temp;
				aSet.doTransition(aSetIdx, m->mNodeId, computeLikelihoodWalker(m, aSet, aSetIdx), x);
				//for(int i=0; i < N; ++i) aNode->mProb0[i+N*aSetIdx] *= x[i];
				for(int i=0; i < N; ++i) aNode->mProb[aSetIdx][i] *= x[i];
			}
		}
		else
		{
			if(first)
			{
				//if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mProb0+N*aSetIdx, aNode->mOtherTreeProb[idx]+N*aSetIdx, N*sizeof(double));
				if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mProb[aSetIdx], aNode->mOtherTreeProb[idx]+N*aSetIdx, N*sizeof(double));
				//else aSet.doTransition(aSetIdx, m->mNodeId, m->mProb0+N*aSetIdx, aNode->mProb0+N*aSetIdx);
				else aSet.doTransition(aSetIdx, m->mNodeId, m->mProb[aSetIdx], aNode->mProb[aSetIdx]);
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
					//aSet.doTransition(aSetIdx, m->mNodeId, m->mProb0+N*aSetIdx, temp);
					aSet.doTransition(aSetIdx, m->mNodeId, m->mProb[aSetIdx], temp);
					x = temp;
				}
				//for(int i=0; i < N; ++i) aNode->mProb0[i+N*aSetIdx] *= x[i];
				for(int i=0; i < N; ++i) aNode->mProb[aSetIdx][i] *= x[i];
			}
		}
	}

	//return aNode->mProb0+N*aSetIdx;
	return aNode->mProb[aSetIdx];
}
#else
// Compute likelihood with the new "Long Vector" approach
//
void Forest::computeLikelihood(const TransitionMatrixSet& aSet, std::vector<double>& aLikelihoods)
{
	// Initialize variables
    unsigned int num_sets = aSet.size();
    size_t      num_sites = mRoots.size();
    aLikelihoods.resize(num_sets*num_sites, 1.0);

	// For each level of the tree (except the root)
    std::vector< std::vector<ForestNode*> >::reverse_iterator inbl;
    for(inbl=mNodesByLevel.rbegin(); inbl != mNodesByLevel.rend(); ++inbl)
    {
        int len = inbl->size()*num_sets;

#ifdef _MSC_VER
        #pragma omp parallel for default(none) shared(aSet, len, inbl, num_sets, aLikelihoods, num_sites, std::cerr)
#else
        #pragma omp parallel for default(shared)
#endif
        for(int i=0; i < len; ++i)
        {
            // Compute probability vector along this branch (for the given set)
            unsigned int set_idx = i % num_sets;
            unsigned int branch  = (inbl->at(i / num_sets))->mNodeId+1;

            // For each branch, except the root, compute the transition
            aSet.doTransition2(set_idx,
							   branch-1,
							   num_sites,
							   &mProbs[N*num_sites*Nt*branch+N*num_sites*set_idx],
							   &mProbsOut[N*num_sites*Nt*branch+N*num_sites*set_idx]);
        }

        // Compose with the other results for the branches starting from this node
        ForestNode* curr_node = 0;
        std::vector<ForestNode*>::iterator ifn;
        for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
        {
			ForestNode*    parent_node = (*ifn)->mParent;
            unsigned int parent_branch = parent_node->mNodeId+1;
            unsigned int     my_branch = (*ifn)->mNodeId+1;

			// If this is the first visit to the parent copy the result, otherwise do a element by element mulktiplication
            if(parent_node != curr_node)
            {
                curr_node = parent_node;
				memcpy(&mProbs[N*num_sites*Nt*parent_branch], &mProbsOut[N*num_sites*Nt*my_branch], N*num_sites*Nt*sizeof(double));
            }
            else
            {
#ifdef _MSC_VER
        #pragma omp parallel for default(none) shared(parent_branch, my_branch, num_sites)
#else
        #pragma omp parallel for default(shared)
#endif
                for(int i=0; i < (int)(N*num_sites*Nt); ++i)
                {
                    mProbs[N*num_sites*Nt*parent_branch+i] *= mProbsOut[N*num_sites*Nt*my_branch+i];
                }
            }
        }
    }

	// Compute the final likelyhood
    int len = num_sites*num_sets;
#ifdef _MSC_VER
    #pragma omp parallel for default(none) shared(len, num_sets, aLikelihoods, num_sites)
#else
    #pragma omp parallel for default(shared)
#endif
    for(int i=0; i < len; ++i)
    {
        unsigned int set_idx = i % num_sets;
        unsigned int site    = i / num_sets;
        aLikelihoods[set_idx*num_sites+site] = dot(mCodonFrequencies, &mProbs[set_idx*num_sites*N+site*N]);
    }
}
#endif


int Forest::codon64to61(unsigned int aId64) const
{
	if(aId64 > 63 || aId64 == 10 || aId64 == 11 || aId64 == 14) return -1;

	if(aId64 > 14) return aId64-3;
	if(aId64 > 11) return aId64-2;
	return aId64;
}


void Forest::setCodonFrequenciesF3x4(void)
{
	int k, j;

#ifdef CHECK_ALGO
	// Print the table of codon counts
	for(k=0; k < 64; ++k)
	{
		int id = codon64to61(k);
		if(id < 0) std::cerr << std::setw(5) << 0;
		else       std::cerr << std::setw(5) << mCodonCount[id];
		if(k % 4 == 3) std::cerr << std::endl;
	}
#endif
	// Compute the 3x4 table
	double fb3x4sg[12];

	memset(fb3x4sg, 0, 12*sizeof(double));

    for(k = 0; k < 64; k++)
    {
		int kk = codon64to61(k);
		if(kk < 0) continue;

        fb3x4sg[0 * 4 + k / 16]      += mCodonCount[kk];
        fb3x4sg[1 * 4 + (k / 4) % 4] += mCodonCount[kk];
        fb3x4sg[2 * 4 + k % 4]       += mCodonCount[kk];
    }

    for(j = 0; j < 3; j++)
    {
        double t = 0;
		for(k=0; k < 4; ++k) t += fb3x4sg[j*4+k];
		for(k=0; k < 4; ++k) fb3x4sg[j*4+k] /= t;
    }

#ifdef CHECK_ALGO
	for(k=0; k < 12; ++k)
	{
		std::cerr << std::fixed << std::setprecision(6) << fb3x4sg[k];
		if(k % 4 == 3) std::cerr << std::endl;
	}
	std::cerr << std::endl;
#endif

	// Compute codon frequency from the 3x4 table
	for(k=0; k < 64; ++k)
	{
		int kk = codon64to61(k);
		if(kk < 0) continue;

		mCodonFrequencies[kk] = fb3x4sg[k / 16] * fb3x4sg[4 + (k / 4) % 4] * fb3x4sg[8 + k % 4];
	}
	double t = 0;
	for(k=0; k < N; ++k) t += mCodonFrequencies[k];
	for(k=0; k < N; ++k) mCodonFrequencies[k] /= t;

	// Support values needed for the eigensolver
	for(k=mNumGoodCodons=0; k < N; ++k)
	{
		mCodonFreqSqrt[k] = sqrt(mCodonFrequencies[k]);
		if(mCodonFrequencies[k] > GOOD_CODON_THRESHOLD)
		{
			mGoodCodon[k] = true;
			++mNumGoodCodons;
		}
		else
		{
			mGoodCodon[k] = false;
		}
	}

#ifdef CHECK_ALGO
	for(k=0; k < 61; ++k)
	{
		std::cerr << std::fixed << std::setprecision(10) << mCodonFrequencies[k] << ' ';
		if(k % 4 == 3) std::cerr << std::endl;
	}
	std::cerr << std::endl;
#endif
}


void Forest::setCodonFrequenciesUnif(void)
{
	for(int k=0; k < N; ++k)
	{
		mCodonFrequencies[k] = 1.0/61.0;
		mCodonFreqSqrt[k] = sqrt(mCodonFrequencies[k]);
		mGoodCodon[k] = true;
	}
	mNumGoodCodons = N;
}

void Forest::mapInternalToBranchIdWalker(const ForestNode* aNode)
{
	std::vector<ForestNode *>::const_iterator in;
	for(in=aNode->mChildrenList.begin(); in != aNode->mChildrenList.end(); ++in)
	{
		ForestNode *m = *in;

		if(m->mInternalNodeId != UINT_MAX) mMapInternalToBranchID[m->mInternalNodeId] = m->mNodeId;

		mapInternalToBranchIdWalker(m);
	}
}


void Forest::addAggressiveReduction(void)
{
	size_t nsites = mRoots.size();

	for(size_t i=0; i < nsites; ++i)
	{
		addAggressiveReductionWalker(&mRoots[i]);
	}
}


void Forest::addAggressiveReductionWalker(ForestNode* aNode)
{
	unsigned int i;
	std::vector<ForestNode *>::iterator in;
	for(in=aNode->mChildrenList.begin(), i=0; in != aNode->mChildrenList.end(); ++in, ++i)
	{
		ForestNode *m = *in;

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


#ifdef NEW_LIKELIHOOD
void Forest::prepareNewReduction(ForestNode* aNode)
{
	size_t nsites = mRoots.size();

	if(aNode)
	{
		std::vector<ForestNode *>::iterator icl;
		for(icl=aNode->mChildrenList.begin(); icl != aNode->mChildrenList.end(); ++icl)
		{
			if((*icl)->mOwnTree == aNode->mOwnTree)
			{
				mNodePresent[(*icl)->mNodeId * nsites + aNode->mOwnTree] = Forest::SITE_EXISTS;
				prepareNewReduction(*icl);
			}
			else
			{
				mNodePresent[(*icl)->mNodeId * nsites + aNode->mOwnTree] = (*icl)->mOwnTree;
			}
		}
	}
	else
	{
		// Initialize the intermediate list
		mNodePresent.assign(mNumBranches*nsites, Forest::SITE_NOT_EXISTS);

		// Visit each site tree
		for(size_t i=0; i < nsites; ++i) prepareNewReduction(&mRoots[i]);

		// Print the sequence of level visits
		std::vector< std::vector<unsigned int> >::iterator inbl;
		for(inbl=mBranchByLevel.begin(); inbl != mBranchByLevel.end(); ++inbl)
		{
			std::vector<unsigned int>::iterator ifn;
			for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
			{
				unsigned int branch_idx = (*ifn);

				size_t begin_idx = 0;
				size_t end_idx   = nsites;
				for(; begin_idx < nsites; ++begin_idx)
				{
					int x = mNodePresent[branch_idx*nsites+begin_idx];
					if(x == Forest::SITE_EXISTS) break;
				}
				if(begin_idx == nsites) throw FastCodeMLFatal("No SITE_EXISTS in mNodePresent");
				for(; end_idx > begin_idx; --end_idx)
				{
					int x = mNodePresent[branch_idx*nsites+end_idx-1];
					if(x == Forest::SITE_EXISTS) break;
				}

				// Count the good elements
				unsigned int cnt = 0;
				for(unsigned int k=begin_idx; k < end_idx; ++k) if(mNodePresent[branch_idx*nsites+k] == Forest::SITE_EXISTS) ++cnt;

				std::cerr << std::setw(2) << branch_idx+1 << ": " << std::setw(4) << begin_idx << '-' << std::setw(4) << end_idx-1 << " (" << cnt << ")" << std::endl;
			}
		}

#if 0
		// Print the sequence of level visits
		std::cerr << std::endl;
		unsigned int level = 1;
		std::vector< std::vector<unsigned int> >::iterator inbl;
		for(inbl=mBranchByLevel.begin(); inbl != mBranchByLevel.end(); ++inbl, ++level)
		{
			std::cerr << level << ": ";

			std::vector<unsigned int>::iterator ifn;
			for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
			{
				std::cerr << (*ifn) << ' ';
			}

			std::cerr << std::endl;
		}

		// TEST
		std::cerr << std::endl;
		for(size_t j=0; j < mNumBranches; ++j)
		{
			bool now = false;
			int first_idx = -1;
			std::cerr << "Branch " << j+1 << std::endl;
			for(size_t k = 0; k < nsites; ++k)
			{
				int x = mNodePresent[j*nsites+k];
				if(x == Forest::SITE_NOT_EXISTS)  std::cerr << '-';
				else if(x == Forest::SITE_EXISTS) std::cerr << 'x';
				else                              std::cerr << x;
				std::cerr << ' ';
			}
			std::cerr << std::endl << std::endl;
		}
#endif
	}
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
		if(aCheckId && aNode->mNodeId != aNodeId) std::cerr << "[site: " << std::setw(5) << aSite << "] aNodeId mismatch " << aNode->mNodeId << " should be: " << aNodeId << std::endl;

		unsigned int id = (aNodeId == UINT_MAX) ? 0 : aNodeId+1;

		std::vector<ForestNode *>::const_iterator icl;
		for(icl=aNode->mChildrenList.begin(); icl != aNode->mChildrenList.end(); ++icl)
		{
			if((*icl)->mOwnTree == aNode->mOwnTree)
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

