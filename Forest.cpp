
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include "Forest.h"
#include "ForestNode.h"
#include "Exceptions.h"
#include "MathSupport.h"

void Forest::loadTreeAndGenes(const PhyloTree& aTree, const Genes& aGenes, bool aIgnoreFreq)
{
	// Check coherence between tree and genes
	checkCoherence(aTree, aGenes);

	// Number of branches of one tree
	mNumBranches = aTree.getNumBranches();

	// Count the number of unique sites
	unsigned int nsites = aGenes.getNumSites();
	const unsigned int* mult = aGenes.getSiteMultiplicity();
	unsigned int num_leaves = 0;

	// Initialize the count of codon types
	memset(mCodonCount, 0, N*sizeof(unsigned int));

	// Clone tree inside the forest
	mRoots.reserve(nsites);
	mRoots.resize(nsites);
	for(unsigned int j=0; j < nsites; ++j)
	{
		// Create a copy of the tree
		aTree.cloneTree(&mRoots[j], j);

		// Create a list of pointers to leaves
		std::vector<ForestNode*> leaves;
		leaves.clear();
		mRoots[j].pushLeaf(leaves);
		num_leaves = leaves.size();

		// Add codon code to leaves
		std::vector<ForestNode*>::const_iterator il;
		for(il=leaves.begin(); il != leaves.end(); ++il)
		{
			int codon = aGenes.getCodonIdx((*il)->mNodeName, j);
			(*il)->mCodons.push_back(codon);

			// Set leaves probability vector
			memset((*il)->mProb, 0, N*Nt*sizeof(double));
			for(int k=0; k < Nt; ++k) (*il)->mProb[codon+k*N] = 1.0;

			// Count codons
			mCodonCount[codon] += mult[j];
		}

		// Combine going up to the root
		mRoots[j].gatherCodons();
	}

	// Set the number of internal branches
	mNumInternalBranches = mNumBranches - num_leaves;

	// Set the site multeplicity
	mSiteMultiplicity.reserve(nsites);
	for(unsigned int i=0; i < nsites; ++i)
	{
		double d = (double)mult[i];
		mSiteMultiplicity.push_back(d);
	}

	// Set the codon frequencies and related values needed for the eigensolver
	if(aIgnoreFreq)
		setCodonFrequenciesUnif();
	else
		setCodonFrequenciesF3x4();

	// Set the mapping from internal branch number to branch number
	mapInternalToBranchIdWalker(&mRoots[0]);
}


void Forest::reduceSubtrees(void)
{
	unsigned int i, j;
	unsigned int nsites = mRoots.size();

	// Try to merge equal subtrees
	for(i=0; i < nsites-1; ++i)
	{
		for(j=i+1; j < nsites; ++j)
		{
			reduceSubtreesWalker(&mRoots[i], &mRoots[j]);
		}
	}
}


void Forest::reduceSubtreesWalker(ForestNode* aRoot1, ForestNode* aRoot2)
{
	std::vector<ForestNode*>::iterator irn1;
	std::vector<ForestNode*>::iterator irn2;
	unsigned int i;
	unsigned nc = aRoot1->mChildrenList.size();
	for(i=0; i < nc; ++i)
	{
		// If one of the two has been already reduced, do nothing
		if(!aRoot1->mChildSameTree[i] || !aRoot2->mChildSameTree[i]) continue;

		// Check if same subtree
		std::vector<int>::const_iterator ig1;
		std::vector<int>::const_iterator ig2;
		bool are_equal = true;
		for(ig1=aRoot1->mChildrenList[i]->mCodons.begin(), ig2=aRoot2->mChildrenList[i]->mCodons.begin();
			ig1 != aRoot1->mChildrenList[i]->mCodons.end();
			++ig1, ++ig2)
		{
			if(*ig1 != *ig2)
			{
				are_equal = false;
				break;
			}
		}

		if(are_equal)
		{
			aRoot2->mChildSameTree[i] = false;
			delete aRoot2->mChildrenList[i];
			aRoot2->mChildrenList[i] = aRoot1->mChildrenList[i];
		}
	}

	// Recurse
	for(i=0; i < nc; ++i)
	{
		// If one of the two has been already reduced, do nothing
		if(!aRoot1->mChildSameTree[i] || !aRoot2->mChildSameTree[i]) continue;

		reduceSubtreesWalker(aRoot1->mChildrenList[i], aRoot2->mChildrenList[i]);
	}
}


void Forest::groupByDependency(bool aForceSerial)
{
	unsigned int i, j;
	unsigned int nsites = mRoots.size();
	if(aForceSerial)
	{
		std::vector<unsigned int> v;

		for(i=0; i < nsites; ++i) v.push_back(i);
		mDependenciesClasses.push_back(v);

		return;
	}

	std::vector< std::set<unsigned int> > dependencies;
	std::vector<bool> done;
	done.resize(nsites, false);
	std::vector<bool> prev;
	prev.resize(nsites, false);

	// Collect dependencies
	for(i=0; i < nsites; ++i)
	{
		std::set<unsigned int> dep;
		dep.clear();
		groupByDependencyWalker(&mRoots[i], dep);
		dependencies.push_back(dep);
	}

	// Trees without dependencies
	std::vector<unsigned int> v;
	std::vector< std::set<unsigned int> >::iterator is;
	for(is=dependencies.begin(),i=0; is != dependencies.end(); ++is,++i)
	{
		if(is->empty())
		{
			v.push_back(i);
			done[i] = true;
			prev[i] = true;
		}
	}
	mDependenciesClasses.push_back(v);
	if(mVerbose >= 1) std::cerr << "Trees in class  0: " << std::setw(3) << v.size() << std::endl;

	for(j=1;; ++j)
	{
		v.clear();
		bool all_done = true;
		for(is=dependencies.begin(),i=0; is != dependencies.end(); ++is,++i)
		{
			if(done[i]) continue;

			all_done = false;
			bool all = true;
			std::set<unsigned int>::iterator iv;
			for(iv=is->begin(); iv != is->end(); ++iv)
			{
				if(!prev[*iv])
				{
					all = false;
					break;
				}
			}
			if(all)
			{
				v.push_back(i);
				done[i] = true;
			}
		}
		if(all_done) break;
		mDependenciesClasses.push_back(v);
		prev = done;
		if(mVerbose >= 1) std::cerr << "Trees in class " << std::setw(2) << j << ": " << std::setw(3) << v.size() << std::endl;
	}
}


void Forest::groupByDependencyWalker(ForestNode* aNode, std::set<unsigned int>& aDependency)
{
	unsigned int i;
	unsigned nc = aNode->mChildrenList.size();
	for(i=0; i < nc; ++i)
	{
		if(aNode->mChildSameTree[i])
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
	unsigned int i;

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
			aObj.mRoots[i].print(aOut);
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
		 exportForestWalker(&(*ifn), node_from, node_to, branch_length);
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
	for(unsigned int i=0; i < node_to.size(); ++i)
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
				net << "   label \"" << iv->second << "\"\n";
				net << "   type 1\n";
			}
			net << "]\n";
		}

		for(unsigned int i=0; i < node_from.size(); ++i)
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
		aLength.push_back((*ifn)->mBranchLength);

		if(your_tree_id == my_tree_id) exportForestWalker(*ifn, aNodeFrom, aNodeTo, aLength);
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
	else       aTimes[aNode->mNodeId] = aNode->mBranchLength;

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
		aNode->mBranchLength = aTimes[aNode->mNodeId];

		for(ifnp=aNode->mChildrenList.begin(); ifnp != aNode->mChildrenList.end(); ++ifnp)
		{
			setLengthsFromTimes(aTimes, *ifnp);
		}
	}
}

void Forest::computeLikelihood(const TransitionMatrixSet& aSet, unsigned int aSetIdx, std::vector<double>& aLikelihood)
{
#if 0
	unsigned int num_sites = mRoots.size();
	aLikelihood.reserve(num_sites);
	//aLikelihood.resize(num_sites, 1.0);
	for(unsigned int site=0; site < num_sites; ++site)
	{
		// Compute likelihood array at the root of one tree
		double* g = computeLikelihoodWalker(&mRoots[site], aSet, aSetIdx);
		//aLikelihood[site] = dot(mCodonFrequencies, g);
		aLikelihood.push_back(dot(mCodonFrequencies, g));
	}
#else
	unsigned int num_sites = mRoots.size();
	aLikelihood.resize(num_sites, 1.0);

	std::vector< std::vector<unsigned int> >::iterator ivs;
	for(ivs=mDependenciesClasses.begin(); ivs != mDependenciesClasses.end(); ++ivs)
	{
		int len = ivs->size();
#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(len, ivs, aLikelihood, aSet, aSetIdx)
#else
		#pragma omp parallel for default(none) shared(len, ivs, aSetIdx)
#endif
		for(int i=0; i < len; ++i)
		{
			// Compute likelihood array at the root of one tree
			unsigned int idx = ivs->at(i);
			double* g = computeLikelihoodWalker(&mRoots[idx], aSet, aSetIdx);
			aLikelihood[idx] = dot(mCodonFrequencies, g);
		}
	}
#endif
}


void Forest::computeLikelihood(const TransitionMatrixSet& aSet, std::vector<double>& aLikelihoods)
{
	unsigned int num_sets = aSet.size();
	unsigned int num_sites = mRoots.size();
	aLikelihoods.resize(num_sets*num_sites, 1.0);

	std::vector< std::vector<unsigned int> >::iterator ivs;
	for(ivs=mDependenciesClasses.begin(); ivs != mDependenciesClasses.end(); ++ivs)
	{
		int len = ivs->size()*num_sets;

#ifdef _MSC_VER
		#pragma omp parallel for default(none) shared(aSet, len, ivs, num_sets, num_sites, aLikelihoods)
#else
		#pragma omp parallel for default(none) shared(len, ivs, num_sets, num_sites)
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
		if(aNode->mChildSameTree[idx])
		{
			if(first)
			{
				aSet.doTransition(aSetIdx, m->mNodeId, computeLikelihoodWalker(m, aSet, aSetIdx), aNode->mProb+N*aSetIdx);
				if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mOtherTreeProb[idx]+N*aSetIdx, aNode->mProb+N*aSetIdx, N*sizeof(double));
				first = false;
			}
			else
			{
				double temp[N];
				double* x = aNode->mOtherTreeProb[idx] ? aNode->mOtherTreeProb[idx]+N*aSetIdx : temp;
				aSet.doTransition(aSetIdx, m->mNodeId, computeLikelihoodWalker(m, aSet, aSetIdx), x);
				for(int i=0; i < N; ++i) aNode->mProb[i+N*aSetIdx] *= x[i];
			}
		}
		else
		{
			if(first)
			{
				if(aNode->mOtherTreeProb[idx]) memcpy(aNode->mProb+N*aSetIdx, aNode->mOtherTreeProb[idx]+N*aSetIdx, N*sizeof(double));
				else aSet.doTransition(aSetIdx, m->mNodeId, m->mProb+N*aSetIdx, aNode->mProb+N*aSetIdx);
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
					aSet.doTransition(aSetIdx, m->mNodeId, m->mProb+N*aSetIdx, temp);
					x = temp;
				}
				for(int i=0; i < N; ++i) aNode->mProb[i+N*aSetIdx] *= x[i];
			}
		}
	}

	return aNode->mProb+N*aSetIdx;
}


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

#if 0
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

#if 0
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

#if 0
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
	unsigned int nsites = mRoots.size();

	for(unsigned int i=0; i < nsites; ++i)
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

		if(aNode->mChildSameTree[i])
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
