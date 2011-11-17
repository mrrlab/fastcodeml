
#ifndef FOREST_H
#define FOREST_H

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include "PhyloTree.h"
#include "Genes.h"
#include "ForestNode.h"
#include "TransitionMatrix.h"
#include "TransitionMatrixSet.h"
#include "MatrixSize.h"
#ifdef NEW_LIKELIHOOD
#include "FatVectorTransform.h"
#endif
#include "AlignedAllocator.h"
#include "CodonFrequencies.h"
#include "Types.h"

/// The phylogenetic tree's forest.
/// This class encapsulates the forest of phylogenetic tree that will be used for computing the tree's maximum likelihood
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-02-23 (initial version)
///     @version 1.0
///
class Forest
{
public:
	/// Constructor
	///
	/// @param[in] aVerbose The verbosity level
	///
	explicit Forest(unsigned int aVerbose=0)
		: mNumSites(0), mCodonFreq(0), mNumBranches(0), mVerbose(aVerbose), mNumInternalBranches(0), mMarkedInternalBranch(UINT_MAX)
	{}

	/// Destructor
	///
	~Forest()
	{
		mRoots.clear();
		mNodeNames.clear();
		mBranchLengths.clear();
		mProbs.clear();
		mSiteMultiplicity.clear();		
		mTableInternalToBranchID.clear();
		mDependenciesClasses.clear();	
#ifdef NEW_LIKELIHOOD
		mProbsOut.clear();
		mNodesByLevel.clear();
#endif
	}
	
	/// Build the forest and reduces the subtrees
	///
	/// @param[in] aTree The phylogenetic tree
	/// @param[in] aGenes The corresponding genes
	/// @param[in] aIgnoreFreq Ignore the codon frequencies from file and set them all to 1/61
	///
	void loadTreeAndGenes(const PhyloTree& aTree, const Genes& aGenes, bool aIgnoreFreq=false);

	/// Print the class statistics as: cout << r;
	///
	/// @param[in] aOut Output stream
	/// @param[in] aObj The object to be printed
	///
	/// @return The output stream
	///
	friend std::ostream& operator<< (std::ostream& aOut, const Forest& aObj);

	/// Reduce common subtrees on the whole forest
	///
	/// @param[in] aNoTipPruning If set the branches going to leaves are not pruned
	///
	void reduceSubtrees(bool aNoTipPruning=false);

#ifndef NEW_LIKELIHOOD
	/// Add more aggressive subtree reduction
	///
	/// @param[in] aNode The tree node from which the walker should start (no argument starts from the root)
	///
	void addAggressiveReduction(ForestNode* aNode=0);
#endif

	/// Remove all work data used for reduction
	///
	/// @param[in] aNode The node from which to start. Pass zero to start from the root of all the trees in the forest.
	///
	void cleanReductionWorkingData(ForestNode* aNode=0);

	/// Group trees by dependencies.
	/// First group contains trees with no dependencies.
	/// Second group contains trees that depends only on trees of the first group.
	/// Third group contains trees that depends on first and second groups. And so on.
	///
	/// @param[in] aForceSerial Don't group so the execution is serial
	///
	void groupByDependency(bool aForceSerial);

	/// Compute the log likelihood of the forest given the set of precomputed matrices.
	/// If NEW_LIKELIHOOD is defined, this routine adopts the experimental "Long Vector" approach.
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	///
	void computeLikelihood(const TransitionMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods);

	/// Export the forest in GML format
	///
	/// @param[in] aFilename The filename to be written
	/// @param[in] aCounter Value to substitute \%d or \@d in filename (it is a printf format)
	///
	void exportForest(const char* aFilename, unsigned int aCounter=0) const;

	/// Return the total number of branches
	///
	/// @return The totaal number of branches
	///
	size_t getNumBranches(void) const {return mNumBranches;}

	/// Return the number of internal branches (i.e. the ones that do not connect to leaves)
	///
	/// @return The number of internal branches
	///
	size_t getNumInternalBranches(void) const {return mNumInternalBranches;}

	/// Get the number of sites
	///
	/// @return The number of sites
	///
	size_t getNumSites(void) const {return mNumSites;}

	/// Get the marked internal branch
	///
	/// @return The internal branch index of the branch marked in the tree file. UINT_MAX otherwise.
	///
	size_t getMarkedInternalBranch(void) const {return mMarkedInternalBranch;}

	/// Get site multeplicity values
	///
	/// @return The array of site molteplicities
	///
	const double* getSiteMultiplicity(void) const {return &mSiteMultiplicity[0];}

	/// Set the times (i.e. the branch lengths) from the values read from the tree file
	///
	/// @param[out] aTimes The array with all the tree times
	/// @param[in] aNode The node from which to start (if zero starts from the root)
	///
	void setTimesFromLengths(std::vector<double>& aTimes, const ForestNode* aNode=0) const;

	/// Set the times (i.e. the branch lengths) on the tree from the values read from the times array
	///
	/// @param[in] aTimes The array with all the tree times
	/// @param[out] aNode The node from which to start (if zero starts from the root)
	///
	void setLengthsFromTimes(const std::vector<double>& aTimes, ForestNode* aNode=0);

	/// Change the internal branch identifier for the foreground branch into the corresponding internal branch index.
	///
	/// @param[in] aFgBranch Number of the foreground branch
	///
	/// @return The node index corresponding to the foreground branch
	///
	unsigned int adjustFgBranchIdx(unsigned int aFgBranch) const {return mTableInternalToBranchID[aFgBranch];}

	/// Access the global list of node names.
	///
	/// @return A reference to the list of node names.
	///
	const std::vector<std::string>& getNodeNames(void) const {return mNodeNames;}

#ifdef NEW_LIKELIHOOD
	/// Analyze the forest to prepare the operation to be done to restore the contiguity to the grouped vector approach.
	///
	/// @param[in] aNode The node from which to start. If null then starts from all the trees' roots.
	///
	void prepareNewReduction(ForestNode* aNode=0);

	/// Prepare the data for a forest that has not been reduced
	///
	void prepareNewReductionNoReuse(void);
#endif

#ifdef CHECK_ALGO
	/// Check the forest structure for obvious mistakes (useful only during development)
	///
	/// @param[in] aCheckId If true checks also the node id's (cannot be done after subtree pruning)
	/// @param[in] aNode The node from which to start. If null then start from the roots
	/// @param[in] aSite The site number of the corresponding aNode
	/// @param[in] aNodeId The aNode id
	///
	/// @return The next node ID
	///
	unsigned int checkForest(bool aCheckId=false, const ForestNode* aNode=0, unsigned int aSite=0, unsigned int aNodeId=0) const;
#endif


private:
	/// Reduce the common subtree between two (sub)trees
	///
	/// @param[in] aNode The subtree to be tested (i.e. if it exists in both trees)
	/// @param[in] aNodeDependent The dependent tree (i.e. it could point to subtrees of aNode)
	/// @param[in] aNoTipPruning If set the branches going to leaves are not pruned
	///
	void reduceSubtreesWalker(ForestNode* aNode, ForestNode* aNodeDependent, bool aNoTipPruning);

	/// Check coherence between tree and genes.
	///
	/// @param[in] aTree The phylogenetic tree
	/// @param[in] aGenes The corresponding genes
	///
	/// @exception FastCodeMLFatal Throw exception if the species do not match
	///
	void checkCoherence(const PhyloTree& aTree, const Genes& aGenes) const;

	/// Walker for the exporter
	///
	///	@param[in] aNode The node from which to start
	///	@param[in] aBranchLengths List of all branch lengths
	/// @param[out] aNodeFrom List of starting nodes
	/// @param[out] aNodeTo List of ending nodes
	/// @param[out] aLength Resulting branch lengths to label branches in exported tree
	///
	void exportForestWalker(const ForestNode* aNode,
							const std::vector<double>& aBranchLengths,
							std::vector< std::pair<int, int> >& aNodeFrom,
							std::vector< std::pair<int, int> >& aNodeTo,
							std::vector<double>& aLength) const;

#ifndef NEW_LIKELIHOOD
	/// Walker for the computation of tree likelihood
	///
	/// @param[in] aNode
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	///
	/// @return The vector of codons probabilities at the aNode node
	///
	double* computeLikelihoodWalker(ForestNode* aNode, const TransitionMatrixSet& aSet, unsigned int aSetIdx);
#endif

	/// Walk the tree to fill the mMapInternalToBranchID map.
	///
	///	@param[in] aNode The node from which to start
	/// @param[out] aMapInternalToBranchID Maps internal branch id to branch id
	///
	void mapInternalToBranchIdWalker(const ForestNode* aNode, std::map<unsigned int, unsigned int>& aMapInternalToBranchID);


private:
	///'mNumSites, mCodonFreq, mNumBranches'
	size_t					mNumSites;					///< Number of sites
	const double*			mCodonFreq;					///< Experimental codon frequencies
	size_t					mNumBranches;				///< Total number of branches of the original tree
	std::vector<ForestNode>	mRoots;						///< The roots of the forest's trees. Its length is the number of valid sites
	std::vector<double>		mSiteMultiplicity;			///< Multiplicity of the valid sites
	unsigned int			mVerbose;					///< If greather than zero prints more info
	size_t					mNumInternalBranches;		///< Total number of branches of the original tree
	std::vector<unsigned int>
							mTableInternalToBranchID;	///< Map from internal branch number to branch number
	std::vector< std::vector<unsigned int> >
							mDependenciesClasses;		///< The groups of dependencies between trees

	/// Here are global data that will be removed from the various (site) trees
	std::vector<std::string>
							mNodeNames;					///< List of node names. Zero is the root, then its first child and so on
	std::vector<double>		mBranchLengths;				///< List of branch lengths (read from file or stored here to be exported in the tree file)
	size_t					mMarkedInternalBranch;		///< Number of the internal branch as marked in the tree file

#ifdef NEW_LIKELIHOOD

	/// New loglikelihood computation support
		
	/// The mProbs and mProbsOut layout
	///
	/// [site0][site1][site2]...  [site0][site1][site2]...               each is VECTOR_SLOT bytes long (for which only the first N are significant)
	/// [  set 0                 ][  set 1                 ]...          there are 4 (Nt) sets
	/// [    node 0                                             ]...
	///
	/// site_index = node*(Nt*NumSites*VECTOR_SLOT) + set*(NumSites*VECTOR_SLOT) + site*(VECTOR_SLOT)
	///
	CacheAlignedDoubleVector	mProbs;						///< The concatenation of all the probability vectors for all the nodes and all the classes
	CacheAlignedDoubleVector	mProbsOut;					///< mProbs after multiplication by exp(Qt)
	std::vector< std::vector<ForestNode*> >
								mNodesByLevel;				///< Each level contains a list of pointers to nodes at this level. List start from the root.
	FatVectorTransform			mFatVectorTransform;		///< Compute and manage the transformations to pack the "long vector" based on subtree pruning
#else
	/// Unified array for each branch probability vector
	CacheAlignedDoubleVector	mProbs;						///< The concatenation of all the probability vectors for all the nodes and all the classes
#endif
	std::vector< std::vector<unsigned int> > mTreeDependencies;		///< mTreeDependencies[tj] = <t1 t2 t3> means: tj can be done after: t1 t2 t3
	std::vector< std::vector<unsigned int> > mTreeRevDependencies;	///< mTreeRevDependencies[tj] = <t1 t2 t3> means: tj should be ready before: t1 t2 t3
};

#endif

