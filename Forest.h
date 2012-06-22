
#ifndef FOREST_H
#define FOREST_H

#include <fstream>
#include <vector>
#include <utility>
#include "PhyloTree.h"
#include "Genes.h"
#include "ForestNode.h"
#include "TransitionMatrix.h"
#include "ProbabilityMatrixSet.h"
#include "MatrixSize.h"
#ifdef NEW_LIKELIHOOD
#include "FatVectorTransform.h"
#endif
#include "CodonFrequencies.h"
#include "Types.h"
#include "ForestExport.h"
#include "DAGScheduler.h"

/// Global scaling factor
static const double GLOBAL_SCALING_FACTOR = 1.0e3;

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
							: mNumSites(0), mCodonFreq(0), mNumBranches(0), mVerbose(aVerbose), mNumInternalBranches(0), mMarkedInternalBranch(UINT_MAX) {}

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
#ifdef NEW_LIKELIHOOD
		mProbsOut.clear();
		mNodesByLevel.clear();
#endif
#ifdef NON_RECURSIVE_VISIT
		mVisitTree.clear();
		mVisitTreeParents.clear();
#endif
	}
	
	/// Build the forest and reduces the subtrees
	///
	/// @param[in] aTree The phylogenetic tree
	/// @param[in] aGenes The corresponding genes
	/// @param[in] aCodonFrequencyModel Model to be used to compute the codon empirical frequencies.
	///
	void loadTreeAndGenes(const PhyloTree& aTree,
						  const Genes& aGenes,
						  CodonFrequencies::CodonFrequencyModel aCodonFrequencyModel);

	/// Print the class statistics as: cout << r;
	///
	/// @param[in] aOut Output stream
	/// @param[in] aForest The forest to be printed
	///
	/// @return The output stream
	///
	friend std::ostream& operator<< (std::ostream& aOut, const Forest& aForest);

	/// Reduce common subtrees on the whole forest
	///
	void reduceSubtrees(void);

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

	/// Group trees by dependencies try to balance their parallel execution.
	/// First group contains trees with no dependencies.
	/// Second group contains trees that depend only on trees of the first group.
	/// Third group contains trees that depend on first and second groups. And so on.
	///
	/// @param[in] aForceSerial Don't group so the execution is serial
	///
	void prepareDependencies(bool aForceSerial);

#if !defined(NON_RECURSIVE_VISIT) && !defined(NEW_LIKELIHOOD)
	/// Compute likelihood visiting the trees in a non-recursive way
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	/// @param[in] aHyp The hypothesis to be computed (H0: 0; H1: 1)
	///
	void computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int aHyp);
#endif

#ifdef NON_RECURSIVE_VISIT
	/// Prepare the list of threading pointers for non-recursive trees visit
	///
	void prepareNonRecursiveVisit(void);

	/// Compute likelihood visiting the trees in a non-recursive way
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	/// @param[in] aHyp The hypothesis to be computed (H0: 0; H1: 1)
	///
	void computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int aHyp);
#endif

#ifdef NEW_LIKELIHOOD
	/// Compute the log likelihood of the forest given the set of precomputed matrices.
	/// If NEW_LIKELIHOOD is defined, this routine adopts the experimental "Long Vector" approach.
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	/// @param[in] aHyp The hypothesis to be computed (H0: 0; H1: 1) (currently ignored)
	///
	void computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int /*aHyp*/);
#endif

	/// Export the forest as graph file
	///
	friend class ForestExport;

	/// Return the total number of branches
	///
	/// @return The total number of branches
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

	/// Get site multeplicity values.
	///
	/// @return Reference to the array of site multiplicities
	///
	const std::vector<double>& getSiteMultiplicity(void) const {return mSiteMultiplicity;}

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
	unsigned int adjustFgBranchIdx(size_t aFgBranch) const {return mTableInternalToBranchID[aFgBranch];}

	/// Access the global list of node names.
	///
	/// @return A reference to the list of node names.
	///
	const std::vector<std::string>& getNodeNames(void) const {return mNodeNames;}

#ifdef NEW_LIKELIHOOD
	/// All the preparatory steps needed for the Fat Vector approach.
	///
	void postLoad(void);

	/// Analyze the forest to prepare the operation to be done to restore the contiguity to the grouped vector approach.
	///
	/// @param[in] aNode The node from which to start. If null then starts from all the trees' roots.
	///
	void prepareNewReduction(ForestNode* aNode=0);

	/// Prepare the data for a forest that has not been reduced
	///
	void prepareNewReductionNoReuse(void);
#endif

	/// Load the forest into a DAG
	///
	/// @param[in] aDAG The DAG structure to be loaded
	///
	void loadForestIntoDAG(DAGScheduler& aDAG) const;

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
	/// Group trees by dependencies.
	/// First group contains trees with no dependencies.
	/// Second group contains trees that depend only on trees of the first group.
	/// Third group contains trees that depend on first and second groups. And so on.
	///
	/// @param[in] aForceSerial Don't group so the execution is serial
	///
	void groupByDependency(bool aForceSerial);

	/// Balance the groups so they have a number of elemnent multiple of the number of available threads.
	///
	/// @param[in] aForceSerial Don't group so the execution is serial
	/// @param[in] aHyp The hypothesis to be computed (H0: 0; H1: 1)
	/// @param[in] aGreedy Try to move as much as possible
	///
	/// @return True if the grouping changed
	///
	bool balanceDependenciesClassesAndTrees(bool aForceSerial, int aHyp, bool aGreedy);

	/// Print the size of each class
	///
	void printDependenciesClassesAndTrees(void);

	/// Measure the number of branches to be computed at each site
	///
	/// @param[out] aEffort Effort per site
	///
	void measureEffort(std::vector<unsigned int>& aEffort);

	/// For each group print the total effort per thread. using the new Site/Class structure
	///
	/// @param[in] aEffort Effort per site
	/// @param[in] aHyp The hypothesis to consider (could be 0 or 1)
	///
	void printEffortByGroup(const std::vector<unsigned int>& aEffort, unsigned int aHyp);

	/// For each group print the total effort per thread.
	///
	/// @param[in] aEffort Effort per site
	/// @param[in] aHyp The hypothesis to consider (could be 0 or 1)
	///
	/// @return The sum of all maxima velues per class (a crude approximation of the runtime value)
	///
	unsigned int totalEffort(const std::vector<unsigned int>& aEffort, unsigned int aHyp);

	/// Balance effort inside each class.
	///
	/// @param[in] aEffort Effort per site
	/// @param[in] aHyp The hypothesis to consider (could be 0 or 1)
	///
	void balanceEffort(const std::vector<unsigned int>& aEffort, unsigned int aHyp);

	/// Reduce the common subtree between two (sub)trees
	///
	/// @param[in] aNode The subtree to be tested (i.e. if it exists in both trees)
	/// @param[in] aNodeDependent The dependent tree (i.e. it could point to subtrees of aNode)
	///
	void reduceSubtreesWalker(ForestNode* aNode, ForestNode* aNodeDependent);

#if !defined(NON_RECURSIVE_VISIT) && !defined(NEW_LIKELIHOOD)
	/// Walker for the computation of tree likelihood
	///
	/// @param[in] aNode The node from which the visit should start
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	///
	/// @return The vector of codons probabilities at the aNode node
	///
	double* computeLikelihoodsWalkerTC(ForestNode* aNode, const ProbabilityMatrixSet& aSet, unsigned int aSetIdx);
#endif

#ifdef NON_RECURSIVE_VISIT
	/// Walker to prepare the non recursive visit list
	///
	/// @param[in] aNode The current node to be visited
	/// @param[in] aParentNode Parent node for aNode
	/// @param[in] aSite The current site
	/// @param[in,out] aVisitList The list of nodes to be visited in the order of visit.
	/// @param[in,out] aParentList The corresponding parent nodes
	///
	void prepareNonRecursiveVisitWalker(ForestNode* aNode,
										ForestNode* aParentNode,
										unsigned int aSite, 
										std::vector<ForestNode*>& aVisitList, 
										std::vector<ForestNode*>& aParentList);

	/// Walker for the computation of tree likelihood
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	/// @param[in] aSiteIdx The site under computation
	///
	void computeLikelihoodsWalkerNR(const ProbabilityMatrixSet& aSet, unsigned int aSetIdx, unsigned int aSiteIdx);
#endif

	/// Walk the tree to fill the mMapInternalToBranchID map.
	///
	///	@param[in] aNode The node from which to start
	/// @param[out] aMapInternalToBranchID Maps internal branch id to branch id
	///
	void mapInternalToBranchIdWalker(const ForestNode* aNode, std::map<unsigned int, unsigned int>& aMapInternalToBranchID);

	/// Walker for the loadForestIntoDAGroutine.
	///
	/// @param[in,out] aDAG The DAG to be built.
	/// @param[in] aNode The current node of the forest to be visited.
	///
	void loadForestIntoDAGWalker(DAGScheduler& aDAG, const ForestNode* aNode) const;

private:
	size_t					mNumSites;					///< Number of sites
	const double*			mCodonFreq;					///< Experimental codon frequencies
	size_t					mNumBranches;				///< Total number of branches of the original tree
	std::vector<ForestNode>	mRoots;						///< The roots of the forest's trees. Its length is the number of valid sites
	std::vector<double>		mSiteMultiplicity;			///< Multiplicity of the valid sites
	unsigned int			mVerbose;					///< If greather than zero prints more info
	size_t					mNumInternalBranches;		///< Total number of branches of the original tree
	std::vector<unsigned int>
							mTableInternalToBranchID;	///< Map from internal branch number to branch number
	typedef std::vector< std::vector<std::pair<unsigned int, unsigned int> > >
							ListDependencies;			///< List (each list depends on the previous) of list (sites to be executed in parallel) of pairs (site, site class)
	ListDependencies mDependenciesClassesAndTrees[2];	///< The groups of dependencies between trees (The two entries are for the two hypothesis)

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
	std::vector< std::vector<unsigned int> >	mTreeDependencies;		///< mTreeDependencies[tj] = [t1 t2 t3] means: tj can be done after: t1 t2 t3
	std::vector< std::vector<unsigned int> >	mTreeRevDependencies;	///< mTreeRevDependencies[tj] = [t1 t2 t3] means: tj should be ready before: t1 t2 t3

#ifdef NON_RECURSIVE_VISIT
	std::vector< std::vector<ForestNode*> >		mVisitTree;				///< List of pointers to tree nodes (a list per site) in the non-recursive visit order
	std::vector< std::vector<ForestNode*> >		mVisitTreeParents;		///< List of parent pointers for the corresponding nodes in the mVisitTree
#endif
};

#endif

