
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
#include "FatVectorTransform.h"

/// If codon probability is greater than this value, the codon is marked as "good codon".
///
const double GOOD_CODON_THRESHOLD = 1e-100;

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
	Forest(unsigned int aVerbose=0) : mCodonFrequencies(N, 1./(double)N), mCodonFreqSqrt(N, 1./sqrt((double)N)), mCodonCount(N, 0)
	{
		mVerbose = aVerbose;
		mNumBranches = 0;
		mNumInternalBranches = 0;
		//memset(mCodonCount, 0, N*sizeof(unsigned int));
		mMarkedInternalBranch = UINT_MAX;
		mNumGoodCodons = 0;
		mNumSites = 0;
		//for(int i=0; i < N; ++i) mCodonFrequencies[i] = 1./N;
		//for(int i=0; i < N; ++i) mCodonFreqSqrt[i] = 1./sqrt((double)N);
		for(int i=0; i < N; ++i) mGoodCodon[i] = true;
	}

	/// Destructor
	///
	~Forest()
	{
		mRoots.clear();
		mNodeNames.clear();
		mBranchLengths.clear();
		mProbs.clear();
		mSiteMultiplicity.clear();		
		mMapInternalToBranchID.clear();	
		mDependenciesClasses.clear();	
#ifdef NEW_LIKELIHOOD
		mProbsOut.clear();
		mNodesByLevel.clear();
#endif
		mCodonFrequencies.clear();
		mCodonFreqSqrt.clear();
		mCodonCount.clear();
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
	void reduceSubtrees(void);

#ifndef NEW_LIKELIHOOD
	/// Add more aggressive subtree reduction
	///
	void addAggressiveReduction(void);
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
	void computeLikelihood(const TransitionMatrixSet& aSet, std::vector<double>& aLikelihoods);

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
	size_t getNumSites(void) const {return mRoots.size();}

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

	/// Return codon frequencies
	///
	/// @return The pointer to the codon frequency array (length: 61)
	///
	const double* getCodonFrequencies(void) const {return &mCodonFrequencies[0];}

	/// Return the array of square roots of codon frequencies.
	///
	/// @return The pointer to the sqrt codon frequency array
	///
	const double* getSqrtCodonFrequencies(void) const {return &mCodonFreqSqrt[0];}

	/// Return an indicator array marking codons whose frequency is over GOOD_CODON_THRESHOLD
	///
	/// @return The indicator array (true if the corresponding codon frequency is above GOOD_CODON_THRESHOLD)
	///
	const bool* getGoodCodonFrequencies(void) const {return mGoodCodon;}

	/// Return the count of codons whose frequency is over the threshold.
	///
	/// @return The cound of good codons.
	///
	unsigned int numGoodCodonFrequencies(void) const {return mNumGoodCodons;}

	/// Change the internal branch identifier for the foreground branch into the corresponding internal branch index.
	///
	/// @param[in] aFgBranch Number of the foreground branch
	///
	/// @return The node index corresponding to the foreground branch
	///
	unsigned int adjustFgBranchIdx(unsigned int aFgBranch) const {return mMapInternalToBranchID.find(aFgBranch)->second;}

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
	///
	void reduceSubtreesWalker(ForestNode* aNode, ForestNode* aNodeDependent);

#ifndef NEW_LIKELIHOOD
	/// Add aggresssive reduction to common subtree between two trees already identified
	///
	/// @param[in] aNode The tree node from which the walker should start
	///
	void addAggressiveReductionWalker(ForestNode* aNode);
#endif

	/// Check dependencies of a tree with other trees
	///
	/// @param[in] aNode The tree node from which the walker should start
	/// @param[out] aDependency The id of the trees on which this tree depends
	///
	void groupByDependencyWalker(ForestNode* aNode, std::set<unsigned int>& aDependency);

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

	/// Change the index into the full list of codons (64) into the non-stop codons list used here (61)
	///
	/// @param[in] aId64 The index in the range 0..63 (i.e. from TTT to GGG)
	///
	/// @return The reduced index (range 0..60) or -1 if aId64 is out of range or represents a stop codon.
	///
	int codon64to61(unsigned int aId64) const;

	/// Compute the codon frequency using the F3x4 model
	///
	void setCodonFrequenciesF3x4(void);

	/// Set all codon frequencies to 1/61
	///
	void setCodonFrequenciesUnif(void);

	/// Walk the tree to fill the mMapInternalToBranchID map.
	///
	///	@param[in] aNode The node from which to start
	///
	void mapInternalToBranchIdWalker(const ForestNode* aNode);


private:
	std::vector<ForestNode>	mRoots;						///< The roots of the forest's trees. Its length is the number of valid sites
	std::vector<double>		mSiteMultiplicity;			///< Multiplicity of the valid sites
	unsigned int			mVerbose;					///< If greather than zero prints more info
	size_t					mNumBranches;				///< Total number of branches of the original tree
	size_t					mNumInternalBranches;		///< Total number of branches of the original tree
	//double					mCodonFrequencies[N];		///< Experimental codon frequencies
	std::vector<double>		mCodonFrequencies;			///< Experimental codon frequencies
	//double					mCodonFreqSqrt[N];			///< Square Root of experimental codon frequencies
	std::vector<double>		mCodonFreqSqrt;				///< Square Root of experimental codon frequencies
	bool					mGoodCodon[N];				///< True if the corresponding codon frequency is not small
	//std::vector<bool>		mGoodCodon;					///< True if the corresponding codon frequency is not small
	unsigned int			mNumGoodCodons;				///< Number of codons whose frequency is not zero
	//unsigned int			mCodonCount[N];				///< Count of codon of each type
	std::vector<unsigned int>
							mCodonCount;				///< Count of codon of each type
	std::map<unsigned int, unsigned int>
							mMapInternalToBranchID;		///< Map from internal branch number to branch number
	std::vector< std::vector<unsigned int> >
							mDependenciesClasses;		///< The groups of dependencies between trees
	size_t					mNumSites;					///< Number of sites

	/// Here are global data that will be removed from the various (site) trees
	std::vector<std::string>
							mNodeNames;					///< List of node names. Zero is the root, then its first child and so on
	std::vector<double>		mBranchLengths;				///< List of branch lengths (read from file or stored here to be exported in the tree file)
	size_t					mMarkedInternalBranch;		///< Number of the internal branch as marked in the tree file

#ifdef NEW_LIKELIHOOD

	/// New loglikelihood computation support
		
	/// The mProbs and mProbsOut layout
	///
	/// [site0][site1][site2]...  [site0][site1][site2]...               each is 61 bytes long
	/// [ set 0                  ][ set 1                  ]...          there are 4 (Nt) sets
	/// [   node 0                                             ]...
	///
	/// site_index = node*(Nt*NumSites*N)+set*(NumSites*N)+site*(N)
	///
	std::vector<double>		mProbs;						///< The concatenation of all the probability vectors for all the nodes and all the classes
	std::vector<double>		mProbsOut;					///< mProbs after multiplication by exp(Qt)
	std::vector< std::vector<ForestNode*> >
							mNodesByLevel;				///< Each level contains a list of pointers to nodes at this level. List start from the root.
	FatVectorTransform		mFatVectorTransform;		///< Compute and manage the transformations to pack the "long vector" based on subtree pruning
#else
	/// Unified array for each branch probability vector
	std::vector<double>		mProbs;						///< The concatenation of all the probability vectors for all the nodes and all the classes
#endif
};

#endif

