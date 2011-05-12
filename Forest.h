
#ifndef FOREST_H
#define FOREST_H

#include <iostream>
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

/// If codon probability over this value, the codon is marked as "good codon".
///
const double GOOD_CODON_THRESHOLD = 1e-100;

/// The phylogenetic tree's forest.
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
	Forest(unsigned int aVerbose=0)
	{
		mVerbose = aVerbose;
		mNumBranches = 0;
		mNumInternalBranches = 0;
		memset(mCodonCount, 0, N*sizeof(unsigned int));
	}

	/// Destructor
	///
	~Forest()
	{
		mRoots.clear();
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

	/// Add more aggressive subtree reduction
	///
	void addAggressiveReduction(void);

	/// Group trees by dependencies.
	/// First group contains trees with no dependencies.
	/// Second group contains trees that depends only on trees of the first group.
	/// Third group contains trees that depends on first and second groups. And so on.
	///
	/// @param[in] aForceSerial Don't group so the execution is serial
	///
	void groupByDependency(bool aForceSerial);

	/// Compute the log likelihood of the tree given the set of precomputed matrices.
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	/// @param[out] aLikelihood Values of the codon probabilities at the tree root
	///
	void computeLikelihood(const TransitionMatrixSet& aSet, unsigned int aSetIdx, std::vector<double>& aLikelihood);

	/// Compute the log likelihood of the tree given the set of precomputed matrices.
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

	/// Return the number of internal branches (i.e. the ones that do not connect to leaves)
	///
	/// @return The number of internal branches
	///
	unsigned int getNumInternalBranches(void) const {return mNumInternalBranches;}

	/// Get the number of sites
	///
	/// @return The number of sites
	///
	unsigned int getNumSites(void) const {return mRoots.size();}

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
	const double* getCodonFrequencies(void) const {return mCodonFrequencies;}

	/// Return the array of square roots of codon frequencies.
	///
	/// @return The pointer to the sqrt codon frequency array
	///
	const double* getSqrtCodonFrequencies(void) const {return mCodonFreqSqrt;}

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


private:
	/// Reduce the common subtree between two trees
	///
	/// @param[in] aRoot1 The first tree
	/// @param[in] aRoot2 The second tree
	///
	void reduceSubtreesWalker(ForestNode* aRoot1, ForestNode* aRoot2);

	/// Add aggresssive reduction to common subtree between two trees already identified
	///
	/// @param[in] aNode The tree node from which the walker should start
	///
	void addAggressiveReductionWalker(ForestNode* aNode);

	/// Check dependencies of a tree with other trees
	///
	/// @param[in] aNode The tree node from which the walker should start
	/// @param[out] aDependency The id of the trees on which this tree depends
	///
	void groupByDependencyWalker(ForestNode* aNode, std::set<unsigned int>& aDependency);

	/// Check coherence between tree and genes.
	/// If the species do not match, throw a FastCodeMLFatal exception
	///
	/// @param[in] aTree The phylogenetic tree
	/// @param[in] aGenes The corresponding genes
	///
	void checkCoherence(const PhyloTree& aTree, const Genes& aGenes) const;

	/// Walker for the exporter
	///
	///	@param[in] aNode The node from which to start
	/// @param[out] aNodeFrom List of starting nodes
	/// @param[out] aNodeTo List of ending nodes
	/// @param[out] aLength Resulting branch lengths to label branches in exported tree
	///
	void exportForestWalker(const ForestNode* aNode,
							std::vector< std::pair<int, int> >& aNodeFrom,
							std::vector< std::pair<int, int> >& aNodeTo,
							std::vector<double>& aLength) const;

	/// Walker for the computation of tree likelihood
	///
	/// @param[in] aNode
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	///
	/// @return The vector of codons probabilities at the aNode node
	///
	double* computeLikelihoodWalker(ForestNode* aNode, const TransitionMatrixSet& aSet, unsigned int aSetIdx);

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

	/// Set the codon frequencies all to 1/61
	///
	void setCodonFrequenciesUnif(void);

	/// Walk the tree to fill the mMapInternalToBranchID map.
	///
	///	@param[in] aNode The node from which to start
	///
	void mapInternalToBranchIdWalker(const ForestNode* aNode);


private:
	std::vector<ForestNode>	mRoots;					///< The roots of the forest's trees. Its length is the number of valid sites
	std::vector<double>		mSiteMultiplicity;		///< Multiplicity of the valid sites
	unsigned int			mVerbose;				///< If greather than zero prints more info
	unsigned int			mNumBranches;			///< Total number of branches of the original tree
	unsigned int			mNumInternalBranches;	///< Total number of branches of the original tree
	double					mCodonFrequencies[N];	///< Experimental codon frequencies
	double					mCodonFreqSqrt[N];		///< Square Root of experimental codon frequencies
	bool					mGoodCodon[N];			///< True if the corresponding codon frequency is not small
	unsigned int			mNumGoodCodons;			///< Number of codons whose frequency is not zero
	unsigned int			mCodonCount[N];			///< Count of codon of each type
	std::map<unsigned int, unsigned int>
							mMapInternalToBranchID;	///< Map from internal branch number to branch number
	std::vector< std::vector<unsigned int> >
							mDependenciesClasses;	///< The groups of dependencies between trees


};

#endif

