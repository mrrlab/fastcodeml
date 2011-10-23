
#ifndef FATVECTORTRANSFORM_H
#define FATVECTORTRANSFORM_H

#include <vector>
#include <utility>
#include "ForestNode.h"

/// Manipolations on the per-branch probability vector array.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-09-31 (initial version)
///     @version 1.0
///
class FatVectorTransform
{
public:
	/// Constructor.
	///
	FatVectorTransform() : mNumBranches(0), mNumSites(0), mNoTransformations(true) {}

	/// Destructor.
	///
	~FatVectorTransform()
	{
		mNodeStatus.clear();
		mLimits.clear();
		mCopyCmds.clear();
		mReuseCmds.clear();
		mFirstForLevel.clear();
		mBranchByLevel.clear();
	}

	/// Create the dependency list between levels in the tree
	///
	/// @param[in] aNodesByLevel List of lists of pointers to nodes one for each level
	///
	void setBranchDependencies(const std::vector< std::vector<ForestNode*> >& aNodesByLevel);

	/// Initialize the class instance
	///
	/// @param[in] aNumBranches Number of branches
	/// @param[in] aNumSites Number of sites
	///
	inline void initNodeStatus(unsigned int aNumBranches, unsigned int aNumSites)
	{
		mNumBranches = aNumBranches;
		mNumSites = aNumSites;
		mNodeStatus.assign(aNumBranches*aNumSites, FatVectorTransform::SITE_NOT_EXISTS);
		mLimits.assign(aNumBranches, std::make_pair(0, aNumSites));
		mCopyCmds.clear();
		mReuseCmds.clear();
		mNoTransformations = false;
	}
	
	/// Initialize the class instance so can be used if no subtree pruning is present
	///
	/// @param[in] aNumBranches Number of branches
	/// @param[in] aNumSites Number of sites
	///
	inline void initNodeStatusMinimal(unsigned int aNumBranches, unsigned int aNumSites)
	{
		mNumBranches = aNumBranches;
		mNumSites = aNumSites;
		mLimits.assign(aNumBranches, std::make_pair(0, aNumSites));
		mCopyCmds.clear();
		mReuseCmds.clear();
		mNoTransformations = true;
	}

	/// Set the corresponding node as existing for the given site
	///
	/// @param[in] aBranch Branch for which the node has been verified as existing
	/// @param[in] aSite Site for which the node has been verified as existing
	///
	inline void setNodeExists(unsigned int aBranch, unsigned int aSite)
	{
		mNodeStatus[aBranch * mNumSites + aSite] = FatVectorTransform::SITE_EXISTS;
	}

	/// Set the corresponding node as taking its value from another site
	///
	/// @param[in] aBranch Branch for which the node has been checked
	/// @param[in] aSite Site for which the node has been checked
	/// @param[in] aReusedSite Site from which the node takes its value
	///
	inline void setNodeReuses(unsigned int aBranch, unsigned int aSite, unsigned int aReusedSite)
	{
		mNodeStatus[aBranch * mNumSites + aSite] = aReusedSite;
	}

	/// Prints (on stderr) for each branch the first and last valid positions and the valid entries in this range).
	///
	void printCountGoodElements(void) const;

	/// Prints (on stderr) the visit sequence of branches.
	///
	void printBranchVisitSequence(void) const;

	/// Prints (on sterr) for each branch and each site if it is valid, if it is not present and if takes the value from another site
	///
	void printNodeStatus(void) const;

	/// Compute the commands needed to compact the various fat vectors (matrices, one for each branch) of probability vectors
	///
	void compactMatrix(void);

	/// Print the lists of generated commands
	///
	void printCommands(void) const;

	/// Get the first index to be used for computation
	///
	/// @param[in] aBranch Specify which branch should be returned.
	///
	/// @return The starting index in the fat vector
	///
	inline unsigned int getLowerIndex(unsigned int aBranch) const {return mLimits[aBranch].first;}

	/// Get the number of items to be used for computation
	///
	/// @param[in] aBranch Specify which branch should be returned.
	///
	/// @return The count of sites to be used
	///
	inline unsigned int getCount(unsigned int aBranch) const {return mLimits[aBranch].second;}

	/// Compact the fat vector at the leaves.
	///
	/// @param[in,out] aProbs The fat probability vector that will be changed at the leaves level
	///
	void preCompactLeaves(std::vector<double>& aProbs);

	/// Compact the fat vector at a certain level in the tree (TBD)
	///
	void postCompact(std::vector<double>& aStepResults, std::vector<double>& aProbs, unsigned int aLevel, unsigned int aNumSets);


private:
	unsigned int			mNumBranches;			///< The number of branches.
	unsigned int			mNumSites;				///< The number of valid sites.
	std::vector<int>		mNodeStatus;			///< For each (Branch, Site) (idx = branch*NumSites+site) the values are:
													///< - SITE_EXISTS if the corresponding: (Branch, Site) exists
													///< - SITE_NOT_EXISTS if doesn't exist
													///< - The site number from which the value is taken (reused)
	enum {
		SITE_EXISTS     = -2,						///< The position (Branch, Site) in mNodePresent exists
		SITE_NOT_EXISTS = -1,						///< The position (Branch, Site) in mNodePresent refers to a not existend node
		SITE_FIRST_NUM  =  0						///< if greather or equal to this value the position contains the index from which the value should be copied
	};

	/// Representation of a range to be copied
	struct Range
	{
		Range(unsigned int aFrom, unsigned int aTo) {from = aFrom; to = aTo; cnt = 1;}

		unsigned int from;		///< Starting index from which to copy
		unsigned int to;		///< Starting index to which the values should be copied
		unsigned int cnt;		///< How many items (if zero, skip this entry)

		//bool operator<(Range& rhs) { return from < rhs.from; } ///< This is needed for sorting
	};
	struct RangeNoCnt
	{
		RangeNoCnt(unsigned int aFrom, unsigned int aTo) {from = aFrom; to = aTo;}

		unsigned int from;		///< Starting index from which to copy
		unsigned int to;		///< Starting index to which the values should be copied
	};
	typedef std::vector<Range> VectorOfRanges;
	typedef std::vector<VectorOfRanges> VectorOfVectorOfRanges;
	typedef std::vector<RangeNoCnt> VectorOfRangesNoCnt;
	typedef std::vector<VectorOfRangesNoCnt> VectorOfVectorOfRangesNoCnt;
	typedef std::vector< std::pair<unsigned int, unsigned int> > VectorOfPairs;

	VectorOfPairs						mLimits;				///< Lower index and total count for each branch
	VectorOfVectorOfRanges				mCopyCmds;				///< Ranges to be copied to fill the holes (one list for each branch)
	VectorOfVectorOfRangesNoCnt			mReuseCmds;				///< Ranges to be reused copying the computed value (one list for each branch)
	std::vector<bool>					mFirstForLevel;			///< One entry for branch set to true if it is the first entry for its level
	bool								mNoTransformations;		///< If set no transformation will take place (corresponds to no tree prune case)
	std::vector< std::vector<unsigned int> >
										mBranchByLevel;			///< Each level contains a list of branch numbers at this level. List start from the leaves.
	std::vector<unsigned int>			mParentNode;			///< One entry for branch set to the parent node index
};

#endif

