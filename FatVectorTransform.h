
#ifndef FATVECTORTRANSFORM_H
#define FATVECTORTRANSFORM_H

#include <vector>
#include <utility>

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
	FatVectorTransform() : mNumBranches(0), mNumSites(0) {}

	/// Destructor.
	///
	~FatVectorTransform()
	{
		mNodeStatus.clear();
		mLimits.clear();
		mCopyCmds.clear();
		mReuseCmds.clear();
	}

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
		mLimits.assign(aNumBranches, std::make_pair(0, aNumSites-1));
		mCopyCmds.clear();
		mReuseCmds.clear();
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
		mLimits.assign(aNumBranches, std::make_pair(0, aNumSites-1));
		mCopyCmds.clear();
		mReuseCmds.clear();
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
	/// @param[in] aBranchByLevel List (one for each level starting with the one without dependencies) of lists (the parallel branches)
	///
	void printBranchVisitSequence(const std::vector< std::vector<unsigned int> >& aBranchByLevel) const;

	/// Prints (on sterr) for each branch and each site if it is valid, if it is not present and if takes the value from another site
	///
	void printNodeStatus(void) const;

	/// Compute the commands needed to compact the various lists of commands
	///
	void compactMatrix(void);

	/// Print the lists of generated commands
	///
	void printCommands(void) const;

	/// Get the first index to be used for computation
	///
	/// @param[in] aBranch Specify the value for which branch should be returned.
	///
	/// @return The starting index
	///
	inline unsigned int getLowerIndex(unsigned int aBranch) const {return mLimits[aBranch].first;}

	/// Get the last index to be used for computation
	///
	/// @param[in] aBranch Specify the value for which branch should be returned.
	///
	/// @return The ending index
	///
	inline unsigned int getUpperIndex(unsigned int aBranch) const {return mLimits[aBranch].second;}

	void preCompact(std::vector<double>& aProbs, unsigned int aBranch, unsigned int aSet);


private:
	unsigned int			mNumBranches;			///< The number of branches
	unsigned int			mNumSites;				///< The number of valid sites. The values are:
	std::vector<int>		mNodeStatus;			///< - -2 if the correstponding: Branch -> Site exists
													///< -1 if doesn't exist
													///< The site number from which the value is taken
	enum {
		SITE_EXISTS     = -2,						///< The position (Branch, Site) in mNodePresent exists
		SITE_NOT_EXISTS = -1,						///< The position (Branch, Site) in mNodePresent refers to a not existend node
		SITE_FIRST_NUM  =  0						///< if greather or equal to this the position contains the index from which the value should be copied
	};

	/// Representation of a range to be copied
	struct Range
	{
		Range(unsigned int aFrom, unsigned int aTo, unsigned int aCnt=1) {from = aFrom; to = aTo; cnt = aCnt;}

		unsigned int from;		///< Starting index from which to copy
		unsigned int to;		///< Starting index to which the values should be copied
		unsigned int cnt;		///< How many items (if zero, skip this entry)

		//bool operator<(Range& rhs) { return from < rhs.from; } ///< This is needed for sorting
	};

	typedef std::vector<Range> VectorOfRanges;
	typedef std::vector< std::vector<Range> > VectorOfVectorOfRanges;
	typedef std::vector< std::pair<unsigned int, unsigned int> > VectorOfPairs;

	VectorOfPairs						mLimits;				///< Lower and upper limits for each branch
	VectorOfVectorOfRanges				mCopyCmds;				///< Ranges to be copied to fill the holes
	VectorOfVectorOfRanges				mReuseCmds;				///< Ranges to be reused copying the computed value
};

#endif

