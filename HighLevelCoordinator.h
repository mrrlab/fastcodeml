
#ifndef HIGHLEVELCOORDINATOR_H
#define HIGHLEVELCOORDINATOR_H

#include <vector>
#include "Forest.h"

/// The rank of the master job
///
static const int MASTER_JOB = 0;

/// Coordinator class for high level parallization.
/// This class encapsulates MPI usage to parallelize FastCodeML above the maximizer level
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-11-22 (initial version)
///     @version 1.0
///

class HighLevelCoordinator
{
public:
	/// Constructor
	///
	/// @param[in,out] aRgc Pointer to the number of arguments
	/// @param[in,out] aRgv Pointer to the arguments' list
	///
	HighLevelCoordinator(int* aRgc, char*** aRgv);

	/// Destructor
	///
	~HighLevelCoordinator();

	/// Starts the high level parallelization of the FastCodeMP application
	///
	/// @param[in,out] aForest The filled forest
	/// @param[in] aSeed The random number generator seed 
	/// @param[in] aVerbose The verbose level
	/// @param[in] aNoMaximization If true non likelihood maximization takes place
	/// @param[in] aTimesFromFile If true the initial branch length values are retrieved from the tree file
	/// @param[in] aOptimizationAlg The maximization algorithm to be used
	///
	/// @return True if the execution can go parallel at this level.
	///
	bool startWork(Forest& aForest, unsigned int aSeed, unsigned int aVerbose=0, bool aNoMaximization=false, bool aTimesFromFile=true, unsigned int aOptimizationAlgo=0);

	/// Is this process the master one?
	///
	/// @return True if this is the master process
	///
	bool isMaster(void) const {return mRank == MASTER_JOB;}

	/// Return the number of MPI processes
	///
	/// @return The number of MPI processes
	///
	int  numJobs(void) const {return mSize;}


private:
	/// The master coordination job
	///
	void doMaster(void);

	/// The worker high level loop
	///
	/// @param[in,out] aForest The filled forest
	/// @param[in] aSeed The random number generator seed 
	/// @param[in] aVerbose The verbose level
	/// @param[in] aNoMaximization If true non likelihood maximization takes place
	/// @param[in] aTimesFromFile If true the initial branch length values are retrieved from the tree file
	/// @param[in] aOptimizationAlg The maximization algorithm to be used
	///
	void doWorker(Forest& aForest, unsigned int aSeed, bool aNoMaximization, bool aTimesFromFile, unsigned int aOptimizationAlgo);


private:
	unsigned int		mVerbose;				///< The verbose level
	int					mRank;					///< Rank of the current process (Master as rank == MASTER_JOB)
	int					mSize;					///< Number of MPI processes
	unsigned int		mNumInternalBranches;	///< Number of internal branches (i.e. the ones that can be foreground branch)

	struct WorkTable;
	WorkTable*			mWorkTable;				///< Management of the work list
};




#endif
