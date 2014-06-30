
#ifndef HIGHLEVELCOORDINATOR_H
#define HIGHLEVELCOORDINATOR_H

#include "Forest.h"
#include "CmdLine.h"
#include "WriteResults.h"

/// The rank of the master job
///
static const int MASTER_JOB = 0;

/// Coordinator class for high level parallelization.
/// This class encapsulates MPI usage to parallelize FastCodeML above the maximizer level.
///
/// @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
/// @date 2011-11-22 (initial version)
/// @version 1.1
///

class HighLevelCoordinator
{
public:
	/// Constructor.
	///
	/// @param[in,out] aRgc Pointer to the number of arguments
	/// @param[in,out] aRgv Pointer to the arguments' list
	///
	/// @exception FastCodeMLFatal MPI Failed to initialize
	/// @exception FastCodeMLSuccess To terminate unused worker processes
	///
	HighLevelCoordinator(int* aRgc, char*** aRgv);

	/// Destructor.
	///
	~HighLevelCoordinator();

	/// Starts the high level parallelization of the FastCodeML application
	///
	/// @param[in,out] aForests The filled forests
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	/// @return True if the execution can go parallel at this level.
	///
	bool startWork(std::vector<Forest*> aForests, const CmdLine& aCmdLine);

	/// Is this process the master one?
	///
	/// @return True if this is the master process
	///
	bool isMaster(void) const {return mRank == MASTER_JOB;}

	/// Return the number of MPI processes.
	///
	/// @return The number of MPI processes
	///
	int  numJobs(void) const {return mSize;}

	/// Return the current process rank.
	///
	/// @return The current process rank.
	///
	int getRank(void) const {return mRank;}

private:
	/// The master coordination job
	///
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	/// @exception FastCodeMLFatal Invalid job request found
	///
	void doMaster(const CmdLine& aCmdLine);

	/// The worker high level loop
	///
	/// @param[in,out] aForest The filled forest
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	void doWorker(Forest& aForest, const CmdLine& aCmdLine);


private:
	unsigned int		    mVerbose;				///< The verbose level
	int					    mRank;					///< Rank of the current process (Master has rank == MASTER_JOB)
	int					    mSize;					///< Number of MPI processes

	struct WorkTable;
	std::vector<WorkTable*> mWorkTables;            ///< Management of the work list for each forest.

    std::map<int, int>      mWorkerForestIndexMap;  ///< Lookup worker -> forest / worktable

    /// Register forest - adds a work table for master to track this forest.
	///
	/// @param[in] aForest  The filled forest
	/// @param[in] aCmdLine The command line object
	void registerForest(Forest *aForest, const CmdLine &aCmdLine);

    /// Helper - Worker to forest lookup.
	///
	/// @param[in]  aRank
	/// @return     The index of the forest in aForests as well as the work table that aRank is working on
	int getForestIndexGivenWorker(int aRank) const;

    /// Initializer - Worker -> forest map.
	///
    /// @param[in]     aForests  The forests we are analyzing
    /// @param[in]     aCmdLine  The command line object
	void initWorkerForestIndexMap(
           std::vector<Forest*> aForests,
           const CmdLine &aCmdLine);

    /// Initializer - Worker -> forest map when cmd.mMultipleMode == T and
    /// we have enough processors to treat each branch for H1 separately
    /// @param[in]     aForests    The forests we are analyzing
    /// @param[in,out] aCurWorker  current worker count
    void initWorkerForestIndexMapEnoughProcs(
           std::vector<Forest*> aForests,
           size_t *aCurWorker);

    /// Initializer - Worker -> forest map when cmd.mMultipleMode == T and
    /// we don't have enough processors to treat each branch for H1 separately
    /// @param[in]     aNumWorkers number of worker processes (total)
    /// @param[in]     aForestSize number of forests we are analyzing
    /// @param[in,out] aCurWorker  current worker count
    /// @exception FastCodeMLFatal For less processors available than we have tree/gene pairs
    void initWorkerForestIndexMapTooFewProcs(
           size_t aNumWorkers,
           size_t aForestSize,
           size_t *aCurWorker);

    /// Initializer - Worker -> forest map when cmd.mMultipleMode == T and
    /// we have leftover processes to allocate
    /// @param[in]     aNumWorkers number of worker processes (total)
    /// @param[in]     aForests    Pointers to the forests we are analyzing
    /// @param[in,out] aCurWorker  current worker count
    void initRemainingWorkerForestIndex(
            size_t aNumWorkers,
            std::vector<Forest*> aForests,
            size_t *aCurWorker);

	/// Prints out the results for a given work table
    ///
	/// @param[in]  workTable The pointer to the work table to print
	void printWorkTableResults(WorkTable *aWorkTable) const;

};


#endif
