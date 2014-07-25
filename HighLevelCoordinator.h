
#ifndef HIGHLEVELCOORDINATOR_H
#define HIGHLEVELCOORDINATOR_H

#include <vector>

#include "ForestGroup.h"
#include "CmdLine.h"
#include "WriteResults.h"
#include "WorkTable.h"



/// Coordinator class for high level parallelization.
/// This class encapsulates MPI usage to parallelize FastCodeML above the maximizer level.
///
/// @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
/// @date 2011-11-22 (initial version)
/// @version 1.1
///

enum JobRequestType
{
	REQ_ANNOUNCE_WORKER,		///< Worker asking for a new job to execute
	REQ_SITE_RESULT,
};
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
	/// @param[in,out] aForestGroup The filled forest group. This class takes control of the ptr.
	/// @param[in]     aCmdLine     The parameters from the command line of the main program
	///
	/// @return True if the execution can go parallel at this level.
	///
	bool startWork(ForestGroup *aForestGroup, const CmdLine& aCmdLine);

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
	static int getRank(void) {return mRank;}

	static int setRank(int aRank) { mRank = aRank; }

    static int getSize(void) {return mSize;}

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
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	void doWorker(const CmdLine& aCmdLine);

    /// Debug function to print the number of processors and forests we have to treat
    ///
    /// @param[in] aCmdLine                  The parameters from the command line of the main program
    /// @param[in] aTotalNumInternalBranches The total number of internal branches to be analyzed
    void checkProcCount(const CmdLine &aCmdLine, size_t aTotalNumInternalBranches) const;

public:
    static int num_threads;

private:
	unsigned int		    mVerbose;				///< The verbose level
	static int					    mRank;			///< Rank of the current process (Master has rank == MASTER_JOB)
	static int					    mSize;			///< Number of MPI processes


    ForestGroup            *mForestGroup;           ///< The group of forests we are analyzing

	/// Prints out the results for a given work table
    ///
	/// @param[in]  aWorkTable The pointer to the work table to print
	void printWorkTableResults(WorkTable *aWorkTable) const;

	/// Traces the MPI calls by outputing the job we sent
	///
	/// @param[in]     aJobType     The job type to print
	/// @param[in]     aBranch      The branch number
	/// @param[in]     aWorker      The worker number
    /// @param[in]     aForestIndex The index of the forest in forest group object
    /// @param[in-out] aOut         The stream to print to
    void printMPITrace(int aJobType, int aBranch, int aWorker, int aForestIndex, std::ostream &aOut) const;
};

template <typename T>
void outputVector(const std::vector<T> &v)
{
    for (int ii= 0; ii< v.size(); ii++)
    {
        std::cout << v[ii];
        if (ii != v.size() - 1)
        {
            std::cout << std::endl;
        }
    }
}



#endif
