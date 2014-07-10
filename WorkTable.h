#ifndef WORKTABLE_H
#define WORKTABLE_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <queue>
#include "Job.h"
#include "ForestGroup.h"
#include "WriteResults.h"

/// Table of work to be done in MPI mode.
///
/// This class represents a work table, alternately, a list of jobs to execute by the
/// MPI processes.
///
/// @author Kareem Ali - UNIL
/// @date 2014-07-03 (initial version)
/// @version 1.2
///

/// Flag for invalid process rank
static const int INVALID_RANK = -1;

/// Scaling value for transforming probabilities into integers for transmission
static const double PROB_SCALING = 1.0e9;

/// Results for one branch
///
struct ResultSet
{
    double				mLnl[2];				///< Likelihood values for H0 and H1
    std::vector<double> mHxVariables[2];		///< Variables for H0 and H1
    std::vector<int>    mPositiveSelSites;		///< Sites (if any) under positive selection
    std::vector<double>	mPositiveSelProbs;		///< Corresponding probabilities
    bool                mSkipped;               ///< The results for this branch were not computed

    /// Default constructor.
    ///
    ResultSet()
    {
        mLnl[0] = mLnl[1] = -DBL_MAX;
        mSkipped=true;
    }
};

class WorkTable
{
public:
	/// Constructor
	///
	/// @param[in] aForestGroup The forest group we are working on
	/// @param[in] aCmdLine     The command line parameters
	///
	explicit WorkTable(const ForestGroup &aForestGroup, const CmdLine &aCmdLine);

	/// Get the next job to be executed. If no more jobs then set aJob to shutdown and return false
	///
	/// @param[out] aJob The job request: [0] is set to the kind of job (JOB_H0, JOB_H1, JOB_BEB, JOB_SHUTDOWN);
	///									  [1] to the fg branch number (or zero for JOB_SHUTDOWN);
	///                                   [2] zero or the number of variables for a JOB_BEB or JOB_H0 requests
	///                                   [3] the index of the forest in the forest group object [0, ... number of forests - 1]
	///
	/// @return True if a job has been assigned, false if the job is JOB_SHUTDOWN
	///
	bool getNextJob(int* aJob);

    /// Record the results of one of the hypothesis tests (H0 or H1)
    ///
    /// @param[in] aJob     the Job these results related to
    /// @param[in] aResults The results of the job
    /// @param[in] aWorker  The id of the worker who did the job
    /// @param[in] aVerbose Verbosity - controls if we output the lnl value to cout
    ///
    void recordResultHyp(const Job &aJob, const std::vector<double> &aResults, int aWorker, unsigned int aVerbose);

    /// Record the results of the BEB
    ///
    /// @param[in] aJob     the Job these results related to
    /// @param[in] aResults The results of the job
    /// @param[in] aWorker  The id of the worker who did the job
    /// @param[in] aVerbose Verbosity - controls if we output the lnl value to cout
    ///
    void recordResultBEB(const Job &aJob, const std::vector<int> &aResults, int aWorker, unsigned int aVerbose);

	/// Mark the job assigned to the aRank worker as finished
	///
	/// @param[in] aRank The rank of the current worker
	///
	/// @return The identifier of the finished job (it is branch*JOBS_PER_BRANCH+job_type)
	///
	/// @exception FastCodeMLFatal If job not found
	///
	void queueNextJob(const Job &aJob, unsigned int aVerbose, const bool aInitFromH1);

	/// Print completed branches.
	/// This routine should be called after markJobFinished().
	///
	/// @param[in] aIdx  The identifier of the finished job (it is branch*JOBS_PER_BRANCH+job_type) as returned by markJobFinished().
	///
	void printFinishedBranch(const Job &aJob) const;

	/// Print the optimized variables.
	///
	/// @param[in] aBranch The branch to be printed.
	/// @param[in] aHyp The hypothesis results to be printed
	/// @param[in] aOut The stream on which the print should be done.
	///
	void printVariables(size_t aBranch, unsigned int aHyp, std::ostream& aOut,
        const ResultSet &aResultSet) const;

    /// Get a pointer to the result set
    ///
    /// @param[in] aForestIndex    Point to the forest in the forest group
    /// @return                    Pointer to the vector of result sets
	std::vector<ResultSet>* getResultSetPtr(size_t aForestIndex);

    /// Print the result set to the given stream
    ///
    /// @param[in-out] aOut       The output stream
    /// @param[in]     aCmdLine   The command line object
    void printResults(std::ostream &aOut, const CmdLine &aCmdLine) const;

private:
    std::queue<Job> mJobQueue;                          ///<The queued jobs to be executed

    std::map<size_t, std::vector<ResultSet> > mResults; ///<The results (forest index -> result set)

    std::map<size_t, WriteResults> mWriteResults;       ///<The write results objects (forest index > write result set)

    bool mStaticQueue;                                  ///<This variable is true when aHypothesisToDo is populated

    /// Add a forest to the work table
    ///
    /// @param[in] aForest      Reference to the forest
    /// @param[in] aForestIndex The index of the forest in the forest group
    /// @param[in] aCmdLine     The command line object
    ///
    void addForest(const Forest &aForest, size_t aForestIndex, const CmdLine &aCmdLine);

    /// Print a result set to the given output
    ///
    /// @param[in-out] aOut       The output stream
    /// @param[in]     aResultSet The result set object
    ///
    void printResultSet(std::ostream &aOut, const std::vector<ResultSet> &aResultSet) const;

	/// Default ctor is blocked
	WorkTable();

	/// No equality operator
	WorkTable& operator=(const WorkTable &that);

	/// No copying is allowed
	WorkTable(const WorkTable &that);
};

#endif
