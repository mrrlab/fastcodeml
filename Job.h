#ifndef TASK_H
#define TASK_H

#include <stdlib.h>
#include <ostream>

/// Representation of a single MPI Job.
///
/// This class represents a single MPI job, as used in the high level coordinator.
/// In versions 1.2.0 and previous, these were stored in vectors in the high level coordinator.
///
/// The creation of this class stems from the requirement to allow for the -mult mode,
/// and to streamline extensibility of the MPI coordinator.
///
/// @author Kareem Ali - UNIL
/// @date 2014-07-03 (initial version)
/// @version 1.2
///

/// Job type enum. The Worktable routines depend on the order and values of this enum.
enum JobType
{
    INVALID_JOB_TYPE = -1,  ///< Invalid job
	JOB_H0			 = 0,	///< H0 computation
	JOB_H1			 = 1,	///< H1 computation
	JOB_BEB			 = 2,	///< BEB computation (done only if H0 and H1 already done and passing the likelihood ratio test)
	JOB_SHUTDOWN	 = 3,	///< Shutdown the worker
	MAX_JOB_TYPE     = 4    ///< Largest job type
};

/// Invalid forest flag
static const int INVALID_FOREST = -1;

/// Invalid branch flag
static const int INVALID_BRANCH = -1;

struct Job
{
public:
    /// ctor for Job. Ints are used as the data type for consistency with MPI & invalid indicators
	/// @param[in] aForestIndex   the index of the forest in ForestGroup
	/// @param[in] aBranch        the branch this job is regarding
	/// @param[in] aJobType       the type of job this object represents
	/// @param[in] aSizeOtherData the size of the other data we need to pass to the worker thread
	///
	/// @return True if a job has been assigned, false if the job is JOB_SHUTDOWN
	///
    explicit Job(
      int aForestIndex, int aBranch, JobType aJobType, int aSizeOtherData);

    /// Setter for size of other data to send to worker
	/// @param[in] aSizeOtherData the size of the other data to be sent to the work
    void setSizeOtherData(int aSizeOtherData) { mSizeOtherData = aSizeOtherData; }
    /// Getter for size of other data to send to worker
	/// @return   the size of the other data to be sent to the work
    int getSizeOtherData() const {return mSizeOtherData; }

    /// Setter for the job type of the process executing the job
	/// @param[in] aJobType the the job type of the job
    void setJobType(JobType aJobType) { mJobType = aJobType; }
    /// Getter for the job type
	/// @return    the job type
    JobType getJobType() const {return mJobType; }

    /// Getter for the index of the related forest in forest group
	/// @return    the index of the related forest in forest gorup
    int getForestIndex() const {return mForestIndex; }
    /// Getter for the branch
	/// @return    the branch of the forest we are treating with this job
    int getBranch() const {return mBranch; }

    friend std::ostream& operator<<(std::ostream &aOut, const Job &aJob);

    Job& operator=(const Job &that);
    Job(const Job &that);

private:
    int       mForestIndex;    ///<Index of the forest this job is for in the ForestGroup object.
    int       mBranch;         ///<Foreground branch this work item is for, in the context of the forest above
    JobType   mJobType;        ///<The job type of this job
    int       mSizeOtherData;  ///<The size of the other job data, to be sent in the job definition to the worker

    Job();                    ///<Default construction is not allowed.
};
#endif
