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
    unsigned int site;

    unsigned int set_idx;
    double val;

    Job() : site(0), set_idx(0) {}

};
#endif
