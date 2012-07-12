
#ifdef USE_MPI

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <mpi.h>
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic pop
#endif

#include <iostream>
#include <iomanip>
#ifndef VTRACE
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

#include "HighLevelCoordinator.h"
#include "BranchSiteModel.h"
#include "Exceptions.h"
#include "BayesTest.h"
#include "VerbosityLevels.h"

enum JobStatus
{
	JOB_WAITING,			///< Not yet assigned
	JOB_ASSIGNED,			///< Assigned but not finished
	JOB_COMPLETED,			///< Job completed
	JOB_SKIP				///< Job should not execute due to preconditions
};

enum JobType
{
	JOB_H0			= 0,	///< H0 computation
	JOB_H1			= 1,	///< H1 computation
	JOB_BEB			= 2,	///< BEB computation (done only if H0 and H1 already done and passing the likelihood ratio test)
	JOB_SHUTDOWN	= 3		///< Shutdown the worker
};

enum MessageType
{
	MSG_WORK_REQUEST,		///< Worker asking for a new job to execute
	MSG_NEW_JOB,			///< New job from the master
	MSG_GET_RESULTS			///< Get step results from worker
};

enum JobRequestType
{
	REQ_ANNOUNCE_WORKER,		///< Worker asking for a new job to execute
	REQ_HX_RESULT,				///< The master gets the results of a H0 or H1 job
	REQ_BEB_RESULT				///< The master gets the results of a BEB job
};

/// Scaling value for transforming probabilities into integers for transmission
static const double PROB_SCALING = 1e9;

/// For each internal branch there are three jobs to be executed: H0, H1 and BEB
static const int JOBS_PER_BRANCH = 3;


/// Table of work to be done and intermediate results.
///
struct HighLevelCoordinator::WorkTable
{
	/// Results for one branch
	///
	struct ResultSet
	{
		double				mLnl[2];				///< Likelihood values for H0 and H1
		std::vector<double> mHxVariables[2];		///< Variables for H0 and H1
		std::vector<int>    mPositiveSelSites;		///< Sites (if any) under positive selection
		std::vector<double>	mPositiveSelProbs;		///< Corresponding probabilities
	};

	/// Constructor
	///
	/// @param[in] aNumInternalBranches Number of internal branches that can be marked as foreground branch.
	///
	WorkTable(size_t aNumInternalBranches) :
					mNumInternalBranches(aNumInternalBranches),
					mJobStatus(aNumInternalBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mWorkList(aNumInternalBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mResults(aNumInternalBranches) {}

	size_t					mNumInternalBranches;	///< Number of internal branches that can be marked as foreground branch.
	std::vector<int>		mJobStatus;				///< Corresponding step status
	std::vector<int>		mWorkList;				///< Who is doing this step
	std::vector<ResultSet>	mResults;				///< Results for each branch

	/// Get the next job to be executed. If no more jobs then set aJob to shutdown and return false
	///
	/// @param[out] aJob The job request: [0] is set to the kind of job (JOB_H0, JOB_H1, JOB_BEB, JOB_SHUTDOWN) and [1] to the branch number (or zero for JOB_SHUTDOWN)
	/// @param[in] aRank The current worker rank
	///
	/// @return True if a job has been assigned, false if the job is JOB_SHUTDOWN
	///
	bool getNextJob(int* aJob, int aRank);

	/// Mark the job assigned to the aRank worker as finished
	///
	/// @param[in] aRank The rank of the current worker
	///
	/// @return The identifier of the finished job (it is branch*JOBS_PER_BRANCH+job_type)
	///
	/// @exception FastCodeMLFatal If job not found
	///
	int markJobFinished(int aRank);

	/// Check that all jobs have been processed
	///
	/// @exception FastCodeMLFatalNoMsg If jobs still pending
	///
	void checkAllJobsDone(void) const;
};


bool HighLevelCoordinator::WorkTable::getNextJob(int* aJob, int aRank)
{
	// Assign all H0 jobs, then all H1 jobs
	for(int kind=JOB_H0; kind <= JOB_H1; ++kind)
	{
		for(size_t branch=0; branch < mNumInternalBranches; ++branch)
		{
			size_t idx = branch*JOBS_PER_BRANCH+kind;
			if(mJobStatus[idx] == JOB_WAITING)
			{
				aJob[0] = kind;
				aJob[1] = static_cast<int>(branch);
				mJobStatus[idx] = JOB_ASSIGNED;
				mWorkList[idx]  = aRank;
				return true;
			}
		}
	}

	// Assign BEB jobs if possible
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		// The BEB job should be pending
		if(mJobStatus[branch*JOBS_PER_BRANCH+JOB_BEB] != JOB_WAITING) continue;

		// The corresponding H0 and H1 jobs should be completed
		if(mJobStatus[branch*JOBS_PER_BRANCH+JOB_H0] != JOB_COMPLETED || mJobStatus[branch*JOBS_PER_BRANCH+JOB_H1] != JOB_COMPLETED) continue;

		// If the previous results do not pass the LRT, then skip the BEB computation
		if(!BranchSiteModel::performLRT(mResults[branch].mLnl[0], mResults[branch].mLnl[1]))
		{
			mJobStatus[branch*JOBS_PER_BRANCH+JOB_BEB] = JOB_SKIP;
			continue;
		}

		// Assign the BEB job
		aJob[0] = JOB_BEB;
		aJob[1] = static_cast<int>(branch);

		// Mark it as assigned
		size_t idx = branch*JOBS_PER_BRANCH+JOB_BEB;
		mJobStatus[idx] = JOB_ASSIGNED;
		mWorkList[idx]  = aRank;
		return true;
	}

	// If no job available, shutdown this worker
	aJob[0] = JOB_SHUTDOWN;
	aJob[1] = 0;
	return false;
}


int HighLevelCoordinator::WorkTable::markJobFinished(int aRank)
{
	for(unsigned int i=0; i < mNumInternalBranches*JOBS_PER_BRANCH; ++i)
	{
		if(mJobStatus[i] == JOB_ASSIGNED && mWorkList[i] == aRank)
		{
			mJobStatus[i] = JOB_COMPLETED;
			return i;
		}
	}

	throw FastCodeMLFatal("No job found in markJobFinished");
}



void HighLevelCoordinator::WorkTable::checkAllJobsDone(void) const
{
	bool any_error = false;
	for(unsigned int i=0; i < mJobStatus.size(); ++i)
	{
		if(mJobStatus[i] != JOB_COMPLETED && mJobStatus[i] != JOB_SKIP)
		{
			std::cerr << "Job kind: " << (i%JOBS_PER_BRANCH) << " for branch " << (i/JOBS_PER_BRANCH) << " still in status: " << mJobStatus[i] << std::endl;
			any_error = true;
		}
	}
	
	if(any_error) throw FastCodeMLFatal();
}


HighLevelCoordinator::HighLevelCoordinator(int* aRgc, char*** aRgv) : mVerbose(0), mRank(MASTER_JOB), mSize(0), mNumInternalBranches(0), mWorkTable(0)
{
#ifdef _OPENMP
	const int requested = MPI_THREAD_SERIALIZED;
#else
    const int requested = MPI_THREAD_SINGLE;
#endif
	int provided;
	MPI_Init_thread(aRgc, aRgv, requested, &provided);
	if(provided != requested)
	{
		mSize = 0;		// Don't support threads, don't execute under MPI
	}
	else
	{
		MPI_Comm_size(MPI_COMM_WORLD, &mSize);
		MPI_Comm_rank(MPI_COMM_WORLD, &mRank);
	}

	// If there are too few MPI processes, terminate the unused workers
	if(mSize < 3 && mRank  != MASTER_JOB)
	{
		MPI_Finalize();
		throw FastCodeMLSuccess();
	}
}


HighLevelCoordinator::~HighLevelCoordinator()
{
	delete mWorkTable;
	MPI_Finalize();
}


bool HighLevelCoordinator::startWork(Forest& aForest, const CmdLine& aCmdLine, unsigned int aVerbose)
{
	// You need more than 2 MPI process to take advantage of it. Otherwise run as single process, OpenMP only.
	if(mSize < 3) return false;

	// Start the jobs
	if(mRank == MASTER_JOB)
	{
		// Initialize structures
		mVerbose = aVerbose;
		mNumInternalBranches = aForest.getNumInternalBranches();

		// Check if the number of worker is ok
		if(mSize > static_cast<int>(2*mNumInternalBranches+1) && mVerbose >= VERBOSE_INFO_OUTPUT) std::cerr << "Too many MPI jobs. " << mSize-1-2*mNumInternalBranches << " of them will not be used" << std::endl;

		// In the master initialize the work table
		delete mWorkTable;
		mWorkTable = new WorkTable(mNumInternalBranches);

		// In the master process initialize the master
		doMaster();
	}
	else
	{
		// Start a worker
		doWorker(aForest, aCmdLine);
	}

	// All done
	return true;
}


void HighLevelCoordinator::doMaster(void)
{
	// Push work to free workers
	unsigned int num_workers = 0;

	// Prepare variables to hold results from workers
	std::vector<double> results_double;
	std::vector<int>    results_integer;

	for(;;)
	{
		// Wait for a request of work packet
		int job_request[2];
		MPI_Status status;
		MPI_Recv((void*)&job_request, 2, MPI_INTEGER, MPI_ANY_SOURCE, MSG_WORK_REQUEST, MPI_COMM_WORLD, &status);
		int worker = status.MPI_SOURCE;

		// Act on the request (job_request[0] values are from the JobRequestType enum, [1] is the response length)
		switch(job_request[0])
		{
		case REQ_ANNOUNCE_WORKER:
			// This is an initial request for work
			++num_workers;
			break;

		case REQ_HX_RESULT:
			{
			// Get the variables and last the loglikelihood value
			results_double.resize(job_request[1]);
			MPI_Recv((void*)&results_double[0], job_request[1], MPI_DOUBLE, worker, MSG_GET_RESULTS, MPI_COMM_WORLD, &status);

			// Mark the step as done (and compute branch and hypothesis)
			int idx    = mWorkTable->markJobFinished(worker);
			int branch = idx / JOBS_PER_BRANCH;
			int h      = idx % JOBS_PER_BRANCH;

			// Save all results (lnl + all variables)
			double lnl = results_double[job_request[1]-1];
			mWorkTable->mResults[branch].mLnl[h] = lnl;
			mWorkTable->mResults[branch].mHxVariables[h].assign(results_double.begin(), results_double.end()-1);

			// Output a status message
			if(mVerbose >= VERBOSE_MORE_DEBUG) std::cerr << std::fixed << std::setprecision(8) << "Lnl: " << lnl << " for H" << h << " from worker " << worker << std::endl;
			}
			break;

		case REQ_BEB_RESULT:
			{
			// Get results
			if(job_request[1] > 0)
			{
				results_integer.resize(job_request[1]);
				MPI_Recv((void*)&results_integer[0], job_request[1], MPI_INTEGER, worker, MSG_GET_RESULTS, MPI_COMM_WORLD, &status);
			}
			
			// Mark the step as done (and compute branch)
			int idx    = mWorkTable->markJobFinished(worker);
			int branch = idx / JOBS_PER_BRANCH;

			// Save all results (positive selection sites and corresponding probability)
			mWorkTable->mResults[branch].mPositiveSelSites.clear();
			mWorkTable->mResults[branch].mPositiveSelProbs.clear();
			for(int i=0; i < job_request[1]/2; ++i)
			{
				int site = results_integer[2*i+0];
				mWorkTable->mResults[branch].mPositiveSelSites.push_back(site);
				double prob = static_cast<double>(results_integer[2*i+1])/static_cast<double>(PROB_SCALING);
				mWorkTable->mResults[branch].mPositiveSelProbs.push_back(prob);
			}

			// Output a status message
			if(mVerbose >= VERBOSE_MORE_DEBUG) std::cerr << "BEB num of results: " << job_request[1]/2 << " from worker " << worker << std::endl;
			}
			break;
		}

		// Send work packet or shutdown request (job[1] is the fg branch)
		int job[2];
		mWorkTable->getNextJob(job, worker);
		MPI_Send((void*)job, 2, MPI_INTEGER, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
		if(mVerbose >= VERBOSE_MORE_DEBUG)
		{
			switch(job[0])
			{
			case JOB_H0:
				std::cerr << "Sent H0 [" << job[1] << "] to "  << worker << std::endl;
				break;
			case JOB_H1:
				std::cerr << "Sent H1 [" << job[1] << "] to "  << worker << std::endl;
				break;
			case JOB_BEB:
				std::cerr << "Sent BEB [" << job[1] << "] to " << worker << std::endl;
				break;
			case JOB_SHUTDOWN:
				std::cerr << "Sent SHUTDOWN to "               << worker << std::endl;
				break;
			default:
				std::cerr << "Sent " << job[0] << " [" << job[1] << "] to " <<  worker << std::endl;
				break;
			}
		}

		// If no more jobs
		if(job[0] == JOB_SHUTDOWN)
		{
			--num_workers;
			if(mVerbose >= VERBOSE_MORE_DEBUG) std::cerr << "Workers remaining " << num_workers << std::endl;
			if(num_workers == 0) break;
		}
	}

	// Verify all jobs have been done
	mWorkTable->checkAllJobsDone();

	// Print likelihoods
	if(mVerbose < VERBOSE_ONLY_RESULTS) return;
	std::cerr << std::endl;
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		std::cerr << "Branch: "   << std::fixed << std::setw(3) << branch <<
					 "  Lnl H0: " << std::setw(12) << std::setprecision(15) << mWorkTable->mResults[branch].mLnl[0] << 
					 "  Lnl H1: " << std::setw(12) << std::setprecision(15) << mWorkTable->mResults[branch].mLnl[1] << std::endl;
	}

	// Check if at least one site is under positive selection
	bool has_positive_selection_sites = false;
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		if(!mWorkTable->mResults[branch].mPositiveSelSites.empty())
		{
			has_positive_selection_sites = true;
			break;
		}
	}

	// If there are print the site and corresponding probability
	if(has_positive_selection_sites)
	{
		std::cerr << std::endl << "Positive selection sites" << std::endl;
		for(size_t branch=0; branch < mNumInternalBranches; ++branch)
		{
			if(mWorkTable->mResults[branch].mPositiveSelSites.empty()) continue;

			std::cerr << "Branch: "   << std::fixed << std::setw(3) << branch << std::endl;
			for(size_t pss=0; pss < mWorkTable->mResults[branch].mPositiveSelSites.size(); ++pss)
			{
				std::cerr << std::setw(5) << mWorkTable->mResults[branch].mPositiveSelSites[pss] <<
							 std::fixed << std::setw(12) << std::setprecision(6) << mWorkTable->mResults[branch].mPositiveSelProbs[pss] << std::endl;
			}
		}
	}
}


void HighLevelCoordinator::doWorker(Forest& aForest, const CmdLine& aCmdLine)
{
	// Initialize the two hypothesis
	BranchSiteModelNullHyp h0(aForest, aCmdLine);
	BranchSiteModelAltHyp  h1(aForest, aCmdLine);

	// This value signals that this is the first work request
	int job_request[2] = {REQ_ANNOUNCE_WORKER, 0};

	// Variables for communication between master and workers
	std::vector<double> values_double;
	std::vector<int>    values_integer;
	MPI_Status          status;

	for(;;)
	{
		// Signal I'm ready for work
		MPI_Request request;
		MPI_Isend((void*)&job_request, 2, MPI_INTEGER, MASTER_JOB, MSG_WORK_REQUEST, MPI_COMM_WORLD, &request);

		// If needed (i.e. it is not a worker announcement), send step results
		switch(job_request[0])
		{
		case REQ_HX_RESULT:
			MPI_Send((void*)&values_double[0], job_request[1], MPI_DOUBLE, MASTER_JOB, MSG_GET_RESULTS, MPI_COMM_WORLD);
			break;

		case REQ_BEB_RESULT:
			if(job_request[1]) MPI_Send((void*)&values_integer[0], job_request[1], MPI_INTEGER, MASTER_JOB, MSG_GET_RESULTS, MPI_COMM_WORLD);
			break;
		}

		// Receive the job to execute or the shutdown request
		int job[2];
		MPI_Recv((void*)job, 2, MPI_INTEGER, MASTER_JOB, MSG_NEW_JOB, MPI_COMM_WORLD, &status);

		switch(job[0])
		{
		case JOB_SHUTDOWN:
			return;

		case JOB_H0:
			{
			// Compute H0
			if(aCmdLine.mInitFromParams)		h0.initFromTreeAndParams();
			else if(aCmdLine.mTimesFromFile)	h0.initFromTree();
			double lnl = h0(job[1]);

			// Assemble the results to be passed to the master
			h0.getVariables(values_double);
			values_double.push_back(lnl);
			job_request[0] = REQ_HX_RESULT;
			job_request[1] = static_cast<int>(values_double.size());
			}
			break;

		case JOB_H1:
			{
			// Compute H1
			if(aCmdLine.mInitFromParams)		h1.initFromTreeAndParams();
			else if(aCmdLine.mTimesFromFile)	h1.initFromTree();
			double lnl = h1(job[1]);

			// Assemble the results to be passed to the master
			h1.getVariables(values_double);
			values_double.push_back(lnl);
			job_request[0] = REQ_HX_RESULT;
			job_request[1] = static_cast<int>(values_double.size());
			}
			break;

		case JOB_BEB:
			{
			// Compute the BEB
			BayesTest bt(aForest.getNumSites());
			bt.computeBEB();

			// Extract the results
			std::vector<unsigned int> positive_sel_sites;
			std::vector<double>       positive_sel_sites_prob;
			bt.extractPositiveSelSites(positive_sel_sites, positive_sel_sites_prob);
			size_t num_sites = positive_sel_sites.size();

			// Assemble the results
			job_request[0] = REQ_BEB_RESULT;
			job_request[1] = 2*static_cast<int>(num_sites);
			if(num_sites)
			{
				// Assemble consecutive pairs (site, probability) into an array of integers
				values_integer.clear();
				for(size_t i=0; i < num_sites; ++i)
				{
					values_integer.push_back(positive_sel_sites[i]);
					int v = static_cast<int>(positive_sel_sites_prob[i]*PROB_SCALING+0.5);
					values_integer.push_back(v);
				}
			}
			}
			break;
		}
	}
}

#endif

