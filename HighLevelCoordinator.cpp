
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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "HighLevelCoordinator.h"
#include "BranchSiteModel.h"
#include "Exceptions.h"
#include "BayesTest.h"

/// For each internal branch there are three jobs to be executed: H0, H1 and BEB
static const int JOBS_PER_BRANCH = 3;

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
	REQ_HX_RESULT,				///< New job from the master
	REQ_BEB_RESULT				///< New job from the master
};

// Scaling value for transforming probabilities into integers for transmission
static const double PROB_SCALING = 1e9;


/// Table of work to be done and intermediate results
///
struct HighLevelCoordinator::WorkTable
{
	/// Constructor
	///
	/// @param[in] aNumInternalBranches Number of internal branches that can be marked as foreground branch.
	///
	WorkTable(unsigned int aNumInternalBranches) :
					mNumInternalBranches(aNumInternalBranches),
					mJobStatus(aNumInternalBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mWorkList(aNumInternalBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mResults(aNumInternalBranches*JOBS_PER_BRANCH, 0.0) {}

	unsigned int		mNumInternalBranches;	///< Number of internal branches that can be marked as foreground branch.
	std::vector<int>	mJobStatus;				///< Corresponding step status
	std::vector<int>	mWorkList;				///< Who is doing this step
	std::vector<double> mResults;				///< LnL results from H0 and H1 (here the BEB column contains the number of positive selection sites)

	/// Get the next job to be executed. If no more jobs then set aJob to shutdown and return false
	///
	/// @param[out] aJob The job request: [0] is set to the kind of job (JOB_H0, JOB_H1, JOB_BEB, JOB_SHUTDOWN) and [1] to the branch number (or zero for JOB_SHUTDOWN)
	/// @param [in] aRank The current worker rank
	///
	/// @return True if a job has been assigned, false if the job is JOB_SHUTDOWN
	///
	bool getNextJob(int* aJob, int aRank);

	/// Mark the job assigned to the aRank worker as finished
	///
	/// @param[in] aRank The rank of the current worker
	///
	/// @return The identifier of the finished job
	///
	/// @exception FastCodeMLFatal If job not found
	///
	int markJobFinished(int aRank);

	/// Get the job type (JOB_H0, JOB_H1, JOB_BEB) for the given job identifier.
	///
	/// @param[in] aIdx Job identifier returned by markJobFinished()
	///
	/// @return The job type
	///
	int getJobType(int aIdx) const {return aIdx % JOBS_PER_BRANCH;}

	/// Check that all jobs have been processed
	///
	/// @exception FastCodeMLFatalNoMsg If jobs still pending
	///
	void checkAllJobsDone(void) const;
};


bool HighLevelCoordinator::WorkTable::getNextJob(int* aJob, int aRank)
{
	// Assign all H0 jobs, then all H1 jobs
	for(unsigned int kind=JOB_H0; kind <= JOB_H1; ++kind)
	{
		for(unsigned int branch=0; branch < mNumInternalBranches; ++branch)
		{
			unsigned int idx = branch*JOBS_PER_BRANCH+kind;
			if(mJobStatus[idx] == JOB_WAITING)
			{
				aJob[0] = kind;
				aJob[1] = branch;
				mJobStatus[idx] = JOB_ASSIGNED;
				mWorkList[idx]  = aRank;
				return true;
			}
		}
	}

	// Assign BEB jobs if possible
	for(unsigned int branch=0; branch < mNumInternalBranches; ++branch)
	{
		// The BEB job should be pending
		if(mJobStatus[branch*JOBS_PER_BRANCH+JOB_BEB] != JOB_WAITING) continue;

		// The corresponding H0 and H1 jobs should be completed
		if(mJobStatus[branch*JOBS_PER_BRANCH+JOB_H0] != JOB_COMPLETED || mJobStatus[branch*JOBS_PER_BRANCH+JOB_H1] != JOB_COMPLETED) continue;

		// If the previous results do not pass the LRT, then skip the BEB computation
		if(!BranchSiteModel::performLRT(mResults[branch*JOBS_PER_BRANCH+JOB_H0], mResults[branch*JOBS_PER_BRANCH+JOB_H1]))
		{
			mJobStatus[branch*JOBS_PER_BRANCH+JOB_BEB] = JOB_SKIP;
			continue;
		}

		// Assign the BEB job
		aJob[0] = JOB_BEB;
		aJob[1] = branch;

		// Mark it as assigned
		unsigned int idx = branch*JOBS_PER_BRANCH+JOB_BEB;
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
	
	if(any_error) throw FastCodeMLFatalNoMsg();
}


HighLevelCoordinator::HighLevelCoordinator(int* aRgc, char*** aRgv) : mVerbose(0), mRank(MASTER_JOB), mSize(0), mNumInternalBranches(0), mWorkTable(0)
{
	const int requested = MPI_THREAD_SERIALIZED;
	int provided;
	MPI_Init_thread(aRgc, aRgv, requested, &provided);
	if(provided != requested)
	{
		mSize = 0;
		// Don't support threads, don't execute under MPI
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


bool HighLevelCoordinator::startWork(Forest& aForest, unsigned int aSeed, unsigned int aVerbose, bool aNoMaximization,
									 bool aTimesFromFile, bool aInitFromConst, unsigned int aOptimizationAlgo, double aDeltaValueForGradient)
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
		if(mSize > static_cast<int>(2*mNumInternalBranches+1) && mVerbose >= 1) std::cerr << "Too many MPI jobs. " << mSize-1-2*mNumInternalBranches << " of them will not be used" << std::endl;

		// In the master initialize the work table
		mWorkTable = new WorkTable(mNumInternalBranches);

		// In the master process initialize the master or master+one worker
		doMaster();
	}
	else
	{
		doWorker(aForest, aSeed, aNoMaximization, aTimesFromFile, aInitFromConst, aOptimizationAlgo, aDeltaValueForGradient);
	}

	// All done
	return true;
}


void HighLevelCoordinator::doMaster(void)
{
	// Push work to free workers
	unsigned int num_workers = 0;

	std::vector<double> resultsDouble;
	std::vector<int> resultsInteger;
	double lnl;

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
			resultsDouble.resize(job_request[1]);
			MPI_Recv((void*)&resultsDouble[0], job_request[1], MPI_DOUBLE, worker, MSG_GET_RESULTS, MPI_COMM_WORLD, &status);
			lnl = resultsDouble[job_request[1]-1];

			// Mark the step as done and save the results (only lnl for now)
			int idx = mWorkTable->markJobFinished(worker);
			int job_type = mWorkTable->getJobType(idx);
			if(mVerbose >= 1) std::cerr << std::fixed << std::setprecision(8) << "Lnl: " << lnl << " for task: " << job_type << " from worker " << worker << std::endl;
			mWorkTable->mResults[idx] = lnl;
			}
			break;

		case REQ_BEB_RESULT:
			{
			// Get results
			if(job_request[1] > 0)
			{
				resultsInteger.resize(job_request[1]);
				MPI_Recv((void*)&resultsInteger[0], job_request[1], MPI_INTEGER, worker, MSG_GET_RESULTS, MPI_COMM_WORLD, &status);
			}
			
			// Look how to code and visualize the results of BEB
			int idx = mWorkTable->markJobFinished(worker);
			if(mVerbose >= 1) std::cerr << "BEB num of results: " << job_request[1]/2 << " from worker " << worker << std::endl;
			mWorkTable->mResults[idx] = static_cast<double>(job_request[1]/2);
			}
			break;
		}

		// Send work packet or shutdown request (job[1] is the fg branch)
		int job[2];
		mWorkTable->getNextJob(job, worker);
		MPI_Send((void*)job, 2, MPI_INTEGER, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
		if(mVerbose >= 1)
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
			if(mVerbose >= 1) std::cerr << "Workers remaining " << num_workers << std::endl;
			if(num_workers == 0) break;
		}
	}

	// Verify all jobs have been done
	mWorkTable->checkAllJobsDone();

	// Print results
	std::cerr << std::endl;
	for(unsigned int branch=0; branch < mNumInternalBranches; ++branch)
	{
		std::cerr << "Branch: "   << std::fixed << std::setw(3) << branch <<
					 "  Lnl H0: " << std::setw(12) << std::setprecision(8) << mWorkTable->mResults[branch*JOBS_PER_BRANCH+JOB_H0] << 
					 "  Lnl H1: " << std::setw(12) << std::setprecision(8) << mWorkTable->mResults[branch*JOBS_PER_BRANCH+JOB_H1];
		int ns = static_cast<int>(mWorkTable->mResults[branch*JOBS_PER_BRANCH+JOB_BEB]);
		if(ns > 0) 	std::cerr << "  Pos. Sel. Sites: "   << std::setw(3) << ns;
		std::cerr << std::endl;
	}
}

void HighLevelCoordinator::doWorker(Forest& aForest, unsigned int aSeed, bool aNoMaximization, bool aTimesFromFile, bool aInitFromConst,
									unsigned int aOptimizationAlgo, double aDeltaValueForGradient)
{
	BranchSiteModelNullHyp h0(aForest, aSeed, aNoMaximization, false, aOptimizationAlgo, aDeltaValueForGradient);
	BranchSiteModelAltHyp  h1(aForest, aSeed, aNoMaximization, false, aOptimizationAlgo, aDeltaValueForGradient);

	// This value signals that this is the first work request
	int job_request[2] = {REQ_ANNOUNCE_WORKER, 0};
	std::vector<double> valuesDouble;
	std::vector<int> valuesInteger;
	MPI_Status status;

	for(;;)
	{
		// Signal I'm ready for work
		MPI_Request request;
		MPI_Isend((void*)&job_request, 2, MPI_INTEGER, MASTER_JOB, MSG_WORK_REQUEST, MPI_COMM_WORLD, &request);

		// If needed, send step results
		switch(job_request[0])
		{
		case REQ_HX_RESULT:
			MPI_Send((void*)&valuesDouble[0], job_request[1], MPI_DOUBLE, MASTER_JOB, MSG_GET_RESULTS, MPI_COMM_WORLD);
			break;

		case REQ_BEB_RESULT:
			if(job_request[1]) MPI_Send((void*)&valuesInteger[0], job_request[1], MPI_INTEGER, MASTER_JOB, MSG_GET_RESULTS, MPI_COMM_WORLD);
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
			if(aTimesFromFile) h0.initFromTree();
			else if(aInitFromConst) h0.initFromTreeAndFixed();
			double lnl = h0(job[1]);

			// Assemble the results to be passed to the master
			h0.getVariables(valuesDouble);
			valuesDouble.push_back(lnl);
			job_request[0] = REQ_HX_RESULT;
			job_request[1] = valuesDouble.size();
			}
			break;

		case JOB_H1:
			{
			// Compute H1
			if(aTimesFromFile) h1.initFromTree();
			else if(aInitFromConst) h1.initFromTreeAndFixed();
			double lnl = h1(job[1]);

			// Assemble the results to be passed to the master
			h1.getVariables(valuesDouble);
			valuesDouble.push_back(lnl);
			job_request[0] = REQ_HX_RESULT;
			job_request[1] = valuesDouble.size();
			}
			break;

		case JOB_BEB:
			{
			// Compute the BEB
			BayesTest bt(aForest.getNumSites());
			bt.computeBEB();

			// Extract the results
			std::vector<unsigned int> aPositiveSelSites;
			std::vector<double> aPositiveSelSitesProb;
			bt.extractPositiveSelSites(aPositiveSelSites, aPositiveSelSitesProb);
			unsigned int num_sites = aPositiveSelSites.size();

			// Assemble the results
			job_request[0] = REQ_BEB_RESULT;
			job_request[1] = 2*num_sites;
			if(num_sites)
			{
				// Assemble consecutive pairs (site, probability) into an array of integers
				valuesInteger.clear();
				for(unsigned int i=0; i < num_sites; ++i)
				{
					valuesInteger.push_back(aPositiveSelSites[i]);
					int v = static_cast<int>(aPositiveSelSitesProb[i]*PROB_SCALING+0.5);
					valuesInteger.push_back(v);
				}
			}
			}
			break;
		}
	}
}


#endif

