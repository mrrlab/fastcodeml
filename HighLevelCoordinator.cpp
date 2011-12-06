
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
	MSG_NEW_JOB				///< New job from the master
};

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

	unsigned int mNumInternalBranches;	///< Number of internal branches that can be marked as foreground branch.
	std::vector<int> mJobStatus;		///< Corresponding step status
	std::vector<int> mWorkList;			///< Who is doing this step
	std::vector<double> mResults;		///< LnL results from H0 and H1 (here the BEB column is here to simplify processing)

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

		// If the previous results passes the LRT, then compute the BEB
		// LRT test: -2*(lnl0-lnl1) > 3.84
		if(mResults[branch*JOBS_PER_BRANCH+JOB_H1] - mResults[branch*JOBS_PER_BRANCH+JOB_H0] <= 1.920729)
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
#ifdef USE_THREAD_MPI
	//const int requested = MPI_THREAD_MULTIPLE;
	const int requested = MPI_THREAD_SINGLE;
	//const int requested = MPI_THREAD_SERIALIZED;
#else
	const int requested = MPI_THREAD_SERIALIZED;
#endif
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


bool HighLevelCoordinator::startWork(Forest& aForest, unsigned int aSeed, unsigned int aVerbose, bool aNoMaximization, bool aTimesFromFile, unsigned int aOptimizationAlgo)
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
		if(mSize > (int)(2*mNumInternalBranches+1) && mVerbose >= 1) std::cerr << "Too many MPI jobs. " << mSize-1-2*mNumInternalBranches << " of them will not be used" << std::endl;

		// In the master initialize the work table
		mWorkTable = new WorkTable(mNumInternalBranches);

		// In the master process initialize the master or master+one worker
#if defined(_OPENMP) && defined(USE_THREAD_MPI)
		unsigned int nthreads = omp_get_max_threads();
		if(nthreads < 2)
		{
			doMaster();
		}
		else
		{
			omp_set_nested(1); 

			#pragma omp parallel
			{
				#pragma omp sections nowait
				{
					#pragma omp section
					{
						doMaster();
					}
					#pragma omp section
					{
						if(nthreads > 2) omp_set_num_threads(nthreads-1);
						doWorker(aForest, aSeed, aNoMaximization, aTimesFromFile, aOptimizationAlgo);
					}
				}
			}
		}
#else
		doMaster();
#endif
	}
	else
	{
		doWorker(aForest, aSeed, aNoMaximization, aTimesFromFile, aOptimizationAlgo);
	}

	// All done
	return true;
}

void HighLevelCoordinator::doMaster(void)
{
	// Push work to free workers
	unsigned int num_workers = 0;
	for(;;)
	{
		// Wait for a request of work packet (or the result of the previous work packet)
		double lnl;
		MPI_Status status;
		MPI_Recv((void*)&lnl, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MSG_WORK_REQUEST, MPI_COMM_WORLD, &status);
		int worker = status.MPI_SOURCE;

		// If message contains the result of a previous step, collect the results
		if(lnl < DBL_MAX)
		{
			int idx = mWorkTable->markJobFinished(worker);
			int job_type = mWorkTable->getJobType(idx);
			if(mVerbose >= 1) std::cerr << "Lnl: " << lnl << " for task: " << job_type << " from worker " << worker << std::endl;
			mWorkTable->mResults[idx] = lnl;
		}
		else
		{
			// This is an initial request for work
			++num_workers;
		}

		// Send work packet or shutdown request
		int job[2];
		mWorkTable->getNextJob(job, worker);
		MPI_Send((void*)job, 2, MPI_INTEGER, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
		if(mVerbose >= 1)
		{
			switch(job[0])
			{
			case JOB_H0:
				std::cerr << "Sent H0 [" << job[1] << "] to " <<  worker << std::endl;
				break;
			case JOB_H1:
				std::cerr << "Sent H1 [" << job[1] << "] to " <<  worker << std::endl;
				break;
			case JOB_BEB:
				std::cerr << "Sent BEB [" << job[1] << "] to " <<  worker << std::endl;
				break;
			case JOB_SHUTDOWN:
				std::cerr << "Sent SHUTDOWN to " <<  worker << std::endl;
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
			continue;
		}

		// Send the work packet content
		//world.send(k, tag_data, s.unit_cells[job]);
		//world.send(k, tag_data, s.num_atoms[job]);
		//world.send(k, tag_data, s.coords[job]);
		//world.send(k, tag_data, cutoff);
	}

	// Verify all jobs have been done
	mWorkTable->checkAllJobsDone();

	// Print results
	std::cerr << std::endl;
	for(unsigned int branch=0; branch < mNumInternalBranches; ++branch)
	{
		std::cerr << "Branch: "  << std::setw(3)  << branch <<
					"  Lnl H0: " << std::setw(12) << std::setprecision(8) << mWorkTable->mResults[branch*JOBS_PER_BRANCH+JOB_H0] << 
					"  Lnl H1: " << std::setw(12) << std::setprecision(8) << mWorkTable->mResults[branch*JOBS_PER_BRANCH+JOB_H1] << std::endl;
	}
}

void HighLevelCoordinator::doWorker(Forest& aForest, unsigned int aSeed, bool aNoMaximization, bool aTimesFromFile, unsigned int aOptimizationAlgo)
{
	// This value signals that this is the first work request
	double lnl = DBL_MAX;
	for(;;)
	{
		// Signal I'm ready for work
		MPI_Request request;
		MPI_Isend((void*)&lnl, 1, MPI_DOUBLE, MASTER_JOB, MSG_WORK_REQUEST, MPI_COMM_WORLD, &request);

		// Receive the job number or the shutdown request
		int job[2];
		MPI_Status status;
		MPI_Recv((void*)job, 2, MPI_INTEGER, MASTER_JOB, MSG_NEW_JOB, MPI_COMM_WORLD, &status);

		switch(job[0])
		{
		case JOB_SHUTDOWN:
			return;

		case JOB_H0:
			{
			BranchSiteModelNullHyp h0(aForest.getNumBranches(), aForest.getNumSites(), aSeed);
			lnl = h0.computeModel(aForest, job[1], aNoMaximization, aTimesFromFile, false, aOptimizationAlgo);
			}
			break;

		case JOB_H1:
			{
			BranchSiteModelAltHyp h1(aForest.getNumBranches(), aForest.getNumSites(), aSeed);
			const double* starting_values = 0; // for now
			lnl = h1.computeModel(aForest, job[1], aNoMaximization, aTimesFromFile, false, starting_values, aOptimizationAlgo);
			}
			break;

		case JOB_BEB:
			{
			BayesTest bt(aForest.getNumSites());
			bt.computeBEB();
			std::vector<unsigned int> aPositiveSelSites;
			std::vector<double> aPositiveSelSitesProb;
			unsigned int num_sites = bt.extractPositiveSelSites(aPositiveSelSites, aPositiveSelSitesProb);
			if(num_sites)
			{
				// Return the two vectors to master
			}
			lnl = 0;
			}
			break;
		}

		// Receive all the work data
		//std::vector<float> unit_cell;
		//world.recv(master_process, tag_data, unit_cell);
		//std::vector<int> num_atoms;
		//world.recv(master_process, tag_data, num_atoms);
		//std::vector<float> coords;
		//world.recv(master_process, tag_data, coords);
		//float cutoff;
		//world.recv(master_process, tag_data, cutoff);
	}
}


#endif

