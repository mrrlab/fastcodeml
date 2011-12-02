
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

static const int JOBS_PER_BRANCH = 3;
static const int MASTER_JOB = 0;

enum JobStatus
{
	JOB_WAITING,
	JOB_ASSIGNED,
	JOB_COMPLETED,
	JOB_SKIP
};

enum JobType
{
	JOB_H0=0,
	JOB_H1=1,
	JOB_BEB=2,
	JOB_SHUTDOWN=3
};

enum MessageType
{
	MSG_WORK_REQUEST,
	MSG_FINISH_AND_WORK_REQUEST,
	MSG_NEW_JOB
};


struct HighLevelCoordinator::WorkTable
{
	WorkTable(unsigned int aNumBranches) :
					mNumBranches(aNumBranches),
					mJobStatus(aNumBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mWorkList(aNumBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mResults(aNumBranches*JOBS_PER_BRANCH, 0.0) {}

	unsigned int mNumBranches;
	std::vector<int> mJobStatus;
	std::vector<int> mWorkList;
	std::vector<double> mResults;

	bool getNextJob(int* aJob, int aRank);
	int markJobFinished(int aRank);
	void checkAllJobsDone(void) const;
};


bool HighLevelCoordinator::WorkTable::getNextJob(int* aJob, int aRank)
{
	// Assign all H0 jobs, then all H1 jobs
	for(unsigned int kind=JOB_H0; kind <= JOB_H1; ++kind)
	{
		for(unsigned int branch=0; branch < mNumBranches; ++branch)
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
	for(unsigned int branch=0; branch < mNumBranches; ++branch)
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
		unsigned int idx = branch*JOBS_PER_BRANCH+JOB_BEB;
		aJob[0] = JOB_BEB;
		aJob[1] = branch;
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
	for(unsigned int i=0; i < mNumBranches*JOBS_PER_BRANCH; ++i)
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


HighLevelCoordinator::HighLevelCoordinator(int* aArgc, char*** aArgv) : mVerbose(0), mRank(MASTER_JOB), mSize(0), mNumBranches(0), mWorkTable(0)
{
#ifdef USE_THREAD_MPI
	const int requested = MPI_THREAD_MULTIPLE;
#else
	const int requested = MPI_THREAD_FUNNELED;
#endif
	int provided;
	MPI_Init_thread(aArgc, aArgv, requested, &provided);
	if(provided != requested) throw FastCodeMLFatal("Requested threading level not provided by MPI library");

	MPI_Comm_size(MPI_COMM_WORLD, &mSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mRank);	
}


HighLevelCoordinator::~HighLevelCoordinator()
{
	MPI_Finalize();
	delete mWorkTable;
}


bool HighLevelCoordinator::startWork(Forest& aForest, unsigned int aSeed, unsigned int aVerbose, bool aNoMaximization, bool aTimesFromFile, unsigned int aOptimizationAlgo)
{
	// If only one MPI process, return
	if(mSize <= 1) return false;

	// Initialize structures
	mVerbose = aVerbose;
	mNumBranches = aForest.getNumInternalBranches();

	// Start the jobs
	if(mRank == MASTER_JOB)
	{
		// In the master initialize the work table
		mWorkTable = new WorkTable(mNumBranches);

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
			if(mVerbose >= 1) std::cerr << "Lnl: " << lnl << " from worker " << worker << std::endl;
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
		if(mVerbose >= 1) std::cerr << "Sent [" << job[0] << ' ' << job[1] << "] to " <<  worker << std::endl;

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
	for(unsigned int branch=0; branch < mNumBranches; ++branch)
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
			bt.printPositiveSelSites(job[1]);
			}
			lnl = 0;
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

