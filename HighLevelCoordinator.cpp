
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

static const int INVALID_RANK = -1;

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

		/// Default constructor.
		///
		ResultSet()
		{
			mLnl[0] = mLnl[1] = -DBL_MAX;
		}
	};

	size_t					mNumInternalBranches;	///< Number of internal branches that can be marked as foreground branch.
	std::vector<int>		mJobStatus;				///< Corresponding step status
	std::vector<int>		mWorkList;				///< Who is doing this step
	std::vector<ResultSet>	mResults;				///< Results for each branchh

	/// Constructor
	///
	/// @param[in] aNumInternalBranches Number of internal branches that can be marked as foreground branch.
	///
	explicit WorkTable(size_t aNumInternalBranches) :
					mNumInternalBranches(aNumInternalBranches),
					mJobStatus(aNumInternalBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mWorkList(aNumInternalBranches*JOBS_PER_BRANCH, JOB_WAITING),
					mResults(aNumInternalBranches) {}

	/// Get the next job to be executed. If no more jobs then set aJob to shutdown and return false
	///
	/// @param[out] aJob The job request: [0] is set to the kind of job (JOB_H0, JOB_H1, JOB_BEB, JOB_SHUTDOWN);
	///									  [1] to the fg branch number (or zero for JOB_SHUTDOWN);
	///                                   [2] zero or the number of variables for a JOB_BEB request
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
	/// @exception FastCodeMLFatal If found jobs still pending
	///
	void checkAllJobsDone(void) const;
};


bool HighLevelCoordinator::WorkTable::getNextJob(int* aJob, int aRank)
{
	// Assign all H1 jobs
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		size_t idx = branch*JOBS_PER_BRANCH+JOB_H1;
		if(mJobStatus[idx] == JOB_WAITING)
		{
			aJob[0] = JOB_H1;
			aJob[1] = static_cast<int>(branch);
			aJob[2] = 0; // No additional data sent
			mJobStatus[idx] = JOB_ASSIGNED;
			mWorkList[idx]  = aRank;
			return true;
		}
	}

	// Then assign all H0 jobs with corresponding H1 already completed
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		size_t idx = branch*JOBS_PER_BRANCH+JOB_H0;
		if(mJobStatus[idx] == JOB_WAITING && mJobStatus[branch*JOBS_PER_BRANCH+JOB_H1] == JOB_COMPLETED)
		{
			aJob[0] = JOB_H0;
			aJob[1] = static_cast<int>(branch);
			aJob[2] = 1; // Send the lnl of the corresponding H1 step
			mJobStatus[idx] = JOB_ASSIGNED;
			mWorkList[idx]  = aRank;
			return true;
		}
	}

	// Then assign all remaining H0 jobs
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		size_t idx = branch*JOBS_PER_BRANCH+JOB_H0;
		if(mJobStatus[idx] == JOB_WAITING)
		{
			aJob[0] = JOB_H0;
			aJob[1] = static_cast<int>(branch);
			aJob[2] = 0; // No additional data sent
			mJobStatus[idx] = JOB_ASSIGNED;
			mWorkList[idx]  = aRank;
			return true;
		}
	}

	// Assign BEB jobs if possible
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		// The BEB job should be pending
		if(mJobStatus[branch*JOBS_PER_BRANCH+JOB_BEB] != JOB_WAITING) continue;

		// The corresponding H0 and H1 jobs should be completed
		if(mJobStatus[branch*JOBS_PER_BRANCH+JOB_H0] != JOB_COMPLETED || mJobStatus[branch*JOBS_PER_BRANCH+JOB_H1] != JOB_COMPLETED) continue;

		// If the previous results do not pass the LRT or H0 has been interrupted, then skip the BEB computation
		if(mResults[branch].mLnl[0] == DBL_MAX || !BranchSiteModel::performLRT(mResults[branch].mLnl[0], mResults[branch].mLnl[1]))
		{
			mJobStatus[branch*JOBS_PER_BRANCH+JOB_BEB] = JOB_SKIP;
			continue;
		}

		// Assign the BEB job
		aJob[0] = JOB_BEB;
		aJob[1] = static_cast<int>(branch);
		aJob[2] = static_cast<int>(mResults[branch].mHxVariables[1].size());

		// Mark it as assigned
		size_t idx = branch*JOBS_PER_BRANCH+JOB_BEB;
		mJobStatus[idx] = JOB_ASSIGNED;
		mWorkList[idx]  = aRank;
		return true;
	}

	// If no job available, shutdown this worker
	aJob[0] = JOB_SHUTDOWN;
	aJob[1] = 0;
	aJob[2] = 0;
	return false;
}


int HighLevelCoordinator::WorkTable::markJobFinished(int aRank)
{
	for(size_t i=0; i < mNumInternalBranches*JOBS_PER_BRANCH; ++i)
	{
		if(mJobStatus[i] == JOB_ASSIGNED && mWorkList[i] == aRank)
		{
			mJobStatus[i] = JOB_COMPLETED;
			return static_cast<int>(i);
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


HighLevelCoordinator::HighLevelCoordinator(int* aRgc, char*** aRgv) : mVerbose(0), mRank(INVALID_RANK), mSize(0), mNumInternalBranches(0), mWorkTable(NULL)
{
#ifdef _OPENMP
#ifdef VTRACE
    const int requested = MPI_THREAD_SINGLE;
#else
	const int requested = (omp_get_max_threads() <= 1) ? MPI_THREAD_SINGLE : MPI_THREAD_FUNNELED; // Change to MPI_THREAD_SERIALIZED if master process do more
#endif
#else
    const int requested = MPI_THREAD_SINGLE;
#endif
	int provided = MPI_THREAD_SINGLE;
	int mpi_status = MPI_Init_thread(aRgc, aRgv, requested, &provided);
	if(mpi_status != MPI_SUCCESS)
	{
		throw FastCodeMLFatal("MPI Failed to initalize");
	}
	else if(requested > MPI_THREAD_SINGLE && provided < MPI_THREAD_FUNNELED)
	{
		// Don't support threads. Disable them
#ifdef _OPENMP
		omp_set_num_threads(1);
#endif
	}

	// Get num of MPI processes and own rank
	MPI_Comm_size(MPI_COMM_WORLD, &mSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mRank);

	// If there are too few MPI processes, terminate the unusable workers
	if(mSize < 3 && mRank != MASTER_JOB)
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

bool HighLevelCoordinator::startWork(Forest& aForest, const CmdLine& aCmdLine)
{
	// You need more than 2 MPI process to take advantage of it. Otherwise run as single process, OpenMP only.
	if(mSize < 3) return false;

	// Start the jobs
	if(mRank == MASTER_JOB)
	{
		// If users set the fg branch tell them it is ignored
		if((aCmdLine.mBranchFromFile || aCmdLine.mBranch != UINT_MAX) && mVerbose >= VERBOSE_INFO_OUTPUT)
			std::cerr << "Cannot specify fg branch if run under MPI. Ignoring." << std::endl;

		// If the user asks for one hypothesis only
		if(aCmdLine.mComputeHypothesis < 2 && mVerbose >= VERBOSE_INFO_OUTPUT)
			std::cerr << "Cannot compute only one hypothesis under MPI. Ignoring." << std::endl;

		// Initialize structures
		mVerbose = aCmdLine.mVerboseLevel;
		mNumInternalBranches = aForest.getNumInternalBranches();

		// Check if the number of worker is ok
		if(mSize > static_cast<int>(2*mNumInternalBranches+1) && mVerbose >= VERBOSE_INFO_OUTPUT)
			std::cerr << "Too many MPI jobs: " << mSize-1-2*mNumInternalBranches << " of them will not be used." << std::endl;

		// In the master initialize the work table
		delete mWorkTable;
		mWorkTable = new WorkTable(mNumInternalBranches);

		// Prepare the results file output
		WriteResults output_results(aCmdLine.mResultsFile);

		// In the master process initialize the master
		doMaster(output_results);
	}
	else
	{
		// Start a worker
		doWorker(aForest, aCmdLine);
	}

	// All done
	return true;
}


void HighLevelCoordinator::doMaster(WriteResults& aOutputResults)
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
		MPI_Recv(static_cast<void*>(job_request), 2, MPI_INTEGER, MPI_ANY_SOURCE, MSG_WORK_REQUEST, MPI_COMM_WORLD, &status);
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
			results_double.resize(static_cast<size_t>(job_request[1]));
			MPI_Recv(static_cast<void*>(&results_double[0]), job_request[1], MPI_DOUBLE, worker, MSG_GET_RESULTS, MPI_COMM_WORLD, &status);

			// Mark the step as done (and compute branch and hypothesis)
			int idx    = mWorkTable->markJobFinished(worker);
			int branch = idx / JOBS_PER_BRANCH;
			int h      = idx % JOBS_PER_BRANCH;

			// Save all results (lnl + all variables)
			// For H1 there are also the two scale values at the end of variables for BEB computation
			double lnl = results_double[job_request[1]-1];
			mWorkTable->mResults[branch].mLnl[h] = lnl;
			mWorkTable->mResults[branch].mHxVariables[h].assign(results_double.begin(), results_double.end()-1);

			// Save for the results file (if has been computed)
			if(lnl < DBL_MAX) aOutputResults.saveLnL(static_cast<size_t>(branch), lnl, h);

			// Output a status message
			if(mVerbose >= VERBOSE_MPI_TRACE) std::cerr << std::fixed << std::setprecision(8) << "Lnl: " << lnl << " for H" << h << " from worker " << worker << std::endl;
			}
			break;

		case REQ_BEB_RESULT:
			{
			// Get results
			if(job_request[1] > 0)
			{
				results_integer.resize(static_cast<size_t>(job_request[1]));
				MPI_Recv(static_cast<void*>(&results_integer[0]), job_request[1], MPI_INTEGER, worker, MSG_GET_RESULTS, MPI_COMM_WORLD, &status);
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

			// If there are sites under positive selection, save them for the results file
			if(job_request[1] > 0)
			{
				std::vector<unsigned int> sites;
				for(int i=0; i < job_request[1]/2; ++i)
				{
					unsigned int u = static_cast<unsigned int>(results_integer[2*i+0]);
					sites.push_back(u);
				}
				aOutputResults.savePositiveSelSites(static_cast<size_t>(branch), sites, mWorkTable->mResults[branch].mPositiveSelProbs);
			}

			// Output a status message
			if(mVerbose >= VERBOSE_MPI_TRACE) std::cerr << "BEB num of results: " << job_request[1]/2 << " from worker " << worker << std::endl;
			}
			break;

		default:
			throw FastCodeMLFatal("Invalid job request in doMaster");
		}

		// Send work packet or shutdown request (job[1] is the fg branch, job[2] the length of the additional data)
		int job[3];
		mWorkTable->getNextJob(job, worker);
		MPI_Send(static_cast<void*>(job), 3, MPI_INTEGER, worker, MSG_NEW_JOB, MPI_COMM_WORLD);

		// For BEB send the variables from H1; for H0 send the lnl value from corresponding H1
		if(job[2] > 0)
		{
			double* v;
			switch(job[0])
			{
			case JOB_BEB:
				v = &(mWorkTable->mResults[job[1]].mHxVariables[1][0]);
				MPI_Send(static_cast<void*>(v), job[2], MPI_DOUBLE, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
				break;

			case JOB_H0:
				v = &(mWorkTable->mResults[job[1]].mLnl[1]);
				MPI_Send(static_cast<void*>(v), 1, MPI_DOUBLE, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
				break;
			}
		}

		// Trace the messages
		if(mVerbose >= VERBOSE_MPI_TRACE)
		{
			switch(job[0])
			{
			case JOB_H0:
				std::cerr << "Sent H0 [branch " << job[1] << "] to "  << worker << std::endl;
				break;
			case JOB_H1:
				std::cerr << "Sent H1 [branch " << job[1] << "] to "  << worker << std::endl;
				break;
			case JOB_BEB:
				std::cerr << "Sent BEB [branch " << job[1] << "] to " << worker << std::endl;
				break;
			case JOB_SHUTDOWN:
				std::cerr << "Sent SHUTDOWN to "                      << worker << std::endl;
				break;
			default:
				std::cerr << "Sent " << job[0] << " [branch " << job[1] << "] to " <<  worker << std::endl;
				break;
			}
		}

		// If no more jobs
		if(job[0] == JOB_SHUTDOWN)
		{
			--num_workers;
			if(mVerbose >= VERBOSE_MPI_TRACE) std::cerr << "Workers remaining " << num_workers << std::endl;
			if(num_workers == 0) break;
		}
	}

	// Verify all jobs have been done
	mWorkTable->checkAllJobsDone();

	// Save results in the results file for later processing
	aOutputResults.outputResults();

	// Print likelihoods
	if(mVerbose < VERBOSE_ONLY_RESULTS) return;
	std::cerr << std::endl;
	for(size_t branch=0; branch < mNumInternalBranches; ++branch)
	{
		std::cerr << "Branch: "   << std::fixed << std::setw(3) << branch;
		if(mWorkTable->mResults[branch].mLnl[0] == DBL_MAX)
		{
			std::cerr << "  Lnl H0: " << std::setw(12) << "NA";
		}
		else
		{
			std::cerr << "  Lnl H0: " << std::setw(12) << std::setprecision(15) << mWorkTable->mResults[branch].mLnl[0];
		}
		std::cerr << "  Lnl H1: " << std::setw(12) << std::setprecision(15) << mWorkTable->mResults[branch].mLnl[1] << std::endl;
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
			WorkTable::ResultSet& branch_results = mWorkTable->mResults[branch];
			if(branch_results.mPositiveSelSites.empty()) continue;

			std::cerr << "Branch: "   << std::fixed << std::setw(3) << branch << std::endl;
			for(size_t pss=0; pss < branch_results.mPositiveSelSites.size(); ++pss)
			{
				// Get probability
				double prob = branch_results.mPositiveSelProbs[pss];

				// Set significance
				const char* sig;
				if(prob > TWO_STARS_PROB)     sig = "**";
				else if(prob > ONE_STAR_PROB) sig = "*";
				else                          sig = "";

				// Adjust the site number because it starts from 1 and not zero
				std::cerr << std::setw(6) << branch_results.mPositiveSelSites[pss] + 1 <<
							 std::fixed << std::setprecision(6) << prob << sig << std::endl;
			}
		}
	}
}


void HighLevelCoordinator::doWorker(Forest& aForest, const CmdLine& aCmdLine)
{
	// Initialize the two hypothesis
	BranchSiteModelNullHyp h0(aForest, aCmdLine);
	BranchSiteModelAltHyp  h1(aForest, aCmdLine);
		
	// Initialize the BEB (no verbose at all)
	BayesTest beb(aForest, 0, aCmdLine.mDoNotReduceForest);

	// This value signals that this is the first work request
	int job_request[2] = {REQ_ANNOUNCE_WORKER, 0};

	// Variables for communication between master and workers
	std::vector<double> values_double;
	std::vector<int>    values_integer;
	MPI_Status          status;

	for(;;)
	{
		// Signal that I'm ready for work
		MPI_Request request;
		MPI_Isend(static_cast<void*>(job_request), 2, MPI_INTEGER, MASTER_JOB, MSG_WORK_REQUEST, MPI_COMM_WORLD, &request);

		// If needed (i.e. it is not a worker announcement), send step results
		switch(job_request[0])
		{
		case REQ_HX_RESULT:
			MPI_Send(static_cast<void*>(&values_double[0]), job_request[1], MPI_DOUBLE, MASTER_JOB, MSG_GET_RESULTS, MPI_COMM_WORLD);
			break;

		case REQ_BEB_RESULT:
			if(job_request[1]) MPI_Send(static_cast<void*>(&values_integer[0]), job_request[1], MPI_INTEGER, MASTER_JOB, MSG_GET_RESULTS, MPI_COMM_WORLD);
			break;
		}

		// Receive the job to execute or the shutdown request (job[0] the request; [1] fg branch; [2] optional number of variables)
		int job[3];
		MPI_Recv(static_cast<void*>(job), 3, MPI_INTEGER, MASTER_JOB, MSG_NEW_JOB, MPI_COMM_WORLD, &status);

		// If there is additional data
		if(job[2] > 0)
		{
			values_double.resize(static_cast<size_t>(job[2]));
			MPI_Recv(static_cast<void*>(&values_double[0]), job[2], MPI_DOUBLE, MASTER_JOB, MSG_NEW_JOB, MPI_COMM_WORLD, &status);
		}

		// Do the work
		switch(job[0])
		{
		case JOB_SHUTDOWN:
			return;

		case JOB_H0:
			{
			// Initialize maximizer
			if(aCmdLine.mInitFromParams)		h0.initFromTreeAndParams();
			else if(aCmdLine.mTimesFromFile)	h0.initFromTree();

			// Get the lnl from the corresponding H1 step if any
			double threshold = (job[2] > 0) ? values_double[0]-THRESHOLD_FOR_LRT : 0.;

			// Compute H0
			double lnl = h0(static_cast<size_t>(job[1]), aCmdLine.mStopIfNotLRT && job[2] > 0, threshold);
	
			// Assemble the results to be passed to the master
			h0.getVariables(values_double);
			values_double.push_back(lnl);
			job_request[0] = REQ_HX_RESULT;
			job_request[1] = static_cast<int>(values_double.size());
			}
			break;

		case JOB_H1:
			{
			// Initialize maximizer
			if(aCmdLine.mInitFromParams)		h1.initFromTreeAndParams();
			else if(aCmdLine.mTimesFromFile)	h1.initFromTree();

			// Compute H1
			double lnl = h1(static_cast<size_t>(job[1]));

			// Assemble the results to be passed to the master (variables, lnl and scales for BEB)
			h1.getVariables(values_double);
			std::vector<double> scales(2);
			h1.getScales(scales);
			values_double.push_back(scales[0]); // bg scale
			values_double.push_back(scales[1]); // fg scale
			values_double.push_back(lnl);
			job_request[0] = REQ_HX_RESULT;
			job_request[1] = static_cast<int>(values_double.size());
			}
			break;

		case JOB_BEB:
			{
			// Get the scale values
			std::vector<double> scales(2);
			scales.assign(values_double.end()-2, values_double.end());

			// Compute the BEB with the vars are taken from the master
			beb.computeBEB(values_double, static_cast<size_t>(job[1]), scales);

			// Extract the results
			std::vector<unsigned int> positive_sel_sites;
			std::vector<double>       positive_sel_sites_prob;
			beb.extractPositiveSelSites(positive_sel_sites, positive_sel_sites_prob);
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

