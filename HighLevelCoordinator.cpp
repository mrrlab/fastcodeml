
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
#include <limits>
#include <algorithm>
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

HighLevelCoordinator::HighLevelCoordinator(int* aRgc, char*** aRgv) : mVerbose(0), mRank(INVALID_RANK), mSize(0), mWorkTable(NULL), mForestGroup(NULL)
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
		throw FastCodeMLFatal("MPI Failed to initialize");
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
    if (mWorkTable != NULL)
    {
        delete(mWorkTable);
        mWorkTable = NULL;
    }

    if (mForestGroup != NULL)
    {
        delete(mForestGroup);
        mForestGroup = NULL;
    }

	MPI_Finalize();
}

bool HighLevelCoordinator::startWork(
     ForestGroup *aForestGroup,
     const CmdLine& aCmdLine)
{
	// You need more than 2 MPI process to take advantage of it. Otherwise run as single process, OpenMP only.
	if(mSize < 3) return false;

    mForestGroup = aForestGroup;

    // Set verbosity level.
	mVerbose = aCmdLine.mVerboseLevel;

	// Start the jobs
	if(mRank == MASTER_JOB)
	{
	    // Initialize the work table so the master can track jobs done
	    mWorkTable = new WorkTable(*aForestGroup, aCmdLine);

        // Output a debug statement based on number of processors in MPI, if requested
        checkProcCount(aCmdLine, mForestGroup->getTotalNumInternalBranches(aCmdLine));

		// In the master process initialize the master
		doMaster(aCmdLine);
	}
	else
	{
		// Start a worker
		doWorker(aCmdLine);
	}

	// All done
	return true;
}

void HighLevelCoordinator::doMaster(const CmdLine& aCmdLine)
{
	// Push work to free workers
	unsigned int num_workers = 0;

	// Prepare variables to hold results from workers
	std::vector<double> results_double;
	std::vector<int>    results_integer;

	for(;;)
	{
		// Wait for a request of work packet
		int job_request[5];
		MPI_Status status;
		MPI_Recv(static_cast<void*>(job_request), 5, MPI_INTEGER, MPI_ANY_SOURCE, MSG_WORK_REQUEST, MPI_COMM_WORLD, &status);
		int worker = status.MPI_SOURCE;

        // The current job (see below for definition of entries.)
        Job cur_job(job_request[2], job_request[3], (JobType)job_request[4], 0);

		// Act on the request (job_request[0] values are from the JobRequestType enum, [1] is the response length,
        // [2] is the forest index, [3] is the branch index, [4] is the job type)
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

			// Record the results
			mWorkTable->recordResultHyp(cur_job, results_double, worker, mVerbose);

			// Mark the step as done (and compute branch and hypothesis)
			mWorkTable->queueNextJob(cur_job, mVerbose, aCmdLine.mInitH0fromH1);
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

            // Record the results
			mWorkTable->recordResultBEB(cur_job, results_integer, worker, mVerbose);

			// Mark the step as done (and compute branch)
            mWorkTable->queueNextJob(cur_job, mVerbose, aCmdLine.mInitH0fromH1);
			}
			break;

		default:
			throw FastCodeMLFatal("Invalid job request in doMaster");
		}

		// Send work packet or shutdown request
		// (job[0] is the step to be done, job[1] is the fg branch, job[2] the length of the additional data, job[3] is the forest index in forest group)
		int job[4];
		mWorkTable->getNextJob(job);
		MPI_Send(static_cast<void*>(job), 4, MPI_INTEGER, worker, MSG_NEW_JOB, MPI_COMM_WORLD);

		// For BEB send the variables from H1; for H0 send the lnl value from corresponding H1
		if(job[2] > 0)
		{
		    std::vector<ResultSet> &result_set = *(mWorkTable->getResultSetPtr(cur_job.getForestIndex()));
			double* v = NULL;
			switch(job[0])
			{
			case JOB_BEB:
				v = &(result_set[job[1]].mHxVariables[1][0]);
				MPI_Send(static_cast<void*>(v), job[2], MPI_DOUBLE, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
				break;

			case JOB_H0:
				if(job[2] == 1)
				{
					v = &(result_set[job[1]].mLnl[1]);
				}
				else
				{
					results_double.assign(result_set[job[1]].mHxVariables[1].begin(), result_set[job[1]].mHxVariables[1].end());
					results_double.push_back(result_set[job[1]].mLnl[1]);
					v = &results_double[0];
				}
				MPI_Send(static_cast<void*>(v), job[2], MPI_DOUBLE, worker, MSG_NEW_JOB, MPI_COMM_WORLD);
				break;
			}
		}

		// Trace the messages
		if(mVerbose >= VERBOSE_MPI_TRACE)
		{
            printMPITrace(job[0], job[1], worker, job[3], std::cout);
		}

		// If no more jobs
		if(job[0] == JOB_SHUTDOWN)
		{
			--num_workers;
			if(mVerbose >= VERBOSE_MPI_TRACE) std::cout << "Workers remaining: " << num_workers << std::endl;
			if(num_workers == 0) break;
		}
	}

    std::ostringstream results;
    mWorkTable->printResults(results, aCmdLine);
    if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cout << results.str() << std::endl;

    WriteResults::outputResultsToFile(aCmdLine.mResultsFile, results.str());
} // doMaster

void HighLevelCoordinator::doWorker(const CmdLine& aCmdLine)
{
    // The current forest index
    int forest_index = INVALID_FOREST;
    int branch       = INVALID_BRANCH;
    int job_type     = (int) INVALID_JOB_TYPE;

	// This value signals that this is the first work request
	// job_request[0] = type of request, job_request[1] = the size of the results variable
	// job_request[2] = the forest this relates to, job_request[3] = the branch in the forest this relates to
    // job_request[4] - the job type (h0, h1, beb).
	int job_request[5] = {REQ_ANNOUNCE_WORKER, 0, forest_index, branch, job_type};

	// Variables for communication between master and workers
	std::vector<double> values_double;
	std::vector<int>    values_integer;
	MPI_Status          status;

	for(;;)
	{
		// Signal that I'm ready for work
		MPI_Request request;
		MPI_Isend(static_cast<void*>(job_request), 5, MPI_INTEGER, MASTER_JOB, MSG_WORK_REQUEST, MPI_COMM_WORLD, &request);

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

		// Receive the job to execute or the shutdown request
		// (job[0] the request; [1] fg branch; [2] optional number of variables, [3] the forest index)
		int job[4];
		MPI_Recv(static_cast<void*>(job), 4, MPI_INTEGER, MASTER_JOB, MSG_NEW_JOB, MPI_COMM_WORLD, &status);

		// If there is additional data
		if(job[2] > 0)
		{
			values_double.resize(static_cast<size_t>(job[2]));
			MPI_Recv(static_cast<void*>(&values_double[0]), job[2], MPI_DOUBLE, MASTER_JOB, MSG_NEW_JOB, MPI_COMM_WORLD, &status);
		}

        // Translate the job to the respective items
        forest_index = job[3];
        branch       = job[1];
        job_type     = job[0];

        // Translate the items to be appended in the job request. Master uses these to log the results
        job_request[2] = forest_index;
        job_request[3] = branch;
        job_request[4] = job_type;

		// Do the work
		switch(job[0])
		{
		case JOB_SHUTDOWN:
			return;

		case JOB_H0:
			{
            BranchSiteModelNullHyp &h0 = *(mForestGroup->getNullHypothesisTest(forest_index));

			// Initialize maximizer
			if(aCmdLine.mInitH0fromH1 && job[2] > 1) h0.initFromResult(values_double, static_cast<unsigned int>(values_double.size())-1u);
			else
			{
				if(aCmdLine.mInitFromParams)		h0.initFromParams();
				if(aCmdLine.mBranchLengthsFromFile)	h0.initFromTree();
			}

			// Get the lnl from the corresponding H1 step if any
			double threshold = (job[2] > 0) ? values_double.back()-THRESHOLD_FOR_LRT : 0.;

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
            BranchSiteModelAltHyp &h1 = *(mForestGroup->getAltHypothesisTest(forest_index));

			// Initialize maximizer
			if(aCmdLine.mInitFromParams)		h1.initFromParams();
			if(aCmdLine.mBranchLengthsFromFile)	h1.initFromTree();

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
            BayesTest &beb = *(mForestGroup->getBayesTest(forest_index));

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

void HighLevelCoordinator::checkProcCount(const CmdLine &aCmdLine, size_t aTotalNumInternalBranches) const
{
    if (mVerbose < VERBOSE_INFO_OUTPUT)
        return; // done here.

    int nb = static_cast<int>(aTotalNumInternalBranches);
    int jobs = (aCmdLine.mComputeHypothesis < 2) ? nb+1 : nb*2+1;

    if (!aCmdLine.mMultipleMode)
    {
        // Compute how many MPI processes needed (master + 1 or 2 proc per branch)
        int surplus = mSize - jobs;

        // Show if there are too many or too few processes
        if(surplus > 0)
            std::cout << "Too many MPI jobs: " << surplus << " of them will not be used." << std::endl;
        else if(surplus < 0)
            std::cout << "For top performances " << -surplus << " more MPI jobs needed." << std::endl;
    }
    else
    {
        double leverage = (double)mSize / jobs;

        // Show if the leverage < 1 (1 proc : many jobs ) or > 1 (1: proc 0-1 job)
        if (leverage > 1.)
            std::cout << "Leverage (num processors / num jobs) " << leverage << " is greater than 1 - we have 1 MPI processor per job." << std::endl;
        else
            std::cout << "Leverage  (num processors / num jobs)  " << leverage << " is less than 1. For best performance increase number of jobs by "
                << jobs - mSize << "." << std::endl;
    }
} // registerForest

void HighLevelCoordinator::printMPITrace(int aJobType, int aBranch, int aWorker, int aForestIndex,
       std::ostream &aOut) const
{
    switch(aJobType)
    {
        case JOB_H0:
            aOut << "Sent H0 [branch " << aBranch << "] to "  << aWorker << " for forest " << aForestIndex <<  std::endl;
            break;
        case JOB_H1:
            aOut << "Sent H1 [branch " << aBranch << "] to "  << aWorker  << " for forest " << aForestIndex <<  std::endl;
            break;
        case JOB_BEB:
            aOut << "Sent BEB [branch " << aBranch << "] to " << aWorker  << " for forest " << aForestIndex <<  std::endl;
            break;
        case JOB_SHUTDOWN:
            aOut << "Sent SHUTDOWN to "                      << aWorker  << std::endl;
            break;
        default:
            aOut << "Sent " << aJobType << " [branch " << aBranch << "] to "  << " for forest " << aForestIndex <<  std::endl;
            break;
    }
}

#endif
