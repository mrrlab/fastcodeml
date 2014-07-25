
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

int HighLevelCoordinator::mRank = -1;
int HighLevelCoordinator::mSize = -1;

int HighLevelCoordinator::num_threads = 1;

HighLevelCoordinator::HighLevelCoordinator(int* aRgc, char*** aRgv) : mVerbose(0),  mForestGroup(NULL)
{
    mSize=0;
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
    std::cout << "Hello from MASTER " << HighLevelCoordinator::getRank() << std::endl;

    BranchSiteModelAltHyp &h1 = *(mForestGroup->getAltHypothesisTest(0));

    // do branch 1
    std::cout << "Calculating H1"<< std::endl;
    double lnl = h1(0);

    std::cout << "master done: lnl is: " << lnl << std::endl;

    // if im the master send out the data!
    // param size

    std::vector<double> vars;
    h1.getVariables(vars);
    int var_size = vars.size();

    std::vector<double> data(var_size, -999.0);
    double *v = NULL;

    v = &data[0];
    //std::cout << "Master sending : " << std::setw(4) << aVar.size() << std::endl;


    for (int ii=1; ii <HighLevelCoordinator::getSize(); ii++)
    {
        std::cout << "sending kill to : " << ii << " of " << HighLevelCoordinator::getSize() << " with ";
        for (int jj = 0 ; jj < data.size(); jj++)
        {
            std::cout << data[jj];
        }
        std::cout << std::endl;
        MPI_Send(static_cast<void*>(v), data.size(), MPI_DOUBLE, ii, MSG_KICK_OFF, MPI_COMM_WORLD);
    }


} // doMaster

void HighLevelCoordinator::doWorker(const CmdLine& aCmdLine)
{

    // do branch 1
    //std::cout << "Calculating H1"<< std::endl;
    //double lnl = h1(0);


    //computeExpoential(pmatrix);

    // receive the values
    BranchSiteModelAltHyp &h1 = *(mForestGroup->getAltHypothesisTest(0));

    h1.init(0);

    std::vector<double> vars;
    h1.getVariables(vars);
    int var_size = vars.size();

    while(1)
    {
        MPI_Status          status;
        std::vector<double> values(var_size, -999.0);


        MPI_Recv(static_cast<void*>(&values[0]), var_size, MPI_DOUBLE, MASTER_JOB, MSG_KICK_OFF, MPI_COMM_WORLD, &status);


        if(values[0] < -998.0)
        {
            if (mRank == 2)
            {
            std::cout << "Worker no " << mRank << " breaking. " << values[0] <<std::endl;

            }
            break;
        }


        if(mRank == 2)
        {
           //std::cout << "kicking off do work " << mRank << " from source " << status.MPI_SOURCE << " for value ";
           //for (int ii = 0; ii < values.size(); ii++)
           //{
           //    std::cout << values[ii] << ";";
          // }
           std::cout << std::endl;
        }
        double lnl = h1.computeLikelihood(values, false);

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
