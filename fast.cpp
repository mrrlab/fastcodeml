/// @mainpage FastCodeML
///
/// @section intro_sect Introduction
///
/// FastCodeML is a rewrite of CodeML based directly on the pseudocode document.
/// It incorporates various parallelization strategies to be able to exploit modern HPC machines architecture.
/// For this reason there are various parts of the code that can be selected at compile time or run time to experiment with various, possible solutions.
///
/// @section contacts_sect Contacts
///
/// Contact us if you want more information on the project, want to collaborate or suggest new ideas.
///
///- Ing. <a href="mailto:mvalle@cscs.ch">Mario Valle</a> - Swiss National Supercomputing Centre (CSCS) - Switzerland
///- The HP2C <a href="mailto:selectome@hp2c.ch">Selectome</a> Project Group - Mainly based in University of Lausanne - Switzerland
///

#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include "CmdLine.h"
#include "Newick.h"
#include "Phylip.h"
#include "BayesTest.h"
#include "Forest.h"
#include "ForestGroup.h"
#include "Exceptions.h"
#include "ParseParameters.h"
#include "VerbosityLevels.h"
#include "BranchSiteModel.h"
#include "WriteResults.h"

#ifndef VTRACE
#ifdef _OPENMP
#include <omp.h>
#endif
#endif
#ifdef USE_MKL_VML
#include <mkl_vml.h>
#endif
#include "Timer.h"
#ifdef USE_MPI
#include "HighLevelCoordinator.h"
#endif

/// Main program for FastCodeML.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.1
///
///	@param[in] aRgc Number of command line parameters
/// @param[in] aRgv Command line parameters
///
const char* version="1.2.1";

int main(int aRgc, char **aRgv)
{
	try
	{

#ifdef USE_MKL_VML
	// If used, intitialize the MKL VML library
	vmlSetMode(VML_HA|VML_DOUBLE_CONSISTENT);
#endif

    int num_jobs = 0;

#ifdef USE_MPI
	// Start the high level parallel executor (based on MPI)
	HighLevelCoordinator hlc(&aRgc, &aRgv);
	num_jobs = hlc.numJobs();
#endif

	// Parse the command line
	CmdLine cmd;
	cmd.parseCmdLine(aRgc, aRgv);

	// Adjust and report the number of threads that will be used
#ifdef _OPENMP
	int num_threads = omp_get_max_threads();
    //std::cout<<"max num of thr: "<< num_threads <<std::endl;

     if((cmd.mNumThreads >=1)&&(cmd.mNumThreads <= (unsigned int)num_threads))
        num_threads = cmd.mNumThreads;
    // std::cout<<"num of thr: "<< num_threads <<std::endl;

     omp_set_num_threads(num_threads);
         HighLevelCoordinator::num_threads = num_threads;

#else
	cmd.mNumThreads=1;
	int num_threads = 1;
	cmd.mForceSerial = true;
	    HighLevelCoordinator::num_threads = num_threads;

#endif


    std::ostringstream header;
    header <<std::endl<<"------------------"<< std::endl<<"FastCodeML V"<<version<<std::endl<<"------------------" << std::endl;
#ifdef USE_MPI
	// Shutdown messages from all MPI processes except the master
	if(!hlc.isMaster()) cmd.mVerboseLevel = VERBOSE_NONE;
	if(hlc.isMaster())  std::cout << header.str() <<std::endl;
#else
    std::cout << header.str() << std::endl;
#endif

	// Write out command line parameters (if not quiet i.e. if verbose level > 0)
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) cmd.printCmdLine(num_threads, num_jobs);

	// Initialize the random number generator (0 means it is not set on the command line)
#ifdef USE_MPI
	// Insure that each MPI process starts with a different seed
	if(cmd.mSeed == 0) cmd.mSeed = static_cast<unsigned int>(time(NULL)) + static_cast<unsigned int>(hlc.getRank()) * 1000;
#else
	if(cmd.mSeed == 0) cmd.mSeed = static_cast<unsigned int>(time(NULL));
#endif
	srand(cmd.mSeed);

	// Verify the optimizer algorithm selected on the command line
	if(!cmd.mNoMaximization) BranchSiteModel::verifyOptimizerAlgo(cmd.mOptimizationAlgo);

	// Start a timer (to measure serial part over parallel one)
	Timer timer;
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) timer.start();

    ForestGroup *forestGroup = new ForestGroup();
    forestGroup->initForests(cmd);

	// Get the time needed by data preprocessing
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {timer.stop(); std::cout << std::endl << "TIMER (preprocessing) ncores: " << std::setw(2) << num_threads << " time: " << timer.get() << std::endl;}

#ifdef USE_MPI
	// Distribute the work. If run under MPI then finish, else return to the standard execution flow
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) timer.start();
	bool has_run_under_MPI = hlc.startWork(forestGroup, cmd);

	// If executed under MPI report the time spent, otherwise stop the timer so it can be restarted around the serial execution
	if(has_run_under_MPI)
	{
		if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {timer.stop(); std::cout << std::endl << "TIMER (processing) ncores: " << std::setw(2) << num_threads*(hlc.numJobs()-1)+1 << " time: " << timer.get() << std::endl;}
		return 0;
	}
	else
	{
		timer.stop();
	}
#endif

	// Start timing parallel part
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) timer.start();

    std::vector<Forest*> forests(forestGroup->getForests());

    std::ostringstream results;
    for (size_t ii = 0; ii < forests.size(); ii++)
    {
        results << forestGroup->solveForest(*(forests[ii]), ii, cmd);
    }

	// Get the time needed by the parallel part
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {timer.stop(); std::cout << std::endl << "TIMER (processing) ncores: " << std::setw(2) << num_threads << " time: " << timer.get() << std::endl;}

    WriteResults::outputResultsToFile(cmd.mResultsFile, results.str());

    // Clean up
    if (forestGroup != NULL) delete(forestGroup);
    forestGroup = NULL;

	////////////////////////////////////////////////////////////////////
	// Catch all exceptions
	}
	catch(const FastCodeMLSuccess&)
	{
		return 0;
	}
	catch(const FastCodeMLFatal& e)
	{
		// If a message associated (i.e. no empty string), display it
		if(e.what()[0]) std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	}
	catch(const FastCodeMLMemoryError& e)
	{
		std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	}
	catch(const std::bad_alloc& e)
	{
		std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	}
	catch(...)
	{
		std::cout << std::endl << "Default exception caught." << std::endl << std::endl;
		return 1;
	}

	return 0;
}

/// @page cppstd_page C++ Coding Standard
/// Here are collected few rules for coding this project.
///
/// @section cnames_sect Class names
/// Class names are CamelCase with first letter uppercase.
///
/// Ex: %PhyloTree
///
/// @section cmeth_sect Class methods
/// Class methods names are CamelCase with the first letter lowercase.
/// Only very common and specific names should be all lowercase, like read, clean, size.
///
/// Ex: testFillQ
///
/// @section cdatamemb_sect Class data members
/// Class member variables names start with 'm' followed by CamelCase name.
///
/// Ex: mFgBranch
///
/// @section carg_sect Function arguments
/// Function arguments names start with 'a' followed by CamelCase name.
///
/// Ex: aFgBranch
///
/// @section const_sect Constants and enumeration
/// Constants and enumerations are all uppercase with words separated by '_'.
/// The first letters specify the kind of constant (like: STS_ for status, OPT_ for option value).
///
/// Ex: STS_CANT_OPEN
///
/// @section stack_sect Temporary variables
/// All the other variables are all lower case with parts separated by '_'.
///
/// Ex: branch_list
///
/// @section misc_sect Miscellaneous rules
/// In case of error main should return 1.
///
/// Array sizes and corresponding indexes should be size_t. The remaining counters should be unsigned int.
///
/// The null pointer should be written as NULL, not 0 to make clear its purpose.
///

/**
@page cmd_page Command Line Switches
Here is a quick list of the valid command line switches for FastCodeML.

The input `tree_file` is in %Newick format with the file containing only one tree. The `alignment_file` instead is in %Phylip format.

@verbatim

Usage:
    FastCodeML [options] tree_file alignment_file

-d  --debug  -v  --verbose (required argument)
        Verbosity level (0: none; 1: results only; 2: normal info; 3: MPI trace; 4: more debug) (default: 1)

-q  --quiet (no argument)
        No messages except results

-?  -h  --help (no argument)
        This help

-s  --seed (required argument)
        Random number generator seed (0 < seed < 1000000000)

-b  --branch (required argument)
        Do only this branch as foreground branch (count from 0)

-bs  --branch-start (required argument)
        Start computing from this branch as foreground one (count from 0) (default: first one)

-be  --branch-end (required argument)
        End computing at this branch as foreground one (count from 0) (default: last one)

-i  --ignore-freq (no argument)
        Ignore computed codon frequency and set all of them to 1/61

-e  --export (required argument)
        Export forest in GML format (if %03d or @03d is present, one is created for each fg branch)

-nr  --no-reduce (no argument)
        Do not reduce forest by merging common subtrees

-l  --lengths-from-file  --times-from-file (no argument)
        Initial branch lengths from tree file

-o  --initial-step (no argument)
        Only the initial step is performed (no maximization)

-c  --export-comp-times (required argument)
        Export the computed times from H0 if 0, H1 if 1, otherwise the one read in the phylo tree

-r  --trace (no argument)
        Trace the maximization run

-nt  --number-of-threads (required argument)
        Number of threads (1 for non parallel execution)

-bf  --branch-from-file (no argument)
        Do only the branch marked in the file as foreground branch

-hy  --only-hyp (required argument)
        Compute only H0 if 0, H1 if 1

-i1  --init-from-h1 (no argument)
        Start H0 optimization from H1 results

-m  --maximizer (required argument)
        Optimizer algorithm (0:LBFGS, 1:VAR1, 2:VAR2, 3:SLSQP, 11:BOBYQA, 22:FromCodeML, 99:MLSL_LDS) (default: 22)

-sd  --small-diff (required argument)
        Delta used in gradient computation (default: 1.49e-8)

-p  --init-param (required argument)
        Pass initialization parameter in the form: P=value (P: w0, k, p0, p1, w2)

-ic  --init-default (no argument)
        Start from default parameter values and times from tree file

-x  --extra-debug (required argument)
        Extra debug parameter (zero disables it)

-re  --relative-error (required argument)
        Relative error where to stop maximization (default: 1e-3)

-ou  --output (required argument)
        Write results formatted to this file

-cl  --clean-data (no argument)
        Remove ambiguous or missing sites from the MSA (default: no)

-ps  --no-pre-stop (no argument)
        Don't stop H0 maximization even if it cannot satisfy LRT (default: stop)

-mi  --max-iterations (required argument)
        Maximum number of iterations for the maximizer (default: 10000)

-bl  --branch-lengths-fixed (no argument)
        The length of the brances is fixed

@endverbatim
*/

/// @page vampir_page Using Vampir for profiling
/// On Linux we use VampirTrace to collect profile data and Vampir to display the results (http://www.vampir.eu/).
///
/// Before invoking CMAKE define CXX=vtCC
///
/// Define CMAKE_BUILD_TYPE as: RelWithDebInfo
///
/// Run CMAKE and configure.
///
/// Then define CMAKE_CXX_FLAGS as: -vt:mt -vt:noopari
/// If you want to trace also the OpenMP calls then change it to: -vt:mt -vt:preprocess -DVTRACE
///
/// Then proceed as usual to build the executable.
///
/// Before running the executable, define the following environment variables:
///
///     export VT_BUFFER_SIZE=512M
///     export VT_MAX_FLUSHES=0
///     export VT_SYNCH_FLUSH=yes
///     export VT_GPUTRACE=no
///     export VT_UNIFY=no
///
/// Due to a VampirTrace bug, at the end of the execution, run the vtunify executable by itself.
///
/// Now you can analyze the results by running vampir on the *.otf file generated.
///
