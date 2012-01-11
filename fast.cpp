/// @mainpage FastCodeML
///
/// @section intro_sect Introduction
/// 
/// FastCodeML is a rewrite of CodeML based directly on the pseudocode document.
/// It incorporate various parallelization strategies to be able to exploit modern HPC machines architecture.
/// For this reason there are various part of the code that can be selected at compile time or run time to experiment with various, possible solutions.
///
/// @section contacts_sect Contacts
///
/// Contact us if you want more information on the project, want to collaborate or suggest new ideas.
/// 
///     - Ing. <a href="mailto:mvalle@cscs.ch">Mario Valle</a> - Swiss National Supercomputing Centre (CSCS) - Switzerland
///     - The HP2C <a href="mailto:selectome@hp2c.ch">Selectome</a> Project Group - Mainly based in University of Lausanne - Switzerland
///

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <cfloat>
#include "PhyloTree.h"
#include "CmdLine.h"
#include "Genes.h"
#include "MatrixSize.h"
#include "BayesTest.h"
#include "Forest.h"
#include "Exceptions.h"
#include "BranchSiteModel.h"
#ifdef _OPENMP
#include <omp.h>
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
///     @version 1.0
///
///	@param[in] ac Number of command line parameters
/// @param[in] av Command line parameters
///
int main(int ac, char **av)
{
	try
	{
	// If used, intitialize the MKL VML library
#ifdef USE_MKL_VML
	vmlSetMode(VML_HA|VML_DOUBLE_CONSISTENT);
#endif

#ifdef USE_MPI
	// Start the high level parallel executor (based on MPI)
	HighLevelCoordinator hlc(&ac, &av);
#endif

	// Parse the command line
	CmdLine cmd;
	cmd.parseCmdLine(ac, av);

	// Adjust and report the number of threads that will be used
#ifdef _OPENMP
	unsigned int num_threads = omp_get_max_threads();
	if(num_threads < 2 || cmd.mForceSerial)
	{
		cmd.mForceSerial = true;
		num_threads = 1;
		omp_set_num_threads(1);
	}
#else
	cmd.mForceSerial = true;
	unsigned int num_threads = 1;
#endif

#ifdef USE_MPI
	// Shutdown messages from all MPI processes except the master
	int verbose_level = cmd.mVerboseLevel;
	if(!hlc.isMaster()) cmd.mVerboseLevel = 0;
#endif

	// Write out command line parameters (if not quiet i.e. if verbose level > 0)
	if(cmd.mVerboseLevel > 0)
	{
													std::cerr << std::endl;
													std::cerr << "Tree file:     " << cmd.mTreeFile << std::endl;
													std::cerr << "Gene file:     " << cmd.mGeneFile << std::endl;
													std::cerr << "Verbose level: " << cmd.mVerboseLevel << std::endl;
		if(cmd.mSeed)								std::cerr << "Seed:          " << cmd.mSeed << std::endl;
		if(cmd.mBranchFromFile)						std::cerr << "Branch:        From tree file" << std::endl;
		else if(cmd.mBranch != UINT_MAX)			std::cerr << "Branch:        " << cmd.mBranch << std::endl;
		if(cmd.mIgnoreFreq)							std::cerr << "Codon freq.:   Ignore" << std::endl;
		if(cmd.mDoNotReduceForest)					std::cerr << "Reduce forest: Do not reduce" << std::endl;
		else if(cmd.mNoAggressiveStep)				std::cerr << "Reduce forest: Normal" << std::endl;
		else										std::cerr << "Reduce forest: Aggressive" << std::endl;
		if(cmd.mTimesFromFile)						std::cerr << "Times:         From tree file" << std::endl;
		else if(cmd.mInitFromConst)					std::cerr << "Times:         From tree file (rest hardcoded)" << std::endl;
		else if(cmd.mInitH1fromH0)					std::cerr << "Times:         From H0" << std::endl;
		if(cmd.mNoMaximization)						std::cerr << "Maximization:  No" << std::endl;
		if(cmd.mExportComputedTimes != UINT_MAX)	std::cerr << "Graph times:   From H" << cmd.mExportComputedTimes << std::endl;
		if(cmd.mTrace)								std::cerr << "Trace:         On" << std::endl;
		if(cmd.mGraphFile)							std::cerr << "Graph file:    " << cmd.mGraphFile << std::endl;
													std::cerr << "Optimizer:     " << cmd.mOptimizationAlgo << std::endl;
		if(cmd.mDeltaValueForGradient > 0.0)		std::cerr << "Delta value:   " << cmd.mDeltaValueForGradient << std::endl;

#ifdef _OPENMP
		if(num_threads > 1)
		{
													std::cerr << "Num. threads:  " << num_threads << std::endl
		                                                      << "Num. cores:    " << omp_get_num_procs() << std::endl;
		}
		else
		{
													std::cerr << "Num. threads:  1 serial" << std::endl
		                                                      << "Num. cores:    1"  << std::endl;
		}
#else
													std::cerr << "Num. threads:  1 serial" << std::endl
		                                                      << "Num. cores:    1"  << std::endl;
#endif
#ifdef USE_MPI
		if(hlc.numJobs() > 1)						std::cerr << "Num. MPI proc: 1 (master) + " << hlc.numJobs()-1 << " (workers)" << std::endl;
		else										std::cerr << "Num. MPI proc: Insufficient, single task execution" << std::endl;
#endif
													std::cerr << "Compiled with: ";
#ifdef _OPENMP
													std::cerr << "USE_OPENMP ";
#endif
#ifdef USE_MPI
													std::cerr << "USE_MPI ";
#endif
#ifdef NEW_LIKELIHOOD
													std::cerr << "NEW_LIKELIHOOD ";
#endif
#ifdef USE_ORIGINAL_PROPORTIONS
													std::cerr << "USE_ORIGINAL_PROPORTIONS ";
#endif
#ifdef USE_LAPACK
													std::cerr << "USE_LAPACK ";
#endif
#ifdef USE_MKL_VML
													std::cerr << "USE_MKL_VML";
#endif
													std::cerr << std::endl << std::endl;
	}

	// Initialize the random number generator (0 means it is not set on the command line)
	if(cmd.mSeed == 0) cmd.mSeed = (unsigned int)time(NULL);
	srand(cmd.mSeed);

	// Start a timer (to measure serial part over parallel one)
	Timer timer;
	if(cmd.mVerboseLevel >= 1) timer.start();

	// Load the genes
	Genes g(cmd.mVerboseLevel);
	g.loadGenesFile(cmd.mGeneFile);

	// Load the phylogenetic tree
	PhyloTree t(cmd.mVerboseLevel);
	t.loadTreeFile(cmd.mTreeFile);

	// Create and load the forest
	Forest forest(cmd.mVerboseLevel);
	forest.loadTreeAndGenes(t, g, cmd.mIgnoreFreq);

#ifdef CHECK_ALGO
	// Check if forest is in shape
	forest.checkForest(true);
#endif

	// Remove the genes and the phylotree objects not used anymore
	t.clear();
	g.clear();

	// Reduce the forest merging common subtrees. Add also more reduction, then clean the no more useful data.
	if(!cmd.mDoNotReduceForest)
	{
		forest.reduceSubtrees();
#ifndef NEW_LIKELIHOOD
		if(!cmd.mNoAggressiveStep) forest.addAggressiveReduction();
#endif
		forest.cleanReductionWorkingData();		
#ifdef NEW_LIKELIHOOD
		forest.prepareNewReduction();
#endif
#ifdef CHECK_ALGO
		// Recheck the forest
		forest.checkForest();
#endif
	}
#ifdef NEW_LIKELIHOOD
	else
	{
		forest.prepareNewReductionNoReuse();
	}
#endif

	// Subdivide the trees in groups based on dependencies
	forest.prepareDependencies(cmd.mForceSerial || cmd.mDoNotReduceForest);

	// Get the time needed by data preprocessing
	if(cmd.mVerboseLevel >= 1) {timer.stop(); std::cerr << std::endl << "TIMER (preprocessing) ncores: " << std::setw(2) << num_threads << " time: " << std::setprecision(3) << timer.get() << std::endl;}

#ifdef USE_MPI
	// Distribute the work. If run under MPI then finish, else return to the standard execution flow
	if(cmd.mVerboseLevel >= 1) timer.start();
	bool sts = hlc.startWork(forest, cmd.mSeed, verbose_level, cmd.mNoMaximization, cmd.mTimesFromFile, cmd.mInitFromConst, cmd.mOptimizationAlgo, cmd.mDeltaValueForGradient);

	// If executed under MPI report the time spent, otherwise stop the timer so it can be restarted around the serial execution
	if(sts)
	{
		if(cmd.mVerboseLevel >= 1) {timer.stop(); std::cerr << std::endl << "TIMER (processing) ncores: " << std::setw(2) << num_threads*(hlc.numJobs()-1)+1 << " time: " << std::setprecision(3) << timer.get() << std::endl;}
		return 0;
	}
	else
	{
		timer.stop();
	}
#endif

	// Compute the range of branches to compute
	size_t branch_start, branch_end;
	const size_t num_branches  = forest.getNumInternalBranches();
	const size_t marked_branch = forest.getMarkedInternalBranch();
	if(cmd.mBranchFromFile && marked_branch < num_branches)
	{
		// Branch from file and valid
		branch_start = marked_branch;
		branch_end   = marked_branch+1;
	}
	else if(cmd.mBranch < num_branches)
	{
		// Branch explicitely requested on the command line and valid
		branch_start = cmd.mBranch;
		branch_end   = cmd.mBranch+1;
	}
	else
	{
		// Do all internal branches
		branch_start = 0;
		branch_end   = num_branches;
	}

	// Print few statistics
	if(cmd.mVerboseLevel >= 1) std::cerr << forest;

	// Start timing parallel part
	if(cmd.mVerboseLevel >= 1) timer.start();

	// Initialize the models
	BranchSiteModelNullHyp h0(forest, cmd.mSeed, cmd.mNoMaximization, cmd.mTrace, cmd.mOptimizationAlgo, cmd.mDeltaValueForGradient);
	BranchSiteModelAltHyp  h1(forest, cmd.mSeed, cmd.mNoMaximization, cmd.mTrace, cmd.mOptimizationAlgo, cmd.mDeltaValueForGradient);

	// For all requested internal branches
	for(size_t fg_branch=branch_start; fg_branch < branch_end; ++fg_branch)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Doing branch " << fg_branch << std::endl;

		// Compute the null model maximum loglikelihood
		double lnl0 = 0;
		if(cmd.mComputeHypothesis != 1)
		{
			if(cmd.mTimesFromFile) h0.initFromTree();
			else if(cmd.mInitFromConst) h0.initFromTreeAndFixed();
			
			lnl0 = h0(fg_branch);
		}

		// Compute the alternate model maximum loglikelihood
		double lnl1 = 0;
		if(cmd.mComputeHypothesis != 0)
		{
			if(cmd.mInitH1fromH0)
			{
				std::vector<double> starting_values;
				h0.getVariables(starting_values);
				h1.initFromResult(starting_values);
			}
			else if(cmd.mTimesFromFile) h1.initFromTree();
			else if(cmd.mInitFromConst) h1.initFromTreeAndFixed();

			lnl1 = h1(fg_branch);
		}

		if(cmd.mVerboseLevel >= 1)
		{
			std::cerr << std::endl;
			if(cmd.mComputeHypothesis != 1)
			{
				std::cerr << "LnL0: ";
				std::cerr << std::setprecision(6) << lnl0;
				std::cerr << " ";
			}
			if(cmd.mComputeHypothesis != 0)
			{
				std::cerr << "LnL1: ";
				std::cerr << std::setprecision(6) << lnl1;
			}
			std::cerr << std::endl;
		}

		// If requested set the time in the forest and export to a graph visualization tool
		if(cmd.mGraphFile)
		{
			switch(cmd.mExportComputedTimes)
			{
			case 0:
				h0.saveComputedTimes();
				break;

			case 1:
				h1.saveComputedTimes();
				break;
			}
			forest.exportForest(cmd.mGraphFile, fg_branch);
		}

		// If the run passes the LRT, then compute the BEB
		if(cmd.mComputeHypothesis > 1 && BranchSiteModel::performLRT(lnl0, lnl1))
		{
			BayesTest bt(forest.getNumSites());
			bt.computeBEB();
			bt.printPositiveSelSites(fg_branch);
		}
	}

	// Get the time needed by the parallel part
	if(cmd.mVerboseLevel >= 1) {timer.stop(); std::cerr << std::endl << "TIMER (processing) ncores: " << std::setw(2) << num_threads << " time: " << std::setprecision(3) << timer.get() << std::endl;}

	////////////////////////////////////////////////////////////////////
	// Catch all exceptions
	}
	catch(FastCodeMLSuccess&)
	{
		return 0;
	}
	catch(FastCodeMLFatalNoMsg&)
	{
		// Don't print throw messages for which the message has been already printed
		return 1;
	}
	catch(FastCodeMLFatal& e)
	{
		std::cerr << std::endl << e.what() << std::endl;
		return 1;
	}
	catch(...)
	{
		std::cerr << std::endl << "Default exception" << std::endl;
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
/// Ex: PhyloTree
///
/// @section cmemb_sect Class methods
/// Class methods names are CamelCase with the first letter lowercase.
///
/// Ex: testFillQ
///
/// @section cmemb_sect Class data members
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
///
/// @page cmd_page Command Line Switches
///
/// Here is a quick list of the valid command line switches for FastCodeML.
///
/// @verbatim
/*
    Usage:
        FastCodeML [options] tree_file gene_file

    -d  --debug  -v  --verbose (required argument)
            Verbosity level

    -q  --quiet (no argument)
            No messages except results

    -?  -h  --help (no argument)
            This help

    -s  --seed (required argument)
            Random number generator seed

    -b  --branch (required argument)
            Do only this branch as foreground branch

    -i  --ignore-freq (no argument)
            Ignore computed codon frequency and set all to 1/61

    -e  --export (required argument)
            Export forest in GML format (if %03d or @03d is present,
            one is created for each fg branch)

    -nr  --no-reduce (no argument)
            Do not reduce forest by merging common subtrees

    -l  --times-from-file (no argument)
            Initial branch lengths from tree file

    -o  --initial-step (no argument)
            Only the initial step is performed (no maximization)

    -c  --export-comp-times (required argument)
            Export the computed times from H0 if 0, H1 if 1, otherwise
            the one read in the phylo tree

    -r  --trace (no argument)
            Trace the maximization run

    -na  --no-aggressive (no argument)
            Don't apply aggressive forest reduction

    -np  --no-parallel (no argument)
            Don't use parallel execution

    -bf  --branch-from-file (no argument)
            Do only the branch marked in the file as foreground branch

    -hy  --only-hyp (required argument)
            Compute only H0 if 0, H1 if 1

    -i0  --init-from-h0 (no argument)
            Start H1 optimization from H0 results

    -m  --maximizer (required argument)
            Optimizer algorithm (0:LBFGS, 1:VAR1, 2:VAR2, 3:SLSQP, 11:BOBYQA, 22:FromCodeML)

    -ic  --init-from-const (no argument)
           Initial branch lengths from tree file and the rest from
           hardcoded constants

    -sd  --small-diff (required argument)
           Delta used in gradient computation

@endverbatim
*/
