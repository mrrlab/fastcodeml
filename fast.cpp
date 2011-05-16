/// @mainpage FastCodeML
///
/// @section intro_sect Introduction
/// 
/// FastCodeML is a rewrite of CodeML based directly on the pseudocode document and will be the base
/// on which the parallel (HPC) version will be created.
/// Still missing computeBEB and computeNEB routines, but the rest should be functional.
///
/// @section contacts_sect Contacts
///
/// Contact us if you want more information on the project, want to collaborate or suggest new ideas.
/// 
///     - Ing. <a href="mailto:mvalle@cscs.ch">Mario Valle</a> - Swiss National Supercomputing Centre (CSCS) - Switzerland
///     - The HP2C <a href="mailto:selectome@hp2c.ch">Selectome</a> Project Group
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
#include "Timer.h"

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

	// Parse the command line
	CmdLine cmd;
	cmd.parseCmdLine(ac, av);

	// Write out command line parameters (if not quiet i.e. if verbose level > 0)
	if(cmd.mVerboseLevel > 0)
	{
		std::cerr << std::endl;
													std::cerr << "Tree file:     " << cmd.mTreeFile << std::endl;
													std::cerr << "Gene file:     " << cmd.mGeneFile << std::endl;
													std::cerr << "Verbose level: " << cmd.mVerboseLevel << std::endl;
		if(cmd.mSeed)								std::cerr << "Seed:          " << cmd.mSeed << std::endl;
		if(cmd.mBranch != UINT_MAX)					std::cerr << "Branch:        " << cmd.mBranch << std::endl;
		if(cmd.mIgnoreFreq)							std::cerr << "Codon freq.:   Ignore" << std::endl;
		if(cmd.mDoNotReduceForest)					std::cerr << "Reduce forest: Do not reduce" << std::endl;
		else if(cmd.mNoAggressiveStep)				std::cerr << "Reduce forest: Normal" << std::endl;
		else                    					std::cerr << "Reduce forest: Aggressive" << std::endl;
		if(cmd.mTimesFromFile)						std::cerr << "Times:         From tree file" << std::endl;
		if(cmd.mNoMaximization)						std::cerr << "Maximization:  No" << std::endl;
		if(cmd.mExportComputedTimes != UINT_MAX)	std::cerr << "Graph times:   From H" << cmd.mExportComputedTimes << std::endl;
		if(cmd.mTrace)								std::cerr << "Trace:         On" << std::endl;
		if(cmd.mGraphFile)							std::cerr << "Graph file:    " << cmd.mGraphFile << std::endl;
#ifdef _OPENMP
		if(!cmd.mForceSerial)						std::cerr << "Num. threads:  " << omp_get_max_threads() << std::endl;
#endif
		std::cerr << std::endl;
	}

	// Initialize the random number generator (0 means it is not set on the command line)
	if(cmd.mSeed == 0) cmd.mSeed = (unsigned int)time(NULL);
	srand(cmd.mSeed);

	// Start a timer (to measure serial part over parallel one)
	Timer timer;
	if(cmd.mVerboseLevel >= 2) timer.start();

	// Load the genes
	Genes g(cmd.mVerboseLevel);
	g.loadGenesFile(cmd.mGeneFile);

	// Load the phylogenetic tree
	PhyloTree t(cmd.mVerboseLevel);
	t.loadTree(cmd.mTreeFile);

	// Create and load the forest
	Forest forest(cmd.mVerboseLevel);
	forest.loadTreeAndGenes(t, g, cmd.mIgnoreFreq);

	// Reduce the forest merging common subtrees. Add also more reduction
	if(!cmd.mDoNotReduceForest)
	{
		forest.reduceSubtrees();
		if(!cmd.mNoAggressiveStep) forest.addAggressiveReduction();
	}

	// Subdivide the trees in groups based on dependencies
	forest.groupByDependency(cmd.mForceSerial);

	// Get the time needed by the serial part
	if(cmd.mVerboseLevel >= 2) {timer.stop(); std::cerr << "TIMER: " << timer.get() << std::endl;}

	// Print few statistics
	if(cmd.mVerboseLevel >= 2) std::cerr << forest;

	// Start timing parallel part
	if(cmd.mVerboseLevel >= 2) timer.start();

	// For all internal branches (only if a single branch has not been requested i.e. mBranch == UINT_MAX)
	size_t num_branches = forest.getNumInternalBranches();
	unsigned int branch_start = (cmd.mBranch < num_branches) ? cmd.mBranch   : 0;
	unsigned int branch_end   = (cmd.mBranch < num_branches) ? cmd.mBranch+1 : num_branches;
	for(unsigned int fg_branch=branch_start; fg_branch < branch_end; ++fg_branch)
	{
		if(cmd.mVerboseLevel >= 1) std::cerr << std::endl << "Doing branch " << fg_branch << std::endl;

		// Compute the null model maximum loglikelihood
		BranchSiteModelNullHyp h0(t.getNumBranches(), cmd.mSeed);
		double lnl0 = h0.computeModel(forest, fg_branch, cmd.mNoMaximization, cmd.mTimesFromFile, cmd.mTrace);

		// Compute the alternate model maximum loglikelihood
		BranchSiteModelAltHyp  h1(t.getNumBranches(), cmd.mSeed);
		double lnl1 = h1.computeModel(forest, fg_branch, cmd.mNoMaximization, cmd.mTimesFromFile, cmd.mTrace);

		if(cmd.mVerboseLevel >= 1)
		{
			std::cerr << std::endl;
			std::cerr << "LnL0: ";
			std::cerr << lnl0;
			std::cerr << " LnL1: ";
			std::cerr << lnl1;
			std::cerr << std::endl;
		}

		// If requested set the time in the forest and export to a graph visualization tool
		if(cmd.mGraphFile)
		{
			switch(cmd.mExportComputedTimes)
			{
			case 0:
				h0.saveComputedTimes(forest);
				break;

			case 1:
				h1.saveComputedTimes(forest);
				break;
			}
			forest.exportForest(cmd.mGraphFile, fg_branch);
		}

		// If the run passes the LRT, then compute the BEB
		// LRT test: -2*(lnl0-lnl1) > 3.84
		if(lnl1 - lnl0 > 1.920729)
		{
			BayesTest bt(g.getNumSites());
			bt.computeBEB();
			bt.printPositiveSelSites(fg_branch);
		}
	}

	// Get the time needed by the parallel part
	if(cmd.mVerboseLevel >= 2) {timer.stop(); std::cerr << "TIMER: " << timer.get() << std::endl;}

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
/// The first letters specify the kind of constant (like: STS_ status, OPT_ option value).
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
/// Counters should be unsigned int.
///

/// @page cmd_page Command Line switches
///
/// @verbatim
///
/// Usage:
///     FastCodeML [options] tree_file gene_file
/// 
/// -d  --debug  -v  --verbose (required argument)
///         Verbosity level
/// 
/// -q  --quiet (no argument)
///         No messages except results
/// 
/// -?  -h  --help (no argument)
///         This help
/// 
/// -s  --seed (required argument)
///         Random number generator seed
/// 
/// -b  --branch (required argument)
///         Do only this branch as foreground branch
/// 
/// -i  --ignore-freq (no argument)
///         Ignore computed codon frequency and set all to 1/61
/// 
/// -e  --export (required argument)
///         Export forest in GML format (if \%03d or \@03d is present, one is created for each fg branch)
/// 
/// -nr  --no-reduce (no argument)
///         Do not reduce forest by merging common subtrees
/// 
/// -l  --times-from-file (no argument)
///         Initial branch lengths from tree file
/// 
/// -o  --initial-step (no argument)
///         Only the initial step is performed (no maximization)
/// 
/// -c  --export-comp-times (required argument)
///         Export the computed times from H0 if 0, H1 if 1, otherwise the one read in the phylo tree
/// 
/// -r  --trace (no argument)
///         Trace the maximization run
/// 
/// -na  --no-aggressive (no argument)
///         Don't apply aggressive forest reduction
/// 
/// -np  --no-parallel (no argument)
///         Don't use parallel execution
/// 
/// @endverbatim
