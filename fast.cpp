/// test
/// @mainpage FastCodeML
///
/// @section intro_sect Introduction
///
/// FastCodeML is a rewrite of CodeML based directly on the pseudocode document.
/// It incorporates various parallelization strategies to be able to exploit
/// modern HPC machines architecture.
/// For this reason there are various parts of the code that can be selected at
/// compile time or run time to experiment with various, possible solutions.
///
/// @section contacts_sect Contacts
///
/// Contact us if you want more information on the project, want to collaborate
/// or suggest new ideas.
///
///- Ing. <a href="mailto:mvalle@cscs.ch">Mario Valle</a> - Swiss National
///Supercomputing Centre (CSCS) - Switzerland
///- The HP2C <a href="mailto:selectome@hp2c.ch">Selectome</a> Project Group -
///Mainly based in University of Lausanne - Switzerland
///

#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include "CmdLine.h"
#include "Newick.h"
#include "Phylip.h"
#include "BayesTest.h"
#include "Forest.h"
#include "Exceptions.h"
#include "BranchSiteModel.h"
#include "ParseParameters.h"
#include "VerbosityLevels.h"
#include "WriteResults.h"
#include "MathSupport.h"

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
///		@author Mario Valle - Swiss National Supercomputing Centre
///(CSCS)
///		@date 2010-12-22 (initial version)
///		@version 1.1
///
///	@param[in] aRgc Number of command line parameters
/// @param[in] aRgv Command line parameters
///
const char *version = "1.3.0";

int main(int aRgc, char **aRgv) {

	Timer timer_app;

	timer_app.start();

	try {

#ifdef USE_MKL_VML
		// If used, intitialize the MKL VML library
		vmlSetMode(VML_HA | VML_DOUBLE_CONSISTENT);
#endif

#ifdef USE_MPI
		// Start the high level parallel executor (based on MPI)
		HighLevelCoordinator hlc(&aRgc, &aRgv);
#endif

		// Parse the command line
		CmdLine cmd;
		cmd.parseCmdLine(aRgc, aRgv);

// Adjust and report the number of threads that will be used
#ifdef _OPENMP
		int num_threads = omp_get_max_threads();
		//std::cout<<"max num of thr: "<< num_threads <<std::endl;

		if ((cmd.mNumThreads >= 1) &&
				(cmd.mNumThreads <= (unsigned int)num_threads))
		num_threads = cmd.mNumThreads;
		//else
		//num_threads = num_threads/2 + 1;// to prevent possible false sharing when it is set to maximum
		// std::cout<<"num of thr: "<< num_threads <<std::endl;

		omp_set_num_threads(num_threads);
		//std::cout<<"current num of thr: "<< num_threads <<", command line num of thr: "<<cmd.mNumThreads<<std::endl;
		/*if (num_threads < 2)
		 cmd.mForceSerial = true;
		 else
		 cmd.mForceSerial = false;*/
#else
		cmd.mNumThreads = 1;
		int num_threads = 1;
		cmd.mForceSerial = true;
#endif

		/*#ifdef _OPENMP
		 int num_threads = omp_get_max_threads();
		 if(num_threads < 2 || cmd.mForceSerial)
		 {
		 cmd.mForceSerial = true;
		 num_threads = 1;
		 omp_set_num_threads(1);
		 }
		 #else
		 cmd.mForceSerial = true;
		 int num_threads = 1;
		 #endif*/

#ifdef USE_MPI
		// Shutdown messages from all MPI processes except the master
		if (!hlc.isMaster())
		cmd.mVerboseLevel = VERBOSE_NONE;
#endif

		//	  std::cout <<std::endl<<"------------------"<< std::endl<<"FastCodeML
		//V"<<version<<std::endl<<"------------------"<<std::endl;
		// Write out command line parameters (if not quiet i.e. if verbose level >
		// 0)
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {

			std::cout << "------------------------------------" << std::endl;
			std::cout << "FastCodeML V" << version << std::endl;
			std::cout << "------------------------------------" << std::endl;
			std::cout << std::endl;
			std::cout << "Tree file:	  " << cmd.mTreeFile << std::endl;
			std::cout << "Gene file:	  " << cmd.mGeneFile << std::endl;
			std::cout << "Verbose level:  " << cmd.mVerboseLevel << " ("
					<< decodeVerboseLevel(cmd.mVerboseLevel) << ')'
					<< std::endl;
			if (cmd.mSeed)
				std::cout << "Seed:			  " << cmd.mSeed << std::endl;
			if (cmd.mBranchFromFile)
				std::cout << "Branch:		  From tree file" << std::endl;
			else if (cmd.mBranchAll)
				std::cout << "FG Branches:	  All (internals + leaves) "
						<< std::endl;
			// else if(cmd.mBranchStart != UINT_MAX && cmd.mBranchStart ==
			// cmd.mBranchEnd)
			//											std::cout
			//<< "Branch:		  " << cmd.mBranchStart << std::endl;
			// else if(cmd.mBranchStart != UINT_MAX && cmd.mBranchEnd == UINT_MAX)
			//											std::cout
			//<< "Branches:		  " << cmd.mBranchStart << "-end" << std::endl;
			// else if(cmd.mBranchStart != UINT_MAX && cmd.mBranchEnd != UINT_MAX)
			//											std::cout
			//<< "Branches:		  " << cmd.mBranchStart << '-' <<
			//cmd.mBranchEnd << std::endl;
			if (!cmd.mStopIfNotLRT)
				std::cout << "H0 pre stop:	  No" << std::endl;
			if (cmd.mIgnoreFreq)
				std::cout << "Codon freq.:	  Ignore" << std::endl;
			if (cmd.mDoNotReduceForest)
				std::cout << "Reduce forest:  Do not reduce" << std::endl;
			else
				std::cout << "Reduce forest:  Aggressive" << std::endl;
			if (cmd.mInitH0fromH1)
				std::cout << "Starting val.:  From H1" << std::endl;
			else if (cmd.mInitFromParams && cmd.mBranchLengthsFromFile)
				std::cout
						<< "Starting val.:  Times from tree file and params from "
								"const (see below)" << std::endl;
			else if (cmd.mInitFromParams)
				std::cout << "Starting val.:  Params from const (see below)"
						<< std::endl;
			else if (cmd.mBranchLengthsFromFile)
				std::cout << "Starting val.:  Times from tree file"
						<< std::endl;
			if (cmd.mNoMaximization)
				std::cout << "Maximization:	  No" << std::endl;
			if (cmd.mTrace)
				std::cout << "Trace:		  On" << std::endl;
			if (cmd.mCleanData)
				std::cout << "Clean data:	  On" << std::endl;
			else
				std::cout << "Clean data:	  Off" << std::endl;
			if (cmd.mGraphFile)
				std::cout << "Graph file:	  " << cmd.mGraphFile << std::endl;
			if (cmd.mGraphFile && cmd.mExportComputedTimes != UINT_MAX)
				std::cout << "Graph times:	  From H" << cmd.mExportComputedTimes
						<< std::endl;
			if (!cmd.mNoMaximization)
				std::cout << "Optimizer:	  " << cmd.mOptimizationAlgo
						<< std::endl;
			if (cmd.mMaxIterations != MAX_ITERATIONS)
				std::cout << "Max iterations: " << cmd.mMaxIterations
						<< std::endl;
			if (cmd.mDeltaValueForGradient > 0.0)
				std::cout << "Delta value:	  " << cmd.mDeltaValueForGradient
						<< std::endl;
			std::cout << "Relative error: " << cmd.mRelativeError << std::endl;
			if (cmd.mResultsFile)
				std::cout << "Results file:	  " << cmd.mResultsFile
						<< std::endl;
			if (cmd.mNumThreads)
				std::cout << "Number of threads: " << cmd.mNumThreads
						<< std::endl;
			if (cmd.mFixedBranchLength)
				std::cout << "Branch lengths are fixed" << std::endl;
#ifdef _OPENMP
			if (num_threads > 1) {
				std::cout << "Current num. threads:	  " << num_threads << std::endl
				<< "Total num. cores:	  " << omp_get_num_procs() << std::endl;
			} else
#endif
			{
				std::cout << "Num. threads:	  1 serial" << std::endl
						<< "Num. cores:	  1" << std::endl;
			}
#ifdef USE_MPI
			if (hlc.numJobs() > 2)
			std::cout << "Num. MPI proc:  1 (master) + " << hlc.numJobs() - 1
			<< " (workers)" << std::endl;
			else
			std::cout << "Num. MPI proc:  Insufficient, single task execution"
			<< std::endl;
#endif
			std::cout << "Compiled with:  ";
#ifdef _OPENMP
			std::cout << "USE_OPENMP ";
#endif
#ifdef USE_MPI
			std::cout << "USE_MPI ";
#endif
#ifdef USE_CPV_SCALING
			std::cout << "USE_CPV_SCALING ";
#endif
#ifdef NEW_LIKELIHOOD
			std::cout << "NEW_LIKELIHOOD ";
#endif
#ifdef NON_RECURSIVE_VISIT
			std::cout << "NON_RECURSIVE_VISIT ";
#endif
#ifdef USE_DAG
			std::cout << "USE_DAG ";
#endif
#ifdef USE_ORIGINAL_PROPORTIONS
			std::cout << "USE_ORIGINAL_PROPORTIONS ";
#endif
#ifdef USE_LAPACK
			std::cout << "USE_LAPACK ";
#endif
#ifdef USE_MKL_VML
			std::cout << "USE_MKL_VML";
#endif
			std::cout << std::endl << std::endl;
			if (cmd.mInitFromParams) {
				std::cout << "Param initial values:" << std::endl << std::endl
						<< ParseParameters::getInstance();
			}
		}

// Initialize the random number generator (0 means it is not set on the command
// line)
#ifdef USE_MPI
		// Insure that each MPI process starts with a different seed
		if (cmd.mSeed == 0)
		cmd.mSeed = static_cast<unsigned int>(time(NULL)) +
		static_cast<unsigned int>(hlc.getRank()) * 1000;
#else
		if (cmd.mSeed == 0)
			cmd.mSeed = static_cast<unsigned int>(time(NULL));
#endif
		// srand(cmd.mSeed); // fastcodeml seed
		SetSeedCodeml(cmd.mSeed, 0); // codeml seed is 1

		// Verify the optimizer algorithm selected on the command line
		if (!cmd.mNoMaximization)
			BranchSiteModel::verifyOptimizerAlgo(cmd.mOptimizationAlgo);

		// Start a timer (to measure serial part over parallel one)
		Timer timer;
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT)
			timer.start();

		// Create the forest
		Forest forest(cmd.mVerboseLevel);

		// Enclose file loading into a block so temporary structures could be
		// deleted when no more needed
		//{
		// Load the multiple sequence alignment (MSA)
		Phylip msa(cmd.mVerboseLevel);
		msa.readFile(cmd.mGeneFile, cmd.mCleanData);

		// Load the phylogenetic tree
		Newick tree(cmd.mVerboseLevel);
		tree.readFile(cmd.mTreeFile);

		// Check coherence between the two files
		msa.checkNameCoherence(tree.getSpecies());

		// Check root and unrooting if tree is rooted
		tree.checkRootBranches();

		// If times from file then check for null branch lengths for any leaf
		if (cmd.mBranchLengthsFromFile) {
			int zero_on_leaf_cnt = 0;
			int zero_on_int_cnt = 0;
			tree.countNullBranchLengths(zero_on_leaf_cnt, zero_on_int_cnt);

			if (zero_on_leaf_cnt > 0 || zero_on_int_cnt > 0) {
				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
					std::cout
							<< "Found null or missing branch length in tree file: on "
							<< zero_on_leaf_cnt << " leave(s) and on "
							<< zero_on_int_cnt << " internal branch(es)."
							<< std::endl;
				}
			}
		}

		// Print the tree with the numbering of internal branches
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT)
			tree.printTreeAnnotated(std::cout);

		// Load the forest
		forest.loadTreeAndGenes(tree, msa,
				cmd.mIgnoreFreq ?
						CodonFrequencies::CODON_FREQ_MODEL_UNIF :
						CodonFrequencies::CODON_FREQ_MODEL_F3X4);

		// Reduce the forest merging common subtrees. Add also more reduction, then
		// clean the no more useful data.
		if (!cmd.mDoNotReduceForest) {
			// bool sts = forest.reduceSubtrees(cmd.mNumReductionBlocks);
			forest.reduceSubtrees();

#ifndef NEW_LIKELIHOOD
			forest.addAggressiveReduction();
#endif
			forest.cleanReductionWorkingData();
#ifdef NEW_LIKELIHOOD
			forest.prepareNewReduction();
#endif
		}
#ifdef NEW_LIKELIHOOD
		else {
			forest.prepareNewReductionNoReuse();
		}
#endif

#ifdef NON_RECURSIVE_VISIT
		// Prepare the pointers to visit the trees without recursion
		forest.prepareNonRecursiveVisit();
#endif

// Subdivide the trees in groups based on dependencies
// forest.prepareDependencies(cmd.mForceSerial || cmd.mDoNotReduceForest);

#ifdef USE_DAG
		// Load the forest into a DAG
		forest.loadForestIntoDAG(Nt);
#endif

		// Get the time needed by data preprocessing
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {
			timer.stop();
			std::cout << std::endl << "TIMER (preprocessing) ncores: "
					<< std::setw(2) << num_threads << " time: " << timer.get()
					<< std::endl;
		}

		// Print few statistics
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT)
			std::cout << forest;

#ifdef USE_MPI
		// Distribute the work. If run under MPI then finish, else return to the
		// standard execution flow
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT)
		timer.start();
		bool has_run_under_MPI = hlc.startWork(forest, cmd);

		// If executed under MPI report the time spent, otherwise stop the timer so
		// it can be restarted around the serial execution
		if (has_run_under_MPI) {
			if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {
				timer.stop();
				std::cout << std::endl
				<< "TIMER (processing) ncores: " << std::setw(2)
				<< num_threads * (hlc.numJobs() - 1) + 1
				<< " time: " << timer.get() << std::endl;
			}
			return 0;
		} else {
			timer.stop();
		}
#endif

		// Compute the range of branches to mark as foreground
		size_t branch_start, branch_end;

		std::set<int> fg_set;     // to save a list of fg branches from the
								  // getBranchRange function
		std::set<int> ib_set;    // to save a list of internal branches from the
								 // getBranchRange function
		std::vector<double> mVar; // to save optimization variables

		forest.getBranchRange(cmd, branch_start, branch_end, fg_set, ib_set); // fgset is added to save a list of fg branches

		// for (std::set<int>::iterator it=ib_set.begin(); it!=ib_set.end(); ++it)
		// std::cout << " " << *it << ",";

		// Start timing parallel part
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT)
			timer.start();

		double lnl0, lnl1 = 0.;

		if (!fg_set.empty()) // in case of marked fg branches (one or multiple fg)
		{
			//initial the output results object
			WriteResultsMfg output_results_mfg(cmd.mResultsFile);
			
			if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
				std::cout << std::endl
						<< "Doing foreground branch(es) from tree file "
						<< std::endl;
			if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
				std::cout << "-------------------------------------------"
						<< std::endl;
			if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
				std::cout << "Doing foreground branch(es) ";
				for (std::set<int>::iterator it = fg_set.begin();
						it != fg_set.end(); ++it)
					std::cout << " " << *it << " ";
				std::cout << std::endl;
			}
			// Initialize the models
			MfgBranchSiteModelNullHyp h0(forest, cmd);
			MfgBranchSiteModelAltHyp h1(forest, cmd);

			// Initialize the test
			MfgBayesTest beb(forest, cmd.mVerboseLevel, cmd.mDoNotReduceForest);

			// Compute the alternate model maximum loglikelihood
			double lnl1 = 0.;
			if (cmd.mComputeHypothesis != 0) {
				if (cmd.mInitFromParams)
					h1.initFromParams();
				if (cmd.mBranchLengthsFromFile)
					h1.initFromTree();

				lnl1 = h1(fg_set);
				// h1.saveComputedTimes();

				// std::cout << "lnl1 = " << lnl1 << std::endl;
				// Save the value for formatted output
				 output_results_mfg.saveLnL(fg_set, lnl1, 1);
			}

			// Compute the null model maximum loglikelihood
			double lnl0 = 0.;
			if (cmd.mComputeHypothesis != 1) {
				if (cmd.mInitH0fromH1)
					h0.initFromResult(h1.getVariables());
				else {
					if (cmd.mInitFromParams)
						h0.initFromParams();
					if (cmd.mBranchLengthsFromFile)
						h0.initFromTree();
				}

				lnl0 = h0(fg_set,
						cmd.mStopIfNotLRT && cmd.mComputeHypothesis != 0,
						lnl1 - THRESHOLD_FOR_LRT);

				// std::cout << "lnl0 = " << lnl0 << std::endl;
				// Save the value for formatted output (only if has not be forced to
				// stop)
				 if(lnl0 < DBL_MAX) output_results_mfg.saveLnL(fg_set, lnl0, 0);
			}

			if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
				std::cout << std::endl;
				if (cmd.mComputeHypothesis != 1) {
					std::cout << "LnL0: ";
					if (lnl0 == std::numeric_limits<double>::infinity())
						std::cout << "**Invalid result**";
					else if (lnl0 < DBL_MAX)
						std::cout << std::setprecision(15) << std::fixed
								<< lnl0;
					else
						std::cout << "(Doesn't pass LRT, skipping)";
					std::cout << " Function calls: " << h0.getNumEvaluations()
							<< "	  ";
					std::cout << std::endl << std::endl;
					if (lnl0 != std::numeric_limits<double>::infinity()) {
						std::string s0 = h0.printFinalVars(std::cout);
						// std::cout<<"EDW0: "<< s0 <<std::endl;
						output_results_mfg.saveParameters(fg_set, s0, 0);
					}
					std::cout << std::endl;
				}
				if (cmd.mComputeHypothesis != 0) {
					std::cout << "LnL1: ";
					if (lnl1 == std::numeric_limits<double>::infinity())
						std::cout << "**Invalid result**";
					else
						std::cout << std::setprecision(15) << std::fixed
								<< lnl1;
					std::cout << " Function calls: " << h1.getNumEvaluations();
					std::cout << std::endl << std::endl;
					if (lnl1 != std::numeric_limits<double>::infinity()) {
						std::string s1 = h1.printFinalVars(std::cout);
						// std::cout<<"EDW1: "<< s1 <<std::endl;
					        output_results_mfg.saveParameters(fg_set, s1, 1);
					}
					std::cout << std::endl;
				}
				if (cmd.mComputeHypothesis > 1) {
					if (lnl0 == std::numeric_limits<double>::infinity()
							|| lnl1 == std::numeric_limits<double>::infinity())
						std::cout << "LRT: **Invalid result**";
					else if (lnl0 < DBL_MAX)
						std::cout << "LRT: " << std::setprecision(15)
								<< std::fixed << lnl1 - lnl0 << "	 (threshold: "
								<< std::setprecision(15) << std::fixed
								<< THRESHOLD_FOR_LRT << ')';
					else
						std::cout << "LRT: < " << std::setprecision(15)
								<< std::fixed << THRESHOLD_FOR_LRT;
					std::cout << std::endl;
				}
			}

			// If requested set the time in the forest and export to a graph
			// visualization tool
			if (cmd.mGraphFile) {
				switch (cmd.mExportComputedTimes) {
				case 0:
					h0.saveComputedTimes();
					break;

				case 1:
					h1.saveComputedTimes();
					break;

				default:
					break;
				}

				// Use the forest export class
				ForestExport fe(forest);
				fe.exportForest(cmd.mGraphFile, 0);
			}

			// If the two hypothesis are computed, H0 has not been stopped and the run
			// passes the LRT, then compute the BEB
			if (cmd.mComputeHypothesis > 1 && lnl0 < DBL_MAX
					&& BranchSiteModel::performLRT(lnl0, lnl1)) {
				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
					std::cout << std::endl
							<< "LRT is significant. Computing sites under positive "
									"selection ... " << std::endl;

				// Get the scale values from the latest optimized h1.
				std::vector<double> scales(2);
				h1.getScales(scales);

				// Run the BEB test
				// if(cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) std::cout << std::endl
				// << "LRT is significant. Computing sites under positive selection ...
				// " << std::endl ;
				beb.computeBEB(h1.getVariables(), fg_set, scales);

				// Output the sites under positive selection (if any)
				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
					beb.printPositiveSelSites(fg_set);

				// Get the sites under positive selection for printing in the results
				// file (if defined)

				if(output_results_mfg.isWriteResultsEnabled())
				 {
				 std::vector<unsigned int> positive_sel_sites;
				 std::vector<double> positive_sel_sites_prob;
				 beb.extractPositiveSelSites(positive_sel_sites,
				 positive_sel_sites_prob);

				 //if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
				 //std::cout << std::endl
				 //<< "Positively selected sites and their probabilities : ";
				 //std::cout << std::endl;
				 //for (std::vector<unsigned int>::iterator it =
				 //positive_sel_sites.begin();
				 //it != positive_sel_sites.end(); ++it)
				 //std::cout << " " << *it << ",";
				 //std::cout << std::endl;
				 //for (std::vector<double>::iterator it =
				 //positive_sel_sites_prob.begin();
				 //it != positive_sel_sites_prob.end(); ++it)
				 //std::cout << " " << *it << ",";
				 //std::cout << std::endl;
				 //}

				 output_results_mfg.savePositiveSelSites(fg_set, positive_sel_sites,
				 positive_sel_sites_prob);
				 }
			}

			// if branches are fixed

			if (cmd.mFixedBranchLength)

			{
				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
					std::cout << std::endl << "Final ";
					tree.printTreeAnnotated(std::cout, NULL, 0, true);
					std::cout << std::endl;
					tree.printTreeAnnotated(std::cout, NULL, 0, true, false);
				}
			}

			else

			{
				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
					if (cmd.mComputeHypothesis == 0)
						mVar = h0.getVariables();
					else
						mVar = h1.getVariables();

					if (cmd.mComputeHypothesis == 0)
						std::cout << std::endl << "H0 Final ";
					else
						std::cout << std::endl << "H1 Final ";
					tree.printTreeAnnotatedWithEstLens(std::cout, NULL, 0, true,
							&mVar);
					std::cout << std::endl;
					tree.printTreeAnnotatedWithEstLens(std::cout, NULL, 0, true,
							&mVar, false);
				}
				// tree.printTreeAnnotatedWithEstLens(std::cout, NULL, 0, true,
				// &h1.getVariables());
				// std :: cout << std::endl;
			}

			/*if(cmd.mInitFromParams)			h0.initFromParams();
			 if(cmd.mBranchLengthsFromFile)	h0.initFromTree();

			 lnl0 = h0(fg_set, cmd.mStopIfNotLRT && cmd.mComputeHypothesis != 0,
			 0-THRESHOLD_FOR_LRT);
			 std::cout << "lnl0 (multiple fg) = " << lnl0 << std::endl;

			 if(cmd.mInitFromParams)			h1.initFromParams();
			 if(cmd.mBranchLengthsFromFile)	h1.initFromTree();

			 lnl1 = h1(fg_set);
			 std::cout << "lnl1 (multiple fg) = " << lnl1 << std::endl;*/

			timer_app.stop();
			if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
				std::cout << std::endl << "Time used: " << timer_app.get() / 60000
									<< "m:" << (timer_app.get() / 1000) % 60 << "s" << std::endl;
				std::cout << "Cores used: " << num_threads << std::endl;
			}
			
			// Output the results
		
			 output_results_mfg.outputResults(fg_set);
		
			return 0;
		}

		// Else for all requested internal branches

		// Initialize the output results file (if the argument is null, no file is
		// created)
		WriteResults output_results(cmd.mResultsFile);
		

		// Initialize the models
		BranchSiteModelNullHyp h0(forest, cmd);
		BranchSiteModelAltHyp h1(forest, cmd);

		// Initialize the test
		BayesTest beb(forest, cmd.mVerboseLevel, cmd.mDoNotReduceForest);

		if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
			if (cmd.mBranchAll)
				std::cout << std::endl << "Doing for all branches (internal + leaf) as foreground"
						<< std::endl;
			else
				std::cout << std::endl << "Doing for all internal branches as foreground"
						<< std::endl;
			std::cout << "------------------------------------" << std::endl;
		}

		for (size_t fg_branch = branch_start; fg_branch <= branch_end;
				++fg_branch) {

			if (cmd.mBranchAll
					or (!cmd.mBranchAll
							and ib_set.find(fg_branch) != ib_set.end()))

					{

				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
					std::cout << "Doing foreground branch " << fg_branch
							<< std::endl;

				// Compute the alternate model maximum loglikelihood
				double lnl1 = 0.;
				if (cmd.mComputeHypothesis != 0) {
					if (cmd.mInitFromParams)
						h1.initFromParams();
					if (cmd.mBranchLengthsFromFile)
						h1.initFromTree();

					lnl1 = h1(fg_branch);
					// h1.saveComputedTimes();
					// h1.mBranches

					// Save the value for formatted output
					 output_results.saveLnL(fg_branch, lnl1, 1);
				}

				// Compute the null model maximum loglikelihood
				double lnl0 = 0.;
				if (cmd.mComputeHypothesis != 1) {
					if (cmd.mInitH0fromH1)
						h0.initFromResult(h1.getVariables());
					else {
						if (cmd.mInitFromParams)
							h0.initFromParams();
						if (cmd.mBranchLengthsFromFile)
							h0.initFromTree();
					}

					lnl0 = h0(fg_branch,
							cmd.mStopIfNotLRT && cmd.mComputeHypothesis != 0,
							lnl1 - THRESHOLD_FOR_LRT);

					// Save the value for formatted output (only if has not be forced to
					// stop)
					 if(lnl0 < DBL_MAX) output_results.saveLnL(fg_branch, lnl0, 0);
				}

				if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
					std::cout << std::endl;
					if (cmd.mComputeHypothesis != 1) {
						std::cout << "LnL0: ";
						if (lnl0 == std::numeric_limits<double>::infinity())
							std::cout << "**Invalid result**";
						else if (lnl0 < DBL_MAX)
							std::cout << std::setprecision(15) << std::fixed
									<< lnl0;
						else
							std::cout << "(Doesn't pass LRT, skipping)";
						std::cout << " Function calls: "
								<< h0.getNumEvaluations() << "	  ";
						std::cout << std::endl << std::endl;
						if (lnl0 != std::numeric_limits<double>::infinity()) {
							std::string s0 = h0.printFinalVars(std::cout);
							output_results.saveParameters(fg_branch, s0, 0);
						}
						std::cout << std::endl;
					}
					if (cmd.mComputeHypothesis != 0) {
						std::cout << "LnL1: ";
						if (lnl1 == std::numeric_limits<double>::infinity())
							std::cout << "**Invalid result**";
						else
							std::cout << std::setprecision(15) << std::fixed
									<< lnl1;
						std::cout << " Function calls: "
								<< h1.getNumEvaluations();
						std::cout << std::endl << std::endl;
						if (lnl1 != std::numeric_limits<double>::infinity()) {
							std::string s1 = h1.printFinalVars(std::cout);
							output_results.saveParameters(fg_branch, s1, 1);
						}
						std::cout << std::endl;
					}
					if (cmd.mComputeHypothesis > 1) {
						if (lnl0 == std::numeric_limits<double>::infinity()
								|| lnl1
										== std::numeric_limits<double>::infinity())
							std::cout << "LRT: **Invalid result**";
						else if (lnl0 < DBL_MAX)
							std::cout << "LRT: " << std::setprecision(15)
									<< std::fixed << lnl1 - lnl0
									<< "	 (threshold: " << std::setprecision(15)
									<< std::fixed << THRESHOLD_FOR_LRT << ')';
						else
							std::cout << "LRT: < " << std::setprecision(15)
									<< std::fixed << THRESHOLD_FOR_LRT;
						std::cout << std::endl;
					}
				}

				// If requested set the time in the forest and export to a graph
				// visualization tool
				if (cmd.mGraphFile) {
					switch (cmd.mExportComputedTimes) {
					case 0:
						h0.saveComputedTimes();
						break; 

					case 1:
						h1.saveComputedTimes();
						break;

					default:
						break;
					}

					// Use the forest export class
					ForestExport fe(forest);
					fe.exportForest(cmd.mGraphFile, fg_branch);
				}

				// If the two hypothesis are computed, H0 has not been stopped and the
				// run passes the LRT, then compute the BEB
				if (cmd.mComputeHypothesis > 1 && lnl0 < DBL_MAX
						&& BranchSiteModel::performLRT(lnl0, lnl1)) {
					if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
						std::cout << std::endl
								<< "LRT is significant. Computing sites under positive "
										"selection ... " << std::endl;

					// Get the scale values from the latest optimized h1.
					std::vector<double> scales(2);
					h1.getScales(scales);

					// Run the BEB test
					// if(cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) std::cout <<
					// std::endl << "LRT is significant. Computing sites under positive
					// selection ... " << std::endl ;
					beb.computeBEB(h1.getVariables(), fg_branch, scales);

					// Output the sites under positive selection (if any)
					if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
						beb.printPositiveSelSites(fg_branch);

					// Get the sites under positive selection for printing in the results
					// file (if defined)
					if (output_results.isWriteResultsEnabled()) {
						std::vector<unsigned int> positive_sel_sites;
						std::vector<double> positive_sel_sites_prob;
						beb.extractPositiveSelSites(positive_sel_sites,
								positive_sel_sites_prob);
						output_results.savePositiveSelSites(fg_branch,
								positive_sel_sites, positive_sel_sites_prob);
					}
				}

				if (cmd.mFixedBranchLength)

				{
					if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
						std::cout << std::endl << "Final ";
						tree.printTreeAnnotated(std::cout, NULL, 0, true);
						std::cout << std::endl;
						tree.printTreeAnnotated(std::cout, NULL, 0, true, false);
						std::cout << std::endl;
					}
				}

				else {
					// std :: cout << std::endl;
					if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
						if (cmd.mComputeHypothesis == 0)
							mVar = h0.getVariables();
						else
							mVar = h1.getVariables();

						if (cmd.mComputeHypothesis == 0)
							std::cout << std::endl << "H0 Final ";
						else
							std::cout << std::endl << "H1 Final ";
						tree.printTreeAnnotatedWithEstLens(std::cout, NULL, 0,
								true, &mVar);
						std::cout << std::endl;
						tree.printTreeAnnotatedWithEstLens(std::cout, NULL, 0,
														true, &mVar, false);
						std::cout << std::endl;
					}
				}
			}
		}

		// Get the time needed by the parallel part
		if (cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {
			timer.stop();
			std::cout << std::endl << "TIMER (processing) ncores: "
					<< std::setw(2) << num_threads << " time: " << timer.get()
					<< std::endl;
		}

		timer_app.stop();
		if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) {
			std::cout << std::endl << "Time used: " << timer_app.get() / 60000
					<< "m:" << (timer_app.get() / 1000) % 60 << "s" << std::endl;
			std::cout << "Cores used: " << num_threads << std::endl;
		}

		// Output the results
		output_results.outputResults(ib_set,cmd.mBranchAll);

		// Catch all exceptions
	} catch (const FastCodeMLSuccess &) {
		return 0;
	} catch (const FastCodeMLFatal &e) {
		// If a message associated (i.e. no empty string), display it
		if (e.what()[0])
			std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	} catch (const FastCodeMLMemoryError &e) {
		std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	} catch (const std::bad_alloc &e) {
		std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	} catch (...) {
		std::cout << std::endl << "Default exception caught." << std::endl
				<< std::endl;
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
/// Only very common and specific names should be all lowercase, like read,
/// clean, size.
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
/// The first letters specify the kind of constant (like: STS_ for status, OPT_
/// for option value).
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
/// Array sizes and corresponding indexes should be size_t. The remaining
/// counters should be unsigned int.
///
/// The null pointer should be written as NULL, not 0 to make clear its purpose.
///
/// @page vampir_page Using Vampir for profiling
/// On Linux we use VampirTrace to collect profile data and Vampir to display
/// the results (http://www.vampir.eu/).
///
/// Before invoking CMAKE define CXX=vtCC
///
/// Define CMAKE_BUILD_TYPE as: RelWithDebInfo
///
/// Run CMAKE and configure.
///
/// Then define CMAKE_CXX_FLAGS as: -vt:mt -vt:noopari
/// If you want to trace also the OpenMP calls then change it to: -vt:mt
/// -vt:preprocess -DVTRACE
///
/// Then proceed as usual to build the executable.
///
/// Before running the executable, define the following environment variables:
///
///		export VT_BUFFER_SIZE=512M
///		export VT_MAX_FLUSHES=0
///		export VT_SYNCH_FLUSH=yes
///		export VT_GPUTRACE=no
///		export VT_UNIFY=no
///
/// Due to a VampirTrace bug, at the end of the execution, run the vtunify
/// executable by itself.
///
/// Now you can analyze the results by running vampir on the *.otf file
/// generated.
///
