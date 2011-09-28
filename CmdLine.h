
#ifndef CMDLINE_H
#define CMDLINE_H

#include "simpleopt/SimpleOpt.h"

/// Parse the command line flags.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
///
class CmdLine
{
public:
	/// Constructor.
	///
	CmdLine();
	
	/// Parse the command line. On error throw an exception.
	///
	/// @param[in] aCnt Number of command line parameters (from main argc)
	/// @param[in] aVal Array of command line parameters (from main argv)
	///
	void parseCmdLine(int aCnt, char **aVal);


public:
	unsigned int	mVerboseLevel;			///< Verbosity level. 0: no messages; 1: basic messages; 2: messages useful for debugging; 3: really annoying
	unsigned int	mSeed;					///< Random number generator seed (0 means not set from command line)
	unsigned int	mBranch;				///< Do only this branch. The numbering starts at 0 (UINT_MAX means run all branches)
	unsigned int	mExportComputedTimes;	///< If 0 or 1 the times exported are the computed ones in H0 or H1, otherwise (UINT_MAX) the one read from file
	const char*		mTreeFile;				///< Newick tree file name
	const char*		mGeneFile;				///< %Genes file name
	const char*		mGraphFile;				///< If not null export the forest to this file in GML format to be visualized using R igraph package or yEd editor
	bool			mIgnoreFreq;			///< Ignore the computed codon frequencies and set them all to 1/61
	bool			mDoNotReduceForest;		///< If true do not reduce the forest merging common subtrees
	bool			mTimesFromFile;			///< The initial value of the branch lengths is taken from the phylo tree file
	bool			mNoMaximization;		///< Only the first step of the likelihood maximization is taken
	bool			mTrace;					///< Trace the optimization steps
	bool			mNoAggressiveStep;		///< Do not apply aggressive common subtree reduction
	bool			mForceSerial;			///< Disable all parallelism
	bool			mBranchFromFile;		///< Read the foreground branch to use from the Phylo Tree file (it is marked as #1)
	unsigned int	mComputeHypothesis;		///< If set to 0 compute only H0, if set to 1 compute only H1, otherwise compute both
	bool			mInitH1fromH0;			///< If set starts the H1 computation from the H0 result
	unsigned int	mOptimizationAlgo;		///< Select the optimization algorithm to use


private:

	/// Return the text corresponding to an error code.
	///
	/// @param[in] aOptParser The command line parser object
	///
	/// @return The human readable error message
	///
	const char *getLastErrorText(CSimpleOpt& aOptParser);

	/// Print the help about the parameters.
	///
	/// @param[in] aParserOptions The table of options definitions
	///
	void showHelp(const CSimpleOpt::SOption *aParserOptions);
};

#endif

