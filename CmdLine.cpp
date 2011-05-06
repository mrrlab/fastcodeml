
#include <iostream>
#include <climits>
#include "CmdLine.h"
#include "simpleopt/SimpleOpt.h"
#include "Exceptions.h"
#ifdef _OPENMP
#include <omp.h>
#endif

const char *CmdLine::getLastErrorText(CSimpleOpt& aOptParser)
{
    switch(aOptParser.LastError())
	{
    case SO_SUCCESS:            return "Success";
    case SO_OPT_INVALID:        return "Unrecognized option";
    case SO_OPT_MULTIPLE:       return "Option matched multiple strings";
    case SO_ARG_INVALID:        return "Option does not accept argument";
    case SO_ARG_INVALID_TYPE:   return "Invalid argument format";
    case SO_ARG_MISSING:        return "Required argument is missing";
    case SO_ARG_INVALID_DATA:   return "Invalid argument data";
    default:					return "Unknown error";
    }
}


void CmdLine::showHelp(const CSimpleOpt::SOption *aParserOptions)
{
	int i, j, cnt;

	// Count entries and create an indicator array
	for(cnt=0; aParserOptions[cnt].pszArg != NULL; ++cnt) {}
	bool *done = new bool[cnt];
	for(i=0; i < cnt; ++i) done[i] = false;

	// For each different option
	for(i=0; i < cnt; ++i)
	{
		if(done[i]) continue;
		done[i] = true;

		std::cerr << aParserOptions[i].pszArg;
		for(j=i+1; j < cnt; ++j)
		{
			if(done[j] || aParserOptions[j].nId != aParserOptions[i].nId) continue;
			done[j] = true;
			std::cerr << "  " << aParserOptions[j].pszArg;
		}

		// Translate the kind of argument
		const char* type = "";
		switch(aParserOptions[i].nArgType)
		{
		case SO_NONE:   
			type = "(no argument)";
			break;

		case SO_REQ_SEP:
		case SO_REQ_CMB:
			type = "(required argument)";
			break;

		case SO_OPT:
			type = "(optional argument)";
			break;

		case SO_MULTI:
			type = "(multiple arguments)";
			break;
		}

		std::cerr << " " << type << std::endl;
		std::cerr << "        " << aParserOptions[i].pszHelp << std::endl << std::endl;
	}

	delete [] done;
}


CmdLine::CmdLine()
{
	// Initialize options values
	mVerboseLevel			= 1;
	mGeneFile				= 0;
	mTreeFile				= 0;
	mSeed					= 0;
	mBranch					= UINT_MAX;
	mGraphFile				= 0;
	mIgnoreFreq				= false;
	mDoNotReduceForest		= false;
	mTimesFromFile			= false;
	mNoMaximization			= false;
	mExportComputedTimes	= UINT_MAX;
	mTrace					= false;
	mNoAggressiveStep		= false;
	mForceSerial			= false;
}


void CmdLine::parseCmdLine(int aCnt, char **aVal)
{
	// Setup the command line parser
	enum {
		OPT_VERBOSE,
		OPT_QUIET,
		OPT_HELP,
		OPT_SEED,
		OPT_BRANCH,
		OPT_IGNORE_FREQ,
		OPT_EXPORT,
		OPT_NOT_REDUCE,
		OPT_TIMES_FROM_FILE,
		OPT_ONE_STEP,
		OPT_COMP_TIMES,
		OPT_TRACE,
		OPT_NO_AGGRESSIVE,
		OPT_FORCE_SERIAL
	};

	CSimpleOpt::SOption parser_options[] = {
		{ OPT_VERBOSE,			"-d",					SO_REQ_SEP, "Verbosity level" },
		{ OPT_VERBOSE,			"--debug",				SO_REQ_SEP, "" },
		{ OPT_VERBOSE,			"-v",					SO_REQ_SEP, "" },
		{ OPT_VERBOSE,			"--verbose",			SO_REQ_SEP, "" },
		{ OPT_QUIET,			"-q",					SO_NONE,    "No messages except results" },
		{ OPT_QUIET,			"--quiet",				SO_NONE,    "" },
		{ OPT_HELP,				"-?",					SO_NONE,    "This help" },
		{ OPT_HELP,				"-h",					SO_NONE,    "" },
		{ OPT_HELP,				"--help",				SO_NONE,    "" },
		{ OPT_SEED,				"-s",					SO_REQ_SEP, "Random number generator seed" },
		{ OPT_SEED,				"--seed",				SO_REQ_SEP, "" },
		{ OPT_BRANCH,			"-b",					SO_REQ_SEP, "Do only this branch as foreground branch" },
		{ OPT_BRANCH,			"--branch",				SO_REQ_SEP, "" },
		{ OPT_IGNORE_FREQ,		"-i",					SO_NONE,	"Ignore computed codon frequency and set all to 1/61" },
		{ OPT_IGNORE_FREQ,		"--ignore-freq",		SO_NONE,	"" },
		{ OPT_EXPORT,			"-e",					SO_REQ_SEP,	"Export forest in GML format (if %03d or @03d is present, one is created for each fg branch)" },
		{ OPT_EXPORT,			"--export",				SO_REQ_SEP,	"" },
		{ OPT_NOT_REDUCE,		"-nr",					SO_NONE,	"Do not reduce forest by merging common subtrees" },
		{ OPT_NOT_REDUCE,		"--no-reduce",			SO_NONE,	"" },
		{ OPT_TIMES_FROM_FILE,	"-l",					SO_NONE,	"Initial branch lengths from tree file" },
		{ OPT_TIMES_FROM_FILE,	"--times-from-file",	SO_NONE,	"" },
		{ OPT_ONE_STEP,			"-o",					SO_NONE,	"Only the initial step is performed (no maximization)" },
		{ OPT_ONE_STEP,			"--initial-step",		SO_NONE,	"" },
		{ OPT_COMP_TIMES,		"-c",					SO_REQ_SEP,	"Export the computed times from H0 if 0, H1 if 1, otherwise the one read in the phylo tree" },
		{ OPT_COMP_TIMES,		"--export-comp-times",	SO_REQ_SEP,	"" },
		{ OPT_TRACE,			"-r",					SO_NONE,	"Trace the maximization run" },
		{ OPT_TRACE,			"--trace",				SO_NONE,	"" },
		{ OPT_NO_AGGRESSIVE,	"-na",					SO_NONE,	"Don't apply aggressive forest reduction" },
		{ OPT_NO_AGGRESSIVE,	"--no-aggressive",		SO_NONE,	"" },
		{ OPT_FORCE_SERIAL,		"-np",					SO_NONE,	"Don't use parallel execution" },
		{ OPT_FORCE_SERIAL,		"--no-parallel",		SO_NONE,	"" },
		SO_END_OF_OPTIONS
	};

	// Setup the usage string
	const char* usage_msg = "FastCodeML [options] tree_file gene_file";

    // Declare our options parser, pass in the arguments from main as well as our array of valid options.
    CSimpleOpt args(aCnt, aVal, parser_options, SO_O_NOSLASH);

    // While there are arguments left to process
    while(args.Next())
	{
		if(args.LastError() == SO_OPT_INVALID)
		{
			std::cerr << "Error: " << getLastErrorText(args) << ": " << args.OptionText() << std::endl;
			throw FastCodeMLFatalNoMsg();
        }
        if(args.LastError() != SO_SUCCESS)
		{
            std::cerr << "Error: " << getLastErrorText(args) << " for: " << args.OptionText() << std::endl;
			throw FastCodeMLFatalNoMsg();
        }

		switch(args.OptionId())
		{
		case OPT_VERBOSE:
			mVerboseLevel = atoi(args.OptionArg());
			break;

		case OPT_QUIET:
			mVerboseLevel = 0;
			break;

		case OPT_HELP:
			std::cerr << "Usage:" << std::endl;
			std::cerr << "    " << usage_msg << std::endl << std::endl;
			showHelp(parser_options);
			throw FastCodeMLSuccess();

		case OPT_SEED:
			mSeed = atoi(args.OptionArg());
			break;

		case OPT_BRANCH:
			mBranch = atoi(args.OptionArg());
			break;

		case OPT_IGNORE_FREQ:
			mIgnoreFreq = true;

		case OPT_EXPORT:
			mGraphFile = args.OptionArg();
			break;

		case OPT_NOT_REDUCE:
			mDoNotReduceForest = true;
			break;

		case OPT_TIMES_FROM_FILE:
			mTimesFromFile = true;
			break;

		case OPT_ONE_STEP:
			mNoMaximization = true;
			break;

		case OPT_COMP_TIMES:
			mExportComputedTimes = atoi(args.OptionArg());
			if(mExportComputedTimes > 1) mExportComputedTimes = UINT_MAX;
			break;

		case OPT_TRACE:
			mTrace = true;
			break;

		case OPT_NO_AGGRESSIVE:
			mNoAggressiveStep = true;
			break;

		case OPT_FORCE_SERIAL:
			mForceSerial = true;
			break;
		}
	}
	
	// Parse the file arguments
	switch(args.FileCount())
	{
	case 0:
		std::cerr << "Missing NEWICK TREE file" << std::endl;
		// Falltrough
	case 1:

		std::cerr << "Missing GENE file" << std::endl << std::endl;
		std::cerr << "Usage:" << std::endl;
		std::cerr << "    " << usage_msg << std::endl << std::endl;
		showHelp(parser_options);
		throw FastCodeMLFatalNoMsg();

	default:
		mTreeFile = args.File(0);
		mGeneFile = args.File(1);
		break;
	}
	
	// Some final checks and sets
	if(!mGraphFile) mExportComputedTimes = UINT_MAX;
	if(mDoNotReduceForest) mNoAggressiveStep = true;
#ifdef _OPENMP
	if(omp_get_max_threads() < 2) mForceSerial = true;
	if(mForceSerial) omp_set_num_threads(1);
#else
	mForceSerial = true;
#endif
}

