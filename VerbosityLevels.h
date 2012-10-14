
#ifndef VERBOSITYLEVELS_H
#define VERBOSITYLEVELS_H

/// Verbosity levels as set by the command line
///
enum VerbosityLevelsEnum
{
	VERBOSE_NONE			= 0,	///< No output at all
	VERBOSE_ONLY_RESULTS	= 1,	///< Output only results
	VERBOSE_INFO_OUTPUT		= 2,	///< Standard info output
	VERBOSE_MORE_DEBUG		= 3,	///< Mild debug messages
	VERBOSE_DSTRUCT_DUMP	= 4		///< Dump internal data structures
};

/// Decode the verbosity level.
///
/// @param[in] aLevel The verbosity level
///
/// @return The decoded level as string
///
extern const char* decodeVerboseLevel(unsigned int aLevel);

#endif

