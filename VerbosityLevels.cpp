
#include "VerbosityLevels.h"

const char* decodeVerboseLevel(int aLevel)
{
	switch(aLevel)
	{
	case VERBOSE_NONE:			return "No output at all";
	case VERBOSE_ONLY_RESULTS:	return "Output only results";
	case VERBOSE_INFO_OUTPUT:	return "Standard info output";
	case VERBOSE_MORE_DEBUG:	return "Mild debug messages";
	case VERBOSE_DSTRUCT_DUMP:	return "Dump internal structures";
	default:					return "Which level is this?";
	}
}

