
#ifndef FASTCODEMLEXCEPTIONS_H
#define FASTCODEMLEXCEPTIONS_H

#include <stdexcept>
#include <string>
#include <sstream>

/// Fatal error in FastCodeML.
/// The message explains the reason
///
class FastCodeMLFatal : public std::runtime_error
{
public:
	/// Constructor.
	/// No message because it has already been printed.
	///
	FastCodeMLFatal(void) : runtime_error("")
	{}

	/// Constructor.
	///
	/// @param[in] aMessage The message to be printed before termination
	///
	FastCodeMLFatal(const char *aMessage) : runtime_error(aMessage)
	{}

	/// Constructor.
	///
	/// @param[in] aMessage The message to be printed before termination (has been formatted printing into a std::ostringstream).
	///
	FastCodeMLFatal(const std::ostringstream& aMessage) : runtime_error(aMessage.str().c_str())
	{}
};


/// Early successful termination exception.
/// It does not signal a fatal error, but something like having printed the help.
///
class FastCodeMLSuccess : public std::exception
{
public:
	/// Constructor.
	///
	FastCodeMLSuccess() : exception()
	{}
};


#endif

