
#ifndef FASTCODEMLEXCEPTIONS_H
#define FASTCODEMLEXCEPTIONS_H

#include <stdexcept>


/// Fatal error in FastCodeML.
/// The message explains the reason
///
class FastCodeMLFatal : public std::runtime_error
{
public:
	/// Constructor.
	///
	/// @param[in] aMessage The message to be printed before termination
	///
	FastCodeMLFatal(const char *aMessage) : runtime_error(aMessage)
	{}
};

/// Fatal error in FastCodeML.
/// The message explains the reason but it is not printed because a more detailed message has already been output.
///
class FastCodeMLFatalNoMsg : public FastCodeMLFatal
{
public:
	/// Constructor
	///
	FastCodeMLFatalNoMsg() : FastCodeMLFatal("")
	{}
};

/// Early termination message.
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

