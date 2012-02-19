
#ifndef TIMER_H
#define TIMER_H

#ifdef _MSC_VER

#include <time.h>

/// Simple timer.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-08-31 (initial version)
///     @version 1.0
///
class Timer
{
public:
	/// Constructor
	///
	Timer() : mDelta(0.) {}

	/// Start the timer
	///
	void start(void)
	{
		time(&mStartTime);
	}

	/// Stop the timer
	/// @return The elapsed time in milliseconds
	///
	double stop(void)
	{
		time_t end_time;
		time(&end_time);

		mDelta = difftime(end_time, mStartTime)*1000;
		return mDelta;
	}

	/// Return the elapsed time (after a start/stop cycle)
	/// @return The elapsed time in milliseconds
	///
	double get(void) const
	{
		return mDelta;
	}

private:
	time_t mStartTime;			///< The start time
	double mDelta;				///< The elapsed time in milliseconds
};

#else

#include <sys/time.h> // gettimeofday

/// Simple timer
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-08-31 (initial version)
///     @version 1.0
///
class Timer
{
public:
	/// Constructor
	///
	Timer() : mDelta(0.) {start();}

	/// Start the timer
	///
	void start(void)
	{
		gettimeofday(&mStartTime, NULL);
	}

	/// Stop the timer
	/// @return The elapsed time in milliseconds
	///
	double stop(void)
	{
		struct timeval end_time;
		gettimeofday(&end_time, NULL);

		//mDelta = (double) end_time.tv_sec - mStartTime.tv_sec;
		//mDelta += (end_time.tv_usec - mStartTime.tv_usec) * 1e-6;

		mDelta  = static_cast<double>(end_time.tv_sec)*1e6+end_time.tv_usec;
		mDelta -= static_cast<double>(mStartTime.tv_sec)*1e6+mStartTime.tv_usec;
		mDelta /= 1000.;

		return mDelta;
	}

	/// Return the elapsed time (after a start/stop cycle)
	/// @return The elapsed time in milliseconds
	///
	double get(void) const
	{
		return mDelta;
	}

private:
	struct timeval mStartTime;	///< The start time
	double           mDelta;		///< The elapsed time in milliseconds
};

#endif
#endif
