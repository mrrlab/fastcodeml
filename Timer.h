
#ifndef TIMER_H
#define TIMER_H

#include <time.h>

#ifdef WIN32

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
	Timer() : delta(0.) {}

	/// Start the timer
	///
	void start(void)
	{
		time(&start_time);
	}

	/// Stop the timer
	/// @return The elapsed time in milliseconds
	///
	double stop(void)
	{
		time_t end_time;
		time(&end_time);

		delta = difftime(end_time, start_time)*1000;
		return delta;
	}

	/// Return the elapsed time (after a start/stop cycle)
	/// @return The elapsed time in milliseconds
	///
	double get(void) const
	{
		return delta;
	}

private:
	time_t start_time;			///< The start time
	double delta;				///< The elapsed time in milliseconds
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
	Timer() : delta(0.) {}

	/// Start the timer
	///
	void start(void)
	{
		gettimeofday(&start_time, NULL);
	}

	/// Stop the timer
	/// @return The elapsed time in milliseconds
	///
	double stop(void)
	{
		struct timeval end_time;
		gettimeofday(&end_time, NULL);

		//delta = (double) end_time.tv_sec - start_time.tv_sec;
		//delta += (end_time.tv_usec - start_time.tv_usec) * 1e-6;

		delta = (double)end_time.tv_sec*1e6+end_time.tv_usec;
		delta -= (double)start_time.tv_sec*1e6+start_time.tv_usec;
		delta /= 1000;

		return delta;
	}

	/// Return the elapsed time (after a start/stop cycle)
	/// @return The elapsed time in milliseconds
	///
	double get(void) const
	{
		return delta;
	}

private:
	struct timeval start_time;	///< The start time
	double delta;				///< The elapsed time in milliseconds
};

#endif
#endif



