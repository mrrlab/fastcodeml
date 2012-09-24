

#ifndef WRITERESULTS_H
#define WRITERESULTS_H


/// Write all results on a given file.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-24 (initial version)
///     @version 1.0
///
///
class WriteResults
{
public:
	/// Constructor
	///
	explicit WriteResults(const char* aFilename) : mFilename(aFilename) {}

	void saveResults(const char* aFilename);

private:
	const char*	mFilename;
};

#endif

