
#ifndef HIGHLEVELCOORDINATOR_H
#define HIGHLEVELCOORDINATOR_H

#include <vector>
#include "Forest.h"

class HighLevelCoordinator
{
public:
	HighLevelCoordinator(int* aArgc, char*** aArgv);
	~HighLevelCoordinator();

	bool startWork(Forest& aForest, unsigned int aSeed, unsigned int aVerbose=0, bool aNoMaximization=false, bool aTimesFromFile=true, unsigned int aOptimizationAlgo=0);
	bool isMaster(void) const {return mRank == 0;}
	int  numJobs(void) const {return mSize;}


private:
	void doMaster(void);
	void doWorker(Forest& aForest, unsigned int aSeed, bool aNoMaximization, bool aTimesFromFile, unsigned int aOptimizationAlgo);


private:
	unsigned int		mVerbose;
	int					mRank;
	int					mSize;
	unsigned int		mNumBranches;

	struct WorkTable;
	WorkTable*			mWorkTable;
};




#endif
