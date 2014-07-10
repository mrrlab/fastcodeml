#include "Job.h"

#include <iomanip>

Job::Job(
      int aForestIndex, int aBranch, JobType aJobType, int aSizeOtherData)
    : mForestIndex(aForestIndex), mBranch(aBranch), mJobType(aJobType),
      mSizeOtherData(aSizeOtherData)
{
    // Nothing to do.
}

Job& Job::operator=(const Job &that)
{
    if (this == &that)
    {
        return (*this);
    }
    this->mForestIndex    = that.mForestIndex;
    this->mBranch         = that.mBranch;
    this->mJobType        = that.mJobType;
    this->mSizeOtherData  = that.mSizeOtherData;

    return(*this);
}

Job::Job(const Job &that)
    : mForestIndex(that.mForestIndex), mBranch(that.mBranch), mJobType(that.mJobType),
      mSizeOtherData(that.mSizeOtherData)
{
    // Nothing to do.
}

std::ostream&
operator<<(std::ostream &os, const Job &aJob)
{
    os << std::setw(20) << "Forest index: "    << aJob.mForestIndex   << ", "
       << std::setw(20) << "Branch index: "    << aJob.mBranch        << ", "
       << std::setw(20) << "Job type: "        << aJob.mJobType       << ", "
       << std::setw(20) << "Size other data: " << aJob.mSizeOtherData;
    return(os);
}
