
#include <iostream>
#include <iomanip>
//#include <algorithm>
#include <cstring>
#include "FatVectorTransform.h"
#include "MatrixSize.h"

void FatVectorTransform::printCountGoodElements(void) const
{
	std::cerr << std::endl;
	for(unsigned int b=0; b < mNumBranches; ++b)
	{
		size_t begin_idx = 0;
		for(; begin_idx < mNumSites; ++begin_idx)
		{
			int x = mNodeStatus[b*mNumSites+begin_idx];
			if(x == FatVectorTransform::SITE_EXISTS) break;
		}
		if(begin_idx == mNumSites)
		{
			std::cerr << "No SITE_EXISTS in mNodePresent at branch: " << b+1 << std::endl;
			return;
		}

		size_t end_idx = mNumSites;
		for(; end_idx > begin_idx; --end_idx)
		{
			int x = mNodeStatus[b*mNumSites+end_idx-1];
			if(x == FatVectorTransform::SITE_EXISTS) break;
		}

		// Count the good elements
		unsigned int cnt = 0;
		for(unsigned int k=begin_idx; k < end_idx; ++k) if(mNodeStatus[b*mNumSites+k] == FatVectorTransform::SITE_EXISTS) ++cnt;

		std::cerr << std::setw(2) << b+1 << ": " << std::setw(4) << begin_idx << '-' << std::setw(4) << end_idx-1 << " (" << cnt << ")" << std::endl;
	}
}


void FatVectorTransform::printBranchVisitSequence(const std::vector< std::vector<unsigned int> >& aBranchByLevel) const
{
	std::cerr << std::endl;
	unsigned int level = 1;
	std::vector< std::vector<unsigned int> >::const_iterator inbl;
	for(inbl=aBranchByLevel.begin(); inbl != aBranchByLevel.end(); ++inbl, ++level)
	{
		std::cerr << level << ": ";

		std::vector<unsigned int>::const_iterator ifn;
		for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
		{
			std::cerr << (*ifn)+1 << ' ';
		}

		std::cerr << std::endl;
	}
}


void FatVectorTransform::printNodeStatus(void) const
{
	std::cerr << std::endl;
	for(size_t j=0; j < mNumBranches; ++j)
	{
		std::cerr << "Branch " << j+1 << std::endl;
		bool is_num = false;
		for(size_t k = 0; k < mNumSites; ++k)
		{
			int x = mNodeStatus[j*mNumSites+k];
			if(x == FatVectorTransform::SITE_NOT_EXISTS)  {std::cerr << '-'; is_num = false;}
			else if(x == FatVectorTransform::SITE_EXISTS) {std::cerr << 'x'; is_num = false;}
			else                                          {if(is_num) std::cerr << ' '; std::cerr << x; is_num = true;}
		}
		std::cerr << std::endl << std::endl;
	}
}


void FatVectorTransform::compactMatrix(void)
{
	// For each branch
	unsigned int b;
	for(b=0; b < mNumBranches; ++b)
	{
		// Compute the index of the first valid site
		size_t begin_idx = 0;
		for(; begin_idx < mNumSites; ++begin_idx)
		{
			if(mNodeStatus[b*mNumSites+begin_idx] == FatVectorTransform::SITE_EXISTS) break;
		}
		if(begin_idx == mNumSites)
		{
			std::cerr << "No SITE_EXISTS in mNodePresent at branch: " << b+1 << std::endl;
			return;
		}

		// Compute the last valid site (actually it points one after)
		size_t end_idx   = mNumSites;
		for(; end_idx > begin_idx; --end_idx)
		{
			if(mNodeStatus[b*mNumSites+end_idx-1] == FatVectorTransform::SITE_EXISTS) break;
		}

		// Get the compaction moves
		VectorOfRanges cmds;
		for(unsigned int site_to=end_idx-1; site_to >= begin_idx; --site_to)
		{
			// Select the first hole (from right)
			if(mNodeStatus[b*mNumSites+site_to] == FatVectorTransform::SITE_EXISTS) continue;

			// From left find the first valid entry
			unsigned int site_from = begin_idx;

			// Save the move command
			cmds.push_back(Range(site_from, site_to));

			// Update the left limit
			for(++begin_idx; begin_idx < site_to; ++begin_idx)
			{
				// Select the first valid site (from left)
				if(mNodeStatus[b*mNumSites+begin_idx] == FatVectorTransform::SITE_EXISTS) break;
			}
		}
		mCopyCmds.push_back(cmds);

		// Save the new limits
		mLimits[b] = std::make_pair(begin_idx, end_idx-1);

		// Compute the reuse of another value moves
		cmds.clear();
		for(unsigned int k=0; k < mNumSites; ++k)
		{
			// Select a reuse pointer
			int x = mNodeStatus[b*mNumSites+k];
			if(x >= FatVectorTransform::SITE_FIRST_NUM)
			{
				cmds.push_back(Range(x, k));
			}
		}
		mReuseCmds.push_back(cmds);
	}

	// Remove the node status array no more needed
	mNodeStatus.clear();

	// Try to combine contiguous ranges
	for(b=0; b < mNumBranches; ++b)
	{
		unsigned int nc = mCopyCmds[b].size();

		// Start with two valid
		for(unsigned int i=0; i < nc-1;)
		{
			if(mCopyCmds[b][i].from+1 == mCopyCmds[b][i+1].from && mCopyCmds[b][i].to == mCopyCmds[b][i+1].to+1)
			{
				unsigned int j=i+1;
				for(; j < nc-1; ++j)
				{
					if(mCopyCmds[b][j].from+1 != mCopyCmds[b][j+1].from || mCopyCmds[b][j].to != mCopyCmds[b][j+1].to+1) break;
				}

				// Update the command list
				mCopyCmds[b][i].cnt = j-i+1;
				for(unsigned int k=i+1; k <= j; ++k) mCopyCmds[b][k].cnt = 0;

				i = j+1;
			}
			else
			{
				++i;
			}
		}
	}
}

void FatVectorTransform::printCommands(void) const
{
	for(unsigned int b=0; b < mNumBranches; ++b)
	{
		std::cerr << std::endl << "*** Branch " << b+1 << std::endl;

		VectorOfRanges::const_iterator icc;
		for(icc=mCopyCmds[b].begin(); icc != mCopyCmds[b].end(); ++icc)
		{
			if(icc->cnt > 0)
				std::cerr << "C " << std::setw(4) << icc->from << '-' << std::setw(4) << icc->to << " (" << icc->cnt << ")" << std::endl;
		}

		for(icc=mReuseCmds[b].begin(); icc != mReuseCmds[b].end(); ++icc)
		{
			if(icc->cnt == 1)
				std::cerr << "R " << std::setw(4) << icc->from << '-' << std::setw(4) << icc->to << std::endl;
			else
				std::cerr << "R " << std::setw(4) << icc->from << '-' << std::setw(4) << icc->to << " (" << icc->cnt << ")" << std::endl;
		}

		std::cerr << "L " << std::setw(4) << mLimits[b].first << '-' << std::setw(4) << mLimits[b].second << std::endl;
	}
}


void FatVectorTransform::preCompact(std::vector<double>& aProbs, unsigned int aBranch, unsigned int aSet)
{
	VectorOfRanges::const_iterator icc;
	for(icc=mCopyCmds[aBranch].begin(); icc != mCopyCmds[aBranch].end(); ++icc)
	{
		if(icc->cnt > 0)
		{
			memcpy(&aProbs[N*mNumSites*Nt*aBranch+N*mNumSites*aSet+N*icc->to],
				   &aProbs[N*mNumSites*Nt*aBranch+N*mNumSites*aSet+N*icc->from],
				   N*icc->cnt*sizeof(double));
		}
	}
}
