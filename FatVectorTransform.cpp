
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <set>
#include "FatVectorTransform.h"
#include "MatrixSize.h"
#include "Exceptions.h"

void FatVectorTransform::setBranchDependencies(const std::vector< std::vector<ForestNode*> >& aNodesByLevel)
{
	// Push only the branch id's (and compute num branches). The root is not pushed!
	mNumBranches = 0;
	mBranchByLevel.clear();
	std::vector< std::vector<ForestNode*> >::const_reverse_iterator inbl;
	for(inbl=aNodesByLevel.rbegin(); inbl != aNodesByLevel.rend(); ++inbl)
	{
		std::vector<unsigned int> v;
		std::vector<ForestNode*>::const_iterator ifn;
		for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
		{
			v.push_back((*ifn)->mBranchId);
			++mNumBranches;
		}
		mBranchByLevel.push_back(v);
	}

	// Mark the first branch for nodes at the level below
	mFirstForLevel.assign(mNumBranches, false);
    ForestNode* curr_node = 0;
	for(inbl=aNodesByLevel.rbegin(); inbl != aNodesByLevel.rend(); ++inbl)
    {
		std::vector<ForestNode*>::const_iterator ifn;
		for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
		{
			ForestNode* parent_node = (*ifn)->mParent;

			// If this is the first visit to the parent copy the result, otherwise do a element by element multiplication
			if(parent_node != curr_node)
			{
				curr_node = parent_node;
				mFirstForLevel[(*ifn)->mBranchId] = true;
			}
		}
	}

	// Discover the parent of each branch
	mParentNode.resize(mNumBranches);
	for(inbl=aNodesByLevel.rbegin(); inbl != aNodesByLevel.rend(); ++inbl)
	{
		std::vector<ForestNode*>::const_iterator ifn;
		for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
		{
			mParentNode[(*ifn)->mBranchId] = (*ifn)->mParent->mBranchId+1; // The parent node is the node from which the branch originate
		}
	}
}


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
			char msg[128];
			sprintf(msg, "No SITE_EXISTS in mNodePresent at branch: %u", b);
			throw FastCodeMLFatal(msg);
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

		std::cerr << std::setw(2) << b << ": " << std::setw(4) << begin_idx << '-' << std::setw(4) << end_idx-1 << " (" << cnt << ")" << std::endl;
	}
}


void FatVectorTransform::printBranchVisitSequence(void) const
{
	std::cerr << std::endl << "Branch at level" << std::endl;
	unsigned int level = 1;
	std::vector< std::vector<unsigned int> >::const_iterator inbl;
	for(inbl=mBranchByLevel.begin(); inbl != mBranchByLevel.end(); ++inbl, ++level)
	{
		std::cerr << level << ": ";

		std::vector<unsigned int>::const_iterator ifn;
		for(ifn=inbl->begin(); ifn != inbl->end(); ++ifn)
		{
			std::cerr << (*ifn) << ' ';
		}

		std::cerr << std::endl;
	}

	std::cerr << std::endl << "Parent node for branch" << std::endl;
	for(unsigned int i=0; i < mNumBranches; ++i)
	{
		std::cerr << std::setw(2) << i << " -> " << std::setw(2) << mParentNode[i] << std::endl;
	}
}



void FatVectorTransform::printNodeStatus(void) const
{
	std::cerr << std::endl;
	for(unsigned int b=0; b < mNumBranches; ++b)
	{
		std::cerr << "Branch " << b << std::endl;
		bool is_num = false;
		for(unsigned int k = 0; k < mNumSites; ++k)
		{
			int x = mNodeStatus[b*mNumSites+k];
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
			char msg[128];
			sprintf(msg, "No SITE_EXISTS in mNodePresent at branch: %u", b);
			throw FastCodeMLFatal(msg);
		}

		// Compute the last valid site (actually it points one after)
		size_t end_idx = mNumSites;
		for(; end_idx > begin_idx; --end_idx)
		{
			if(mNodeStatus[b*mNumSites+end_idx-1] == FatVectorTransform::SITE_EXISTS) break;
		}

		// Get the compaction moves
		VectorOfRanges cmds;
		for(int site_to=end_idx-1; site_to >= (int)begin_idx; --site_to)
		{
			// Select the first hole (from right)
			if(mNodeStatus[b*mNumSites+site_to] == FatVectorTransform::SITE_EXISTS) continue;

			// From left find the first valid entry
			unsigned int site_from = begin_idx;

			// Save the move command
			cmds.push_back(Range(site_from, site_to));

			// Update the left limit
			for(++begin_idx; begin_idx < (unsigned int)site_to; ++begin_idx)
			{
				// Select the first valid site (from left)
				if(mNodeStatus[b*mNumSites+begin_idx] == FatVectorTransform::SITE_EXISTS) break;
			}
		}
		mCopyCmds.push_back(cmds);

		// Save the new start index and count
		mLimits[b] = std::make_pair(begin_idx, end_idx-begin_idx);

		// Compute the reuse of another value moves
		VectorOfRangesNoCnt reuse;
		for(unsigned int k=0; k < mNumSites; ++k)
		{
			// Select a reuse pointer
			int x = mNodeStatus[b*mNumSites+k];
			if(x >= FatVectorTransform::SITE_FIRST_NUM)
			{
				reuse.push_back(RangeNoCnt(x, k));
			}
		}
		mReuseCmds.push_back(reuse);
	}

	// Remove the node status array no more needed
	mNodeStatus.clear();

	// Try to combine contiguous ranges
	for(b=0; b < mNumBranches; ++b)
	{
		unsigned int nc = mCopyCmds[b].size();
		if(nc < 2) continue;

		// Start with two valid
		for(unsigned int i=0; i < nc-1;)
		{
			if(mCopyCmds[b][i].from+1 == mCopyCmds[b][i+1].from && mCopyCmds[b][i].to == mCopyCmds[b][i+1].to+1)
			{
				// Try to extend the range to other with the same ordering
				unsigned int j=i+1;
				for(; j < nc-1; ++j)
				{
					if(mCopyCmds[b][j].from+1 != mCopyCmds[b][j+1].from || mCopyCmds[b][j].to != mCopyCmds[b][j+1].to+1) break;
				}

				// Update the command list
				// Example: (100, 10, 1) (101, 9, 1) --> (100, 9, 2) (101, 9, 0)
				mCopyCmds[b][i].cnt = j-i+1;
				mCopyCmds[b][i].to  = mCopyCmds[b][j].to;
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
		std::cerr << std::endl << "*** Branch " << b << std::endl;

		VectorOfRanges::const_iterator icc;
		for(icc=mCopyCmds[b].begin(); icc != mCopyCmds[b].end(); ++icc)
		{
			if(icc->cnt == 1)
				std::cerr << "C " << std::setw(4) << icc->from << " - " << std::setw(4) << icc->to << std::endl;
			else if(icc->cnt > 1)
				std::cerr << "C " << std::setw(4) << icc->from << " - " << std::setw(4) << icc->to << " (" << icc->cnt << ")" << std::endl;
		}

		VectorOfRangesNoCnt::const_iterator icr;
		for(icr=mReuseCmds[b].begin(); icr != mReuseCmds[b].end(); ++icr)
		{
			std::cerr << "R " << std::setw(4) << icr->from << " - " << std::setw(4) << icr->to << std::endl;
		}

		std::cerr << "L   from: " << mLimits[b].first << " cnt: " << mLimits[b].second << std::endl;
	}
}


void FatVectorTransform::preCompactLeaves(CacheAlignedDoubleVector& aProbs)
{
	// If the forest has not been reduced do nothing (this should not happens)
	if(mNoTransformations) return;

	// Find all the nodes that are parent of other nodes (ie. they are not leaves)
	std::set<unsigned int> non_leaves;
	non_leaves.insert(mParentNode.begin(), mParentNode.end());

	// Make a list of leaf nodes
	std::vector<unsigned int> leaves;
	for(unsigned int node=1; node <= mNumBranches; ++node)
	{
		// Check if the node is a leaf, if not skip it
		if(non_leaves.find(node) == non_leaves.end()) leaves.push_back(node);
	}

	// For all leaves and all sets
	int len = leaves.size()*Nt;
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(len, leaves, aProbs)
#else
	#pragma omp parallel for default(shared)
#endif
	for(int i=0; i < len; ++i)
	{
		unsigned int node    = leaves[i / Nt];
		unsigned int set_idx = i % Nt;
		unsigned int start   = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites);

		// Do all the copies as requested
		VectorOfRanges::const_iterator icc;
		for(icc=mCopyCmds[node-1].begin(); icc != mCopyCmds[node-1].end(); ++icc)
		{
			if(icc->cnt == 1)
			{
				memcpy(&aProbs[start+VECTOR_SLOT*icc->to], &aProbs[start+VECTOR_SLOT*icc->from], N*sizeof(double));
			}
			else if(icc->cnt > 1)
			{
				memcpy(&aProbs[start+VECTOR_SLOT*icc->to], &aProbs[start+VECTOR_SLOT*icc->from], (VECTOR_SLOT*icc->cnt-(VECTOR_SLOT-N))*sizeof(double));
			}
		}
	}
}


void FatVectorTransform::postCompact(CacheAlignedDoubleVector& aStepResults, CacheAlignedDoubleVector& aProbs, unsigned int aLevel, unsigned int aNumSets)
{
	int nsns = VECTOR_SLOT*mNumSites*aNumSets;
	if(mNoTransformations)
	{
		std::vector<unsigned int>::const_iterator ibl;
		for(ibl=mBranchByLevel[aLevel].begin(); ibl != mBranchByLevel[aLevel].end(); ++ibl)
		{
			unsigned int   my_branch = *ibl;
			unsigned int parent_node = mParentNode[my_branch];
			unsigned int     my_node = my_branch + 1;

			if(mFirstForLevel[my_branch])
			{
				memcpy(&aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node],
					   &aStepResults[VECTOR_SLOT*mNumSites*Nt*my_node],
					   (VECTOR_SLOT*mNumSites*aNumSets-(VECTOR_SLOT-N))*sizeof(double));
			}
			else
			{
#ifdef _MSC_VER
				#pragma omp parallel for default(none) shared(parent_node, my_node, aNumSets, aProbs, aStepResults, nsns)
#else
				#pragma omp parallel for default(shared)
#endif
                for(int i=0; i < nsns; ++i)
                {
                    aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node+i] *= aStepResults[VECTOR_SLOT*mNumSites*Nt*my_node+i];
                }
			}
		}
	}
	else
	{
		// For all the branches just processed
		std::vector<unsigned int>::const_iterator ibl;
		for(ibl=mBranchByLevel[aLevel].begin(); ibl != mBranchByLevel[aLevel].end(); ++ibl)
		{
			unsigned int branch      = *ibl;
			unsigned int node        = branch + 1;
			unsigned int parent_node = mParentNode[branch];

			// Reverse all copies (copy back the values copied in the previous step to fill holes)
			VectorOfRanges::const_iterator icc;
			for(icc=mCopyCmds[branch].begin(); icc != mCopyCmds[branch].end(); ++icc)
			{
				if(icc->cnt == 1)
				{
					for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
					{
						unsigned int from_idx = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites+icc->from);
						unsigned int to_idx   = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites+icc->to);

						memcpy(&aStepResults[from_idx],
							   &aStepResults[to_idx],
							   N*sizeof(double));
					}
				}
				else if(icc->cnt > 1)
				{
					for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
					{
						unsigned int from_idx = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites+icc->from);
						unsigned int to_idx   = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites+icc->to);

						memcpy(&aStepResults[from_idx],
							   &aStepResults[to_idx],
							   (VECTOR_SLOT*icc->cnt-(VECTOR_SLOT-N))*sizeof(double));
					}
				}
			}

			// Reuse values 
			VectorOfRangesNoCnt::const_iterator icr;
			for(icr=mReuseCmds[branch].begin(); icr != mReuseCmds[branch].end(); ++icr)
			{
				for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
				{
					unsigned int from_idx = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites+icr->from);
					unsigned int to_idx   = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites+icr->to);

					memcpy(&aStepResults[to_idx],
						   &aStepResults[from_idx],
						   N*sizeof(double));
				}
			}

			if(mFirstForLevel[branch])
			{
				memcpy(&aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node], &aStepResults[VECTOR_SLOT*mNumSites*Nt*node], VECTOR_SLOT*mNumSites*aNumSets*sizeof(double));
			}
			else
			{
#ifdef _MSC_VER
				#pragma omp parallel for default(none) shared(parent_node, node, aNumSets, aProbs, aStepResults, nsns)
#else
				#pragma omp parallel for default(shared)
#endif
                for(int i=0; i < nsns; ++i)
                {
                    aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node+i] *= aStepResults[VECTOR_SLOT*mNumSites*Nt*node+i];
                }
			}

			// Copy for the next branch (if this branch does not lead to the root)
			if(parent_node)
			{
				// Do all the copies as requested
				for(icc=mCopyCmds[parent_node-1].begin(); icc != mCopyCmds[parent_node-1].end(); ++icc)
				{
					if(icc->cnt == 1)
					{
						for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
						{
							unsigned int from_idx = VECTOR_SLOT*(mNumSites*Nt*parent_node+set_idx*mNumSites+icc->from);
							unsigned int to_idx   = VECTOR_SLOT*(mNumSites*Nt*parent_node+set_idx*mNumSites+icc->to);

							memcpy(&aProbs[to_idx],
								   &aProbs[from_idx],
								   N*sizeof(double));
						}
					}
					if(icc->cnt > 1)
					{
						for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
						{
							unsigned int from_idx = VECTOR_SLOT*(mNumSites*Nt*parent_node+set_idx*mNumSites+icc->from);
							unsigned int to_idx   = VECTOR_SLOT*(mNumSites*Nt*parent_node+set_idx*mNumSites+icc->to);

							memcpy(&aProbs[to_idx],
								   &aProbs[from_idx],
								   (VECTOR_SLOT*icc->cnt-(VECTOR_SLOT-N))*sizeof(double));
						}
					}
				}
			}
		}
	}
}
