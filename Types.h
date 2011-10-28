
#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include "AlignedAllocator.h"
#include "MatrixSize.h"

typedef std::vector<double, AlignedAllocator<double, CACHE_LINE_ALIGN> > CacheAlignedDoubleVector;
//typedef std::vector<double> CacheAlignedDoubleVector;

#endif

