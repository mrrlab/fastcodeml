// The following headers are required for all allocators.
#include <stddef.h>  // Required for size_t and ptrdiff_t and NULL
#include <stdexcept> // Required for std::length_error

// The following headers contain stuff that AlignedAllocator uses.
#include <stdlib.h>  // For malloc() and free()

#include "AlignedMalloc.h"
#include "AlignedAllocator.h"

#ifndef _MSC_VER
#include <stdint.h>  // for uintptr_t
#endif

// Alignment must be power of 2 (1,2,4,8,16...)
void* alignedMalloc(size_t size, size_t alignment)
{
    uintptr_t r = reinterpret_cast<uintptr_t>(malloc(size + --alignment + sizeof(uintptr_t)));
    uintptr_t t = r + sizeof(uintptr_t);
    uintptr_t o = (t + alignment) & ~static_cast<uintptr_t>(alignment);
    if(!r) return NULL;
    reinterpret_cast<uintptr_t*>(o)[-1] = r;
    return reinterpret_cast<void*>(o);
}

void alignedFree(void* p)
{
    if(!p) return;
    free(reinterpret_cast<void*>(reinterpret_cast<uintptr_t*>(p)[-1]));
}


#if 0

// The following headers contain stuff that main() uses.
#include <iostream>  // For std::cerr
#include <ostream>   // For std::endl
#include <vector>    // For std::vector

int main()
{
    using namespace std;

    cerr << "Constructing l:" << endl;

    vector<double, AlignedAllocator<double, 8> > l;
	l.reserve(10);
    cerr << endl << "l.push_back(1729):" << endl;

    l.push_back(1729.);

    cerr << endl << "l.push_back(2161):" << endl;

    l.push_back(2161.);

    cerr << endl;
	double* p = &l[0];
	int x = reinterpret_cast<int>(p);
	cerr << "Aligned on 16: " << x%16 << endl;
	cerr << "Aligned on 8:  " << x%8 << endl;
	cerr << "Aligned on 4:  " << x%4 << endl;
    cerr << endl;

    for (vector<double, AlignedAllocator<double, 8> >::const_iterator i = l.begin(); i != l.end(); ++i) {
        cerr << "Element: " << *i << endl;
    }

    cerr << endl << "Destroying l:" << endl;
}
#endif

