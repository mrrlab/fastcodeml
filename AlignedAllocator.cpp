// The following headers are required for all allocators.
#include <stddef.h>  // Required for size_t and ptrdiff_t and NULL
//#include <new>       // Required for placement new and std::bad_alloc
#include <stdexcept> // Required for std::length_error

// The following headers contain stuff that AlignedAllocator uses.
#include <stdlib.h>  // For malloc() and free()
//#include <iostream>  // For std::cout
//#include <ostream>   // For std::endl

// The following headers contain stuff that main() uses.
//#include <vector>      // For std::vector

#include "AlignedMalloc.h"
#include "AlignedAllocator.h"

#ifndef _MSC_VER
#include <stdint.h>  // for uintptr_t
#endif

// Alignment must be power of 2 (1,2,4,8,16...)
void* alignedMalloc(size_t size, size_t alignment)
{
    uintptr_t r = (uintptr_t)malloc(size + --alignment + sizeof(uintptr_t));
    uintptr_t t = r + sizeof(uintptr_t);
    uintptr_t o = (t + alignment) & ~(uintptr_t)alignment;
    if(!r) return NULL;
    ((uintptr_t*)o)[-1] = r;
    return (void*)o;
}

void alignedFree(void* p)
{
    if(!p) return;
    free((void*)(((uintptr_t*)p)[-1]));
}


#if 0
int main()
{
    using namespace std;

    cout << "Constructing l:" << endl;

    vector<double, AlignedAllocator<double, 8> > l;
	l.reserve(10);
    cout << endl << "l.push_back(1729):" << endl;

    l.push_back(1729.);

    cout << endl << "l.push_back(2161):" << endl;

    l.push_back(2161.);

    cout << endl;
	double* p = &l[0];
	int x = (int)p;
	cout << "Aligned on 16: " << x%16 << endl;
	cout << "Aligned on 8:  " << x%8 << endl;
	cout << "Aligned on 4:  " << x%4 << endl;
    cout << endl;

    for (vector<double, AlignedAllocator<double, 8> >::const_iterator i = l.begin(); i != l.end(); ++i) {
        cout << "Element: " << *i << endl;
    }

    cout << endl << "Destroying l:" << endl;
}
#endif

