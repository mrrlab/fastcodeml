
#ifndef ALIGNEDMALLOC_H
#define ALIGNEDMALLOC_H

// Alignment must be power of 2 (1,2,4,8,16...)
extern void* alignedMalloc(size_t size, size_t alignment);
extern void alignedFree(void* p);



#endif
