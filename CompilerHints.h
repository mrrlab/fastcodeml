
#ifndef COMPILERHINTS_H
#define COMPILERHINTS_H


#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
//  Intel
#define ALIGN64 __declspec(align(64))
//#define PURE    __declspec(const)

#elif defined(__GNUC__)
//  GNU C++
#define ALIGN64 __attribute__ ((aligned (64)))
//#define PURE    __attribute__ ((pure))

#elif defined(_MSC_VER)
// Microsoft Visual C++
#define ALIGN64 __declspec(align(64))
//#define PURE

#elif defined(__PGI)
//  PGI C++
#define ALIGN64 __attribute__ ((aligned (64)))
//#define PURE    __attribute__ ((pure))

#else
#warning "Unknown compiler detected"
#define ALIGN64
//#define PURE

#endif


#endif

