
#ifndef MATRIXSIZE_H
#define MATRIXSIZE_H
/// The matrix side (64 possible three-letters codons less 3 STOP codons)
///
static const int N = 61;

/// Number of codon classes (4)
///
//static const int Nc = 4;

/// Number of possible codons (64)
///
static const int N64 = 64;

/// Max number of possible tree traversals for each cycle (4)
///
static const int Nt = 4;

/// Slot size for a matrix. It should be equal or larger than N*N. Filler based on measurement.
static const int MATRIX_SLOT = 61*61+3;

/// Slot size for a vector. It should be equal or larger than N. Size based on measurements.
#ifdef NEW_LIKELIHOOD
static const int VECTOR_SLOT = 61;
#else
static const int VECTOR_SLOT = 66;
#endif

/// Alignment to avoid cache line false sharing
static const int CACHE_LINE_ALIGN = 64;

#endif

