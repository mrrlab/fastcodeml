
#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H

#include "MatrixSize.h"

#ifdef USE_LAPACK
#include "blas.h"
#endif

/// Dot product for vectors of variable length.
///
/// @param[in] aV1 First vector
/// @param[in] aV2 Second vector
/// @param[in] aCnt Vectors lenght
///
/// @return The dot product
///
inline double dot(const double* aV1, const double* aV2, unsigned int aCnt)
{
#ifdef USE_LAPACK
	int n = aCnt;
	return ddot_(&n, aV1, &I1, aV2, &I1);
#else
	double tot = 0.;
	for(unsigned int i=0; i < aCnt; ++i) tot += aV1[i]*aV2[i];
	return tot;
#endif
}


/// Dot product specialized for 61 elements vectors.
///
/// @param[in] aV1 First vector
/// @param[in] aV2 Second vector
///
/// @return The dot product
///
inline double dot(const double* aV1, const double* aV2)
{
#ifdef USE_LAPACK
	return ddot_(&N, aV1, &I1, aV2, &I1);
#else
	double tot = 0.;
	for(int i=0; i < N; ++i) tot += aV1[i]*aV2[i];
	return tot;
#endif
}


#endif

