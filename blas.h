
#ifndef BLAS_H
#define BLAS_H

#ifdef __cplusplus
extern "C" { 
#endif
///  DGEMV  performs one of the matrix-vector operations
///
///     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
///
///  where alpha and beta are scalars, x and y are vectors and A is an
///  m by n matrix.
///
///  @param[in] trans On entry should be "T" to transpose the "C" input A matrix
///
///  @param[in] m On entry, m specifies the number of rows of the matrix a. m must be at least zero.
///
///  @param[in] n On entry, n specifies the number of columns of the matrix a. n must be at least zero.
///
///  @param[in] alpha On entry, ALPHA specifies the scalar alpha.
///
///  @param[in] a Before entry, the leading m by n part of the array a must contain the matrix of coefficients.
///
///  @param[in] lda  On entry, lda specifies the first dimension of a as declared
///           in the calling (sub) program. lda must be at least max( 1, m ).
///
///  @param[in] x Before entry, the incremented array X must contain the vector x.
///
///  @param[in] incx On entry, incx specifies the increment for the elements of x. incx must not be zero.
///
///  @param[in] beta On entry, beta specifies the scalar beta. When beta is
///           supplied as zero then Y need not be set on input.
///
///  @param[in,out] y Before entry with beta non-zero, the incremented array y
///           must contain the vector y. On exit, y is overwritten by the
///           updated vector y.
///
///  @param[in] incy On entry, incy specifies the increment for the elements of y. incy must not be zero.
///
void dgemv_(const char *trans,
				const int *m,
				const int *n,
				const double *alpha,
				const double *a,
				const int *lda,
				const double *x,
				const int *incx,
				const double *beta,
				double *y,
				const int *incy);

//void dgemm_(const char *transa,
//				const char *transb,
//				const int *m,
//				const int *n,
//				const int *k,
//				const double *alpha,
//				double *a,
//				const int *lda,
//				double *b,
//				const int *ldb,
//				const double *beta,
//				double *c,
//				const int *ldc);

/// Forms the dot product of two vectors.
///
/// @param[in] n Vector length
/// @param[in] dx The first vector
/// @param[in] incx The stride of the vector dx
/// @param[in] dy The second vector
/// @param[in] incy The stride of the vector dy
///
/// @return The dot product of dx and dy
///
double ddot_(const int *n,
				const double *dx,
				const int *incx,
				const double *dy,
				const int *incy);

//double dnrm2_(const int *n,
//				const double *dx,
//				const int *incx);
//
//void dscal_(const int *n,
//				const double *alpha,
//				double *dx,
//				const int *incx);

#ifdef __cplusplus
}
#endif

/// Constants used by blas and lapack routines.
///
const	int		I0 = 0;		///< Integer zero.
const	int		I1 = 1;		///< Integer one.
const	double	D0 = 0.;	///< Double zero.
const	double	D1 = 1.;	///< Double one.

#endif

