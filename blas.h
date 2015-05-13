
#ifndef BLAS_H
#define BLAS_H

//#define USE_REFERENCE_BLAS

// These mappings are for ACML 5.0.1 64 bits
#ifdef _MSC_VER
#ifdef USE_REFERENCE_BLAS
#define ddot_   _ddot_
#define dcopy_	_dcopy_
#define dscal_  _dscal_
#define dsymv_  _dsymv_
#define dsymm_  _dsymm_
#define dsyevr_ _dsyevr_
#define dsyevd_ _dsyevd_
#define dsyrk_  _dsyrk_
#define dnrm2_  _dnrm2_
#define daxpy_  _daxpy_
#define dgemm_	_dgemm_
#else
#define ddot_   DDOT
#define dcopy_	DCOPY
#define dscal_  DSCAL
#define dsymv_  DSYMV
#define dsymm_  DSYMM
#define dsyevr_ DSYEVR
#define dsyevd_ DSYEVD
#define dsyrk_  DSYRK
#define dnrm2_  DNRM2
#define daxpy_  DAXPY
#define dgemv_	DGEMV
#define dgemm_	DGEMM
#endif
#endif

#ifdef _BGQ_COMPILER
#define dsyevr_ dsyevr
#define ddot_  ddot
#define dnrm2_ dnrm2
#define dcopy_ dcopy
#define dscal_ dscal
#define dsyrk_ dsyrk
#define dsymv_ dsymv
#define daxpy_ daxpy
#define dgemv_ dgemv
#define dgemm_ dgemm
#endif


#ifdef __cplusplus
extern "C" { 
#endif

/// DSYRK  performs one of the symmetric rank k operations
/// C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C,
/// where alpha and beta are scalars, C is an n by n symmetric matrix
/// and A is an n by k matrix in the first case and a k by n matrix
/// in the second case.
///
/// @param[in] uplo On entry, UPLO specifies whether the upper or lower
/// triangular part of the array C is to be referenced as
/// follows:
/// UPLO = 'U' or 'u' Only the upper triangular part of C
/// is to be referenced.
/// UPLO = 'L' or 'l' Only the lower triangular part of C
/// is to be referenced.
///
/// @param[in] trans On entry, TRANS specifies the operation to be performed as follows:
/// TRANS = 'N' or 'n' C := alpha*A*A**T + beta*C.
/// TRANS = 'T' or 't' C := alpha*A**T*A + beta*C.
/// TRANS = 'C' or 'c' C := alpha*A**T*A + beta*C.
///
/// @param[in] n On entry, N specifies the order of the matrix C. N must be
/// at least zero.
///
/// @param[in] k On entry with TRANS = 'N' or 'n', K specifies the number
/// of columns of the matrix A, and on entry with
/// TRANS = 'T' or 't' or 'C' or 'c', K specifies the number
/// of rows of the matrix A. K must be at least zero.
///
/// @param[in] alpha On entry, ALPHA specifies the scalar alpha.
///
/// @param[in] a Array of DIMENSION ( LDA, ka ), where ka is
/// k when TRANS = 'N' or 'n', and is n otherwise.
/// Before entry with TRANS = 'N' or 'n', the leading n by k
/// part of the array A must contain the matrix A, otherwise
/// the leading k by n part of the array A must contain the
/// matrix A.
///
/// @param[in] lda On entry, LDA specifies the first dimension of A as declared
/// in the calling (sub) program. When TRANS = 'N' or 'n'
/// then LDA must be at least max( 1, n ), otherwise LDA must be at least max( 1, k ).
///
/// @param[in] beta On entry, BETA specifies the scalar beta.
///
/// @param[in,out] c Array of DIMENSION ( LDC, n ).
/// Before entry with UPLO = 'U' or 'u', the leading n by n
/// upper triangular part of the array C must contain the upper
/// triangular part of the symmetric matrix and the strictly
/// lower triangular part of C is not referenced. On exit, the
/// upper triangular part of the array C is overwritten by the
/// upper triangular part of the updated matrix.
/// Before entry with UPLO = 'L' or 'l', the leading n by n
/// lower triangular part of the array C must contain the lower
/// triangular part of the symmetric matrix and the strictly
/// upper triangular part of C is not referenced. On exit, the
/// lower triangular part of the array C is overwritten by the
/// lower triangular part of the updated matrix.
///
/// @param[in] ldc On entry, LDC specifies the first dimension of C as declared
/// in the calling (sub) program. LDC must be at least max(1, n).
///
void dsyrk_(const char *uplo,
				const char *trans,
				const int *n,
				const int *k,
				const double *alpha,
				const double *a,
				const int *lda,
				const double *beta,
				double *c,
				const int *ldc);

///  DSYMV  performs the matrix-vector operation: y := alpha*A*x + beta*y
///  where alpha and beta are scalars, x and y are n element vectors and
///  A is an n by n symmetric matrix.
///
///  @param[in] uplo On entry, UPLO specifies whether the upper or lower
///           triangular part of the array A is to be referenced as
///           follows:
///              UPLO = 'U' or 'u'   Only the upper triangular part of A
///                                  is to be referenced.
///              UPLO = 'L' or 'l'   Only the lower triangular part of A
///                                  is to be referenced.
///
///  @param[in] n On entry, N specifies the order of the matrix A.
///           N must be at least zero.
///
///  @param[in] alpha On entry, ALPHA specifies the scalar alpha.
///
///  @param[in] a Array of DIMENSION ( LDA, n ).
///           Before entry with  UPLO = 'U' or 'u', the leading n by n
///           upper triangular part of the array A must contain the upper
///           triangular part of the symmetric matrix and the strictly
///           lower triangular part of A is not referenced.
///           Before entry with UPLO = 'L' or 'l', the leading n by n
///           lower triangular part of the array A must contain the lower
///           triangular part of the symmetric matrix and the strictly
///           upper triangular part of A is not referenced.
///
///  @param[in] lda On entry, LDA specifies the first dimension of A as declared
///           in the calling (sub) program. LDA must be at least
///           max( 1, n ).
///
///  @param[in] x Array of dimension at least
///           ( 1 + ( n - 1 )*abs( INCX ) ).
///           Before entry, the incremented array X must contain the n
///           element vector x.
///
///  @param[in] incx On entry, INCX specifies the increment for the elements of
///           X. INCX must not be zero.
///
///  @param[in] beta On entry, BETA specifies the scalar beta. When BETA is
///           supplied as zero then Y need not be set on input.
///
///  @param[in,out] y Array of dimension at least
///           ( 1 + ( n - 1 )*abs( INCY ) ).
///           Before entry, the incremented array Y must contain the n
///           element vector y. On exit, Y is overwritten by the updated
///           vector y.
///
///  @param[in] incy On entry, INCY specifies the increment for the elements of
///           Y. INCY must not be zero.
///
void dsymv_(const char *uplo,
				const int *n,
				const double *alpha,
				const double *a,
				const int *lda,
				const double *x,
				const int *incx,
				const double *beta,
				double *y,
				const int *incy);


/// DSYMM  performs one of the matrix-matrix operations
///
///  C := alpha*A*B + beta*C  or  C := alpha*B*A + beta*C
///
///  where alpha and beta are scalars, A is a symmetric matrix and B and
///  C are m by n matrices.
///
///  Arguments
///  ==========
///
///  @param[in] side
///           On entry,  SIDE  specifies whether  the  symmetric matrix  A
///           appears on the  left or right  in the  operation as follows:
///              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
///              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
///
///  @param[in] uplo
///           On entry, UPLO specifies whether the upper or lower
///           triangular part of the symmetric matrix  A is to be
///           referenced as follows:
///           UPLO = 'U' or 'u'  Only the upper triangular part of the symmetric matrix is to be referenced.
///           UPLO = 'L' or 'l'  Only the lower triangular part of the symmetric matrix is to be referenced.
///
///  @param[in] m
///           On entry, M specifies the number of rows of the matrix  C.
///           M must be at least zero.
///
///  @param[in] n
///           On entry, N specifies the number of columns of the matrix C.
///           N must be at least zero.
///
///  @param[in] alpha
///           On entry, ALPHA specifies the scalar alpha.
///
///  @param[in] a
///           Array of DIMENSION ( LDA, ka ), where ka is
///           m when  SIDE = 'L' or 'l'  and is  n otherwise.
///           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
///           the array  A  must contain the  symmetric matrix,  such that
///           when  UPLO = 'U' or 'u', the leading m by m upper triangular
///           part of the array  A  must contain the upper triangular part
///           of the  symmetric matrix and the  strictly  lower triangular
///           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
///           the leading  m by m  lower triangular part  of the  array  A
///           must  contain  the  lower triangular part  of the  symmetric
///           matrix and the strictly upper triangular part of  A  is not
///           referenced.
///           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
///           the array  A  must contain the  symmetric matrix,  such that
///           when  UPLO = 'U' or 'u', the leading n by n upper triangular
///           part of the array  A  must contain the upper triangular part
///           of the  symmetric matrix and the  strictly  lower triangular
///           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
///           the leading  n by n  lower triangular part  of the  array  A
///           must  contain  the  lower triangular part  of the  symmetric
///           matrix and the  strictly upper triangular part of  A  is not
///           referenced.
///
///  @param[in] lda
///           On entry, LDA specifies the first dimension of A as declared
///           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
///           LDA must be at least  max( 1, m ), otherwise  LDA must be at
///           least  max( 1, n ).
///
///  @param[in] b
///           Array of DIMENSION ( LDB, n ).
///           Before entry, the leading  m by n part of the array  B  must
///           contain the matrix B.
///
///  @param[in] ldb
///           On entry, LDB specifies the first dimension of B as declared
///           in  the  calling  (sub)  program.   LDB  must  be  at  least
///           max( 1, m ).
///
///  @param[in] beta
///           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
///           supplied as zero then C need not be set on input.
///
///  @param[in, out] c
///           Array of DIMENSION ( LDC, n ).
///           Before entry, the leading  m by n  part of the array  C must
///           contain the matrix  C,  except when  beta  is zero, in which
///           case C need not be set on entry.
///           On exit, the array  C  is overwritten by the  m by n updated
///           matrix.
///
///  @param[in] ldc
///           On entry, LDC specifies the first dimension of C as declared
///           in  the  calling  (sub)  program.   LDC  must  be  at  least
///           max( 1, m ).
///
void dsymm_(const char *side,
				const char *uplo,
				const int *m,
				const int *n,
				const double *alpha,
				const double *a,
				const int *lda,
				const double *b,
				const int *ldb,
				const double *beta,
				double *c,
				const int *ldc);


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


/// Computes the Euclidean norm of a vector.
/// 
/// DNRM2 computes the Euclidean (L2) norm of a double precision real vector
///
/// @param[in]      n  Number of elements in the operand vector.
/// @param[in]      dx  Array of dimension (n-1) * |incx| + 1. Array x contains the operand vector.
/// @param[in]      incx  Increment between elements of x.
/// 
/// @return Resulting Euclidean norm.
/// 
double dnrm2_(const int *n,
				const double *dx,
				const int *incx);


/// Scales a double precision vector.
/// DSCAL scales a double precision vector with a double precision scalar.
/// DSCAL scales the vector x of length n and increment incx  by  the  constant alpha.
///
/// This routine performs the following vector operation:  x <-- alpha x
/// where alpha is a double precision scalar, and x is a double precision vector.
///
/// @param[in] n  Number of elements in the vector. If n <= 0, this routine returns without computation.
///
/// @param[in]  alpha   The scaling value.
///
/// @param[in,out]  x Array of dimension (n-1) * |incx| + 1. Vector to be scaled.
///
/// @param[in]  incx Increment between elements of x. If incx = 0, the results will be unpredictable.
///
void dscal_(const int *n,
			const double *alpha,
			double *x,
			const int *incx);



/// copy a vector x into another location y.
///
/// @param[in] n Vector length
/// @param[in] dx The first vector
/// @param[in] incx The stride of the vector dx
/// @param[in] dy The second vector
/// @param[in] incy The stride of the vector dy
///
void dcopy_(const int *n,
				const double *dx,
				const int *incx,
				double *dy,
				const int *incy);



/// Scale by alpha a vector x and add it to a vector y 
///
/// This routine performs the following vector operation:  y <-- alpha x + y
/// where alpha is a double precision scalar, x and y are double precision vectors.
///
/// @param[in] n  Number of elements in the vector. If n <= 0, this routine returns without computation.
///
/// @param[in]  alpha   The scaling value.
///
/// @param[in]  x Array of dimension (n-1) * |incx| + 1.
/// @param[in]  incx Increment between elements of x. If incx = 0, the results will be unpredictable.
///
/// @param[in,out]  y Array of dimension (n-1) * |incy| + 1. Vector which contains the result.
/// @param[in]  incy Increment between elements of x. If incy = 0, the results will be unpredictable.
///

void daxpy_(const int *n,
			const double *alpha,
			const double *x,
			const int *incx,
			double *y,
			const int *incy);


/// performs one of the matrix-vector operations
/// y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
///
/// where alpha and beta are scalars, x and y are vectors and A is an
/// m by n matrix.
///
/// @param[in] trans: use the transpose of the matrix A if set to 'T', 't', 'C' or 'c'
///					  use the matrix A if set to 'N' or 'n'
///
/// @param[in] m Number of rows in A.
/// @param[in] n Number of columns in A.
///
/// @param[in]  alpha   The scaling value.
/// @param[in]  A   The m by n matrix.
/// @param[in]  LDA On entry, LDA specifies the first dimension of A as declared
///       		in the calling (sub) program. LDA must be at least
///				max( 1, m )..
///
/// @param[in]  x Array of dimension (n(or m)-1) * |incx| + 1 when trans='N' (or 'T')
/// @param[in]  incx Increment between elements of x. If incx = 0, the results will be unpredictable.
///
/// @param[in]  beta   The scaling value.
/// @param[in,out]  y Array of dimension (m(or n)-1) * |incy| + 1 when trans='N' (or 'T'). Vector which contains the result.
/// @param[in]  incy Increment between elements of x. If incy = 0, the results will be unpredictable.
///
void dgemv_(const char *trans,
			const int *m,
			const int *n,
			const double *alpha,
			const double *A,
			const int *LDA,
			const double *x,
			const int *incx,
			const double *beta,
			double *y,
			const int *incy); 	




/// dgemm_ performs one of the matrix-matrix operations
///    C := alpha*op( A )*op( B ) + beta*C,
/// where  op( X ) is one of
///    op( X ) = X   or   op( X ) = X**T,
///
/// alpha and beta are scalars, and A, B and C are matrices, with op( A )
/// an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
///
///  Arguments:
///  ==========
///
/// @param[in] TRANSA
///          TRANSA is CHARACTER*1
///           On entry, TRANSA specifies the form of op( A ) to be used in
///           the matrix multiplication as follows:
///              TRANSA = 'N' or 'n',  op( A ) = A.
///              TRANSA = 'T' or 't',  op( A ) = A**T.
///              TRANSA = 'C' or 'c',  op( A ) = A**T.
///
/// @param[in] TRANSB
///          TRANSB is CHARACTER*1
///           On entry, TRANSB specifies the form of op( B ) to be used in
///           the matrix multiplication as follows:
///              TRANSB = 'N' or 'n',  op( B ) = B.
///              TRANSB = 'T' or 't',  op( B ) = B**T.
///              TRANSB = 'C' or 'c',  op( B ) = B**T.
///
/// @param[in] M
///          M is INTEGER
///           On entry,  M  specifies  the number  of rows  of the  matrix
///           op( A )  and of the  matrix  C.  M  must  be at least  zero.
///
/// @param[in] N
///          N is INTEGER
///           On entry,  N  specifies the number  of columns of the matrix
///           op( B ) and the number of columns of the matrix C. N must be
///           at least zero.
///
/// @param[in] K
///          K is INTEGER
///           On entry,  K  specifies  the number of columns of the matrix
///           op( A ) and the number of rows of the matrix op( B ). K must
///           be at least  zero.
///
/// @param[in] ALPHA
///          ALPHA is DOUBLE PRECISION.
///           On entry, ALPHA specifies the scalar alpha.
///
/// @param[in] A
///          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
///           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
///           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
///           part of the array  A  must contain the matrix  A,  otherwise
///           the leading  k by m  part of the array  A  must contain  the
///           matrix A.
///
/// @param[in] LDA
///          LDA is INTEGER
///           On entry, LDA specifies the first dimension of A as declared
///           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
///           LDA must be at least  max( 1, m ), otherwise  LDA must be at
///           least  max( 1, k ).
///
/// @param[in] B
///          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
///           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
///           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
///           part of the array  B  must contain the matrix  B,  otherwise
///           the leading  n by k  part of the array  B  must contain  the
///           matrix B.
///
/// @param[in] LDB
///          LDB is INTEGER
///           On entry, LDB specifies the first dimension of B as declared
///           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
///           LDB must be at least  max( 1, k ), otherwise  LDB must be at
///           least  max( 1, n ).
///
/// @param[in] BETA
///          BETA is DOUBLE PRECISION.
///           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
///           supplied as zero then C need not be set on input.
///
/// @param[in,out] C
///          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
///           Before entry, the leading  m by n  part of the array  C must
///           contain the matrix  C,  except when  beta  is zero, in which
///           case C need not be set on entry.
///           On exit, the array  C  is overwritten by the  m by n  matrix
///           ( alpha*op( A )*op( B ) + beta*C ).
///
/// @param[in] LDC
///          LDC is INTEGER
///           On entry, LDC specifies the first dimension of C as declared
///           in  the  calling  (sub)  program.   LDC  must  be  at  least
///           max( 1, m ).
///
void dgemm_(const char *transA
	,const char *transB
	,const int *M
	,const int *N
	,const int *K
	,const double *ALPHA
	,const double *A
	,const int *LDA
	,const double *B
	,const int *LDB
	,const double *BETA
	,double *C
	,const int *LDC);


#ifdef __cplusplus
}
#endif

// Constants used by blas and lapack routines
static const int	I0 = 0;		///< Integer zero.
static const int	I1 = 1;		///< Integer one.
static const int	I2 = 2;		///< Integer two.
static const double	D0 = 0.;	///< Float double zero.
static const double	D1 = 1.;	///< Float double one.
static const double minus_one(-1.);	///< Float double -1

#endif

