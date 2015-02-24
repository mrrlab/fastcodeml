
#ifndef LAPACK_H
#define LAPACK_H

///  DSYEVR computes selected eigenvalues and, optionally, eigenvectors
///  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
///  selected by specifying either a range of values or a range of
///  indices for the desired eigenvalues.
///
///
///  @param[in] jobz    (input) CHARACTER*1
///          = 'N':  Compute eigenvalues only;
///          = 'V':  Compute eigenvalues and eigenvectors.
///
///  @param[in] range   (input) CHARACTER*1
///          = 'A': all eigenvalues will be found.
///          = 'V': all eigenvalues in the half-open interval (VL,VU]
///                 will be found.
///          = 'I': the IL-th through IU-th eigenvalues will be found.
///
///  @param[in] uplo    (input) CHARACTER*1
///          = 'U':  Upper triangle of A is stored;
///          = 'L':  Lower triangle of A is stored.
///
///  @param[in] n       (input) INTEGER
///          The order of the matrix A.  N >= 0.
///
///  @param[in,out] a       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
///          On entry, the symmetric matrix A.  If UPLO = 'U', the
///          leading N-by-N upper triangular part of A contains the
///          upper triangular part of the matrix A.  If UPLO = 'L',
///          the leading N-by-N lower triangular part of A contains
///          the lower triangular part of the matrix A.
///          On exit, the lower triangle (if UPLO='L') or the upper
///          triangle (if UPLO='U') of A, including the diagonal, is
///          destroyed.
///
///  @param[in] lda     (input) INTEGER
///          The leading dimension of the array A.  LDA >= max(1,N).
///
///  @param[in] vl      (input) DOUBLE PRECISION
///  @param[in] vu      (input) DOUBLE PRECISION
///          If RANGE='V', the lower and upper bounds of the interval to
///          be searched for eigenvalues. VL < VU.
///          Not referenced if RANGE = 'A' or 'I'.
///
///  @param[in] il      (input) INTEGER
///  @param[in] iu      (input) INTEGER
///          If RANGE='I', the indices (in ascending order) of the
///          smallest and largest eigenvalues to be returned.
///          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
///          Not referenced if RANGE = 'A' or 'V'.
///
///  @param[in] abstol  (input) DOUBLE PRECISION
///          The absolute error tolerance for the eigenvalues.
///          An approximate eigenvalue is accepted as converged
///          when it is determined to lie in an interval [a,b]
///          of width less than or equal to ABSTOL + EPS *   max( |a|,|b| ),
///          where EPS is the machine precision.  If ABSTOL is less than
///          or equal to zero, then  EPS*|T|  will be used in its place,
///          where |T| is the 1-norm of the tridiagonal matrix obtained
///          by reducing A to tridiagonal form.
///          See "Computing Small Singular Values of Bidiagonal Matrices
///          with Guaranteed High Relative Accuracy," by Demmel and
///          Kahan, LAPACK Working Note #3.
///          If high relative accuracy is important, set ABSTOL to
///          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
///          eigenvalues are computed to high relative accuracy when
///          possible in future releases.  The current code does not
///          make any guarantees about high relative accuracy, but
///          future releases will. See J. Barlow and J. Demmel,
///          "Computing Accurate Eigensystems of Scaled Diagonally
///          Dominant Matrices", LAPACK Working Note #7, for a discussion
///          of which matrices define their eigenvalues to high relative
///          accuracy.
///
///  @param[out] m       (output) INTEGER
///          The total number of eigenvalues found.  0 <= M <= N.
///          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
///
///  @param[out] w       (output) DOUBLE PRECISION array, dimension (N)
///          The first M elements contain the selected eigenvalues in
///          ascending order.
///
///  @param[out] z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
///          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
///          contain the orthonormal eigenvectors of the matrix A
///          corresponding to the selected eigenvalues, with the i-th
///          column of Z holding the eigenvector associated with W(i).
///          If JOBZ = 'N', then Z is not referenced.
///          Note: the user must ensure that at least max(1,M) columns are
///          supplied in the array Z; if RANGE = 'V', the exact value of M
///          is not known in advance and an upper bound must be used.
///          Supplying N columns is always safe.
///
///  @param[in] ldz     (input) INTEGER
///          The leading dimension of the array Z.  LDZ >= 1, and if
///          JOBZ = 'V', LDZ >= max(1,N).
///
///  @param[out] isuppz  (output) INTEGER array, dimension ( 2*max(1,M) )
///          The support of the eigenvectors in Z, i.e., the indices
///          indicating the nonzero elements in Z. The i-th eigenvector
///          is nonzero only in elements ISUPPZ( 2*i-1 ) through
///          ISUPPZ( 2*i ). Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
///
///  @param[out] work    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
///          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
///
///  @param[in] lwork   (input) INTEGER
///          The dimension of the array WORK.  LWORK >= max(1,26*N).
///          For optimal efficiency, LWORK >= (NB+6)*N,
///          where NB is the max of the blocksize for DSYTRD and DORMTR
///          returned by ILAENV.
///          If LWORK = -1, then a workspace query is assumed; the routine
///          only calculates the optimal size of the WORK array, returns
///          this value as the first entry of the WORK array, and no error
///          message related to LWORK is issued by XERBLA.
///
///  @param[out] iwork   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
///          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
///
///  @param[in] liwork  (input) INTEGER
///          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
///          If LIWORK = -1, then a workspace query is assumed; the
///          routine only calculates the optimal size of the IWORK array,
///          returns this value as the first entry of the IWORK array, and
///          no error message related to LIWORK is issued by XERBLA.
///
///  @param[out] info 
///          = 0:  successful exit
///          < 0:  if INFO = -i, the i-th argument had an illegal value
///          > 0:  Internal error
///
extern "C" void dsyevr_(const char *jobz,
                        const char *range,
                        const char *uplo,
                        const int *n,
                        double *a,
                        const int *lda,
                        const double *vl,
                        const double *vu,
                        const int *il,
                        const int *iu,
                        const double *abstol,
                        int *m,
                        double *w,
                        double *z,
                        const int *ldz,
                        int *isuppz,
                        double *work,
                        const int *lwork,
                        int *iwork,
                        const int *liwork,
                        int *info);



///   DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
///	  real symmetric matrix A. If eigenvectors are desired, it uses a
///	  divide and conquer algorithm.
///	
///	  The divide and conquer algorithm makes very mild assumptions about
///	  floating point arithmetic. It will work on machines with a guard
///	  digit in add/subtract, or on those binary machines without guard
///	  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
///	  Cray-2. It could conceivably fail on hexadecimal or decimal machines
///	  without guard digits, but we know of none.
///	
///	  Because of large use of BLAS of level 3, DSYEVD needs N**2 more
///	  workspace than DSYEVX.
///	
///	  Arguments
///	  =========
///	
///	  @param[in] jobz    (input) CHARACTER*1
///	          = 'N':  Compute eigenvalues only;
///	          = 'V':  Compute eigenvalues and eigenvectors.
///	
///	  @param[in] uplo    (input) CHARACTER*1
///	          = 'U':  Upper triangle of A is stored;
///	          = 'L':  Lower triangle of A is stored.
///	
///	  @param[in] n       (input) INTEGER
///	          The order of the matrix A.  N >= 0.
///	
///	  @param[in,out] a       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
///	          On entry, the symmetric matrix A.  If UPLO = 'U', the
///	          leading N-by-N upper triangular part of A contains the
///	          upper triangular part of the matrix A.  If UPLO = 'L',
///	          the leading N-by-N lower triangular part of A contains
///	          the lower triangular part of the matrix A.
///	          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
///	          orthonormal eigenvectors of the matrix A.
///	          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
///	          or the upper triangle (if UPLO='U') of A, including the
///	          diagonal, is destroyed.
///	
///	  @param[in] lda     (input) INTEGER
///	          The leading dimension of the array A.  LDA >= max(1,N).
///	
///	  @param[out] w       (output) DOUBLE PRECISION array, dimension (N)
///	          If INFO = 0, the eigenvalues in ascending order.
///	
///	  @param[out] work    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
///	          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
///	
///	  @param[in] lwork   (input) INTEGER
///	          The dimension of the array WORK.
///	          If N <= 1,               LWORK must be at least 1.
///	          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
///	          If JOBZ = 'V' and N > 1, LWORK must be at least 1 + 6*N + 2*N**2.
///	          If LWORK = -1, then a workspace query is assumed; the routine
///	          only calculates the optimal sizes of the WORK and IWORK
///	          arrays, returns these values as the first entries of the WORK
///	          and IWORK arrays, and no error message related to LWORK or
///	          LIWORK is issued by XERBLA.
///	
///	  @param[out] iwork   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
///	          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
///	
///	  @param[in] liwork  (input) INTEGER
///	          The dimension of the array IWORK.
///	          If N <= 1,                LIWORK must be at least 1.
///	          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
///	          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
///	          If LIWORK = -1, then a workspace query is assumed; the
///	          routine only calculates the optimal sizes of the WORK and
///	          IWORK arrays, returns these values as the first entries of
///	          the WORK and IWORK arrays, and no error message related to
///	          LWORK or LIWORK is issued by XERBLA.
///	
///	  @param[out] info    (output) INTEGER
///	          = 0:  successful exit
///	          < 0:  if INFO = -i, the i-th argument had an illegal value
///	          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
///	                to converge; i off-diagonal elements of an intermediate
///	                tridiagonal form did not converge to zero;
///	                if INFO = i and JOBZ = 'V', then the algorithm failed
///	                to compute an eigenvalue while working on the submatrix
///	                lying in rows and columns INFO/(N+1) through
///	                mod(INFO,N+1).
///	
extern "C" void dsyevd_(const char *jobz,
                        const char *uplo,
                        const int *n,
                        double *a,
                        const int *lda,
                        double *w,
                        double *work,
                        const int *lwork,
                        int *iwork,
                        const int *liwork,
                        int *info);
                        
                        

/// The pdgeqrf routine forms the QR factorization of a general m-by-n distributed matrix 
/// sub(A)= A(ia:ia+m-1,ja:ja+n-1) as A=Q*R
///
/// @params[in] m  (global) INTEGER. The number of rows in the submatrix sub(Q) (m ≥ 0).
///
/// @params[in] n  (global) INTEGER. The number of columns in the submatrix sub(Q) (m ≥ n ≥ 0).
///
/// @params[in,out] a  (local)  DOUBLE PRECISION 
///			in:	   Pointer into the local memory to an array of local dimension (lld_a, LOCc(ja+n-1)). 
///				   Contains the local pieces of the distributed matrix sub(A) to be factored.
///			out:   The elements on and above the diagonal of sub(A) contain the min(m,n)-by-n upper 
///				   trapezoidal matrix R (R is upper triangular if m ≥ n); the elements below the diagonal, 
///				   with the array tau, represent the orthogonal/unitary matrix Q as a product of elementary 
///				   reflectors (see Application Notes below).
///
/// @params[in] ia
/// @params[in] ja (global) INTEGER. 
/// 			The row and column indices in the global array a indicating the first row and the first 
///				column of the submatrix A(ia:ia+m-1,ja:ja+n-1), respectively.
///
/// @params[in] desca (global and local) INTEGER array, dimension (dlen_). 
///				The array descriptor for the distributed matrix A.
///
/// @params[out] tau (local) DOUBLE PRECISION Array, DIMENSION LOCc(ja+k-1).
///				Contains the scalar factor tau of elementary reflectors. tau is tied to the distributed matrix A.
///
/// @params[in,out] work (local) DOUBLE PRECISION
///			    Workspace array of dimension of lwork.
///				On exit, work(1) contains the minimum value of lwork required for optimum performance.
///
/// @params[in] lwork (local or global) INTEGER, dimension of work.
///				Must be at least lwork ≥ nb_a*(nqa0 + mpa0 + nb_a), where
///			    iroffa = mod(ia-1, mb_a), icoffa = mod(ja-1, nb_a),
///				iarow = indxg2p(ia, mb_a, MYROW, rsrc_a, NPROW),
///				iacol = indxg2p(ja, nb_a, MYCOL, csrc_a, NPCOL),
///				mpa0 = numroc(m+iroffa, mb_a, MYROW, iarow, NPROW),
///				nqa0 = numroc(n+icoffa, nb_a, MYCOL, iacol, NPCOL);
///				indxg2p and numroc are ScaLAPACK tool functions; MYROW, MYCOL, NPROW and NPCOL can be 
///				determined by calling the subroutine blacs_gridinfo.
///				If lwork = -1, then lwork is global input and a workspace query is assumed; 
///				the routine only calculates the minimum and optimal size for all work arrays. 
///				Each of these values is returned in the first entry of the corresponding work array, 
///				and no error message is issued by pxerbla.
///
/// @params[out] info (global) INTEGER.
///					= 0: the execution is successful.
///					< 0: if the i-th argument is an array and the j-entry had an illegal value, 
///					then info = - (i* 100+j), if the i-th argument is a scalar and had an illegal value, 
///					then info = -i.
///
extern "C" void pdgeqrf(const int *m
					   ,const int *n
					   ,double *A
					   ,const int *IA
					   ,const int *JA
					   ,const int *desca
					   ,double *tau
					   ,double *work
					   ,const int *lwork
					   ,int *info);



/// The pdorgqr routine generates the whole or part of m-by-n real distributed matrix Q denoting 
/// A(ia:ia+m-1, ja:ja+n-1) with orthonormal columns, which is defined as the first n columns of 
/// a product of k elementary reflectors of order m Q= H(1)*H(2)*...*H(k) as returned by p?geqrf.
///
/// @params[in] m  (global) INTEGER. The number of rows in the submatrix sub(Q) (m ≥ 0).
///
/// @params[in] n  (global) INTEGER. The number of columns in the submatrix sub(Q) (m ≥ n ≥ 0).
///
/// @params[in] k  (global) INTEGER. The number of elementary reflectors whose product defines 
/// 			   					 the matrix Q (n ≥ k ≥ 0).
///
/// @params[in,out] a  (local)  DOUBLE PRECISION 
///			in:	   Pointer into the local memory to an array of local dimension (lld_a, LOCc(ja+n-1)). 
///				   The j-th column must contain the vector which defines the elementary reflector 
///				   H(j), ja≤j≤ja +k-1, as returned by p?geqrf in the k columns of its distributed 
///				   matrix argument A(ia:*, ja:ja+k-1).
///			out:   Contains the local pieces of the m-by-n distributed matrix Q.
///
/// @params[in] ia
/// @params[in] ja (global) INTEGER. 
/// 			The row and column indices in the global array a indicating the first row and the first 
///				column of the submatrix A(ia:ia+m-1,ja:ja+n-1), respectively.
///
/// @params[in] desca (global and local) INTEGER array, dimension (dlen_). 
///				The array descriptor for the distributed matrix A.
///
/// @params[in] tau (local) DOUBLE PRECISION Array, DIMENSION LOCc(ja+k-1).
///				Contains the scalar factor tau (j) of elementary reflectors H(j) as returned by p?geqrf. 
///				tau is tied to the distributed matrix A.
///
/// @params[in,out] work (local) DOUBLE PRECISION
///			    Workspace array of dimension of lwork.
///				On exit, work(1) contains the minimum value of lwork required for optimum performance.
///
/// @params[in] lwork (local or global) INTEGER, dimension of work.
///				Must be at least lwork ≥ nb_a*(nqa0 + mpa0 + nb_a), where
///			    iroffa = mod(ia-1, mb_a), icoffa = mod(ja-1, nb_a),
///				iarow = indxg2p(ia, mb_a, MYROW, rsrc_a, NPROW),
///				iacol = indxg2p(ja, nb_a, MYCOL, csrc_a, NPCOL),
///				mpa0 = numroc(m+iroffa, mb_a, MYROW, iarow, NPROW),
///				nqa0 = numroc(n+icoffa, nb_a, MYCOL, iacol, NPCOL);
///				indxg2p and numroc are ScaLAPACK tool functions; MYROW, MYCOL, NPROW and NPCOL can be 
///				determined by calling the subroutine blacs_gridinfo.
///				If lwork = -1, then lwork is global input and a workspace query is assumed; 
///				the routine only calculates the minimum and optimal size for all work arrays. 
///				Each of these values is returned in the first entry of the corresponding work array, 
///				and no error message is issued by pxerbla.
///
/// @params[out] info (global) INTEGER.
///					= 0: the execution is successful.
///					< 0: if the i-th argument is an array and the j-entry had an illegal value, 
///					then info = - (i* 100+j), if the i-th argument is a scalar and had an illegal value, 
///					then info = -i.
///
extern "C" void pdorgqr(const int *m
					   ,const int *n
					   ,const int *k
					   ,double *A
					   ,const int *IA
					   ,const int *JA
					   ,const int *desca
					   ,double *tau
					   ,double *work
					   ,const int *lwork
					   ,int *info);


#endif

