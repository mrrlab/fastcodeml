
#ifndef BOXCQP_H
#define BOXCQP_H

#include <vector>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif


/// BOXCQP class
/// solves a box constrained quadratic program using the BOXCQP algorithm
/// problems of type
/// 			find min 1/2 x^t B x + d^t x 
///				      x
///				s.t. a <= x <= b
/// see http://cs.uoi.gr/~lagaris/papers/PREPRINTS/boxcqp.pdf
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-30 (initial version)
///     @version 1.1
///

// uncomment this to solve reduced linear systems, at the cost of copying the 
// values sequentially
#define USE_SUBMATRIX_QP

class BOXCQP
{
public:
	/// Constructor
	///
	/// @param[in] aN The size of the problem
	///	@param[in] aLowerBound the vector containing the lower bounds
	///	@param[in] aUpperBound the vector containing the upper bounds
	///
	BOXCQP(const int& aN
		  ,double* aLowerBound
		  ,double* aUpperBound)
		:mN(aN)
		,ma(aLowerBound)
		,mb(aUpperBound)
	{
		alocateMemory();
	}
	
	
	/// solveQP
	/// solve the quadratic program using the BOXCQP algorithm
	///
	/// @param[in]  B the matrix B
	/// @param[in]  d the d vector
	/// @param[in]	LDA the leading dimension of matrix B
	/// @param[out] x the solution vector
	/// @param[out] aSolutionOnBorder true if the solution is on the bounds, false otherwise
	///
	void solveQP(const double *B, const double *d, const int *LDA, double *x, bool *aSolutionOnBorder);
	
	
private:
	
	/// alocateMemory
	/// alocate the space for work
	///
	void alocateMemory(void);
	
	/// updateSets
	/// update the vector representing the sets
	///
	/// @param[in] ax The current solution vector
	///
	void updateSets(double *ax);
	
private:
	
	/// activeSet
	/// defines the sets in which each variable is
	///
	enum ActiveSet
	{
		LSET,	//< set {i: xi < ai , or xi = a1 and lambda_i >= 0}	
		USET,	//< set {i: xi > bi , or xi = b1 and mu_i >= 0}
		SSET	//< set {i: ai < xi < bi , or xi = a1 and lambdai < 0 , or xi = b1 and mu_i < 0}
	};

private:

	int 					mN;				//< size of the problem
	
	std::vector<double>		mSpace;			//< workspace
	double*					mRHS;			//< right hand side vector used for the linear system 
	double*					mLHS;			//< left hand side matrix used for the linear system
	
	double*					mB;				//< B matrix
	double*					md;				//< d vector
	double*					mLambda;		//< Lagrange variables, lower constraints
	double*					mMu;			//< Lagrange variables, upper constraints
	
	std::vector<ActiveSet>	mSets;			//< sets in which each variable belongs
#ifdef USE_SUBMATRIX_QP
	std::vector<int>		mListLset;		//< list containing the L set
	std::vector<int>		mListUset;		//< list containing the U set
	std::vector<int>		mListSset;		//< list containing the S set
#else
	double*					mx_known;		//<	helper for solving linear system
	double*					mMu_known;		//<	helper for solving linear system
	double*					mLambda_known;	//<	helper for solving linear system
#endif //USE_SUBMATRIX_QP

	double*					ma;				//< lower bounds (constraints)
	double*					mb;				//< upper bounds (constraints)
};

#endif // BOXCQP_H

