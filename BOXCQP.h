
#ifndef BOXCQP_H
#define BOXCQP_H

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


/// BOXCQP class
/// solves a box constrained convex quadratic program using the BOXCQP algorithm
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
		  ,double* aUpperBound) :
			mN(aN),
			mSpace(),
			mRHS(NULL),
			mLHS(NULL),
			mB(NULL),
			mLambda(NULL),
			mMu(NULL),
			mSets(),
			mListLset(),
			mListUset(),
			mListSset(),
			mDWorkSpace(),
			mIWorkspace(),
			mSubmatrixFact(NULL),
			mDiagScaling(NULL),
			mSubSolution(NULL),
			mLowerBounds(aLowerBound),
			mUpperBounds(aUpperBound)
	{
		allocateMemory();
	}
	
	
	/// solveQP
	/// solve positive definite quadratic program using the BOXCQP algorithm. Returns true if converged.
	///
	/// @param[in]		aB the matrix B
	/// @param[in]		aD the d vector
	/// @param[in]		aLDB the leading dimension of matrix B
	/// @param[out]		aX the solution vector
	/// @param[out]		aSolutionOnBorder true if the solution is on the bounds, false otherwise
	/// @param[in,out]	aUnconstrainedDirection	Solution for the unconstrained problem if not NULL.
	///
	bool solveQP(const double *aB, const double *aD, const int *aLDB, double *aX, bool *aSolutionOnBorder, double *aUnconstrainedDirection = NULL);
	
	
private:
	
	/// allocateMemory
	/// allocate the space for work
	///
	void allocateMemory(void);
	
	/// updateSets
	/// update the vector representing the sets
	///
	/// @param[in] aX The current solution vector
	///
	void updateSets(double *aX);
	
private:
	
	/// ActiveSet
	/// defines the sets in which each variable is
	///
	enum ActiveSet
	{
		LSET,	///< set {i: xi < ai , or xi = ai and lambda_i >= 0}	
		USET,	///< set {i: xi > bi , or xi = bi and mu_i >= 0}
		SSET	///< set {i: ai < xi < bi , or xi = ai and lambda_i < 0 , or xi = bi and mu_i < 0}
	};

private:

	int 					mN;				///< size of the problem
	
	std::vector<double>		mSpace;			///< workspace
	double*					mRHS;			///< right hand side vector used for the linear system 
	double*					mLHS;			///< left hand side matrix used for the linear system
	
	double*					mB;				///< B matrix
	double*					mLambda;		///< Lagrange variables, lower constraints
	double*					mMu;			///< Lagrange variables, upper constraints
	
	std::vector<ActiveSet>	mSets;			///< sets in which each variable belongs

	std::vector<int>		mListLset;		///< list containing the L set
	std::vector<int>		mListUset;		///< list containing the U set
	std::vector<int>		mListSset;		///< list containing the S set
	
	double*					mSubmatrixFact;	///< submatrix (mLHS) factorized
	double*					mDiagScaling;	///< diagonal scaling for reducing the condition number of the submatrix
	double*					mSubSolution;	///< solution of subsystems
	
	std::vector<double>		mDWorkSpace;	///< double workspace for lapack
	std::vector<int>		mIWorkspace;	///< int workspace for lapack

	double*					mLowerBounds;	///< lower bounds (constraints)
	double*					mUpperBounds;	///< upper bounds (constraints)
};

#endif // BOXCQP_H

