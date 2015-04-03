
#ifndef OPTSQP_H
#define OPTSQP_H


#include <cstdio>
#include <vector>
#include <memory>
#include "BranchSiteModel.h"



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
	/// @param[out] x the solution vector
	///
	void solveQP(const double *B, const double *d, double *x);
	
	
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
	double*					mRHS;			//< right hand side matrix used for the linear system 
	double*					mLHS;			//< left hand side vector used for the linear system
	double*					mx_known;		//<	helper for solving linear system
	double*					mMu_known;		//<	helper for solving linear system
	double*					mLambda_known;	//<	helper for solving linear system
	
	double*					mB;				//< B matrix
	double*					md;				//< d vector
	double*					mLambda;		//< Lagrange variables, lower constraints
	double*					mMu;			//< Lagrange variables, upper constraints
	std::vector<ActiveSet>	mSets;			//< sets in which each variable belongs
	double*					ma;				//< lower bounds (constraints)
	double*					mb;				//< upper bounds (constraints)
};



// uncomment to use strong wolfe conditions as a stopping criterion for the line search
// comment it to use only the first Wolfe condition
//#define STRONG_WOLFE_LINE_SEARCH

/// OptSQP class.
/// sequential quadratic programming optimizer
/// see http://www.neos-guide.org/content/sequential-quadratic-programming
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-30 (initial version)
///     @version 1.1
///

class OptSQP
{
public:
	/// Constructor
	/// 
	/// @param[in] aModel				The pointer to the hypothesis class that will be used
	/// @param[in] aTrace				Trace or not the optimizer progress
	/// @param[in] aVerbose				The verbose level
	/// @param[in] aLowerBound			Lower limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aUpperBound			Upper limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aAbsoluteError		Absolute error to stop computation
	/// @param[in] aStopIfBigger		If true stop computation as soon as value is over aThreshold
	/// @param[in] aThreshold			The threshold at which the maximization should be stopped
	/// @param[in] aMaxIterations		Maximum number of iterations for the maximization
	/// @param[in] aNumTimes			Number of branches
	/// 
	OptSQP(BranchSiteModel* aModel
		  ,bool aTrace
		  ,unsigned int aVerbose
		  ,const std::vector<double>& aLowerBound
		  ,const std::vector<double>& aUpperBound
		  ,double aAbsoluteError
		  ,bool aStopIfBigger
		  ,double aThreshold
		  ,int aMaxIterations
		  ,int aNumTimes) 
		:mModel(aModel)
		,mTrace(aTrace)
		,mTraceFun(aTrace)
		,mLowerBound(aLowerBound)
		,mUpperBound(aUpperBound)
		,mAbsoluteError(aAbsoluteError)
		,mVerbose(aVerbose)
		,mStopIfBigger(aStopIfBigger)
		,mThreshold(-aThreshold)
		,mMaxIterations(aMaxIterations)
		,mNumTimes(aNumTimes)
		,mN(0)
		,mStep(0)
		{}
	
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars);
	

private:

	/// alocateMemory
	/// alocate the space for storage
	///
	void alocateMemory(void);
	
	
	/// SQPminimizer
	/// performs a sequential quadratic approximation of the function to estimate its minimum
	///	uses quadratic programming to solve local constrained subproblems 
	/// 
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void SQPminimizer(double *f, double *x);
	
	/// evaluateFunction
	///	evaluates the function at point x
	///
	/// @param[in]	x the point at which one evaluates the function
	/// @param[in]	aTrace The trace 
	///
	/// @return the function value
	///
	/// @exception nlopt::forced_stop To force halt the maximization because LRT is already not satisfied
	///
	double evaluateFunction(const double *x, bool aTrace);
	
	/// evaluateFunctionForLineSearch
	///	evaluates the function at point x + alpha*mP
	///
	/// @param[in]	x the point x
	/// @param[in]	alpha the step length
	///
	/// @return the function value
	///
	/// @exception nlopt::forced_stop To force halt the maximization because LRT is already not satisfied
	///
	double evaluateFunctionForLineSearch(const double* x, double alpha);
	
	/// computeGradient
	///	compute the gradient at point x using finite differences aproximation
	///
	/// @param[in]	x the point at which one evaluates the gradient
	/// @param[in]	f0 The value at point x
	/// @param[out]	aGrad The gradient 
	///
	void computeGradient(const double *x, double f0, double *aGrad);
	
	/// BFGSupdate
	/// perform the BFGS formula to update the hessian matrix.
	/// uses a modified version of BFGS to keep the matrix positive definite
	///
	void BFGSupdate(void);
	
	/// lineSearch
	/// perform a line search in th mP direction
	/// see http://pages.cs.wisc.edu/~ferris/cs730/chap3.pdf
	///
	/// @param[in,out] aalpha 	in: initial guess of step length
	///							out: step length
	/// @param[in,out] x		in: the original position
	///							out: if success, the new position
	/// @param[in,out] f		in: the original function value
	///							out: if success, the new value
	///
	void lineSearch(double *aalpha, double *x, double *f);
	
#ifdef STRONG_WOLFE_LINE_SEARCH
	/// zoom
	///	subroutine used in the lineSearch for the strong Wolfe conditions
	/// see http://pages.cs.wisc.edu/~ferris/cs730/chap3.pdf pp.60-61 for
	/// more informations
	///
	/// @param[in]	low			Low value of the step length
	/// @param[in]	high		High value of the step length
	/// @param[in]	x			The current position
	/// @param[in]	phi_0		The function value at point x
	/// @param[in]	phi_0_prime	The derivative in the mP direction at point x
	/// @param[in]	c1			Constant for the strong Wolfe conditions
	/// @param[in]	c2			Constant for the strong Wolfe conditions
	///
	/// @return the optimal step length
	///
	double zoom(double low, double high, const double *x, double const& phi_0, double const& phi_0_prime, const double& c1, const double& c2);
#endif

private:
	
	//double (f*)(unsigned, const double* x, double* grad, void* f_data); ///< Function to minimize
	
	int 						mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	std::vector<double>			mSpace;				///< Work and storage space
	std::vector<double> 		mXEvaluator;		///< Workspace for function evaluations
	
	double*						mGradient;			///< Gradient of the function. mN components
	double*						mHessian;			///< positive definite hessian approximation using BFGS; mN*mN components
	
	double*						mP;					///< search direction
	
	double*						mSk;				///< position change, i.e. mSk = xk - xk-1
	double*						mYk;				///< gradient change, i.e. mYk = dfk - dfk-1
	
	double*						mXPrev;				///< previous position
	double*						mGradPrev;			///< previous gradient
	
	double*						mWorkSpaceVect;		///< workspace. size of a mN vector.
	double*						mWorkSpaceMat;		///< workspace. size of a mN by mN matrix
	
	int							mStep;				///< current step	
	int							mNumCallZoom;		///< step of the zoom iteration in line search
	
	std::auto_ptr<BOXCQP>		mQPsolver;			///< box constrained quadratic program solver

	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	bool						mTraceFun;			///< If a trace has been selected for the inner function computeLikelihood()
	const std::vector<double>&	mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	const std::vector<double>&	mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	double						mAbsoluteError;		///< The absolute error at which the computation stops
	unsigned int				mVerbose;			///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;		///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;			///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
	int							mMaxIterations;		///< Maximum number of iterations for the maximization
	int							mNumTimes;			///< Number of branches in the optimizer
};

#endif // OPTSQP_H
