
#ifndef OPTALTERNATORSQP_H
#define OPTALTERNATORSQP_H


#include <cstdio>
#include <vector>
#include <memory>
#include "BranchSiteModel.h"
#include "OptSQP.h"


/// OptAlternatorSQP class.
/// sequential quadratic programming optimizer
/// see http://www.neos-guide.org/content/sequential-quadratic-programming
/// alternating between optimization on subspaces of the full space
/// this is motivated by the higher cost of evaluating the likelihood when
/// the matrix parameters (w0, w2, kappa, v1, v2) change
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-30 (initial version)
///     @version 1.1
///

class OptAlternatorSQP
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
	OptAlternatorSQP(BranchSiteModel* aModel
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

	/// SearchSpace
	/// defines the space to use by the function evaluator
	///
	enum SearchSpace
	{
		SPACE_FULL,				///< full space
		SPACE_BRANCHES_ONLY,	///< only branchlengths
		SPACE_EXTRA_ONLY		///< only parameters {v0, v1, w0, k, w2}
	};

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


private:
	
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


#endif // OPTALTERNATORSQP_H
