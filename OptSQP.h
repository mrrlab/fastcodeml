
#ifndef OPTSQP_H
#define OPTSQP_H

#include <vector>
#include <memory>
#include "BranchSiteModel.h"
#include "BOXCQP.h"


/// OptSQP class.
/// sequential quadratic programming optimizer
/// see http://djvuru.512.com1.ru:8073/WWW/e7e02357929ed3ac5afcd17cac4f44de.pdf
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
	/// @param[in] aModel					The pointer to the hypothesis class that will be used
	/// @param[in] aTrace					Trace or not the optimizer progress
	/// @param[in] aVerbose					The verbose level
	/// @param[in] aLowerBound				Lower limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aUpperBound				Upper limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aAbsoluteError			Absolute error to stop computation
	/// @param[in] aHimmelblauTermination	Use Himmelblau stopping criterion
	/// @param[in] aStopIfBigger			If true stop computation as soon as value is over aThreshold
	/// @param[in] aThreshold				The threshold at which the maximization should be stopped
	/// @param[in] aMaxIterations			Maximum number of iterations for the maximization
	/// @param[in] aNumTimes				Number of branches
	/// 
	OptSQP(BranchSiteModel* aModel
		  ,bool aTrace
		  ,unsigned int aVerbose
		  ,const std::vector<double>& aLowerBound
		  ,const std::vector<double>& aUpperBound
		  ,double aAbsoluteError
		  ,bool aHimmelblauTermination
		  ,bool aStopIfBigger
		  ,double aThreshold
		  ,int aMaxIterations
		  ,int aNumTimes) 
		:mN(0)
		,mStep(0)
		,mModel(aModel)
		,mTrace(aTrace)
		,mTraceFun(aTrace)
		,mLowerBound(aLowerBound)
		,mUpperBound(aUpperBound)
		,mAbsoluteError(aAbsoluteError)
		,mHimmelblauTermination(aHimmelblauTermination)
		,mVerbose(aVerbose)
		,mStopIfBigger(aStopIfBigger)
		,mThreshold(-aThreshold)
		,mMaxIterations(aMaxIterations)
		,mNumTimes(aNumTimes)
		,mSizeVect(0)
		,mGradient(NULL)
		,mHessian(NULL)
		,mP(NULL)
		,mSk(NULL)
		,mYk(NULL)
		,mXPrev(NULL)
		,mGradPrev(NULL)
		,mWorkSpaceVect(NULL)
		,mWorkSpaceMat(NULL)
		{}
	
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars);
	

private:

	/// allocateMemory
	/// allocate the space for storage
	///
	void allocateMemory(void);

	/// SQPminimizer
	/// performs a sequential quadratic approximation of the function to estimate its minimum
	///	uses quadratic programming to solve local constrained subproblems 
	/// 
	/// @param[out] aF The minimized function value
	/// @param[in,out] aX The variables to be optimized
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void SQPminimizer(double *aF, double *aX);
	
	/// evaluateFunction
	///	evaluates the function at point x
	///
	/// @param[in]	aX the point at which one evaluates the function
	/// @param[in]	aTrace The trace 
	///
	/// @return the function value
	///
	/// @exception nlopt::forced_stop To force halt the maximization because LRT is already not satisfied
	///
	double evaluateFunction(const double *aX, bool aTrace);
	
	/// evaluateFunctionForLineSearch
	///	evaluates the function at point x + alpha*mP
	///
	/// @param[in]	aX the point x
	/// @param[in]	aAlpha the step length
	///
	/// @return the function value
	///
	/// @exception FastCodeMLEarlyStopLRT To force halt the maximization because LRT is already not satisfied
	///
	double evaluateFunctionForLineSearch(const double* aX, double aAlpha);
	
	/// computeGradient
	///	compute the gradient at point x using finite differences aproximation
	///
	/// @param[in]	aX the point at which one evaluates the gradient
	/// @param[in]	aF0 The value at point x
	/// @param[out]	aGrad The gradient 
	///
	void computeGradient(const double *aX, double aF0, double *aGrad);
	
	/// hessianInitialization
	///
	/// initialize the hessian matrix
	///
	void hessianInitialization(void);
	
	/// BFGSupdate
	/// performs the BFGS formula to update the hessian matrix.
	/// uses a modified version of BFGS to keep the matrix positive definite
	///
	void BFGSupdate(void);
	
	/// activeSetUpdate
	/// updates the activeSet and sets counters so the gradient is not always calculated 
	/// when it is not needed 
	///
	/// @param[in,out] aX The current variables. will be clipped to the boundaries
	/// @param[in] aTolerance The tolerance for a variable to be in the active set
	///
	void activeSetUpdate(double *aX, const double aTolerance);
	
	/// lineSearch
	/// perform a line search in the mP direction
	///
	/// see http://djvuru.512.com1.ru:8073/WWW/e7e02357929ed3ac5afcd17cac4f44de.pdf, 
	/// chap3 pp.59-60 for more informations on the line search algorithm using strong
	/// wolfe condition.
	///
	/// @param[in,out] aAlpha 	in: initial guess of step length
	///							out: step length
	/// @param[in,out] aX		in: the original position
	///							out: if success, the new position
	/// @param[in,out] aF		in: the original function value
	///							out: if success, the new value
	///
	void lineSearch(double *aAlpha, double *aX, double *aF);
	
	/// zoom
	/// used in the linesearch function to "zoom" in an interval [alo, ahi]
	///
	/// see http://djvuru.512.com1.ru:8073/WWW/e7e02357929ed3ac5afcd17cac4f44de.pdf, 
	/// chap3 pp.59-60 for more informations on the line search algorithm
	///
	/// @param[in] aAlo	The lower bound of the interval
	/// @param[in] aAhi	The upper bound of the interval
	/// @param[in] aX	The previous position
	/// @param[in] aPhi0		The value phi(0) = f(x + 0.mP)
	/// @param[in] aPhi0Prime	The derivative of phi (with phi(a) = f(x+a.mP) at point a=0.
	/// @param[in] aPhiLo		The value of the function at point alo
	/// @param[in] c1	The first wolfe variable
	/// @param[in] c2	The second wolfe variable
	///
	/// @return	The (approximate) optimal value of a in the interval [aAlo, aAhi]
	///
	double zoom(double aAlo, double aAhi, const double *aX, const double& aPhi0, const double& aPhi0Prime, const double& aPhiLo, const double& c1, const double& c2);

private:
	
	int 						mN;						///< Number of unknown parameters
	size_t						mSizeVect;				///< Size in memory of a mN vector
	
	std::vector<double>			mSpace;					///< Work and storage space
	std::vector<double> 		mXEvaluator;			///< Workspace for function evaluations
	
	double*						mGradient;				///< Gradient of the function. mN components
	double*						mHessian;				///< positive definite hessian approximation using BFGS; mN*mN components
	
	double*						mP;						///< search direction

	double*						mSk;					///< position change, i.e. mSk = xk+1 - xk
	double*						mYk;					///< gradient change, i.e. mYk = gk+1 - gk

	double*						mXPrev;					///< previous position
	double*						mGradPrev;				///< previous gradient
	
	std::vector<int>			mActiveSet;				///< active set used to reduce the gradient computaion
	
	double*						mWorkSpaceVect;			///< workspace. size of a mN vector.
	double*						mWorkSpaceMat;			///< workspace. size of a mN by mN matrix
	
	int							mStep;					///< current step	
	
	std::auto_ptr<BOXCQP>		mQPsolver;				///< box constrained quadratic program solver

	BranchSiteModel*			mModel;					///< The model for which the optimization should be computed
	bool						mTrace;					///< If a trace has been selected
	bool						mTraceFun;				///< If a trace has been selected for the inner function computeLikelihood()
	
	std::vector<double>			mLowerBound;			///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	std::vector<double>			mUpperBound;			///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	
	double						mAbsoluteError;			///< The absolute error at which the computation stops
	bool						mHimmelblauTermination;	///< Use a Himmelblau stopping criterion instead of function reduction only
	unsigned int				mVerbose;				///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;			///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;				///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
	int							mMaxIterations;			///< Maximum number of iterations for the maximization
	int							mNumTimes;				///< Number of branches in the optimizer
};

#endif // OPTSQP_H
