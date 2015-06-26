
#ifndef OPTSQPSR1_H
#define OPTSQPSR1_H

#include <vector>
#include <memory>
#include "BranchSiteModel.h"
#include "BOXCQP.h"

// uncomment to rescale the variables before the optimization process
//#define SCALE_OPT_VARIABLES_SR1


/// OptSQPSR1 class.
/// sequential quadratic programming optimizer with an SR1 update instead of the standard 
/// BFGS update. It should be more appropriate for constrained optimization as the hessian
/// approximation can converge to the hessian, while the BFGS update does not guarantee this
/// result. The drawback of SR1 update is a non positive definite matrix.
/// see http://www.meregold.com/abstract.html
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-04-23 (initial version)
///     @version 1.1
///

class OptSQPSR1
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
	OptSQPSR1(BranchSiteModel* aModel
		  	 ,bool aTrace
		  	 ,unsigned int aVerbose
		  	 ,std::vector<double> aLowerBound
		  	 ,std::vector<double> aUpperBound
		  	 ,double aAbsoluteError
			 ,bool aStopIfBigger
			 ,double aThreshold
			 ,int aMaxIterations
			 ,int aNumTimes) 
		:mN(0)
		,mStep(0)
		,mModel(aModel)
		,mTrace(aTrace)
		,mTraceFun(aTrace)
		,mLowerBoundUnscaled(aLowerBound)
		,mUpperBoundUnscaled(aUpperBound)
		,mLowerBound(aLowerBound)
		,mUpperBound(aUpperBound)
		,mAbsoluteError(aAbsoluteError)
		,mVerbose(aVerbose)
		,mStopIfBigger(aStopIfBigger)
		,mThreshold(-aThreshold)
		,mMaxIterations(aMaxIterations)
		,mNumTimes(aNumTimes)
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
	
#ifdef SCALE_OPT_VARIABLES_SR1
	/// scaleVariables
	/// scale the variables of the problem linearly so it is more adapted to the optimization method
	///
	/// @param[in,out] x The variables to be scaled
	/// 
	void scaleVariables(double *x);
	
	/// unscaleVariables
	/// unscale the variables of the problem linearly to recover the "true" values of the variables
	///
	/// @param[in,out] x The variables to be unscaled
	/// 
	void unscaleVariables(double *x);
#endif // SCALE_OPT_VARIABLES_SR1

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
	/// @exception FastCodeMLEarlyStopLRT To force halt the maximization because LRT is already not satisfied
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
	
	/// SR1update
	/// performs the a symmetric rank one (SR1) update to approximate the hessian matrix
	///
	void SR1update(void);
	
	/// activeSetUpdate
	/// updates the activeSet and sets counters so the gradient is not always calculated 
	/// when it is not needed 
	///
	/// @param[in,out] aX The current variables. will be clipped to the boundaries
	/// @param[in] aTolerance The tolerance for a variable to be in the active set
	///
	void activeSetUpdate(double *aX, const double aTolerance);
	
	/// updateFixedVariables
	/// update the set of fixed variables based on the current position
	///
	/// @param[in] aX The current position
	///
	void updateFixedVariables(const double *aX);
	
	/// spectralProjectedGradientIteration
	/// perform an iteration of the spectral projected gradient
	/// This is used to escape from "bad" faces
	///
	/// @param[in,out] aX The position to update
	/// @param[in,out] aF The corresponding function value
	///
	void spectralProjectedGradientIteration(double *aX, double *aF);

	/// computeSearchDirection
	///
	/// compute a search direction (descent direction) before line search
	///
	/// @param[in] aX				The current solution
	/// @param[in] aLocalLowerBound	The lower bound for the search direction
	/// @param[in] aLocalUpperBound	The upper bound for the search direction
	///
	void computeSearchDirection(const double *aX, const double *aLocalLowerBound, const double *aLocalUpperBound);
	
	/// lineSearch
	/// perform a line search in the mP direction satisfying the Armijo condition
	///
	/// @param[in,out] aX The current position (is updated during the function)
	/// @param[in,out] aF The current function value (is updated during the function)
	///
	void lineSearch(double *aX, double *aF);
	
	/// backtracking
	/// perform a backtrace for line search algorithm
	///
	/// @param[in,out]	aX	Current position
	/// @param[in,out]	aF	Current function value
	/// @param[in]		aAlpha	line search position
	/// @param[in]		aPhi0	function value at x
	///	@param[in]		aPhi0_prime	line derivative at alpha=0 	
	///
	void backtrackingLineSearch(double *aX, double *aF, const double aAlpha, const double aPhi0, const double aPhi0_prime);
	
	/// extrapolatingLineSearch
	/// perform an extrapolation for the line search algorithm
	///
	/// @param[in,out]	aX	Current position
	/// @param[in,out]	aF	Current function value
	/// @param[in]		aAlpha		The line position
	/// @param[in]		aAlphaMax	The max line position
	void extrapolatingLineSearch(double *aX, double *aF, const double aAlpha, const double aAlphaMax);

private:
		
	int 						mN;					///< Number of unknown parameters
	size_t						mSizeVect;			///< Size in memory of a mN vector
	
	std::vector<double>			mSpace;				///< Work and storage space
	std::vector<double>			mXEvaluator;		///< Workspace for function evaluations
	
	double*						mGradient;			///< Gradient of the function. mN components
	double*						mProjectedGradient;	///< current projected gradient
	double*						mHessian;			///< hessian approximation using SR1 update; mN*mN components
	
	double*						mP;					///< search direction (convex part)
	
	double*						mSk;				///< position change, i.e. mSk = xk - xk-1
	double*						mYk;				///< gradient change, i.e. mYk = dfk - dfk-1
	
	double*						mXPrev;				///< previous position
	double*						mGradPrev;			///< previous gradient
	
	std::vector<int>			mActiveSet;			///< active set used to reduce the gradient computaion
	std::vector<int>			mFixedVariables;	///< set of fixed variables (on the boundaries)
	
	double*						mWorkSpaceVect;		///< workspace. size of a mN vector.
	double*						mWorkSpaceMat;		///< workspace. size of a mN by mN matrix
	
	int							mStep;				///< current step	
	
	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	bool						mTraceFun;			///< If a trace has been selected for the inner function computeLikelihood()
	
	const std::vector<double>&	mLowerBoundUnscaled;	///< original lower bounds, before scaling	
	const std::vector<double>&	mUpperBoundUnscaled;	///< original upper bounds, before scaling
	
	std::vector<double>			mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	std::vector<double>			mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	
	double						mAbsoluteError;		///< The absolute error at which the computation stops
	unsigned int				mVerbose;			///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;		///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;			///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
	int							mMaxIterations;		///< Maximum number of iterations for the maximization
	int							mNumTimes;			///< Number of branches in the optimizer
};

#endif // OPTSQPSR1_H
