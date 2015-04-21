
#ifndef OPT_TRUST_REGION_H
#define OPT_TRUST_REGION_H

#include <cstdio>
#include <vector>
#include <memory>
#include "BranchSiteModel.h"
#include "BOXCQP.h"


// uncomment to rescale the variables before the optimization process
#define SCALE_OPT_TRUST_REGION_VARIABLES


/// OptTrustRegion class.
/// trust region optimizer
/// see http://www.optimization-online.org/DB_FILE/2010/06/2643.pdf
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-04-21 (initial version)
///     @version 1.1
///

class OptTrustRegion
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
	OptTrustRegion(BranchSiteModel* aModel
				  ,bool aTrace
				  ,unsigned int aVerbose
				  ,std::vector<double> aLowerBound
				  ,std::vector<double> aUpperBound
				  ,double aAbsoluteError
				  ,bool aStopIfBigger
				  ,double aThreshold
				  ,int aMaxIterations
				  ,int aNumTimes) 
		:mModel(aModel)
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
	
#ifdef SCALE_OPT_TRUST_REGION_VARIABLES
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
#endif // SCALE_OPT_TRUST_REGION_VARIABLES

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
	
	/// computeRatio
	///	compute the improvement ratio between the function value and the model
	///
	/// @param[in]	f0	The function value at origin point x0
	/// @param[in]	dx	The step for which we want to measure the improvement
	/// @param[in]	f	The function value at point x0+dx
	///
	/// @return the improvement ratio
	///
	double computeRatio(double f0, double *dx, double f);
	
	/// computeGradient
	///	compute the gradient at point x using finite differences aproximation
	///
	/// @param[in]	x the point at which one evaluates the gradient
	/// @param[in]	f0 The value at point x
	/// @param[out]	aGrad The gradient 
	///
	void computeGradient(const double *x, double f0, double *aGrad);
	
	/// BFGSupdate
	/// performs the BFGS formula to update the hessian matrix.
	/// uses a modified version of BFGS to keep the matrix positive definite
	///
	void BFGSupdate(void);
	

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
	
	std::vector<int>			mActiveSet;			///< active set used to reduce the gradient computaion
	
	double*						mWorkSpaceVect;		///< workspace. size of a mN vector.
	double*						mWorkSpaceMat;		///< workspace. size of a mN by mN matrix
	
	int							mStep;				///< current step	
	int							mNumCallZoom;		///< step of the zoom iteration in line search
	
	std::auto_ptr<BOXCQP>		mQPsolver;			///< box constrained quadratic program solver

	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	bool						mTraceFun;			///< If a trace has been selected for the inner function computeLikelihood()
	
	const std::vector<double>	mLowerBoundUnscaled;	///< original lower bounds, before scaling	
	const std::vector<double>	mUpperBoundUnscaled;	///< original upper bounds, before scaling
	
	std::vector<double>			mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	std::vector<double>			mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	
	double						mAbsoluteError;		///< The absolute error at which the computation stops
	unsigned int				mVerbose;			///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;		///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;			///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
	int							mMaxIterations;		///< Maximum number of iterations for the maximization
	int							mNumTimes;			///< Number of branches in the optimizer
};

#endif // OPT_TRUST_REGION_H
