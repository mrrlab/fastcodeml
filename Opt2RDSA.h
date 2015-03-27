
#ifndef OPT2RDSA_H
#define OPT2RDSA_H


#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"


/// Opt2RDSA class.
/// 2RDSA optimizer, stochastic gradient based
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-10 (initial version)
///     @version 1.1
///


class Opt2RDSA
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
	/// @param[in] aNumTimes			Number of branch lengths to optimize
	///
	Opt2RDSA(BranchSiteModel* aModel
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
		,mN(0)
		,mStep(0)
		,eta(0.25)
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

	/// Opt2RDSAminimizer
	/// performs the 2RDSA algorithm to estimate the minimum of the function
	/// 
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	///
	/// @return Optimization status (-1 check convergence; 0 success; 1 fail)
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	int Opt2RDSAminimizer(double *f, double *x);
	
	
	/// alocateMemory
	/// alocate all the needed space to store the gradient, hessian matrix and workspace
	///
	void alocateMemory();
	
	
	/// updateParameters
	///	update the parameters aN, deltaN as a function of the step mStep
	///
	void updateParameters();
	
	
	/// computeGradient
	/// computes estimators of the gradient and the hessian matrix of the log likelihood function at point x
	/// 
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aVars Variables to be optimized
	///
	void computeGradientAndHessian(double aPointValue, const double *aVars);
	
	
	/// performOneStage
	/// update the variable x and its value f using a second order stochastic algorithm
	/// 
	/// @param[in,out] f The value of the function at point x
	/// @param[in,out] x Variables to be optimized
	///	
	void performOneStage(double *f, double *x);
	
	
	/// evaluateFunction
	/// evaluate the likelihood at point x
	///
	/// @param[in] x The point at which we want to evaluate the function
	/// @param[in] trace The trace
	///
	/// @return the function value
	///
	double evaluateFunction(double *x, bool trace);
	
	/// evaluateFunction
	/// evaluate the likelihood at point x
	///
	/// @param[in] x The point at which we want to evaluate the function
	/// @param[in] trace The trace
	///
	/// @return the function value
	///
	double evaluateFunction(std::vector<double> &x, bool trace);
	
private:
	
	//double (f*)(unsigned, const double* x, double* grad, void* f_data); ///< Function to minimize
	
	int 						mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	std::vector<double>			mSpace;				///< Workspace
	std::vector<int>			IPIV;				///< Workspace used by lapack	
	std::vector<double> 		x_;					///< Workspace for function evaluation
	std::vector<double>			x_unscaled;			///< Variable unscaled, used to evaluate the function
	int							packSize;			///< Number of elements in the packed matrix
	
	int							mStep;				///< current step	
	double*						mGradient;			///< current gradient estimate
	double*						Hn;					///< current hessian estimate
	
	double						deltaN;				///< Parameter for the gradient approximation
	double						aN;					///< parameter for the update at each iteration
	double						eta;				///< parameter for the random number generation 
	

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
	int 						mNumTimes;			///< Number of branch lengths to optimize
};

#endif // OPTSESOP_H
