
#ifndef OPTALTERNATOR_H
#define OPTALTERNATOR_H


#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"



/// OptAlternator class.
/// Optimize by alternating the space of search
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-27 (initial version)
///     @version 1.1
///


class OptAlternator
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
	OptAlternator(BranchSiteModel* aModel
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
		,mSearchState(SPACE_EXTRA_ONLY)
		{}
	
	/// maximizeFunction
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] 	aVars The variables to be optimized
	/// @return 		The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars);
	

private:

	/// AlternatorMinimizer
	/// minimize the function by alternating the search space 
	/// 
	/// @param[out] 	f The minimized function value
	/// @param[in,out]	x The variables to be optimized
	///
	///	@exception 		FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void AlternatorMinimizer(double *f, double *x);
	
	
	
private:

	/// evaluateFunction
	/// evaluate the likelihood in space described by the current mSearchState.
	/// aVar and aGrad should have length of this search space
	///
	/// @param[in]		aVar	The point where we evaluate the likelihood
	/// @param[in,out] 	aGrad	The gradient if it is not set to NULL pointer. 
	///							NULL otherwise. 
	/// @return 				The function value
	///
	double evaluateFunction(const double *aVar, double *aGrad);
	
	/// Wrapper to be passed to the nLopt optimizer
	///
	/// @param[in] 	x 		Variables to be optimized, size of the current space
	/// @param[out] grad 	Gradient values
	/// @param[in] 	data 	Opaque pointer containing the function to be passed to the optimizer
	///
	/// @return 			The evaluated function
	///
	static double subspaceEvaluatorWrapper(const std::vector<double>& x, std::vector<double>& grad, void *data)
    {
    	return (reinterpret_cast<OptAlternator*>(data))->evaluateFunction(&x[0], &grad[0]);
	}
	
private:
	
	/// SearchSpace
	/// defines the space to use by the function evaluator
	///
	enum SearchSpace
	{
		SPACE_FULL,				///< full space
		SPACE_BRANCHES_ONLY,	///< only branchlengths
		SPACE_EXTRA_ONLY		///< only parameters in {v0, v1, w0, k, w2}
	};
	
	/// getSpaceProperties
	/// gives the indices of variables of current space (bounds)
	///
	/// @param[out]	idFirstVar	The index of first variable of the space
	/// @param[out]	idLastVar	The index of last variable of the space
	///
	void getSpaceProperties(int& idFirstVar, int& idLastVar) const;

private:
	
	int 						mN;					///< Number of unknown parameters (total)
	int							mNextra;			///< Number of parameters without branchlengths, can be 4 (hyp0) or 5 (hyp1)
	int							mNumTimes;			///< Number of branches in the optimizer
		
	size_t						size_vect;			///< Size in memory of a mN vector
	
	std::vector<double>			mSpace;				///< Work and storage space
	std::vector<double> 		x_;					///< Workspace for function evaluation
		

	int							mStep;				///< current step 	
	int							mMaxIterations;		///< Maximum number of iterations for the maximization
	SearchSpace					mSearchState;		///< current search space

	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	const std::vector<double>&	mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	const std::vector<double>&	mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	double						mAbsoluteError;		///< The absolute error at which the computation stops
	unsigned int				mVerbose;			///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;		///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;			///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
};

#endif // OPTALTERNATOR_H

