
#ifndef CDOSOPTIMIZER_H
#define CDOSOPTIMIZER_H

#define USE_CODE_ML_LINE_SEARCH

#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"

/// CDOSOptimizer class.
/// derivative free optimizer taken from 
/// 'Universal derivative-free optimization method with
///  quadratic convergence', Sergey N. Moiseev, 2011 
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-02-23 (initial version)
///     @version 1.1
///

class CDOSOptimizer
{
public:
	/// Constructor
	/// 
	/// @param[in] aModel				The pointer to the hypothesis class that will be used
	/// @param[in] aTrace				Trace or not the optimizer progress
	/// @param[in] aVerbose				The verbose level
	/// @param[in] aLowerBound			Lower limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aUpperBound			Upper limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aRelativeError		Relative error to stop computation
	/// @param[in] aStopIfBigger		If true stop computation as soon as value is over aThreshold
	/// @param[in] aThreshold			The threshold at which the maximization should be stopped
	/// @param[in] aMaxIterations		Maximum number of iterations for the maximization
	/// 
	CDOSOptimizer(BranchSiteModel* aModel
				 ,bool aTrace
				 ,unsigned int aVerbose
				 ,const std::vector<double>& aLowerBound
				 ,const std::vector<double>& aUpperBound
				 ,double aRelativeError
				 ,bool aStopIfBigger
				 ,double aThreshold
				 ,int aMaxIterations) 
		:mModel(aModel)
		,mTrace(aTrace)
		,mTraceFun(aTrace)
		,mLowerBound(aLowerBound)
		,mUpperBound(aUpperBound)
		,mRelativeError(aRelativeError)
		,mVerbose(aVerbose)
		,mStopIfBigger(aStopIfBigger)
		,mThreshold(-aThreshold)
		,mMaxIterations(aMaxIterations)
		,mN(0)
		,mLambda(0.01)
		{}
	
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars);
	
private:

	/// CDOS minimizer
	/// 
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	///
	/// @return Optimization status (-1 check convergence; 0 success; 1 fail)
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	int CDOSminimizer(double *f, double *x);
	
	
	
	/// PerformLineSearch
	/// Perform a line search in the p direction
	/// 
	/// @params[in,out] y The variable to be optimized
	/// @params[in] p The direction of search
	/// @params[in] step The (initial) step value, copy passing 
	///					so it is not modified
	/// @params[out] fy the value of the function
	///
	void PerformLineSearch(double *y, double *p, double step, double *fy);
	
	
#ifndef USE_CODE_ML_LINE_SEARCH
	/// LineSearch
	/// Perform a line search in the p direction
	/// 
	/// @params[in,out] y The variable to be optimized
	/// @params[in] p The direction of search
	/// @params[in] step The (initial) step value, copy passing 
	///					so it is not modified
	/// @params[out] fy the value of the function
	///
	void LineSearch(double *y, double *p, double step, double *fy);
	
	/// isFeasiblePoint
	/// check if the point satisfies the constraints
	///
	/// @params[in] x the point we consider
	/// @return true if the point satisfies the constrains, false otherwise
	/// 
	bool isFeasiblePoint(double *x) const;
	
	/// reduce_step
	///	reduces the step value when needed, change its size if numIter is too high
	///
	/// @params[in,out] step The step size to reduce or change size
	/// @params[in] numIter The number of iterations we performed
	///
	void reduce_step(double& step, int const& numIter) const;
#endif
	
	/// QRdecomposition
	/// perform a QR decomposition of the first i vectors u.
	/// store the result in qi (last orthonormal vector)
	///
	/// @params[in] width The number of columns in U to be decomposed
	///
	void QRdecomposition(int width);
	
	/// initSearchDirections
	/// initialize the vectors 'u' in each direction as the canonical form
	/// i.e. u1 = (1,0,0,...,0)^T
	///		 u2 = (0,1,0,...,0)^T 
	///		 ...
	///
	void initSearchDirections();
	
	/// alocateWorkspace
	/// alocate the workspace for the QR decomposition
	/// 
	void alocateWorkspace();	
	
	/// computeGradient
	/// Compute the gradient of the function to minimize.
	/// Only used to bootstrap the optimizer
	///
	/// @params[in] f0 The value of the function at the point x0
	/// @params[in] x0 The point at which we compute the gradient
	/// @params[out] g The value of the gradient
	///
	void computeGradient(double f0, double *x0, double *g);
	
#ifdef USE_CODE_ML_LINE_SEARCH
	/// Compute the function moving along p starting from x0 by a percentage t.
	///
	/// @param[in] t Percentage move along p
	/// @param[in] x0 Starting point
	/// @param[in] p Search line vector
	/// @param[out] x The position on which the function should be evaluated
	/// @param[in] n Number of coordinates
	///
	/// @return The function value computed at point x
	///
	double fun_LineSearch(double t, const double x0[], const double p[], double x[], int n);

	/// Linear search using quadratic interpolation from x0[] in the direction of p[].
	/// The formula used is:
    ///                x = x0 + a*p        a ~(0,limit)
	///
	/// Adapted from: Wolfe M. A.  1978.  Numerical methods for unconstrained
    /// optimization: An introduction.  Van Nostrand Reinhold Company, New York. pp. 62-73.
    ///
	/// @param[in,out] f Contains f(x0) for input and f(x) for output
	/// @param[in] x0 Starting point for the search
	/// @param[in] p Search line vector
	/// @param[in] step Is used to find the bracket and is increased or reduced as necessary, and is not terribly important.
	/// @param[in] limit Limit the range of search between 0 and this value
	/// @param[in] e (Unknown)
	/// @param[out] space Workspace
	/// @param[in] iround Iteration number just for reporting
	/// @param[in] n Number of coordinates
	///
	/// @return The value of a as in: x = x0 + a*p  a ~(0,limit)
	///
	double LineSearch2(double *f, const double x0[], const double p[], double step, double limit, double e, double space[], int iround, int n);
#endif
private:
	
	//double (f*)(unsigned, const double* x, double* grad, void* f_data); ///< Function to minimize
	
	int 						mN;					///< Number of unknown parameters
	std::vector<double>			mU;					///< Search directions
	std::vector<double>			mQ;					///< Orthogonal directions (work space for the QR decomposition)
	std::vector<double>			mqi;				///< Last orthogonal vector, used to update the search directions
	std::vector<double>			mSpace;				///< Workspace for the QR decomposition; also used to store the tau variable in QR decomposition
	int							mlwork;				///< Size of the workspace needed for the QR decomposition
	double 						mLambda;			///< Step variable
	
	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	bool						mTraceFun;			///< If a trace has been selected for the inner function computeLikelihood()
	const std::vector<double>&	mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	const std::vector<double>&	mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	double						mRelativeError;		///< The relative error at which the computation stops
	unsigned int				mVerbose;			///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;		///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;			///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
	int							mMaxIterations;		///< Maximum number of iterations for the maximization
	
};

#endif // CDOSOPTIMIZER_H
