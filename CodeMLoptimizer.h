
#ifndef CODEMLOPTIMIZER_H
#define CODEMLOPTIMIZER_H

#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"

/// Minimizer from CodeML.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-01-11 (initial version)
///     @version 1.0
///
class Ming2
{
public:
	/// Constructor
	///
	/// @param[in] aModel				The pointer to the hypothesis class that will be used
	/// @param[in] aTrace				Trace or not the optimizer progress
	/// @param[in] aLowerBound			Lower limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aUpperBound			Upper limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aDeltaForGradient	Delta used to compue the gradient
	///
	Ming2(BranchSiteModel* aModel, bool aTrace, const std::vector<double>& aLowerBound, const std::vector<double>& aUpperBound, double aDeltaForGradient) :
			mModel(aModel),
			mTrace(aTrace),
			mTraceFun(false),
			mLowerBound(aLowerBound),
			mUpperBound(aUpperBound),
			mDeltaForGradient(aDeltaForGradient),
			mFgBranch(0),
			mNoisy(0),
			mAlwaysCenter(false),
			mIround(0),
			mSIZEp(0.) {}

	/// Do the minimization of: aModel->computeLikelihood(x, n, mFgBranch, mTraceFun);
	///
	/// @param[in] aFgBranch The foreground branch number
	/// @param[in,out] aVars The variables that should be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double minimizeFunction(unsigned int aFgBranch, std::vector<double>& aVars);

private:
	/// The original ming2 minimizer.
	/// Few parameters dropped
	///
	/// @param[in] fout File descriptor on which the trace of the variables is written
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	/// @param[in] xl Lower limits for x
	/// @param[in] xu Upper limits for x
	/// @param[out] space Working space
	/// @param[out] ispace Working space (integer)
	/// @param[in] e Error (?)
	/// @param[in] n Number of variables
	///
	/// @return Optimization status (-1 check convergence; 0 success; 1 fail)
	///
	int ming2(FILE *fout, double *f, double x[], const double xl[], const double xu[], double space[], int ispace[], double e, int n);

	/// Compute the gradient at point x
	///
	/// @param[in] n Number of variables
	/// @param[in] x The point on which the gradient should be computed
	/// @param[in] f0 The function value in x
	/// @param[out] g The computed gradient
	/// @param[out] space Workspace
	/// @param[in] xmark 0: central; 1: upper; -1: down
	///
	void gradientB(int n, const double x[], double f0, double g[], double space[], const int xmark[]) const;

	/// Himmelblau termination rule.
	///
	/// @param[in] x0 (Unknown)
	/// @param[in] x1 (Unknown)
	/// @param[in] f0 (Unknown)
	/// @param[in] f1 (Unknown)
	/// @param[in] e1 (Unknown)
	/// @param[in] e2 (Unknown)
	/// @param[in] n Number of variables in x0 and x1
	///
	/// @return true for stop, false otherwise.
	///
	bool H_end(const double x0[], const double x1[], double f0, double f1, double e1, double e2, int n) const;

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
	/// @param[in] n Number of coordinates
	///
	/// @return The value of a as in: x = x0 + a*p  a ~(0,limit)
	///
	double LineSearch2(double *f, const double x0[], const double p[], double step, double limit, double e, double space[], int n);


private:
	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	bool						mTraceFun;			///< If a trace has been selected for the inner function computeLikelihood()
	const std::vector<double>&	mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	const std::vector<double>&	mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed
	double						mDeltaForGradient;	///< This is the original Small_Diff value
	unsigned int				mFgBranch;			///< The selected foreground branch

private:
	/// The following variables are from the original code
	bool						mAlwaysCenter;		///< From the original code
	int							mNoisy;				///< How much rubbish on the screen. Valid values: 0,1,2,3,9
	int							mIround;			///< From the original code
	double						mSIZEp;				///< From the original code
};

#endif
