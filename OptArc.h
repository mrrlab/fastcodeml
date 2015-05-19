
#ifndef OPTARC_H
#define OPTARC_H

#include <vector>
#include "BranchSiteModel.h"


/// OptArc class.
/// uses arc search for the SR1 hessian updtate
/// see http://web.stanford.edu/group/SOL/dissertations/nwh-thesis-prepress.pdf
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-05-06 (initial version)
///     @version 1.1
///

class OptArc
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
	OptArc(BranchSiteModel* aModel
		  ,bool aTrace
		  ,unsigned int aVerbose
		  ,const std::vector<double>& aLowerBound
		  ,const std::vector<double>& aUpperBound
		  ,double aAbsoluteError
		  ,bool aStopIfBigger
		  ,double aThreshold
		  ,int aMaxIterations
		  ,int aNumTimes) 
		:mN(0)
		,mNs(2)
		,mStep(0)
		,mModel(aModel)
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
	/// alocate the space for storage
	///
	void allocateMemory(void);


	/// ArcMinimizer
	/// Optimization using arc search with SR1 matrix update for hessian approximation
	/// 
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void ArcMinimizer(double *f, double *x);
	
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
	
	/// evaluateFunctionForArcSearch
	///	evaluates the function at point x + Gamma(alpha)
	///
	/// @param[in]	x the point x
	/// @param[in]	alpha the step length
	///
	/// @return the function value
	///
	/// @exception FastCodeMLEarlyStopLRT To force halt the maximization because LRT is already not satisfied
	///
	double evaluateFunctionForArcSearch(const double* x, double alpha);
	
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
	
	/// computeSubspaceArcSearch
	/// computes the subspace for arc search and the orthonormal matrix Q
	///
	void computeSubspaceArcSearch(const double *x);
	
	/// arcSearch
	/// perform an arc search along the arc calculated as in
	/// http://web.stanford.edu/group/SOL/dissertations/nwh-thesis-prepress.pdf
	///
	/// @param[in,out] aalpha in: The initial step
	///						  out: The step chosen by the arc search
	/// @param[in,out] x	  in: The current position
	///						  out: The resulting position, i.e. x+Gamma(aalpha), Gamma being the arc
	/// @param[in,out] f	  in: The function value at x
	///						  out: The function value after the arc search
	///
	void arcSearch(double *aalpha, double *x, double *f);
	
#if 0
	/// zoom
	/// used in the linesearch function to "zoom" in an interval [alo, ahi]
	///
	/// see http://djvuru.512.com1.ru:8073/WWW/e7e02357929ed3ac5afcd17cac4f44de.pdf, 
	/// chap3 pp.59-60 for more informations on the line search algorithm
	///
	/// @param[in] alo	The lower bound of the interval
	/// @param[in] ahi	The upper bound of the interval
	/// @param[in] x	The previous position
	/// @param[in] phi_0		The value phi(0) = f(x + 0.mP)
	/// @param[in] phi_0_prime	The derivative of phi (with phi(a) = f(x+a.mP) at point a=0.
	/// @param[in] phi_lo		The value of the function at point alo
	/// @param[in] c1	The first wolfe variable
	/// @param[in] c2	The second wolfe variable
	///
	/// @return	The (approximate) optimal value of a in the interval [alo, ahi]
	///
	double zoom(double alo, double ahi, double *x, const double& phi_0, const double& phi_0_prime, const double& phi_lo, const double& c1, const double& c2);
#endif 

	/// findMaxStep
	/// finds the maximum step a such that l <= x+Gamma(a) <= u
	///
	/// @param[in] x The current position
	/// @param[out] amax The maximum step size 
	///
	void findMaxStep(const double *x, double *amax);
	
	/// updateActiveSet
	/// updates the active set of position x.
	/// The active set is the set where x is on a bound and the gradient tries to 
	/// "pull" x out of the domain
	///
	/// @param[in] x The current position
	///
	void updateActiveSet(const double *x);
	
	/// projectActiveSet
	/// sets all variables of a vector to 0 if they are in the current active set
	///
	/// @param[in, out] aVect the vector to project
	///
	void projectActiveSet(double *aVect);
	
	/// projectedDirection
	/// computes the projected direction at point x
	///
	/// @param[in]  x The current position
	/// @param[in, out] p in: The direction to project, out: The projected direction
	///
	void projectedDirection(const double *x, double *p);
private:
		
	int 						mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	std::vector<double>			mSpace;				///< Work and storage space
	std::vector<double>			mXEvaluator;		///< Workspace for function evaluations
	
	double*						mGradient;			///< Gradient of the function. mN components
	double*						mHessian;			///< hessian approximation using SR1 update; mN*mN components
	
	double*						mP;					///< search direction
	
	double*						mSk;				///< position change, i.e. mSk = xk - xk-1
	double*						mYk;				///< gradient change, i.e. mYk = dfk - dfk-1
	
	double*						mXPrev;				///< previous position
	double*						mGradPrev;			///< previous gradient
	
	int							mNs;				///< size of subspace S in which we perform the arc search (tipically 2 or 4)
	std::vector<double>			mArcSpace;			///< memory space for the arc search
	double*						mQ;					///< orthogonal matrix of size mN x mNs
	double*						mU;					///< eigenValues matrix of the matrix QT H Q (size mNs x mNs)
	double*						mV;					///< eigenvalues of the matrix QT H Q (size mNs)
	double						mLambdaMin;	
	std::vector<int>			mActiveSet;			///< active set used to reduce the gradient computaion
	
	double*						mWorkSpaceVect;		///< workspace. size of a mN vector.
	double*						mWorkSpaceMat;		///< workspace. size of a mN by mN matrix
	
	int							mStep;				///< current step	
	
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

#endif // OPTARC_H
