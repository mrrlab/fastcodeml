
#ifndef OPTSESOP_H
#define OPTSESOP_H


#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"

/// OptSESOP class.
/// SESOP optimizer, gradient based
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-02-27 (initial version)
///     @version 1.1
///


// if use the BOBYQA optimizer with bound constraints
//#define USE_BOUND_CONSTRAINTS

// if use the SLSQP optimizer with linear constraints (maybe more relevant)
#define USE_TRANSFORMED_CONSTRAINTS

class OptSESOP
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
	OptSESOP(BranchSiteModel* aModel
			,bool aTrace
			,unsigned int aVerbose
			,const std::vector<double>& aLowerBound
			,const std::vector<double>& aUpperBound
			,double aRelativeError
			,bool aStopIfBigger
			,double aThreshold
			,int aMaxIterations
			,int aNumTimes) 
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
		,mNumTimes(aNumTimes)
		,mN(0)
		,mStep(0)
		,mOmega(1.)
		{}
	
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars);
	
public:

#ifdef USE_TRANSFORMED_CONSTRAINTS
	// structure only used to store data for the constraints
	struct data_constraint
	{
		OptSESOP*	sesop;		///< Pointer on the instance containing all the informations such as the bounds of the constraints 	
		int 		line;		///< index of the bound
		int 		bound_type; ///< type of bound, 0 meaning lower bound and 1 meaning upper bound
	};
#endif // USE_TRANSFORMED_CONSTRAINTS

private:

	/// SESOPminimizer
	/// performs the SESOP algorithm to estimate the minimum of the function
	/// 
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	///
	/// @return Optimization status (-1 check convergence; 0 success; 1 fail)
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	int SESOPminimizer(double *f, double *x);
	
	
	/// Computes the function.
	/// Note that the computation of the gradient is not yet implemented
	///
	/// @param[in] aVarsAlpha Variables to be optimized
	/// @param[out] aGrad Gradient values
	///
	/// @return The evaluated function
	///
	/// @exception nlopt::forced_stop To force halt the maximization because LRT is already not satisfied
	///
	double operator()(const std::vector<double>& aVarsAlpha, std::vector<double>& aGrad);
	
	
	/// Wrapper to be passed to the nLopt optimizer
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	/// @param[in] aData Opaque pointer containing the function to be passed to the optimizer
	///
	/// @return The evaluated function
	///
	static double subspaceEvaluatorWrapper(const std::vector<double> &x, std::vector<double> &grad, void *data)
    {
    	return (*reinterpret_cast<OptSESOP*>(data))(x, grad);
	}
	
	/// updateOmega
	/// Update the omega weight at stage k
	///
	void updateOmega();
	
	/// alocateMemory
	/// alocate all the needed space to store the directions of search, as well as the workspace
	///
	void alocateMemory();
	
	
	/// updateDMatrix
	/// Updates the D matrix by calling all the required subroutines
	/// Scales the columns of the matrix
	///
	void updateDMatrix();
	
#ifdef USE_TRANSFORMED_CONSTRAINTS
	/// Wrapper to be passed to the nLopt optimizer for the constraints
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values of the constraint
	/// @param[in] aData Opaque pointer containing the function to be passed to the optimizer as well as what constraint we evaluate
	///
	/// @return The evaluated function
	///
	static double myconstraintWrapper(const std::vector<double> &alpha, std::vector<double> &grad, void *data)
	{
		unsigned n = alpha.size(); 
		data_constraint *data_ = (data_constraint*)(data);
		
		return (*reinterpret_cast<OptSESOP*>(data_->sesop))(n, alpha, grad, data);
	}
	
	
	/// operator()
	/// Implements the contraints for the problem reduced at the subspace
	/// constraints have the form x0_i + (D*alpha)_i - u_i <= 0
	/// 						  l_i - x0_i - (D*alpha)_i <= 0
	///
	/// @params[in] alpha The subspace variable to optimize
	/// @params[out] grad The gradient of the constraint
	/// @params[in] data The line i and the upper(1)/lower bound(0) b 
	///				(put a int* containing two ints i and b)
	///
	double operator()(unsigned n, const std::vector<double> &alpha, std::vector<double> &grad, void *data); 
	
#endif // if USE_TRANSFORMED_CONSTRAINTS

#ifdef USE_BOUND_CONSTRAINTS
	/// updateBoundsAndAlpha
	/// set the lower and upper bounds for alpha and initialize alpha
	///
	void updateBoundsAndAlpha();
	
	/// subspaceLowerBoundIsInSpace
	/// verify if the lower bound of the subspace is contained in the full space
	///
	/// @return True if contained in the space
	///			False otherwise
	///
	bool subspaceLowerBoundIsInSpace();
	
	/// subspaceUpperBoundIsInSpace
	/// verify if the upper bound of the subspace is contained in the full space
	///
	/// @return True if contained in the space
	///			False otherwise
	///
	bool subspaceUpperBoundIsInSpace();
#endif // if USE_BOUND_CONSTRAINTS

	/// computeGradient
	/// computes the gradient of the log likelihood function at point x
	/// 
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values computed
	///
	void computeGradient(double aPointValue, const double *aVars, double* aGrad) const;
	
	
	/// computeGradientSubspace
	/// computes the gradient of the log likelihood function at point alpha in the subspace
	/// 
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aAlpha Variables to be optimized
	/// @param[out] aGrad Gradient values computed
	///
	void computeGradientSubspace(double aPointValue, const std::vector<double>& aAlpha, std::vector<double>& aGrad);
	
	
private:
	
	//double (f*)(unsigned, const double* x, double* grad, void* f_data); ///< Function to minimize
	
	int 						mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	int 						mSubspaceStorage;	///< storage size of the subspace directions
	std::vector<double>			mSpace;				///< Workspace
	std::vector<double> 		x_;					///< Workspace for function evaluation
	
	int							mStep;				///< current step	
	double*						mGradient;			///< current gradient
	double*						mGradient_times;	///< search direction with only the components other than w0, w2 and kappa
	double*						md1;				///< Nemirovski direction 1 (=xk - x0)
	double*						md2;				///< Nemirovski direction 2 (=\sum_{i=1}^k omega_i gradient(x_i))
	double						mOmega;				///< used to compute Nemirovski direction 2
	
	int 						mM;					///< size of the current subspace
	std::vector<double>			alpha;				///< step to optimize in the current subspace
	double						*mD;				///< D matrix, represent the current subspace. Only a pointer on Workspace mSpace
#ifdef USE_TRANSFORMED_CONSTRAINTS
	std::vector<data_constraint>			data_constraints;	///< Data used to compute the constraints
#endif // USE_TRANSFORMED_CONSTRAINTS
#ifdef USE_BOUND_CONSTRAINTS
	std::vector<double>			mLowerBoundSubspace;///< Lower bounds for the subspace optimization
	std::vector<double>			mUpperBoundSubspace;///< Upper bounds for the subspace optimization
#endif // USE_BOUND_CONSTRAINTS

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
	int							mNumTimes;			///< Number of branches in the optimizer
};

#endif // OPTSESOP_H
