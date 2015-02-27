
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
		,mStep(0)
		,mOmega(1.)
		,s1(3)
		,s2(3)
		{}
	
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars);
	
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
	/// TODO: throw exceptions
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
	
	/// saveDirection
	/// save the last direction dir and reorder the pointers
	///
	/// @params[in] dir pointer on the new direction
	///
	void saveDirection(double *dir);
	
	/// saveGradient
	/// save the last gradient grad and reorder the pointers
	///
	/// @params[in] grad pointer on the new grad
	///
	void saveGradient(double *grad);
	
	/// updateDMatrix
	/// Updates the D matrix by calling all the required subroutines
	/// Scales the columns of the matrix
	///
	void updateDMatrix();
	
	
	/// computeGradient
	/// computes the gradient of the log likelihood function at point x
	/// 
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values computed
	///
	void computeGradient(double aPointValue, const double *aVars, double* aGrad) const;
	
private:
	
	//double (f*)(unsigned, const double* x, double* grad, void* f_data); ///< Function to minimize
	
	int 						mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	int 						mSubspaceStorage;	///< storage size of the subspace directions
	std::vector<double>			mSpace;				///< Workspace
	
	int							mStep;				///< current step	
	int 						s1;					///< number of previous directions to store
	int 						s2;					///< number of previous gradients to store
	double*						mGradient;			///< current gradient
	std::vector<double*>		mGradient_prev;		///< previous gradients
	std::vector<double*>		mDirection_prev;	///< previous directions
	double*						md1;				///< Nemirovski direction 1 (=xk - x0)
	double*						md2;				///< Nemirovski direction 2 (=\sum_{i=1}^k omega_i gradient(x_i))
	double						mOmega;				///< used to compute Nemirovski direction 2
	
	int 						mM;					///< size of the current subspace
	std::vector<double>			alpha;				///< step to optimize in the current subspace
	double						*mD;				///< D matrix, represent the current subspace. Only a pointer on Workspace mSpace
	std::vector<double>			mLowerBoundSubspace;///< Lower bounds for the subspace optimization
	std::vector<double>			mUpperBoundSubspace;///< Upper bounds for the subspace optimization
	
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

#endif // OPTSESOP_H
