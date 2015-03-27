
#ifndef OPTNES_H
#define OPTNES_H


#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
typedef boost::random::mt19937 RNGType;

/// OptNES class.
/// Natural Evolution Strategy optimizer
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-26 (initial version)
///     @version 1.1
///


class OptNES
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
	OptNES(BranchSiteModel* aModel
		  ,bool aTrace
		  ,unsigned int aVerbose
		  ,const std::vector<double>& aLowerBound
		  ,const std::vector<double>& aUpperBound
		  ,double aAbsoluteError
		  ,bool aStopIfBigger
		  ,double aThreshold
		  ,int aMaxIterations
		  ,int aNumTimes
		  ,unsigned int aSeed) 
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
		,mNumTimes(aNumTimes)
		,mN(0)
		,mStep(0)
		,rng( aSeed )
		{}
	
	/// Compute the maximum of computeLikelihood()
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value
	///
	double maximizeFunction(std::vector<double>& aVars, unsigned int popSize = 10);
	

private:

	/// NESminimizer
	/// performs the natural evolution algorithm to estimate the minimum of the function
	/// 
	/// @param[out] f The minimized function value
	/// @param[in,out] x The variables to be optimized
	///
	/// @return Optimization status (-1 check convergence; 0 success; 1 fail)
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	int NESminimizer(double *f, double *x);
	
	/// alocateMemory
	/// alocate the space for storage
	///
	void alocateMemory(void);
	
	/// initializeDistribution
	/// initialize the parameters mu and sigma of each variable
	///
	void initializeDistribution(void);
	
	/// generatePopulation
	/// generate the mLambda points according to the current distribution
	///
	void generatePopulation(void);
	
	/// computeGradients
	/// compute the gradients of the expected value of the fitness 
	/// gradients w.r.t. mu and sigma
	///
	void computeGradients(void);
	
	/// updateParameters
	/// update the parameters of the distribution
	///
	void updateParameters(void);
	
	/// sortFitness
	/// sort the points mS according to the fitness values
	///
	void sortFitness(void);
	
private:

	/// CompareFitness
	/// utility class used to sort the fitness without changing the vectors.
	/// only used to compare the fitness according to their index
	///
	class CompareFitness
	{
		public:
		/// CompareFitness 
		/// constructor
		/// @param[in] aFitness The fitness reference
		///
		CompareFitness(const double* aFitness)
		:mFitness(aFitness){}
		
		/// operator()
		/// used for sorting according to mFitness
		/// 
		/// @param[in] i the first index value to compare
		/// @param[in] j the second index value to compare
		/// 
		bool operator()(int i, int j) const {return mFitness[i] < mFitness[j];}
		
		private:
		const double* mFitness;
	};
	
private:
	
	//double (f*)(unsigned, const double* x, double* grad, void* f_data); ///< Function to minimize
	
	int 						mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	int							mLambda;			///< population size
	
	std::vector<double>			mSpace;				///< Work and storage space
	std::vector<double> 		x_;					///< Workspace for function evaluation
		
	
	RNGType 									rng;	///< Random generator
	boost::random::normal_distribution<double>	Norm;	///< normal distribution N(0,1) centered in 0 with standard deviation 1
	
	double*						mMu;				///< mean parameters for the distribution
	double*						mSigma;				///< standard deviation parameters for the distribution
	
	double*						mGradMu;			///< gradient of the mean w.r.t. mu parameters
	double*						mGradSigma;			///< gradient of the mean w.r.t. sigma parameters
	
	double*						mPopFitness;		///< fitness values of the population points
	double*						mPopPos;			///< points in the current population
	
	double*						mUtility;			///< utility of the fitness; used to make the algorithm more robust
	double*						mS;					///< random variables (population size) according the current distribution
	
	std::vector<int> 			mPermutation;		///< vector used for sorting (see sortFitness())
	
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

#endif // OPTNES_H
