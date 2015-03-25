
#ifndef BOOTSTRAP_RANDOM_H
#define BOOTSTRAP_RANDOM_H


#include <cstdio>
#include <vector>
#include "BranchSiteModel.h"

// boost random generation
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
typedef boost::random::mt19937 RNGType;


// uncomment this to use the evolution strategy algorithm bootstrap
#define BOOTSTRAP_ES


/// BootstrapType
/// type of bootstrap to use
///
enum BootstrapType
{
	ONLY_RANDOM_TRIES			= 1,	///< tries multiple initializations from distributions, take the best
	RANDOM_TRIES_SEPARATE_VARS	= 2,	///< tries by changing variables one at the time and take the best
	EVOLUTION_STRATEGY			= 3		///< use an evolution strategy (metaheuristics)
};


/// BootstrapRandom class.
/// bootstrap the optimization using the distributions of the variables
///
///     @author Lucas Amoudruz - EPFL.
///     @date 2015-03-23 (initial version)
///     @version 1.1
///

class BootstrapRandom	
{
public:
	/// Constructor
	/// 
	/// @param[in] aModel				The pointer to the hypothesis class that will be used
	/// @param[in] aTrace				Trace or not the optimizer progress
	/// @param[in] aVerbose				The verbose level
	/// @param[in] aLowerBound			Lower limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aUpperBound			Upper limit of the variables to constrain the interval on which the optimum should be computed
	/// @param[in] aStopIfBigger		If true stop computation as soon as value is over aThreshold
	/// @param[in] aThreshold			The threshold at which the maximization should be stopped
	/// @param[in] aMaxIterations		Maximum number of iterations for the maximization
	/// @param[in] aNumTimes			Number of branch lengths to optimize. Set 0 if fixed branch lenghts.
	/// @param[in] aSeed				seed for the random generation
	///
	BootstrapRandom(BranchSiteModel* aModel
				   ,bool aTrace
				   ,unsigned int aVerbose
				   ,const std::vector<double>& aLowerBound
				   ,const std::vector<double>& aUpperBound
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
		,mVerbose(aVerbose)
		,mStopIfBigger(aStopIfBigger)
		,mThreshold(-aThreshold)
		,mMaxIterations(aMaxIterations)
		,mNumTimes(aNumTimes)
		,rng(aSeed)
		,gamma_dist_T(0.5031126, 0.1844347)
		,normal_dist_T(0.058, 0.03)
		//,exp_dist_p0(4.605203)
		//,exp_dist_p1(5.807218)
		,exp_dist_v0(9.441686)
		,gamma_dist_v1(0.8764469, 0.126567)
		,beta_dist_w0(1.638631, 21.841174)
		,gamma_dist_k(7.547445, 0.5789037)
		,gamma_dist_w2(0.209740957, 274.537247)
		,mN(0)
		{}
	
	/// Compute the maximum of computeLikelihood() over several tries
	///
	/// @param[in,out] aVars The variables to be optimized
	///
	/// @return The maximum loglikelihood value over the tries
	///
	double bootstrap(std::vector<double>& aVars);
	

private:
	
	/// evaluateLikelihood
	/// evaluate the likelihood value of the model at point x
	///
	/// @param[in] x the point at which we evaluate the function
	///
	/// @return the function value
	///
	double evaluateLikelihood(double *x);
	double evaluateLikelihood(const std::vector<double> &x);


private:
	
	/// alocateMemory
	/// alocate the needed space to store the current best variables and the workspace
	///
	void alocateMemory(void);
	
	
	/// generateRandom
	/// generates a random number according to the distribution of the ith variable
	/// distributions are based on data computed previously
	///
	/// @param[in] i the index of the variable to generate
	///
	/// @return the random value
	///
	double generateRandom(unsigned int i);
	
	
	/// bootstrapRandomly
	/// performs a bootstrap algorithm using the distribution of the variables x 
	/// pick up random points in the space according the distributions of variables
	/// 
	/// @param[out] f The function value at x(out)
	/// @param[in,out] x The variables to be optimized
	/// @param[in] numTries The number of tests we perform
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void bootstrapRandomly(double *f, double *x, unsigned int numTries);	
	
	
	/// bootstrapEachDirectionRandomly
	/// performs a bootstrap algorithm using the distribution of the variables x 
	/// pick up random points in each direction according the distributions of variable
	/// 
	/// @param[out] f The function value at x(out)
	/// @param[in,out] x The variables to be optimized
	/// @param[in] numGlobal The number of numTries in bootstrapRandomly called at the beginning
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void bootstrapEachDirectionRandomly(double *f, double *x, int numGlobal);	
	
	
	
	
#ifdef BOOTSTRAP_ES

	/// bootstrapEvolutionStrategy
	///	use evolution strategy to bootstrap the optimization, i.e. help to find a good starting point
	/// 
	/// @param[out] f The function value at x(out)
	/// @param[out] x The variables to be optimized; only output, will be initialized in the routine
	/// @param[in] maxNumGenerations maximum number of generations
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void bootstrapEvolutionStrategy(double *f, double *x, int maxNumGenerations);
	
private:
	
	int							mPopSize;			///< Population size of the ES algorithm
	std::vector<double> 		mGASpace;			///< space used to store the informations for the ES algorithm
	double*						mPopPos;			///< array containing the mPopSize positions of the individuals of the population 
	double*						mPopFitness;		///< Values of the fitness function (likelihood) of each individual
	
#endif // BOOTSTRAP_ES
	
	
private:
	
	RNGType 										rng;			///< Uniform andom number generator
	
	boost::random::gamma_distribution<double> 		gamma_dist_T;	///< distribution of the branchLengths(mixture)
	boost::random::normal_distribution<double>		normal_dist_T;	///< distribution of the branchLengths(mixture)
    boost::random::exponential_distribution<double> exp_dist_v0;	///< distribution of 1 - v0
    boost::random::gamma_distribution<double> 		gamma_dist_v1;	///< distribution of 1 - v1
	boost::random::beta_distribution<double> 		beta_dist_w0;	///< distribution of w0
    boost::random::gamma_distribution<double> 		gamma_dist_k;	///< distribution of kappa
    boost::random::gamma_distribution<double> 		gamma_dist_w2;	///< distribution of w2 - 1
	
	
private:
		
	int							mN;					///< Number of unknown parameters
	size_t						size_vect;			///< Size in memory of a mN vector
	
	std::vector<double>			x_;					///< Copy of the variables used to evaluate the function
	std::vector<double>			mSpace;				///< work space

	BranchSiteModel*			mModel;				///< The model for which the optimization should be computed
	bool						mTrace;				///< If a trace has been selected
	bool						mTraceFun;			///< If a trace has been selected for the inner function computeLikelihood()
	const std::vector<double>&	mLowerBound;		///< Lower limit of the variables to constrain the interval on which the optimum should be computed
	const std::vector<double>&	mUpperBound;		///< Upper limit of the variables to constrain the interval on which the optimum should be computed

	unsigned int				mVerbose;			///< The verbose flag from the BranchSiteModel class
	bool						mStopIfBigger;		///< When true stop if lnL is bigger than mThreshold
	double						mThreshold;			///< Threshold for the early stop of optimization if LRT non satisfied (the value is stored with sign changed)
	unsigned int				mMaxIterations;		///< Maximum number of iterations for the maximization
	unsigned int				mNumTimes;			///< Number of branches in the optimizer
};

#endif // BOOTSTRAP_RANDOM_H
