
#ifndef BOOTSTRAP_RANDOM_H
#define BOOTSTRAP_RANDOM_H

#include <vector>
#include "BranchSiteModel.h"

// boost random generation
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
typedef boost::random::mt19937 RNGType;


// uncomment this to use the genetic algorithm bootstrap
#define BOOTSTRAP_GA

#ifndef BOOTSTRAP_GA
// uncomment this to use the PSO bootstrap
#define BOOTSTRAP_PSO
#endif //BOOTSTRAP_GA

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
	/// @param[in] aInitStatus			initialization status (what has been initialized and how) (see InitVarStatus enum type in BranchSiteModel.h)
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
				   ,unsigned int aSeed
				   ,unsigned int aInitStatus) 
		:mUnifRandNumGenerator(aSeed)
		,mGammaDistT(0.5031126, 0.1844347)
		,mExpDistV0(9.441686)
		,mGammaDistV1(0.8764469, 0.126567)
		,mBetaDistW0(1.638631, 21.841174)
		,mGammaDistK(7.547445, 0.5789037)
		,mGammaDistW2(0.209741, 0.5)
		,mN(0)
		,mSizeVect(0)
		,mModel(aModel)
		,mTrace(aTrace)
		,mTraceFun(aTrace)
		,mLowerBound(aLowerBound)
		,mUpperBound(aUpperBound)
		,mVerbose(aVerbose)
		,mStopIfBigger(aStopIfBigger)
		,mThreshold(-aThreshold)
		,mMaxIterations(aMaxIterations)
		,mNumTimes(aNumTimes)
		,mInitStatus(aInitStatus)
		,mIndexBegin(0)
		,mIndexEnd(0)
#ifdef BOOTSTRAP_GA
		,mPopSize(0)
		,mPopPos(NULL)
		,mPopFitness(NULL)
#endif //BOOTSTRAP_GA
#ifdef BOOTSTRAP_PSO
		,mPopSize(0)
		,mWorkSpace(NULL)
		,mPositions(NULL)
		,mFitnesses(NULL)
		,mBestPositions(NULL)
		,mBestFitnesses(NULL)
		,mVelocities(NULL)
#endif // BOOTSTRAP_PSO
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
	/// @param[in] aX the point at which we evaluate the function
	///
	/// @return the function value
	///
	double evaluateLikelihood(const double *aX);
	
	/// evaluateLikelihood
	/// evaluate the likelihood value of the model at point x (std::vector version)
	///
	/// @param[in] aX the point at which we evaluate the function
	///
	/// @return the function value
	///
	double evaluateLikelihood(const std::vector<double> &aX);


private:
	
	/// allocateMemory
	/// allocate the needed space to store the current best variables and the workspace
	///
	void allocateMemory(void);
	
	
	/// generateRandom
	/// generates a random number according to the distribution of the ith variable
	/// distributions are based on data computed previously
	///
	/// @param[in] aIndexVariable the index of the variable to generate
	///
	/// @return the random value
	///
	double generateRandom(unsigned int aIndexVariable);
	
	
#ifdef BOOTSTRAP_GA

	/// bootstrapGeneticAlgorithm
	///	use genetic algorithm to bootstrap the optimization, i.e. help to find a good starting point
	/// 
	/// @param[out] aF The function value at x(out)
	/// @param[out] aX The variables to be optimized; only output, will be initialized in the routine
	/// @param[in] aMaxNumGenerations maximum number of generations
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void bootstrapGeneticAlgorithm(double *aF, double *aX, int aMaxNumGenerations);
	
private:
	
	int							mPopSize;			///< Population size of the ES algorithm
	std::vector<double> 		mGASpace;			///< space used to store the informations for the ES algorithm
	double*						mPopPos;			///< array containing the mPopSize positions of the individuals of the population 
	double*						mPopFitness;		///< Values of the fitness function (likelihood) of each individual
	
#endif // BOOTSTRAP_GA
#ifdef BOOTSTRAP_PSO

	/// bootstrapParticlSwarm
	///	use particle swarm optimization to  bootstrap the optimization, i.e. help to find a good starting point
	/// 
	/// @param[out] aF The function value at x(out)
	/// @param[out] aX The variables to be optimized; only output, will be initialized in the routine
	/// @param[in] aMaxNumGenerations maximum number of generations
	///
	///	@exception FastCodeMLEarlyStopLRT If the optimization has been stopped in advance because LRT is not satisfied
	///
	void bootstrapParticlSwarm(double *aF, double *aX, int aMaxNumGenerations);
	
	private:
	
	int						mPopSize;				///< Population size of the ES algorithm
	std::vector<double> 	mPSOSpace;				///< workspace + space to store positions/velocities of the particles
	double*					mWorkSpace;				///< work space
	double*					mPositions;				///< Positions of the particles
	double*					mFitnesses;				///< corresponding log-likelihood values
	double*					mBestPositions;			///< best position of each particle
	double*					mBestFitnesses;			///< corresponding log-likelihood values
	double*					mVelocities;			///< velocity of each particle

#endif // BOOTSTRAP_PSO
	
	
private:
	
	RNGType 										mUnifRandNumGenerator;	///< Uniform random number generator
	
	boost::random::gamma_distribution<double> 		mGammaDistT;	///< distribution of the branchLengths(mixture)
    boost::random::exponential_distribution<double> mExpDistV0;		///< distribution of 1 - v0
    boost::random::gamma_distribution<double> 		mGammaDistV1;	///< distribution of 1 - v1
	boost::random::beta_distribution<double> 		mBetaDistW0;	///< distribution of w0
    boost::random::gamma_distribution<double> 		mGammaDistK;	///< distribution of kappa
    boost::random::gamma_distribution<double> 		mGammaDistW2;	///< distribution of w2 - 1
	
	
private:
		
	int							mN;					///< Number of unknown parameters
	size_t						mSizeVect;			///< Size in memory of a mN vector
	
	std::vector<double>			mVarsCopy;			///< Copy of the variables used to evaluate the function
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
	
	unsigned int				mInitStatus;		///< Which variables have been initialized and how
	unsigned int				mIndexBegin;		///< index of the first variable to bootstrap
	unsigned int				mIndexEnd;			///< index of the last variable to bootstrap
};

#endif // BOOTSTRAP_RANDOM_H
