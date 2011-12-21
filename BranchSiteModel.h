
#ifndef BRANCHSITEMODEL_H
#define BRANCHSITEMODEL_H

#include <cstring>
#include <vector>

#include "TransitionMatrix.h"
#include "TransitionMatrixSet.h"
#include "Forest.h"

/// Common routines for the Hypothesis test.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-23 (initial version)
///     @version 1.0
///
///
class BranchSiteModel
{
public:
	/// Constructor.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aNumBranches Number of tree branches
	/// @param[in] aNumSites Number of sites
	/// @param[in] aSeed Random number generator seed
	/// @param[in] aNumVariables Number of extra variables (k, w0, w2, p0, p1)
	/// @param[in] aOnlyInitialStep Compute only the first step, no optimization involved
	/// @param[in] aTimesFromTree Takes the time from the input tree
	/// @param[in] aTrace If set print a trace of the maximization process
	/// @param[in] aOptAlgo Maximization algorithm to be used
	///
	BranchSiteModel(Forest& aForest,
					size_t aNumBranches,
					size_t aNumSites,
					unsigned int aSeed,
					unsigned int aNumVariables,
					bool aOnlyInitialStep,
					bool aTimesFromTree,
					bool aTrace,
					unsigned int aOptAlgo=0)
		: mForest(aForest),
		  mNumTimes(aNumBranches),
		  mNumVariables(aNumVariables),
		  mMaxLnL(-DBL_MAX),
		  mNumEvaluations(0),
		  mOnlyInitialStep(aOnlyInitialStep),
		  mTimesFromTree(aTimesFromTree),
		  mTrace(aTrace),
		  mOptAlgo(aOptAlgo),
		  mSeed(aSeed)
	{
		mVar.resize(mNumTimes+mNumVariables);
		mLowerBound.reserve(mNumTimes+mNumVariables);
		mUpperBound.reserve(mNumTimes+mNumVariables);
		mLikelihoods.resize(Nt*aNumSites);
	}

	/// Destructor.
	///
	virtual ~BranchSiteModel() {}

	/// Set the times on the tree from the variables
	///
	void saveComputedTimes(void) const
	{
		mForest.setLengthsFromTimes(mVar);
	}

	/// Formatted print of the maximizer variables array
	///
	/// @param[in] aVars The variables array to be printed
	/// @param[in] aLnl The likelihood value to be printed
	///
	void printVar(const std::vector<double>& aVars, double aLnl=DBL_MAX) const;

	/// Compute the maximum likelihood for the given forest
	///
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	///
	/// @return The maximum Likelihood value
	///
	double maximizeLikelihood(size_t aFgBranch);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	virtual double computeLikelihood(unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace) =0;

	/// Get variable values
	///
	/// @return Pointer to the set of optimized variables
	///
	const double* getStartingValues(void) const {return &mVar[0];}

	/// Get variable values
	///
	/// @param[out] aVariables Vector that will be filled with the variables
	///
	void getVariables(std::vector<double>& aVariables) const {aVariables = mVar;}

	enum {
		OPTIM_LD_LBFGS		= 0,	///< Same optimizer as CodeML
		OPTIM_LD_VAR1		= 1,	///< Shifted limited-memory variable-metric rank-1 method
		OPTIM_LD_VAR2		= 2,	///< Shifted limited-memory variable-metric rank-2 method
		OPTIM_LD_TNEWTON	= 3,	///< Preconditioned inexact truncated Newton algorithm with restart

		OPTIM_LN_BOBYQA		= 11,	///< One gradient free optimizer
		OPTIM_LN_COBYLA		= 12,	///< Another gradient free optimizer

		OPTIM_MLSL_LDS		= 99	///< A global optimizer
	};

protected:
	/// Compute the four site proportions from the two values in the optimization variables
	///
	/// @param[in] aV0 The first optimization variables
	/// @param[in] aV1 The second optimization variables
	/// @param[out] aProportions The four proportions output
	///
	void getProportions(double aV0, double aV1, double* aProportions) const
	{
#ifdef USE_ORIGINAL_PROPORTIONS
		aProportions[0] = exp(aV0);
		aProportions[1] = exp(aV1);
		double tot = aProportions[0] + aProportions[1] + 1;
		aProportions[0] /= tot;
		aProportions[1] /= tot;
		tot = aProportions[0] + aProportions[1];

		aProportions[2] = (1. - tot)*aProportions[0]/tot;
		aProportions[3] = (1. - tot)*aProportions[1]/tot;
#else
		aProportions[0] = aV0*aV1;
		aProportions[1] = aV0*(1-aV1);
		aProportions[2] = (1-aV0)*aV1;
		aProportions[3] = (1-aV0)*(1-aV1);
#endif
	}

	/// Check if two values are sufficiently different
	///
	/// @param[in] aFirst First number to compare
	/// @param[in] aSecond Second term to compare
	///
	/// @return True if the two parameters differs more than TOL
	///
	static bool isDifferent(double aFirst, double aSecond)
	{
		static const double TOL = 1e-7;
		const double diff = aFirst - aSecond;
		return (diff > TOL || diff < -TOL);
	}

protected:
	Forest&						mForest;			///< The forest to be used
	unsigned int				mNumTimes;			///< Number of branch lengths
	unsigned int				mNumVariables;		///< The number of extra variables (4 for H0 and 5 for H1)
	std::vector<double>			mVar;				///< Variable to optimize (first the branch lengths then the remaining variables)
	std::vector<double>			mLowerBound;		///< Lower limits for the variables to be optimized
	std::vector<double>			mUpperBound;		///< Upper limits for the variables to be optimized
	double						mProportions[4];	///< The four proportions
	double						mMaxLnL;			///< Maximum value of LnL found during optimization
	unsigned int				mNumEvaluations;	///< Counter of the likelihood function evaluations
	CacheAlignedDoubleVector	mLikelihoods;		///< Computed likelihoods at the root of all trees. Defined here to make it aligned.
	bool						mOnlyInitialStep;	///< Only the initial step is executed, no optimization
	bool						mTimesFromTree;		///< Read the initial times from the tree
	bool						mTrace;				///< Enable maximization tracing
	unsigned int				mOptAlgo;			///< Optimization algorithm to use

private:
	unsigned int				mSeed;				///< Random number generator seed to be passed to the optimizer
};

/// Null Hypothesis test.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-23 (initial version)
///     @version 1.0
///
///
class BranchSiteModelNullHyp : public BranchSiteModel
{
public:
	/// Constructor.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aSeed Random number generator seed
	/// @param[in] aOnlyInitialStep If true no optimization is done, only the initial step is run
	/// @param[in] aTimesFromTree Initial times are from the file plus fixed values for the other variables
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aOptAlgo The optimization algorithm to use
	///
	BranchSiteModelNullHyp(Forest& aForest, unsigned int aSeed, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace, unsigned int aOptAlgo=0)
		: BranchSiteModel(aForest, aForest.getNumBranches(), aForest.getNumSites(), aSeed, 4, aOnlyInitialStep, aTimesFromTree, aTrace, aOptAlgo), mSet(aForest.getNumBranches(), 3), mPrevK(DBL_MAX), mPrevOmega0(DBL_MAX) {}

	/// Compute the null hypothesis log likelihood.
	///
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	///
	/// @return The log likelihood for the null hypothesis
	///
	double operator()(size_t aFgBranch);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	double computeLikelihood(unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace);


private:
	TransitionMatrix 	mQw0;			///< Q matrix for the omega0 case
	TransitionMatrix 	mQ1;			///< Q matrix for the omega1 == 1 case
	TransitionMatrixSet mSet;			///< Set of matrices used for the tree visits
	double				mPrevK;			///< Previous k value used to compute matrices
	double				mPrevOmega0;	///< Previous w0 value used to compute matrices
	double				mScaleQw0;		///< Scale value for Qw0
	double				mScaleQ1;		///< Scale value for Q1
};



/// Alternate Hypothesis test.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-23 (initial version)
///     @version 1.0
///
///
class BranchSiteModelAltHyp : public BranchSiteModel
{
public:
	/// Constructor.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aSeed Random number generator seed
	/// @param[in] aOnlyInitialStep If true no optimization is done, only the initial step is run
	/// @param[in] aTimesFromTree Initial times are from the file plus fixed values for the other variables
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aOptAlgo The optimization algorithm to use
	///
	BranchSiteModelAltHyp(Forest& aForest, unsigned int aSeed, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace, unsigned int aOptAlgo=0)
		: BranchSiteModel(aForest, aForest.getNumBranches(), aForest.getNumSites(), aSeed, 5, aOnlyInitialStep, aTimesFromTree, aTrace, aOptAlgo), mSet(aForest.getNumBranches(), 4), mPrevK(DBL_MAX), mPrevOmega0(DBL_MAX), mPrevOmega2(DBL_MAX) {}

	/// Compute the alternative hypothesis log likelihood.
	///
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aInitFromH0 If not null uses these results from H0 to initalize H1 values
	///
	/// @return The log likelihood for the alternative hypothesis
	///
	double operator()(size_t aFgBranch, const double* aInitFromH0=0);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	double computeLikelihood(unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace);


private:
	TransitionMatrix    mQw0;			///< Q matrix for the omega0 case
	TransitionMatrix    mQw2;			///< Q matrix for the omega2 case
	TransitionMatrix    mQ1;  			///< Q matrix for the omega1 == 1 case
	TransitionMatrixSet mSet;			///< Set of matrices used for the tree visits
	double				mPrevK;			///< Previous k value used to compute matrices
	double				mPrevOmega0;	///< Previous w0 value used to compute matrices
	double				mPrevOmega2;	///< Previous w2 value used to compute matrices
	double				mScaleQw0;		///< Scale value for Qw0
	double				mScaleQw2;		///< Scale value for Qw2
	double				mScaleQ1;		///< Scale value for Q1
};


#endif
