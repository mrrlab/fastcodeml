
#ifndef BRANCHSITEMODEL_H
#define BRANCHSITEMODEL_H

#include <cstring>
#include <vector>
#include <cstdlib>

#include "TransitionMatrix.h"
#include "ProbabilityMatrixSet.h"
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
	/// @param[in] aTrace If set print a trace of the maximization process
	/// @param[in] aOptAlgo Maximization algorithm to be used
	/// @param[in] aDeltaValueForGradient The variable increment to compute gradient
	///
	BranchSiteModel(Forest& aForest,
					size_t aNumBranches,
					size_t aNumSites,
					unsigned int aSeed,
					unsigned int aNumVariables,
					bool aOnlyInitialStep,
					bool aTrace,
					unsigned int aOptAlgo,
					double aDeltaValueForGradient)
		: mForest(aForest),
		  mNumTimes(aNumBranches),
		  mNumVariables(aNumVariables),
		  mMaxLnL(-DBL_MAX),
		  mNumEvaluations(0),
		  mOnlyInitialStep(aOnlyInitialStep),
		  mTrace(aTrace),
		  mOptAlgo(aOptAlgo),
		  mInitType(INIT_TYPE_NONE),
		  mDeltaForGradient((aDeltaValueForGradient > 0.0) ? aDeltaValueForGradient : sqrt(DBL_EPSILON)),
		  mSeed(aSeed)
	{
		mVar.resize(mNumTimes+mNumVariables);
		mLikelihoods.resize(Nt*aNumSites);
		setLimits(mNumTimes, mNumVariables);
	}

	/// Destructor.
	///
	virtual ~BranchSiteModel() {}

	/// Set the times on the tree from the variables
	///
	void saveComputedTimes(void) const {mForest.setLengthsFromTimes(mVar);}

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

	/// Compute the likelihood for the given forest and the given set of parameters.
	///
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	virtual double computeLikelihood(const std::vector<double>& aVar, bool aTrace) =0;

	/// Compute the likelihood for the given forest and the given set of parameters.
	/// This version is for use inside Ming2 minimizer
	///
	/// @param[in] aVar The optimizer variables
	/// @param[in] aVarLen The optimizer variables array length
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	double computeLikelihood(double* aVar, int aVarLen, bool aTrace)
	{
		std::vector<double> x(aVar, aVar+aVarLen);
		return computeLikelihood(x, aTrace);
	}

	/// Get variable values
	///
	/// @param[out] aVariables Vector that will be filled with the variables
	///
	void getVariables(std::vector<double>& aVariables) const {aVariables = mVar;}

	/// Get the number of function evaluation.
	///
	/// @return The number of likelihood function calls
	///
	unsigned int getNumEvaluations(void) const {return mNumEvaluations;}

	/// Perform the Likelihood Ratio Test.
	/// LRT test: -2*(lnl0-lnl1) > chisq(.95, df=1)
	///
	/// @param[in] aLnL0 Max LogLikelihood for H0
	/// @param[in] aLnL1 Max LogLikelihood for H1
	///
	///	@return True if the test is passed
	///
	static bool performLRT(double aLnL0, double aLnL1) {return (aLnL1 - aLnL0) > 1.92072941;}

	/// Valid values (on the command line) for the optimization algorithm
	enum OptimAlgoIdentifier
	{
		OPTIM_LD_LBFGS		= 0,	///< Same optimizer as CodeML
		OPTIM_LD_VAR1		= 1,	///< Shifted limited-memory variable-metric rank-1 method
		OPTIM_LD_VAR2		= 2,	///< Shifted limited-memory variable-metric rank-2 method
		OPTIM_LD_SLSQP		= 3,	///< Sequential quadratic programming (SQP) algorithm

		OPTIM_LN_BOBYQA		= 11,	///< Derivative-free bound-constrained optimization using an iteratively constructed quadratic approximation for the objective function

		OPTIM_LD_MING2		= 22,	///< The optimizer extracted from CodeML

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
		aProportions[1] = aV0*(1.-aV1);
		aProportions[2] = (1.-aV0)*aV1;
		aProportions[3] = (1.-aV0)*(1.-aV1);
#endif
	}

	/// Initialize variables to be optimized.
	/// It uses mInitType to know what has been already initialized by initFromTree() or initFromResult()
	///
	void initVariables(void);

private:
	/// Set upper and lower limits for the maximization domain
	///
	/// @param[in] aNumTimes Number of times (ie branch lengths)
	/// @param[in] aNumVariables Number of other variables (4 for H0, 5 for H1)
	///
	void setLimits(unsigned int aNumTimes, unsigned int aNumVariables);

	/// Generate a double random number between 0 and 1
	///
	/// @return The random number
	///
	static inline double randFrom0to1(void) {return static_cast<double>(rand())/static_cast<double>(RAND_MAX);}


public:
	/// Initialize the times from the input phylogenetic tree
	///
	void initFromTree(void);

	/// Initialize the times from the input phylogenetic tree and set the other values to hardcoded constants.
	/// This routine could be used in place of initFromTree()
	///
	void initFromTreeAndFixed(void);

	/// Initialize the times from the input phylogenetic tree and set the other values to hardcoded constants with P0=1 and P1=0.
	/// This routine could be used in place of initFromTree()
	///
	void initFromTreeAndFixedP0(void);

	/// Initialize variables from a previous optimization result
	///
	/// @param[in] aPreviousResult A previous result from another model (obtained by getVariables())
	/// @param[in] aValidLen If set gives how many values to take from aPreviousResult
	///
	void initFromResult(const std::vector<double>& aPreviousResult, unsigned int aValidLen=0);

protected:
	/// Valid values for the mInitType variable depicting from where the variables have been initialized.
	enum InitVarStatus
	{
		INIT_TYPE_NONE,			///< All variables to optimize should be initialized
		INIT_TYPE_TIMES,		///< All variables to optimize should be initialized except times
		INIT_TYPE_RES_4,		///< All variables to optimize should be initialized except times, w0, k, v1, v2
		INIT_TYPE_RES_5			///< All variables to optimize have been initialized
	};

private:
	/// Disabled assignment operator to avoid warning on Windows
	///
	/// @fn BranchSiteModel& operator=(const BranchSiteModel& aObj)
	///
	/// @param[in] aObj The object to be assigned
	///
	/// @return The object receiving the assignment
	///
	BranchSiteModel& operator=(const BranchSiteModel& /*aObj*/) {return *this;}


protected:
	Forest&						mForest;			///< The forest to be used
	unsigned int				mNumTimes;			///< Number of branch lengths values
	unsigned int				mNumVariables;		///< The number of extra variables (4 for H0 and 5 for H1)
	std::vector<double>			mVar;				///< Variable to optimize (first the branch lengths then the remaining variables)
	std::vector<double>			mLowerBound;		///< Lower limits for the variables to be optimized
	std::vector<double>			mUpperBound;		///< Upper limits for the variables to be optimized
	double						mProportions[4];	///< The four proportions
	double						mMaxLnL;			///< Maximum value of LnL found during optimization
	unsigned int				mNumEvaluations;	///< Counter of the likelihood function evaluations
	CacheAlignedDoubleVector	mLikelihoods;		///< Computed likelihoods at the root of all trees. Defined here to make it aligned.
	bool						mOnlyInitialStep;	///< Only the initial step is executed, no optimization
	bool						mTrace;				///< Enable maximization tracing
	unsigned int				mOptAlgo;			///< Optimization algorithm to use
	InitVarStatus				mInitType;			///< From where the variables have been initialized
	double						mDeltaForGradient;	///< Value used to change the variables to compute gradient

private:
	unsigned int				mSeed;				///< Random number generator seed to be used also by the optimizer
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
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aOptAlgo The optimization algorithm to use
	/// @param[in] aDeltaValueForGradient The variable increment to compute gradient
	///
	BranchSiteModelNullHyp(Forest& aForest, unsigned int aSeed, bool aOnlyInitialStep, bool aTrace, unsigned int aOptAlgo=0, double aDeltaValueForGradient=0.0)
		: BranchSiteModel(aForest, aForest.getNumBranches(), aForest.getNumSites(), aSeed, 4, aOnlyInitialStep, aTrace, aOptAlgo, aDeltaValueForGradient), mSet(aForest.getNumBranches(), 3), mPrevK(DBL_MAX), mPrevOmega0(DBL_MAX) {}

	/// Compute the null hypothesis log likelihood.
	///
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	///
	/// @return The log likelihood under the null hypothesis
	///
	double operator()(size_t aFgBranch);

	/// Compute the likelihood for the given forest and the given set of parameters.
	///
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	double computeLikelihood(const std::vector<double>& aVar, bool aTrace);

private:
	/// Disabled assignment operator to avoid warning on Windows
	///
	/// @fn BranchSiteModelNullHyp& operator=(const BranchSiteModelNullHyp& aObj)
	///
	/// @param[in] aObj The object to be assigned
	///
	/// @return The object receiving the assignment
	///
	BranchSiteModelNullHyp& operator=(const BranchSiteModelNullHyp& /*aObj*/) {return *this;}

private:
	TransitionMatrix 	mQw0;			///< Q matrix for the omega0 case
	TransitionMatrix 	mQ1;			///< Q matrix for the omega1 == 1 case
	ProbabilityMatrixSet mSet;			///< Set of matrices used for the tree visits
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
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aOptAlgo The optimization algorithm to use
	/// @param[in] aDeltaValueForGradient The variable increment to compute gradient
	///
	BranchSiteModelAltHyp(Forest& aForest, unsigned int aSeed, bool aOnlyInitialStep, bool aTrace, unsigned int aOptAlgo=0, double aDeltaValueForGradient=0.0)
		: BranchSiteModel(aForest, aForest.getNumBranches(), aForest.getNumSites(), aSeed, 5, aOnlyInitialStep, aTrace, aOptAlgo, aDeltaValueForGradient), mSet(aForest.getNumBranches(), 4), mPrevK(DBL_MAX), mPrevOmega0(DBL_MAX), mPrevOmega2(DBL_MAX) {}

	/// Compute the alternative hypothesis log likelihood.
	///
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	///
	/// @return The log likelihood under the alternative hypothesis
	///
	double operator()(size_t aFgBranch);

	/// Compute the likelihood for the given forest and the given set of parameters.
	///
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace If set visualize the best result so far
	///
	/// @return The maximum Likelihood value
	///
	double computeLikelihood(const std::vector<double>& aVar, bool aTrace);

private:
	/// Disabled assignment operator to avoid warning on Windows
	///
	/// @fn BranchSiteModelAltHyp& operator=(const BranchSiteModelAltHyp& aObj)
	///
	/// @param[in] aObj The object to be assigned
	///
	/// @return The object receiving the assignment
	///
	BranchSiteModelAltHyp& operator=(const BranchSiteModelAltHyp& /*aObj*/) {return *this;}


private:
	TransitionMatrix    mQw0;			///< Q matrix for the omega0 case
	TransitionMatrix    mQw2;			///< Q matrix for the omega2 case
	TransitionMatrix    mQ1;  			///< Q matrix for the omega1 == 1 case
	ProbabilityMatrixSet mSet;			///< Set of matrices used for the tree visits
	double				mPrevK;			///< Previous k value used to compute matrices
	double				mPrevOmega0;	///< Previous w0 value used to compute matrices
	double				mPrevOmega2;	///< Previous w2 value used to compute matrices
	double				mScaleQw0;		///< Scale value for Qw0
	double				mScaleQw2;		///< Scale value for Qw2
	double				mScaleQ1;		///< Scale value for Q1
};


#endif
