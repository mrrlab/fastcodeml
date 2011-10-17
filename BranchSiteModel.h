
#ifndef BRANCHSITEMODEL_H
#define BRANCHSITEMODEL_H

#include <cstring>
#include <vector>

#include "TransitionMatrix.h"
#include "TransitionMatrixSet.h"
#include "Forest.h"

// Uncomment to use the original CodeML proportion definition
//#define USE_ORIGINAL_PROPORTIONS


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
	/// @param[in] aNumBranches Number of tree branches
	/// @param[in] aNumVariables Number of extra variables (k, w0, w2, p0, p1)
	/// @param[in] aSeed Random number generator seed
	///
	BranchSiteModel(size_t aNumBranches, unsigned int aNumVariables, unsigned int aSeed)
	{
		mNumTimes     = aNumBranches;
		mNumVariables = aNumVariables;
		mVar.reserve(mNumTimes+mNumVariables);
		mVar.resize(mNumTimes+mNumVariables);
		mLowerBound.reserve(mNumTimes+mNumVariables);
		mUpperBound.reserve(mNumTimes+mNumVariables);
		mSeed = aSeed;
		mNumEvaluations = 0;
		mMaxLnL = -DBL_MAX;
	}

	/// Destructor.
	///
	virtual ~BranchSiteModel() {}

	/// Set the times on the tree from the variables
	///
	/// @param[in,out] aForest The forest to be updated
	///
	void saveComputedTimes(Forest& aForest) const
	{
		aForest.setLengthsFromTimes(mVar);
	}

	/// Formatted print of the maximizer variables array
	///
	/// @param[in] aVars The variables array to be printed
	///
	void printVar(const std::vector<double>& aVars) const;

	/// Compute the maximum likelihood for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aOnlyInitialStep If set do not maximize, compute only the starting point
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aOptAlgo The optimization method to be used
	///
	/// @return The maximum Likelihood value
	///
	double maximizeLikelihood(Forest& aForest, size_t aFgBranch, bool aOnlyInitialStep, bool aTrace, unsigned int aOptAlgo=0);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace Set to trace the cycle result
	///
	/// @return The maximum Likelihood value
	///
	virtual double oneCycleMaximizer(Forest& aForest, size_t aFgBranch, const std::vector<double>& aVar, bool aTrace) =0;

	const double* getStartingValues(void) const {return &mVar[0];}

	enum {
		OPTIM_LD_LBFGS		= 0,
		OPTIM_LN_BOBYQA		= 1,
		OPTIM_MLSL_LDS		= 99
	};

protected:
	/// Compute the four site proportions from the two values in the optimization variables
	///
	/// @param[in] aV0 The first optimization variables
	/// @param[in] aV1 The second optimization variables
	/// @param[out] aProportions The four proportions output
	///
	inline void getProportions(double aV0, double aV1, double* aProportions) const
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
	inline static bool isDifferent(double aFirst, double aSecond)
	{
		static const double TOL = 1e-7;
		double diff = aFirst - aSecond;
		return (diff > TOL || diff < -TOL);
	}

protected:
	size_t				mNumTimes;			///< Number of branch lengths
	unsigned int		mNumVariables;		///< The number of extra variables (4 for H0 and 5 for H1)
	std::vector<double>	mVar;				///< Variable to optimize (first the times then the remaining variables)
	std::vector<double>	mLowerBound;		///< Lower limits for the variables to be optimized
	std::vector<double>	mUpperBound;		///< Upper limits for the variables to be optimized
	double				mProportions[4];	///< The four proportions
	double				mCodonFreq[N];		///< %Codon frequencies
	double				mMaxLnL;			///< Maximum value of LnL found during optimization
	unsigned int		mNumEvaluations;	///< Counter of the likelihood function evaluations

private:
	unsigned int		mSeed;				///< Random number generator seed to be passed to the optimizer
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
	/// @param[in] aNumBranches Number of tree branches
	/// @param[in] aSeed Random number generator seed
	///
	BranchSiteModelNullHyp(size_t aNumBranches, unsigned int aSeed)
		: BranchSiteModel(aNumBranches, 4, aSeed), mSet(aNumBranches, 3), mPrevK(DBL_MAX), mPrevOmega0(DBL_MAX) {}

	/// Compute the null hypothesis log likelihood.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aOnlyInitialStep If true no optimization is done, only the initial step is run
	/// @param[in] aTimesFromTree Initial times are from the file plus fixed values for the other variables
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aOptAlgo The optimization algorithm to use
	///
	/// @return The log likelihood for the null hypothesis
	///
	double computeModel(Forest& aForest, size_t aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace, unsigned int aOptAlgo=0);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace Set to trace the cycle result
	///
	/// @return The maximum Likelihood value
	///
	double oneCycleMaximizer(Forest& aForest, size_t aFgBranch, const std::vector<double>& aVar, bool aTrace);


private:
	TransitionMatrix    mQw0;			///< Q matrix for the omega0 case
	TransitionMatrix    mQ1;			///< Q matrix for the omega1 == 1 case
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
	/// @param[in] aNumBranches Number of tree branches
	/// @param[in] aSeed Random number generator seed
	///
	BranchSiteModelAltHyp(size_t aNumBranches, unsigned int aSeed)
		: BranchSiteModel(aNumBranches, 5, aSeed), mSet(aNumBranches, 4), mPrevK(DBL_MAX), mPrevOmega0(DBL_MAX), mPrevOmega2(DBL_MAX) {}

	/// Compute the alternative hypothesis log likelihood.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aOnlyInitialStep If true no optimization is done, only the initial step is run
	/// @param[in] aTimesFromTree Initial times are from the file plus fixed values for the other variables
	/// @param[in] aTrace If set the maximization is traced
	/// @param[in] aInitFromH0 If not null uses these results from H0 to initalize H1 values
	/// @param[in] aOptAlgo The optimization algorithm to use
	///
	/// @return The log likelihood for the alternative hypothesis
	///
	double computeModel(Forest& aForest, size_t aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace, const double* aInitFromH0, unsigned int aOptAlgo=0);


	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace Set to trace the cycle result
	///
	/// @return The maximum Likelihood value
	///
	double oneCycleMaximizer(Forest& aForest, size_t aFgBranch, const std::vector<double>& aVar, bool aTrace);


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
