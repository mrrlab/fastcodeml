
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
	/// @param[in] aNumBranches Number of tree branches
	/// @param[in] aNumVariables Number of extra variables (k, w0, w2, p0, p1)
	/// @param[in] aSeed Random number generator seed
	///
	BranchSiteModel(unsigned int aNumBranches, unsigned int aNumVariables, unsigned int aSeed)
	{
		mNumTimes     = aNumBranches;
		mNumVariables = aNumVariables;
		mVar.reserve(mNumTimes+mNumVariables);
		mVar.resize(mNumTimes+mNumVariables);
		mLowerBound.reserve(mNumTimes+mNumVariables);
		mLowerBound.resize(mNumTimes+mNumVariables);
		mUpperBound.reserve(mNumTimes+mNumVariables);
		mUpperBound.resize(mNumTimes+mNumVariables);
		mSeed = aSeed;
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

	/// Compute the four site proportions from the two values in the optimization variables
	///
	/// @param[in] aV0 The first optimization variables
	/// @param[in] aV1 The second optimization variables
	/// @param[out] aProportions The four proportions output (size: 4)
	///
	void getProportions(double aV0, double aV1, double* aProportions) const;

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
	///
	/// @return The maximum Likelihood value
	///
	double maximizeLikelihood(Forest& aForest, unsigned int aFgBranch, bool aOnlyInitialStep, bool aTrace);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace Set to trace the cycle result
	///
	/// @return The maximum Likelihood value
	///
	virtual double oneCycleMaximizer(Forest& aForest, unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace) =0;

protected:
	unsigned int		mNumTimes;			///< Number of branch lengths
	unsigned int		mNumVariables;		///< The number of extra variables (4 for H0 and 5 for H1)
	std::vector<double>	mVar;				///< Variable to optimize (first the times then the remaining variables)
	std::vector<double>	mLowerBound;		///< Lower limits for the variables to be optimized
	std::vector<double>	mUpperBound;		///< Upper limits for the variables to be optimized
	double				mProportions[4];	///< The four proportions
	double				mCodonFreq[N];		///< %Codon frequencies
	std::vector<double>	mLnLsite;			///< Loglik values for each site
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
	BranchSiteModelNullHyp(unsigned int aNumBranches, unsigned int aSeed)
		: BranchSiteModel(aNumBranches, 4, aSeed), mSet(aNumBranches, 3) {}

	/// Compute the null hypothesis log likelihood.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aOnlyInitialStep If true no optimization is done, only the initial step is run
	/// @param[in] aTimesFromTree Initial times are from the file plus fixed values for the other variables
	/// @param[in] aTrace If set the maximization is traced
	///
	/// @return The log likelihood for the null hypothesis
	///
	double computeModel(Forest& aForest, unsigned int aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace Set to trace the cycle result
	///
	/// @return The maximum Likelihood value
	///
	double oneCycleMaximizer(Forest& aForest, unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace);


private:
	TransitionMatrix    mQw0;	///< Q matrix for the omega0 case
	TransitionMatrix    mQ1;	///< Q matrix for the omega1 == 1 case
	TransitionMatrixSet mSet;	///< Set of matrices used for the tree visits
};


class BranchSiteModelAltHyp : public BranchSiteModel
{
public:
	/// Constructor.
	///
	/// @param[in] aNumBranches Number of tree branches
	/// @param[in] aSeed Random number generator seed
	///
	BranchSiteModelAltHyp(unsigned int aNumBranches, unsigned int aSeed)
		: BranchSiteModel(aNumBranches, 5, aSeed), mSet(aNumBranches, 4) {}

	/// Compute the alternative hypothesis log likelihood.
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aOnlyInitialStep If true no optimization is done, only the initial step is run
	/// @param[in] aTimesFromTree Initial times are from the file plus fixed values for the other variables
	/// @param[in] aTrace If set the maximization is traced
	///
	/// @return The log likelihood for the alternative hypothesis
	///
	double computeModel(Forest& aForest, unsigned int aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace);

	/// Compute one iteration of the maximum likelihood computation for the given forest
	///
	/// @param[in] aForest The forest for which the maximum likelihood should be computed
	/// @param[in] aFgBranch The number of the internal branch to be marked as foreground
	/// @param[in] aVar The optimizer variables
	/// @param[in] aTrace Set to trace the cycle result
	///
	/// @return The maximum Likelihood value
	///
	double oneCycleMaximizer(Forest& aForest, unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace);


private:
	TransitionMatrix    mQw0;	///< Q matrix for the omega0 case
	TransitionMatrix    mQw2;	///< Q matrix for the omega2 case
	TransitionMatrix    mQ1;  	///< Q matrix for the omega1 == 1 case
	TransitionMatrixSet mSet;	///< Set of matrices used for the tree visits
};

#endif
