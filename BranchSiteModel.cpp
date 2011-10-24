
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include <memory>

#include "nlopt.hpp"

#include "BranchSiteModel.h"
#include "MathSupport.h"
#include "Exceptions.h"

/// How much a variable should be changed to compute gradient
static const double SMALL_DIFFERENCE = sqrt(DBL_EPSILON);

// Starting value for the maximum likelihood
// Beware: -HUGE_VAL is too low and on Linux it is converted to -Inf (with subsequent NLopt crash)
static const double VERY_LOW_LIKELIHOOD = -1e14;

void BranchSiteModel::printVar(const std::vector<double>& aVars) const
{
	std::vector<double>::const_iterator ix;
	int k;
	double v0 = 0;
	for(ix=aVars.begin(),k = -(int)mNumTimes; ix != aVars.end(); ++ix,++k)
	{
		switch(k)
		{
		case 0:
			std::cerr << std::endl;
			std::cerr << std::setprecision(4) <<   "w0: " << *ix;
			break;
		case 1:
			std::cerr << std::setprecision(4) << "  k: " << *ix;
			break;
		case 2:
			std::cerr << std::setprecision(4) << "  v0: " << *ix;
			v0 = *ix;
			break;
		case 3:
			std::cerr << std::setprecision(4) << "  v1: " << *ix;
			{
				double p[4];
				getProportions(v0, *ix, p);
				std::cerr << "  [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]";
			}
			break;
		case 4:
			std::cerr << std::setprecision(4) << "  w2: " << *ix;
			break;
		default:
			std::cerr << std::setprecision(7) << *ix << ' ';
			break;
		}
	}
	std::cerr << std::endl;
}


double BranchSiteModelNullHyp::computeModel(Forest& aForest, size_t aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace, unsigned int aOptAlgo)
{
	unsigned int i;

	// Initialize the variables to be optimized
	if(aTimesFromTree)
	{
		// Initialize branch lengths from the phylo tree
		aForest.setTimesFromLengths(mVar);

		// Initialization as in CodeML (seems)
		mVar[mNumTimes+0] = 0.235087;											// w0
		mVar[mNumTimes+1] = 0.4;												// k

#ifdef USE_ORIGINAL_PROPORTIONS
		mVar[mNumTimes+2] = 1.04885;											// x0 -> p0
		mVar[mNumTimes+3] = 0.12437;											// x1 -> p1
#else
		mVar[mNumTimes+2] = 0.813836;											// p0+p1
		mVar[mNumTimes+3] = 0.7260434;											// p0/(p0+p1)
#endif
	}
	else
	{
		for(i=0; i < mNumTimes; ++i) mVar[i] = rand()/(double)RAND_MAX*.1+0.01;	// T
		mVar[mNumTimes+0] = rand()/(double)RAND_MAX*0.8 + 0.1;					// w0
		mVar[mNumTimes+1] = 2.0;												// k
#ifdef USE_ORIGINAL_PROPORTIONS
		mVar[mNumTimes+2] = 1.0 + 0.2 * rand()/(double)RAND_MAX;				// x0 -> p0
		mVar[mNumTimes+3] = 0.2*rand()/(double)RAND_MAX;						// x1 -> p1
#else
		mVar[mNumTimes+2] = rand()/(double)RAND_MAX;							// p0+p1
		mVar[mNumTimes+3] =	rand()/(double)RAND_MAX;							// p0/(p0+p1)
#endif
	}

	// Set lower constrains
	mLowerBound.assign(mNumTimes, 4e-6);	// T
	mLowerBound.push_back(1e-6);			// w0
	mLowerBound.push_back(0.0001);			// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mLowerBound.push_back(-99.0);			// x0 -> p0
	mLowerBound.push_back(-99.0);			// x1 -> p1
#else
	mLowerBound.push_back(0.0);				// p0+p1
	mLowerBound.push_back(0.0);				// p0/(p0+p1)
#endif

#if 0
	for(i=0; i < mNumTimes; ++i) mLowerBound[i] = 4e-6;		// T
	mLowerBound[mNumTimes+0] = 1e-6;						// w0
	mLowerBound[mNumTimes+1] = 0.0001;						// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mLowerBound[mNumTimes+2] = -99;							// x0 -> p0
	mLowerBound[mNumTimes+3] = -99;							// x1 -> p1
#else
	mLowerBound[mNumTimes+2] = 0.0;							// p0+p1
	mLowerBound[mNumTimes+3] = 0.0;							// p0/(p0+p1)
#endif
#endif

	// Set upper constrains
	mUpperBound.assign(mNumTimes, 50.0);	// T
	mUpperBound.push_back(1.0);				// w0
	mUpperBound.push_back(20.0);			// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mUpperBound.push_back(99.0);			// x0 -> p0
	mUpperBound.push_back(99.0);			// x1 -> p1
#else
	mUpperBound.push_back(1.0);				// p0+p1
	mUpperBound.push_back(1.0);				// p0/(p0+p1)
#endif

#if 0
	for(i=0; i < mNumTimes; ++i) mUpperBound[i] = 50.0;		// T
	mUpperBound[mNumTimes+0] = 1.0;							// w0
	mUpperBound[mNumTimes+1] = 20.0;						// k (in the old code is 999)
#ifdef USE_ORIGINAL_PROPORTIONS
	mUpperBound[mNumTimes+2] = 99;							// x0 -> p0
	mUpperBound[mNumTimes+3] = 99;							// x1 -> p1
#else
	mUpperBound[mNumTimes+2] = 1.0;							// p0+p1
	mUpperBound[mNumTimes+3] = 1.0;							// p0/(p0+p1)
#endif
#endif

	// Check the initial values are inside the domain
	for(i=0; i < mNumTimes+4; ++i)
	{
		if(mVar[i] < mLowerBound[i]) mVar[i] = mLowerBound[i];
		if(mVar[i] > mUpperBound[i]) mVar[i] = mUpperBound[i];
	}

	// Run the optimizer
	return maximizeLikelihood(aForest, aFgBranch, aOnlyInitialStep, aTrace, aOptAlgo);
}


double BranchSiteModelAltHyp::computeModel(Forest& aForest, size_t aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace, const double* aInitFromH0, unsigned int aOptAlgo)
{
	unsigned int i;

	// Initialize the variables to be optimized from the H0 values
	if(aInitFromH0)
	{
		mVar.assign(aInitFromH0, aInitFromH0+mNumTimes+4);
		mVar.push_back(1.001);
	}
	else 
	{
		// Initialize from the tree (plus fixed values guessed from CodeML code)
		if(aTimesFromTree)
		{
			// Initialize branch lengths from the phylo tree
			aForest.setTimesFromLengths(mVar);

			// Initialization as in CodeML (seems)
			mVar[mNumTimes+0] = 0.235087;											// w0
			mVar[mNumTimes+1] = 0.4;												// k
			mVar[mNumTimes+4] = 1.14833;											// w2

#ifdef USE_ORIGINAL_PROPORTIONS
			mVar[mNumTimes+2] = 1.04885;											// x0 -> p0
			mVar[mNumTimes+3] = 0.12437;											// x1 -> p1
#else
			mVar[mNumTimes+2] = 0.813836;											// p0+p1
			mVar[mNumTimes+3] = 0.7260434;											// p0/(p0+p1)
#endif
		}
		else
		{
			// Initialize the variables to be optimized from random values
			for(i=0; i < mNumTimes; ++i) mVar[i] = rand()/(double)RAND_MAX*.1+0.01;		// T
			mVar[mNumTimes+0] = rand()/(double)RAND_MAX*0.8 + 0.1;						// w0
			mVar[mNumTimes+1] = 2.0;													// k
#ifdef USE_ORIGINAL_PROPORTIONS
			mVar[mNumTimes+2] = 1.0 + 0.2 * rand()/(double)RAND_MAX;					// x0 -> p0
			mVar[mNumTimes+3] = 0.2*rand()/(double)RAND_MAX;							// x1 -> p1
#else
			mVar[mNumTimes+2] = rand()/(double)RAND_MAX;								// p0+p1
			mVar[mNumTimes+3] = rand()/(double)RAND_MAX;								// p0/(p0+p1)
#endif
			mVar[mNumTimes+4] = 1.001;													// w2
		}
	}

	// Set lower constrains
	mLowerBound.assign(mNumTimes, 4e-6);	// T
	mLowerBound.push_back(1e-6);			// w0
	mLowerBound.push_back(0.0001);			// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mLowerBound.push_back(-99.0);			// x0 -> p0
	mLowerBound.push_back(-99.0);			// x1 -> p1
#else
	mLowerBound.push_back(0.0);				// p0+p1
	mLowerBound.push_back(0.0);				// p0/(p0+p1)
#endif
	mLowerBound.push_back(1.0);				// w2

#if 0
	// Set lower constrains
	for(i=0; i < mNumTimes; ++i) mLowerBound[i] = 4e-6;		// T

	mLowerBound[mNumTimes+0] = 1e-6;						// w0
	mLowerBound[mNumTimes+1] = 0.0001;						// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mLowerBound[mNumTimes+2] = -99;							// x0 -> p0
	mLowerBound[mNumTimes+3] = -99;							// x1 -> p1
#else
	mLowerBound[mNumTimes+2] = 0.0;							// p0+p1
	mLowerBound[mNumTimes+3] = 0.0;							// p0/(p0+p1)
#endif
	mLowerBound[mNumTimes+4] = 1.0;							// w2
#endif

	// Set upper constrains
	mUpperBound.assign(mNumTimes, 50.0);	// T
	mUpperBound.push_back(1.0);				// w0
	mUpperBound.push_back(20.0);			// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mUpperBound.push_back(99.0);			// x0 -> p0
	mUpperBound.push_back(99.0);			// x1 -> p1
#else
	mUpperBound.push_back(1.0);				// p0+p1
	mUpperBound.push_back(1.0);				// p0/(p0+p1)
#endif
	mUpperBound.push_back(999.0);			// w2 (in the old code is 999)

#if 0
	// Set upper constrains
	for(i=0; i < mNumTimes; ++i) mUpperBound[i] = 50.0;		// T
	mUpperBound[mNumTimes+0] = 1.0;							// w0
	mUpperBound[mNumTimes+1] = 20.0;						// k (in the old code is 999)
#ifdef USE_ORIGINAL_PROPORTIONS
	mUpperBound[mNumTimes+2] = 99;							// x0 -> p0
	mUpperBound[mNumTimes+3] = 99;							// x1 -> p1
#else
	mUpperBound[mNumTimes+2] = 1.0;							// p0+p1
	mUpperBound[mNumTimes+3] = 1.0;							// p0/(p0+p1)
#endif
	mUpperBound[mNumTimes+4] = 999.0;						// w2 (in the old code is 999)
#endif

	// Check the initial values are inside the domain
	for(i=0; i < mNumTimes+5; ++i)
	{
		if(mVar[i] < mLowerBound[i]) mVar[i] = mLowerBound[i];
		if(mVar[i] > mUpperBound[i]) mVar[i] = mUpperBound[i];
	}

	// Run the optimizer
	return maximizeLikelihood(aForest, aFgBranch, aOnlyInitialStep, aTrace, aOptAlgo);
}


double BranchSiteModelNullHyp::oneCycleMaximizer(Forest& aForest, size_t aFgBranch, const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Check if steps can be skipped
	bool changed_w0 = BranchSiteModel::isDifferent(aVar[mNumTimes+0], mPrevOmega0);
	bool changed_k  = BranchSiteModel::isDifferent(aVar[mNumTimes+1], mPrevK);
	if(changed_w0) mPrevOmega0 = aVar[mNumTimes+0];
	if(changed_k)  mPrevK      = aVar[mNumTimes+1];

	// Fill the matrices and compute their eigen decomposition. Not worth the effort to parallelize this section.
	if(changed_w0 || changed_k)
	{
		mScaleQw0 = mQw0.fillQ(aVar[mNumTimes+0], aVar[mNumTimes+1]);
		mQw0.eigenQREV();
	}
	if(changed_k)
	{
		mScaleQ1  = mQ1.fillQ(                    aVar[mNumTimes+1]);
		mQ1.eigenQREV();
	}

	// Compute all proportions
	getProportions(aVar[mNumTimes+2], aVar[mNumTimes+3], mProportions);

	// Compute the scale values
	double fg_scale = mProportions[0]*mScaleQw0 +
					  mProportions[1]*mScaleQ1  +
					  mProportions[2]*mScaleQ1  +
					  mProportions[3]*mScaleQ1;

	//double bg_scale = mProportions[0]*scale_qw0 +
	//				  mProportions[1]*scale_q1  +
	//				  mProportions[2]*scale_qw0 +
	//				  mProportions[3]*scale_q1;
	double bg_scale = 1./(mProportions[0]+ mProportions[1])*(mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1);

	// Fill the Transition Matrix sets
	mSet.computeMatrixSetH0(mQw0, mQ1, bg_scale, fg_scale, aForest.adjustFgBranchIdx(aFgBranch), aVar);

	// Compute likelihoods
	std::vector<double> likelihoods;
	aForest.computeLikelihood(mSet, likelihoods);

	// For all (valid) sites. Don't parallelize: time increase and the results are errant
	size_t num_sites = aForest.getNumSites();
	const double* mult = aForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		//double p = mProportions[0]*likelihoods[0*num_sites+site] +
		//		   (mProportions[1]+mProportions[3])*likelihoods[1*num_sites+site] +
		//		   mProportions[2]*likelihoods[2*num_sites+site];
		double p = likelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= mProportions[0];
		double x = likelihoods[1*num_sites+site];
		if(x > 0) p += (mProportions[1]+mProportions[3])*x;
		x = likelihoods[2*num_sites+site];
		if(x > 0) p += mProportions[2]*x;

		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];

		//if(p <= 0) std::cerr << std::setw(4) << site << ' ' << std::setw(10) << p << ' ' << std::setw(10) << x << ' ' << mMaxLnL << std::endl;
		//if(p <= 0)
		//{
		//	std::cerr << std::setw(4) << site << ' ';
		//	std::cerr << std::setw(14) << likelihoods[0*num_sites+site] << ' ';
		//	std::cerr << std::setw(14) << likelihoods[1*num_sites+site] << ' ';
		//	std::cerr << std::setw(14) << likelihoods[2*num_sites+site] << std::endl;
		//}
	}

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		std::cerr << std::endl << lnl << std::endl;
		printVar(aVar);
	}
	//std::cerr << lnl << std::endl;

	return lnl;
}

	
double BranchSiteModelAltHyp::oneCycleMaximizer(Forest& aForest, size_t aFgBranch, const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Check if steps can be skipped
	bool changed_w0 = BranchSiteModel::isDifferent(aVar[mNumTimes+0], mPrevOmega0);
	bool changed_w2 = BranchSiteModel::isDifferent(aVar[mNumTimes+4], mPrevOmega2);
	bool changed_k  = BranchSiteModel::isDifferent(aVar[mNumTimes+1], mPrevK);
	if(changed_w0) mPrevOmega0 = aVar[mNumTimes+0];
	if(changed_w2) mPrevOmega2 = aVar[mNumTimes+4];
	if(changed_k)  mPrevK      = aVar[mNumTimes+1];

	// Fill the matrices and compute their eigen decomposition. Not worth the effort to parallelize this section.
	if(changed_w0 || changed_k)
	{
		mScaleQw0 = mQw0.fillQ(aVar[mNumTimes+0], aVar[mNumTimes+1]);
		mQw0.eigenQREV();
	}
	if(changed_w2 || changed_k)
	{
		mScaleQw2 = mQw2.fillQ(aVar[mNumTimes+4], aVar[mNumTimes+1]);
		mQw2.eigenQREV();
	}
	if(changed_k)
	{
		mScaleQ1  = mQ1.fillQ(                    aVar[mNumTimes+1]);
		mQ1.eigenQREV();
	}

	// Compute all proportions
	getProportions(aVar[mNumTimes+2], aVar[mNumTimes+3], mProportions);

	// Compute the scale values
	double fg_scale = mProportions[0]*mScaleQw0 +
					  mProportions[1]*mScaleQ1  +
					  mProportions[2]*mScaleQw2 +
					  mProportions[3]*mScaleQw2;
	double bg_scale = 1./(mProportions[0]+ mProportions[1])*(mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1);

	// Fill the Transition Matrix sets
	mSet.computeMatrixSetH1(mQw0, mQ1, mQw2, bg_scale, fg_scale, aForest.adjustFgBranchIdx(aFgBranch), aVar);

	// Compute likelihoods
	std::vector<double> likelihoods;
	aForest.computeLikelihood(mSet, likelihoods);

	// For all (valid) sites. Don't parallelize: time increase and the results are errant
	size_t num_sites = aForest.getNumSites();
	const double* mult = aForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		//double p = mProportions[0]*likelihoods[0*num_sites+site] +
		//		     mProportions[1]*likelihoods[1*num_sites+site] +
		//		     mProportions[2]*likelihoods[2*num_sites+site] +
		//		     mProportions[3]*likelihoods[3*num_sites+site];
		//
		double p = likelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= mProportions[0];
		double x = likelihoods[1*num_sites+site];
		if(x > 0) p += mProportions[1]*x;
		x = likelihoods[2*num_sites+site];
		if(x > 0) p += mProportions[2]*x;
		x = likelihoods[3*num_sites+site];
		if(x > 0) p += mProportions[3]*x;

		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];
	}

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		std::cerr << std::endl << lnl << std::endl;
		printVar(aVar);
	}

	return lnl;
}

/// Adapter class to pass the routine to the optimizer.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-01-13 (initial version)
///     @version 1.0
///
///
class MaximizerFunction
{
public:
	/// Constructor.
	/// Saves the parameters for the routine.
	///
	/// @param[in] aModel The model to be evaluated
	/// @param[in] aForest The phylogenetic tree forest
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aTrace If set the optimization progress is traced
	/// @param[in] aUpper Upper limit for the variables (to constrain the gradient computation)
	///
	MaximizerFunction(BranchSiteModel* aModel, Forest* aForest, unsigned int aFgBranch, bool aTrace, std::vector<double>& aUpper)
					: mModel(aModel), mForest(aForest), mFgBranch(aFgBranch), mTrace(aTrace), mUpper(aUpper) {}

	/// Functor.
	/// It computes the function and the gradient if needed.
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	/// @return The evaluated function
	///
	double operator()(const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		double f0 = mModel->oneCycleMaximizer(*mForest, mFgBranch, aVars, mTrace);

		if(!aGrad.empty())
		{
			gradient(f0, aVars, aGrad);
		}
		return f0;
	}

	/// Compute the function gradient
	///
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values computed
	///
	void gradient(double aPointValue, const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		std::vector<double> x = aVars;
		unsigned int vs = aVars.size();
		for(unsigned int i=0; i < vs; ++i)
		{
			double v = aVars[i];
			double eh = SMALL_DIFFERENCE * (fabs(v)+1.);
			
			x[i] += eh;
			if(x[i] >= mUpper[i]) {x[i] -= 2*eh; eh = -eh;}

			double f1 = mModel->oneCycleMaximizer(*mForest, mFgBranch, x, false);

			aGrad[i] = (f1-aPointValue)/eh;

			x[i] = v;
		}
	}

	/// Wrapper to be passed to the optimizer
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	/// @param[in] aData Opaque pointer containing the function to be passed to the optimizer
	///
	/// @return The evaluated function
	///
	static double wrap(const std::vector<double>& aVars, std::vector<double>& aGrad, void* aData)
	{
		return (*reinterpret_cast<MaximizerFunction*>(aData))(aVars, aGrad);
	}

private:
	BranchSiteModel*	mModel;		///< Pointer to the model to be evaluated
	Forest*				mForest;	///< Pointer to the forest
	unsigned int		mFgBranch;	///< Branch number of the foreground branch
	bool				mTrace;		///< If set traces the optimization progresses
	std::vector<double>	mUpper;		///< Upper limit of the variables to constrain the interval on which the gradient should be computed
};


double BranchSiteModel::maximizeLikelihood(Forest& aForest, size_t aFgBranch, bool aOnlyInitialStep, bool aTrace, unsigned int aOptAlgo)
{
	// Print starting values
	if(aTrace)
	{
		std::cerr << std::endl;
		std::cerr << "*****************************************" << std::endl;
		std::cerr << "*** Starting" << std::endl;
		printVar(mVar);
		std::cerr << "*** Upper" << std::endl;
		printVar(mUpperBound);
		std::cerr << "*** Lower" << std::endl;
		printVar(mLowerBound);
		std::cerr << std::endl;
	}

	// Initialize the maximum value found and the function evaluations counter
	mMaxLnL = VERY_LOW_LIKELIHOOD;
	mNumEvaluations = 0;

	// If only the initial step is requested, do it and return
	if(aOnlyInitialStep) return oneCycleMaximizer(aForest, aFgBranch, mVar, aTrace);

	// Select the maximizer algorithm
	std::auto_ptr<nlopt::opt> opt;
	switch(aOptAlgo)
	{
	case OPTIM_LD_LBFGS:
		opt.reset(new nlopt::opt(nlopt::LD_LBFGS,    mNumTimes+mNumVariables));
		break;

	case OPTIM_LN_BOBYQA:
		opt.reset(new nlopt::opt(nlopt::LN_BOBYQA,   mNumTimes+mNumVariables));
		break;

	case OPTIM_LN_COBYLA:
		opt.reset(new nlopt::opt(nlopt::LN_COBYLA,   mNumTimes+mNumVariables));
		break;

	case OPTIM_MLSL_LDS:
		opt.reset(new nlopt::opt(nlopt::G_MLSL_LDS,  mNumTimes+mNumVariables));
		{
		// For global optimization put a timeout of one hour
		opt->set_maxtime(60*60);

		// This algorithm requires a local optimizer, add it
		nlopt::opt local_opt(nlopt::LN_BOBYQA, mNumTimes+mNumVariables);
		opt->set_local_optimizer(local_opt);
		}
		break;

	//	opt = new nlopt::opt(nlopt::GN_DIRECT_L, mNumTimes+mNumVariables);
	//	opt = new nlopt::opt(nlopt::GN_ISRES,    mNumTimes+mNumVariables);
	//	opt = new nlopt::opt(nlopt::LN_COBYLA,   mNumTimes+mNumVariables);
	//	opt = new nlopt::opt(nlopt::LN_BOBYQA,   mNumTimes+mNumVariables);
	//	opt = new nlopt::opt(nlopt::LN_SBPLX,    mNumTimes+mNumVariables);
	//	opt = new nlopt::opt(nlopt::LD_MMA,      mNumTimes+mNumVariables);
	//	opt = new nlopt::opt(nlopt::LD_SLSQP,    mNumTimes+mNumVariables);

	default:
		throw FastCodeMLFatal("Invalid optimization algorithm");
	}

	// Initialize bounds and termination criteria
	try
	{
		opt->set_lower_bounds(mLowerBound);
		opt->set_upper_bounds(mUpperBound);
    	opt->set_ftol_abs(1e-4);
		nlopt::srand((unsigned long)mSeed);
	}
	catch(std::exception& e)
	{
		std::cerr << "Exception during inizialization: " << e.what() << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	// Optimize the function
	double maxl = 0;
	try
	{
		MaximizerFunction compute(this, &aForest, aFgBranch, aTrace, mUpperBound);

		opt->set_max_objective(MaximizerFunction::wrap, &compute);

		nlopt::result result = opt->optimize(mVar, maxl);

		// Print the final optimum value
		if(aTrace)
		{
			std::cerr << std::endl << "Function invocations:       " << mNumEvaluations << std::endl;
			switch(result)
			{
			case nlopt::SUCCESS:
				break;

			case nlopt::STOPVAL_REACHED:
				std::cerr << "Optimization stopped because stopval was reached." << std::endl;
				break;

			case nlopt::FTOL_REACHED:
				std::cerr << "Optimization stopped because ftol_rel or ftol_abs was reached." << std::endl;
				break;

			case nlopt::XTOL_REACHED:
				std::cerr << "Optimization stopped because xtol_rel or xtol_abs was reached." << std::endl;
				break;

			case nlopt::MAXEVAL_REACHED:
				std::cerr << "Optimization stopped because maxeval was reached." << std::endl;
				break;

			case nlopt::MAXTIME_REACHED:
				std::cerr << "Optimization stopped because maxtime was reached." << std::endl;	
				break;

			default:
				std::cerr << "Other reason: " << result << std::endl;	
				break;
			}
			std::cerr << "Final log-likelihood value: " << maxl << std::endl;
			printVar(mVar);
		}
	}
	catch(std::exception& e)
	{
		std::cerr << "Exception in computation: " << e.what() << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	return maxl;
}

