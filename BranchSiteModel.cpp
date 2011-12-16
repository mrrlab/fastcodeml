
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


void BranchSiteModel::printVar(const std::vector<double>& aVars, double aLnl) const
{
	// Write the data with an uniform precision
	std::streamsize prec = std::cerr.precision(7);
	std::cerr.setf(std::ios::fixed);

	// Write the LnL value (if set)
	if(aLnl != DBL_MAX) std::cerr << std::endl << aLnl << std::endl;

	// Print all the variables
	std::vector<double>::const_iterator ix;
	int k;
	double v0 = 0;
	for(ix=aVars.begin(),k = -(int)mNumTimes; ix != aVars.end(); ++ix,++k)
	{
		switch(k)
		{
		case 0:
			std::cerr << std::endl;
			std::cerr <<   "w0: " << *ix;
			break;
		case 1:
			std::cerr <<  "  k: " << *ix;
			break;
		case 2:
			std::cerr << "  v0: " << *ix;
			v0 = *ix;
			break;
		case 3:
			std::cerr << "  v1: " << *ix;
			{
				double p[4];
				getProportions(v0, *ix, p);
				std::cerr << "  [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]";
			}
			break;
		case 4:
			std::cerr << "  w2: " << *ix;
			break;
		default:
			std::cerr << *ix << ' ';
			break;
		}
	}
	std::cerr << std::endl;
	std::cerr.precision(prec);
}


double BranchSiteModelNullHyp::operator()(size_t aFgBranch)
{
	unsigned int i;

	// Initialize the variables to be optimized
	if(mTimesFromTree)
	{
		// Initialize branch lengths from the phylo tree
		mForest.setTimesFromLengths(mVar);

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

	// Check the initial values are inside the domain
	for(i=0; i < mNumTimes+4; ++i)
	{
		if(mVar[i] < mLowerBound[i]) mVar[i] = mLowerBound[i];
		if(mVar[i] > mUpperBound[i]) mVar[i] = mUpperBound[i];
	}

	// Initialize the variables used to avoid unneeded recomputing
	mPrevK      = DBL_MAX;
	mPrevOmega0 = DBL_MAX;

	// Run the optimizer
	return maximizeLikelihood(aFgBranch);
}


double BranchSiteModelAltHyp::operator()(size_t aFgBranch, const double* aInitFromH0)
{
	unsigned int i;

	// Initialize the variables to be optimized from the H0 values
	if(aInitFromH0)
	{
		mVar.assign(aInitFromH0, aInitFromH0+mNumTimes+4);
		mVar.push_back(1.001);
	}
	else if(mTimesFromTree)
	{
		// Initialize branch lengths from the phylo tree
		mForest.setTimesFromLengths(mVar);

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
		mVar[mNumTimes+4] = 1.14833;											// w2
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

	// Check the initial values are inside the domain
	for(i=0; i < mNumTimes+5; ++i)
	{
		if(mVar[i] < mLowerBound[i]) mVar[i] = mLowerBound[i];
		if(mVar[i] > mUpperBound[i]) mVar[i] = mUpperBound[i];
	}

	// Initialize the variables used to avoid unneeded recomputing
	mPrevK      = DBL_MAX;
	mPrevOmega0 = DBL_MAX;
	mPrevOmega2 = DBL_MAX;

	// Run the optimizer
	return maximizeLikelihood(aFgBranch);
}


double BranchSiteModelNullHyp::computeLikelihood(unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Check if steps can be skipped
	const bool changed_w0 = BranchSiteModel::isDifferent(aVar[mNumTimes+0], mPrevOmega0);
	const bool changed_k  = BranchSiteModel::isDifferent(aVar[mNumTimes+1], mPrevK);
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
	const double fg_scale = mProportions[0]*mScaleQw0 +
							mProportions[1]*mScaleQ1  +
							mProportions[2]*mScaleQ1  +
							mProportions[3]*mScaleQ1;

	const double bg_scale = 1./(mProportions[0]+ mProportions[1])*(mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1);

	// Fill the Transition Matrix sets
	mSet.computeMatrixSetH0(mQw0, mQ1, bg_scale, fg_scale, mForest.adjustFgBranchIdx(aFgBranch), aVar);

	// Compute likelihoods
	mForest.computeLikelihoods(mSet, mLikelihoods);

	// For all (valid) sites. Don't parallelize: time increase and the results are errant
	const size_t num_sites = mForest.getNumSites();
	const double* mult = mForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		//double p = mProportions[0]*mLikelihoods[0*num_sites+site] +
		//		   (mProportions[1]+mProportions[3])*mLikelihoods[1*num_sites+site] +
		//		   mProportions[2]*mLikelihoods[2*num_sites+site];
		double p = mLikelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= mProportions[0];
		double x = mLikelihoods[1*num_sites+site];
		if(x > 0) p += (mProportions[1]+mProportions[3])*x;
		x = mLikelihoods[2*num_sites+site];
		if(x > 0) p += mProportions[2]*x;

		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];

		// DBG
		//if(p <= 0)
		//{
		//	std::cerr << std::setw(4) << site << ' ';
		//	std::cerr << std::setw(14) << mLikelihoods[0*num_sites+site] << ' ';
		//	std::cerr << std::setw(14) << mLikelihoods[1*num_sites+site] << ' ';
		//	std::cerr << std::setw(14) << mLikelihoods[2*num_sites+site] << std::endl;
		//}
	}

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
	}
	//std::cerr << lnl << std::endl;

	return lnl;
}

	
double BranchSiteModelAltHyp::computeLikelihood(unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Check if steps can be skipped
	const bool changed_w0 = BranchSiteModel::isDifferent(aVar[mNumTimes+0], mPrevOmega0);
	const bool changed_w2 = BranchSiteModel::isDifferent(aVar[mNumTimes+4], mPrevOmega2);
	const bool changed_k  = BranchSiteModel::isDifferent(aVar[mNumTimes+1], mPrevK);
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
	const double fg_scale = mProportions[0]*mScaleQw0 +
							mProportions[1]*mScaleQ1  +
							mProportions[2]*mScaleQw2 +
							mProportions[3]*mScaleQw2;

	const double bg_scale = 1./(mProportions[0]+ mProportions[1])*(mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1);

	// Fill the Transition Matrix sets
	mSet.computeMatrixSetH1(mQw0, mQ1, mQw2, bg_scale, fg_scale, mForest.adjustFgBranchIdx(aFgBranch), aVar);

	// Compute likelihoods
	mForest.computeLikelihoods(mSet, mLikelihoods);

	// For all (valid) sites. Don't parallelize: time increase and the results are errant
	const size_t num_sites = mForest.getNumSites();
	const double* mult = mForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		//double p = mProportions[0]*mLikelihoods[0*num_sites+site] +
		//		     mProportions[1]*mLikelihoods[1*num_sites+site] +
		//		     mProportions[2]*mLikelihoods[2*num_sites+site] +
		//		     mProportions[3]*mLikelihoods[3*num_sites+site];
		//
		double p = mLikelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= mProportions[0];
		double x = mLikelihoods[1*num_sites+site];
		if(x > 0) p += mProportions[1]*x;
		x = mLikelihoods[2*num_sites+site];
		if(x > 0) p += mProportions[2]*x;
		x = mLikelihoods[3*num_sites+site];
		if(x > 0) p += mProportions[3]*x;

		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];
	}

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
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
	/// @param[in] aFgBranch The identifier for the branch marked as foreground branch
	/// @param[in] aTrace If set the optimization progress is traced
	/// @param[in] aUpper Upper limit for the variables (to constrain the gradient computation)
	///
	MaximizerFunction(BranchSiteModel* aModel, unsigned int aFgBranch, bool aTrace, std::vector<double>& aUpper)
					: mModel(aModel), mFgBranch(aFgBranch), mTrace(aTrace), mUpper(aUpper) {}

	/// Functor.
	/// It computes the function and the gradient if needed.
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	/// @return The evaluated function
	///
	double operator()(const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		double f0 = mModel->computeLikelihood(mFgBranch, aVars, mTrace);

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
			const double v = aVars[i];
			double eh = SMALL_DIFFERENCE * (fabs(v)+1.);
			
			x[i] += eh;
			if(x[i] >= mUpper[i]) {x[i] -= 2*eh; eh = -eh;}

			const double f1 = mModel->computeLikelihood(mFgBranch, x, false);

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
	unsigned int		mFgBranch;	///< Branch number of the foreground branch
	bool				mTrace;		///< If set traces the optimization progresses
	std::vector<double>	mUpper;		///< Upper limit of the variables to constrain the interval on which the gradient should be computed
};


double BranchSiteModel::maximizeLikelihood(size_t aFgBranch)
{
	// Print starting values
	if(mTrace)
	{
		std::cerr << std::endl;
		std::cerr << "*****************************************" << std::endl;
		std::cerr << "*** Starting branch " << aFgBranch << std::endl;
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
	if(mOnlyInitialStep) return computeLikelihood(aFgBranch, mVar, mTrace);

	// Select the maximizer algorithm
	std::auto_ptr<nlopt::opt> opt;
	switch(mOptAlgo)
	{
	case OPTIM_LD_LBFGS:
		opt.reset(new nlopt::opt(nlopt::LD_LBFGS,   mNumTimes+mNumVariables));
		break;

	case OPTIM_LN_BOBYQA:
		opt.reset(new nlopt::opt(nlopt::LN_BOBYQA,  mNumTimes+mNumVariables));
		break;

	case OPTIM_LN_COBYLA:
		opt.reset(new nlopt::opt(nlopt::LN_COBYLA,  mNumTimes+mNumVariables));
		break;

	case OPTIM_MLSL_LDS:
		opt.reset(new nlopt::opt(nlopt::G_MLSL_LDS, mNumTimes+mNumVariables));
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
		MaximizerFunction compute(this, aFgBranch, mTrace, mUpperBound);

		opt->set_max_objective(MaximizerFunction::wrap, &compute);

		nlopt::result result = opt->optimize(mVar, maxl);

		// Print the final optimum value
		if(mTrace)
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

