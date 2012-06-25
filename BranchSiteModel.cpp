
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <memory>

#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4267) // warning C4267: 'argument' : conversion from 'size_t' to 'unsigned int', possible loss of data
#endif

#include "nlopt.hpp"

#ifdef _MSC_VER
    #pragma warning(pop)
#endif

#include "BranchSiteModel.h"
#include "MathSupport.h"
#include "Exceptions.h"
#include "CodeMLoptimizer.h"
#include "ParseParameters.h"

/// Starting value for the computed maximum likelihood.
/// Beware: -HUGE_VAL is too low and on Linux it is converted to -Inf (with subsequent NLopt crash)
static const double VERY_LOW_LIKELIHOOD = -1e14;


void BranchSiteModel::printVar(const std::vector<double>& aVars, double aLnl, std::ostream& aOut) const
{
	// Write the data with an uniform precision
	std::streamsize prec = aOut.precision(7);
	aOut.setf(std::ios::fixed, std::ios::floatfield);

	// Write the LnL value (if set)
	if(aLnl < DBL_MAX) aOut << std::endl << aLnl << std::endl;

	// Print all variables formatted to be readable
	int k;
	double v0 = 0;
	std::vector<double>::const_iterator ix = aVars.begin();
	const std::vector<double>::const_iterator end = aVars.end();
	for(k = -static_cast<int>(mNumTimes); ix != end; ++ix,++k)
	{
		switch(k)
		{
		case 0:
			aOut << std::endl;
			aOut <<   "w0: " << *ix;
			break;
		case 1:
			aOut <<  "  k: " << *ix;
			break;
		case 2:
			aOut << "  v0: " << *ix;
			v0 = *ix;
			break;
		case 3:
			aOut << "  v1: " << *ix;
			{
				double p[4];
				getProportions(v0, *ix, p);
				aOut << "  [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]";
			}
			break;
		case 4:
			aOut << "  w2: " << *ix;
			break;
		default:
			aOut << *ix << ' ';
			break;
		}
	}
	aOut << std::endl;
	aOut.precision(prec);
}


void BranchSiteModel::setLimits(unsigned int aNumTimes, unsigned int aNumVariables)
{
	// Reserve space
	mLowerBound.reserve(aNumTimes+aNumVariables);	mUpperBound.reserve(aNumTimes+aNumVariables);
	
	// Set lower constrains							// Set upper constrains
	mLowerBound.assign(aNumTimes, 4e-6);			mUpperBound.assign(aNumTimes, 50.0);	// T
	mLowerBound.push_back(1e-6);					mUpperBound.push_back(1.0);				// w0
	mLowerBound.push_back(0.0001);					mUpperBound.push_back(20.0);			// k
#ifdef USE_ORIGINAL_PROPORTIONS
	mLowerBound.push_back(-99.0);					mUpperBound.push_back(99.0);			// x0 -> p0
	mLowerBound.push_back(-99.0);					mUpperBound.push_back(99.0);			// x1 -> p1
#else
	mLowerBound.push_back(0.0);						mUpperBound.push_back(1.0);				// p0+p1
	mLowerBound.push_back(0.0);						mUpperBound.push_back(1.0);				// p0/(p0+p1)
#endif
	if(aNumVariables >= 5)
	{
		mLowerBound.push_back(1.0);					mUpperBound.push_back(999.0);			// w2
	}
}


void BranchSiteModel::initFromTree(void)
{
	// Initialize branch lengths from the phylo tree
	mForest.setTimesFromLengths(mVar);

	// Ask for initialization completion
	mInitType = INIT_TYPE_TIMES;
}


void BranchSiteModel::initFromTreeAndParams(void)
{
	// Initialize branch lengths from the phylo tree
	mForest.setTimesFromLengths(mVar);

	// Get the parameters
	ParseParameters* params = ParseParameters::getInstance();

	// Initialization as in CodeML (seems)
	mVar[mNumTimes+0] = params->getParameter("w0");							// w0
	mVar[mNumTimes+1] = params->getParameter("k");							// k

	double p0 = params->getParameter("p0");
	double p1 = params->getParameter("p1");
#ifdef USE_ORIGINAL_PROPORTIONS
	if(p0 <= 0 || p1 <= 0) throw FastCodeMLFatal("Invalid p0 and p1 values");
	mVar[mNumTimes+2] = log(p0);											// x0 -> p0
	mVar[mNumTimes+3] = log(p1);											// x1 -> p1
#else
	if(p0 < 0 || p1 < 0 || (p0+p1) < 1e-15) throw FastCodeMLFatal("Invalid p0 and p1 values");
	mVar[mNumTimes+2] = p0+p1;												// p0+p1
	mVar[mNumTimes+3] = p0/(p0+p1);											// p0/(p0+p1)
#endif
	if(mNumVariables == 5 && mInitType != INIT_TYPE_RES_5)
	{
		mVar[mNumTimes+4] = params->getParameter("w2");						// w2

		// Ask for initialization completion
		mInitType = INIT_TYPE_RES_5;
	}
	else
	{
		// Ask for initialization completion
		mInitType = INIT_TYPE_RES_4;
	}
}

void BranchSiteModel::initFromResult(const std::vector<double>& aPreviousResult, unsigned int aValidLen)
{
	// Adjust the length to be copied
	if(aValidLen == 0) aValidLen = static_cast<unsigned int>(aPreviousResult.size());

	// Too long, cut. Too short, ignore. 
	if(aValidLen > mNumTimes+mNumVariables) aValidLen = mNumTimes+mNumVariables;
	else if(aValidLen < mNumTimes)
	{
		mInitType = INIT_TYPE_NONE;
		return;
	}
	else if(aValidLen < mNumTimes+4) aValidLen = mNumTimes;

	// Copy the requested values
	mVar.assign(aPreviousResult.begin(), aPreviousResult.begin()+aValidLen);
	mVar.resize(mNumTimes+mNumVariables);

	// Ask for initialization completion
	if(aValidLen == mNumTimes)        mInitType = INIT_TYPE_TIMES;
	else if(aValidLen == mNumTimes+4) mInitType = INIT_TYPE_RES_4;
	else                              mInitType = INIT_TYPE_RES_5;
}


void BranchSiteModel::initVariables(void)
{
	unsigned int i;

	// Initialize time
	if(mInitType == INIT_TYPE_NONE)
	{
		for(i=0; i < mNumTimes; ++i) mVar[i] = randFrom0to1()*.1+0.01;	// T
	}

	// Initialize w0, k, v1, v2
	if(mInitType == INIT_TYPE_TIMES || mInitType == INIT_TYPE_NONE)
	{
		mVar[mNumTimes+0] = randFrom0to1()*0.8 + 0.1;					// w0
		mVar[mNumTimes+1] = 2.0;										// k
#ifdef USE_ORIGINAL_PROPORTIONS
		mVar[mNumTimes+2] = 1.0 + 0.2 * randFrom0to1();					// x0 -> p0
		mVar[mNumTimes+3] = 0.2*randFrom0to1();							// x1 -> p1
#else
		mVar[mNumTimes+2] = randFrom0to1();								// p0+p1
		mVar[mNumTimes+3] = randFrom0to1();								// p0/(p0+p1)
#endif
	}

	// Initialize w2 if needed
	if(mNumVariables == 5 && mInitType != INIT_TYPE_RES_5)
	{
		mVar[mNumTimes+4] = 1.001 + 0.149 * randFrom0to1();				// w2
	}

	// Re-initialize the next time
	mInitType = INIT_TYPE_NONE;

	// Check the initial values to be inside the domain (otherwise clamp them to the domain)
	unsigned int nv = mNumTimes+mNumVariables;
	for(i=0; i < nv; ++i)
	{
		if(mVar[i] < mLowerBound[i]) mVar[i] = mLowerBound[i];
		if(mVar[i] > mUpperBound[i]) mVar[i] = mUpperBound[i];
	}
}


double BranchSiteModelNullHyp::operator()(size_t aFgBranch)
{
	// Initialize the variables to be optimized
	initVariables();

	// Initialize the variables used to avoid unneeded recomputing
	mPrevK      = DBL_MAX;
	mPrevOmega0 = DBL_MAX;

	// Initialize the matrix set
	mSet.initForH0(mForest.adjustFgBranchIdx(aFgBranch));

	// Run the optimizer
	return maximizeLikelihood(aFgBranch);
}


double BranchSiteModelAltHyp::operator()(size_t aFgBranch)
{
	// Initialize the variables to be optimized
	initVariables();

	// Initialize the variables used to avoid unneeded recomputing
	mPrevK      = DBL_MAX;
	mPrevOmega0 = DBL_MAX;
	mPrevOmega2 = DBL_MAX;

	// Initialize the matrix set
	mSet.initForH1(mForest.adjustFgBranchIdx(aFgBranch));

	// Run the optimizer
	return maximizeLikelihood(aFgBranch);
}


double BranchSiteModelNullHyp::computeLikelihood(const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Check if steps can be skipped
	const bool changed_w0 = isDifferent(aVar[mNumTimes+0], mPrevOmega0);
	const bool changed_k  = isDifferent(aVar[mNumTimes+1], mPrevK);
	if(changed_w0) mPrevOmega0 = aVar[mNumTimes+0];
	if(changed_k)  mPrevK      = aVar[mNumTimes+1];

	// Fill the matrices and compute their eigen decomposition. Not worth the effort to parallelize this section.
	if(changed_w0 || changed_k)
	{
		mScaleQw0 = mQw0.fillMatrix(aVar[mNumTimes+0], aVar[mNumTimes+1]);
		mQw0.eigenQREV();
	}
	if(changed_k)
	{
		mScaleQ1  = mQ1.fillMatrix(                    aVar[mNumTimes+1]);
		mQ1.eigenQREV();
	}

	// Compute all proportions
	getProportions(aVar[mNumTimes+2], aVar[mNumTimes+3], mProportions);

	// Compute the scale values
	const double fg_scale = mProportions[0]*mScaleQw0 +
							mProportions[1]*mScaleQ1  +
							mProportions[2]*mScaleQ1  +
							mProportions[3]*mScaleQ1;

	const double bg_scale = (mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1)/(mProportions[0]+mProportions[1]);

	if(mExtraDebug > 0)
	{
		std::cerr << "FG: " << std::setprecision(8) << fg_scale << " BG: " << bg_scale << std::endl;
		std::cerr << "The following is the value printed by CodeML" << std::endl;
		std::cerr << "FG: " << std::setprecision(8) << 1./fg_scale << " BG: " << 1./bg_scale << std::endl;
		std::cerr << "Q0 " << mScaleQw0 << std::endl;
		std::cerr << "Q1 " << mScaleQ1 << std::endl << std::endl;
	}

	// Fill the Transition Matrix Sets
	mSet.computeMatrixSetH0(mQw0, mQ1, changed_w0 || changed_k, bg_scale, fg_scale, aVar);

	// Compute likelihoods
	mForest.computeLikelihoods(mSet, mLikelihoods, 0);

	// Precompute the proportions to be used
	const double p0 = mProportions[0];
	const double p1_p2b = mProportions[1]+mProportions[3];
	const double p2a = mProportions[2];

	// For all (valid) sites. Don't parallelize: time increases and results are errant
	const size_t num_sites = mForest.getNumSites();
	const std::vector<double>& mult = mForest.getSiteMultiplicity();
	double lnl = 0;
	double scale = 0;
	for(size_t site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		// double p = mProportions[0]*mLikelihoods[0*num_sites+site] +
		//		     (mProportions[1]+mProportions[3])*mLikelihoods[1*num_sites+site] +
		//		      mProportions[2]*mLikelihoods[2*num_sites+site];
		double p = mLikelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= p0;
		double x = mLikelihoods[1*num_sites+site];
		if(x > 0) p += p1_p2b*x;
		x = mLikelihoods[2*num_sites+site];
		if(x > 0) p += p2a*x;

		//x = (p > 0) ? log(p)-(mForest.getNumBranches()-mForest.getNumInternalBranches())*log(GLOBAL_SCALING_FACTOR) : mMaxLnL-100000;
		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];
		scale += mult[site]*(mForest.getNumBranches()-mForest.getNumInternalBranches());

		if(mExtraDebug > 1)
		{
			std::cerr << std::setw(4) << site << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[0*num_sites+site] << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[1*num_sites+site] << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[2*num_sites+site] << " -> ";
			std::cerr << std::fixed << std::setw(14) << x*mult[site] << std::endl;
		}
	}
	lnl -= scale*log(GLOBAL_SCALING_FACTOR);

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
	}

	return lnl;
}

	
double BranchSiteModelAltHyp::computeLikelihood(const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Check if steps can be skipped
	const bool changed_w0 = isDifferent(aVar[mNumTimes+0], mPrevOmega0);
	const bool changed_w2 = isDifferent(aVar[mNumTimes+4], mPrevOmega2);
	const bool changed_k  = isDifferent(aVar[mNumTimes+1], mPrevK);
	if(changed_w0) mPrevOmega0 = aVar[mNumTimes+0];
	if(changed_w2) mPrevOmega2 = aVar[mNumTimes+4];
	if(changed_k)  mPrevK      = aVar[mNumTimes+1];

	// Fill the matrices and compute their eigen decomposition. Not worth the effort to parallelize this section.
	if(changed_w0 || changed_k)
	{
		mScaleQw0 = mQw0.fillMatrix(aVar[mNumTimes+0], aVar[mNumTimes+1]);
		mQw0.eigenQREV();
	}
	if(changed_w2 || changed_k)
	{
		mScaleQw2 = mQw2.fillMatrix(aVar[mNumTimes+4], aVar[mNumTimes+1]);
		mQw2.eigenQREV();
	}
	if(changed_k)
	{
		mScaleQ1  = mQ1.fillMatrix(                    aVar[mNumTimes+1]);
		mQ1.eigenQREV();
	}

	// Compute all proportions
	getProportions(aVar[mNumTimes+2], aVar[mNumTimes+3], mProportions);

	// Compute the scale values
	const double fg_scale = mProportions[0]*mScaleQw0 +
							mProportions[1]*mScaleQ1  +
							mProportions[2]*mScaleQw2 +
							mProportions[3]*mScaleQw2;

	const double bg_scale = (mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1)/(mProportions[0]+mProportions[1]);

	if(mExtraDebug > 0)
	{
		std::cerr << "FG: " << std::setprecision(8) << fg_scale << " BG: " << bg_scale << std::endl;
		std::cerr << "The following is the value printed by CodeML" << std::endl;
		std::cerr << "FG: " << std::setprecision(8) << 1./fg_scale << " BG: " << 1./bg_scale << std::endl;
		std::cerr << "Q0 " << mScaleQw0 << std::endl;
		std::cerr << "Q1 " << mScaleQ1 << std::endl;
		std::cerr << "Q2 " << mScaleQw2 << std::endl << std::endl;
	}

	// Fill the Transition Matrix Sets
	mSet.computeMatrixSetH1(mQw0, mQ1, mQw2, changed_w2 || changed_k, bg_scale, fg_scale, aVar);

	// Compute likelihoods
	mForest.computeLikelihoods(mSet, mLikelihoods, 1);

	// Precompute the proportions to be used
	const double p0  = mProportions[0];
	const double p1  = mProportions[1];
	const double p2a = mProportions[2];
	const double p2b = mProportions[3];

	// For all (valid) sites. Don't parallelize: time increase and the results are errant
	const size_t num_sites = mForest.getNumSites();
	const std::vector<double>& mult = mForest.getSiteMultiplicity();
	double lnl = 0;
	double scale = 0;
	for(size_t site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		//double p = mProportions[0]*mLikelihoods[0*num_sites+site] +
		//		     mProportions[1]*mLikelihoods[1*num_sites+site] +
		//		     mProportions[2]*mLikelihoods[2*num_sites+site] +
		//		     mProportions[3]*mLikelihoods[3*num_sites+site];
		//
		double p = mLikelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= p0;
		double x = mLikelihoods[1*num_sites+site];
		if(x > 0) p += p1*x;
		x = mLikelihoods[2*num_sites+site];
		if(x > 0) p += p2a*x;
		x = mLikelihoods[3*num_sites+site];
		if(x > 0) p += p2b*x;

		//x = (p > 0) ? log(p)-(mForest.getNumBranches()-mForest.getNumInternalBranches())*log(GLOBAL_SCALING_FACTOR) : mMaxLnL-100000;
		//lnl += x*mult[site];
		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];
		scale += mult[site]*(mForest.getNumBranches()-mForest.getNumInternalBranches());

		if(mExtraDebug > 1)
		{
			std::cerr << std::setw(4) << site << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[0*num_sites+site] << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[1*num_sites+site] << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[2*num_sites+site] << ' ';
			std::cerr << std::scientific << std::setw(14) << mLikelihoods[3*num_sites+site] << " -> ";
			std::cerr << std::fixed << std::setw(14) << x*mult[site] << std::endl;
		}
	}
	lnl -= scale*log(GLOBAL_SCALING_FACTOR);

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
	/// @param[in] aTrace If set the optimization progress is traced
	/// @param[in] aUpper Upper limit for the variables (to constrain the gradient computation)
	/// @param[in] aDeltaForGradient The variable increment to compute gradient
	///
	MaximizerFunction(BranchSiteModel* aModel,
					  bool aTrace,
					  const std::vector<double>& aUpper,
					  double aDeltaForGradient)
					: mModel(aModel), mTrace(aTrace), mUpper(aUpper), mDeltaForGradient(aDeltaForGradient) {}

	/// Functor.
	/// It computes the function and the gradient if needed.
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	/// @return The evaluated function
	///
	double operator()(const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		double f0 = mModel->computeLikelihood(aVars, mTrace);

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
		const size_t vs = aVars.size();
		for(size_t i=0; i < vs; ++i)
		{
			const double v = aVars[i];

#ifdef USE_ORIGINAL_PROPORTIONS
			double eh = mDeltaForGradient * (fabs(v)+1.);
#else
			double eh = mDeltaForGradient * (v+1.);
#endif

			x[i] += eh;
			if(x[i] >= mUpper[i]) {x[i] -= 2*eh; eh = -eh;}

			const double f1 = mModel->computeLikelihood(x, false);

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
	BranchSiteModel*	mModel;				///< Pointer to the model to be evaluated
	bool				mTrace;				///< If set traces the optimization progresses
	std::vector<double>	mUpper;				///< Upper limit of the variables to constrain the interval on which the gradient should be computed
	double				mDeltaForGradient;	///< The variable increment to compute gradient
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
	if(mOnlyInitialStep) return computeLikelihood(mVar, mTrace);

	// Special case for the CodeML optimizer
	if(mOptAlgo == OPTIM_LD_MING2)
	{
		try
		{
			// Create the optimizer
			Ming2 optim(this, mTrace, mLowerBound, mUpperBound, mDeltaForGradient);

			// Do the maximization
			double maxl = optim.minimizeFunction(mVar);
			
			if(mTrace)
			{
				std::cerr << std::endl << "Function invocations:       " << mNumEvaluations << std::endl;
				std::cerr <<              "Final log-likelihood value: " << maxl << std::endl;
				printVar(mVar);
			}
			return maxl;
		}
		catch(std::exception& e)
		{
			std::ostringstream o;
			o << "Exception in computation: " << e.what() << std::endl;
			throw FastCodeMLFatal(o);
		}
	}

	// Select the maximizer algorithm (the listed ones works and are reasonably fast for FastCodeML)
	std::auto_ptr<nlopt::opt> opt;
	switch(mOptAlgo)
	{
	case OPTIM_LD_LBFGS:
		opt.reset(new nlopt::opt(nlopt::LD_LBFGS,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LD_VAR1:
		opt.reset(new nlopt::opt(nlopt::LD_VAR1,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LD_VAR2:
		opt.reset(new nlopt::opt(nlopt::LD_VAR2,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LD_SLSQP:
		opt.reset(new nlopt::opt(nlopt::LD_SLSQP,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LN_BOBYQA:
		opt.reset(new nlopt::opt(nlopt::LN_BOBYQA,  mNumTimes+mNumVariables));
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

	default:
		throw FastCodeMLFatal("Invalid optimization algorithm identifier on the command line");
	}

	// Initialize bounds and termination criteria
	try
	{
		opt->set_lower_bounds(mLowerBound);
		opt->set_upper_bounds(mUpperBound);
    	opt->set_ftol_abs(1e-4);
		nlopt::srand(static_cast<unsigned long>(mSeed));
	}
	catch(std::exception& e)
	{
		std::ostringstream o;
		o << "Exception during inizialization: " << e.what() << std::endl;
		throw FastCodeMLFatal(o);
	}

	// Optimize the function
	double maxl = 0;
	try
	{
		MaximizerFunction compute(this, mTrace, mUpperBound, mDeltaForGradient);

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
		std::ostringstream o;
		o << "Exception in computation: " << e.what() << std::endl;
		throw FastCodeMLFatal(o);
	}

	return maxl;
}

