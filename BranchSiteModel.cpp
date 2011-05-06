
#include <iostream>
#include <iomanip>
#include <cfloat>

#ifdef USE_OPTIMIZER
#include "nlopt.hpp"
#endif

/// Uncomment to use the original CodeML proportion definition
//#define USE_ORIGINAL_PROPORTIONS

/// Uncomment if you want to use the old likelihood method
//#define OLD_LIKELIHOOD

#include "BranchSiteModel.h"
#include "MathSupport.h"
#include "Exceptions.h"

/// How much a variable should be changed to compute gradient
static const double SMALL_DIFFERENCE = 1e-6;

void BranchSiteModel::getProportions(double aV0, double aV1, double* aProportions) const
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


double BranchSiteModelNullHyp::computeModel(Forest& aForest, unsigned int aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace)
{
	unsigned int i;

	// Initialize the variables to be optimized
	if(aTimesFromTree)
	{
		aForest.setTimesFromLengths(mVar);
#if 0
		// Initialization as in Octave
		mVar[mNumTimes+0] = 0.062748;											// w0
		mVar[mNumTimes+1] = 2.0;												// k
#else
		// Initialization as in CodeML (seems)
		mVar[mNumTimes+0] = 0.235087;											// w0
		mVar[mNumTimes+1] = 0.4;												// k
#endif

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

	// Allocate the site likelihood array
	mLnLsite.reserve(aForest.getNumSites());
	mLnLsite.resize(aForest.getNumSites());

	// Run the optimizer
	return maximizeLikelihood(aForest, aFgBranch, aOnlyInitialStep, aTrace);
}


double BranchSiteModelAltHyp::computeModel(Forest& aForest, unsigned int aFgBranch, bool aOnlyInitialStep, bool aTimesFromTree, bool aTrace)
{
	unsigned int i;

	// Initialize the variables to be optimized
	if(aTimesFromTree)
	{
		aForest.setTimesFromLengths(mVar);
#if 0
		// Initialization as in Octave
		mVar[mNumTimes+0] = 0.062748;											// w0
		mVar[mNumTimes+1] = 2.0;												// k
		mVar[mNumTimes+4] = 2.5;												// w2
#else
		// Initialization as in CodeML (seems)
		mVar[mNumTimes+0] = 0.235087;											// w0
		mVar[mNumTimes+1] = 0.4;												// k
		mVar[mNumTimes+4] = 1.14833;											// w2
#endif

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
		// Initialize the variables to be optimized
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

	// Allocate the site likelihood array
	mLnLsite.reserve(aForest.getNumSites());
	mLnLsite.resize(aForest.getNumSites());

	// Run the optimizer
	return maximizeLikelihood(aForest, aFgBranch, aOnlyInitialStep, aTrace);
}

	
double BranchSiteModelNullHyp::oneCycleMaximizer(Forest& aForest, unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Compute all proportions
	getProportions(aVar[mNumTimes+2], aVar[mNumTimes+3], mProportions);

	// Compute the matrices
	double scale_qw0 = mQw0.fillQ(aVar[mNumTimes+0], aVar[mNumTimes+1], aForest.getCodonFrequencies());
	double scale_q1  = mQ1.fillQ(                    aVar[mNumTimes+1], aForest.getCodonFrequencies());

	// Compute the scale values
	double fg_scale = mProportions[0]*scale_qw0 +
					  mProportions[1]*scale_q1  +
					  mProportions[2]*scale_q1  +
					  mProportions[3]*scale_q1;

	//double bg_scale = mProportions[0]*scale_qw0 +
	//				  mProportions[1]*scale_q1  +
	//				  mProportions[2]*scale_qw0 +
	//				  mProportions[3]*scale_q1;
	double bg_scale = 1./(mProportions[0]+ mProportions[1])*(mProportions[0]*scale_qw0+mProportions[1]*scale_q1);

	// Compute eigen decomposition of the matrices
	mQw0.eigenQREV(aForest.numGoodCodonFrequencies(), aForest.getSqrtCodonFrequencies(), aForest.getGoodCodonFrequencies());
	mQ1.eigenQREV(aForest.numGoodCodonFrequencies(), aForest.getSqrtCodonFrequencies(), aForest.getGoodCodonFrequencies());

	// Fill the Transition Matrix sets
	mSet.computeMatrixSetH0(mQw0, mQ1, bg_scale, fg_scale, aForest.adjustFgBranchIdx(aFgBranch), aVar);

#ifdef OLD_LIKELIHOOD
	// Compute the likelihood values on the forest
	std::vector<double> like0, like1, like2;
	aForest.computeLikelihood(mSet, 0, like0);
	aForest.computeLikelihood(mSet, 1, like1);
	aForest.computeLikelihood(mSet, 2, like2);

	// For all (valid) sites
	unsigned int num_sites = aForest.getNumSites();
	const double* mult = aForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		double p = mProportions[0]*like0[site] +
                  (mProportions[1]+mProportions[3])*like1[site] +
				   mProportions[2]*like2[site];

		if(p <= 0) return mMaxLnL-100000; // To avoid invalid maxima
		double x = log(p);

		mLnLsite[site] = x;
		lnl += x*mult[site];
	}
#else
	std::vector<double> likelihoods;
	aForest.computeLikelihood(mSet, likelihoods);

	// For all (valid) sites
	unsigned int num_sites = aForest.getNumSites();
	const double* mult = aForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		double p = mProportions[0]*likelihoods[0*num_sites+site] +
				   (mProportions[1]+mProportions[3])*likelihoods[1*num_sites+site] +
				   mProportions[2]*likelihoods[2*num_sites+site];

		if(p <= 0) return mMaxLnL-100000; // To avoid invalid maxima
		double x = log(p);

		mLnLsite[site] = x;
		lnl += x*mult[site];
	}
#endif

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		std::cerr << std::endl << lnl << std::endl;
		printVar(aVar);
	}

	return lnl;
}

	
double BranchSiteModelAltHyp::oneCycleMaximizer(Forest& aForest, unsigned int aFgBranch, const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Compute all proportions
	getProportions(aVar[mNumTimes+2], aVar[mNumTimes+3], mProportions);

	// Compute the matrices
	double scale_qw0 = mQw0.fillQ(aVar[mNumTimes+0], aVar[mNumTimes+1], aForest.getCodonFrequencies());
	double scale_qw2 = mQw2.fillQ(aVar[mNumTimes+4], aVar[mNumTimes+1], aForest.getCodonFrequencies());
	double scale_q1  = mQ1.fillQ(                    aVar[mNumTimes+1], aForest.getCodonFrequencies());

	// Compute the scale values
	double fg_scale = mProportions[0]*scale_qw0 +
					  mProportions[1]*scale_q1  +
					  mProportions[2]*scale_qw2 +
					  mProportions[3]*scale_qw2;

	//double bg_scale = mProportions[0]*scale_qw0 +
	//			      mProportions[1]*scale_q1  +
	//				  mProportions[2]*scale_qw0 +
	//				  mProportions[3]*scale_q1;
	double bg_scale = 1./(mProportions[0]+ mProportions[1])*(mProportions[0]*scale_qw0+mProportions[1]*scale_q1);

	// Compute eigen decomposition of the matrices
	mQw0.eigenQREV(aForest.numGoodCodonFrequencies(), aForest.getSqrtCodonFrequencies(), aForest.getGoodCodonFrequencies());
	mQw2.eigenQREV(aForest.numGoodCodonFrequencies(), aForest.getSqrtCodonFrequencies(), aForest.getGoodCodonFrequencies());
	mQ1.eigenQREV(aForest.numGoodCodonFrequencies(),  aForest.getSqrtCodonFrequencies(), aForest.getGoodCodonFrequencies());

	// Fill the Transition Matrix sets
	mSet.computeMatrixSetH1(mQw0, mQ1, mQw2, bg_scale, fg_scale, aForest.adjustFgBranchIdx(aFgBranch), aVar);

#ifdef OLD_LIKELIHOOD
	// Compute the likelihood values on the forest
	std::vector<double> like0, like1, like2, like3;
	aForest.computeLikelihood(mSet, 0, like0);
	aForest.computeLikelihood(mSet, 1, like1);
	aForest.computeLikelihood(mSet, 2, like2);
	aForest.computeLikelihood(mSet, 3, like3);

	// For all sites
	unsigned int num_sites = aForest.getNumSites();
	const double* mult = aForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		double p = mProportions[0]*like0[site] +
				   mProportions[1]*like1[site] +
				   mProportions[2]*like2[site] +
				   mProportions[3]*like3[site];

		if(p <= 0) return mMaxLnL-100000; // To avoid invalid maxima
		double x = log(p);

		mLnLsite[site] = x;
		lnl += x*mult[site];
	}
#else
	std::vector<double> likelihoods;
	aForest.computeLikelihood(mSet, likelihoods);

	// For all sites
	unsigned int num_sites = aForest.getNumSites();
	const double* mult = aForest.getSiteMultiplicity();
	double lnl = 0;
	for(unsigned int site=0; site < num_sites; ++site)
	{
		double p = mProportions[0]*likelihoods[0*num_sites+site] +
				   mProportions[1]*likelihoods[1*num_sites+site] +
				   mProportions[2]*likelihoods[2*num_sites+site] +
				   mProportions[3]*likelihoods[3*num_sites+site];

		if(p <= 0) return mMaxLnL-100000; // To avoid invalid maxima
		double x = log(p);

		mLnLsite[site] = x;
		lnl += x*mult[site];
	}
#endif

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		std::cerr << std::endl << lnl << std::endl;
		printVar(aVar);
	}

	return lnl;
}

#ifdef USE_OPTIMIZER
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
		for(unsigned int i=0; i < aVars.size(); ++i)
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
#endif

double BranchSiteModel::maximizeLikelihood(Forest& aForest, unsigned int aFgBranch, bool aOnlyInitialStep, bool aTrace)
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

#ifdef USE_OPTIMIZER
	// Initialize the maximum value found and the function evaluations counter
	mMaxLnL = -HUGE_VAL;
	mNumEvaluations = 0;

	// If only the initial step is requested, do it and return
	if(aOnlyInitialStep) return oneCycleMaximizer(aForest, aFgBranch, mVar, aTrace);

	// Select the optimization algorithm
//	nlopt::opt opt(nlopt::GN_DIRECT_L, mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::GN_ISRES,    mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::LN_COBYLA,   mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::LN_BOBYQA,   mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::LN_SBPLX,    mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::G_MLSL_LDS,  mNumTimes+mNumVariables);

	nlopt::opt opt(nlopt::LD_LBFGS,    mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::LD_MMA,      mNumTimes+mNumVariables);
//	nlopt::opt opt(nlopt::LD_SLSQP,    mNumTimes+mNumVariables);

	// Initialize bounds and termination criteria
	try
	{
		opt.set_lower_bounds(mLowerBound);
		opt.set_upper_bounds(mUpperBound);
    	opt.set_ftol_abs(1e-4);
		nlopt::srand(mSeed);

		// If the algorithm requires a local optimizer, then add it
		if(opt.get_algorithm() == nlopt::G_MLSL_LDS)
		{
			// For global optimization put a timeout of one hour
			opt.set_maxtime(60*60);

			nlopt::opt local_opt(nlopt::LN_BOBYQA, mNumTimes+mNumVariables);
			opt.set_local_optimizer(local_opt);
		}
	}
	catch(std::exception& e)
	{
		std::cerr << "Exception during inizialization: " << e.what() << std::endl;
	}

	// Optimize the function
	double maxl = 0;
	try
	{
		MaximizerFunction compute(this, &aForest, aFgBranch, aTrace, mUpperBound);

		opt.set_max_objective(MaximizerFunction::wrap, &compute);

		nlopt::result result = opt.optimize(mVar, maxl);

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
#else
	// If no maximizer available return only the first step result
	mMaxLnL = -HUGE_VAL;
	return oneCycleMaximizer(aForest, aFgBranch, mVar, aTrace);
#endif
}

