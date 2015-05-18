
#include "OptNES.h"
#include "MathSupport.h"
#include "lapack.h"
#include <algorithm>
#include <iomanip>

// ----------------------------------------------------------------------
//	Class members definition
// ----------------------------------------------------------------------
double OptNES::maximizeFunction(std::vector<double>& aVars, unsigned int popSize)
{
	mN = static_cast<int>(aVars.size());
	mLambda = popSize;
	
	if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "Starting OptNES with " << mN << " variables and population of size " << mLambda << std::endl;
	
	alocateMemory();
	
	
	double maxl = 1000000;
	/*int success = */ NESminimizer(&maxl, &aVars[0]);
	return maxl;
}

// ----------------------------------------------------------------------
int OptNES::NESminimizer(double *f, double *x)
{
	// initialize distribution parameters
	initializeDistribution();
	
	// initialize the utilities
	
	double scale( 0.0 ), extra_term( 1.0/double(mLambda) );
	
	#pragma omp parallel for reduction(+:scale)
	for(int k=0; k<mLambda; ++k)
	{
		double kthValue = log(double(mLambda)*0.5 + 1.0) - log(double(k+1));
		kthValue = (kthValue > 0.) ? kthValue : 0.0 ;
		
		mUtility[k] = kthValue;
		scale += kthValue;
	}
	scale = 1./scale;
	
	#pragma omp parallel for
	for(int k=0; k<mLambda; ++k)
	{
		mUtility[k] = mUtility[k]*scale - extra_term;
	}	
	
	std::cout << "Utility:\n";
	for(int k=0; k<mLambda; ++k)
		std::cout << std::setprecision(10) << mUtility[k] << std::endl;
	
	
	// main loop
	
	for(mStep = 0; mStep < 10; ++mStep)
	{
	
		if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		{
			std::cout << "Parameters:\n--Mu:\n";
			for(int i=0; i<mN; ++i)
				std::cout << std::setprecision(10) << mMu[i] << std::endl;
			std::cout << "\n--Sigma:\n";
			for(int i=0; i<mN; ++i)
				std::cout << std::setprecision(10) << mSigma[i] << std::endl;
		}
		
		generatePopulation();
		sortFitness();
		
		// save current best result
		int bestId = mPermutation[0];
		
		if(*f > mPopFitness[bestId])
		{
			*f = mPopFitness[bestId];
			memcpy(x, &mPopPos[bestId*mN], size_vect);
		}
		
		computeGradients();
		
		updateParameters();
	}
	return 0;
}

// ----------------------------------------------------------------------
void OptNES::alocateMemory(void)
{
	size_vect = mN*sizeof(double);
	
	x_.resize(mN);	
	mSpace.resize( mN + 4*mN + mLambda*(2*mN+2) );
	
	mMu			= &mSpace[mN];
	mSigma		= mMu + mN;
	mGradMu		= mSigma + mN;
	mGradSigma	= mGradMu + mN;
	
	mPopFitness	= mGradSigma + mN;
	mUtility	= mPopFitness + mLambda;
	mPopPos		= mUtility + mLambda;
	mS			= mPopPos + (mLambda*mN);
	
	mPermutation.resize(mLambda);
}


// ----------------------------------------------------------------------
void OptNES::initializeDistribution(void)
{
	double shape, rate, a, b;
	
	double mf = 3.0;
	double sf = 10.0;
	
	// branch lengths:
	shape = 0.5031126;
	rate  = 0.1844347;
	#pragma omp parallel for
	for(int k=0; k<mNumTimes; ++k)
	{
		mMu[k] 		= mf * shape*rate;
		mSigma[k]	= sf * shape*rate*rate;
	}
	
	// v0
	rate = 9.441686;
	mMu[mNumTimes+0] 	= 1.0 - mf / rate;
	mSigma[mNumTimes+0] = sf / (rate*rate);
	
	// v1
	shape = 0.8764469;
	rate  = 0.126567;
	mMu[mNumTimes+1]	= 1. - mf * shape*rate;
	mSigma[mNumTimes+1] = sf * shape*rate*rate;
	
	// w0
	a = 1.638631;
	b = 21.841174;
	mMu[mNumTimes+2] 	= mf * a / (a + b);
	mSigma[mNumTimes+2] = sf * a*b / ( square(a+b)*(a+b+1.0) );
						
	// kappa	
	shape = 7.547445;
	rate  = 0.5789037;
	mMu[mNumTimes+3] 	= shape*rate;
	mSigma[mNumTimes+3] = shape*rate*rate;
	
	// w2
	if(mNumTimes + 4 > mN)
	{
		std::cout<<"w2"<<std::endl;
		shape = 0.209741;
		rate  = 274.5373;
		mMu[mNumTimes+4] 	= 1.0 + shape*rate;
		mSigma[mNumTimes+4] = shape*rate*rate;
	}
}

// ----------------------------------------------------------------------
void OptNES::generatePopulation(void)
{
	// generate randomly the positions
	
	for(size_t k( 0 ); k<mLambda; ++k)
	{
		for(size_t i( 0 ); i<mN; ++i)
		{
			// generate a normal distributed number while it is not in the bounds
			double si, xi;
			do{
				si = Norm( rng );
				xi = mMu[i] + mSigma[i]*si;
			}while( xi < mLowerBound[i] || mUpperBound[i] < xi );
			
			mS[k*mN+i] = si;
			mPopPos[k*mN+i] = xi;
		}
	}
	
	// compute the fitness of each position
	for(size_t k( 0 ); k<mLambda; ++k)
	{
		memcpy(&x_[0], &mPopPos[k*mN], size_vect);
		mPopFitness[k] = -mModel->computeLikelihood(x_, mTrace);
	}
}


// ----------------------------------------------------------------------
void OptNES::sortFitness(void)
{	
	#pragma omp parallel for
	for(int k=0; k<mLambda; ++k) {mPermutation[k] = k;}
	
	std::sort(mPermutation.begin(), mPermutation.end(), CompareFitness(mPopFitness));
}

// ----------------------------------------------------------------------
void OptNES::computeGradients(void)
{
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
	{
		double partialValueMu(0.);
		double partialValueSigma(0.);
		//#pragma omp parallel for reduction(+:partialValueMu) reduction(+:partialValueSigma)
		for(size_t k(0); k<mLambda; ++k)
		{
			int idPos = mPermutation[k];
			double sk = mPopPos[idPos*mN + i];
			
			partialValueMu += mUtility[k] * sk;
			partialValueSigma += mUtility[k] * (sk*sk - 1.0);
		}
		mGradMu[i] = partialValueMu;
		mGradSigma[i] = partialValueSigma;
	}
}
	
// ----------------------------------------------------------------------
void OptNES::updateParameters(void)
{
	// learning rates
	
	double eta_mu = 2.0;
	double half_eta_sigma = 0.5;
	
	// effect of the natural gradient
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
	{
		mMu[i] += eta_mu * mSigma[i]*mGradMu[i];
		mSigma[i] *= exp( half_eta_sigma * mGradSigma[i] );		
	}
	
	// effect of the bounds: require (mu +/- alpha*sigma) in the bounds
	const double alpha = 3.;
	#pragma omp parallel for
	for(int i=0; i<mN; ++i)
	{
		double sd = mSigma[i];
		double mn = mMu[i];
		double lb = mLowerBound[i];
		double ub = mUpperBound[i];
		
		if(sd > 0.2*(ub-lb))
			sd = (ub-lb)/10.0;
		
		// mean should be in the bounds
		if( mn < lb )
			mMu[i] = lb + sd;
		mn = mMu[i];
		if( mn > ub )
			mMu[i] = lb - sd;
		mn = mMu[i];
		
		// avoid too large standard deviations that would generate too
		// many variables outside the bounds
		if( (mn - alpha*sd) < lb )
			mSigma[i] = (mn-lb)/alpha;
		sd = mSigma[i];
		if( (mn + alpha*sd) > ub )
			mSigma[i] = (ub-mn)/alpha;
	}
}



