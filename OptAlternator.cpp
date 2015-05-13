
#include "OptAlternator.h"
#include "nlopt.hpp"


// ----------------------------------------------------------------------
double OptAlternator::maximizeFunction(std::vector<double>& aVars)
{
	// alocate space
	
	mN = static_cast<int>(aVars.size());
	mNextra = mN - mNumTimes;
	
	size_vect = mN*sizeof(double);
	
	x_.resize(mN);
	mSpace.resize(mN);
	
	
	double maxLnL( -1e14 );
	AlternatorMinimizer(&maxLnL, &aVars[0]);
	return - maxLnL;
}


// ----------------------------------------------------------------------
void OptAlternator::AlternatorMinimizer(double *f, double *x)
{
	memcpy(&x_[0], x, size_vect);
	
	mSearchState = SPACE_FULL;
	int sizeSpace = mN;
	int idFirstVar, idLastVar;
	
	// nLopt optimizer
	std::auto_ptr<nlopt::opt> opt;
	
	// vector to pass to the nLopt optimizer
	std::vector<double> vars;
	std::vector<double> upBounds;
	std::vector<double> loBounds;
	vars.reserve(mN);
	upBounds.reserve(mN);
	loBounds.reserve(mN);
	
	
	// relative error of the step
	double absoluteTempError;
	double f_prev, f_curr;
	f_curr = *f;
	bool stopCriterionReached( false );
	
	mSearchState = SPACE_EXTRA_ONLY;
	
	while(!stopCriterionReached)
	{
		absoluteTempError = 0.5/(1.0+exp(mStep));
		++mStep;
		f_prev = f_curr;
		
		std::cout << "Starting step " << mStep << std::endl;
		
		// update space variables
		getSpaceProperties(idFirstVar, idLastVar);
		sizeSpace = idLastVar-idFirstVar;
		
		vars.resize(sizeSpace);
		upBounds.resize(sizeSpace);
		loBounds.resize(sizeSpace);
		memcpy(&vars[0], &x[idFirstVar], sizeSpace*sizeof(double));
		memcpy(&loBounds[0], &mLowerBound[idFirstVar], sizeSpace*sizeof(double));
		memcpy(&upBounds[0], &mUpperBound[idFirstVar], sizeSpace*sizeof(double));
		
		// setup optimizer
		opt.reset(new nlopt::opt(nlopt::LD_SLSQP, sizeSpace));
		//opt.reset(new nlopt::opt(nlopt::LD_LBFGS, sizeSpace));
		//opt.reset(new nlopt::opt(nlopt::LN_BOBYQA, sizeSpace));
		opt->set_vector_storage(20);
		
		opt->set_ftol_abs(absoluteTempError);
		
		opt->set_lower_bounds(loBounds);
		opt->set_upper_bounds(upBounds);		
		
		
		opt->set_max_objective(OptAlternator::subspaceEvaluatorWrapper, this);
		
		// optimize
		try
		{
			nlopt::result result = opt->optimize(vars, f_curr);
		}
		catch(const nlopt::forced_stop&)
		{
			if(mTrace) std::cout << "Optimization stopped because LRT not satisfied" << std::endl;
		}
		catch(const nlopt::roundoff_limited&)
		{
			throw FastCodeMLFatal("Exception in computation: Halted because roundoff errors limited progress, equivalent to NLOPT_ROUNDOFF_LIMITED.");
		}
		catch(const std::runtime_error&)
		{
			throw FastCodeMLFatal("Exception in computation: Generic failure, equivalent to NLOPT_FAILURE.");
		}
		catch(const std::invalid_argument&)
		{
			throw FastCodeMLFatal("Exception in computation: Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera), equivalent to NLOPT_INVALID_ARGS.");
		}
		catch(const std::bad_alloc&)
		{
			throw FastCodeMLFatal("Exception in computation: Ran out of memory (a memory allocation failed), equivalent to NLOPT_OUT_OF_MEMORY.");
		}
		catch(const std::exception& e)
		{
			std::ostringstream o;
			o << "Exception in computation: " << e.what();
			throw FastCodeMLFatal(o);
		}
		
		// change the subspace
		switch(mSearchState)
		{
			case SPACE_FULL:
				mSearchState = SPACE_BRANCHES_ONLY;
				break;
			
			case SPACE_BRANCHES_ONLY:
				mSearchState = SPACE_FULL;
				break;
		
			case SPACE_EXTRA_ONLY:
				mSearchState = SPACE_FULL;
				break;		
		}
		
		// verify convergence
		std::cout << " df : " << fabs(f_prev-f_curr) << std::endl;
		stopCriterionReached = fabs(f_prev-f_curr) < mAbsoluteError;
	
	}
	*f = f_curr;
}


// ----------------------------------------------------------------------
double OptAlternator::evaluateFunction(const double *aVar, double *aGrad)
{
	// select the good space
	int imin, i, imax, k; 
	
	getSpaceProperties(imin, imax);
	i = imin;
	
	memcpy(&x_[imin], aVar, (imax-imin)*sizeof(double));
	double f = mModel->computeLikelihood(x_, mTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if(mStopIfBigger && f >= mThreshold) throw nlopt::forced_stop();
	
	// gradient computation
	if(aGrad)
	{	
		volatile double eh;
		double fp;
		
		for(k=0; i<imax; ++k, ++i)
		{
			eh = 1e-6 * ( (aVar[i]>1.0) ? aVar[i] : 1.0 );
			
			if( aVar[k] + eh > mUpperBound[i] )
				eh = -eh;
				
			x_[i] += eh;
			eh = x_[i]-aVar[k];
			
			fp = mModel->computeLikelihood(x_, false);
			aGrad[k] = (fp-f) / eh;
			x_[i] = aVar[k];
		}
	}
		
	return f;
}


// ----------------------------------------------------------------------
void OptAlternator::getSpaceProperties(int& idFirstVar, int& idLastVar) const
{
	switch(mSearchState)
	{
		case SPACE_FULL:
			idFirstVar	= 0;
			idLastVar	= mN;
			break;
			
		case SPACE_BRANCHES_ONLY:
			idFirstVar	= 0;
			idLastVar	= mNumTimes;
			break;
		
		case SPACE_EXTRA_ONLY:
			idFirstVar	= mNumTimes;
			idLastVar	= mN;
			break;		
	}
}


