
#include "OptSESOP.h"
#include "MathSupport.h"
#include "lapack.h"
#include "nlopt.hpp"

double OptSESOP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	alocateMemory();
	
	double maxl = -1000000;
	int success = SESOPminimizer(&maxl, &aVars[0]);
	return -maxl;
}



// ----------------------------------------------------------------------
int OptSESOP::SESOPminimizer(double *f, double *x)
{	
	// Store the initial state
	
	std::vector<double> x_0, x_prev;
	x_0.resize(mN);
	x_prev.resize(mN);
	memcpy(&x_0[0], x, size_vect);
	memcpy(&x_prev[0], x, size_vect);
	
	*f = -mModel->computeLikelihood(x_0, mTrace);
	double f_prev(*f);
	
	
	// nLopt optimizer parameters
	std::auto_ptr<nlopt::opt> opt;
	
	// initialize the matrix D
	computeGradient(*f, x, mGradient);
	mM = 1;
	memcpy(md2, mGradient, size_vect);
		
	// loop until convergene reached
	bool convergence_reached( false );
	while(!convergence_reached)
	{
		if(mVerbose > 2)
			std::cout << "Starting step " << mStep << ":\n";
		// update D
		if(mStep > 0)
		{	
			// compute and store all the required vectors to construct the subspace
			
			// -- gradient related variables
			// mGradient_prev
			saveGradient(mGradient);
			// mGradient
			computeGradient(*f, x, mGradient);
			// md2
			updateOmega();
			daxpy_(&mN, &mOmega, mGradient, &I1, md2, &I1);
			
			// -- md1
			memcpy(md1, x, size_vect);
			daxpy_(&mN, &minus_one, &x_0[0], &I1, md1, &I1);
			
			// -- previous direction
			memcpy(&mSpace[0], x, size_vect);
			daxpy_(&mN, &minus_one, &x_prev[0], &I1, &mSpace[0], &I1);
			saveDirection(&mSpace[0]);			
		}
		updateDMatrix();
		
		memcpy(&x_prev[0], x, size_vect); 
		
		if(mVerbose > 2)
			std::cout << "Matrix D updated, optimizing...\n";
		
		// optimize in the subspace, i.e. find best alpha: f(x+D*alpha) minimized
		memcpy(&mSpace[0], x, size_vect); // we store the current state in the workspace
		
		opt.reset(new nlopt::opt(nlopt::LN_BOBYQA, mM));
		
		// TODO lower and upper bounds/use COBYLA or linear inequalty supported optimizer
		
		opt->set_lower_bounds(mLowerBoundSubspace);
		opt->set_upper_bounds(mUpperBoundSubspace);
		opt->set_ftol_rel(1e-6);
		//nlopt::srand(static_cast<unsigned long>(mSeed));
		opt->set_max_objective(OptSESOP::subspaceEvaluatorWrapper, this);
		
		try
		{
			nlopt::result result = opt->optimize(alpha, *f);
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
		
		
		// update current state
		char trans = 'N';
		dgemv_(&trans, &mN, &mM, &D1, mD, &mM, &alpha[0], &I1, &D1, x, &I1);
		
		if(mVerbose > 2)
			std::cout << "Optimized at this step, obtained likelihood: " << *f << "\n";
		
		// check convergence
		convergence_reached = fabs(*f - f_prev) < mRelativeError;
		if(mStep > mMaxIterations)
			return 1;
		
		f_prev = *f;
		mStep ++;
	}
	return 0;
}

// ----------------------------------------------------------------------
double OptSESOP::operator()(const std::vector<double>& aVarsAlpha, std::vector<double>& aGrad)
{
	std::vector<double> x(mN);
	// TODO!!!!!!! prealocate
	memcpy(&x[0], &mSpace[0], size_vect);
		
	// compute x = x + D*alpha
	char trans = 'N';
	dgemv_(&trans, &mN, &mM, &D1, mD, &mM, &aVarsAlpha[0], &I1, &D1, &x[0], &I1);
	
	return - mModel->computeLikelihood(x, mTrace);
}



// ----------------------------------------------------------------------
void OptSESOP::updateOmega() {0.5 + sqrt(0.25 + square(mOmega));}

// ----------------------------------------------------------------------
void OptSESOP::alocateMemory()
{
	size_vect = mN*sizeof(mN);
	
	mSubspaceStorage = mN*(3+s1+s2);
	mSpace.resize(2*mSubspaceStorage+mN);
	
	mGradient = &mSpace[mN];
	md1 = &mSpace[2*mN];
	md2 = &mSpace[3*mN];
	
	mGradient_prev.resize(s1);
	mDirection_prev.resize(s2);
	for(int i(0); i<s1; i++) {mGradient_prev[i] = &mSpace[4*mN+i];}
	for(int i(0); i<s2; i++) {mDirection_prev[i] = &mSpace[(4+s1)*mN+i];}	
	
	
	alpha.reserve(3+s1+s2);
	mM = 1;
	alpha.resize(mM);
	
	mD = &mSpace[mN+mSubspaceStorage];
}

// ----------------------------------------------------------------------
void OptSESOP::saveDirection(double *dir)
{
	if(mStep < s1)
	{
		memcpy(mDirection_prev[mStep], dir, size_vect);
	}
	else
	{
		for(int i(0); i<s1-1; i++)
		{
			mDirection_prev[i] = mDirection_prev[i+1];
		}
		memcpy(mDirection_prev.back(), dir, size_vect);
	}
}

// ----------------------------------------------------------------------
void OptSESOP::saveGradient(double *grad)
{
	if(mStep < s2)
	{
		memcpy(mGradient_prev[mStep], grad, size_vect);
	}
	else
	{
		for(int i(0); i<s2-1; i++)
		{
			mGradient_prev[i] = mGradient_prev[i+1];
		}
		memcpy(mGradient_prev.back(), grad, size_vect);
	}
}
	
// ----------------------------------------------------------------------
void OptSESOP::updateDMatrix()
{
	// local variables
	double norm;
	double inv_norm;
	double *vec;

	// first direction of the subspace: gradient
	memcpy(mD, mGradient, size_vect);
	
	norm = dnrm2_(&mN, mD, &I1);
	inv_norm = 1./norm;
	dscal_(&mN, &inv_norm, mD, &I1);
	
	mM = 1;
	
	// next directions of the subspace: Nemirovski directions
	if(mStep > 1)
	{
		vec = mD+size_vect;
		memcpy(vec, md1, size_vect);
				
		norm = dnrm2_(&mN, vec, &I1);
		inv_norm = 1./norm;
		dscal_(&mN, &inv_norm, vec, &I1);
		
		vec = mD+2*size_vect;
		memcpy(vec, md2, size_vect);
		
		norm = dnrm2_(&N, vec, &I1);
		inv_norm = 1./norm;
		dscal_(&mN, &inv_norm, vec, &I1);
		
		mM = 3;
	}
	// last directions of the subspace: previous gradients and directions
	if(mStep > max2(s1,s2))
	{
		for(int i(0); i<s1; i++)
		{
			vec = mD+(3+i)*size_vect;
			memcpy(vec, mDirection_prev[i], size_vect);
			norm = dnrm2_(&N, vec, &I1);
			inv_norm = 1./norm;
			dscal_(&mN, &inv_norm, vec, &I1);
		}
		for(int i(0); i<s2; i++)
		{
			vec = mD+(3+s1+i)*size_vect;
			memcpy(vec, mGradient_prev[i], size_vect);
			norm = dnrm2_(&N, vec, &I1);
			inv_norm = 1./norm;
			dscal_(&mN, &inv_norm, vec, &I1);
		}
		mM = 3+s1+s2;
	}
}

// ----------------------------------------------------------------------
void OptSESOP::computeGradient(double aPointValue, const double *aVars, double* aGrad) const
{
	std::vector<double> vars_working_copy(mN);
	memcpy(&vars_working_copy[0], aVars, size_vect);
	std::vector<double> delta(mN);

	size_t i = 0;
	for(; i < mN; ++i)
	{
		double eh = 1e-6 * (vars_working_copy[i]+1.);
			
		// If it is going over the upper limit reverse the delta
		vars_working_copy[i] += eh;
		if(vars_working_copy[i] >= mUpperBound[i])
		{
			vars_working_copy[i] -= 2*eh;
			delta[i] = -eh;
		}
		else
		{
			delta[i] = eh;
		}
		
		const double f1 = -mModel->computeLikelihood(vars_working_copy, false);
		aGrad[i] = (f1-aPointValue)/delta[i];
		
		vars_working_copy[i] = aVars[i];
	}
}


