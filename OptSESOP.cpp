
#include "OptSESOP.h"
#include "MathSupport.h"
#include "lapack.h"
#include "nlopt.hpp"

double OptSESOP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	if(mVerbose > 2)
	{
		std::cout << "Size of the problem: N=" << mN << "\n";
	}
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
	memcpy(&x_0[0]   , x, size_vect);
	memcpy(&x_prev[0], x, size_vect);	
	
	*f = mModel->computeLikelihood(x_0, mTrace);
	double f_prev = *f - 2.;
	
	// nLopt optimizer
	std::auto_ptr<nlopt::opt> opt;

	// initialize the matrix D
	computeGradient(*f, x, mGradient);
	
	memcpy(mGradient_times, mGradient, mNumTimes*sizeof(double));
	for(int i(mNumTimes); i<mN; i++)
	{
		mGradient_times[i] = 0.;
	}
	
	
	// mGradient_others
	for(size_t i(0); i<mN; i++)
		mGradient_others[i] = 0.;
	//mGradient_others[mNumTimes+3] = mGradient[mNumTimes+3]; // kappa
	mGradient_others[mNumTimes+1] = mGradient[mNumTimes+1];	// v1
	mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
	
	// Initialize the second Nemirovski direction
	memcpy(md2, mGradient, size_vect);
	
		
	// loop until convergene reached
	bool convergence_reached( false );
	while(!convergence_reached)
	{
		double f_diff = fabs(*f - f_prev);
		f_prev = *f;		
		
		if(mVerbose > 2)
			std::cout << "Starting step " << mStep << ":\n";
		// update D
		if(mStep > 0)
		{	
			// compute and store all the required vectors to construct the subspace
			
			// mGradient
			computeGradient(*f, x, mGradient);
			
			// mGradient_times			
			memcpy(mGradient_times, mGradient, mNumTimes*sizeof(double));			
			
			
			// randomly choose the search directions for the matrix parameters:
			for(size_t i(mNumTimes); i<mN; i++)
			{
				mGradient_others[i] = 0.;
			}
			double rand_tmp = 1. - (1. - exp(-f_diff)) * randFrom0to1();
			if(rand_tmp < 0.3)
			{
				mGradient_others[mNumTimes+1] = mGradient[mNumTimes+1];	// v1
				mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
			}
			else if(rand_tmp < 0.5)
			{	
				mGradient_others[mNumTimes+3] = mGradient[mNumTimes+3];	// kappa
				mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
			}
			else
			{
				mGradient_others[mNumTimes+3] = 1.; // kappa
			}
			
			
			/*
			//mGradient_others[mNumTimes+3] = mGradient[mNumTimes+3]; // kappa
			if(f_diff > 5.*mRelativeError)
			{
				mGradient_others[mNumTimes+1] = mGradient[mNumTimes+1];	// v1
				mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
			}
			else
			{
				mGradient_others[mNumTimes+1] = 0.; // v1
				mGradient_others[mNumTimes+2] = 1.; // w0
			}
			*/
			
			
			// md2
			updateOmega();
			daxpy_(&mN, &mOmega, mGradient, &I1, md2, &I1);	
		}
		updateDMatrix();
		
		
		if(mVerbose > 3)
		{
			std::cout << "\n\n --------------- Matrix D: ------------- \n\n";
			for(int j(0); j<mM; j++)
			{
				memcpy(&x_[0], mD+j*mN, size_vect);
				std::cout << "\nColumn " << j << " of the D matrix:\n";
				mModel->printVar(x_, 0);
			}
		}
		
		
		memcpy(&x_prev[0], x, size_vect); 
		
		if(mVerbose > 2)
		{
			std::cout << "Matrix D updated, optimizing...\n";
		}
		
		// optimize in the subspace, i.e. find best alpha s.t. f(x+D*alpha) is minimized
		// we store the current state in the workspace
		memcpy(&mSpace[0], x, size_vect); 
		
#ifdef USE_TRANSFORMED_CONSTRAINTS
		opt.reset(new nlopt::opt(nlopt::LD_SLSQP, mM));
		opt->set_vector_storage(20);
		
		//opt.reset(new nlopt::opt(nlopt::LN_COBYLA, mM));
		
		// add the constraints
		for(int i(0); i<mN; i++)
		{
			data_constraints[2*i  ].sesop = this;
			data_constraints[2*i  ].line = i;
			data_constraints[2*i  ].bound_type = 0;
			data_constraints[2*i+1].sesop = this;
			data_constraints[2*i+1].line = i;
			data_constraints[2*i+1].bound_type = 1;
		}
		for(int constraint_id(0); constraint_id<2*mN; constraint_id++)
		{
			opt->add_inequality_constraint(OptSESOP::myconstraintWrapper, (void*) &data_constraints[constraint_id], 1e-8);
		}
		
		if(mVerbose > 2)
		{
			std::cout << "Inequality constraints added\n";
		}
		
		// initialize alpha to 0 so we start at point x
		for(int i(0); i<mM; i++)
			alpha[i] = 0.;
		
#endif // USE_TRANSFORMED_CONSTRAINTS

#ifdef USE_BOUND_CONSTRAINTS
		opt.reset(new nlopt::opt(nlopt::LN_BOBYQA, mM));
		
		// lower and upper bounds/use COBYLA or linear inequalty supported optimizer
		updateBoundsAndAlpha();
		
		opt->set_lower_bounds(mLowerBoundSubspace);
		opt->set_upper_bounds(mUpperBoundSubspace);
		
		if(mVerbose > 2)
			std::cout << "Bounds set...\n";
#endif // USE_BOUND_CONSTRAINTS

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
		
#ifdef USE_TRANSFORMED_CONSTRAINTS
		if(mVerbose > 2)
			std::cout << "Removing inequality constraints...\n";
		opt->remove_inequality_constraints();
#endif // USE_TRANSFORMED_CONSTRAINTS

		
		
		if(mVerbose > 2)
			std::cout << "Updating current state...\n";
			
		
		// update current state
		char trans = 'N';
		memcpy(x, &x_prev[0], size_vect);
		dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &alpha[0], &I1, &D1, x, &I1);
		
		if(mVerbose > 2)
		{
			std::cout << "Obtained alpha:\n";
			for(int i(0); i<mM; i++)
				std::cout << alpha[i] << " ";
			
			std::cout << "\n\nOptimized at this step, obtained likelihood: ";
			memcpy(&x_[0], x, size_vect);
			mModel->printVar(x_, *f);
			
		}
		
		// test to make the times move...
		
		size_t iter_move = 0;
		size_t num_iter_move = size_t(3. / (1 + square(2.*float(mStep)/float(mN))));
		for(; iter_move < num_iter_move; iter_move ++)
		{
			memcpy(&x_[0], x, size_vect);
			for(size_t i(0); i<mNumTimes; i++)
			{
				//! TEST
				x_[i]+=1e-1*sqrt(randFrom0to1())*mGradient[i]*square(x[i]);				
			
				if(x_[i] < mLowerBound[i])
				{
					x_[i] = 0.;
				}
			}
			double f_test = mModel->computeLikelihood(x_, mTrace);
			if(mVerbose > 2)
				std::cout << "Previous f: " << *f << ", new one: " << f_test << std::endl;
			if (f_test > *f)
			{
				*f = f_test;
				memcpy(x, &x_[0], size_vect);
			}
		}
		
		
		// check convergence
		convergence_reached = f_diff < mRelativeError;
		//					&& square(dnrm2_(&mN, mGradient, &I1))/mN < mRelativeError;
		
		if(mStep > 100)
		{
			if(mVerbose > 2)
			{
				std::cout << "SESOP aborted, number of iterations exceeded.\n";
			}
			return 1;
		}
		
		mStep ++;
	}
		
	return 0;
}

// ----------------------------------------------------------------------
double OptSESOP::operator()(const std::vector<double>& aVarsAlpha, std::vector<double>& aGrad)
{
	memcpy(&x_[0], &mSpace[0], size_vect);
	// compute x = x + D*alpha
	char trans = 'N';
	dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &aVarsAlpha[0], &I1, &D1, &x_[0], &I1);
	
	double pointValue = mModel->computeLikelihood(x_, mTrace);
	
	
	if(mVerbose > 2)
	{
		std::cout << "Computed the point, coords:\n";
		mModel->printVar(x_, pointValue);
	}
	
	if(!aGrad.empty())
	{
		computeGradientSubspace(pointValue, aVarsAlpha, aGrad);
	}
	
	// Stop optimization if value is greater or equal to threshold
	if(mStopIfBigger && pointValue >= mThreshold) throw nlopt::forced_stop();
	
	return pointValue;
}



// ----------------------------------------------------------------------
void OptSESOP::updateOmega() {0.5 + sqrt(0.25 + square(mOmega));}

// ----------------------------------------------------------------------
void OptSESOP::alocateMemory()
{
	int max_size_alpha = 4;
	size_vect = mN*sizeof(double);
	
	mSubspaceStorage = mN*max_size_alpha;
	mSpace.resize(2*mSubspaceStorage + mN);
	x_.resize(mN);
	
	mGradient = &mSpace[mN];
	mGradient_times = &mSpace[2*mN];
	mGradient_others = &mSpace[3*mN];
	md2 = &mSpace[4*mN];	
	
	alpha.reserve(max_size_alpha);
	mM = 3;
	alpha.resize(mM);
	
#ifdef USE_BOUND_CONSTRAINTS
	mLowerBoundSubspace.reserve(max_size_alpha);
	mUpperBoundSubspace.reserve(max_size_alpha);
	mLowerBoundSubspace.resize(mM);
	mUpperBoundSubspace.resize(mM);
#endif //USE_BOUND_CONSTRAINTS
	
#ifdef USE_TRANSFORMED_CONSTRAINTS
	data_constraints.resize(2*mN);
#endif // USE_TRANSFORMED_CONSTRAINTS
	
	mD = &mSpace[mN+mSubspaceStorage];
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
	
	// second direction of the subspace: gradient_times
	vec = &mD[mN];
	memcpy(vec, mGradient_times, size_vect);
	
	norm = dnrm2_(&mN, vec, &I1);
	inv_norm = 1./norm;
	dscal_(&mN, &inv_norm, vec, &I1);
	
	// third direction
	vec = &mD[2*mN];
	memcpy(vec, mGradient_others, size_vect);
				
	norm = dnrm2_(&mN, vec, &I1);
	inv_norm = 1./norm;
	dscal_(&mN, &inv_norm, vec, &I1);
	
	mM = 3;
	alpha.resize(mM);
	
	// next direction of the subspace: Nemirovski direction (2);
	// only available at second step
	if(mStep > 0)
	{
		vec = &mD[3*mN];
				
		norm = dnrm2_(&N, md2, &I1);
		if( norm < 1e-3 )
		{
			mM = 3;
		}
		else
		{
			mM = 4;
			memcpy(vec, md2, size_vect);
			inv_norm = 1./norm;
			dscal_(&mN, &inv_norm, vec, &I1);
		}
		alpha.resize(mM);		
		
#ifdef USE_BOUND_CONSTRAINTS
		mLowerBoundSubspace.resize(mM);
		mUpperBoundSubspace.resize(mM);
#endif // USE_BOUND_CONSTRAINTS
	}
}

#ifdef USE_TRANSFORMED_CONSTRAINTS
double OptSESOP::operator()(unsigned n, const std::vector<double> &alpha, std::vector<double> &grad, void *data)
{
	data_constraint *data_ = (data_constraint*)(data);
	int i = data_->line;
	int bound_type = data_->bound_type;
	
	// compute the dot product between alpha and the ith row of D
	double Dalpha_i = ddot_(&mM, &alpha[0], &I1, &mD[i], &mN);
	
	
	if(bound_type == 0) // lower bound constraint
	{
		if(!grad.empty())
		{
			dcopy_(&mM, &mD[i], &mN, &grad[0], &I1);	
			dscal_(&mM, &minus_one, &grad[0], &I1);
		}
		return mLowerBound[i] - mSpace[i] - Dalpha_i;
	}
	else if(bound_type == 1) // upper bound constraint
	{
		if(!grad.empty())
		{
			dcopy_(&mM, &mD[i], &mN, &grad[0], &I1);
		}
		return mSpace[i] + Dalpha_i - mUpperBound[i];
	}
	else
	{
		throw FastCodeMLFatal("Exception in computation: Constraint evaluation: Wrong data[1] (corresponding to the bound type, low or up): should be either 0 or 1.\n");
		return 0.;
	}
}
#endif // if USE_TRANSFORMED_CONSTRAINTS

#ifdef USE_BOUND_CONSTRAINTS
// ----------------------------------------------------------------------
void OptSESOP::updateBoundsAndAlpha()
{
	// initialize mLowerBoud to -0.1 and mUpperBound to 0.5
	// and then reduce it if needed
	// TODO: add 'prefered' directions 
	for (int i(0); i<mM; i++)
	{
		mLowerBoundSubspace[i] = -0.001;
		mUpperBoundSubspace[i] = 0.1;
	}
	
	double factor = 0.4;
	int iter = 0;
	int max_iter = 30;
	
	while( !(subspaceLowerBoundIsInSpace() && subspaceUpperBoundIsInSpace()) )
	{
		// reduce the Lower bound
		dscal_(&mN, &factor, &mLowerBoundSubspace[0], &I1);
		// reduce the Upper bound
		dscal_(&mN, &factor, &mUpperBoundSubspace[0], &I1);
		iter ++;
		if (iter > max_iter)
		{
			factor = 0.;
		}
	}
	
	// initialize alpha
	for (int i(0); i<mM; i++)
	{
		alpha[i] = mLowerBoundSubspace[i] + 0.5*(mUpperBoundSubspace[i]-mLowerBoundSubspace[i]);
	}
}


// ----------------------------------------------------------------------
bool OptSESOP::subspaceLowerBoundIsInSpace()
{
	// compute x = x + D*LowerBound
	memcpy(&mSpace[0], &x_[0], size_vect);
	char trans = 'N';
	dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &mLowerBoundSubspace[0], &I1, &D1, &mSpace[0], &I1);
	
	for(int i(0); i<mN; i++)
	{
		if(mSpace[i] < mLowerBound[i] || mSpace[i] > mUpperBound[i])
		{
			return false;
		}	
	}
	return true;
}

bool OptSESOP::subspaceUpperBoundIsInSpace()
{
	// compute x = x + D*UpperBound
	memcpy(&mSpace[0], &x_[0], size_vect);
	char trans = 'N';
	dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &mUpperBoundSubspace[0], &I1, &D1, &mSpace[0], &I1);
	
	for(int i(0); i<mN; i++)
	{
		if(mSpace[i] < mLowerBound[i] || mSpace[i] > mUpperBound[i])
			return false;
	}
	return true;
}
#endif // if USE_BOUND_CONSTRAINTS


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
		
		const double f1 = mModel->computeLikelihood(vars_working_copy, false);
		aGrad[i] = (f1-aPointValue)/delta[i];
		
		vars_working_copy[i] = aVars[i];
	}
}

// ----------------------------------------------------------------------
void OptSESOP::computeGradientSubspace(double aPointValue, const std::vector<double>& aAlpha, std::vector<double>& aGrad)
{
	std::vector<double> vars_working_copy(mM);
	std::vector<double> empty;
	memcpy(&vars_working_copy[0], &aAlpha[0], mM*sizeof(double));

	
	double alpha_norm = dnrm2_(&mM, &aAlpha[0], &I1);
	
	if(alpha_norm > 1e-5)
	{
		for(size_t i(0); i < mM; i++)
		{
			double eh = 1e-6;
			vars_working_copy[i] += eh;
		
			// compute x = x + D*alpha
			memcpy(&x_[0], &mSpace[0], size_vect);
			char trans = 'N';
			dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &vars_working_copy[0], &I1, &D1, &x_[0], &I1);
		
			double f1 = mModel->computeLikelihood(x_, false);
			aGrad[i] = (f1-aPointValue)/eh;
		
			vars_working_copy[i] = aAlpha[i];
		}
	}
	else
	{
		// try with a huge approximation, works only if |alpha| << 1
		for(size_t i(0); i<mM; i++)
		{
			aGrad[i] = ddot_(&mN, mGradient, &I1, &mD[i*mN], &I1);
		}
	}

	if(mVerbose > 2)
	{
		std::cout << "\n Subspace Gradient: \n df/dalpha = [ ";
		for(size_t i(0); i<mM; i++)
			std::cout << aGrad[i] << " ";
		std::cout << " ].\n";
	}
}

