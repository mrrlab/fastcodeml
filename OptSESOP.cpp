
#include "OptSESOP.h"
#include "MathSupport.h"
#include "lapack.h"
#include "nlopt.hpp"

double OptSESOP::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	if(mVerbose > 2)
		std::cout << "Size of the problem: N=" << mN << "\n";

	alocateMemory();
	
	double maxl = -1000000;
	int success = SESOPminimizer(&maxl, &aVars[0]);
	return -maxl;
}



// ----------------------------------------------------------------------
int OptSESOP::SESOPminimizer(double *f, double *x)
{	
	// Store the initial state
	
	std::vector<double> x_0;
	x_0.resize(mN);
	memcpy(&x_0[0], x, size_vect);
	memcpy(mXPrev,  x, size_vect);	
	
	*f = mModel->computeLikelihood(x_0, mTrace);
	double f_prev = *f - 2.;
	
	// nLopt optimizer
	std::auto_ptr<nlopt::opt> opt;

	// initialize the matrix D
	computeGradient(*f, x, mGradient);
	
	memcpy(mGradPrev, mGradient, size_vect);
	memcpy(mGradient_times, mGradient, mNumTimes*sizeof(double));
	
	for(int i(mNumTimes); i<mN; i++)
		mGradient_times[i] = 0.;
			
	
	// mGradient_others
	for(size_t i(0); i<mN; i++)
		mGradient_others[i] = 0.;
	mGradient_others[mNumTimes+1] = mGradient[mNumTimes+1];	// v1
	mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
	//mGradient_others[mNumTimes+3] = mGradient[mNumTimes+3]; // kappa
		
	// Initialize the second Nemirovski direction
	memcpy(md2, mGradient, size_vect);
	
	// Initialize the Hessian matrix to identity
	for(size_t i(0); i<mN; i++) {mHdiag[i] = 1.;}
	
	
	// loop until convergene reached
	bool convergence_reached( false );
	while(!convergence_reached)
	{
		double f_diff = fabs(*f - f_prev);
		f_prev = *f;		
		
		double stepTolerance, orderF_diff;
		orderF_diff = log(f_diff)/log(10.);
		
		stepTolerance = max2(mRelativeError, exp(-double(mStep)));
		//stepTolerance = max2(mRelativeError, pow(0.1, 5.-orderF_diff) );
		
		std::cout << "Step " << mStep << ", tolerance " << stepTolerance << ", order of f_diff:" << orderF_diff << ".\n";
		
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
			// TODO: improve with correlations?
			for(size_t i(mNumTimes); i<mN; i++)
			{
				mGradient_others[i] = 0.;
			}
			double rand_tmp = 1. - (1. - exp(-f_diff/mRelativeError)) * randFrom0to1();
			if(rand_tmp < 0.6)
			{
				mGradient_others[mNumTimes+1] = mGradient[mNumTimes+1];	// v1
				mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
			}
			else if(rand_tmp < 0.85)
			{	
				mGradient_others[mNumTimes+3] = mGradient[mNumTimes+3];	// kappa
				mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
			}
			else
			{
				mGradient_others[mNumTimes+3] = 1.; // kappa
			}
			
			// md2
			updateOmega();
			daxpy_(&mN, &mOmega, mGradient, &I1, md2, &I1);	
			
			
			// -- Hessian approximation
			
			//mS
			memcpy(mS, x, size_vect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mS, &I1);
			
			//mY
			memcpy(mY, mGradient, size_vect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mY, &I1);
			
			// Hessian approximation update (diagonal only)
			double sum_s_sq = 0.;
			double sTDs		= 0.;
			#pragma omp parallel for reduction(+:sum_s_sq) reduction(+:sTDs)
			for(size_t i(0); i<mN; ++i)
			{
				// local variable
				double s_sq = square(mS[i]);
				
				sTDs 	 += mHdiag[i]*s_sq;
				sum_s_sq += square(s_sq);
				
				mSpace[i] = s_sq;
			}
			double SY = ddot_(&mN, mS, &I1, mY, &I1);
			double factor = (SY - sTDs) / sum_s_sq;
			
			daxpy_(&mN, &factor, &mSpace[0], &I1, mHdiag, &I1);
		}
		updateDMatrix();
		
		
		if(mVerbose > 3	)
		{
			std::cout << "\n\n --------------- Matrix D: ------------- \n\n";
			for(int j(0); j<mM; ++j)
			{
				memcpy(&x_[0], mD+j*mN, size_vect);
				std::cout << "\nColumn " << j << " of the D matrix:\n";
				mModel->printVar(x_, 0);
			}
		}
		
		
		// save position and gradient
		memcpy(mGradPrev, mGradient , size_vect);
		memcpy(mXPrev	, x			, size_vect);
		
		
		if(mVerbose > 2)
			std::cout << "Matrix D updated, optimizing...\n";
				
		// optimize in the subspace, i.e. find best alpha s.t. f(x+D*alpha) is minimized
		// we store the current state in the workspace
		memcpy(&mSpace[0], x, size_vect); 
		
		opt.reset(new nlopt::opt(nlopt::LD_SLSQP, mM));
		//opt->set_vector_storage(mM);
		
		//opt.reset(new nlopt::opt(nlopt::LN_COBYLA, mM));
		
		// add the constraints
		for(int constraint_id(0); constraint_id<2*mN; constraint_id++)
			opt->add_inequality_constraint(OptSESOP::myconstraintWrapper, (void*) &data_constraints[constraint_id], 1e-6);
		
		if(mVerbose > 2)
			std::cout << "Inequality constraints added\n";
			
		
		// initialize alpha to 0 so we start at point x
		for(int i(0); i<mM; i++)
			alpha[i] = 0.;

		opt->set_ftol_rel(stepTolerance);
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
			//throw FastCodeMLFatal("Exception in computation: Halted because roundoff errors limited progress, equivalent to NLOPT_ROUNDOFF_LIMITED.");
			std::cout << "Esception NLOPT_ROUNDOFF_LIMITED caught. continue in another subspace...\n";
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
		
		opt->remove_inequality_constraints();
					
		
		// update current state
		char trans = 'N';
		dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &alpha[0], &I1, &D1, x, &I1);
		
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			if(x[i] < mLowerBound[i])
				x[i] = mLowerBound[i];
			if(x[i] > mUpperBound[i])
				x[i] = mUpperBound[i]-mRelativeError;
		}
		
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
		//size_t num_iter_move = 1+size_t(4. / (1 + square(2.*float(mStep)/float(mN))));
		size_t num_iter_move = 1+size_t(5. / (1 + square(2.*float(mStep)/float(mN))));
		
		for(; iter_move < num_iter_move; ++iter_move)
		{
			memcpy(&x_[0], x, size_vect);
			for(size_t i(0); i<mNumTimes; ++i)
			{				
				x_[i] += 1e-1 * square(x[i])*square(randFrom0to1()) * mGradient[i];
							
				
				if(x_[i] < mLowerBound[i])
					x_[i] = mLowerBound[i];
			}
			double f_test = mModel->computeLikelihood(x_, mTrace);
			
			if(mVerbose > 2) std::cout << "Step " << mStep << ", previous f: " << *f << ", new one: " << f_test << std::endl;
			
			if (f_test > *f)
			{
				*f = f_test;
				memcpy(x, &x_[0], size_vect);
			}
		}		
		
		// check convergence
		convergence_reached = f_diff < mRelativeError;

		
		if(mStep > mN) // TODO: have to find optimal parameter here
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
void OptSESOP::updateOmega() {mOmega = 0.5 + sqrt(0.25 + square(mOmega));}

// ----------------------------------------------------------------------
void OptSESOP::alocateMemory()
{
	int max_size_alpha = 4;
	size_vect = mN*sizeof(double);
	
	mSubspaceStorage = mN*max_size_alpha;
	mSpace.resize(2*mSubspaceStorage + 6*mN);
	// memory space:
	// mN of workspace
	// mN for gradient
	// mN for gradient_times
	// mN for gradient others
	// mN for md2
	// mSubspaceStorage for matrix D
	// mN for Hessian diagonal
	// mN for mS
	// mN for mXPrev
	// mN for mY
	// mN for mGradPrev
	
	x_.resize(mN);
	
	mGradient = &mSpace[mN];
	mGradient_times 	= mGradient			+ mN;
	mGradient_others 	= mGradient_times	+ mN;
	md2 				= mGradient_others	+ mN;	
	
	mD = &mSpace[mN+mSubspaceStorage];
	
	mHdiag 		= mD		+ mSubspaceStorage;
	mS			= mHdiag	+ mN;
	mXPrev		= mS		+ mN;
	mY			= mXPrev	+ mN;
	mGradPrev	= mY		+ mN;
	
	mM = 3;
	
	alpha.reserve(max_size_alpha);
	alpha.resize(mM);
	
	data_constraints.resize(2*mN);
	
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i)
	{
		data_constraints[2*i  ].sesop = this;
		data_constraints[2*i  ].line = i;
		data_constraints[2*i  ].bound_type = 0;
		data_constraints[2*i+1].sesop = this;
		data_constraints[2*i+1].line = i;
		data_constraints[2*i+1].bound_type = 1;
	}
	
	
}

	
// ----------------------------------------------------------------------
void OptSESOP::updateDMatrix()
{
	// local variables
	double *vec;

	// first direction of the subspace: modified gradient
	memcpy(mD, mGradient, size_vect);
	
#if 1
	#pragma omp parallel for	
	for(size_t i(0); i<mN; ++i) 
	{
		if(fabs(mHdiag[i]) > 1e-5)
			mD[i] /= mHdiag[i];
	}
#endif
	normalizeVector(&mN, mD);
	
	// second direction of the subspace: gradient_times
	vec = &mD[mN];
	memcpy(vec, mGradient_times, size_vect);
	normalizeVector(&mN, vec);	
	
	// third direction
	vec = &mD[2*mN];
	memcpy(vec, mGradient_others, size_vect);
	normalizeVector(&mN, vec);
	
	mM = 3;
	
	// next direction of the subspace: Nemirovski direction (2);
	// only available at second step
	if(mStep > 0)
	{
		vec = &mD[3*mN];
		memcpy(vec, md2, size_vect);
		if( normalizeVector(&mN, vec) > 1e-5 )
			mM = 4;
	}
	alpha.resize(mM);
}


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

// ----------------------------------------------------------------------
void OptSESOP::computeGradient(double aPointValue, const double *aVars, double* aGrad)
{
	memcpy(&x_[0], aVars, size_vect);
	double eh, delta;
	
	for(size_t i(0); i < mN; ++i)
	{
		eh = 1e-7 * max2(abs(x_[i]), 1.);
			
		// If it is going over the upper limit reverse the delta
		x_[i] += eh;
		if(x_[i] >= mUpperBound[i])
		{
			x_[i] -= 2*eh;
			delta = -eh;
		}
		else
		{
			delta = eh;
		}
		
		const double f1 = mModel->computeLikelihood(x_, false);
		aGrad[i] = (f1-aPointValue)/delta;
		
		x_[i] = aVars[i];
	}
}

// ----------------------------------------------------------------------
void OptSESOP::computeGradientSubspace(double aPointValue, const std::vector<double>& aAlpha, std::vector<double>& aGrad)
{
	std::vector<double> vars_working_copy(mM);
	memcpy(&vars_working_copy[0], &aAlpha[0], mM*sizeof(double));

	
	double alpha_norm = dnrm2_(&mM, &aAlpha[0], &I1);
	
	if(alpha_norm > 1e-5)
	{
		for(size_t i(0); i < mM; i++)
		{
			double eh = 1e-7 * max2(aAlpha[i], 1.);
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
			aGrad[i] = ddot_(&mN, mGradient, &I1, &mD[i*mN], &I1);
	}

	if(mVerbose > 2)
	{
		std::cout << "\n Subspace Gradient: \n df/dalpha = [ ";
		for(size_t i(0); i<mM; i++)
			std::cout << aGrad[i] << " ";
		std::cout << "].\n";
	}
}

