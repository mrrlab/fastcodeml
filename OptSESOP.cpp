
#include "OptSESOP.h"
#include "MathSupport.h"
#include "lapack.h"
#include "nlopt.hpp"

// ----------------------------------------------------------------------
//	Utility functions
// ----------------------------------------------------------------------

#ifdef SESOP_CORR_SELECTION
static size_t binarySearch(double aValue, double* aVect, size_t aLeft, size_t aRight)
{
	if(aLeft >= aRight-1)
		return aLeft;

	size_t mid = (aLeft + aRight)>>1;
	
	if(aValue > aVect[mid]) 
		return binarySearch(aValue, aVect, mid, aRight);
	else
		return binarySearch(aValue, aVect, aLeft, mid);
}
#endif

// ----------------------------------------------------------------------
//	Class members definition
// ----------------------------------------------------------------------
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

	*f = mModel->computeLikelihood(x_0, mTrace);
	double f_prev = *f - 2.;
	
	// nLopt optimizer
	std::auto_ptr<nlopt::opt> opt;

	// initialize the matrix D
	computeGradient(*f, x, mGradient);
	
	// first take all the branch lengths 
	memcpy(mGradient_times, mGradient, mNumTimes*sizeof(double));
	for(size_t i(mNumTimes); i<mN; ++i) mGradient_times[i] = 0.;
	
#ifdef SESOP_HESSIAN_APPROX
	memcpy(mXPrev,  x, size_vect);	
	memcpy(mGradPrev, mGradient, size_vect);
	
	// Initialize the Hessian matrix to identity
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i) mHdiag[i] = 1.;
#endif	
	
	// mGradient_others
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i) mGradient_others[i] = 0.;
		
	mGradient_others[mNumTimes+1] = mGradient[mNumTimes+1];	// v1
	mGradient_others[mNumTimes+2] = mGradient[mNumTimes+2]; // w0
	//mGradient_others[mNumTimes+3] = mGradient[mNumTimes+3]; // kappa
		
	// Initialize the second Nemirovski direction
	memcpy(md2, mGradient, size_vect);
	
	
	// loop until convergene reached
	// local variables
	double f_diff, stepTolerance;
	size_t numIterBeforeConverged( 0 ),
		   maxNumIterBeforeConverged( 1 );
	
	bool convergence_reached( false );
	while(!convergence_reached)
	{
		f_diff = fabs(*f - f_prev);
		f_prev = *f;
		
		// adaptative tolerance so we don't spend too much time at each iteration
		stepTolerance = max2(mAbsoluteError, exp(-double(mStep)));
		
		if(mVerbose > 2)
			std::cout << "Starting step " << mStep << ":\n";
			
		// compute/update search directions
		if(mStep > 0)
		{	
			// compute and store all the required vectors to construct the subspace
			
			// mGradient
			computeGradient(*f, x, mGradient);
			
			// mGradient_times			
#ifdef SESOP_CORR_SELECTION
			// Correlation matrix
			updateCorr();
			// select variables
			selectCorrVariables();			
#else
			memcpy(mGradient_times, mGradient, mNumTimes*sizeof(double));	
#endif
			// gradient other
			for(size_t i(mNumTimes); i<mN; i++)
			{
				mGradient_others[i] = 0.;
			}
			double rand_tmp = randFrom0to1();
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
			

#ifdef SESOP_HESSIAN_APPROX			
			
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
			
#endif
		}
		// copy and normalize all the search directions in matrix D
		updateDMatrix();
		
		
		if(mVerbose > 2	)
		{
			std::cout << "\n\n --------------- Matrix D: ------------- \n\n";
			for(int j(0); j<mM; ++j)
			{
				memcpy(&x_[0], mD+j*mN, size_vect);
				std::cout << "\nColumn " << j << " of the D matrix:\n";
				mModel->printVar(x_, 0);
			}
		}
		
#ifdef SESOP_HESSIAN_APPROX
		// save position and gradient
		memcpy(mGradPrev, mGradient , size_vect);
		memcpy(mXPrev	, x			, size_vect);
#endif
		
		if(mVerbose > 2)
			std::cout << "Matrix D updated.\n";
				
		// optimize in the subspace, i.e. find best alpha s.t. f(x+D*alpha) is minimized
		// we store the current state in the workspace
		memcpy(&mSpace[0], x, size_vect); 
		
#ifdef SESOP_USE_COBYLA
		opt.reset(new nlopt::opt(nlopt::LN_COBYLA, mM));
#else
		opt.reset(new nlopt::opt(nlopt::LD_SLSQP, mM));
		//opt->set_vector_storage(mM);
#endif	
		// add the constraints
		for(size_t constraint_id(0); constraint_id<2*mN; ++constraint_id)
			opt->add_inequality_constraint(OptSESOP::myconstraintWrapper, &data_constraints[constraint_id], 1e-6);
		
		if(mVerbose > 2)
			std::cout << "Inequality constraints added\n";
		
		
		// initialize alpha to 0 so the optimizer starts at current point x
		for(size_t i(0); i<mM; ++i)	alpha[i] = 0.;
		
		try
		{
			opt->set_ftol_rel(stepTolerance);
			//opt->set_ftol_abs(stepTolerance);
			opt->set_max_objective(OptSESOP::subspaceEvaluatorWrapper, this);
			
#ifdef SESOP_USE_COBYLA
			// compute a step smaller and smaller as mStep grows so we optimize in a smaller subregion
			mInitStepCobyla.resize(mM);
			if(mStep == 0)
			{
				opt->get_initial_step(alpha, mInitStepCobyla);
				mInitStepLenghth0 = dnrm2_(&mM, &mInitStepCobyla[0], &I1);
			}
			else
			{
				double step_length, gamma, minStepPossible;
				gamma = 1e-2;
				double *p = mGradient; // initial step direction
				step_length = mInitStepLenghth0 * gamma/(gamma+double(mStep));
				
				// compute the minumum step length we can do
				minStepPossible = step_length + 1.;
				//#pragma omp paralel for reduction(min:minStepPossible)
				for(size_t i(0); i<mN; ++i)
				{
					double minTmp;
					if(fabs(p[i]) > 1e-5)
					{
						if(p[i] > 0.)
						{
							minTmp = (mUpperBound[i]-x[i])/p[i];
							if(minStepPossible > minTmp)
								minStepPossible = minTmp;
							std::cout << "debug: minTmp " << minTmp << " xi " << x[i] << " ui "<< mUpperBound[i] << " p1 " << p[i] << std::endl;
						}
						else
						{
							minTmp = (mLowerBound[i]-x[i])/p[i];
							if(minStepPossible > minTmp)
								minStepPossible = minTmp;
							std::cout << "debug: minTmp " << minTmp << " xi " << x[i] << " li "<< mLowerBound[i] << " p1 " << p[i] << std::endl;
						}
					}
				}
				if( step_length > minStepPossible )
					step_length = 0.5 * minStepPossible;
			
				if(step_length < 0.)
				{
					for(size_t i(0); i<mM; ++i) mInitStepCobyla[i] = 0.;
					mInitStepCobyla[0] = step_length;
					std::cout << "Size of the first step: " << step_length << ", array of size " << mInitStepCobyla.size() << "\n";
					opt->set_initial_step(mInitStepCobyla);
				}
			}
#endif
			if(mVerbose > 2)
				std::cout << "Optimizer set, optimizing...\n";
			
			
			nlopt::result result = opt->optimize(alpha, *f);
			opt->remove_inequality_constraints();
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
		
		// update current state
		char trans = 'N';
		dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &alpha[0], &I1, &D1, x, &I1);
		
		#pragma omp parallel for
		for(size_t i(0); i<mN; ++i)
		{
			if(x[i] < mLowerBound[i])
				x[i] = mLowerBound[i];
			if(x[i] > mUpperBound[i])
				x[i] = mUpperBound[i]-mAbsoluteError;
		}
		
		if(mVerbose > 2)
		{
			std::cout << "Obtained alpha:\n";
			for(int i(0); i<mM; i++)
				std::cout << alpha[i] << " ";			
		}
		
#if 1
		// test to make the times move...
		size_t iter_move = 0;
		//size_t num_iter_move = 1+size_t(4. / (1 + square(2.*float(mStep)/float(mN))));
		size_t num_iter_move = 1+size_t(6. / (1 + square(2.*float(mStep)/float(mN))));
		
		for(; iter_move < num_iter_move; ++iter_move)
		{
			memcpy(&x_[0], x, size_vect);
			for(size_t i(0); i<mNumTimes; ++i)
			{				
				x_[i] += 1e-1 * square(x[i])*square(randFrom0to1()) * mGradient[i];
							
				if(x_[i] < mLowerBound[i])
					x_[i] = mLowerBound[i];
				
				if(x_[i] > mUpperBound[i])
					x_[i] = mUpperBound[i]-1e-3;
			}
			double f_test = mModel->computeLikelihood(x_, mTrace);
			
			if(mVerbose > 2) std::cout << "Step " << mStep << ", previous f: " << *f << ", new one: " << f_test << std::endl;
			
			if (f_test > *f)
			{
				*f = f_test;
				memcpy(x, &x_[0], size_vect);
			}
		}		
		
#endif
		if(mVerbose > 2)
		{
			std::cout << "\n\nOptimized at this step, obtained likelihood: ";
			memcpy(&x_[0], x, size_vect);
			mModel->printVar(x_, *f);
		}

		// check convergence
		
		if(f_diff < mAbsoluteError)
			++numIterBeforeConverged;
		else
			numIterBeforeConverged = 0;
		
		convergence_reached = numIterBeforeConverged >= maxNumIterBeforeConverged;

		
		if(mStep > mN)
		{
			if(mVerbose > 2)
				std::cout << "SESOP aborted, number of iterations exceeded.\n";
			return 1;
		}
		
		mStep ++;
	}
		
	return 0;
}

// ----------------------------------------------------------------------
double OptSESOP::eValuateFunctionSubspace(const std::vector<double>& aVarsAlpha, std::vector<double>& aGrad)
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
void OptSESOP::updateOmega(void) {mOmega = 0.5 + sqrt(0.25 + square(mOmega));}

// ----------------------------------------------------------------------
void OptSESOP::alocateMemory(void)
{
	int max_size_alpha = 4;
	mSubspaceStorage = mN*max_size_alpha;
	size_t size_mSpace(mN + 2*mSubspaceStorage);
	size_vect = mN*sizeof(double);
	
#ifdef SESOP_HESSIAN_APPROX
	size_mSpace += 5*mN;
#endif	

	mSpace.resize(size_mSpace);
	// --- memory space:
	// mN of workspace
	// mN for gradient
	// mN for gradient_times
	// mN for gradient others
	// mN for md2
	// mSubspaceStorage for matrix D
	// --- following only ifdef SESOP_HESSIAN_APPROX
	// mN for Hessian diagonal
	// mN for mS
	// mN for mXPrev
	// mN for mY
	// mN for mGradPrev
	// ---
	
	x_.resize(mN);
	
	mGradient = &mSpace[mN];
	mGradient_times 	= mGradient			+ mN;
	mGradient_others 	= mGradient_times	+ mN;
	md2 				= mGradient_others	+ mN;	
	
	mD = &mSpace[mN+mSubspaceStorage];

#ifdef SESOP_HESSIAN_APPROX
	mHdiag 		= mD		+ mSubspaceStorage;
	mS			= mHdiag	+ mN;
	mXPrev		= mS		+ mN;
	mY			= mXPrev	+ mN;
	mGradPrev	= mY		+ mN;
#endif


#ifdef SESOP_CORR_SELECTION
	// --- memory space
	// size_corr for the correlation matrix (don't store the diagonal and the lower triangle part)
	// size_corr for the prefix sum
	// size_corr for the prefix sum workspace
	// ---
	size_corr = mN*(mN-1)>>1;
 #ifdef SESOP_CORR_MEM_OPTIM
	mCorrContainer.resize(size_corr);
	mCorr				= &mCorrContainer[0];
 #else
	mCorrContainer.resize(3*size_corr);
	mCorr				= &mCorrContainer[0];
	mCorrPrefSumm 		= mCorr			+ size_corr;
	mPrefSumWorkspace	= mCorrPrefSumm	+ size_corr;
 #endif
	#pragma omp parallel for
	for(size_t k(0); k<size_corr; ++k) mCorr[k] = 0.;
#endif
	mM = 3;
	
	alpha.reserve(max_size_alpha);
	alpha.resize(mM);
	
#ifdef SESOP_USE_COBYLA
	mInitStepCobyla.reserve(max_size_alpha);
	mInitStepCobyla.resize(mM);
#endif
	
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
void OptSESOP::updateDMatrix(void)
{
	// local variables
	double *vec;

	// first direction of the subspace: (modified) gradient
	memcpy(mD, mGradient, size_vect);
	
#ifdef SESOP_HESSIAN_APPROX
	#pragma omp parallel for	
	for(size_t i(0); i<mN; ++i) 
	{
		if(fabs(mHdiag[i]) > 1e-3)
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
		// only keep this direction if it is not too small
		if( normalizeVector(&mN, vec) > 1e-5 )
			mM = 4;
	}
	alpha.resize(mM);
}

// ----------------------------------------------------------------------
#ifdef SESOP_CORR_SELECTION
void OptSESOP::updateCorr(void)
{
	// update the correlation matrix
	double an = 1./(1.+double(mStep));
	
	#pragma omp parallel for
	for(size_t i(1); i<mN; ++i)
	{
		for(size_t j(0); j<i; ++j)
		{
			// index of the correlation between ith and jth variables
			size_t k = (i*(i-1))>>1 + j;
			mCorr[k] = (1.-an)*mCorr[k] + an*fabs(mGradient[i]*mGradient[j]);
		}
	}
		
 #ifdef SESOP_CORR_MEM_OPTIM

 #else	
	// compute the prefix sum with openMP
	memcpy(mCorrPrefSumm, mCorr, size_corr*sizeof(double));
	
	for(size_t step(0); step<log2(size_corr); ++step)
	{
		size_t i;
		#pragma omp parallel private(i)
		{
			#pragma omp for
			for(i = 1<<step; i<size_corr; ++i)
				mPrefSumWorkspace[i] = mCorrPrefSumm[i] + mCorrPrefSumm[i-(1<<step)];				
			#pragma omp barrier 
			#pragma omp for
			for(i = 1<<step; i<size_corr; ++i)
				mCorrPrefSumm[i] = mPrefSumWorkspace[i];
			#pragma omp barrier 
		}
	}
 #endif
}


// ----------------------------------------------------------------------
void OptSESOP::selectCorrVariables(void)
{
	double *columnCorrVar = mGradient_times;
	// randomly choose the search directions for the matrix parameters using the correlation:
	
	#pragma omp parallel for
	for(size_t i(0); i<mN; ++i) columnCorrVar[i] = 0.;
		
	size_t numVarSelected(sqrt(mN));
	
	#pragma omp parallel for
	for(size_t i(0); i<numVarSelected; ++i)
	{
		size_t id1, id2, k;
 #ifdef SESOP_CORR_MEM_OPTIM
		size_t rand_stencil_selector = 1+randFrom0to1()*5;
		size_t rand_selector = rand_stencil_selector + randFrom0to1() * (size_corr-1-rand_stencil_selector);
		k = rand_selector;
		for(size_t walker(0); walker<rand_stencil_selector; ++walker)
		{
			if(mCorr[rand_selector-walker] > mCorr[k])
				k = rand_selector-walker;
			if(mCorr[rand_selector+walker] > mCorr[k])
				k = rand_selector+walker;
		}
 #else
		double rand_selector = randFrom0to1() * mCorrPrefSumm[size_corr-1];
		
		// Find the index k in the correlation matrix using prefix sum.
		// this is the good way to select k but it is costly in terms 
		// of computational time and memory (requires to store the prefix-
		// sum array + a working space)
		k = binarySearch(rand_selector, mCorrPrefSumm, 0, size_corr);
 #endif		
		// Find the indices id1 and id2 of the variables
		// 
		// this is done as we have k = i(i-1)/2 + j with j<i
		// so i^2-i <= 2*k < i^2+i
		
		id1 = size_t(round(sqrt(2*k)));
		if(id1*id1+id1 <= 2*k) ++id1;
		id2 = k - ((id1-1)*id1>>1);

		std::cout << "Selected variables " << id1 << " and " << id2 << "\n";
		
		// add the variables gradient to the search direction
		columnCorrVar[id1] = mGradient[id1];
		columnCorrVar[id2] = mGradient[id2];
	}
}
#endif

// ----------------------------------------------------------------------
double OptSESOP::eValuateConstraintsSubspace(const std::vector<double> &alpha, std::vector<double> &grad, void *data)
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
#ifdef SESOP_STOCHASTIC_GRADIENT
		if( i>=mNumTimes || randFrom0to1() <= (mStep == 0 ? 1. : 0.4 + 0.6*randFrom0to1())) //TODO
		{
#endif
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
#ifdef SESOP_STOCHASTIC_GRADIENT
		}
		else
		{
			if( i==0 )
			{
				aGrad[i] = 0.5*randFrom0to1()(mGradient[i+1] + mGradient[i+2]);
			}
			else
			{
				aGrad[i] = 0.5*(mGradient[i-1] + mGradient[i+1]);
			}
		}
#endif
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
		volatile double da;
		double eh, fp;		
		char trans = 'N';
		
		for(size_t i(0); i < mM; ++i)
		{
			eh = 1e-6 * (1. + randFrom0to1());// * max2(fabs(aAlpha[i]), 1.);
			
			// take only points between x and x+alpha so we have more chance to be in the 
			if(aAlpha[i] > 0.) eh = -eh; 
			
			// compute in the forward finite differences
			vars_working_copy[i] += eh;
			da = vars_working_copy[i]-aAlpha[i];
			
			// compute x + D*(alpha+da)
			memcpy(&x_[0], &mSpace[0], size_vect);
			dgemv_(&trans, &mN, &mM, &D1, mD, &mN, &vars_working_copy[0], &I1, &D1, &x_[0], &I1);
		
			fp = mModel->computeLikelihood(x_, false);
			aGrad[i] = (fp-aPointValue)/da;
		}
	}
	else
	{
		// try with a huge approximation, works only if |alpha| << 1
		for(size_t i(0); i<mM; ++i)
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

