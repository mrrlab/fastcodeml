
#include "Opt2RDSA.h"
#include "MathSupport.h"
#include "lapack.h"


inline double mapKappa_(const double &k) 
{
	double alpha = 0.04;
	double x_lim = pow(alpha, 2./7.);
	
	if(k < x_lim)
		return alpha*sqrt(k);
	else
		return  square(square(k));
}

inline double mapKappaInverse_(const double &k)
{
	double alpha = 0.04;
	double x_lim = pow(alpha, 1./7.);
	
	if(k < x_lim)
		return square(k)/alpha;
	else
		return  pow(k, 0.25);
}



double Opt2RDSA::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	if(mVerbose > 2)
	{
		std::cout << "Size of the problem: N=" << mN << "\n";
	}
	alocateMemory();
	
	double maxl = -1000000;
	int success = Opt2RDSAminimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------

int Opt2RDSA::Opt2RDSAminimizer(double *f, double *x)
{
	//memcpy(&x_[0], x, size_vect);
	//*f = mModel->computeLikelihood(x_, mTrace);
	
	
	// scale x:
	for(size_t i(0); i<mN; i++)
	{
		x[i] = (x[i] - mLowerBound[i]) / (mUpperBound[i] - mLowerBound[i]);
	}
	//x[mNumTimes+3] =  mapKappaInverse_(x[mNumTimes+3] - mLowerBound[mNumTimes+3]) / (mUpperBound[mNumTimes+3] - mLowerBound[mNumTimes+3]);
	
	*f = evaluateFunction(x, mTrace);
	
	// convergence parametters
	bool convergenceReached(false);
	int num_good_iter(0);
	double f_prev(*f);
	
	std::vector<double> y(mN);
	double fy;
	
	
	while(!convergenceReached)
	{
		f_prev = *f;
		
		if(mVerbose > 2)
			std::cout << "Starting step " << mStep << ":\n"; 
			
		updateParameters();
		
		if(mVerbose > 2)
			std::cout << "Parameters updated!\n"; 
		
		computeGradientAndHessian(*f, x);
		
		if(mVerbose > 2)
			std::cout << "Gradient and hessian updated!\n";
		
		
		memcpy(&y[0], x, size_vect);
		fy = *f;
		
		performOneStage(&fy, &y[0]);
		
		if(fy > *f)
		{
			*f = fy;
			memcpy(x, &y[0], size_vect);
		}
		
		memcpy(&y[0], x, size_vect);
		//y[mNumTimes+3] += mGradient[mNumTimes+3]/fabs(mGradient[mNumTimes+3])*sqrt(randFrom0to1())*1e-2;
		if(y[mNumTimes+3] < 0.)
			y[mNumTimes+3] = 0.;
		if(y[mNumTimes+3] >1.)
			y[mNumTimes+3] = 1.;
		
		fy = evaluateFunction(y, false);
		
		if(fy > *f)
		{
			*f = fy;
			memcpy(x, &y[0], size_vect);
		}
		
		
		if(mVerbose > 2)
		{
			std::cout << "State updated! Obtained values (scaled!):\n";
			memcpy(&x_[0], x, size_vect);
			mModel->printVar(x_, *f);
		}
		
		
		mStep ++;
		
		if(fabs(*f-f_prev)<mAbsoluteError)
			num_good_iter++;
		else
			num_good_iter=0;
			
		convergenceReached = (num_good_iter > 20);
	}
	
	// unscale x:
	for(size_t i(0); i<mN; i++)
	{
		x[i] = mLowerBound[i] + x[i] * (mUpperBound[i] - mLowerBound[i]);
	}
	//x[mNumTimes+3] =  mLowerBound[mNumTimes+3] + mapKappaInverse_(x[mNumTimes+3] * (mUpperBound[mNumTimes+3] - mLowerBound[mNumTimes+3]));
}


// ----------------------------------------------------------------------


void Opt2RDSA::alocateMemory()
{
	size_vect = mN*sizeof(double);
	packSize = mN*(mN+1)/2;
	
	// workspace
	// organised as follows:
	// mN elements for workspace
	// mN elements for gradient
	// packSize elements for hessian
	// packSize elements for modified positive definite hessian
	// packSize elements for the workspace to solve the linear system
	mSpace.resize(2*mN + 3*packSize);
	
	// variable to get the permutaions, acts as a workspace
	IPIV.resize(mN);
	
	// vectors used to evaluate the function
	x_.resize(mN);
	x_unscaled.resize(mN);
	
	// make the gradient and Hessian variables point on a region of workspace
	mGradient = &mSpace[mN];
	Hn = mGradient+mN;
	
	// initialize the matrix to zero
	for(size_t i(0); i<packSize; i++)
		Hn[i] = 0.;
}


// ----------------------------------------------------------------------

void Opt2RDSA::updateParameters()
{
	aN = 1. / ( double(mN) + pow(double(mStep+1), 0.602));
	deltaN = 3.8 / pow(double(mStep+1), 0.101);
}


// ----------------------------------------------------------------------

void Opt2RDSA::computeGradientAndHessian(double aPointValue, const double *aVars)
{
	size_t numAv(1);
	double one_over_NumAv = 1./double(numAv);
	
	double alpha_1, alpha_2;
	alpha_1 = double(mStep) / (1.+double(mStep));
	alpha_2 = 1. / (1.+double(mStep));
	
	// set the gradient to 0
	
	for(size_t i(0); i<mN; i++)
		mGradient[i] = 0.;
		
		
	// set the hessian to alpha1*Hn-1 so we can easily add the mean hessian
	dscal_(&packSize, &alpha_1, Hn, &I1);
	
	for(size_t avStep(0); avStep<numAv; avStep++)
	{
		// generate the random varables in the workspace
		double *dN = Hn+packSize;
		double *tmp_grad = &mSpace[0];
	
		bool isPositive;
		for(size_t i(0); i<mN; i++)
		{
			dN[i] = eta * (2.0*randFrom0to1() - 1.0);
		
			// treat boundary problems
			isPositive = (dN[i] >= 0.);
		
			if(aVars[i] - deltaN*fabs(dN[i]) < 0.)
			{
				dN[i] = (0.-aVars[i]) / deltaN;
				if(isPositive)
					dN[i] = -dN[i] - mAbsoluteError;
				else
					dN[i] += mAbsoluteError;
			}
		
			if(aVars[i] + deltaN*fabs(dN[i]) > 1.)
			{
				dN[i] = (1.-aVars[i]) / deltaN;
				if(isPositive)
					dN[i] -= mAbsoluteError;
				else
					dN[i] = -dN[i] + mAbsoluteError;
			}
		}
	
		double ynp, ynm;
		double deltaNm(-deltaN);
	
		// compute yn+
		memcpy(&x_[0], aVars, size_vect); 
		daxpy_(&mN, &deltaN, dN, &I1, &x_[0], &I1);
		ynp = evaluateFunction(x_, false);
	
		// compute yn-
		memcpy(&x_[0], aVars, size_vect); 
		daxpy_(&mN, &deltaNm, dN, &I1, &x_[0], &I1);
		ynm = evaluateFunction(x_, false);
	
		// compute the gradient estimator
		double factor_grad = 1.5 * (ynp-ynm) / (square(eta)*deltaN);
		memcpy(tmp_grad, dN, size_vect);
		dscal_(&mN, &factor_grad, tmp_grad, &I1);
		
		// add contribution of this to the mean gradient
		daxpy_(&mN, &one_over_NumAv, tmp_grad, &I1, mGradient, &I1);
	
	
	
		// Update the Hessian estimator
		double eta_sq_3 = square(eta)/3.;
		double factor_hessian = one_over_NumAv*4.5 * (ynp+ynm-2.*aPointValue) / (square(square(eta)*deltaN));
	
	
		size_t diagId = 0;
		for(size_t i(0); i<mN; i++)
		{
			// diagonal terms
			Hn[diagId] += alpha_2*factor_hessian*2.5 * (square(dN[i]) - eta_sq_3);
		
			for(size_t j(0); j<i; j++)
			{
				Hn[diagId+j+1] += alpha_2*factor_hessian*dN[j]*dN[i];
			}
			diagId += i+2;
		}
	}
}


// ----------------------------------------------------------------------

void Opt2RDSA::performOneStage(double *f, double *x)
{
	// compute the variation delta_x: x_n = x_{n-1} + aN*delta_x
	double *delta_x(&mSpace[0]);
	
	// -- solve the linear system Hn * delta_x = mGradient
	char Uplo = 'U';
	
	int INFO, N_sq;
	N_sq = square(mN);
	
	double *A = Hn + packSize;
//	double *work = Hn + 2*packSize;
	memcpy(A, Hn, mN*size_vect);
	memcpy(delta_x, mGradient, size_vect);

	// change a bit A so it becomes "more" positive definite
	size_t diagId = 0;
	for(size_t i(0); i<mN; i++)
	{
		A[diagId] += deltaN;
		diagId += i+2;
	}
	
	
	// solve the system
	
	if(mVerbose > 2)
		std::cout << "Solving the linear system...\n"; 
	
	dspsv(&Uplo, &mN, &I1, A, &IPIV[0], delta_x, &mN, &INFO);

#ifdef USE_FULL_MATRIX_REPRESENTATION
	dsysv(&Uplo, &mN, &I1, A, &mN
		 ,&IPIV[0]
		 ,delta_x, &mN
		 ,work, &N_sq, &INFO);
#endif

	// check error messages
	if(INFO < 0)
		std::cout << "Argument " << INFO << "Has an illegal value.\n"; 
	else if(INFO > 0)
		std::cout << "The solution could not be computed, singular value occured.\n"; 
	
	if(mVerbose > 2)
		std::cout << "Done!\n"; 
	
	// update the state x
	double minus_aN = -aN;
	daxpy_(&mN, &minus_aN, delta_x, &I1, x, &I1);
	
	// verify the boundary conditions // TODO: better way...
	for(size_t i(0); i<mN; i++)
	{
		if(x[i] < 0.)
			x[i] = mAbsoluteError;
		if(x[i] > 1.)
			x[i] = 1. - mAbsoluteError;
	}
	
	if(mVerbose > 2)
		std::cout << "x updated!\n"; 
	
	if(mVerbose > 2)
		std::cout << "Computing f...\n"; 
	
	*f = evaluateFunction(x, mTrace);
	
	if(mVerbose > 2)
		std::cout << "f updated!\n"; 
}


double Opt2RDSA::evaluateFunction(double *x, bool trace)
{
	memcpy(&x_[0], x, size_vect);
	return evaluateFunction(x_, trace);
}

double Opt2RDSA::evaluateFunction(std::vector<double> &x, bool trace)
{
	for(size_t i(0); i<mN; i++)
	{
		x_unscaled[i] = mLowerBound[i] + x[i]*(mUpperBound[i]-mLowerBound[i]);
	}
	//x_unscaled[mNumTimes+3] =  mLowerBound[mNumTimes+3] + mapKappaInverse_(x[mNumTimes+3] * (mUpperBound[mNumTimes+3] - mLowerBound[mNumTimes+3]));
	
	return mModel->computeLikelihood(x_unscaled, trace);
}

