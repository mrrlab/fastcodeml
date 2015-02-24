#include "CDOSOPtimizer.h"
#include "blas.h"
#include "lapack.h"

inline static double square(double a) {return a*a;}

// ----------------------------------------------------------------
double CDOSOptimizer::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	double lnL(0.0);
	CDOSminimizer(&lnL, &aVars[0]);
	return -lnL; 
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------
int CDOSOptimizer::CDOSminimizer(double *f, double x[])
{
	const double minus_one(-1);
	alocateWorkspace();
	InitSearchDirections();
	
	// stage 1
	// first conjugate direction is the anti-gradient normalized
	std::vector<double> gradient;
	gradient.resize(mN, 0.);
	// TODO compute the gradient
	
#ifdef USE_LAPACK
	double scale = -1./dnrm2_(&mN, &gradient[0], &I1);

	dscal_(&mN, &scale, &gradient[0], &I1);
#else
	double scale = 0.;
	for(int i=0; i < mN; ++i) scale += gradient[i]*gradient[i];
	scale = -1./sqrt(scale);
	for(int i=0; i < mN; ++i) gradient[i] *= scale;
#endif
	
	memcpy(&mU[0], &gradient[0], mN*sizeof(double));
	
	LineSearch(x, mU[0]);
	
	// stage 2
	double mLambdaS (0.62*mLambda);
	std::vector<double> y;
	y.resize(mN);
	
	for(int i(1); i<mN; i++)
	{
		// QR decomposition to find the next orthogonal shift direction mqi
		QRdecomposition(i);
		
		// perform y = x+lamndaS*mqi using blas
		memcpy(&y[0], &x[0], mN*sizeof(double));
		daxpy_(&mN, &mlambdaS, &y[0], &I1, &mqi[0]);
		
		for(int j(0); j<i; j++)
		{
			LineSearch(&y[0], &mU[j*mN]);
		}
		// update ui
		// TODO
	}
	
	// stage 3 (for non quadratic functions only, which is here the case)
	
	// compute the step length
	//daxpy_(&mN, &minus_one,	double *x,	const int *incx, double *y,	const int *incy);
	
	bool stop_condition_reached(false);
	while(!stop_condition_reached)
	{
		//TODO
	}
	
	
	return 0;
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------
void QRdecomposition(int width)
{
	// copy the 'width' first vectors of U into Q
	memcpy(&mQ[0], &mU[0], mN*width*sizeof(double));
	
	// size of the space to use
	const int lwork = mlwork;
	int info;
	
	// use the lapack functions to first decompose mQ
	dgeqrf(&mN, &width, &mQ[0], &width, &mSpace[lwork], &mSpace[0], &lwork, &info);
	//TODO: check errors
	// and then compute the first 'width' vectors of Q
	dorgqr(&mN, &width, &width, &mQ[0], &width, &mSpace[lwork], &mSpace[0], &lwork, &info);
	//TODO: check errors
	
	// copy the interesting vector in mqi
	memcpy(&mqi[0], &mQ[mN*width], mN*sizeof(double));
}


// ----------------------------------------------------------------
void CDOSOptimizer::initSearchDirections()
{
	mU.resize(mN*mN, 0.0);
	for(size_t i(0); i<mN; i++)
	{
		mU[i*(mN+1)] = 1.0;
	}
	mQ.resize(mN*mN);
	mqi.resize(mN);
}

// ----------------------------------------------------------------
void CDOSOptimizer::alocateWorkspace()
{
	mlwork = mN*mN;
	mSpace.resize(mlwork+mN);
}
