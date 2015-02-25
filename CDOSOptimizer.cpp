
#include "CDOSOptimizer.h"
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
int CDOSOptimizer::CDOSminimizer(double *f, double *x)
{
	// local variables
	double fx, fy, lambda_i;
	std::vector<double> x_prev(mN);

	alocateWorkspace();
	initSearchDirections();
	
	// ------------------------------ stage 1 ------------------------------
	
	// first conjugate direction is the anti-gradient normalized
	std::vector<double> gradient;
	gradient.resize(mN);
	
	fy = -mModel->computeLikelihood(x, mN, mTraceFun);
	fx = fy;
	
	// compute the gradient at point x
	computeGradient(fx, x, &gradient[0]);
	
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
	
	LineSearch(x, &mU[0], mLambda, &fy);
	
	// ------------------------------ stage 2 ------------------------------
	
	double lambdaS (0.62*mLambda);
	std::vector<double> y;
	y.resize(mN);
	
	for(int i(1); i<mN; i++)
	{
		// QR decomposition to find the next orthogonal shift direction mqi
		QRdecomposition(i);
		
		// perform y = x+lamndaS*mqi using blas
		memcpy(&y[0], &x[0], mN*sizeof(double));
		daxpy_(&mN, &lambdaS, &y[0], &I1, &mqi[0], &I1);
		
		// look for a new point y (better than x)
		for(int j(0); j<i-1; j++)
		{	
			// line search in uj direction
			// update y if better point
			LineSearch(&y[0], &mU[j*mN], lambdaS, &fy);
		}
		
		// --- update ui
			double searchDirection = (fx<fy) ? -1. : 1.;
			// store y-x in ui
			memcpy(&mU[mN*i], &y[0], mN*sizeof(double));
			daxpy_(&mN, &minus_one, &x[0], &I1, &mU[mN*i], &I1);
			// scale it
			double scale = searchDirection/dnrm2_(&mN, &mU[mN*i], &I1);
			dscal_(&mN, &scale, &mU[mN*i], &I1);
		// --
		
		// update x
		memcpy(&x_prev[0], &x[0], mN*sizeof(double));
		if(fy < fx)
		{
			memcpy(&x[0], &y[0], mN*sizeof(double));
			fx = fy;
		}
	}
	
	// ------------------------------ stage 3 ------------------------------
	// (for non quadratic functions only, which is here the case)
	
	
	// compute the step value lambda_i 
	daxpy_(&mN, &minus_one,	&x[0], &I1, &x_prev[0],	&I1);
	
	lambda_i = 0.3 * dnrm2_(&mN, &x_prev[0], &I1) + 0.091 * mLambda;
	lambdaS = 0.62 * lambda_i;
	
	int num_iter(0);
	bool stop_condition_reached(false);
	while(!stop_condition_reached)
	{
		// update the search directions
		// reverse the column order (is it needed??)
		for(int column(0); column<mN/2; column++)
		{
			swap_content(&mU[column*mN], &mU[(mN-column-1)*mN], mN);
		}
		// compute the new direction
		QRdecomposition(mN-1);
		
		// change the sign if necessary
		double d1 = distance(&mqi[0], &mU[0], &mSpace[0], mN);
		double d2 = distance(&mqi[0], &mU[0], &mSpace[0], mN, -1.);
		if (d1 > d2) {dscal_(&mN, &minus_one, &mqi[0], &I1);}
		
		// perform y = x+lamndaS*mqi using blas
		memcpy(&y[0], &x[0], mN*sizeof(double));
		daxpy_(&mN, &lambdaS, &y[0], &I1, &mqi[0], &I1);
		
		// -- shift the conjugate directions
		memcpy(&mSpace[0], &mU[0], mN*sizeof(double));
		for(int column(0); column < mN-1; column++)
			memcpy(&mU[column*mN], &mU[(column+1)*mN], mN*sizeof(double));
		memcpy(&mU[(mN-1)*mN], &mSpace[0], mN*sizeof(double));
		// --
	
		for(int j(0); j<mN-1; j++)
		{
			// line search in uj direction
			// update y if better point
			LineSearch(&y[0], &mU[j*mN], 3.*lambdaS, &fy);
		}
		
		// --- update Un
			double searchDirection = (fx<fy) ? -1. : 1.;
			// store y-x in ui
			memcpy(&mU[mN*(mN-1)], &y[0], mN*sizeof(double));
			daxpy_(&mN, &minus_one, &x[0], &I1, &mU[mN*(mN-1)], &I1);
			// scale it
			double scale = searchDirection/dnrm2_(&mN, &mU[mN*(mN-1)], &I1);
			dscal_(&mN, &scale, &mU[mN*(mN-1)], &I1);
		// ---
		
		// update x
		memcpy(&x_prev[0], &x[0], mN*sizeof(double));
		if(fy < fx)
		{
			// TODO: use lambdaS and N_exit
			if( fx - fy < mRelativeError )
				stop_condition_reached = true;
				
			memcpy(&x[0], &y[0], mN*sizeof(double));
			fx = fy;	
		}
		
		// update the step value lambda_i 
		daxpy_(&mN, &minus_one,	&x[0], &I1, &x_prev[0],	&I1);
		lambda_i = 0.3 * dnrm2_(&mN, &x_prev[0], &I1) + 0.091 * mLambda;
		lambdaS = 0.62 * lambda_i;
	
	
		num_iter++;
		*f = fx;
		if (num_iter > mMaxIterations) 
		{
			stop_condition_reached = true;
			return 1;
		}	
	}
	
	
	return 0;
}


// ----------------------------------------------------------------
void CDOSOptimizer::LineSearch(double *y, double *p, double step, double *fy)
{
	// Local variables
	// we use the workspace to store the three required points
	std::vector<double> F;
	std::vector<double*> X;
	
	double *x, f_tmp;
	x  = &mSpace[0];
	
	double step_size(step);
	
	// stopping criterion variables
	int iter( 0 );
	int MaxIter( 20 );
	bool keep_searching(true);
	
	
	while(keep_searching)
	{
		iter ++;
		if(iter > MaxIter)
		{
			keep_searching = false;
		}
		
		// try a new point (stored at the begining of the workspace)
		memcpy(x, &y[0], mN*sizeof(double));
		daxpy_(&mN, &step_size,	&p[0], &I1, x,	&I1);
		
		if(isFeasiblePoint(&mSpace[0]))
		{
			f_tmp = -mModel->computeLikelihood(x, mN, mTraceFun);
			if (f_tmp < *fy)
			{
				// update y
				*fy = f_tmp;
				memcpy(&y[0], x, mN*sizeof(double));
				
				// update the step value
				step_size *= 2.;
				
				// keep track of the point
				F.push_back(f_tmp);
				X.push_back(&mSpace[(F.size()%mN + 1)*mN]);
				memcpy(X.back(), x, mN*sizeof(double));
			}
			else
			{
				if(X.size() >= 3)
					keep_searching = false;
				else
					reduce_step(step_size, iter);
			}
		}
		else
		{
			reduce_step(step_size, iter);
		}
	}
	
	if (X.size() < 3)
		return;

	// compute the quadratic interpolation
	
	double d, d12, d13, d23;
	double f1, f2, f3, *x1, *x2, *x3;
	
	f3 = F.back(); F.pop_back();
	f2 = F.back(); F.pop_back();
	f1 = F.back(); F.pop_back();
	
	x3 = X.back(); X.pop_back();
	x2 = X.back(); X.pop_back();
	x1 = X.back(); X.pop_back();
	
	d12 = distance(x1, x2, &mSpace[0], mN);
	d13 = distance(x1, x3, &mSpace[0], mN);
	d23 = distance(x2, x3, &mSpace[0], mN);
	
	// we need d13 > d12 && d13 > d23
	double f_swapper, *x_swapper;
	if(d13 < d23)
	{
		if(d23 < d12) // d13 < d12 < d23
		{
			x_swapper = x1;
			f_swapper = f1;
			x1 = x2;
			f1 = f2;
			x2 = x_swapper;
			f2 = f_swapper;
		}
		else // d13 < d23 < d12
		{
			x_swapper = x2;
			f_swapper = f2;
			x2 = x3;
			f2 = f3;
			x3 = x_swapper;
			f3 = f_swapper;
		}
	}
	else if(d13 < d12) // d23 < d13 < d12
	{
		x_swapper = x2;
		f_swapper = f2;
		x2 = x3;
		f2 = f3;
		x3 = x_swapper;
		f3 = f_swapper;
	}
	
	// compute the position of the interpolation
	d = 0.5 * ( square(d12)*(f1-f3) + square(d13)*(f2-f1) ) / ( f2*d13 - f3*d12 + f1*(d12-d13) );
	
	// set y = x1 + d*p
	memcpy(y, x1, mN*sizeof(double));
	daxpy_(&mN, &d, p, &I1, y, &I1);
	
	// update fy
	*fy = -mModel->computeLikelihood(y, mN, mTraceFun);
}


// ----------------------------------------------------------------
bool CDOSOptimizer::isFeasiblePoint(double *x) const
{
	// check all the bounds
	for(int i(0); i<mN; i++)
	{
		if(mLowerBound[i] > x[i] || mUpperBound[i] < x[i])
		{
			return false;
		}
	}
	return true;
}

// ----------------------------------------------------------------
void CDOSOptimizer::reduce_step(double& step, int const& numIter) const
{
	if(numIter <= 6)
		step /= 1.1;
	else if(numIter <= 8)
		step /= 1.2;
/*	else if(numIter <= 10)
		step /= 1.5;
	else if(numIter <= 16)
		step *= 0.5;
	else if(numIter <= 20)
		step *= 0.2;
	else if(numIter <= 40)
		step *= 0.1;
	else if(numIter <= 50)
		step *= 0.01; */
	else
		step = -5.*step;
}



// ----------------------------------------------------------------
// ----------------------------------------------------------------
void CDOSOptimizer::QRdecomposition(int width)
{
	// copy the 'width' first vectors of U into Q
	memcpy(&mQ[0], &mU[0], mN*width*sizeof(double));
	
	// size of the space to use
	const int lwork = mlwork;
	int info;
	
	// use the lapack functions to first decompose mQ
	dgeqrf(&mN, &width, &mQ[0], &mN, &mSpace[lwork], &mSpace[0], &lwork, &info);
	//TODO: check errors
	// and then compute the first 'width' vectors of Q
	dorgqr(&mN, &width, &width, &mQ[0], &mN, &mSpace[lwork], &mSpace[0], &lwork, &info);
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
	if (mN <= 5)
		mSpace.resize(6*mN);
	else
		mSpace.resize(mlwork+mN);
}

// ----------------------------------------------------------------
void CDOSOptimizer::computeGradient(double f0, double *x0, double *g)
{
	double delta( sqrt( DBL_MIN ) );
	double di, *x;
	x = &mSpace[0];
	
	for(int i(0); i < mN; i++)
	{
		memcpy(x, x0, mN*sizeof(double));
		
		di = delta * (1.+x0[i]);
		x[i] += di;
		
		if(x[i] > mUpperBound[i])
		{
			di = -di;
			x[i] += 2.*di;
		}
		
		g[i] = (-mModel->computeLikelihood(x, mN, false) - f0) / di;
	}	
}

