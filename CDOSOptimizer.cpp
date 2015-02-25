
#include "CDOSOptimizer.h"
#include "blas.h"
#include "lapack.h"


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
	alocateWorkspace();
	initSearchDirections();
	
	// local variables
	double fx, fy, lambda_i;
	std::vector<double> y(mN);
	std::vector<double> x_prev(mN);
	std::vector<double> traceVar(mN);
	
	// ------------------------------ stage 1 ------------------------------
	
	
	// first conjugate direction is the anti-gradient normalized
	std::vector<double> gradient;
	gradient.resize(mN);
	
	fy = -mModel->computeLikelihood(x, mN, mTraceFun);
	fx = fy;
	
	/*
	// compute the gradient at point x
	computeGradient(fx, x, &gradient[0]);
	*/
	
	/*
	// compute the quasi gradient using linear search
	memcpy(&y[0], x, mN*sizeof(double));
	double delta_search;
	for(int i(0); i<mN; i++)
	{
		delta_search = PerformLineSearch(&y[0], &mU[i*mN], mLambda, &fy);
		//gradient[i] = fabs(delta_search)>1e-6 ? (fy-fx)/delta_search : 0.;
		gradient[i] = delta_search>0 ? fy-fx : fx-fy;
		memcpy(&y[0], x, mN*sizeof(double));
		fy = fx;
	}
	*/
	
	// compute the quasi gradient
	double sign;
	for(int i(0); i<mN; i++)
	{
		memcpy(&y[0], x, mN*sizeof(double));
		y[i] += mLambda;
		sign = 1.;
		if(y[i] > mUpperBound[i])
		{
			y[i] -= 2.*mLambda;
			sign = -1.;
		}
		if(y[i] < mLowerBound[i])
		{
			// TODO lambda too high error
			return -2;
		}
		fy = -mModel->computeLikelihood(&y[0], mN, false);
		gradient[i] = (fy - fx)*sign;
	}
	
	
#ifdef USE_LAPACK
	double scale = -1./dnrm2_(&mN, &gradient[0], &I1);
	dscal_(&mN, &scale, &gradient[0], &I1);
#else
	double scale = 0.;
	for(int i=0; i < mN; ++i) scale += square(gradient[i]);
	scale = -1./sqrt(scale);
	for(int i=0; i < mN; ++i) gradient[i] *= scale;
#endif
	
	memcpy(&mU[0], &gradient[0], mN*sizeof(double));
	
	PerformLineSearch(x, &mU[0], mLambda, &fx);
	
	if (mVerbose > 2)
	{
		std::cout << "Lambda 0 : " << mLambda << '\n';
		memcpy(&traceVar[0], x, mN*sizeof(double));
		std::cout << "Stage 1 done\n";
		mModel->printVar(traceVar, -fx);
	}
	
	
	
	// ------------------------------ stage 2 ------------------------------
	
	double lambdaS (0.62*mLambda);
	
	for(int i(1); i<mN; i++)
	{
		
		// QR decomposition to find the next orthogonal shift direction mqi
		QRdecomposition(i);
		
		
		if(mVerbose > 2)
		{
			std::cout << "Stage 2, iteration i=" << i << "/" << mN << ", with loglikelihood " << -fx << "\n"
					  << "Direction mqi: \n";
			memcpy(&traceVar[0], &mqi[0], mN*sizeof(double));
			mModel->printVar(traceVar, -fx);
			
		}
		
		
		
		// perform y = x+lamndaS*mqi using blas
		memcpy(&y[0], &x[0], mN*sizeof(double));
		daxpy_(&mN, &lambdaS, &mqi[0], &I1, &y[0], &I1);
		
		// look for a new point y (better than x)
		for(int j(0); j<i-1; j++)
		{	
			// line search in uj direction
			// update y if better point
			PerformLineSearch(&y[0], &mU[j*mN], lambdaS, &fy);
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
	
	if (mVerbose > 2)
	{
		memcpy(&traceVar[0], x, mN*sizeof(double));
		std::cout << "Stage 2 done\n";
		mModel->printVar(traceVar, -fx);
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
			PerformLineSearch(&y[0], &mU[j*mN], 3.*lambdaS, &fy);
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


double CDOSOptimizer::PerformLineSearch(double *y, double *p, double step, double *fy)
{
#ifndef USE_CODE_ML_LINE_SEARCH
	return LineSearch(y, p, step, fy);
#else
	
	memcpy(&mSpace[0], y, mN*sizeof(double)); 
	
	double am(0);
	for (int i(0), am = 0.1; i < mN; i++)  	/* max step length */
    {
        if (p[i] > 0 && (mUpperBound[i] - y[i]) / p[i] < am)
        {
            am = (mUpperBound[i] - y[i]) / p[i];
        }
        else if (p[i] < 0 && (mLowerBound[i] - y[i]) / p[i] < am)
        {
            am = (mLowerBound[i] - y[i]) / p[i];
        }
    }
    
    step = am/2.;
	step = LineSearch2(fy, &mSpace[0], p, step, am, 1e-3, &mSpace[mN], 10, mN);
	daxpy_(&mN, &step, p, &I1, y, &I1);
	return 0.; //TODO
		
#endif
}



#ifndef USE_CODE_ML_LINE_SEARCH
// ----------------------------------------------------------------
double CDOSOptimizer::LineSearch(double *y, double *p, double step, double *fy)
{
	// Local variables
	// storage of the current optimal points (we store the steps in dX)
	// and their values in F
	std::vector<double> F;
	std::vector<double> dX;
	
	double *x, f_tmp, *x0;
	x  = &mSpace[0];
	x0 = &mSpace[mN];
	memcpy(x0, y, mN*sizeof(double));
	
	double step_size(step);
	double best_step_size(step);
	
	// stopping criterion variables
	int iter( 0 );
	int MaxIter( 10 );
	bool keep_searching(true);
	
	
	while(keep_searching)
	{
		iter ++;
		if(iter > MaxIter)
		{
			keep_searching = false;
		}
		
		// try a new point (stored at the begining of the workspace)
		// located at x = x0 + step_size*p
		
		memcpy(x, x0, mN*sizeof(double));
		daxpy_(&mN, &step_size,	&p[0], &I1, x,	&I1);
		
		int num_iter_constraint(0);
		while(!isFeasiblePoint(x) && num_iter_constraint < 70)
		{
			reduce_step(step_size, num_iter_constraint);
			if (step_size < 0.)
				step_size = - step;
			num_iter_constraint++;
			memcpy(x, x0, mN*sizeof(double));
			daxpy_(&mN, &step_size,	&p[0], &I1, x,	&I1);
		}
		
		f_tmp = -mModel->computeLikelihood(x, mN, mTraceFun);
		//std::cout << "debug: -- step size: " << step_size << '\n';
		if (f_tmp < *fy)
		{
			//std::cout << "debug: -- f_tmp: " << f_tmp << "\n";
			// update y
			*fy = f_tmp;
			memcpy(y, x, mN*sizeof(double));
			
			// keep track of the point
			F.push_back(f_tmp);
			dX.push_back(step_size);
			
			// update the step value
			best_step_size = step_size;
			step_size *= 2.;
		}
		else
		{
			if(iter == 1)
				step_size = -step;
			else if(dX.size() >= 3)
				keep_searching = false;
			else
				step_size *= 2.;
		}
	}
	
	if (dX.size() < 3)
		return best_step_size;

	// compute the quadratic interpolation
	
	double d, d12, d13, d23, denominator;
	double f1, f2, f3, x1, x2, x3;
	
	f3 = F.back(); F.pop_back();
	f2 = F.back(); F.pop_back();
	f1 = F.back(); F.pop_back();
	
	x3 = dX.back(); dX.pop_back();
	x2 = dX.back(); dX.pop_back();
	x1 = dX.back(); dX.pop_back();
	
	d12 = fabs(x2-x1);
	d13 = fabs(x3-x1);
	d23 = fabs(x3-x2);
	
	std::cout << "debug: -- calculation of d: ----------------------------- \n";
	printf("d12 %7.4f d13 %6.4f d23 %8.4f \n", d12, d13, d23);
	
	// we need d13 > d12 && d13 > d23
	double f_swapper, x_swapper;
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
	
	d12 = fabs(x2-x1);
	d13 = fabs(x3-x1);
	d23 = fabs(x3-x2);
	
	printf("d12 %7.4f d13 %6.4f d23 %8.4f \n", d12, d13, d23);
	
	// compute the position of the interpolation
	
	denominator = f2*d13 - f3*d12 + f1*(d12-d13);
	if( fabs(denominator) > 1e-6 )
	{
		d = 0.5 * ( square(d12)*(f1-f3) + square(d13)*(f2-f1) ) / denominator;
	
		// set y = x1 + d*p = x0 + (d+x1)p
		d += x1;
		memcpy(y, x0, mN*sizeof(double));
		daxpy_(&mN, &d, p, &I1, y, &I1);
	
		// update fy
		*fy = -mModel->computeLikelihood(y, mN, mTraceFun);
		std::cout << "debug: -- obtained d: " << d << " ----------------------------- \n";
		return d;
	}
	return best_step_size;
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
	else if(numIter <= 10)
		step /= 1.5;
	else if(numIter <= 16)
		step *= 0.5;
	else if(numIter <= 20)
		step *= 0.2;
	else if(numIter <= 40)
		step *= 0.1;
	else if(numIter <= 50)
		step *= 0.01;
	else
		step = -1.;
}
#endif



#ifdef USE_CODE_ML_LINE_SEARCH

double CDOSOptimizer::fun_LineSearch(double t, const double x0[], const double p[], double x[], int n)
{
    for(int i=0; i < n; ++i) x[i] = x0[i] + t * p[i];
    //return ((*fun) (x, n));
	return -mModel->computeLikelihood(x, n, mTraceFun);
}



double CDOSOptimizer::LineSearch2(double *f, const double x0[], const double p[], double step, double limit, double e, double space[], int iround, int n)
{
    /* linear search using quadratic interpolation
       from x0[] in the direction of p[],
                    x = x0 + a*p        a ~(0,limit)
       returns (a).    *f: f(x0) for input and f(x) for output

       x0[n] x[n] p[n] space[n]

       adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
       optimization: An introduction.  Van Nostrand Reinhold Company, New York.
       pp. 62-73.
       step is used to find the bracket and is increased or reduced as necessary,
       and is not terribly important.
    */
    int ii = 0, maxround = 10, status, i, nsymb = 0;
    double *x = space, factor = 4, small = 1e-10, smallgapa = 0.2;
    double a0, a1, a2, a3, a4 = -1, a5, a6, f0, f1, f2, f3, f4 = -1, f5, f6;

    /* look for bracket (a1, a2, a3) with function values (f1, f2, f3)
       step length step given, and only in the direction a>=0
    */

    if (mVerbose > 2)
        printf("\n%3d h-m-p %7.4f %6.4f %8.4f ", iround + 1, step, limit, norm(p, n));

    if (step <= 0 || limit < small || step >= limit)
    {
        if (mVerbose > 2)
            printf("\nh-m-p:%20.8e%20.8e%20.8e %12.6f\n", step, limit, norm(p, n), *f);

        return 0;
    }

    a0 = a1 = 0;
    f1 = f0 = *f;
    a2 = a0 + step;
    f2 = fun_LineSearch(a2, x0, p, x, n);

    if (f2 > f1)  		/* reduce step length so the algorithm is decreasing */
    {
        for (;;)
        {
            step /= factor;

            if (step < small)
            {
                return 0;
            }

            a3 = a2;
            f3 = f2;
            a2 = a0 + step;
            f2 = fun_LineSearch(a2, x0, p, x, n);

            if (f2 <= f1)
            {
                break;
            }

            if (mVerbose > 2)   //CMV added correct #if
            {
                putchar('-');
                nsymb++;
            }
        }
    }
    else  			/* step length is too small? */
    {
        for (;;)
        {
            step *= factor;

            if (step > limit)
            {
                step = limit;
            }

            a3 = a0 + step;
            f3 = fun_LineSearch(a3, x0, p, x, n);

            if (f3 >= f2)
            {
                break;
            }

            if (mVerbose > 2)
            {
                putchar('+');
                nsymb++;
            }

            a1 = a2;
            f1 = f2;
            a2 = a3;
            f2 = f3;

            if(step >= limit)
            {
                if(mVerbose > 2) //CMV added correct #if
				{
                    for (; nsymb < 5; nsymb++)
                    {
                        printf(" ");
                    }

                    printf(" %12.6f%3c %6.4f", *f = f3, 'm', a3);
                    //printf(" %12.6f%3c %6.4f %5d", *f = f3, 'm', a3, mNumFunCall);
				}

                *f = f3;
                return a3;
            }
        }
    }

    /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
    for(ii = 0; ii < maxround; ii++)
    {
        /* a4 is the minimum from the parabola over (a1,a2,a3)  */
        a4 = (a2 - a3) * f1 + (a3 - a1) * f2 + (a1 - a2) * f3;

        if (fabs(a4) > 1e-100)
            a4 = ((a2 * a2 - a3 * a3) * f1 + (a3 * a3 - a1 * a1) * f2 + (a1 * a1 - a2 * a2) * f3) / (2 * a4);

        if (a4 > a3 || a4 < a1)  	/* out of range */
        {
            a4 = (a1 + a2) / 2;
            status = 'N';
        }
        else
        {
            if ((a4 <= a2 && a2 - a4 > smallgapa * (a2 - a1))
            || (a4 > a2 && a4 - a2 > smallgapa * (a3 - a2)))
            {
                status = 'Y';
            }
            else
            {
                status = 'C';
            }
        }

        f4 = fun_LineSearch(a4, x0, p, x, n);

        if (mVerbose > 2) //CMV added correct #if
        {
            putchar(status);
        }

        if (fabs(f2 - f4) < e * (1 + fabs(f2)))
        {
            if (mVerbose > 2) //CMV added correct #if
                for (nsymb += ii + 1; nsymb < 5; nsymb++)
                {
                    putchar(' ');
                }
            break;
        }

        /* possible multiple local optima during line search */
        if (mVerbose > 2 && ((a4 < a2 && f4 > f1) || (a4 > a2 && f4 > f3)))   //CMV added correct #if
        {
            printf("\n\na %12.6f %12.6f %12.6f %12.6f", a1, a2, a3, a4);
            printf("\nf %12.6f %12.6f %12.6f %12.6f\n", f1, f2, f3, f4);

            for (a5 = a1; a5 <= a3; a5 += (a3 - a1) / 20)
            {
                printf("\t%.6e ", a5);

                if (n < 5)
                {
                    for(i=0; i < n; ++i) printf("\t%.6f", x0[i] + a5 * p[i]);
                }

                printf("\t%.6f\n", fun_LineSearch(a5, x0, p, x, n));
            }

            puts("Linesearch2 a4: multiple optima?");
        }

        if (a4 <= a2)  		/* fig 2.2.10 */
        {
            if (a2 - a4 > smallgapa * (a2 - a1))
            {
                if (f4 <= f2)
                {
                    a3 = a2;
                    a2 = a4;
                    f3 = f2;
                    f2 = f4;
                }
                else
                {
                    a1 = a4;
                    f1 = f4;
                }
            }
            else
            {
                if (f4 > f2)
                {
                    a5 = (a2 + a3) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 > f2)
                    {
                        a1 = a4;
                        a3 = a5;
                        f1 = f4;
                        f3 = f5;
                    }
                    else
                    {
                        a1 = a2;
                        a2 = a5;
                        f1 = f2;
                        f2 = f5;
                    }
                }
                else
                {
                    a5 = (a1 + a4) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 >= f4)
                    {
                        a3 = a2;
                        a2 = a4;
                        a1 = a5;
                        f3 = f2;
                        f2 = f4;
                        f1 = f5;
                    }
                    else
                    {
                        a6 = (a1 + a5) / 2;
                        f6 = fun_LineSearch(a6, x0, p, x, n);

                        if (f6 > f5)
                        {
                            a1 = a6;
                            a2 = a5;
                            a3 = a4;
                            f1 = f6;
                            f2 = f5;
                            f3 = f4;
                        }
                        else
                        {
                            a2 = a6;
                            a3 = a5;
                            f2 = f6;
                            f3 = f5;
                        }
                    }
                }
            }
        }
        else  		/* fig 2.2.9 */
        {
            if (a4 - a2 > smallgapa * (a3 - a2))
            {
                if (f2 >= f4)
                {
                    a1 = a2;
                    a2 = a4;
                    f1 = f2;
                    f2 = f4;
                }
                else
                {
                    a3 = a4;
                    f3 = f4;
                }
            }
            else
            {
                if (f4 > f2)
                {
                    a5 = (a1 + a2) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 > f2)
                    {
                        a1 = a5;
                        a3 = a4;
                        f1 = f5;
                        f3 = f4;
                    }
                    else
                    {
                        a3 = a2;
                        a2 = a5;
                        f3 = f2;
                        f2 = f5;
                    }
                }
                else
                {
                    a5 = (a3 + a4) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 >= f4)
                    {
                        a1 = a2;
                        a2 = a4;
                        a3 = a5;
                        f1 = f2;
                        f2 = f4;
                        f3 = f5;
                    }
                    else
                    {
                        a6 = (a3 + a5) / 2;
                        f6 = fun_LineSearch(a6, x0, p, x, n);

                        if (f6 > f5)
                        {
                            a1 = a4;
                            a2 = a5;
                            a3 = a6;
                            f1 = f4;
                            f2 = f5;
                            f3 = f6;
                        }
                        else
                        {
                            a1 = a5;
                            a2 = a6;
                            f1 = f5;
                            f2 = f6;
                        }
                    }
                }
            }
        }
    }

    if (f2 > f0 && f4 > f0)
    {
        a4 = 0;
    }

    if (f2 <= f4)
    {
        *f = f2;
        a4 = a2;
    }
    else
    {
        *f = f4;
    }

    if (mVerbose > 2)
    {
        printf(" %12.6f%3d %6.4f", *f, ii, a4);
        //printf(" %12.6f%3d %6.4f %5d", *f, ii, a4, mNumFunCall);
    }

    return a4;
}

#endif

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
	if( info < 0 )
	{
		//TODO: check errors exceptions
		std::cout << "ERROR: dgeqrf: error at column " << -info << "\n";
	}
	
	// and then compute the first 'width' vectors of Q
	dorgqr(&mN, &width, &width, &mQ[0], &mN, &mSpace[lwork], &mSpace[0], &lwork, &info);
	if( info < 0 )
	{
		//TODO: check errors exceptions
		std::cout << "ERROR: dorgqr: error at column " << -info << "\n";
	}
	
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
	double delta( 1e-7 );
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

