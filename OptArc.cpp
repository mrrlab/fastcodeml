
#include "OptArc.h"
#include "blas.h"
#include "lapack.h"
#include <cfloat>
#include <iomanip>
#include <assert.h>


// ----------------------------------------------------------------------
//	Class members definition: OptArc
// ----------------------------------------------------------------------
double OptArc::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	allocateMemory();
	
	double maxl = 1e7;
	ArcMinimizer(&maxl, &aVars[0]);
	return -maxl;
}


// ----------------------------------------------------------------------
void OptArc::allocateMemory(void)
{
	size_vect = mN*sizeof(double);
	
	mXEvaluator.resize(mN);
	mSpace.resize(2*mN*mN + 7*mN);	
	
	mWorkSpaceVect	= &mSpace[0];
	mWorkSpaceMat	= mWorkSpaceVect + mN;
	
	mGradient	= mWorkSpaceMat	+ mN*mN;
	mP			= mGradient		+ mN;
	mHessian	= mP			+ mN;
	
	mSk = mHessian	+ mN*mN;
	mYk = mSk		+ mN;
	
	mXPrev		= mYk		+ mN;
	mGradPrev	= mXPrev	+ mN;
	
	mArcSpace.resize( mNs*mNs + mNs*mN + mNs );
	mQ = &mArcSpace[0];
	mU = mQ+mNs*mN;
	mV = mU+mNs*mNs;
	
	// set active set counters to zero
	mActiveSet.resize(mN, 0);
}


// ----------------------------------------------------------------------
void OptArc::ArcMinimizer(double *f, double *x)
{
	*f = evaluateFunction(x, mTrace);
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
	{
		std::cout << "Solution after Bootstrap:";
		mModel->printVar(mXEvaluator, *f);
	}
	
	// compute current gradient
	computeGradient(x, *f, mGradient);
	
	// initialize hessian matrix to identity
	int N_sq( mN*mN ), diag_stride( mN+1 );
	dcopy_(&N_sq, &D0, &I0, mHessian, &I1);
	dcopy_(&mN, &D1, &I0, mHessian, &diag_stride);
	
	// --- main loop
	bool convergenceReached = false;
	for(mStep = 0; !convergenceReached; ++mStep)
	{
		// save current parameters
		double f_prev = *f;
		memcpy(mGradPrev, mGradient, size_vect);
		memcpy(mXPrev, x, size_vect);
		
		// prepare the arc search
		computeSubspaceArcSearch(x);
		
		// arc search
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Arc Search..." << std::endl;
		double alpha_max = (fabs(mLambdaMin) > 1e-1) ? 1.0/fabs(mLambdaMin) : 1.0;
		double alpha = alpha_max;
		
		// find the maximum step size alpha so it stays in the bounds
		findMaxStep(x, &alpha);
		alpha = min2(alpha, alpha_max);
		std::cout << std::setprecision(14) <<  "alpha_max = " << alpha << std::endl;
		std::cout << "1/vmin = " << 1.0/mLambdaMin << std::endl;
		
		arcSearch(&alpha, x, f);
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
			std::cout << "Step length found:" << alpha << std::endl;
		
		// avoid unsatisfied bounds due to roundoff errors 
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			double x_ = x[i];
			x_ = max2(x_, mLowerBound[i]);
			x_ = min2(x_, mUpperBound[i]);
			x[i] = x_;
		}
		
		if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		{
			std::cout << "New Solution:";
			mModel->printVar(mXEvaluator, *f);
		}		
		
		
		// check convergence
		convergenceReached =   fabs(f_prev - *f) < mAbsoluteError
							|| mStep >= mMaxIterations;
		
		if (!convergenceReached)
		{
			// update the system
			
			computeGradient(x, *f, mGradient);
				
			memcpy(mSk, x, size_vect);
			daxpy_(&mN, &minus_one, mXPrev, &I1, mSk, &I1);
		
			memcpy(mYk, mGradient, size_vect);
			daxpy_(&mN, &minus_one, mGradPrev, &I1, mYk, &I1);
		
			if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
				std::cout << "SR1 update..." << std::endl;
		
			// update the B matrix
			SR1update();
			
			// update the active set
			updateActiveSet(x);
		}
	}
}


// ----------------------------------------------------------------------
double OptArc::evaluateFunction(const double *x, bool aTrace)
{
	memcpy(&mXEvaluator[0], x, size_vect);
	double f = mModel->computeLikelihood(mXEvaluator, aTrace);
	
	// Stop optimization if value is greater or equal to threshold
	if (mStopIfBigger && f >= mThreshold) throw FastCodeMLEarlyStopLRT();
	return -f;
}


// ----------------------------------------------------------------------
double OptArc::evaluateFunctionForArcSearch(const double* x, double alpha)
{
	// compute rho(mV, alpha)
	double *rho = mWorkSpaceVect;
	for (size_t j(0); j<mNs; ++j)
		rho[j] = alpha / (1.0 + alpha*(mV[j]-mLambdaMin));
	
	// introduce perturbation if the gradient and search direction are too similar
	double *g_perturbed = mWorkSpaceMat;
	double *Bp = g_perturbed+mN;
	memcpy(g_perturbed, mGradient, size_vect);
	
	//dscal_(&mN, &minus_one, g_perturbed, &I1);
	//projectedDirection(x, g_perturbed);
	//dscal_(&mN, &minus_one, g_perturbed, &I1);
	projectActiveSet(g_perturbed);
	
	char trans = 'N';
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mGradient, &I1, &D0, Bp, &I1);
	double dp = fabs(ddot_(&mN, mP, &I1, mGradient, &I1));
	double pBp = fabs(ddot_(&mN, mP, &I1, Bp, &I1));
	if (dp > 1e-4*pBp)
	{
		// add perturbation
		daxpy_(&mN, &D1, mP, &I1, g_perturbed, &I1);
	}
	
	// compute gamma(alpha)
	double *w = g_perturbed+mNs;
	double *w_tmp = w+mN;
	trans = 'T';
	dgemv_(&trans, &mN, &mNs, &D1, mQ, &mN, g_perturbed, &I1, &D0, w, &I1);	// multiply by QT
	dgemv_(&trans, &mNs, &mNs, &D1, mU, &mNs, w, &I1, &D0, w_tmp, &I1);		// multiply by UT
	for (int i=0; i<mNs; ++i) {w_tmp[i] *= rho[i];}						// multiply by diagonal matrix rho(V, alpha)
	trans = 'N';
	double *x_arc = mWorkSpaceVect;
	dgemv_(&trans, &mNs, &mNs, &D1, mU, &mNs, w_tmp, &I1, &D0, w, &I1);		// multiply by U
	dgemv_(&trans, &mN, &mNs, &minus_one, mQ, &mN, w, &I1, &D0, x_arc, &I1);// multiply by Q
	
	// compute f(x+gamma(alpha))
	daxpy_(&mN, &D1, x, &I1, x_arc, &I1);
	
#if 1
	// project on the search space
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		double x_ = x_arc[i];
		x_ = max2(x_, mLowerBound[i]);
		x_ = min2(x_, mUpperBound[i]);
		x_arc[i] = x_;
	}
#endif
	return evaluateFunction(x_arc, mTrace);
}


// ----------------------------------------------------------------------
void OptArc::computeGradient(const double *x, double f0, double *aGrad)
{
	volatile double eh;
	double sqrt_eps = sqrt(DBL_EPSILON);
	double f;
	memcpy(&mXEvaluator[0], x, size_vect);
	size_t i;
	double *delta = mWorkSpaceVect;
	
	const double *x_ = x;
	
	// branch lengths
	for(i=0; i<mNumTimes; ++i)
	{
		eh = sqrt_eps * ( 1.0 + x_[i] );
		if( x_[i] + eh > mUpperBound[i] )
			eh = -eh;
		mXEvaluator[i] += eh;
		delta[i] = mXEvaluator[i] - x_[i];
	}

	for(i=0; i<mNumTimes; ++i)
	{
		if(mActiveSet[i] <= 1)
		{
			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-f0)/delta[i];
		}
		// otherwise we don't change it
	}
	
	// other variables
	memcpy(&mXEvaluator[0], x_, size_vect);
	for(; i<mN; ++i)
	{
		if (mActiveSet[i] <= 1)
		{
			eh = sqrt_eps * ( 1.0 + fabs(x_[i]) );
			if ( x_[i] + eh > mUpperBound[i] )
				eh = -eh;
			mXEvaluator[i] += eh;
			eh = mXEvaluator[i] - x_[i];

			f = -mModel->computeLikelihoodForGradient(mXEvaluator, false, i);
			aGrad[i] = (f-f0)/eh;
			mXEvaluator[i] = x_[i];
		}
	}
}


// ----------------------------------------------------------------------
void OptArc::SR1update(void)
{
	// local variables
	double *v, *Bs;
	double vs;
	double *vvT;
	char trans = 'N';
	
	int mN_sq = mN*mN;
	
	// compute vector B*mSk
	Bs = mWorkSpaceVect;
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, mSk, &I1, &D0, Bs, &I1); 	
	
	// compute vector v = -y + Bs
	v = Bs;
	daxpy_(&mN, &minus_one, mYk, &I1, v, &I1);
	
	vs = -ddot_(&mN, v, &I1, mSk, &I1);
	
	// only update if vs is big enough
	if (fabs(vs) > 1e-8)
	{
		double inverse_vs = 1.0/vs;
	
		// compute Matrix v.v^T / vs
		vvT = mWorkSpaceMat;
		#pragma omp parallel for
		for(int i=0; i<mN; ++i)
		{
			double prefactor = v[i] * inverse_vs;
			dcopy_(&mN, v, &I1, &vvT[i*mN], &I1);
			dscal_(&mN, &prefactor, &vvT[i*mN], &I1);
		}
	
		// add the v.v^T / vs contribution if it is not too big compared to B_prev
		
		double Frob_prev = dnrm2_(&mN_sq, mHessian, &I1);
		double Frob_modif = dnrm2_(&mN_sq, vvT, &I1);
		
		if (Frob_modif < 1e8 * (Frob_prev+1.0))
			daxpy_(&mN_sq, &D1, vvT, &I1, mHessian, &I1);
	}
}


// ----------------------------------------------------------------------
void OptArc::arcSearch(double *aalpha, double *x, double *f)
{
	// constants for the Wolfe condition
	double c1 (2e-1);
	
	// compute derivative of arc at alpha=0
	double *QTg = mWorkSpaceMat;
	double *QQTg = mWorkSpaceVect;
	char trans = 'T';
	dgemv_(&trans, &mN, &mNs, &D1, mQ, &mN, mGradient, &I1, &D0, QTg, &I1);
	trans = 'N';
	dgemv_(&trans, &mN, &mNs, &D1, mQ, &mN, QTg, &I1, &D0, QQTg, &I1);
	double phi_0_prime = -ddot_(&mN, QQTg, &I1, mGradient, &I1);
	double phi_0 = *f, phi, phi_prev;
	double a_prev = 0.0;
	double phi_a_prime;
	double a = *aalpha;
	
	if (mVerbose >= VERBOSE_MORE_INFO_OUTPUT)
		std::cout << "phi_0_prime: " << phi_0_prime << std::endl; 
	
	phi_prev = phi_0;
	phi = evaluateFunctionForArcSearch(x, a);
	
	std::cout << "\tDEBUG: phi(" << a << ") = " << phi << std::endl;
	
	double sigma, sigma_bas;
	int maxIterBack, maxIterUp;
	
	// we take a step dependant on the problem size:
	// if the problem is large, we are able to spend more time 
	// to find a better solution.
	// otherwise, we consider that the solution is sufficiently 
	// good and continue. 
	// The time of line search should be small compared to the 
	// gradient computation
	
	maxIterBack = maxIterUp = static_cast<int> (ceil( 3.*log(mN+10.) ));
	sigma_bas 	= pow(1e-3, 1./static_cast<double>(maxIterBack));
	
	
	// begin by a backtrace
	double a_min = a;
	double phi_min = phi_0;
	size_t iter = 0;
	while(phi > phi_0 + phi_0_prime*a*c1 && iter < maxIterBack)
	{
		++iter;
		a_prev = a;
		phi_prev = phi;
		sigma = sigma_bas * (0.9 + 0.2*randFrom0to1());
		a *= sigma;
		phi = evaluateFunctionForArcSearch(x, a);
		std::cout << "\tDEBUG: a = " << a << ", phi = " << phi << std::endl; 
		a_min = (phi<phi_min)?a:a_min;
	}
	
	a = a_min;
	a_prev = a;
	phi = phi_min;
	phi_prev = phi;
	
	// compute the derivative
	double eh = sqrt(DBL_EPSILON);
	if ( a+eh >= 1.0 ) {eh = -eh;}
	phi_a_prime = (phi - evaluateFunctionForArcSearch(x, a+eh))/eh;
	
	iter = 0;
	if (phi_a_prime < 0.0 && fabs(a - *aalpha) > 1e-8)
	{
		double a0 = a_prev;
		while(phi <= phi_prev && iter < maxIterBack)
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			sigma = sigma_bas * (0.55 + 0.3*randFrom0to1());
			a = a + sigma*(a*0.5+a0-a);
			phi = evaluateFunctionForArcSearch(x, a);
		}
		if (phi_prev < phi)
		{
			a = a_prev;
		}
	}
	else
	{
		sigma_bas = 0.7;
		while(phi <= phi_prev && iter < maxIterUp)
		{
			++iter;
			
			a_prev = a;
			phi_prev = phi;
			sigma = sigma_bas * (0.7 + 0.6*randFrom0to1());
			a *= sigma;
			phi = evaluateFunctionForArcSearch(x, a);
		}
		if (phi_prev < phi)
		{
			a = a_prev;
		}
	}
	
	*aalpha = a;
	*f = evaluateFunctionForArcSearch(x, a);
	memcpy(x, &mXEvaluator[0], size_vect);
}


// ----------------------------------------------------------------------
void OptArc::computeSubspaceArcSearch(const double *x)
{
	// compute the minimum eigenValue and its corresponding eigenvector 
	// in order to have the descent direction mP with negative curvature
	char job = 'V';
	char range = 'I';
	char uplo = 'U';
	double *B = mWorkSpaceMat;
	memcpy(B, mHessian, mN*size_vect);
	
	double *vl_not_used = NULL, *vu_not_used = NULL;
	int M;
	double *eigenValues = mWorkSpaceVect;
	int *isupz_not_used = NULL;
	
	int lwork = 26*mN;
	std::vector<double> work(lwork);
	int liwork = 10*mN;
	std::vector<int> iwork(liwork);
	int info;
	
	dsyevr_(&job
		   ,&range
		   ,&uplo
		   ,&mN
		   ,B
		   ,&mN
		   ,vl_not_used
		   ,vu_not_used
		   ,&I1
		   ,&I1
		   ,&mAbsoluteError
		   ,&M
		   ,eigenValues
		   ,mP
		   ,&mN
		   ,isupz_not_used
		   ,&work[0]
		   ,&lwork
		   ,&iwork[0]
		   ,&liwork
		   ,&info);
	
	if (info != 0)
		std::cout << "ERROR : dsyevr_ in OptArc::computeSubspaceArcSearch: info = " << info << std::endl;
	assert(info == 0);
	
	
	mLambdaMin = eigenValues[0];
	std::cout << "min eigenValue of the hessian: " << mLambdaMin << std::endl;
	if (mLambdaMin > 0.0)
	{
		// take the direction mP = - B^-1 g projected
		memcpy(B, mHessian, mN*size_vect);
		memcpy(mP, mGradient, size_vect);
		dgesv(&mN, &I1, B, &mN, &iwork[0], mP, &mN, &info);
		if (info != 0)
			std::cout << "ERROR : dgesv in OptArc::computeSubspaceArcSearch: info = " << info << std::endl;
		dscal_(&mN, &minus_one, mP, &I1);
		// project the search direction
		projectedDirection(x, mP);
	}
	else
	{
		// take the direction such that mP.mGradient < 0 (descent direction)
		if (ddot_(&mN, mP, &I1, mGradient, &I1) > 0.0)
			dscal_(&mN, &minus_one, mP, &I1);
		// project the negative curvature direction
		projectedDirection(x, mP);	
	}
	
	
	// form subspace S
	double *S = mQ;
	memcpy(S, mGradient, size_vect);
	dscal_(&mN, &minus_one, S, &I1);
	memcpy(S+mN, mP, size_vect);
	
	// only consider free variables
	//projectActiveSet(S);
	//projectActiveSet(S+mN);
	projectedDirection(x, S);
	projectedDirection(x, S+mN);
	
	// form the matrix Q from subspace S
	std::vector<double> tau(mNs);
	// workspace querry
	lwork = -1;
	dgeqrf(&mN, &mNs, S, &mN, &tau[0], &work[0], &lwork, &info);
	lwork = static_cast<int>(work[0]);
	work.resize(lwork);
	dgeqrf(&mN, &mNs, S, &mN, &tau[0], &work[0], &lwork, &info);
	if (info != 0)
		std::cout << "ERROR : dgeqrf in OptArc::computeSubspaceArcSearch: info = " << info << std::endl;
	assert(info == 0);
	
	// workspace querry
	lwork = -1;
	dorgqr(&mN, &mNs, &mNs, S, &mN, &tau[0], &work[0], &lwork, &info);
	lwork = static_cast<int>(work[0]);
	work.resize(lwork);
	dorgqr(&mN, &mNs, &mNs, S, &mN, &tau[0], &work[0], &lwork, &info);
	if (info != 0)
		std::cout << "ERROR : dorgqr in OptArc::computeSubspaceArcSearch: info = " << info << std::endl;
	assert(info == 0);
	
	
	// form the matrix QTHQ
	char transA = 'N', transB = 'N'; 
	dgemm_(&transA, &transB, &mN, &mNs, &mN
		  ,&D1, mHessian, &mN, mQ, &mN
		  ,&D0, mWorkSpaceMat, &mN); 
	transA = 'T';
	dgemm_(&transA, &transB, &mNs, &mNs, &mN
		  ,&D1, mQ, &mN, mWorkSpaceMat, &mN
		  ,&D0, mU, &mNs); 
	
	//std::cout << "Matrix QHQ : " << std::setprecision(12) << mU[0] << " " << mU[1] << " " << mU[2] << " " << mU[3] << std::endl;
	
	// compute eigenValues V and eigenVectors U
	// workspace quaerry
	lwork = -1; liwork = -1;
	dsyevd_(&job, &uplo, &mNs, mU, &mNs, mV, &work[0], &lwork, &iwork[0], &liwork, &info);
	lwork = static_cast<int>(work[0]);
	liwork = iwork[0];
	work.resize(lwork);
	iwork.resize(liwork);
	dsyevd_(&job, &uplo, &mNs, mU, &mNs, mV, &work[0], &lwork, &iwork[0], &liwork, &info);
	
	if (info != 0)
		std::cout << "ERROR : dsyevd_ in OptArc::computeSubspaceArcSearch: info = " << info << std::endl;
	assert(info == 0);
	
	mLambdaMin = mV[0];
	//std::cout << "eigenValues of QHQ : " << std::setprecision(12) << mV[0] << " " << mV[1] << std::endl;
	//std::cout << "eigenVectors QHQ : " << std::setprecision(12) << mU[0] << " " << mU[1] << " " << mU[2] << " " << mU[3] << std::endl;
}


// ----------------------------------------------------------------------
void OptArc::findMaxStep(const double *x, double *amax)
{
	double a = 1e16; // take a huge value to allow bigger values than a_max
	
	double *QU = mWorkSpaceMat;
	char trans = 'N'; 
	dgemm_(&trans, &trans, &mN, &mNs, &mNs
		  ,&D1, mQ, &mN, mU, &mNs
		  ,&D0, QU, &mN); 
	
	double *g_perturbed = mWorkSpaceMat;
	double *Bp = g_perturbed+mN;
	
	memcpy(g_perturbed, mGradient, size_vect);
	
	//dscal_(&mN, &minus_one, g_perturbed, &I1);
	//projectedDirection(x, g_perturbed);
	//dscal_(&mN, &minus_one, g_perturbed, &I1);
	projectActiveSet(g_perturbed);
	
	dgemv_(&trans, &mN, &mN, &D1, mHessian, &mN, g_perturbed, &I1, &D0, Bp, &I1);
	
	double dp = fabs(ddot_(&mN, mP, &I1, g_perturbed, &I1));
	double pBp = fabs(ddot_(&mN, g_perturbed, &I1, Bp, &I1));
	if (dp > 1e-4*pBp)
	{
		// add perturbation
		daxpy_(&mN, &D1, mP, &I1, g_perturbed, &I1);
	}
	
	
	double beta1_ = ddot_(&mN, &QU[0] , &I1, g_perturbed, &I1);
	double beta2_ = ddot_(&mN, &QU[mN], &I1, g_perturbed, &I1);
	const double tol_alpha = 1e-6;
	
	for (int i=0; i<mN; ++i)
	{
		// compute the limiting alpha for the ith variable
		double beta1, beta2, gamma, linear_term, const_term;
		double D;
		
		// lower bound
		beta1 = QU[i]*beta1_;
		beta2 = QU[mN+i]*beta2_;
		gamma = mLowerBound[i] - x[i];
		
		linear_term = gamma*(mV[0]+mV[1]) - beta1 - beta2;
		const_term  = gamma*mV[0]*mV[1] - beta1*mV[1] - beta2*mV[0];
		
		D = square(linear_term) - 4.0*gamma*const_term;
		if (D >= 0.0)
		{
			D = sqrt(D);
			double pi1 = (-linear_term + D)/(2.0*gamma);
			double pi2 = -(linear_term + D)/(2.0*gamma);
			double a1 = 1.0/(pi1+mLambdaMin);
			double a2 = 1.0/(pi2+mLambdaMin);
			
			//std::cout << i << " gamma = " << gamma << ", lin = " << linear_term << ", const = " << const_term << std::endl;
			//std::cout << i << " Lower: " << a1 << " " << a2 << std::endl;
			a2 = max2(a1,a2);
			
			if (a2 >= tol_alpha)
				a = min2(a2, a);		
		}		
		
		// upper bound
		beta1 = QU[i]*beta1_;
		beta2 = QU[mN+i]*beta2_;
		gamma = mUpperBound[i] - x[i];
		
		linear_term = gamma*(mV[0]+mV[1]) - beta1 - beta2;
		const_term  = gamma*mV[0]*mV[1] - beta1*mV[1] - beta2*mV[0];
		
		D = square(linear_term) - 4.0*gamma*const_term;
		if (D >= 0.0)
		{
			D = sqrt(D);
			double pi1 = (-linear_term + D)/(2.0*gamma);
			double pi2 = -(linear_term + D)/(2.0*gamma);
			double a1 = 1.0/(pi1+mLambdaMin);
			double a2 = 1.0/(pi2+mLambdaMin);
			
			//std::cout << i << " gamma = " << gamma << ", lin = " << linear_term << ", const = " << const_term << std::endl;
			//std::cout << i << " Upper: " << a1 << " " << a2 << std::endl;
			a2 = max2(a1,a2);
			
			if (a2 >= tol_alpha)
				a = min2(a2, a);
		}
	}
	*amax = a;
}


// ----------------------------------------------------------------------
void OptArc::updateActiveSet(const double *x)
{
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		double x_ = x[i];
		double g_ = mGradient[i];
		double l_ = mLowerBound[i];
		double u_ = mUpperBound[i];
		
		const double tol = 1e-4;
		
		if (x_ < l_+tol && g_ > 0.0)
		{
			mActiveSet[i] = 1;
		}
		else if (x_ > u_-tol && g_ < 0.0)
		{
			mActiveSet[i] = 1;
		}
		else
		{
			mActiveSet[i] = 0;
		}
	}
}


// ----------------------------------------------------------------------
void OptArc::projectActiveSet(double *aVect)
{
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		if (mActiveSet[i] != 0)
			aVect[i] = 0.0;
	}
}

// ----------------------------------------------------------------------
void OptArc::projectedDirection(const double *x, double *p)
{
	// set p <- x+p
	daxpy_(&mN, &D1, x, &I1, p, &I1);
	// projection step
	#pragma omp parallel for
	for (int i=0; i<mN; ++i)
	{
		double pi = p[i];
		pi = max2(pi, mLowerBound[i]);
		pi = min2(pi, mUpperBound[i]);
		p[i] = pi;
	}
	// set p <- P(x+g) - x
	daxpy_(&mN, &minus_one, x, &I1, p, &I1);
}

// ----------------------------------------------------------------------

