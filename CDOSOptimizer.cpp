#include "CDOSOPtimizer.h"

inline static double square(double a) {return a*a;}

// ----------------------------------------------------------------
double CDOSOptimizer::maximizeFunction(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	double lnL(0.0);
	CDOSminimizer(lnL, &aVars[0]);
	return -lnL; 
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------
int CDOSOptimizer::CDOSminimizer(double& f, double x[])
{
	// stage 1
	
	InitSearchDirections();
	// first conjugate direction is the anti-gradient normalized
	std::vector<double> gradient;
	gradient.resize(mN, 0.);
	// TODO compute the gradient and rescale it
	
	LineSearch(x, mU[0]);
	
	// stage 2
	double mLambdaS (0.62*mLambda);
	
	for(int i(1); i<mN; i++)
	{
		// QR decomposition to find the next orthogonal shift direction
		QRdecomposition(i);
		
		// y = x+lamndaS*mqi
		
		for(int j(0); j<i; j++)
		{
			LineSearch(y, mU[j*mN]);
		}
		// update ui
		// TODO
	}
	
	// stage 3
	// TODO
}


// ----------------------------------------------------------------
// ----------------------------------------------------------------
void QRdecomposition(int width)
{
	// copy first the 'width' first vectors of U into Q
	memcpy(mQ, mU, mN*width*sizeof(double));
	
	// local variables
	double s(0.);
	
	for(int j(0); j<mN; j++)
	{
		s = 0.;
		for(int i(j); j<width; j++)	{s += square(mQ[j*mN+i]); }
		s = sqrt(s);
		double ajj( mQ[(j+1)*mN] );
		double dj( ajj>0 ? -s : s );
		double fak( s*(s+fabs(ajj)) );
		mQ[(j+1)*mN] = ajj - dj;
		
		for(int k(j); k<width; k++) {mQ[j*mN+k] /= fak;}
		for(int i(j+1); i<mN; i++)
		{
			s = 0.;
			for(int k(j); k<width; k++) {s += mQ[j*mN+k] * mQ[i*mN+k];}
			for(int k(j); k<width; k++) {mQ[i*mN+k] -= mQ[j*mN+k] * s;}
		}
	}
	// compute now the vector qi (position 'width')
	for(int j(0); j<mN; j++) {mqi[j] = (j==width) ? 1. : 0.;}
	for(int j(0); j<mN; j++)
	{			
		s = (j>=width-1) ? mQ[(width+1)*mN] : 0.;
		for(int k(j); k<width; k++){mqi[k] += mQ[j*mN+k]*s; }
	}
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

