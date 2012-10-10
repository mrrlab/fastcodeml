
#include <iostream>
#include <iomanip>
#include <cmath>
#include "BayesTest.h"
#include "VerbosityLevels.h"
#include "TreeAndSetsDependencies.h"



BayesTest::BayesTest(size_t aNumSites, unsigned int aVerbose) : mSiteClassProb(BEB_DIMS*aNumSites), mNumSites(aNumSites), mVerbose(aVerbose), mPriors(aNumSites*BEB_NUM_CAT)
{
	// Searching range for w0
	const static double w0b0 = 0.;
	const static double w0b1 = 1.;

	// Searching range for w2
	const static double w2b0 = 1.;
	const static double w2b1 = 11.;

	// Initialize the w0 and w2 values to be tested
	for(unsigned int i=0; i < BEB_N1D; i++)
	{
		mPara[0][i] = mPara[1][i] = -1;						// p0 & p1
		mPara[2][i] = w0b0 + (i+0.5)*(w0b1-w0b0)/BEB_N1D;	// w0
		mPara[3][i] = w2b0 + (i+0.5)*(w2b1-w2b0)/BEB_N1D;	// w2
	}
}

double BayesTest::getGridParams(Forest& aForest, const std::vector<double>& aVars, const std::vector<double>& aSiteMultiplicity, size_t aFgBranch)
{
	// Omega for foreground and background branches and the corresponding codon class
	double omega_fg;
	double omega_bg;

	// Get the optimized kappa
	size_t kappa_idx = aVars.size()-2;
	const double kappa = aVars[kappa_idx];

	// Compute the corresponding Q matrices for foreground and background branches
	TransitionMatrix q_fg, q_bg;
	double scale_q_fg, scale_q_bg;

	// Create the dependency list for forest likelihood computation
	TreeAndSetsDependencies dep(aForest, mVerbose);
	dep.computeDependencies(1, true);

	// Initialize the probability list
	ProbabilityMatrixSetBEB beb_set(aForest.getNumBranches());
	beb_set.initializeSet(aForest.adjustFgBranchIdx(aFgBranch));

	// Calculating f(x_h|w) 
	// Order of site classes for iw or f(x_h|w):
	//                     fore   back     #sets
	//   site class 0:      w0     w0        10
	//   site class 1:      w1=1   w1=1       1
	//   site class 2a:     w0     w2       100
	//   site class 2b:     w1=1   w2        10

	for(unsigned int iw=0; iw < BEB_NUM_CAT; ++iw)
	{
		if(iw < BEB_N1D)										// class 0: w0 w0
		{
			omega_fg = omega_bg = mPara[2][iw];
		}
		else if(iw == BEB_N1D)									// class 1: w1 w1
		{
			omega_fg = omega_bg = 1.;
		}
		else if(iw < (BEB_N1D+1+BEB_N1D*BEB_N1D))				// class 2a: w0 w2
		{                                  
            omega_fg = mPara[2][(iw-BEB_N1D-1) / BEB_N1D];
            omega_bg = mPara[3][(iw-BEB_N1D-1) % BEB_N1D];
		}
		else													// class 2b: w1 w2
		{                                                       
            omega_fg = 1.;
            omega_bg = mPara[3][iw-BEB_N1D-1-BEB_N1D*BEB_N1D];
		}

		// Fill the matrices and compute their eigen decomposition.
#ifdef _MSC_VER
		#pragma omp parallel sections default(none) shared(omega_fg, omega_bg, kappa, scale_q_fg, scale_q_bg, q_fg, q_bg)
#else
		#pragma omp parallel sections default(shared)
#endif
		{
			#pragma omp section
			{
				scale_q_fg = q_fg.fillMatrix(omega_fg, kappa);
				q_fg.eigenQREV();
			} 
			#pragma omp section
			{
				scale_q_bg = q_bg.fillMatrix(omega_bg, kappa);
				q_bg.eigenQREV();
			}
		}

		// Fill the matrix set with the new matrices and the times computed before
		beb_set.fillMatrixSet(q_fg, q_bg, scale_q_bg, scale_q_fg, aVars);

		// Compute likelihoods
		CacheAlignedDoubleVector likelihoods(aForest.getNumSites());
		aForest.computeLikelihoods(beb_set, likelihoods, dep.getDependencies());
		memcpy(&mPriors[iw*mNumSites], &likelihoods[0], aForest.getNumSites()*sizeof(double));

		if(mVerbose >= VERBOSE_ONLY_RESULTS)
		{
			std::cerr << std::setw(3) << iw << ' ';
			std::cerr << std::fixed << std::setw(8) << std::setprecision(4) << omega_fg << ' ';
			std::cerr << std::fixed << std::setw(8) << std::setprecision(4) << omega_bg << " = " << std::setprecision(16);
			std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+0] << ' ';        // Print only the first three sites
			std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+1] << ' ';
			std::cerr << std::scientific << std::setw(14) << mPriors[iw*mNumSites+2] << std::endl;
		}
	}

	double scale = 0.;
	for(size_t site=0; site < mNumSites; ++site)
	{
		size_t k;
		double fh = mPriors[site];
		for(k=1; k < BEB_NUM_CAT; ++k)
		{
			if(mPriors[k*mNumSites+site] > fh) fh = mPriors[k*mNumSites+site];
		}
		for(k=0; k < BEB_NUM_CAT; ++k)
		{
			mPriors[k*mNumSites+site] = exp(mPriors[k*mNumSites+site]-fh);
		}
		scale += fh*aSiteMultiplicity[site];
	}

	return scale;
}

void BayesTest::getIndexTernary(double* aProbX, double* aProbY, unsigned int aTriangleIdx)
{
	unsigned int ix = static_cast<unsigned int>(sqrt(static_cast<double>(aTriangleIdx)));
	unsigned int iy = aTriangleIdx - ix*ix;

	*aProbX = (1. + floor(iy/2.)*3. + (iy%2))/(3.*BEB_N1D);
	*aProbY = (1. + (BEB_N1D - 1 - ix)*3. + (iy%2))/(3.*BEB_N1D);
}

void BayesTest::computeBEB(Forest& aForest, const std::vector<double>& aVars, const std::vector<double>& aSiteMultiplicity, size_t aFgBranch)
{
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << std::endl << "Computing BEB for fg branch " << aFgBranch << std::endl;

	// Prepare the prior list
	double scale1 = getGridParams(aForest, aVars, aSiteMultiplicity, aFgBranch);

	// Temporary here till I understand their usage
	int iw[BEB_NGRID*BEB_DIMS];
	double pclassM[BEB_NGRID*BEB_DIMS];

	// Set up im and pclassM, for each igrid and iclassM
	for(unsigned int igrid=0; igrid < BEB_NGRID; ++igrid)
	{
		unsigned int ip[BEB_DIMS];
		unsigned int it = igrid;
		for(int j=BEB_DIMS-1; j >= 0; --j)
		{
			ip[j]= it % BEB_N1D;
			it /= BEB_N1D;
		}

		// In the original code instead of BEB_DIMS there is nclassM. Both values are 4.
		for(unsigned int k=0; k < BEB_DIMS; ++k)
		{
			// Given the point on the grid ip[] and iclassM, this returns iw and pclassM, 
			//   where iw locates the correct f(x_h|w) stored in com.fhK[], and pclassM is 
			//   the proportion of the site class under the model.
			//   The n1d*n1d grid for p0-p1 is mapped onto the ternary graph for p0-p1-p2.  
			//
			//   See get_grid_para_like_AC() for order of iw or site classes.
			//
			//  Parameters are model A: (p0 p1 w0 w2)
			//
			int idx;
			switch(k)
			{
			case 0: idx = ip[2]; break;								/* class 0: w0 */
			case 1: idx = BEB_N1D; break;							/* class 1: w1 */
			case 2: idx = BEB_N1D+1+ip[2]*BEB_N1D+ip[3]; break;		/* class 2a model A: w0 & w2 */
			case 3: idx = BEB_N1D+1+BEB_N1D*BEB_N1D+ip[3]; break;	/* class 2b model A: w1 & w2 */
			}
			iw[igrid*BEB_DIMS+k] = idx;

			// Fill the volume with the probabilities for each class
			double p[3];
			getIndexTernary(&p[0], &p[1], ip[0]*BEB_N1D+ip[1]);
			p[2] = 1.0 - p[0] - p[1];

			pclassM[igrid*BEB_DIMS+k] = (k < 2) ? p[k] : p[2]*p[k-2]/(1.0-p[2]);
		}
	}

	// Calculate marginal prob of data, fX, and postpara[].  scale2 is scale.
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << std::endl << "Calculating f(X), the marginal probability of data." << std::endl;
	double fX = 1.;
	double scale2 = -1e300;
	double lnfXs[BEB_NGRID];

	for(unsigned int j=0; j < BEB_DIMS; ++j)  /* postpara[0-1] for p0p1 ignored */
		for(unsigned int k=0; k < BEB_N1D; ++k) 
			mPostPara[j][k] = 1.;

	double postp0p1[BEB_N1D*BEB_N1D];
	for(unsigned int k=0; k < BEB_N1D*BEB_N1D; k++) postp0p1[k] = 1.;

   	for(unsigned int igrid=0; igrid < BEB_NGRID; ++igrid)
	{
		unsigned int ip[BEB_DIMS];
		unsigned int it = igrid;
		for(int j=BEB_DIMS-1; j >= 0; --j)
		{
			ip[j]= it % BEB_N1D;
			it /= BEB_N1D;
		}

		lnfXs[igrid] = 0.;
		for(unsigned int site=0; site < mNumSites; ++site)
		{
			double fh = 0.;
			for(unsigned int k=0; k < BEB_DIMS; ++k)
				fh += pclassM[igrid*BEB_DIMS+k]*mPriors[iw[igrid*BEB_DIMS+k]*mNumSites+site];
			if(fh < 1e-300)
			{
				if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << "strange: f[" << site << "] = " << fh << " very small." << std::endl;
				continue;
			}
			lnfXs[igrid] += log(fh)*aSiteMultiplicity[site];
		}

		double t = lnfXs[igrid]-scale2;
		if(t > 0)
		{
			/* change scale factor scale2 */
			t = (t<200) ? exp(-t) : 0;
			fX = fX*t+1;
			for(unsigned int j=0; j < BEB_DIMS; ++j) for(unsigned int k=0; k < BEB_N1D; ++k) mPostPara[j][k] *= t;
			for(unsigned int k=0; k < BEB_N1D*BEB_N1D; ++k) postp0p1[k] *= t;

			for(unsigned int j=0; j < BEB_DIMS; ++j) mPostPara[j][ip[j]] ++;
			postp0p1[ip[0]*BEB_N1D+ip[1]] ++;
			scale2 = lnfXs[igrid];
		}
		else if(t > -200)
		{
			t = exp(t);
			fX += t;
			for(unsigned int j=0; j < BEB_DIMS; ++j) mPostPara[j][ip[j]] += t;
			postp0p1[ip[0]*BEB_N1D+ip[1]] += t;
		}
	}

	// Normalize probabilities
	for(unsigned int j=0; j<BEB_DIMS; ++j) for(unsigned int k=0; k < BEB_N1D; ++k) mPostPara[j][k] /= fX;
	for(unsigned int k=0; k<BEB_N1D*BEB_N1D; k++) postp0p1[k] /=fX;

	fX = log(fX)+scale2;

	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << "log(fX) = " << (fX+scale1-BEB_DIMS*log(BEB_N1D*1.))
		                                           << "  Scales = " << scale1 << " " << scale2 << std::endl;

	// Calculate posterior probabilities for sites.  S1 is scale factor */
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cerr << std::endl << "Calculating f(w|X), posterior probs of site classes." << std::endl;

	for(unsigned int site=0; site < mNumSites; ++site)
	{
		scale1 = -1e300;

		for(unsigned int j=0; j < BEB_DIMS; ++j) mSiteClassProb[j*mNumSites+site] = 1.;

   		for(unsigned int igrid=0; igrid < BEB_NGRID; ++igrid)
		{
			//unsigned int ip[BEB_DIMS];
			//unsigned int it = igrid;
			//for(int j=BEB_DIMS-1; j >= 0; --j)
			//{
			//	ip[j]= it % BEB_N1D;
			//	it /= BEB_N1D;
			//}

			double fh = 0.;
			double fhk[BEB_DIMS];
			for(unsigned int k=0; k < BEB_DIMS; ++k) /* duplicated calculation */
			{
				fhk[k] = pclassM[igrid*BEB_DIMS+k]*mPriors[iw[igrid*BEB_DIMS+k]*mNumSites+site];
				fh += fhk[k];
			}

			for(unsigned int iclassM=0; iclassM < BEB_DIMS; ++iclassM)
			{
				fhk[iclassM] /= fh;
				double t = log(fhk[iclassM]) + lnfXs[igrid]; /* t is log of term on grid */
				if(t > scale1 + 50)
				{  /* change scale factor scale1 */
					for(unsigned int j=0; j < BEB_DIMS; ++j)
						mSiteClassProb[j*mNumSites+site] *= exp(scale1-t);
					scale1 = t;
				}
				mSiteClassProb[iclassM*mNumSites+site] += exp(t-scale1);
			}
		}
		for(unsigned int j=0; j<BEB_DIMS; ++j) mSiteClassProb[j*mNumSites+site] *= exp(scale1-fX);

		// For debug
		if(mVerbose >= VERBOSE_ONLY_RESULTS) 
		{
			std::cerr << "Site " << std::setw(4) << site;
			for(unsigned int k=0; k < BEB_DIMS; ++k) std::cerr << std::setw(20) << std::setprecision(12) << mSiteClassProb[k*mNumSites+site];
			std::cerr << std::endl;
		}
	}
}


void BayesTest::printPositiveSelSites(size_t aFgBranch) const
{
	bool print_title = true;

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if is a type 2 site with prob > 50%
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > MIN_PROB)
		{
			// Put a title the firts time
			if(print_title)
			{
				std::cerr << "Printing positive sel sites for branch " << aFgBranch << std::endl;
				print_title = false;
			}
			
			// Set significance
			const char* sig;
			if(prob > TWO_STARS_PROB)     sig = "**";
			else if(prob > ONE_STAR_PROB) sig = "*";
			else                          sig = "";

			// Print site number and probability
			std::cerr << std::setw(6) << site << " " << std::fixed << std::setprecision(6) << prob << sig << std::endl;
		}
	}
}


void BayesTest::extractPositiveSelSites(std::vector<unsigned int>& aPositiveSelSites, std::vector<double>& aPositiveSelSitesProb) const
{
	// Prepare the output vectors
	aPositiveSelSites.clear();
	aPositiveSelSitesProb.clear();

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if it is a type 2 site with prob > minimum cutoff
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > MIN_PROB)
		{
			aPositiveSelSites.push_back(site);
			aPositiveSelSitesProb.push_back(prob);
		}
	}
}
