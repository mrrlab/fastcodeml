
#ifndef BAYESTEST_H
#define BAYESTEST_H

#include <vector>
#include <cstdlib>
#include "BranchSiteModel.h"

/// The minimum value for class 2 sites probability to be a positive selection site.
///
const static double MIN_PROB       = 0.50;
const static double ONE_STAR_PROB  = 0.95;
const static double TWO_STARS_PROB = 0.99;

/// Helper class to compute BEB_N1D^BEB_DIMS at compile time (that is Y^N)
template<unsigned int Y, unsigned int N>
class Pow
{
public:
	static const int value = Y * Pow<Y, N-1>::value;
};
template<unsigned int Y>
class Pow<Y, 1>
{
public:
	static const int value = Y;
};

/// Tests to find the sites under positive selection.
///
///  @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///  @date 2010-12-22 (initial version)
///  @version 1.0
///
class BayesTest
{
public:
	/// Constructor.
	///
	/// @param[in] aNumSites Number of sites
	/// @param[in] aVerbose The verbosity level
	///
	explicit BayesTest(size_t aNumSites, unsigned int aVerbose=0)
						: mSiteClassProb(BEB_DIMS*aNumSites), mNumSites(aNumSites), mVerbose(aVerbose), mPriors(aNumSites*BEB_NUM_CAT) {}

	/// Destructor.
	///
	~BayesTest() {}

	/// Bayes Empirical Bayes (BEB) test.
	///
	/// @param[in] aModel The last computed H1 test
	///
	void computeBEB(Forest& aForest, const std::vector<double>& aVars, const std::vector<double>& aSiteMultiplicity, size_t aFgBranch);
	
	/// Print the sites under positive selection.
	///
	/// @param[in] aFgBranch Identifier of the branch marked as foreground branch
	///
	void printPositiveSelSites(size_t aFgBranch) const;

	/// Extract the sites under positive selection and the corresponding probabilities.
	///
	/// @param[out] aPositiveSelSites Vector of sites under positive selection
	/// @param[out] aPositiveSelSitesProb Corresponding probabilities
	///
	void extractPositiveSelSites(std::vector<unsigned int>& aPositiveSelSites, std::vector<double>& aPositiveSelSitesProb) const;

private:
	/// This sets up the grid (mPara[][]) according to the priors.  
	/// It calculates the probability of data at each site given w: f(f_h|w).  
	/// This is calculated using the branch model (NSsites = 0 model = 2), with 
	/// BayesEB=2 used to force the use of the correct scale factors in GetPMatBranch().
	///
	///@verbatim
	/// Order of site classes for iw or f(x_h|w):
	///                     fore   back     num.sets
	/// Branchsite A (121 sets)
	///   site class 0:      w0     w0        10
	///   site class 1:      w1=1   w1=1       1
	///   site class 2a:     w0     w2       100
	///   site class 2b:     w1=1   w2        10
	///@endverbatim
	///
	double getGridParams(Forest& aForest, const std::vector<double>& aVars, const std::vector<double>& aSiteMultiplicity, size_t aFgBranch);

	///    This gives the indices (ix, iy) and the coordinates (aProbX, aProbY, 1-aProbX-aProbY) for 
	///    the aTriangleIdx-th triangle, with aTriangleIdx from 0, 1, ..., BEB_N1D*BEB_N1D-1.  
	///    The ternary graph (0-1 on each axis) is partitioned into BEB_N1D*BEB_N1D equal-sized triangles.  
	///    In the first row (ix=0), there is one triangle (iy=0);
	///    In the second row (ix=1), there are 3 triangles (iy=0,1,2);
	///    In the i-th row (ix=i), there are 2*i+1 triangles (iy=0,1,...,2*i).
	///
	///    aProbX rises when ix goes up, but aProbY decreases when iy increases.  (aProbX, aProbY) is the 
	///    centroid in the ij-th small triangle.
	///    
	///    aProbX and aProbY each takes on 2*BEB_N1D-1 possible values.
	///
	void getIndexTernary(double* aProbX, double* aProbY, unsigned int aTriangleIdx);


private:
	const static unsigned int BEB_N1D = 10;												///< Number of intervals for w0 and w2
	const static unsigned int BEB_DIMS = 4;												///< Number of codon classes (0, 1, 2a, 2b)
	const static unsigned int BEB_NUM_CAT = BEB_N1D + 1 + BEB_N1D*BEB_N1D + BEB_N1D;	///< Total number of categories for w0 and w2 (it is com.ncatG in codeml.c)
	//const static unsigned int BEB_NGRID = BEB_N1D*BEB_N1D*BEB_N1D*BEB_N1D;			///< Number of points in the grid used to evaluate the integral. It is BEB_N1D^BEB_DIMS
	const static unsigned int BEB_NGRID = Pow<BEB_N1D, BEB_DIMS>::value;				///< Number of points in the grid used to evaluate the integral. It is BEB_N1D^BEB_DIMS

private:
	std::vector<double> mSiteClassProb;					///< Probability of a site to pertain to a given class (one row per class (4 classes), one column per site).
	size_t				mNumSites;						///< Number of sites.
	unsigned int		mVerbose;						///< If greather than zero prints more info
	std::vector<double>	mPriors;						///< Computed priors (each points to a list, one for each site)
};

#endif

