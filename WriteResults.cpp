
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <sstream>
#include "WriteResults.h"
#include <string>

void
WriteResults::outputResultsToFile(
    const char *aFileName,
    const std::string &aResults)
{
	// If no file set, then do nothing
	if(aFileName==0) return;

    if (aResults.length() == 0) return;

	// Open the output file
	// Open the output file
	std::ofstream f_out(aFileName, std::ios_base::trunc | std::ios_base::out);
	if(!f_out.good())
	{
		std::cout << "Cannot create results file <" << aFileName << ">" << std::endl;
		return;
	}

	f_out << aResults;

	f_out.close();
}

void WriteResults::outputResults(void)
{
	std::string results = outputResultsToString();
    WriteResults::outputResultsToFile(mFilename, results);
}

std::string
WriteResults::outputResultsToString() const
{
	std::ostringstream out;

	// Range of branches to be output (for H0 and H1)
	std::map<size_t, double>::const_iterator im;
	std::map<size_t, std::string>::const_iterator ims;

	size_t min_branch = std::numeric_limits<size_t>::max();
	size_t max_branch = 0;
	for(im = mLnL[0].begin(); im != mLnL[0].end(); ++im)
	{
		size_t v = im->first;
		if(v < min_branch) min_branch = v;
		if(v > max_branch) max_branch = v;
	}
	for(im = mLnL[1].begin(); im != mLnL[1].end(); ++im)
	{
		size_t v = im->first;
		if(v < min_branch) min_branch = v;
		if(v > max_branch) max_branch = v;
	}

	// No entries, so do nothing
	if(min_branch == std::numeric_limits<size_t>::max()) return std::string();

	// Write the log-likelihood values (if a value is not present, write NA)
	for(size_t branch = min_branch; branch <= max_branch; ++branch)
	{
		out << "Branch: " << std::setw(4) << branch<<std::endl<< std::endl<< "  LnL0: ";

		// Prints LnL for H0 if present
		im = mLnL[0].find(branch);
		if(im == mLnL[0].end())
		{
			out << std::setw(22) << "NA"<< std::endl<<std::endl;
		}
		else
		{
			out << std::setw(22) << std::setprecision(15) << std::fixed << im->second;

            ims = mParamStr[0].find(branch);

            if (ims != mParamStr[0].end())
			out << std::endl<< std::fixed <<"Branch lengths:" << ims->second << std::endl;
		}
		out << "  LnL1: ";

		// Prints LnL for H1 if present
		im = mLnL[1].find(branch);
		if(im == mLnL[1].end())
		{
			out << std::setw(22) << "NA";
		}
		else
		{
			out << std::setw(22) << std::setprecision(15) << std::fixed << im->second;

            ims = mParamStr[1].find(branch);

            if (ims != mParamStr[1].end())
            out << std::endl<< std::fixed  <<"Branch lengths:" << ims->second << std::endl;
		}
		out << std::endl;
	}

	// Write the positive selection sites (adding one to the site number because they start from 1 and not zero)
	for(size_t branch = min_branch; branch <= max_branch; ++branch)
	{
		std::map<size_t, std::pair<std::vector<unsigned int>, std::vector<double> > >::const_iterator ipss;
		ipss = mPositiveSelSites.find(branch);
		if(ipss != mPositiveSelSites.end())
		{
			const std::vector<unsigned int>& site = ipss->second.first;
			const std::vector<double>& prob       = ipss->second.second;

			const std::vector<size_t>& idx = orderSites(site);

			size_t ns = site.size();
			for(size_t s=0; s < ns; ++s)
			{
				size_t i = idx[s];
				out << "PositiveSelectionSite for branch: " << std::setw(4) << branch;
				out << "  Site: " << std::setw(6) << site[i] + 1 << "  Prob: " << std::setw(9) << std::setprecision(6) << std::fixed << prob[i] << std::endl;
			}
		}
	}

	std::string results(out.str());
	return (results);
} //outputResultsToString

const std::vector<size_t>& WriteResults::orderSites(const std::vector<unsigned int>& aSites) const
{
	// Fill the map with pairs(site, its index)
	std::multimap<unsigned int, size_t> ordered_map;
	size_t end(aSites.size());
	for(size_t i=0; i < end; ++i) ordered_map.insert(std::pair<unsigned int, size_t>(aSites[i], i));

	// Take the list of indices after ordering the list of sites
	mSiteOrder.clear();
	std::multimap<unsigned int, size_t>::const_iterator im(ordered_map.begin());
	std::multimap<unsigned int, size_t>::const_iterator endm(ordered_map.end());
	for(; im != endm; ++im)  mSiteOrder.push_back(im->second);

	return mSiteOrder;
}

void WriteResults::saveLnL(size_t aFgBranch, double aLnL, unsigned int aHypothesis)
{
	// Sanity check
	if(aHypothesis > 1) return;

	// Save the likelihood for later printing
	mLnL[aHypothesis][aFgBranch] = aLnL;
}

void WriteResults::savePositiveSelSites(size_t aFgBranch, const std::vector<unsigned int>& aPositiveSelSites, const std::vector<double>& aPositiveSelSitesProb)
{
	// Save the positive selection sites and corresponding probabilities for later printing
	mPositiveSelSites[aFgBranch] = std::make_pair(aPositiveSelSites, aPositiveSelSitesProb);
}

void WriteResults::saveParameters(size_t aFgBranch, std::string& aParamStr, unsigned int aHypothesis)
{
	// Sanity check
	if(aHypothesis > 1) return;

	// Save the likelihood for later printing
	mParamStr[aHypothesis][aFgBranch] = aParamStr;
}

