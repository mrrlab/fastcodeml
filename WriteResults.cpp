
#include <iostream>
#include <iomanip>
#include <fstream>
#include "WriteResults.h"


void WriteResults::outputResults(bool aOutputToStdout)
{
	// If should not output to screen and no file set, then do nothing
	if(!aOutputToStdout && !mFilename) return;

	// Redirect cout to the file if requested
	std::ofstream out;
	std::streambuf* backup;
	if(!aOutputToStdout)
	{
		out.open(mFilename, std::ios_base::trunc | std::ios_base::out);
		if(out.good())
		{
			backup = std::cout.rdbuf();
			std::cout.rdbuf(out.rdbuf());
		}
		else
		{
			std::cerr << "Cannot create results file <" << mFilename << "> Sending to screen" << std::endl;
			aOutputToStdout = true;
		}
	}

	// Range of branches
	std::map<size_t, double>::const_iterator im = mLnL[0].begin();
	size_t min_branch = im->first;
	size_t max_branch = min_branch;
	for(; im != mLnL[0].end(); ++im)
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

	// Write the log-likelihood values (if a value is not present, write NA)
	size_t branch;
	for(branch = min_branch; branch <= max_branch; ++branch)
	{
		std::cout << "Branch: " << std::setw(4) << branch << "  LnL0: ";

		// Prints LnL for H0 if present
		im = mLnL[0].find(branch);
		if(im == mLnL[0].end())
		{
			std::cout << std::setw(22) << "NA";
		}
		else
		{
			std::cout << std::setw(22) << std::setprecision(15) << std::fixed << im->second;
		}
		std::cout << "  LnL1: ";

		// Prints LnL for H1 if present
		im = mLnL[1].find(branch);
		if(im == mLnL[1].end())
		{
			std::cout << std::setw(22) << "NA";
		}
		else
		{
			std::cout << std::setw(22) << std::setprecision(15) << std::fixed << im->second;
		}
		std::cout << std::endl;
	}

	// Write the positive selection sites
	for(branch = min_branch; branch <= max_branch; ++branch)
	{
		std::map<size_t, std::pair<std::vector<unsigned int>, std::vector<double> > >::const_iterator ipss;
		ipss = mPositiveSelSites.find(branch);
		if(ipss != mPositiveSelSites.end())
		{
			const std::vector<unsigned int>& site = ipss->second.first;
			const std::vector<double>& prob       = ipss->second.second;

			size_t ns = site.size();
			for(size_t s=0; s < ns; ++s)
			{
				std::cout << "PositiveSelectionSite for branch: " << std::setw(4) << branch;
				std::cout << "  Site: " << std::setw(6) << site[s] << "  Prob: " << std::setw(9) << std::setprecision(6) << std::fixed << prob[s] << std::endl;
			}
		}
	}

	// Undo the redirect if setup
	if(!aOutputToStdout)
	{
		std::cout.rdbuf(backup);
		out.close();
	}
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

