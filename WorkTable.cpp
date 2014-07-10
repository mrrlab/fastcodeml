
#include "WorkTable.h"

/// For each internal branch there are three jobs to be executed: H0, H1 and BEB
static const int JOBS_PER_BRANCH = 3;

WorkTable::WorkTable(const ForestGroup &aForestGroup, const CmdLine &aCmdLine)
{
    mStaticQueue = false;

    for (size_t ii = 0; ii < aForestGroup.size(); ii++)
    {
        const Forest *forest = aForestGroup.getForestPtr(ii);
        this->addForest(*forest, ii, aCmdLine);
    }
}

void
WorkTable::addForest(const Forest &aForest, size_t aForestIndex, const CmdLine &aCmdLine)
{
    size_t branch_start(0), branch_end(0);
    bool do_all = aForest.getBranchRange(aCmdLine, branch_start, branch_end);

    size_t num_internal_branches = aForest.getNumInternalBranches();

    unsigned int aHypothesisToDo = aCmdLine.mComputeHypothesis;

    if (aHypothesisToDo == 0 || aHypothesisToDo == 1)
    {
        mStaticQueue = true;
    }

    std::vector<ResultSet> result_set(num_internal_branches);
    for (size_t branch = 0; branch < num_internal_branches; branch++)
    {
        if (branch >= branch_start && branch <= branch_end)
        {
            JobType h;
            // Switch on the command line, note that if aHypothesisToDo is populated, we won't do the BEB because the queue is static.
            switch(aHypothesisToDo)
            {
                case 0: h = JOB_H0; break;
                case 1: h = JOB_H1; break;
                default: h = JOB_H1; break;
            }

            JobType jt             = h;
            size_t size_other_data = 0;
            Job new_job = Job(aForestIndex, branch, jt, size_other_data);
            mJobQueue.push(new_job);
            result_set[branch].mSkipped = false;
        }
    }
    mResults[aForestIndex]      = result_set;
    mWriteResults[aForestIndex] = WriteResults();
}

bool WorkTable::getNextJob(int* aJob)
{
    // Job structure
    // job[0] - the type of job
    // job[1] - send the branch number
    // job[2] - the size of the other data
    // job[3] - the index of the forest in forest group

    // If the queue is non-empty, pop off the next job.
    if (!mJobQueue.empty())
    {
        Job next_job(mJobQueue.front()); // Get the next job
        mJobQueue.pop();                  // Pop it off the queue

        aJob[0] = static_cast<int>(next_job.getJobType());
        aJob[1] = next_job.getBranch();
        aJob[2] = next_job.getSizeOtherData();
        aJob[3] = next_job.getForestIndex();

        return true;
    }
    else
    {
        // We are done, there are no more elements in the queue
        aJob[0] = JOB_SHUTDOWN;
        aJob[1] = 0;
        aJob[2] = 0;
        aJob[3] = 0;
        return false;
    }
}

void WorkTable::recordResultHyp(const Job &aJob, const std::vector<double> &aResults,
    int aWorker, unsigned int aVerbose)
{
    int forest = aJob.getForestIndex();
    int branch = aJob.getBranch();
    int h      = aJob.getJobType();

    // Get the relevant result set
    std::vector<ResultSet> &result_set = mResults[forest];

    // Save all results (lnl + all variables)
    // For H1 there are also the two scale values at the end of variables for BEB computation
    double lnl = aResults[aResults.size() - 1];

    result_set[branch].mLnl[h] = lnl;
    result_set[branch].mHxVariables[h].assign(aResults.begin(), aResults.end()-1);

    if(aVerbose >= VERBOSE_MPI_TRACE)
    {
        if(lnl < DBL_MAX)
            std::cout << std::fixed << std::setprecision(8) << "Lnl: " << lnl << " for H" << h << " - on forest - " << forest << " from worker " << aWorker << std::endl;
        else
            std::cout << std::fixed << "Lnl: NA for H" << h << " - on forest - " << forest << " from worker " << aWorker << std::endl;
    }

    // Save for the results file (if has been computed)
    if(lnl < DBL_MAX) mWriteResults[forest].saveLnL(static_cast<size_t>(branch), lnl, h);
}

void WorkTable::recordResultBEB(const Job &aJob, const std::vector<int> &aResults,
    int aWorker, unsigned int aVerbose)
{
    int forest = aJob.getForestIndex();
    int branch = aJob.getBranch();
    int h      = aJob.getJobType();

    // Get the relevant result set
    std::vector<ResultSet> &result_set = mResults[forest];

    // Save all results (positive selection sites and corresponding probability)
    result_set[branch].mPositiveSelSites.clear();
    result_set[branch].mPositiveSelProbs.clear();
    for(int i=0; i < aResults.size()/2; ++i)
    {
        int site = aResults[2*i+0];
        result_set[branch].mPositiveSelSites.push_back(site);
        double prob = static_cast<double>(aResults[2*i+1])/PROB_SCALING;
        result_set[branch].mPositiveSelProbs.push_back(prob);
    }

    // If there are sites under positive selection, save them for the results file
    if(aResults.size() > 0)
    {
        std::vector<unsigned int> sites;
        for(int i=0; i < aResults.size()/2; ++i)
        {
            unsigned int u = static_cast<unsigned int>(aResults[2*i+0]);
            sites.push_back(u);
        }
        mWriteResults[forest].savePositiveSelSites(static_cast<size_t>(branch), sites, result_set[branch].mPositiveSelProbs);
    }

    // Output a status message
    if(aVerbose >= VERBOSE_MPI_TRACE) std::cout << "BEB num of results: " << aResults.size()/2 << " from worker " << aWorker << std::endl;
}

void WorkTable::queueNextJob(const Job &aJob, unsigned int aVerbose, const bool aInitFromH1)
{
    // Rules
    // 1. When a H1 job completes, we queue the corresponding H0 job
    // 2. When a H0 job completes, we queue the BEB job if it is applicable.
    // 3. When the queue is static (user forced H0, or H1) we won't push anything into the queue

    if (mStaticQueue)
    {
        if(aVerbose >= VERBOSE_INFO_OUTPUT) printFinishedBranch(aJob);
        return; // Done here
    }

    // First invoke copy ctor
    Job new_job(aJob);

    // Get a ref to the results for this forest
    std::vector<ResultSet> &result_set = mResults[aJob.getForestIndex()];
    size_t branch = aJob.getBranch();

    // If we just completed an H1 job, queue the corresponding H0 job
    if (aJob.getJobType() == JOB_H1)
    {
        // Send the lnl of the corresponding H1 step or mResults[branch].mHxVariables[1]
        int other_data = 0;
        if (aInitFromH1){
            other_data = static_cast<int>(result_set[branch].mHxVariables[1].size())+1;
        }else{
            other_data = 1;
        }
        new_job.setJobType(JOB_H0);
        new_job.setSizeOtherData(other_data);
        mJobQueue.push(new_job);
    }

    // If we just finished an H0 job, check to see if we can queue the BEB job and if
    // so, queue it.
    if (aJob.getJobType() == JOB_H0)
    {
		// If the previous results do not pass the LRT or H0 has been interrupted, then skip the BEB computation
		if(result_set[branch].mLnl[0] == DBL_MAX || !BranchSiteModel::performLRT(result_set[branch].mLnl[0], result_set[branch].mLnl[1]))
        {
            if(aVerbose >= VERBOSE_INFO_OUTPUT) printFinishedBranch(aJob);
            return;
        }
        new_job.setJobType(JOB_BEB);
        new_job.setSizeOtherData(static_cast<int>(result_set[branch].mHxVariables[1].size()));

        mJobQueue.push(new_job);
    }

    // If we just finished a beb job, flag the branch as completed
    if (aJob.getJobType() == JOB_BEB)
    {
        if(aVerbose >= VERBOSE_INFO_OUTPUT) printFinishedBranch(aJob);
    }
}

std::vector<ResultSet>* WorkTable::getResultSetPtr(size_t aForestIndex)
{
    return(&mResults[aForestIndex]);
}

void
WorkTable::printResults(std::ostream& aOut, const CmdLine &aCmdLine) const
{
	// Print likelihoods (and variables)
	aOut << std::endl;

	std::map<size_t, std::vector<ResultSet> >::const_iterator cib = mResults.begin();
    std::map<size_t, std::vector<ResultSet> >::const_iterator cif = mResults.end();
    std::map<size_t, std::vector<ResultSet> >::const_iterator cit;
	for (cit = cib; cit != cif; cit++)
    {
        std::string tree_file(aCmdLine.mTreeFiles[cit->first]);
        std::string alignment_file(aCmdLine.mGeneFiles[cit->first]);
        aOut << std::endl << "Tree file :      " <<  tree_file << std::endl
            << "Alignment file : " << alignment_file << std::endl;

        const std::vector<ResultSet> &result_set =cit->second;
        printResultSet(aOut, result_set);
    }

} //printWorkTableResults

void
WorkTable::printResultSet(std::ostream &aOut, const std::vector<ResultSet> &aResultSet) const
{
	for(size_t branch=0; branch < aResultSet.size(); ++branch)
	{
		// Skip branches that were not computed
		if(aResultSet[branch].mSkipped) continue;

		aOut << "Branch: "   << std::fixed << std::setw(3) << branch;
		if(aResultSet[branch].mLnl[0] == std::numeric_limits<double>::infinity())
		{
			aOut << "  Lnl H0: " << std::setw(24) << "Inf";
		}
		else if(aResultSet[branch].mLnl[0] == DBL_MAX)
		{
			aOut << "  Lnl H0: " << std::setw(24) << "NA";
		}
		else
		{
			aOut << "  Lnl H0: " << std::setw(24) << std::setprecision(15) << aResultSet[branch].mLnl[0];
		}
		if(aResultSet[branch].mLnl[1] == std::numeric_limits<double>::infinity())
		{
			aOut << "  Lnl H1: " << std::setw(24) << "Inf";
		}
		else
		{
			aOut << "  Lnl H1: " << std::setw(24) << std::setprecision(15) << aResultSet[branch].mLnl[1];
		}
		if(aResultSet[branch].mLnl[0] == std::numeric_limits<double>::infinity() || aResultSet[branch].mLnl[1] == std::numeric_limits<double>::infinity())
		{
			aOut << "  LRT: " << std::setw(24) << "*Invalid*";
		}
		else if(aResultSet[branch].mLnl[0] < DBL_MAX)
			aOut << "  LRT: " << std::setw(24) << std::setprecision(15) << std::fixed << aResultSet[branch].mLnl[1] - aResultSet[branch].mLnl[0] << "  (threshold: " << std::setprecision(15) << std::fixed << THRESHOLD_FOR_LRT << ')';
		else
			aOut << "  LRT: < " << std::setprecision(15) << std::fixed << THRESHOLD_FOR_LRT;
		aOut << std::endl;
		aOut << std::endl;
		if(aResultSet[branch].mLnl[0] != std::numeric_limits<double>::infinity() && aResultSet[branch].mLnl[0] != DBL_MAX) printVariables(branch, 0, aOut, aResultSet[branch]);
		if(aResultSet[branch].mLnl[1] != std::numeric_limits<double>::infinity()) printVariables(branch, 1, aOut, aResultSet[branch]);
		aOut << std::endl;
	}

	// Check if there are sites under positive selection
	std::vector<size_t> branch_with_pos_selection;
	for(size_t branch=0; branch < aResultSet.size(); ++branch)
	{
		// Skip branches that were not computed
		if(aResultSet[branch].mSkipped) continue;

		if(!aResultSet[branch].mPositiveSelSites.empty())
		{
			branch_with_pos_selection.push_back(branch);
		}
	}

	// If there are print the site and corresponding probability
	if(!branch_with_pos_selection.empty())
	{
		aOut << std::endl << "Positive selection sites" << std::endl;
		std::vector<size_t>::const_iterator ib(branch_with_pos_selection.begin());
		std::vector<size_t>::const_iterator end(branch_with_pos_selection.end());
		for(; ib != end; ++ib)
		{
			const ResultSet& branch_results = aResultSet[*ib];

			// To order the sites
			std::multimap<size_t, size_t> ordered_map;
			std::vector<double> probs;
			size_t current_idx = 0;

			aOut << "Branch: "   << std::fixed << std::setw(3) << *ib << std::endl;
			for(size_t pss=0; pss < branch_results.mPositiveSelSites.size(); ++pss)
			{
				// Get probability
				double prob = branch_results.mPositiveSelProbs[pss];

				// Save site and probability to order output by site
				ordered_map.insert(std::pair<size_t, size_t>(branch_results.mPositiveSelSites[pss], current_idx));
				probs.push_back(prob);
				++current_idx;
			}

			// Print site number and probability after mapping the site number to the original value (and changing numbering so it starts from 1 and not zero)
			std::multimap<size_t, size_t>::const_iterator im(ordered_map.begin());
			std::multimap<size_t, size_t>::const_iterator endm(ordered_map.end());
			for(; im != endm; ++im)
			{
				double prob = probs[im->second];

				// Set significance
				const char* sig;
				if(prob > TWO_STARS_PROB)     sig = "**";
				else if(prob > ONE_STAR_PROB) sig = "*";
				else                          sig = "";

				aOut << std::setw(6) << im->first + 1 << ' ' << std::fixed << std::setprecision(6) << prob << sig << std::endl;
			}
		}
	}
}

void WorkTable::printFinishedBranch(const Job &aJob) const
{
	int branch = aJob.getBranch();
	int forest = aJob.getForestIndex();

    std::cout << "Forest " << std::setw(3) << forest <<  " - Branch " << std::setw(3) << branch << " completed" << std::endl;
}

void WorkTable::printVariables(size_t aBranch, unsigned int aHyp, std::ostream& aOut,
   const ResultSet &aResultSet) const
{
	aOut << "Optimized variables for H" << aHyp << " for fg branch " << aBranch << std::endl;

	// To nicely format num branch lengths per line
	static const unsigned int VARS_PER_LINE = 8;
	unsigned int count_per_line = 0;
	static const std::streamsize VARS_PRECISION = 7;
	static const std::streamsize VARS_WIDTH     = 11;

	// Write the data with uniform precision
	std::streamsize prec = aOut.precision(VARS_PRECISION);
	aOut.setf(std::ios::fixed, std::ios::floatfield);

	// Print all variables formatted to be readable
	int num_times = static_cast<int>(aResultSet.mHxVariables[aHyp].size()) - ((aHyp) ? 7 : 4); // for H1 is 5 + 2 scale factors
	double v0 = 0;
	std::vector<double>::const_iterator ix(aResultSet.mHxVariables[aHyp].begin());
	for(int k = -static_cast<int>(num_times); k < (aHyp ? 5 : 4); ++ix,++k)
	{
		switch(k)
		{
		case 0:
			if(count_per_line) aOut << std::endl;
			v0 = *ix;
			break;
		case 1:
			{
				double p[4];
#ifdef USE_ORIGINAL_PROPORTIONS
				p[0] = exp(v0);
				p[1] = exp(*ix);
				double tot = p[0] + p[1] + 1;
				p[0] /= tot;
				p[1] /= tot;
				tot = p[0] + p[1];

				p[2] = (1. - tot)*p[0]/tot;
				p[3] = (1. - tot)*p[1]/tot;
#else
				p[0] = v0*(*ix);
				p[1] = v0*(1.-(*ix));
				p[2] = (1.-v0)*(*ix);
				p[3] = (1.-v0)*(1.-(*ix));
#endif

				aOut <<   "p0:"  << std::setw(VARS_WIDTH) << p[0];
				aOut << "  p1:"  << std::setw(VARS_WIDTH) << p[1];
				aOut << "  p2a:" << std::setw(VARS_WIDTH) << p[2];
				aOut << "  p2b:" << std::setw(VARS_WIDTH) << p[3];
				aOut << std::endl;
			}
			break;
		case 2:
			aOut << "w0:" << std::setw(VARS_WIDTH) << *ix;
			break;
		case 3:
			aOut << "  k: " << std::setw(VARS_WIDTH) << *ix;
			break;
		case 4:
			aOut << "  w2: " << std::setw(VARS_WIDTH) << *ix;
			break;
		default:
			aOut << std::setw(VARS_WIDTH) << *ix;
			++count_per_line;
			if(count_per_line == VARS_PER_LINE)
			{
				count_per_line = 0;
				aOut << std::endl;
			}
			break;
		}
	}
	aOut << std::endl;
	aOut.precision(prec);
}


