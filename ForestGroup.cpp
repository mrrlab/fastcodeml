#include "ForestGroup.h"
#include <iomanip>
#include "WriteResults.h"
#include "BranchSiteModel.h"
#include "BayesTest.h"

ForestGroup::~ForestGroup()
{
    // Dtor is going to delete all the forests, because we own them.
    for (size_t ii = 0; ii < mForests.size(); ii++)
    {
        if (mForests[ii] != NULL)
        {
            delete(mForests[ii]);
            mForests[ii] = NULL;
        }
        if (mCodonFrequencies[ii] != NULL)
        {
            delete(mCodonFrequencies[ii]);
            mCodonFrequencies[ii] = NULL;
        }
        if (mNullHypModels[ii] != NULL)
        {
            delete(mNullHypModels[ii]);
            mNullHypModels[ii]=NULL;
        }
        if (mAltHypModels[ii] != NULL)
        {
            delete(mAltHypModels[ii]);
            mAltHypModels[ii]=NULL;
        }
        if (mBayesTests[ii] != NULL)
        {
            delete(mBayesTests[ii]);
            mBayesTests[ii]=NULL;
        }
    }
    mForests.clear();
    mNullHypModels.clear();
    mAltHypModels.clear();
    mBayesTests.clear();
}

void
ForestGroup::initForests(const CmdLine &aCmdLine)
{
    // We add all the files that CmdLine picked up for us in the
    // given to us.
    for (size_t ii = 0; ii < aCmdLine.mTreeFiles.size(); ii++)
    {
        this->addForest(
            aCmdLine,
            aCmdLine.mTreeFiles[ii].c_str(),
            aCmdLine.mGeneFiles[ii].c_str());
    }

    if (mForests.size() == 0)
    {
        throw FastCodeMLFatal("Halting - no forests.");
    }
}

void
ForestGroup::addForest(
    const CmdLine &cmd,
    const char *treeFile,
    const char *geneFile)
{
    if (cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
    {
            std::cout << std::setw(15) << "Initializing structures for tree at : " << treeFile << std::endl;
    }
	// Create the forest
	Forest *forest = new Forest(cmd.mVerboseLevel);

    // The codon frequencies object  - loadTreeAndGenes will allocate codon frequencies onto heap, we own the structure after the call.
    CodonFrequencies *codon_frequencies = NULL;

	// Load the multiple sequence alignment (MSA)
	Phylip msa(cmd.mVerboseLevel);
	msa.readFile(geneFile, cmd.mCleanData);

	// Load the phylogenetic tree
	Newick tree(cmd.mVerboseLevel);
	tree.readFile(treeFile);

	// Check coherence between the two files
	msa.checkNameCoherence(tree.getSpecies());

	// Check root
	tree.checkRootBranches();

	// If times from file then check for null branch lengths for any leaf
	if(cmd.mBranchLengthsFromFile)
	{
		int zero_on_leaf_cnt = 0;
		int zero_on_int_cnt  = 0;
		tree.countNullBranchLengths(zero_on_leaf_cnt, zero_on_int_cnt);

		if(zero_on_leaf_cnt > 0 || zero_on_int_cnt > 0)
		{
			std::cout << "Found null or missing branch length in tree file. On leaves: " << zero_on_leaf_cnt << "  on internal branches: " << zero_on_int_cnt << std::endl;
		}

		if(zero_on_leaf_cnt > 0)
		{
			throw FastCodeMLFatal("Null or missing branch length in tree file");
		}
	}

	// Print the tree with the numbering of internal branches
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) tree.printTreeAnnotated(std::cout);

	// Load the forest
	codon_frequencies = forest->loadTreeAndGenes(tree, msa,
        cmd.mIgnoreFreq ? CodonFrequencies::CODON_FREQ_MODEL_UNIF : CodonFrequencies::CODON_FREQ_MODEL_F3X4);
    mCodonFrequencies.push_back(codon_frequencies); // now is owned by this structure.

	// Reduce the forest merging common subtrees. Add also more reduction, then clean the no more useful data.
	if(!cmd.mDoNotReduceForest)
	{
		//bool sts = forest->reduceSubtrees(cmd.mNumReductionBlocks);
		forest->reduceSubtrees();

#ifndef NEW_LIKELIHOOD
		forest->addAggressiveReduction();
#endif
		forest->cleanReductionWorkingData();
#ifdef NEW_LIKELIHOOD
		forest->prepareNewReduction();
#endif
	}
#ifdef NEW_LIKELIHOOD
	else
	{
		forest->prepareNewReductionNoReuse();
	}
#endif

#ifdef NON_RECURSIVE_VISIT
	// Prepare the pointers to visit the trees without recursion
	forest->prepareNonRecursiveVisit();
#endif

	// Subdivide the trees in groups based on dependencies
	//forest->prepareDependencies(cmd.mForceSerial || cmd.mDoNotReduceForest);

#ifdef USE_DAG
	// Load the forest into a DAG
	forest->loadForestIntoDAG(Nt);
#endif

	// Print few statistics
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) std::cout << *forest;

    mForests.push_back(forest);

    // Initialize the two hypothesis
	BranchSiteModelNullHyp *h0 = new BranchSiteModelNullHyp(*forest, mCodonFrequencies.back(), cmd);
	BranchSiteModelAltHyp  *h1 = new BranchSiteModelAltHyp(*forest, mCodonFrequencies.back(), cmd);

	// Initialize the BEB (no verbose at all)
	BayesTest *beb = new BayesTest(*forest, mCodonFrequencies.back(), 0, cmd.mDoNotReduceForest);

	mNullHypModels.push_back(h0);
	mAltHypModels.push_back(h1);
	mBayesTests.push_back(beb);
}

std::vector<Forest*>
ForestGroup::getForests()
{
    return(mForests);
}

std::string
ForestGroup::solveForest(Forest &aForest, size_t aForestIndex, const CmdLine &aCmdLine)
{
    std::ostringstream os;

    os << std::endl << "Solving tree No. " << mCurrentForestIndex  << std::endl
        << "Tree file     - " << aCmdLine.mTreeFiles[aForestIndex] << std::endl
        << "Alignment file- " << aCmdLine.mGeneFiles[aForestIndex] << std::endl;

    mCurrentForestIndex++;

    if (aCmdLine.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
    {
        std::cout << os.str() << std::endl;
    }

    // Initialize the output results file (if the argument is null, no file is created)
	WriteResults output_results;

	// Compute the range of branches to mark as foreground
	size_t branch_start, branch_end;
	aForest.getBranchRange(aCmdLine, branch_start, branch_end);

	// Initialize the models
	BranchSiteModelNullHyp &h0  = *(getNullHypothesisTest(aForestIndex));
	BranchSiteModelAltHyp  &h1  = *(getAltHypothesisTest(aForestIndex));

	// Initialize the test
	BayesTest             &beb  = *(getBayesTest(aForestIndex));

	// For all requested internal branches
	for(size_t fg_branch=branch_start; fg_branch <= branch_end; ++fg_branch)
	{
		if(aCmdLine.mVerboseLevel >= VERBOSE_ONLY_RESULTS) std::cout << std::endl << "Doing branch " << fg_branch << std::endl;

		// Compute the alternate model maximum loglikelihood
		double lnl1 = 0.;
		if(aCmdLine.mComputeHypothesis != 0)
		{
			if(aCmdLine.mInitFromParams)			h1.initFromParams();
			if(aCmdLine.mBranchLengthsFromFile)	h1.initFromTree();

			lnl1 = h1(fg_branch);

			// Save the value for formatted output
			output_results.saveLnL(fg_branch, lnl1, 1);
		}

		// Compute the null model maximum loglikelihood
		double lnl0 = 0.;
		if(aCmdLine.mComputeHypothesis != 1)
		{
			if(aCmdLine.mInitH0fromH1)				h0.initFromResult(h1.getVariables());
			else
			{
				if(aCmdLine.mInitFromParams)			h0.initFromParams();
				if(aCmdLine.mBranchLengthsFromFile)	h0.initFromTree();
			}

			lnl0 = h0(fg_branch, aCmdLine.mStopIfNotLRT && aCmdLine.mComputeHypothesis != 0, lnl1-THRESHOLD_FOR_LRT);

			// Save the value for formatted output (only if has not be forced to stop)
			if(lnl0 < DBL_MAX) output_results.saveLnL(fg_branch, lnl0, 0);
		}

		if(aCmdLine.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
		{
			std::cout << std::endl;
			if(aCmdLine.mComputeHypothesis != 1)
			{
				std::cout << "LnL0: ";
				if(lnl0 == std::numeric_limits<double>::infinity())
					std::cout << "**Invalid result**";
				else if(lnl0 < DBL_MAX)
					std::cout << std::setprecision(15) << std::fixed << lnl0;
				else
					std::cout << "(Doesn't pass LRT, skipping)";
				std::cout << " Function calls: " << h0.getNumEvaluations() << "   ";
				std::cout << std::endl << std::endl;
				if(lnl0 != std::numeric_limits<double>::infinity())
				{
                    std::string s0 = h0.printFinalVars(std::cout);
                    //std::cout<<"EDW0: "<< s0 <<std::endl;
                    output_results.saveParameters(fg_branch, s0, 0);
				}
				std::cout << std::endl;
			}
			if(aCmdLine.mComputeHypothesis != 0)
			{
				std::cout << "LnL1: ";
				if(lnl1 == std::numeric_limits<double>::infinity())
					std::cout << "**Invalid result**";
				else
					std::cout << std::setprecision(15) << std::fixed << lnl1;
				std::cout << " Function calls: " << h1.getNumEvaluations();
				std::cout << std::endl << std::endl;
				if(lnl1 != std::numeric_limits<double>::infinity())
				{
				    std::string s1= h1.printFinalVars(std::cout);
				    //std::cout<<"EDW1: "<< s1 <<std::endl;
				    output_results.saveParameters(fg_branch, s1, 1);
				}
				std::cout << std::endl;
			}
			if(aCmdLine.mComputeHypothesis > 1)
			{
				if(lnl0 == std::numeric_limits<double>::infinity() || lnl1 == std::numeric_limits<double>::infinity())
					std::cout << "LRT: **Invalid result**";
				else if(lnl0 < DBL_MAX)
					std::cout << "LRT: " << std::setprecision(15) << std::fixed << lnl1 - lnl0 << "  (threshold: " << std::setprecision(15) << std::fixed << THRESHOLD_FOR_LRT << ')';
				else
					std::cout << "LRT: < " << std::setprecision(15) << std::fixed << THRESHOLD_FOR_LRT;
				std::cout << std::endl;
			}
		}

		// If requested set the time in the forest and export to a graph visualization tool
		// This switch is not honoured when in multiple mode.
		if(aCmdLine.mGraphFile && !aCmdLine.mMultipleMode)
		{
			switch(aCmdLine.mExportComputedTimes)
			{
			case 0:
				h0.saveComputedTimes();
				break;

			case 1:
				h1.saveComputedTimes();
				break;

			default:
				break;
			}

			// Use the forest export class
			ForestExport fe(aForest);
			fe.exportForest(aCmdLine.mGraphFile, fg_branch);
		}

		// If the two hypothesis are computed, H0 has not been stopped and the run passes the LRT, then compute the BEB
		if(aCmdLine.mComputeHypothesis > 1 && lnl0 < DBL_MAX && BranchSiteModel::performLRT(lnl0, lnl1))
		{
			// Get the scale values from the latest optimized h1.
			std::vector<double> scales(2);
			h1.getScales(scales);

			// Run the BEB test
			beb.computeBEB(h1.getVariables(), fg_branch, scales);

			// Output the sites under positive selection (if any)
			if(aCmdLine.mVerboseLevel >= VERBOSE_ONLY_RESULTS) beb.printPositiveSelSites(fg_branch);

			// Get the sites under positive selection for printing in the results file (if defined)
            std::vector<unsigned int> positive_sel_sites;
            std::vector<double> positive_sel_sites_prob;
            beb.extractPositiveSelSites(positive_sel_sites, positive_sel_sites_prob);
            output_results.savePositiveSelSites(fg_branch, positive_sel_sites, positive_sel_sites_prob);
		}
	}

	// Output the results
	std::string results(output_results.outputResultsToString());
	os << results << std::endl;
	return(os.str());
}

size_t
ForestGroup::getTotalNumInternalBranches(const CmdLine &aCmdLine) const
{
    size_t total_num_internal_branches(0);
    for (size_t ii = 0; ii < mForests.size(); ii++)
    {
        size_t branch_start, branch_end;
        mForests[ii]->getBranchRange(aCmdLine, branch_start, branch_end);

        // The + 1 is because we are inclusive on both ends.
        total_num_internal_branches += (branch_end - branch_start + 1);
    }
    return(total_num_internal_branches);
}
