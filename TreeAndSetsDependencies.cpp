
#include <boost/dynamic_bitset.hpp>
#include "TreeAndSetsDependencies.h"

void TreeAndSetsDependencies::computeDependencies(const Forest& aForest, unsigned int aNumSets, unsigned int aNumThreads)
{
	// Take values from forest
	size_t num_sites = aForest.getNumSites();
	std::vector< std::vector<unsigned int> > tree_dependencies = aForest.getTreeDependencies();

	// Collect the class dependencies
	std::vector< std::vector<unsigned int> > tree_groups_dependencies;

	// If no dependencies
	if(aNumThreads < 2)
	{
		std::vector<unsigned int> v(num_sites);

		for(size_t k=0; k < num_sites; ++k) v[k] = static_cast<unsigned int>(num_sites-k-1); // Remember: prior (could) point to subsequent

		tree_groups_dependencies.push_back(v);
	}
	else
	{
		size_t i, j;

		// Prepare the search of dependencies
		boost::dynamic_bitset<> done(num_sites);	// The sites that has dependencies satisfied in the previous level
		boost::dynamic_bitset<> prev;				// Dependencies till the previous level
		std::vector<unsigned int> v;				// Temporary list of sites

		// Mark trees without dependencies
		// mTreeDependencies[tj] can be done after: t1 t2 t3
		for(i=0; i < num_sites; ++i)
		{
			if(tree_dependencies[i].empty())
			{
				done.set(i);
				v.push_back(static_cast<unsigned int>(i));
			}
		}

		// Prepare the dependency list
		tree_groups_dependencies.push_back(v);
		prev = done;

		// Start to find trees with one, two, ... dependencies
		for(unsigned int numdep=1;; ++numdep)
		{
			v.clear();
			bool all_done = true;
			for(i=0; i < num_sites; ++i)
			{
				// If tree i has been already processed skip it
				if(prev[i]) continue;
				all_done = false;

				size_t nc = tree_dependencies[i].size();
				bool all = true;
				for(j=0; j < nc; ++j) if(!prev[tree_dependencies[i][j]]) {all = false; break;}
				if(all)
				{
					v.push_back(static_cast<unsigned int>(i));
					done.set(i);
				}
			}
			if(all_done) break;
			tree_groups_dependencies.push_back(v);
			prev = done;
		}
	}

	// Number of dependencies classes
	size_t nc = tree_groups_dependencies.size();
	
	// One dependency classes
	std::vector<unsigned int> one_class;

	// Transform the list multiplying the entries by the number of codon classes
	mDependenciesClassesAndTrees.clear();
	for(size_t i=0; i < nc; ++i)
	{
		// Prepare the dependency classe
		one_class.clear();

		// Number of trees in the class
		const size_t nt = tree_groups_dependencies[i].size();
		for(unsigned int set=0; set < aNumSets; ++set)
		{
			for(size_t j=0; j < nt; ++j)
			{
				one_class.push_back(Forest::makePair(tree_groups_dependencies[i][j], set));
			}
		}
		mDependenciesClassesAndTrees.push_back(one_class);
	}
}
