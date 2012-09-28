
#ifndef TREEANDSETSDEPENDENCIES_H
#define TREEANDSETSDEPENDENCIES_H

#include "Forest.h"

/// Create the lists of trees and sets that can be computed concurrently
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-21 (initial version)
///     @version 1.0
///
///
class TreeAndSetsDependencies
{
public:
	/// Constructor.
	///
	/// @param[in] aVerbose The verbosity level
	///
	explicit TreeAndSetsDependencies(unsigned int aVerbose=0) {}

	/// Destructor.
	///
	~TreeAndSetsDependencies() {}

			// Prepare the dependency list for H0: 3 codon classes; H1: 4 codon classes

	void computeDependencies(const Forest& aForest, unsigned int aNumSets, unsigned int aNumThreads);

private:
	typedef std::vector< std::vector<unsigned int> >
							ListDependencies;		///< List (each list depends on the previous) of list (sites to be executed in parallel) of pairs (site, site class) stored as site*4+site_class 
	ListDependencies mDependenciesClassesAndTrees;	///< The groups of dependencies between trees
};


#endif
