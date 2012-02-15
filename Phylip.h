
#ifndef PHYLIP_H
#define PHYLIP_H

#include "Genes.h"


/// The genes of the set of species under analysis.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-02-15 (initial version)
///     @version 1.0
///
///
class Phylip : public Genes
{
public:
	/// Constructor
	///
	/// @param[in] aVerboseLevel The verbosity level
	///
	explicit Phylip(unsigned int aVerboseLevel=0) : Genes(aVerboseLevel) {}

	/// Load the gene file. On error throw exceptions.
	///
	/// @param[in] aFilename The filename containing the genes under analysis
	///
	virtual void loadGenesFile(const char* aFilename);
};

#endif

