
#ifndef GENES_H
#define GENES_H

#include <string>
#include <vector>
#include <map>

/// The genes of the set of species under analysis.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
///
class Genes
{
protected:
	/// Constructor
	///
	/// @param[in] aVerboseLevel The verbosity level
	///
	explicit Genes(unsigned int aVerboseLevel=0) : mVerboseLevel(aVerboseLevel)
	{
		mMapBaseToIdx['t'] = 0;
		mMapBaseToIdx['T'] = 0;
		mMapBaseToIdx['c'] = 1;
		mMapBaseToIdx['C'] = 1;
		mMapBaseToIdx['a'] = 2;
		mMapBaseToIdx['A'] = 2;
		mMapBaseToIdx['g'] = 3;
		mMapBaseToIdx['G'] = 3;
	}

	/// Destructor
	///
	virtual ~Genes();

public:
	/// Load the gene file.
	///
	/// @param[in] aFilename The filename containing the genes under analysis
	///
	/// @exception FastCodeMLFatalNoMsg If cannot open file and other problems
	///
	void readFile(const char* aFilename);

	/// Return the number of valid sites loaded
	///
	/// @return The number of sites
	///
	size_t getNumSites(void) const {return mSiteMultiplicity.size();}

	/// Return the site multiplicity
	///
	/// @return Pointer to the site multiplicity array
	///
	const std::vector<unsigned int>& getSiteMultiplicity(void) const {return mSiteMultiplicity;}

	/// Return the codon index for the one identified by the specie label and the site.
	///
	/// @param[in] aSpecie The specie label.
	/// @param[in] aSite The site index.
	///
	/// @return The codon index or -1 in case of error or invalid codon at the given position.
	///
	int getCodonIdx(std::string aSpecie, size_t aSite) const;

	/// Return the list of species present in the gene file
	///
	/// @param[out] aSpeciesList Will be filled with the list of species
	///
	void getSpecies(std::vector<std::string>& aSpeciesList) const {aSpeciesList = mDnaSpecies;}

	/// Check coherence between tree and genes.
	///
	/// @param[in] aNames The phylogenetic tree species names
	///
	/// @exception FastCodeMLFatal Throw exception if the species do not match
	///
	void checkNameCoherence(const std::vector<std::string>& aNames) const;


protected:
	/// Test if the three letters of the argument represent a valid codon 
	///
	/// @param[in] aCodon String of three letters TCAG repesenting the codon (not null terminated)
	/// @return True if codon is a valid codon and not a stop codon
	///
	bool validCodon(const char* aCodon) const {return idxCodon(aCodon) >= 0;}

	/// Return the index of the given codon
	///
	/// @param[in] aCodon String of three letters TCAG repesenting the codon (not null terminated)
	/// @return Index in the range 0 - 60 (or -1 if the codon is invalid)
	///
	int idxCodon(const char* aCodon) const;

private:
	/// Load the gene file.
	/// This routine is redefined in every derived class to load a specific format.
	///
	/// @param[in] aFilename The filename containing the genes under analysis
	/// @param[out] aSpecies The list of species read
	/// @param[out] aSequences Array of genes, one per specie
	///
	/// @exception FastCodeMLFatalNoMsg On various error conditions
	///
	virtual void loadData(const char* aFilename, std::vector<std::string>& aSpecies, std::vector<std::string>& aSequences) =0;

protected:
	unsigned int						mVerboseLevel;			///< The verbosity level as set in the constructor

private:
	std::vector<std::string>			mDnaSpecies;			///< The list of species labels
	std::vector<std::string>			mDnaGene;				///< The gene DNA basis strings
	std::vector<unsigned int>			mSiteMultiplicity;		///< Site multiplicity (sites with multiplicity of zero has been removed from the site list)
	std::vector<unsigned int>			mMapSiteToDnaGene;		///< Map the site number to the position in mDnaGene
	std::map<std::string, unsigned int> mMapSpecieToDnaGene;	///< Map specie name to position in the gene list mDnaGene
	std::map<char, int>					mMapBaseToIdx;			///< Map DNA base letter (TCAG) to number 0 to 3
};

#endif

