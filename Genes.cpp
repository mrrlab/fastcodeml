
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <limits>
#include <cctype>
#include "Genes.h"
#include "Exceptions.h"
#include "VerbosityLevels.h"


Genes::~Genes()
{
	mDnaSpecies.clear();		
	mDnaGene.clear();			
	mSiteMultiplicity.clear();	
	mMapSiteToDnaGene.clear();	
	mMapSpecieToDnaGene.clear();
	mSitesMappingToOriginal.clear();
	mMapCodonToPosition.clear();
	mCurrentPositions.clear();
}


bool Genes::validCodon(const char* aCodon, bool aRemoveAmbiguous) const
{
	const std::vector<int>& pos = getPositions(aCodon);

	if(pos.empty()) return false;		// Invalid codon
	if(pos.size() == 1) return true;	// Valid, non ambiguous codon
	return !aRemoveAmbiguous;			// Ambiguous codon
}


long long Genes::getCodonIdx(std::string aSpecie, size_t aSite) const
{
	// Find the specie
	const unsigned int idx = mMapSpecieToDnaGene.find(aSpecie)->second;

	// Access its gene
	const char* gene = mDnaGene[idx].c_str();

	// Convert the site to the position on the gene
	unsigned int position_on_gene = mMapSiteToDnaGene[aSite];

	// Get the list of positions for the given codon
	const std::vector<int>& pos = getPositions(gene+3*position_on_gene);

	// Update last decoded set of positions
	mCurrentPositions.assign(pos.begin(), pos.end());

	// Return -1 if invalid, the position or a negative code summarizing all the positions
	if(pos.empty()) return -1;
	if(pos.size() == 1) return pos[0];
	long long code = 0;
	for(size_t i=0; i < pos.size(); ++i) code |= (static_cast<long long>(1) << pos[i]);

	return -code;
}

void Genes::setLeaveProb(double* aLeaveProbVect, double aProb) const
{
	if(mCurrentPositions.empty())
	{
		throw FastCodeMLFatal("Invalid codon found in setLeaveProb.");
	}
	else
	{
		for(size_t i=0; i < mCurrentPositions.size(); ++i) aLeaveProbVect[mCurrentPositions[i]] = aProb;
	}
}

void Genes::updateCodonCount(std::vector<unsigned int>& aCodonCounts, unsigned int aSiteMultiplicity) const
{
	if(mCurrentPositions.empty())
	{
		throw FastCodeMLFatal("Invalid codon found in updateCodonCount.");
	}
	else
	{
		for(size_t i=0; i < mCurrentPositions.size(); ++i) aCodonCounts[mCurrentPositions[i]] += aSiteMultiplicity;
	}
}

bool Genes::compareCodons(const char* aCodon1, const char* aCodon2) const
{
	const std::vector<int>& pos1 = getPositions(aCodon1);
	const std::vector<int>& pos2 = getPositions(aCodon2);

	if(pos1.empty() || pos2.empty()) return false;	// Both should be valid
	if(pos1.size() != pos2.size()) return false;	// They should expand to the same number of positions

	for(size_t i=0; i < pos1.size(); ++i) if(pos1[i] != pos2[i]) return false;	// All positions must be equal
	return true;
}


void Genes::checkNameCoherence(const std::vector<std::string>& aNames) const
{
	// Should at least have the same number of species
	if(aNames.size() != mDnaSpecies.size()) throw FastCodeMLFatal("Different number of species in tree and genes.");

	// Create correspondence between species names
	std::vector<std::string>::const_iterator is1=aNames.begin();
	const std::vector<std::string>::const_iterator end1=aNames.end();
	for(; is1 != end1; ++is1)
	{
		bool found = false;
		std::vector<std::string>::const_iterator is2=mDnaSpecies.begin();
		const std::vector<std::string>::const_iterator end2=mDnaSpecies.end();
		for(; is2 != end2; ++is2)
		{
			if(*is1 == *is2) {found = true; break;}
		}
		if(!found) throw FastCodeMLFatal("Mismatch between species in tree and genes.");
	}
}


void Genes::readFile(const char* aFilename, bool aCleanData)
{
	// Read sequences and the corresponding species
	loadData(aFilename, mDnaSpecies, mDnaGene);

	// Finish postprocessing of the read-in sequences
	size_t i, j;

	// Get the number of basis, codons and species
	size_t nbasis = mDnaGene[0].length();
	size_t ncodons = nbasis/3;
	size_t nspecies = mDnaSpecies.size();

	// Check for too many sites (in Forest site*0+class is coded in a unsigned int)
	if(ncodons >= std::numeric_limits<size_t>::max()/10U)
	{
		std::ostringstream o;
		o << "File \"" << aFilename << "\" has too many basis. Max " << std::numeric_limits<size_t>::max()/10U;
		throw FastCodeMLFatal(o);
	}

	// Inizialize codons multiplicity
	std::vector<unsigned int> codon_multiplicity(ncodons, 1);

	// Remove invalid codons
	for(i=0; i < nspecies; ++i)
    {
		const char *p = mDnaGene[i].c_str();
		for(j=0; j < ncodons; ++j)
		{
			if(!validCodon(&p[3*j], aCleanData)) codon_multiplicity[j] = 0;
		}
	}

	// Check if at least one site remains
	size_t valid_codons = std::count(codon_multiplicity.begin(), codon_multiplicity.end(), 1);
	if(valid_codons == 0) throw FastCodeMLFatal("Not a single valid codon read.");

	// Print statistics
	if(mVerboseLevel >= VERBOSE_INFO_OUTPUT)
	{
		std::cerr << std::endl;
		std::cerr << "Num. species: " << std::setw(6) << nspecies << std::endl;
		std::cerr << "Num. basis:   " << std::setw(6) << nbasis << std::endl;
		std::cerr << "Valid codons: " << std::setw(6) << valid_codons << "/" << ncodons << std::endl;
	}

	// Prepare the mapping from program sites back to original sites
	std::multimap<size_t, size_t> sites_back_mapping;
	mOriginalNumSites = ncodons;

	// Remove duplicated sites
	for(i=0; i < ncodons-1; ++i)
	{
		if(codon_multiplicity[i] == 0) continue;
		for(j=i+1; j < ncodons; ++j)
		{
			if(codon_multiplicity[j] == 0) continue;

			unsigned int k=0;
			for(; k < nspecies; ++k)
			{
				const char *p = mDnaGene[k].c_str();
				if(!compareCodons(p+3*i, p+3*j)) break;
			}
			if(k == nspecies)
			{
				++codon_multiplicity[i];
				codon_multiplicity[j] = 0;

				sites_back_mapping.insert(std::pair<size_t,size_t>(i,j));
			}
		}
	}

	// Compute site multiplicity (remove zero multiplicity sites from codon_multiplicity)
	mSiteMultiplicity = codon_multiplicity;
	std::vector<unsigned int>::iterator pend(std::remove(mSiteMultiplicity.begin(), mSiteMultiplicity.end(), 0));
	mSiteMultiplicity.erase(pend, mSiteMultiplicity.end());

	if(mVerboseLevel >= VERBOSE_INFO_OUTPUT)
	{
		std::cerr << "Sites:        " << std::setw(6) << mSiteMultiplicity.size() << "/" << ncodons << std::endl;
		int multi_codons = static_cast<int>(std::count_if(codon_multiplicity.begin(), codon_multiplicity.end(), std::bind2nd(std::greater<unsigned int>(), 1)));
		std::cerr << "Multi codons: " << std::setw(6) << multi_codons << "/" << ncodons << std::endl;
	}
	if(mVerboseLevel >= VERBOSE_MORE_DEBUG)
	{
		std::cerr << std::endl;
		std::ostream_iterator<unsigned int> out_it(std::cerr, " ");
		std::copy(codon_multiplicity.begin(), codon_multiplicity.end(), out_it);
		std::cerr << std::endl;
	}

	// Prepare the map from reduced site num. (j) to list of corresponding original sites (i).
	for(i=j=0; i < codon_multiplicity.size(); ++i)
	{
		if(codon_multiplicity[i] > 0)
		{
			mSitesMappingToOriginal.insert(std::pair<size_t,size_t>(j,i));

			std::multimap<size_t, size_t>::iterator it;
			std::pair<std::multimap<size_t, size_t>::iterator,std::multimap<size_t, size_t>::iterator> ret;
			ret = sites_back_mapping.equal_range(i);
			for(it=ret.first; it != ret.second; ++it)
			{
				mSitesMappingToOriginal.insert(std::pair<size_t,size_t>(j,it->second));
			}
			++j;
		}
	}

	if(mVerboseLevel >= VERBOSE_MORE_DEBUG)
	{
		std::cerr << std::endl;
		for(j=0; j < mSiteMultiplicity.size(); ++j)
		{
			std::multimap<size_t, size_t>::iterator it;
			std::pair<std::multimap<size_t, size_t>::iterator,std::multimap<size_t, size_t>::iterator> ret;
			ret = mSitesMappingToOriginal.equal_range(j);

			std::cout << "Reduced " << std::setw(4) << j << " maps to original:";
			for(it=ret.first; it != ret.second; ++it)
			{
				std::cout << " " << it->second;
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// Compute map from site to position on mDnaGene
	for(unsigned int position_on_gene = 0; position_on_gene < codon_multiplicity.size(); ++position_on_gene)
	{
		if(codon_multiplicity[position_on_gene] > 0) mMapSiteToDnaGene.push_back(position_on_gene);
	}

	// Map from specie name to position in DnaGene
	unsigned int idx = 0;
	std::vector<std::string>::const_iterator is=mDnaSpecies.begin();
	const std::vector<std::string>::const_iterator end=mDnaSpecies.end();
	for(; is != end; ++is, ++idx) mMapSpecieToDnaGene[*is] = idx;
}

void Genes::initFullCodonMap(void)
{
	int i, j, k;

	// Create the list of codons without ambiguities
	const char* base = "TCAG";
	std::map<std::string, int> codons;

	char codon[4];
	codon[3] = '\0';

	for(i=0; i < 4; ++i)
	{
		for(j=0; j < 4; ++j)
		{
			for(k=0; k < 4; ++k)
			{
				// Skip stop codons
			    if(i == 0)
			    {
					if((j == 2) && (k == 2 || k == 3)) continue;
					else if(j == 3 && k == 2) continue;
			    }

			    // Compute the index
			    int idx = i*4*4+j*4+k;

			    // Adjust for missing stop codons
			    if(idx > 14) idx -= 3;
			    else if(idx >  9) idx -= 2;

				// Create the valid codon
			    codon[0] = base[i];
			    codon[1] = base[j];
			    codon[2] = base[k];

				// Add to the map from codon to position in the CPV
			    codons.insert(std::pair<std::string, int>(std::string(codon), idx));
			}
		}
	}

	// Create the list of codons with ambiguous positions
	int mask[4] = {0x1, 0x4, 0x8, 0x2};
	const char* amb = ".TGKCYSBAWRDMHVN";
	char codona[4];
	codona[3] = '\0';
	codon[3] = '\0';

	for(i=1; i < 16; ++i)
	{
		for(j=1; j < 16; ++j)
		{
			for(k=1; k < 16; ++k)
			{
				// Build one of the possible codons (valid and ambiguous)
				codona[0] = amb[i];					
				codona[1] = amb[j];					
				codona[2] = amb[k];					

				std::vector<int> pos;
				bool valid = false;

				// Translate back to the corresponding non-ambiguous codons
				for(int mi=0; mi < 4; ++mi)
				{
					for(int mj=0; mj < 4; ++mj)
					{
						for(int mk=0; mk < 4; ++mk)
						{
							if((i & mask[mi]) && (j & mask[mj]) && (k & mask[mk]))
							{
								codon[0] = base[mi];
								codon[1] = base[mj];
								codon[2] = base[mk];

								// If valid, record the corresponding position
								std::map<std::string, int>::const_iterator im(codons.find(codon));
								if(im == codons.end()) continue;

								valid = true;
								pos.push_back(im->second);
							}
						}
					}
				}

				// This is a valid, possibly ambiguous codon
				if(valid)
				{
					mMapCodonToPosition.insert(std::pair<std::string, std::vector<int> >(std::string(codona), pos));
				}
			}
		}
	}
}


const std::vector<int>& Genes::getPositions(const char* aCodon) const
{
	// Convert to canonical form (only valid uppercase letters)
	char codon[4];
	if(aCodon[0] == '-') codon[0] = 'N';
	else
	{
		char b = toupper(aCodon[0]);
		codon[0] = (b == 'U') ? 'T' : b;
	}
	if(aCodon[1] == '-') codon[1] = 'N';
	else
	{
		char b = toupper(aCodon[1]);
		codon[1] = (b == 'U') ? 'T' : b;
	}
	if(aCodon[2] == '-') codon[2] = 'N';
	else
	{
		char b = toupper(aCodon[2]);
		codon[2] = (b == 'U') ? 'T' : b;
	}
	codon[3] = '\0';

	// Check if it is in the list of valid codons
	std::map<std::string, std::vector<int> >::const_iterator im(mMapCodonToPosition.find(std::string(codon)));

	// If no, return an empty list, else return the list of corresponding positions
	if(im == mMapCodonToPosition.end()) return mEmptyVector;
	return im->second;
}
