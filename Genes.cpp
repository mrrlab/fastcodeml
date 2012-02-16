
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include "Genes.h"
#include "Exceptions.h"



Genes::~Genes()
{
	clear();
}

void Genes::clear(void)
{
	mDnaSpecies.clear();		
	mDnaGene.clear();			
	mSiteMultiplicity.clear();	
	mMapSiteToDnaGene.clear();	
	mMapSpecieToDnaGene.clear();
}


int Genes::idxCodon(const char* aCodon) const
{
	int i, j, k;
#if 1
	std::map<char, int>::const_iterator im = mMapBaseToIdx.find(aCodon[0]);
	if(im == mMapBaseToIdx.end()) return -1;
	i = im->second;
	if(aCodon[0] == aCodon[1])
	{
		j = i;
	}
	else
	{
		im = mMapBaseToIdx.find(aCodon[1]);
		if(im == mMapBaseToIdx.end()) return -1;
		j = im->second;
	}
	if(aCodon[0] == aCodon[2])
	{
		k = i;
	}
	else if(aCodon[1] == aCodon[2])
	{
		k = j;
	}
	else
	{
		im = mMapBaseToIdx.find(aCodon[2]);
		if(im == mMapBaseToIdx.end()) return -1;
		k = im->second;
	}
#else
	switch(aCodon[0])
	{
	case 'T':
	case 't':
		i = 0;
		break;
	case 'C':
	case 'c':
		i = 1;
		break;
	case 'A':
	case 'a':
		i = 2;
		break;
	case 'G':
	case 'g':
		i = 3;
		break;
	default:
		return -1;
	}

	switch(aCodon[1])
	{
	case 'T':
	case 't':
		j = 0;
		break;
	case 'C':
	case 'c':
		j = 1;
		break;
	case 'A':
	case 'a':
		j = 2;
		break;
	case 'G':
	case 'g':
		j = 3;
		break;
	default:
		return -1;
	}

	switch(aCodon[2])
	{
	case 'T':
	case 't':
		k = 0;
		break;
	case 'C':
	case 'c':
		k = 1;
		break;
	case 'A':
	case 'a':
		k = 2;
		break;
	case 'G':
	case 'g':
		k = 3;
		break;
	default:
		return -1;
	}
#endif
	// Check if it is a stop codon
	if(i == 0)
	{
		if((j == 2) && (k == 2 || k == 3)) return -1;
		else if(j == 3 && k == 2) return -1;
	}

	// Compute the index
	int idx = i*4*4+j*4+k;

	// Adjust for missing stop codons
	if(idx > 14) return idx-3;
	if(idx >  9) return idx-2;
	return idx;
}


int Genes::getCodonIdx(std::string aSpecie, unsigned int aSite) const
{
	// Find the specie
	const unsigned int idx = mMapSpecieToDnaGene.find(aSpecie)->second;

	// Access its gene
	const char* gene = mDnaGene[idx].c_str();

	// Convert the site to the position on the gene
	unsigned int position_on_gene = mMapSiteToDnaGene[aSite];

	// Return the index for the codon
	return idxCodon(gene+3*position_on_gene);
}


void Genes::checkNameCoherence(const std::vector<std::string>& aNames) const
{
	// Should at least have the same number of species
	if(aNames.size() != mDnaSpecies.size()) throw FastCodeMLFatal("Different number of species in tree and genes");

	// Create correspondence between species names
	std::vector<std::string>::const_iterator is1, is2;
	for(is1=aNames.begin(); is1 != aNames.end(); ++is1)
	{
		bool found = false;
		for(is2=mDnaSpecies.begin(); is2 != mDnaSpecies.end(); ++is2)
		{
			if(*is1 == *is2) {found = true; break;}
		}
		if(!found) throw FastCodeMLFatal("Mismatch between species in tree and genes");
	}
}


void Genes::readFile(const char* aFilename)
{
	// Read sequences and the corresponding species
	loadData(aFilename, mDnaSpecies, mDnaGene);

	// Finish postprocessing of the read-in sequences
	unsigned int i, j;

	// Get the number of basis and species
	unsigned int nbasis = mDnaGene[0].length();
	unsigned int nspecies = mDnaSpecies.size();

	// Inizialize codons multiplicity
	std::vector<unsigned int> codon_multiplicity(nbasis/3, 1);

	// Remove invalid codons
	for(i=0; i < nspecies; ++i)
    {
		const char *p = mDnaGene[i].c_str();
		for(j=0; j < nbasis/3; ++j)
		{
			if(!validCodon(&p[3*j])) codon_multiplicity[j] = 0;
		}
	}

	if(mVerboseLevel >= 1)
	{
		std::cerr << std::endl;
		std::cerr << "Num. species: " << std::setw(6) << nspecies << std::endl;
		std::cerr << "Num. basis:   " << std::setw(6) << nbasis << std::endl;
		int valid_codons = std::count(codon_multiplicity.begin(), codon_multiplicity.end(), 1);
		std::cerr << "Valid codons: " << std::setw(6) << valid_codons << "/" << nbasis/3 << std::endl;
	}

	// Remove duplicated sites
	for(i=0; i < nbasis/3-1; ++i)
	{
		if(codon_multiplicity[i] == 0) continue;
		for(j=i+1; j < nbasis/3; ++j)
		{
			if(codon_multiplicity[j] == 0) continue;

			unsigned int k;
			for(k=0; k < nspecies; ++k)
			{
				const char *p = mDnaGene[k].c_str();
				if(p[3*i+0] != p[3*j+0] || p[3*i+1] != p[3*j+1] || p[3*i+2] != p[3*j+2]) break;
			}
			if(k == nspecies)
			{
				++codon_multiplicity[i];
				codon_multiplicity[j] = 0;
			}
		}
	}

	// Compute site multiplicity (remove zero multiplicity sites from codon_multiplicity)
	mSiteMultiplicity = codon_multiplicity;
	std::vector<unsigned int>::iterator pend = std::remove(mSiteMultiplicity.begin(), mSiteMultiplicity.end(), 0);
	mSiteMultiplicity.erase(pend, mSiteMultiplicity.end());

	if(mVerboseLevel >= 1)
	{
		std::cerr << "Sites:        " << std::setw(6) << mSiteMultiplicity.size() << "/" << nbasis/3 << std::endl;
		int multi_codons = std::count_if(codon_multiplicity.begin(), codon_multiplicity.end(), std::bind2nd(std::greater<unsigned int>(), 1));
		std::cerr << "Multi codons: " << std::setw(6) << multi_codons << "/" << nbasis/3 << std::endl;
	}
	if(mVerboseLevel >= 3)
	{
		std::ostream_iterator<unsigned int> out_it(std::cout, " ");
		std::copy(codon_multiplicity.begin(), codon_multiplicity.end(), out_it);
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
	for(; is != mDnaSpecies.end(); ++is, ++idx)
	{
		mMapSpecieToDnaGene[*is] = idx;
	}
}

