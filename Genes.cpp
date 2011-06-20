
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include "Genes.h"
#include "MatrixSize.h"
#include "Exceptions.h"

Genes::Genes(int aVerboseLevel)
{
	mVerboseLevel = aVerboseLevel;
}


Genes::~Genes()
{
}

void Genes::clear(void)
{
	mDnaSpecies.clear();		
	mDnaGene.clear();			
	mCodonMultiplicity.clear();	
	mSiteMultiplicity.clear();	
	mMapSiteToDnaGene.clear();	
	mMapSpecieToDnaGene.clear();
}

void Genes::loadGenesFile(const char* aFilename)
{
	int i;
	unsigned int j;

	std::ifstream in(aFilename);
	if(!in)
	{
		std::cerr << "Cannot open " << aFilename << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

    std::string str;
	if(!getline(in, str) || str.empty())
	{
		in.close();
		std::cerr << "File " << aFilename << " is empty" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	long unsigned int nspecies, nbasis;
	char *endptr;
	const char *next = str.c_str();
	nspecies = strtol(next, &endptr, 10);
	if(endptr == next)
	{
		in.close();
		std::cerr << "File " << aFilename << " is malformed" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	next = endptr;
	nbasis = strtol(next, &endptr, 10);
	if(endptr == next)
	{
		in.close();
		std::cerr << "File " << aFilename << " is malformed" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	if(mVerboseLevel >= 1) std::cerr << std::endl << "Nspecies: " << nspecies << std::endl << "Nbasis:   " << nbasis << std::endl << std::endl;

	// Read and parse the genes
    while(mDnaSpecies.size() < nspecies && getline(in, str))
    {
		// Extract the specie name
        if(str.empty()) continue;
		size_t p1 = str.find_first_not_of(" \t\r");
		if(p1 == std::string::npos) continue;
		size_t p2 = str.find_first_of(" \t", p1);

		std::string s;
		s.assign(str, p1, p2-p1);
		mDnaSpecies.push_back(s);

		// Extract the gene specification
		s.clear();
		for(;;)
		{
			for(;;)
			{
				p1 = str.find_first_not_of(" \t\r", p2);
				if(p1 == std::string::npos) break;

				p2 = str.find_first_of(" \t\r", p1);
				if(p2 == std::string::npos) p2 = str.size();

				s.append(str, p1, p2-p1);
			}
			if(s.size() >= nbasis) break;
			getline(in, str);
			p2 = 0;
		}
        mDnaGene.push_back(s);
	}
	in.close();

	// Check correct number of species loaded
	if(nspecies != mDnaSpecies.size())
	{
		std::cerr << "File " << aFilename << " has number of species mismatch" << std::endl;
		throw FastCodeMLFatalNoMsg();
	}

	// Inizialize codons multiplicity
	mCodonMultiplicity.assign(nbasis/3, 1);

	// Remove invalid codons
	for(i=0; i < (int)nspecies; ++i)
    {
		const char *p = mDnaGene[i].c_str();
		for(j=0; j < nbasis/3; ++j)
		{
			if(!validCodon(&p[3*j])) mCodonMultiplicity[j] = 0;
		}
	}
	
	if(mVerboseLevel >= 1)
	{
		int valid_codons = std::count(mCodonMultiplicity.begin(), mCodonMultiplicity.end(), 1);
		std::cerr << "Valid codons: " << std::setw(6) << valid_codons << "/" << nbasis/3 << std::endl;
	}

	// Remove duplicated sites
	for(i=0; i < (int)nbasis/3-1; ++i)
	{
		if(mCodonMultiplicity[i] == 0) continue;
		for(j=i+1; j < nbasis/3; ++j)
		{
			if(mCodonMultiplicity[j] == 0) continue;

			unsigned int k;
			for(k=0; k < nspecies; ++k)
			{
				const char *p = mDnaGene[k].c_str();
				if(p[3*i+0] != p[3*j+0] || p[3*i+1] != p[3*j+1] || p[3*i+2] != p[3*j+2]) break;
			}
			if(k == nspecies)
			{
				++mCodonMultiplicity[i];
				mCodonMultiplicity[j] = 0;
			}
		}
	}

	// Compute site multiplicity
	mSiteMultiplicity.clear();

    std::vector<int>::const_iterator icm;
    for(icm=mCodonMultiplicity.begin(); icm != mCodonMultiplicity.end(); ++icm)
    {
        if(*icm < 1) continue;

		unsigned int u = *icm;
		mSiteMultiplicity.push_back(u);
    }

	if(mVerboseLevel >= 1)
	{
		std::cerr << "Sites:        " << std::setw(6) << mSiteMultiplicity.size() << "/" << nbasis/3 << std::endl;
		//int multi_codons = 0;
		//for(j=0; j < nbasis/3; ++j) if(mCodonMultiplicity[j] > 1) ++multi_codons;
		int multi_codons = std::count_if(mCodonMultiplicity.begin(), mCodonMultiplicity.end(), std::bind2nd(std::greater<int>(), 1));
		std::cerr << "Multi codons: " << std::setw(6) << multi_codons << "/" << nbasis/3 << std::endl;
	}
	if(mVerboseLevel >= 3)
	{
		for(j=0; j < nbasis/3; ++j) std::cout << " " << mCodonMultiplicity[j];
		std::cout << std::endl;
	}

	// Compute map from site to position on mDnaGene
	for(unsigned int position_on_gene = 0; position_on_gene < mCodonMultiplicity.size(); ++position_on_gene)
	{
		if(mCodonMultiplicity[position_on_gene] > 0) mMapSiteToDnaGene.push_back(position_on_gene);
	}

	// Map from specie name to position in DnaGene
	std::vector<std::string>::const_iterator is;
	unsigned int idx = 0;
	for(is=mDnaSpecies.begin(); is != mDnaSpecies.end(); ++is, ++idx)
	{
		mMapSpecieToDnaGene[*is] = idx;
	}
}


int Genes::idxCodon(const char* aCodon) const
{
	int i, j, k;

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
	unsigned int idx = mMapSpecieToDnaGene.find(aSpecie)->second;

	// Access its gene
	const char* gene = mDnaGene[idx].c_str();

	// Convert the site to the position on the gene
	unsigned int position_on_gene = mMapSiteToDnaGene[aSite];

	// Return the index for the codon
	return idxCodon(gene+3*position_on_gene);
}
