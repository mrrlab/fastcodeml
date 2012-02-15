
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "Phylip.h"
#include "Exceptions.h"


void Phylip::loadGenesFile(const char* aFilename)
{
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

	// Check the number of nucleotides read
	for(unsigned int n=0; n < nspecies; ++n)
	{
		if(mDnaGene[n].length() != nbasis)
		{
			std::cerr << "File " << aFilename << " gene " << n << " has wrong number of nucleotides" << std::endl;
			throw FastCodeMLFatalNoMsg();
		}
	}

	// Finish postprocessing of the read in sequences
	postprocessLoadedGenes();
}

