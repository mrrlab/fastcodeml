
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "Phylip.h"
#include "Exceptions.h"

void Phylip::loadData(const char* aFilename, std::vector<std::string>& aSpecies, std::vector<std::string>& aSequences)
{
	// Open the file and verify it has at least one line
	std::ifstream in(aFilename);
	if(!in)
	{
		std::ostringstream o;
		o << "Cannot open gene file \"" << aFilename << '"';
		throw FastCodeMLFatal(o);
	}

    std::string str;
	if(!getline(in, str) || str.empty())
	{
		in.close();
		std::ostringstream o;
		o << "File \"" << aFilename << "\" is empty";
		throw FastCodeMLFatal(o);
	}

	// From the first line extract number of species and nmber of basis
	long unsigned int nspecies, nbasis;
	char *endptr;
	const char *next = str.c_str();
	nspecies = strtol(next, &endptr, 10);
	if(endptr == next)
	{
		in.close();
		std::ostringstream o;
		o << "File \"" << aFilename << "\" is malformed";
		throw FastCodeMLFatal(o);
	}

	next = endptr;
	nbasis = strtol(next, &endptr, 10);
	if(endptr == next)
	{
		in.close();
		std::ostringstream o;
		o << "File \"" << aFilename << "\" is malformed";
		throw FastCodeMLFatal(o);
	}

	// Read and parse the genes
    while(aSpecies.size() < nspecies && getline(in, str))
    {
		// Extract the specie name
        if(str.empty()) continue;
		size_t p1 = str.find_first_not_of(" \t\r");
		if(p1 == std::string::npos) continue;
		size_t p2 = str.find_first_of(" \t", p1);

		std::string s;
		s.assign(str, p1, p2-p1);
		aSpecies.push_back(s);

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
        aSequences.push_back(s);
	}
	in.close();

	// Check correct number of species loaded
	if(nspecies != aSpecies.size())
	{
		std::ostringstream o;
		o << "File \"" << aFilename << "\" has number of species mismatch";
		throw FastCodeMLFatal(o);
	}

	// Check the number of nucleotides read
	for(unsigned int n=0; n < nspecies; ++n)
	{
		if(aSequences[n].length() != nbasis)
		{
			std::ostringstream o;
			o << "File \"" << aFilename << "\" gene " << n << " has wrong number of nucleotides";
			throw FastCodeMLFatal(o);
		}
	}
	
	// Other sanity checks
	if(nbasis % 3)
	{
		std::ostringstream o;
		o << "File \"" << aFilename << "\" number of basis is not multiple of 3";
		throw FastCodeMLFatal(o);
	}
}

