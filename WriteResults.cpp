
#include <iostream>
#include <iomanip>
#include <fstream>
#include "WriteResults.h"


void WriteResults::saveResults(const char* aFilename)
{
	if(!aFilename) return;

	std::ofstream out(aFilename, std::ios_base::trunc | std::ios_base::out);
	if(!out.good())
	{
	}

	out << "Outputting" << std::endl;

	out.close();
}
