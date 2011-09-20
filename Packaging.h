
#ifndef PACKAGING_H
#define PACKAGING_H

#include <vector>

class Packaging
{
public:
	Packaging();
	~Packaging();
	
private:
	std::vector<int>		mNodePresent;			///< -2 if the correstponding: Branch -> Site exists
														///< -1 if doesn't exist
														///< The site number from which the value is taken
	enum {
		SITE_EXISTS     = -2,							///< The position (Branch, Site) in mNodePresent exists
		SITE_NOT_EXISTS = -1							///< The position (Branch, Site) in mNodePresent refers to a not existend node
	};
};

#endif

