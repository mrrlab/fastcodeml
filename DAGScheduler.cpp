
#include <iostream>
#include <iomanip>
#include "DAGScheduler.h"

void DAGScheduler::clear(void)
{
	mNodes.clear();
	mEdges.clear();
}


void DAGScheduler::loadDependency(const void* aFirst, const void* aDependant)
{
	//std::cerr << std::hex << aFirst << " -> " << std::hex << aDependant << std::endl;

	mNodes.insert(aFirst);
	mNodes.insert(aDependant);

	mEdges.push_back(std::make_pair(aFirst, aDependant));

}


void DAGScheduler::endLoadDependencies(void)
{
}


void* DAGScheduler::getNext(void)
{
	return 0;
}


void DAGScheduler::setDone(const void* aItem)
{
}


void DAGScheduler::resetDependencies(void)
{
}


void DAGScheduler::dumpDAG(std::ostream& aOut) const
{
	std::set<const void*>::const_iterator ind = mNodes.begin();
	for(; ind != mNodes.end(); ++ind)
	{
		size_t v = reinterpret_cast<size_t>(*ind);
		aOut << std::hex << *ind << ' ' << std::setw(6) << v << std::endl;
	}
	aOut << '#' << std::endl;

	std::vector<std::pair<const void*, const void*> >::const_iterator ied = mEdges.begin();
	for(; ied != mEdges.end(); ++ied)
	{
		aOut << ied->first << ' ' << ied->second << std::endl;
	}
}

