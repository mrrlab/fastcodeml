
#ifndef PARSEPARAMETERS_H
#define PARSEPARAMETERS_H

#include <map>
#include <fstream>
#include <string>

/// Parse and accumulate initial values for parameters.
/// New values for the parameters came in the form PPP=VVV
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-02-21 (initial version)
///     @version 1.0
///
class ParseParameters
{
public:
	/// Return a pointer to the singleton instance
	///
	/// @return The pointer to the instance
	///
	static ParseParameters* getInstance(void);

	void addParameter(const char* aParamValuePair);

	double getParameter(const char* aParamName) const;

	/// Print the class statistics as: cout << r;
	///
	/// @param[in] aOut Output stream
	/// @param[in] aForest The forest to be printed
	///
	/// @return The output stream
	///
	friend std::ostream& operator<< (std::ostream& aOut, const ParseParameters* aParamsList);

protected:
	/// Protected constructor
	///
	ParseParameters();


private:
	static ParseParameters*			mInstance;					///< Pointer to the singleton instance
	std::map<std::string, double>	mDictionary;
};

#endif

