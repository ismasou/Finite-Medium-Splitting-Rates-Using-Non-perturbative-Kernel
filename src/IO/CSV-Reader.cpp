#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
/*
 * A class to read data from a csv file.
 */
class CSVReader
{
	std::string fileName;
	std::string delimeter;
	int SkipLines;
 
public:
	CSVReader(std::string filename, std::string delm = ",", int SkipLines_ = 0) :
			fileName(filename), delimeter(delm), SkipLines(SkipLines_)
	{ }
 
	// Function to fetch data from a CSV File
	std::vector<std::vector<double> > getData();
};
 
/*
* Parses through csv file line by line and returns the data
* in vector of vector of strings.
*/
std::vector<std::vector< double > > CSVReader::getData()
{
	std::ifstream file(fileName);
 
	std::vector<std::vector<double > > dataList;
 
	std::string line = "";
	for( int i=0;i<SkipLines;i++){
		getline(file, line);
	}
	// Iterate through each line and split the content using delimeter
	while (getline(file, line))
	{
		std::vector<std::string> vec;
        if(line[0]!='#' && line!=""){
			boost::trim_if(line, boost::is_any_of("\t "));
            boost::algorithm::split(vec, line, boost::is_any_of(delimeter),boost::token_compress_on);
            std::vector<double> doubleVector(vec.size());
            for(int i=0;i<int(vec.size());i++){
				// vec[i].erase(std::remove(vec[i].begin(),vec[i].end(),' '),vec[i].end());
                doubleVector[i]=stod(vec[i]);
                // doubleVector[i]=boost::lexical_cast<double>(vec[i]);
				if(!std::isfinite(doubleVector[i])){
					printf("%s, %g\n", vec[i].c_str(), doubleVector[i]);
				}

            }
            dataList.push_back(doubleVector);            
        }
	}
	// Close the File
	file.close();
 
	return dataList;
}