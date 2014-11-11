/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * FixEPS.cpp
 *
 *  Created on: 20/07/2009
 *      Author: tom
 */

#include "FixEPS.h"
#include "util/exception.h" 
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <boost/algorithm/string.hpp>
 
const double EPS_SCALE = 1;

void CFixEPS::fixEpsBB(const char * szFile)
{  
    /*typedef char                    char_t;
    typedef file_iterator <char_t>  iterator_t;
    typedef scanner<iterator_t>     scanner_t;
    typedef rule <scanner_t>        rule_t;

    rule_t my_rule;

    file_iterator<> first(szFile);

    CHECK (!first, "Unable to open file!\n");

    file_iterator<> last = first.make_end();

    // Define your rule

    parse_info<iterator_t> info = parse(first, last, my_rule);*/
	std::ostringstream newFileContents;
    std::ifstream inFile(szFile);
    while(!inFile.eof())
    {
    	std::string strLine;
    	getline(inFile, strLine);

    	if(strLine.compare(0, 8, "newpath ") == 0) //no clipping
    		continue;

    	if(strLine.compare(0, strlen("%%BoundingBox: "), "%%BoundingBox: ") == 0)
    	{
			std::vector<std::string> strs;
			boost::split(strs, strLine, boost::is_any_of("\t "));

			/*std::cout << strs[1] << std::endl;
			std::cout << strs[2] << std::endl;
			std::cout << strs[3] << std::endl;
			std::cout << strs[4] << std::endl;*/

			double x0 = atof(strs[1].c_str());
			double y0 = atof(strs[2].c_str());
			double x1 = atof(strs[3].c_str());
			double y1 = atof(strs[4].c_str());
			newFileContents << "%%BoundingBox: " << x0-15*EPS_SCALE << ' ' << y0-15*EPS_SCALE << ' ' << x1+15*EPS_SCALE << ' ' << y1+15*EPS_SCALE << '\n';
    	}
    	else
    		newFileContents << strLine << '\n';
    };

    inFile.close();
    std::ofstream outFile(szFile);
    outFile << newFileContents.str();
    outFile.close();
}

