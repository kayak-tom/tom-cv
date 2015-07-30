/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "prfxstream.h"

class redirectCout
{
    std::streambuf * pBackupCout;
    oprfxstream    * pOut;
    std::ostream & cout_or_cerr; //the stream that's redirected

    //char * szSaveLoc;
public:
    redirectCout(const char * szLogfilePath, bool bSpawnGedit, std::ostream & cout_or_cerr = std::cout);
    ~redirectCout();
    
    int numBytesPrinted() const;
    
    void setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in) 
    {
        pOut->setNumBytesPrinted(nMaxNumBytesPrinted_in);
    }    
};
