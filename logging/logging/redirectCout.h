/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "prfxstream.h"

class redirectCout
{
    std::streambuf * pBackupCout;
    oprfxstream    * pOut;

    //char * szSaveLoc;
public:
    redirectCout(const char * szArgv, bool bSpawnGedit);
    ~redirectCout();
    
    int numBytesPrinted() const;
    
    void setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in) 
    {
        pOut->setNumBytesPrinted(nMaxNumBytesPrinted_in);
    }    
};
