/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "redirectCout.h"
#include "prfxstream.h"

redirectCout::redirectCout(const char * szLogfileName, bool bSpawnGedit)
{
    cout << "Redirecting cout...\n";
    cout.flush();
    pOut = new oprfxstream(cout.rdbuf(), szLogfileName, bSpawnGedit);
    //*pOut << "TEST data sent to oprfxstream" << endl;
    pBackupCout = cout.rdbuf();
    cout.rdbuf(pOut->rdbuf());
    //cout << "TEST data sent to cout" << endl;
};

redirectCout::~redirectCout()
{
    pOut->flush();
    cout.rdbuf(pBackupCout);
    delete pOut;
};

int redirectCout::numBytesPrinted() const
{
    return pOut->numBytesPrinted();
}
