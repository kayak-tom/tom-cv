/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "redirectCout.h"
#include "prfxstream.h"
#include <boost/filesystem.hpp>

redirectCout::redirectCout(const char * szLogfileName, bool bSpawnGedit, std::ostream & cout_or_cerr) : cout_or_cerr(cout_or_cerr)
{
    cout_or_cerr << "Redirecting cout...\n";
    cout_or_cerr.flush();
    pOut = new oprfxstream(cout_or_cerr.rdbuf(), szLogfileName, bSpawnGedit);
    //*pOut << "TEST data sent to oprfxstream" << endl;
    pBackupCout = cout_or_cerr.rdbuf();
    cout_or_cerr.rdbuf(pOut->rdbuf());
    //cout << "TEST data sent to cout" << endl;
};

redirectCout::~redirectCout()
{
    pOut->flush();
    cout_or_cerr.rdbuf(pBackupCout);
    delete pOut;
};

int redirectCout::numBytesPrinted() const
{
    return pOut->numBytesPrinted();
}

std::string makeNiceFolderName()
{
    // Setup logging, windows, output directory, etc.

    std::string strFolderName = "logs/";

    const bool bRemoveLogsOnStartup = false;
    if(bRemoveLogsOnStartup) {
        if (boost::filesystem::exists(strFolderName)) {
            try {
                boost::filesystem::remove_all(strFolderName);
            } catch (...) {
                cout << "Error removing " << strFolderName << endl;
            }
        }
    }

    if (!boost::filesystem::exists(strFolderName)) {
        cout << "Creating directory " << strFolderName << endl;
        boost::filesystem::create_directory(strFolderName);
    }

    strFolderName.append(getTimeAndDate(false));
    strFolderName.append("/");

#ifdef _WIN32
    /**************  rdcm: parsing "/" and ":" to underlines "_" for windows *******************/
    char* cStrFolderName =  new char[strFolderName.size()+1];
    strcpy_s(cStrFolderName, strFolderName.size()+1, strFolderName.c_str());
    char* strFolderNameParts = strtok(cStrFolderName,"/:");
    strFolderName = "logs/";
    while(strFolderNameParts!=NULL) {
        strFolderName.append(strFolderNameParts);
        strFolderName.append("_");
        strFolderNameParts = strtok(NULL,"/:");
    }
    strFolderName.append("/");
    /************** END    rdcm: parsing "/" and ":" to underlines "_" for windows ****************/
#endif

    cout << "Creating directory " << strFolderName << endl;
    boost::filesystem::create_directory(strFolderName);

    return strFolderName;
}