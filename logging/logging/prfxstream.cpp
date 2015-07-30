
#include "prfxstream.h"
#include "util/exception.h"
 
#include <time.h>
#ifdef __GNUC__
#define BS '/'
#else
#define BS '\\'
#endif

prfxbuf::prfxbuf(streambuf *sb, const char * szLogfileName, bool bSpawnGedit):
    streambuf(),
    i_sbuf(sb),
    i_cache(EOF), //not too sure what this is doing
    strLine(),
    nNumBytesPrinted(0),
    nMaxNumBytesPrinted(0)//unlimited
{
    string strFile(szLogfileName);

    fd.open(strFile.c_str());

#ifdef __GNUC__
//#ifdef _DEBUG
    if(bSpawnGedit) {
        char szCommand[200];
        //Was: sprintf(szCommand, "gedit %s &", strFile.c_str()); but hangs when child isn't killed on exit
        sprintf(szCommand, "ps x | egrep 'geX?dit' > temp; if [ -s temp ]; then gedit %s; fi; rm temp &", strFile.c_str());
        if(-1==system(szCommand)) THROW("Error spawning stdout display thread");
    }
//        ps x | egrep 'geX?dit' > temp; if [ -s temp ]; then echo 'gedit running'; fi; rm temp
//#endif
#endif
    setp(0, 0);
    setg(0, 0, 0);
}
prfxbuf::~prfxbuf()
{
    fd.flush();
    fd.close();
}

int    prfxbuf::underflow()
{
    return i_cache;
}


int    prfxbuf::uflow()
{
    if (i_cache == EOF) {
        int rc = i_sbuf->sbumpc();
        return rc;
    }
    int rc = i_cache;
    i_cache = EOF;
    return rc;
}

//static std::ofstream logOverflow("logOverflow.log");

int prfxbuf::overflow(int c) //todo tidy up
{
    boost::mutex::scoped_lock scoped_lock(mx);

    //printf("Print %d\n", c);
    nNumBytesPrinted++;
    if(nMaxNumBytesPrinted>0 && nNumBytesPrinted > nMaxNumBytesPrinted) //After 2MB, dump to file
    {
        printf("\n\n######### Printing only to file ##############\n");
        fflush(stdout);
        cout.rdbuf(fd.rdbuf());
        return 0;
    }
    
    
    if (c != EOF) {
        int rc = i_sbuf->sputc((char)c);
        strLine.push_back((char)c);

        if (c == '\n') {
            // Write string to file
            fd.write(strLine.c_str(), (std::streamsize)strLine.size());
            fd.flush();

            strLine.clear();
        }
        return rc;
    }
    return 0;
}

int    prfxbuf::sync()
{
    i_sbuf->pubsync();
    return 0;
}

iprfxstream::iprfxstream(streambuf *sb, const char *prfx):
    istream(new prfxbuf(sb, prfx, false))
{
}

oprfxstream::oprfxstream(streambuf *sb, const char *szLogfileName, bool bSpawnGedit):
    ostream(new prfxbuf(sb, szLogfileName, bSpawnGedit))
{
}

int oprfxstream::numBytesPrinted() const {
    return (int)static_cast<const prfxbuf *>(rdbuf())->numBytesPrinted();
}

void oprfxstream::setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in) 
{
    static_cast<prfxbuf *>(rdbuf())->setNumBytesPrinted(nMaxNumBytesPrinted_in);
}

iprfxstream::~iprfxstream()
{
    delete rdbuf();
}

oprfxstream::~oprfxstream()
{
    delete rdbuf();
}

std::string getTimeAndDate(bool bAugment) //Todo: move
{
#ifdef __GNUC__
    time_t rawtime;

    time ( &rawtime );
    char * dateTimeStr = ctime(&rawtime);
    char acDate[50];
    static int s_int=0;
    s_int++;
    s_int %= 1000;
    memset(acDate, 0, 50);
    char * pos = strstr(dateTimeStr, "200");//2008 or whatever
    if(pos==0) pos = strstr(dateTimeStr, "201"); //2010
    int len = (pos-dateTimeStr)+4;
    strncpy(acDate, dateTimeStr, len);
    if(bAugment)
        sprintf_s(acDate + strlen(acDate), 50, "-%d", s_int);

    int length = strlen(acDate);
    for(int i=0; i < length; i++)
        if(acDate[i] == ' ')
            acDate[i] = '_';
        else if(acDate[i] == ':')
            acDate[i] = '-';

    return acDate;
#else
    char dateStr [128];
    char timeStr [128];
    _strdate_s( dateStr, 9);
    _strtime_s( timeStr, 9);
    char *acDate = new char[50];
    sprintf_s(acDate, 50, "%s_%s", timeStr, dateStr);
    return acDate;
#endif

}
