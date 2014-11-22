/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include "util/exception.h"
#include "stdlib.h"
#include "string.h"
#ifdef __GNUC__ 
#include <unistd.h> 
#endif

#ifndef __GNUC__
#include "windows.h"
#endif
#include <algorithm>

class config
{
    struct ltstr
    {
      bool operator()(const char* s1, const char* s2) const
      {
         int res = strcmp(s1, s2);
         //cout << "comparing " << s1 << " and "<< s2 << " (" << res << ")" << endl << flush;
         return res < 0;
      }
    };

protected:
    typedef std::map<const char*, char*, ltstr> TSzMap;
    typedef std::set<const char*> TAccessedParams;
    TSzMap cfgNameVal;
    TAccessedParams accessedParams;
public:
    const char * getSz(const char * szName) const
    {
        //cout << "Searching for " << szName << " in "; dumpNames();
        if(szName && strlen(szName)>0)
        {
            TSzMap::const_iterator pName = cfgNameVal.find(szName);
            if(pName == cfgNameVal.end())
            {
                //std::cout << "WARNING: '" << szName << "' not in map\n";
                return 0;
            }

            //Flag as accessed
            const_cast<config *>(this)->accessedParams.insert(pName->first);

            //std::cout << "Name: " << pName->first << std::endl;
            //std::cout << "Val: " << pName->second << std::endl;

            return pName->second;
        }
        return 0;
    };

    static void safeSleep(int nMillisecs)
    {
    #ifdef __GNUC__
        usleep(1000*nMillisecs);
    #else
        Sleep(nMillisecs);
    #endif
    }

    static void cfg_strncpy(const char * szOriginal, char * szNew, int n)
    {
#ifdef __GNUC__
        strncpy(szNew, szOriginal, n);
#else
        strncpy_s(szNew, n+1, szOriginal, n);
#endif
    }

    static char * cfg_strndup(const char * szOriginal, int n)
    {
        char * szNew = new char[n+1];
        cfg_strncpy(szOriginal, szNew, n);
/*#ifdef __GNUC__
        strncpy(szNew, szOriginal, n);
#else
        strncpy_s(szNew, n+1, szOriginal, n);
#endif*/
        int nEnd = std::min<int>(n, (int)strlen(szOriginal));
        szNew[nEnd] = 0;
        return szNew;
    }
protected:
    static bool invalid(const char c)
    {
        switch  (c)
        {
            case ' ':
            case '\357'://'UTF-8 byte order marks' apparently
            case '\273':
            case '\277':
            case '\n':
            case '\r':
            case '\t':
                return true;
        }
        return false;
    }

    void load(const char * szName)
    {
        std::ifstream fileSettings(szName);
        while(!fileSettings.is_open())
        {
            std::cout << "Path to config file must be given as a command line argument (" << szName << ")\n";
            std::cout.flush();
            fileSettings.close();
            safeSleep(1000);
            fileSettings.open(szName);
        };
        
        int nLoopCount = 0;
        while(!fileSettings.eof() && nLoopCount < 1000)
        {
            nLoopCount++;
            //cout <<"Not eof-";
            const int nMaxStrLen = 1024;
            char szLine[1024] = "";
            fileSettings.getline(szLine, nMaxStrLen);
            //cout << "Read:" << szLine << "\n";
            //cout.flush();
            char * pcHashPos = strchr(szLine, '#');
            if(pcHashPos) *pcHashPos=0;

            char * pcName = szLine;

            while (invalid(*pcName) ) pcName++;
            
            if(*pcName == '\n' || *pcName == '\0' || strlen(pcName)==0)
                continue; //Includes empty lines
                
            const char * pcEqualsPos = strchr(pcName, '=');
            
            CHECK_P(pcEqualsPos == 0, szLine, "No '=' in config file line");
            
            int nNameLen = (int)(pcEqualsPos - pcName);
            if(nNameLen <= 0 || nNameLen>=nMaxStrLen) //continue; //no string-val pair found
            {
                cout << "Load line " << nLoopCount << " contents=\"" << szLine << "\" from config file failed" << endl;
                THROW("Load line from config file failed")
            }

            char * pcNameEnd = pcName + nNameLen;
            char * pcVal = pcNameEnd + 1;

            while (*(pcNameEnd-1) == ' ' && pcNameEnd>pcName) pcNameEnd--;

            while (*pcVal == ' ') pcVal++;

            int nValLen = (int)strlen(pcVal);
            if(nValLen==0) continue;
            char * pcValEnd = pcVal + nValLen;
            while (invalid(*(pcValEnd-1)) && pcValEnd>pcVal) pcValEnd--;

            if(pcNameEnd>pcName && pcValEnd>pcVal)
            {
                char * szName = cfg_strndup(pcName, (int)(pcNameEnd-pcName));
                char * szVal = cfg_strndup(pcVal, (int)(pcValEnd-pcVal));
                //std::cout << "Loaded param name=" << szName << ", val=" << szVal << std::endl;
                insert(szName, szVal);
            };
        };
        
        CHECK(nLoopCount>=1000, "Error reading config file line-by-line");

        fileSettings.close();
    };
    void insert(const char * szName, char * szVal)
    {
        if(IS_DEBUG) CHECK(!szName || !szVal || strlen(szName)==0 || strlen(szVal)==0 || strlen(szName)>500 || strlen(szVal)>500, "Insert: bad data" );
        if(cfgNameVal.find(szName) == cfgNameVal.end())
            cfgNameVal[szName] = szVal;
        else
        {
            std::cout << "Duplicate parameter in config file:" << szName << "\n";
            delete [] szName;
            delete [] szVal;
            THROW("Duplicate parameter in config file");
        }
    }

public:
    config(bool bAutotest = false)
    {
        if(!bAutotest)
            load("rs.cfg");
        /*char * szName = strdup("SCORE");
        char * szVal = strdup("0");
        insert(szName, szVal);
        cout << "Loaded " << cfgNameVal.size() << " params (should be 0)" << endl;*/
    };
    config(const char * szFileName)
    {
        //if(IS_DEBUG) CHECK(!boost::filesystem::exists(szFileName), "Missing config file");
        load(szFileName);
        std::cout << "Loaded " << cfgNameVal.size() << " params from " << szFileName << std::endl;
    };

    void dumpNames() const
    {
        std::cout << cfgNameVal.size() << ':';
        for(TSzMap::const_iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++)
        {
            //pair<const char*, char*> nv = *pMI;
            std::cout << pMI->first << ',';
        }
        std::cout << std::flush << std::endl;
    };
    
    void checkForUnusedParameters() const
    {
        for(TSzMap::const_iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++)
        {
            std::pair<const char*, char*> nv = *pMI;
            if(accessedParams.find(nv.first) == accessedParams.end())
            {
                std::cout << "Parameter " << nv.first << " in config file not used\n";
                THROW("A parameter in the config file isn't used--it doesn't exist in the 'params' class");
            }
        }
    }

    ~config()
    {
        for(TSzMap::iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++)
        {
            std::pair<const char*, char*> nv = *pMI;

            if(accessedParams.find(nv.first) == accessedParams.end())
                std::cout << "Parameter " << nv.first << " in config file not used\n";

            delete [] nv.first;
            delete [] nv.second;
        }
    };

    void get(const char * szName, int & val) const
    {
        const char * szVal = getSz(szName);
        if(szVal)val=atoi(szVal);
    };
    void get(const char * szName, double & val) const
    {
        const char * szVal = getSz(szName);
        if(szVal)val=atof(szVal);
    };
    void get(const char * szName, char * szValOut) const
    {
        const char * szVal = getSz(szName);
        if(szVal) cfg_strncpy(szVal, szValOut, 500);
    };
};

#define cfgLoad(var) get( #var, var ); cout << #var << "=" << var << '\n';
#define cfgLoadEnum(var) { int nvar; cfg.get( #var, nvar ); var = (var ## s)nvar; cout << #var << "=" << var << '\n'; }


