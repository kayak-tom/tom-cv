#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include "GRCException.h"

using namespace std;

namespace grc
{
class config
{
	struct ltstr
	{
	  bool operator()(const char* s1, const char* s2) const
	  {
	    return strcmp(s1, s2) < 0;
	  }
	};
    typedef map<const char*, char*, ltstr> TSzMap;
    TSzMap cfgNameVal;
public:
    config(const char * szFileName)
    {
        load(szFileName);
    };

    ~config()
    {

        for(TSzMap::iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++)
        {
            pair<const char*, char*> nv = *pMI;
            free((void*)nv.first);
            free(nv.second);
        }
    };
//#ifndef GCC
    char * strndup_int(char *s, int n)
    {
        if (strlen(s) < (size_t)n ) n=strlen(s);

        char *result = (char*)malloc(n + 1);
        if (result == (char*)0)
    	    return (char*)0;

        strncpy(result, s, n);
        result[n] = 0;
        return result;
    }

//#endif

    void load(const char * szName)
    {
        ifstream fileSettings(szName);
        if(!fileSettings.is_open())
        {
            char szErr[500];
#ifdef GCC
            sprintf(szErr, "Config file not found: \"%s\"", szName);
#else
            sprintf_s(szErr, 500, "Config file not found: \"%s\"", szName);
#endif
            throw new GRCException(szErr);
        };

        while(!fileSettings.eof())
        {
        	//cout <<"Not eof-";
            char szLine[1024] = "";
            fileSettings.getline(szLine, 1024);
            int lineLen = strlen(szLine);
            if(szLine[lineLen-1] == '\r' || szLine[lineLen-1] == '\n')  szLine[lineLen-1] = 0;
            if(szLine[lineLen-2] == '\r' || szLine[lineLen-2] == '\n')  szLine[lineLen-2] = 0;

            //cout << "Read:" << szLine << "\n";
            //cout.flush();
            char * pcName = szLine;
            while (*pcName == ' ') pcName++;
            if(*pcName == '\n' || *pcName == '\0' || strlen(pcName)==0)
            	continue;
            int nNameLen = strchr(pcName, '=')-pcName;
            if(nNameLen == 0) continue;

            char * pcNameEnd = pcName + nNameLen;
            char * pcVal = pcNameEnd + 1;

            while (*(pcNameEnd-1) == ' ' && pcNameEnd>pcName) pcNameEnd--;

            while (*pcVal == ' ') pcVal++;

            int nValLen = strlen(pcVal);
            //if(nValLen==0) continue;
            char * pcValEnd = pcVal + nValLen;
            while (*(pcValEnd-1) == ' ' && pcValEnd>pcVal) pcValEnd--; //Todo: allow comments here

            if(pcNameEnd>pcName /*&& pcValEnd>pcVal*/)
            {
                char * szName = strndup_int(pcName, pcNameEnd-pcName );
                char * szVal = strndup_int(pcVal, pcValEnd-pcVal );

                if(cfgNameVal.find(szName) == cfgNameVal.end())
                	cfgNameVal[szName] = szVal;
                else
                {
                    cout << "WARNING: '" << szName << "' duplicated in cfg file\n";
                    free(szName);
                    free(szVal);
                }

            };
        };

        fileSettings.close();
    };

    GRCException * badNameException(const char * szName) const
    {
            char szErr[500];
#ifdef GCC
            sprintf(szErr, "Parameter not found in config file: \"%s\"", szName);
#else
            sprintf_s(szErr, 500, "Parameter not found in config file: \"%s\"", szName);
#endif
            return new GRCException(szErr);
    };

    const char * getSz(const char * szName)
    {
        if(cfgNameVal.find(szName) == cfgNameVal.end())
        {
            throw badNameException(szName);
        }
        return cfgNameVal[szName];
    };

    double getDouble(const char * szName)
    {
        const char * szVal = getSz(szName);
        if(szVal) return atof(szVal);

        throw badNameException(szName);
    };

    int getInt(const char * szName)
    {
        const char * szVal = getSz(szName);
        if(szVal) return atoi(szVal);

        throw badNameException(szName);
    };
};

//#define cfgLoad(var) get( #var, var ); cout << #var << "=" << var << '\n';

}
