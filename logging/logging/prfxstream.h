// <!-*-C++-*- file: prfxstream.h --->
// <!-------------------------------------------------------------------------->
// <! Copyright (C) 1995 Dietmar Kuehl >
// <!>
// <! This file is free software; you can redistribute it and/or modify >
// <! it under the terms of the GNU General Public License as published by >
// <! the Free Software Foundation; either version 2 of the License, or >
// <! (at your option) any later version. >
// <!>
// <! This program is distributed in the hope that it will be useful, >
// <! but WITHOUT ANY WARRANTY; without even the implied warranty of >
// <! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the >
// <! GNU General Public License for more details. >
//----------------------------------------------------------------------------
// <HTML>
// <HEAD>
// <TITLE> prfxstream.h - Header File for a Specialized IOstream </TITLE>
// <H1>Header File for a Specialized IOstream</H1>
// <H2>Reading/Writing files with a per line prefix</H2>
// <HR>
// </HEAD>
#ifndef _PRFXSTREAM_H_
#define _PRFXSTREAM_H_
// </PRE>
//----------------------------------------------------------------------------
#include <iostream>
#include <string.h>
#include<fstream>
#include <boost/thread/mutex.hpp>
// #include <streambuf>
using namespace std;
//#include<iosfwd>
// </PRE> 
class prfxbuf: public streambuf
{
private:
    streambuf *i_sbuf; // the actual streambuf used to read and write chars
//unsigned int i_len; // the length of the prefix
//char *i_prfx; // the prefix
//bool i_newline; // remember whether we are at a new line
    int i_cache; // may cache a read character
    std::string strLine;
    ofstream fd;
    boost::mutex mx;
    size_t nNumBytesPrinted;
    size_t nMaxNumBytesPrinted;
//bool skip_prefix();
protected:
    int overflow(int);
    int underflow();
    int uflow();
    int sync();
public:
    prfxbuf(streambuf *sb, const char * szExeName, bool bSpawnGedit);
    ~prfxbuf();
    size_t numBytesPrinted() const {
        return nNumBytesPrinted;
    }
    
    void setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in) 
    {
        nMaxNumBytesPrinted = nMaxNumBytesPrinted_in;
    }
};
class iprfxstream: public istream
{
public:
    iprfxstream(streambuf *sb, const char *prfx);
    ~iprfxstream();
};
class oprfxstream: public ostream
{
public:
    oprfxstream(streambuf *sb, const char *prfx, bool bSpawnGedit);
    ~oprfxstream();
    int numBytesPrinted() const;
    
    void setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in);
};
std::string getTimeAndDate(bool bAugment);

//Format: ./logs/TIMESTAMP/output.log
std::string makeNiceFolderName();

#endif /* _PRFXSTREAM_H_ */
