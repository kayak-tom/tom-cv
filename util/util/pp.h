#ifndef PP_H
#define PP_H

/*
 * Macros for debugging output.
 * 
 * Turn on debugging with const bool bVerbose = true;
 * 
 * Then use COUT macros to priunt variable names.
 * 
 * Eigen types automatically transposed.
 * 
 * "Strings" are only output once.
 * 
 * ::toString(variable) gives a std::string (e.g. for appending numbers to a string)
 * TO_STRING(variable) gives a std::string with the variable name. std::string("variable=1234")
 * 
 * COUT(variable);
 * 
 * > variable = 1242
 * 
 * COUT2("message", variable);
 * 
 * > message: variable = 1242
 * 
 * */
 
#include <sstream> 
#include <iomanip> //setprecision
#include <iostream>

namespace boost
{
    template<typename T>
    class optional;
}

using std::cout;
using std::endl;
using std::string;
using boost::optional;

template<typename T1, typename T2>
std::ostream & operator<<(std::ostream & s, const std::pair<T1, T2> & info)
{
    s << info.first << ":" << info.second;
    return s;
}

template<class T>
std::ostream & operator<<(std::ostream & s, const optional<T> & info)
{
    s << "Optional";
    if (info)
        s << *info;
    else
        s << '-';
    return s;
}

template<typename T>
void doCout(const char * label, const T & t)
{
    std::ostringstream ss;
    ss << std::setprecision(6);
    ss << t;
    const bool bMultiline = (ss.str().find('\n') != std::string::npos);
    const bool bIsChar = (label[0] == '\"'); // COUT("Hello") syntax

    if(!bIsChar) {
        cout << label;
        cout << " = ";
    }

    if(bMultiline)
        cout << endl;

    cout << std::dec /*Try to stop random hex output*/ << ss.str() << endl;
}

#define COUT(var) if(bVerbose) doCout(#var, (var));
#define COUT2(label, var) if(bVerbose) { cout << label << " "; doCout(#var, (var)); }

#define COUT_CERR(...) { std::ostringstream ss; ss << __VA_ARGS__; cout << ss.str() << endl; std::cerr << ss.str() << endl; }

template<typename T>
std::string toString(const T & t)
{
    std::ostringstream ss;
    ss << std::setprecision(4) << t;
    return ss.str();
}

#define COUTPAIR(var1,var2) if(bVerbose) doCout(#var1 "," #var2, ::toString(var1)+ "," +::toString(var2));

#define TO_STRING(t) (#t "=" + ::toString(t) + " ")

#define PPSTR(str)  (# str " length=" + ::toString(str.length()) + " " + str.substr(0,50) + (str.length() > 50 ? "... " : " "))

#define ALWAYS_VERBOSE const bool bVerbose = true;

#endif // PP_H