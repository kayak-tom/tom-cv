/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * param.h
 *
 *  Created on: 15/10/2009
 *      Author: hilandtom
 */

#ifndef PARAM_H_
#define PARAM_H_

#include "string.h"
#include "stdlib.h"
#include <set>
#include <string>
#include <iostream>
#include "util/exception.h"
#include <boost/noncopyable.hpp>

#ifndef __GNUC__
#  pragma warning(disable:4355) // disable 'this' used in c'tor initialisation warning from msvc
#  pragma warning(push)
#endif

class config;

class CParamHTML : private boost::noncopyable
{
    std::ostream & os;
public:
    CParamHTML(std::ostream & os);

    void addName(const std::string & name);
    void addMMD(double min, double max, double def);
    void addDef(bool def);
    void addDef(const std::string & def);
    void addExplanation(const std::string & explanation);

    ~CParamHTML();
};

class CParamClass;

class CParam : private boost::noncopyable
{
    friend class CParamClass;
    friend class CParamSort;

    //static const int MAX_NAME_LENGTH = 240;
    //const char * NAME, * DESCRIPTION; //Reduce memory use locally (i.e. in cache) TODO: Entire param elsewhere?
    const std::string NAME, DESCRIPTION;
    CParam() { THROW("No default constructor"); } //No default construction either
protected:
    class CParamSort
    {
    public:
       bool operator( ) ( const CParam * a,
                     const CParam * b )
       {
          return a->NAME < b->NAME; //strncmp( b->NAME, a->NAME, CParam::MAX_NAME_LENGTH ) > 0;
       }
    };

    const CParamClass * pParent;

    enum eParamStatus{eUninit, eInit};
    enum eParamUseStatus{eUnused, eUsed};

    eParamStatus initialised;
    eParamUseStatus used;

    void initName(const char * szName);

    virtual void initChildren(config * pConfig);
    virtual void initValue(const char * szVal) = 0;

    bool isUsed() const
    {
        return (used==eUsed);
    }

    enum ePrintAllSettings { ePrintInitUsed, ePrintInitUnused, ePrintUsed, ePrintUnused, ePrintAll };

    bool printMe(ePrintAllSettings eMode) const;

    void setUsed() const
    {
        const_cast<CParam *>(this)->used = eUsed;
    }

    virtual void printAll(ePrintAllSettings eMode) const = 0;
    virtual void printAllCheck(ePrintAllSettings eMode) const;

    virtual void addToHtmlDoc(CParamHTML & doc) const;
    virtual void addValsToHTMLDoc(CParamHTML & doc) const = 0;

public:
    CParam(const char * szName, const char * szDescription, CParamClass * pParent_in);

    void init(config * pConfig)
    {
        if(IS_DEBUG) CHECK(pParent, "This class is not the root of a param tree so should not be initialised directly");
        initChildren(pConfig);
    }

    void name(std::string & addTo) const;

    bool isInit() const
    {
        return (initialised==eInit);
    }

    virtual void printParams(eParamStatus initd) const;
	virtual void printParam(std::ostream & os = std::cout) const = 0;

    virtual void printUseSummary() const;

    //Print out the current config file, with annotations
    void printCfgFile() const;
    void printHtmlDoc(std::ostream & os) const;
};

std::ostream& operator<<(std::ostream& s, const CParam & X);

class CParamClass : public CParam //A param class nust also notify a parent of its existnce
{
    typedef std::set<CParam *, CParam::CParamSort > TChildren;
    TChildren children;
protected:
    virtual void initChildren(config * pConfig);
    virtual void initValue(const char * strName);
    virtual void printAll(ePrintAllSettings eMode) const;
    virtual void printAllCheck(ePrintAllSettings eMode) const;
    virtual void addToHtmlDoc(CParamHTML & doc) const;
    virtual void addValsToHTMLDoc(CParamHTML & doc) const { doc.addDef(std::string()); };
public:
    CParamClass(const char * szName, const char * szDescription, CParamClass * pParent) : CParam(szName, szDescription, pParent)
    {
    }

    void notify(CParam * pChild)
    {
        children.insert(pChild);
    }

    virtual void printParams(eParamStatus initd) const;
    virtual void printParam(std::ostream & os = std::cout) const;
    virtual void printUseSummary() const;
};

class CStringParam : public CParam
{
    const std::string defaultVal; 
    std::string val;
protected:
    virtual void initValue(const char * szVal);
    virtual void printAll(ePrintAllSettings eMode) const;
    virtual void addValsToHTMLDoc(CParamHTML & doc) const { doc.addDef(defaultVal); };
public:
    void operator=(const std::string & newVal) //Should be safe because of initialiser counting
    {
        const bool bChanged = (val!=newVal);
        val=newVal;

        //std::string accumulateName;
        //name(accumulateName);
        if(bChanged)
        {
            std::cout << "Set "; printParam(cout);
        }

        initialised = eInit;

        if(IS_DEBUG) CHECK(!isInit(), "Param initialisation failed");
    }

    CStringParam(const char * szName, const char * defaultVal, const char * szDescription, CParamClass * pParent) : CParam(szName, szDescription, pParent),
        defaultVal(defaultVal), val(defaultVal)
    {
    }

    /*NOT virtual*/operator const std::string & () const { setUsed(); return val; };//Flag as used?? Debug only??

    virtual void printParam(std::ostream & os = std::cout) const;

    const char * asSz() const { setUsed(); return val.c_str(); }
};

template<typename T>
class CNumParam : public CParam
{
    const T minVal, maxVal, defaultVal; 
    T val;
protected:
    virtual void initValue(const char * szVal);
    virtual void printAll(ePrintAllSettings eMode) const;
    virtual void addValsToHTMLDoc(CParamHTML & doc) const;
public:
    void operator=(T newVal) //Should be safe because of initialiser counting
    {
        const bool bChanged = (val!=newVal);
        val=newVal;

        //std::string accumulateName;
        //name(accumulateName);
        if(bChanged) {
        std::cout << "Set "; printParam(); }

        if(newVal < minVal)
            std::cout << "Value too low (min=" << minVal << ")" << std::endl;
        else if(newVal > maxVal)
            std::cout << "Value too high (max=" << maxVal << ")" << std::endl;
        else
        {
            initialised = eInit;
        }
        if(IS_DEBUG) CHECK(!isInit(), "Param initialisation failed");
    }

    CNumParam(const char * szName, const T minVal, const T maxVal, const T defaultVal, const char * szDescription, CParamClass * pParent) : CParam(szName, szDescription, pParent),
        minVal(minVal), maxVal(maxVal), defaultVal(defaultVal), val(defaultVal)
    {
        if(val<minVal || val >maxVal)
        {
            cout << *this << endl;
            THROW("Parameter default value is not within min/max allowed values");
        }
    }

    //NOT virtual
    operator const T () const
    {
        setUsed(); //Flag as used. Debug only??
        if(val<minVal || val >maxVal)
        {
            cout << *this << endl;
            THROW("CNumParam::operator const T(): uninitialised parameter accessed");
        }
        return val;
    };

    virtual void printParam(std::ostream & os = std::cout) const;
};

template<>
void CNumParam<bool>::addValsToHTMLDoc(CParamHTML & doc) const;

template<typename T>
void CNumParam<T>::addValsToHTMLDoc(CParamHTML & doc) const 
{
    doc.addMMD(minVal, maxVal, defaultVal); 
}

template<>
void CNumParam<bool>::printParam(std::ostream & os) const;

template<typename T>
void CNumParam<T>::printParam(std::ostream & os) const
{
    std::string accumulateName;
    name(accumulateName);
    os << accumulateName << "=" << val << std::endl;
}

template<> void CNumParam<bool>::printAll(ePrintAllSettings eMode) const;

template<typename T>
void CNumParam<T>::printAll(CParamClass::ePrintAllSettings) const
{
    std::string accumulateName;
    name(accumulateName);
    std::cout << accumulateName << "=" << val << " # min=" << minVal << ", max=" << maxVal << ", default=" << defaultVal << std::endl;
}

//Parameter that MUST NOT be set from a config file--it will be set dynamically (e.g. image size)
template<typename T>
class CNumParamDerived : public CNumParam<T>
{
protected:
    virtual void initValue(const char *) { THROW( "Value of this parameter is derived and should not be set in config file"); }
    virtual void printAll(CParamClass::ePrintAllSettings) const {}; //Print nothing because should not be set in file
public:
    CNumParamDerived(const char * szName, const T minVal, const T maxVal, const T defaultVal, const char * szDescription, CParamClass * pParent) : CNumParam<T>(szName, minVal, maxVal, defaultVal, szDescription /*"This parameter must be set in code before use (e.g. parameter derived from other parameters)"*/, pParent) {}

    operator const T () const { if(IS_DEBUG) CHECK(!CNumParamDerived<T>::isInit(), "This parameter must be initialised in code (not in config file) before use"); return CNumParam<T>::operator const T(); };//Flag as used?? Debug only??

    using CNumParam<T>::operator=; //Doesn't inherit by default, cos like a c'tor

    virtual void initChildren(config *) {}; //Do not init me from config file (or print warning)
};

#define templateCEnumParam template<typename TEnum, class TfnMembers>
#define TEnumParam CEnumParam<TEnum, TfnMembers>

templateCEnumParam
class CEnumParam : public CParam
{
    const TEnum defaultVal; 
    TEnum val;
    TfnMembers fnMembers;
protected:
    virtual void initValue(const char * szVal);
    virtual void printAll(ePrintAllSettings eMode) const;
    virtual void addValsToHTMLDoc(CParamHTML & doc) const 
    { 
        std::string vals(fnMembers(defaultVal));

        vals += " (";
        for(int i=0;; i++)
        {
            const char * szMember = fnMembers(i);
            if(!szMember) break;

            if(i>0) vals += ", ";
            vals += szMember;
        }
        vals += ")";

        doc.addDef(vals); 
    }
public:
    void operator=(TEnum newVal) //Should be safe because of initialiser counting
    {
        const bool bChanged = (val!=newVal);

        val=newVal;

        //std::string accumulateName;
        //name(accumulateName);
        if(bChanged)
        {
            std::cout << "Set "; printParam();
        }

        initialised = eInit;

        if(IS_DEBUG) CHECK(!isInit(), "Param initialisation failed");
    }

    CEnumParam(const char * szName, const TEnum defaultVal, const char * szDescription, CParamClass * pParent) : CParam(szName, szDescription, pParent),
        defaultVal(defaultVal), val(defaultVal)
    {
    }

    //NOT virtual
    operator const TEnum () const { setUsed(); return val; };
    virtual void printParam(std::ostream & os = std::cout) const;

    static int NUM_OPTIONS()
    {
        TfnMembers fnMembers;
        int nMembers=0;
        while(fnMembers(nMembers)) nMembers++;
        return nMembers;
    }
};

templateCEnumParam
void TEnumParam::printParam(std::ostream & os) const
{
    std::string accumulateName;
    name(accumulateName);
    os << accumulateName << "=" << fnMembers(val) << std::endl;
}

templateCEnumParam
void TEnumParam::printAll(ePrintAllSettings) const
{
    std::string accumulateName;
    name(accumulateName);
    std::cout << accumulateName << "=" << fnMembers(val) << " # Options: ";
    
    for(int i=0;; i++)
    {
        const char * szMember = fnMembers(i);
        if(!szMember) break;

        if(i>0) std::cout << ", ";
        std::cout << szMember;
    }
    
    std::cout << std::endl;
}

templateCEnumParam
void TEnumParam::initValue(const char * szVal)
{
    if(IS_DEBUG) CHECK(!szVal, "Null not allowed");

    int i=0;
    for(;; i++)
    {
        const char * szMember = fnMembers(i);
        if(!szMember) break;

#ifdef __GNUC__
        if(strcasecmp(szMember, szVal) == 0)
#else
        if(_stricmp(szMember, szVal) == 0)
#endif
        {
            *this = (TEnum)i;
            return;
        }
    }
    const int nValInt = atoi(szVal);
    if(nValInt == 0 && szVal[0] != '0')
    {
        printAll(ePrintAll);
        THROW("String to set enum parameter not a valid option");
    }


    if(IS_DEBUG) CHECK(nValInt < 0 || nValInt >= i, "Integer value of enum parameter OOB");

    *this = (TEnum)nValInt;

    std::cout << "Warning: initialised enum param from integer" << std::endl;
}

#define MAKEENUMPARAM1(paramName, e1) _MAKEENUMPARAM(paramName, 1, e1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM2(paramName, e1,e2) _MAKEENUMPARAM(paramName, 2, e1,e2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM3(paramName, e1,e2,e3) _MAKEENUMPARAM(paramName, 3, e1,e2,e3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM4(paramName, e1,e2,e3,e4) _MAKEENUMPARAM(paramName, 4, e1,e2,e3,e4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM5(paramName, e1,e2,e3,e4,e5) _MAKEENUMPARAM(paramName, 5, e1,e2,e3,e4,e5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM6(paramName, e1,e2,e3,e4,e5,e6) _MAKEENUMPARAM(paramName, 6, e1,e2,e3,e4,e5,e6,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM7(paramName, e1,e2,e3,e4,e5,e6,e7) _MAKEENUMPARAM(paramName, 7, e1,e2,e3,e4,e5,e6,e7,0,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM8(paramName, e1,e2,e3,e4,e5,e6,e7,e8) _MAKEENUMPARAM(paramName, 8, e1,e2,e3,e4,e5,e6,e7,e8,0,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM9(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9) _MAKEENUMPARAM(paramName, 9, e1,e2,e3,e4,e5,e6,e7,e8,e9,0,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM10(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10) _MAKEENUMPARAM(paramName, 10, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,0,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM11(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11) _MAKEENUMPARAM(paramName, 11, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,0,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM12(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12) _MAKEENUMPARAM(paramName, 12, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,0,0,0,0,0,0,0,0)
#define MAKEENUMPARAM13(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13) _MAKEENUMPARAM(paramName, 13, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,0,0,0,0,0,0,0)
#define MAKEENUMPARAM14(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14) _MAKEENUMPARAM(paramName, 14, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,0,0,0,0,0,0)
#define MAKEENUMPARAM15(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15) _MAKEENUMPARAM(paramName, 15, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,0,0,0,0,0)
#define MAKEENUMPARAM16(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16) _MAKEENUMPARAM(paramName, 16, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,0,0,0,0)
#define MAKEENUMPARAM17(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17) _MAKEENUMPARAM(paramName, 17, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,0,0,0)
#define MAKEENUMPARAM18(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18) _MAKEENUMPARAM(paramName, 18, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,0,0)
#define MAKEENUMPARAM19(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19) _MAKEENUMPARAM(paramName, 19, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,0)
#define MAKEENUMPARAM20(paramName, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) _MAKEENUMPARAM(paramName, 20, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20)

#define _MAKE_STRINGS1(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,0}
#define _MAKE_STRINGS2(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,0}
#define _MAKE_STRINGS3(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,0}
#define _MAKE_STRINGS4(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,0}
#define _MAKE_STRINGS5(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,0}
#define _MAKE_STRINGS6(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,0}
#define _MAKE_STRINGS7(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,0}
#define _MAKE_STRINGS8(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,0}
#define _MAKE_STRINGS9(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,0}
#define _MAKE_STRINGS10(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,0}
#define _MAKE_STRINGS11(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,0}
#define _MAKE_STRINGS12(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,0}
#define _MAKE_STRINGS13(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,0}
#define _MAKE_STRINGS14(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,0}
#define _MAKE_STRINGS15(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,#e15,0}
#define _MAKE_STRINGS16(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,#e15,#e16,0}
#define _MAKE_STRINGS17(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,#e15,#e16,#e17,0}
#define _MAKE_STRINGS18(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,#e15,#e16,#e17,#e18,0}
#define _MAKE_STRINGS19(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,#e15,#e16,#e17,#e18,#e19,0}
#define _MAKE_STRINGS20(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {#e1,#e2,#e3,#e4,#e5,#e6,#e7,#e8,#e9,#e10,#e11,#e12,#e13,#e14,#e15,#e16,#e17,#e18,#e19,#e20,0}

#define _MAKE_ENUM1(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1}
#define _MAKE_ENUM2(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2}
#define _MAKE_ENUM3(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3}
#define _MAKE_ENUM4(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4}
#define _MAKE_ENUM5(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5}
#define _MAKE_ENUM6(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6}
#define _MAKE_ENUM7(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7}
#define _MAKE_ENUM8(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8}
#define _MAKE_ENUM9(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9}
#define _MAKE_ENUM10(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10}
#define _MAKE_ENUM11(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11}
#define _MAKE_ENUM12(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12}
#define _MAKE_ENUM13(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13}
#define _MAKE_ENUM14(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14}
#define _MAKE_ENUM15(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14,e##e15}
#define _MAKE_ENUM16(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14,e##e15,e##e16}
#define _MAKE_ENUM17(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14,e##e15,e##e16,e##e17}
#define _MAKE_ENUM18(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14,e##e15,e##e16,e##e17,e##e18}
#define _MAKE_ENUM19(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14,e##e15,e##e16,e##e17,e##e18,e##e19}
#define _MAKE_ENUM20(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) {e##e1,e##e2,e##e3,e##e4,e##e5,e##e6,e##e7,e##e8,e##e9,e##e10,e##e11,e##e12,e##e13,e##e14,e##e15,e##e16,e##e17,e##e18,e##e19,e##e20}

#define _MAKEENUMPARAM(paramName, numParams, e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20) \
    private: \
class CGet ## paramName ## Strings { \
    public: \
        const char * operator()(int idx) const \
        { \
            const char * asz[] = _MAKE_STRINGS ## numParams(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20); \
            return asz[idx]; \
        } \
    }; \
    public: \
    enum e ## paramName _MAKE_ENUM ## numParams(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20); \
    CEnumParam<e ## paramName, CGet ## paramName ## Strings > paramName;


#define PARAMCLASS_int(name, nameStr) \
class C ## name ## Params : public CParamClass\
{\
public:\
    C ## name ## Params(const char * szDescription, CParamClass * pParent) : CParamClass(nameStr, szDescription, pParent)

#define PARAMCLASS(name) PARAMCLASS_int(name, #name)
#define WRAPPERPARAMCLASS(name) PARAMCLASS_int(name, 0)

#define PARAM(paramName, TMin, TMax, TDefault, description) , paramName(# paramName, TMin, TMax, TDefault, description, this )
#define PARAMB(paramName, TDefault, description) PARAM(paramName, false, true, TDefault, description)
#define PARAME(paramName, TDefault, description) , paramName(# paramName, e ## TDefault, description, this)
#define PARAMSTR(paramName, description) , paramName(# paramName, "", description, this)
#define PARAMSTR2(paramName, szDefault, description) , paramName(# paramName, szDefault, description, this)

#define CHILDCLASS(name, description) , name( description, this ) 

#define MAKECHILDCLASS(name) C ## name ## Params name;

#define EXPAND_ALL(x) x;

#endif /* PARAM_H_ */
