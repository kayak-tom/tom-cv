/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * param.cpp
 *
 *  Created on: 15/10/2009
 *      Author: tom 
 */

#include "param.h" 
#include "config.h"

using namespace std;

template<> void CNumParam<bool>::printAll(ePrintAllSettings eMode) const { printParam(); };

void CParam::printAllCheck(ePrintAllSettings eMode) const
{
    if(printMe(eMode))
        printAll(eMode);
}

void CParamClass::printAllCheck(ePrintAllSettings eMode) const
{
    printAll(eMode);
}

void CParam::addToHtmlDoc(CParamHTML & doc) const
{
    std::string accumulateName;
    name(accumulateName);
    doc.addName(accumulateName);
    addValsToHTMLDoc(doc);
    doc.addExplanation(DESCRIPTION);
}

template<>
void CNumParam<bool>::addValsToHTMLDoc(CParamHTML & doc) const
{
    doc.addDef(defaultVal);
}

void CParamClass::addToHtmlDoc(CParamHTML & doc) const
{
    std::string accumulateNameBold("<b>");
    {
        std::string accumulateName;
        name(accumulateName);
        accumulateNameBold += accumulateName;
    }
    accumulateNameBold.append("</b>");
    doc.addName(accumulateNameBold);
    addValsToHTMLDoc(doc);

    std::string paramDesc("<b>PARAM CLASS</b> ");
    paramDesc.append(DESCRIPTION);

    doc.addExplanation(paramDesc);

    for(TChildren::const_iterator ppChild = children.begin(); ppChild != children.end(); ppChild++)
    {
        (*ppChild)->addToHtmlDoc(doc);
    }
}

void CParam::printHtmlDoc(ostream & os) const
{
    CParamHTML doc(os);
    addToHtmlDoc(doc);
}


void CParam::printCfgFile() const
{
    std::cout << "Config file used:\n";
    printAll(ePrintInitUsed);
    printAll(ePrintInitUnused);
    printAll(ePrintUsed);
    printAll(ePrintUnused);

    std::cout << std::endl;
}

bool CParam::printMe(ePrintAllSettings eMode) const
{
    switch(eMode)
    {
    case ePrintAll:
        return true;
    case ePrintInitUsed:
        return isInit() && isUsed();
    case ePrintInitUnused:
        return isInit() && !isUsed();
    case ePrintUsed:
        return !isInit() && isUsed();
    case ePrintUnused:
        return !isInit() && !isUsed();
    default:
        THROW("Enum option not handled")
    }
}


CParam::CParam(const char * szName, const char * szDescription, CParamClass * pParent_in) : NAME(szName ? szName : ""), DESCRIPTION(szDescription ? szDescription : ""), pParent(pParent_in), initialised(eUninit), used(eUnused)
{
    if(IS_DEBUG) CHECK(pParent && (!szName || strlen(szName)==0), "CParam::initName: Parent but no name given");

    if(pParent_in)
    {
        pParent_in->notify(this);
        //cout << "Parameter " << *this << " created\n";
    }
    else
    {
        if(IS_DEBUG) CHECK(!szName && (szDescription || pParent_in), "Missing parameter name");
        //if(!szName)
            //cout << "Expected parameter name, found 0..";
        //cout << "Parameter tree rooted at " << *this << " created\n";
    //CAN'T CALL YET    initChildren(pConfig); //Must ALSO call in derived-class c'tor as the class is currently a *CParam* (during construction/destruction)
    }
}
void CParam::printUseSummary() const
{
    if(used == eUnused && initialised == eInit)
    {
        cout << "Initialised but not used: ";
        printParam();
    }
    else if(used == eUsed && initialised == eUninit)
    {
        cout << "Uninitialised and default used: ";
        printParam();
    }
}
void CParamClass::printUseSummary() const
{
    for(TChildren::const_iterator ppChild = children.begin(); ppChild != children.end(); ppChild++)
    {
        (*ppChild)->printUseSummary();
    }
};


void CParam::name(std::string & addTo) const
{
    if(pParent)
    {
        pParent->name(addTo);

        if(addTo.length()>0)
            addTo += '.';
    }

    addTo += NAME;
}

void CParamClass::printParams(eParamStatus initd) const
{
    for(TChildren::const_iterator ppChild = children.begin(); ppChild != children.end(); ppChild++)
    {
        (*ppChild)->printParams(initd);
    }
}
void CParamClass::printAll(ePrintAllSettings eMode) const
{
    for(TChildren::const_iterator ppChild = children.begin(); ppChild != children.end(); ppChild++)
    {
        (*ppChild)->printAllCheck(eMode);
    }
}

void CParam::printParams(eParamStatus initd) const
{
    if(initialised == initd)
    {
        printParam();
    }
}

void CParamClass::printParam(std::ostream & os) const
{
    std::string accumulateName;
    name(accumulateName);
    os << accumulateName << " was initialised" << endl;
}

std::ostream& operator<<(std::ostream& s, const CParam & X)
{
    std::string accumulateName;
    X.name(accumulateName);
    s << accumulateName;
    return s;
}

void CParam::initChildren(config * pConfig)
{
    //Init this parameter value
    std::string accumulateName;
    name(accumulateName);

    const char * szVal = pConfig->getSz(accumulateName.c_str());

    //cout << "Read " << accumulateName << "=" << (szVal ? szVal : "*") << endl;

    if(szVal && strlen(szVal) > 0)
    {
        if(initialised == eInit)
        {
            cout << accumulateName << ": ";
            THROW( "CParam::initChildren: Parameter already initialised")
        }

        initValue(szVal);

        initialised = eInit;
    }
    else
    {
        initialised = eUninit;
    }
}

void CParamClass::initChildren(config * pConfig)
{
    //Init this parameter's children's values
    for(TChildren::iterator ppChild = children.begin(); ppChild != children.end(); ppChild++)
    {
        CParam * pChild = (*ppChild);
        pChild->initChildren(pConfig);
    }
    initialised = eInit; //Can already be init here (from base c'tor), but is a class so doesn't matter.
}

void CParamClass::initValue(const char * szVal)
{
    std::string accumulateName;
    name(accumulateName);
    cout << "Class name: " << accumulateName << ", value = " << szVal << '\n';
    THROW( "Error: Class name in config file")
}

template<> void CNumParam<double>::initValue(const char * szVal)
{
    *this = atof(szVal);
}

template<> void CNumParam<float>::initValue(const char * szVal)
{
    *this = (float)atof(szVal);
}

template<> void CNumParam<int>::initValue(const char * szVal)
{
    *this = atoi(szVal);
}

template<> void CNumParam<bool>::initValue(const char * szVal)
{
    printParam();
    //cout << "Setting val " << szVal << endl;

    std::string strVal(szVal);
    for(int i=0; i<(int)strVal.size();++i)
       strVal[i] = tolower(strVal[i]); //Todo: function

    if(strVal.compare("true") == 0)
        *this = true;
    else if(strVal.compare("on") == 0)
        *this = true;
    else if(strVal.compare("false") == 0)
        *this = false;
    else if(strVal.compare("off") == 0)
        *this = false;
    else if(strVal.compare("t") == 0)
        *this = true;
    else if(strVal.compare("f") == 0)
        *this = false;
    else
        *this = (atoi(szVal) != 0);
}

void CStringParam::printAll(ePrintAllSettings eMode) const { printParam(); }

void CStringParam::printParam(std::ostream & os) const
{
    std::string accumulateName;
    name(accumulateName);
    os << accumulateName << "=\"" << val << '\"' << std::endl;
}

template<>
void CNumParam<bool>::printParam(std::ostream & os) const
{
    std::string accumulateName;
    name(accumulateName);
    os << accumulateName << "=" << (val ? "true" : "false") << std::endl;
}

void CStringParam::initValue(const char * szVal)
{
    CHECK(!szVal, "Null not allowed");
    CHECK_P(szVal[0] != '\"', szVal, "String parameters in config file should be quoted");
    CHECK_P(strlen(szVal) < 2, szVal, "String parameters in config file should be quoted");

    char szValNoQuote[10000];
    //Get rid of quotes
#ifdef __GNUC__
    strncpy(szValNoQuote, szVal+1, strlen(szVal) - 1);
#else
    strncpy_s(szValNoQuote, 10000, szVal+1, strlen(szVal) - 1);
#endif
    //std::cout << szValNoQuote << std::endl;
    szValNoQuote[strlen(szVal) - 2] = 0;
    *this = szValNoQuote;
}


CParamHTML::CParamHTML(std::ostream & os) : os(os)
{
    os << "<html><head><title>Documentation for BoWSLAM, BaySAC, etc. parameters</title></head><body>\n";
    os << "<h2>Documentation for BoWSLAM, BaySAC, etc. parameters</h2>";
    os << "<h3>Generated by PARAMS.printHtmlDoc</h3>"; 
    os << "<p>Add parameters to a config file, each line should read \"NAME=val\", e.g. \"RANSAC.E_INLIER_THRESH_PX=3.5\" (without quotes)</p>";
    os << "<p>Pass path to config file as only command line argument: ./BoWSLAM /path/to/config.cfg</p>"; 
    os << "<p>Image source (directory or attached cam or video file) MUST be set in config file. Calibration data is usally needed too.</p>"; 
    os << "<table width='100%' border='1'>\n";
    os << "<tr><td><b>Name</b></td>";
    os << "<td><b>Min val</b></td>";
    os << "<td><b>Max val</b></td>";
    os << "<td><b>Default val</b></td>";
    os << "<td><b>Notes</b></td></tr>\n";
}

void CParamHTML::addName(const std::string & name)
{
    os << "<tr><td>" << name << "</td>";
}

void CParamHTML::addMMD(double min, double max, double def)
{
    os << "<td>" << min << "</td>";
    os << "<td>";

    const int BIG_INT = 1000000000;
    if(max < BIG_INT)
        os << max;
    else
        os << "[HUGE]";

    os << "</td>";
    os << "<td>";

    if(def < BIG_INT)
        os << def;
    else
        os << "[HUGE]";

    os << "</td>";
}

void CParamHTML::addDef(bool def)
{
    os << "<td>-</td>";
    os << "<td>-</td>";
    os << "<td>" << (def ? "true" : "false") << "</td>";
}

void CParamHTML::addDef(const std::string & def)
{
    os << "<td colspan=3>" << def << "</td>";
}

void CParamHTML::addExplanation(const std::string & explanation)
{
    os << "<td>" << explanation << "</td></tr>\n";
}

CParamHTML::~CParamHTML()
{
    os << "</table>\n";
    os << "</body></html>";
    os.flush();
}
