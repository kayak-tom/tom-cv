/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "params/config.h"
//#define BOOST_SPIRIT_USE_OLD_NAMESPACE
//#include "boost/spirit.hpp"
#include "boost/spirit/include/classic.hpp"
#include "boost/filesystem.hpp"
#include <vector>
#include <set>
#include "string.h"
#include "util/random.h"
#include "util/exception.h"
#include<boost/filesystem/operations.hpp>
#include "time/SpeedTest.h"

using namespace std;
using namespace boost::filesystem;
using namespace boost::spirit::classic;

double run(const bool bOutput, const char* szConfigFile, const char * szOtherConfig);
class paramSet : public config {
    double score, time;
public:
    paramSet(const char * szFileName) : config(szFileName), score(0), time(0) {
        get("SCORE", score);
        score *= 0.75; //only on load
        get("TIME", time);
        time *= 1.02;
    };
    paramSet() : config(true), score(0), time(0) {
    };
    void ageScore() //Drop score with age, to ensure params consistently good rather than coincidentally
    {
        score *= 0.95;
    }
    void saveToFolder(const char * szFolder) const {
        char szFullFN[500];
        sprintf_s(szFullFN, 500, "%s/s=%.2f-t=%.2f.cfg", szFolder, score, time);
        save(szFullFN);
    }
    void save(const char * szCfg) const {
        ofstream file(szCfg);
        for (TSzMap::const_iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++) {
            file << pMI->first << '=' << pMI->second << endl;
        }
        file.close();
    };
    void testSet(const char * szOtherConfig) {
        const char * szCfg = "tune.cfg";
        save(szCfg);

        CStopWatch timer;
        timer.startTimer();
        score = run(false, szCfg, szOtherConfig);
        timer.stopTimer();
        time = fabs(timer.getElapsedTime());

        char szScore[50];
        sprintf_s(szScore, 50, "%f", score);
        addParam(cfg_strndup("SCORE", 10), cfg_strndup(szScore, 50));
        sprintf_s(szScore, 50, "%f", time);
        addParam(cfg_strndup("TIME", 10), cfg_strndup(szScore, 50));

        remove(szCfg);
    };
    void addParam(const char * szName, char * szVal) {
        //delete if already exists
        TSzMap::const_iterator pName = cfgNameVal.find(szName);
        if (pName != cfgNameVal.end()) {
            cfgNameVal.erase(szName);
            free((void*) pName->first);
            free(pName->second);
        }
        insert(szName, szVal);
        pName = cfgNameVal.find(szName);
        CHECK(pName == cfgNameVal.end(), "Insert failed")

    }
    void saveToSS(const char * szFolder) const {
        char szFullFN[500];
        sprintf_s(szFullFN, 500, "%s/params.tsv", szFolder);
        if (!exists(szFullFN)) {
            ofstream file(szFullFN);
            for (TSzMap::const_iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++) {
                file << pMI->first << '\t';
            }
            file << endl;

            file.close();
        }

        ofstream file(szFullFN, ios_base::app);
        char szParamOnly[500];
        for (TSzMap::const_iterator pMI = cfgNameVal.begin(); pMI != cfgNameVal.end(); pMI++) {
            sprintf_s(szParamOnly, 500, "%s", pMI->second);
            char * pcHash = strchr(szParamOnly, '#');
            if (pcHash) *pcHash = 0;
            file << szParamOnly << '\t';
        }
        file << endl;

        file.close();
    };
    double getScore() const {
        return score;
    };
    double getRank() const {
        //if(score <= 20)
        //return score ; //Boost scores for param sets that are ok anyway
        return score / (100. + 0.05 * time);
    };
};
class param {
    enum eType {
        eSz, eDouble, eInt, eDiscrete, eIntFixed, eDoubleFixed
    };
    eType type;
    enum eMethod {
        eChoose1, eChoose2, eAverage
    };
    char * name, * comment, * szValue;

    static const double maxDeviation, mutateRate;
    double dMin, dMax, dNum;
    int nMin, nMax, nNum;
    vector<int> anEnumVals;
    bool bParsedOk;

public:
    param(const char * szLine) : type(eIntFixed), name(0), comment(0), szValue(0), dMin(0), dMax(0), dNum(0), nMin(0), nMax(0), nNum(0), bParsedOk(false) {
        if (parse(szLine, str_p("TIME=")).hit || parse(szLine, str_p("SCORE=")).hit) return;

        //decide what type+setup ranges
        string strName, strVal, strComment;

        rule<const char*> doubleRangeParser = ch_p('[') >> real_p[assign_a(dMin)] >> ch_p(',') >> real_p[assign_a(dMax)] >> ch_p(']')[assign_a(type, eDouble)];
        rule<const char*> intRangeParser = ch_p('[') >> int_p[assign_a(nMin)] >> ch_p(',') >> int_p[assign_a(nMax)] >> ch_p(']')[assign_a(type, eInt)];
        rule<const char*> intSetParser = ch_p('{') >> int_p[push_back_a(anEnumVals)] >> *(ch_p(',') >> int_p[push_back_a(anEnumVals)]) >> ch_p('}')[assign_a(type, eDiscrete)];
        rule<const char*> stringParser = ch_p('"') >> (+~ch_p('"'))[assign_a(strVal)] >> ch_p('"')[assign_a(type, eSz)];
        rule<const char*> intParser = int_p[assign_a(nNum)];
        rule<const char*> doubleParser = real_p[assign_a(dNum)];

        rule<const char*> lineParser = (+~ch_p('='))[assign_a(strName)] >> ch_p('=') >> (intRangeParser || doubleRangeParser || intSetParser || stringParser || intParser || doubleParser) >> *(~ch_p('#')) >> !(ch_p('#') >> (+anychar_p)[assign_a(strComment)]);

        bParsedOk = parse(szLine, lineParser).hit;

        if (bParsedOk) {
            if (type == eIntFixed && dNum != 0) type = eDoubleFixed;

            name = config::cfg_strndup(strName.c_str(), strName.length());
            comment = config::cfg_strndup(strComment.c_str(), strComment.length());
            if (type == eSz)
                szValue = config::cfg_strndup(strVal.c_str(), strVal.length());

            if(IS_DEBUG) CHECK(!name || strlen(name) == 0, "param: Name not initialised")
        } else
            bParsedOk = false;
    };
    ~param() {
        free(szValue);
        free(name);
        free(comment);
    }
    bool isGood() const {
        return bParsedOk;
    };
    const char * getName() const {
        return config::cfg_strndup(name, 1000);
    };
    char * getVal(paramSet & p1, paramSet &p2) const {
        char szVal[2000], szVal2[2000];

        if (type == eSz) {
            if (strchr(szValue, '"'))
                sprintf_s(szVal, 2000, "%s", szValue);
            else
                sprintf_s(szVal, 2000, "\"%s\"", szValue);
        } else if (type == eDoubleFixed) {
            sprintf_s(szVal, 2000, "%f", dNum);
        } else if (type == eIntFixed) {
            sprintf_s(szVal, 2000, "%d", nNum);
        } else if (type == eDouble) {
            double dp1 = 0, dp2 = 0, dNewVal = 0;
            p1.get(name, dp1);
            p2.get(name, dp2);

            dNewVal = mutate(combine(dp1, dp2));
            sprintf_s(szVal, 2000, "%f", dNewVal);
        } else if (type == eInt) {
            int np1 = 0, np2 = 0, nNewVal = 0;
            p1.get(name, np1);
            p2.get(name, np2);

            nNewVal = mutate(combine(np1, np2));
            sprintf_s(szVal, 2000, "%d", nNewVal);
        } else if (type == eDiscrete) {
            int np1 = 0, np2 = 0, nNewVal = 0;
            p1.get(name, np1);
            p2.get(name, np2);

            nNewVal = mutateEnum(combineEnum(np1, np2));
            sprintf_s(szVal, 2000, "%d", nNewVal);
        }

        if (strlen(comment) > 0)
            sprintf_s(szVal2, 2000, "%s # %s", szVal, comment);
        else
            sprintf_s(szVal2, 2000, "%s", szVal);

        //Check it exists, add to child param sets if not
        if (p1.getSz(name) == 0) {
            p1.addParam(config::cfg_strndup(name, 2000), config::cfg_strndup(szVal2, 2000));
        }
        if (p2.getSz(name) == 0) {
            p2.addParam(config::cfg_strndup(name, 2000), config::cfg_strndup(szVal2, 2000));
        }

        return config::cfg_strndup(szVal2, 2000);
    };
    double inRange(double d1) const {
        if (d1 < dMin) return dMin;
        if (d1 > dMax) return dMax;
        return d1;
    };
    double getRandScale() const {
        return CRandom::Uniform(1 - maxDeviation, 1 + maxDeviation);
    };
    double mutate(double d1) const {
        if (CRandom::Uniform(0.0, 1.0) > mutateRate) return d1;

        double dScale = getRandScale();
        return inRange((d1 >= 0 ? 1 : -1) * (0.0 + fabs(d1)) * dScale); //Adding a number tends to grow insensitive parameters
    };
    int mutate(int n1) const {
        if (CRandom::Uniform(0.0, 1.0) > mutateRate) return n1;

        int nDev = 1 + (int) fabs((n1 * maxDeviation));
        return inRange(CRandom::Uniform(n1 - nDev, n1 + nDev + 1));
    };
    int mutateEnum(int e) const {
        if (CRandom::Uniform(0.0, 1.0) > mutateRate) return e;
        return anEnumVals[CRandom::Uniform((int) anEnumVals.size())];
    };
    double inRange(int d1) const {
        if (d1 < nMin) return nMin;
        if (d1 > nMax) return nMax;
        return d1;
    };
    double combine(double d1, double d2) const {
        int rand = CRandom::Uniform(3);
        eMethod method = (eMethod) (rand);

        switch (method) {
            case eChoose1:
                return inRange(d1);
            case eChoose2:
                return inRange(d2);
            case eAverage:
                return inRange(0.5 * (d1 + d2));
        }
        THROW("Combine: value not returned");
    };
    int combineEnum(int n1, int n2) const {
        int rand = CRandom::Uniform(2);
        eMethod method = (eMethod) (rand);

        switch (method) {
            case eChoose1:
                return (n1);
            case eChoose2:
                return (n2);
            case eAverage:
                THROW("Combine: value not handled");
        }
        THROW("Combine: value not returned");
    };
    int combine(int n1, int n2) const {
        int rand = CRandom::Uniform(3);
        eMethod method = (eMethod) (rand);

        switch (method) {
            case eChoose1:
                return inRange(n1);
            case eChoose2:
                return inRange(n2);
            case eAverage:
                return inRange((n1 + n2) / 2);
        }
        THROW("Combine: value not returned");
    };

};
struct paramSetCmp {
    bool operator()(const paramSet* s1, const paramSet * s2) const {
        return s1->getRank() > s2->getRank();
    };
};

typedef multiset<paramSet *, paramSetCmp> TParamPopulation;

typedef vector<paramSet *> TChildParamSets;
class CChildParamSets : public TChildParamSets {
public:
    void test(const char * szOtherParams) {
        for (iterator ppChild = begin(); ppChild < end(); ppChild++)
            (*ppChild)->testSet(szOtherParams);
    }
    void saveToSS(const char * szFolder) const {
        for (const_iterator ppChild = begin(); ppChild < end(); ppChild++)
            (*ppChild)->saveToSS(szFolder);
    }
};
class CParamPopulation : public TParamPopulation {
public:
    ~CParamPopulation() {
        for (iterator ppChild = begin(); ppChild != end(); ppChild++)
            delete *ppChild;
    }
    void addChildren(const CChildParamSets & children) {
        for (CChildParamSets::const_iterator ppChild = children.begin(); ppChild < children.end(); ppChild++)
            insert(*ppChild);
    }

    //Delete weakest in population
    void killWeakest(const int nPopulationSizeTarget) {
        iterator pWorst = end();
        int numToDelete = (int) size() - nPopulationSizeTarget;
        for (int nDeleted = 0; nDeleted < numToDelete; nDeleted++) {
            pWorst--;

            cout << "Deleting child with score " << (*pWorst)->getScore() << endl;
            delete *pWorst;
        }
        erase(pWorst, end());

    }
    void save(const char * szFolder) const {
        //First delete the last lot
        directory_iterator end_itr;
        for (directory_iterator itr(szFolder); itr != end_itr; ++itr) {
            if (is_regular(itr->status())) {
                const path p = itr->path();
                string strImFilename = p.leaf().string(); //Remove .string() if this causes boost errors
                if (strcasestr(strImFilename.c_str(), ".cfg")) {
                    remove(itr->path());
                }
            }
        }
        for (const_iterator pSel = begin(); pSel != end(); pSel++)
            (*pSel)->ageScore();
        for (const_iterator pSel = begin(); pSel != end(); pSel++)
            (*pSel)->saveToFolder(szFolder);
    }
    CParamPopulation(const char * szFolder) {
        //Read all param sets in folder to init population
        directory_iterator end_itr;
        for (directory_iterator itr(szFolder); itr != end_itr; ++itr) {
            if (is_regular(itr->status())) {
                const path p = itr->path();
                string strImFilename = p.leaf().string();
                if (strcasestr(strImFilename.c_str(), ".cfg")) {
                    char szParamPath[2000];
                    sprintf_s(szParamPath, 2000, "%s/%s", szFolder, strImFilename.c_str());
                    insert(new paramSet(szParamPath));
                }
            }
        }

    }
    paramSet *random() {
        int nSelect = CRandom::Uniform((int) size());

        iterator pSel = begin();
        for (int i = 0; i < nSelect; i++)
            pSel++;
        if(IS_DEBUG) CHECK(pSel == end(), "random: Error")
        return *pSel;
    };

};
class paramSettings {
    vector<const param *> params;
public:
    paramSettings(const char * szFolderName) {
        //Parse the param info file:
        char szFileName[200];
        sprintf_s(szFileName, 200, "%s/paramInfo.dat", szFolderName);

        ifstream fileSettings(szFileName);
        CHECK(!fileSettings.is_open(), "paramSettings: Settings file not found")

        while (!fileSettings.eof()) {
            char szLine[1024] = "";
            fileSettings.getline(szLine, 1024);
            if (strlen(szLine)) {
                const param * pParam = new param(szLine);
                if (pParam->isGood())
                    params.push_back(pParam);
                else {
                    cout << szLine << " not parsed\n";
                    delete pParam;
                }
            }
        }
        fileSettings.close();
        CHECK(params.size() == 0, "paramSettings: No settings loaded")
    };
    ~paramSettings() {
        for (vector<const param *>::iterator ppParam = params.begin(); ppParam < params.end(); ppParam++)
            delete *ppParam;
    }
    paramSet * mate(paramSet & p1, paramSet &p2) const {
        paramSet * child = new paramSet();
        for (vector<const param *>::const_iterator ppParam = params.begin(); ppParam < params.end(); ppParam++) {
            const param * pParam = *ppParam;
            const char * szName = pParam->getName();
            char * szVal = pParam->getVal(p1, p2);
            child->addParam(szName, szVal);
            //cout << "adding " << szName << '=' << szVal << ' '; child->dumpNames();
        }

        return child;
    };
    void breed(CChildParamSets & children, CParamPopulation & population, const int nChildren) const {
        if(IS_DEBUG) CHECK(population.size() == 0, "Tune: No population--probably cfg file missing from this folder")
        for (int i = 0; i < nChildren; i++) {
            paramSet * pParent1 = population.random();
            paramSet * pParent2 = population.random();
            children.push_back(mate(*pParent1, *pParent2));
        }
    };

};

const double param::maxDeviation = 0.35;
const double param::mutateRate = 0.2;
int tune(const char * szFolder, const char * szOtherParams) {
    int nPopulation = 8, nChildren = 4;

    paramSettings myParamSettings(szFolder);

    CParamPopulation population(szFolder);
    for (int i = 0; i < 300000; i++) {
        CChildParamSets children;
        myParamSettings.breed(children, population, nChildren);

        children.test(szOtherParams);
        children.saveToSS(szFolder);

        //Insert children into population
        population.addChildren(children);

        //Delete weakest in population
        population.killWeakest(nPopulation);

        population.save(szFolder);

    }
    return 0;
}
