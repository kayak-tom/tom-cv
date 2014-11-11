#pragma once
#ifndef STATS_H
#define STATS_H

#include "dynArray.h"
#include "convert.h"
#include <algorithm>
#include <iostream>

//INEFFICIENT concise way to compute medians, mean, etc.
class CStats
{
    CDynArray<double> data;
	std::string szName;
public:
	CStats() {}
	CStats(const std::string name) : szName(name) {}
	void setName(const char * szName_in) { szName = szName_in; }

    void reset() { data.clear(); }
    void add(double d) { data.push_back(d); };
    double mean() const
    {
        if(data.size()==0)
            return 0;

        double dMean = 0;
        for(CDynArray<double>::const_iterator pd = data.begin(); pd != data.end(); pd++)
            dMean += *pd;
        return dMean / data.size();
    }

    double variance() const
    {
        if(data.size()==0)
            return 0;

        double dMean = mean(), dVar = 0;
        for(CDynArray<double>::const_iterator pd = data.begin(); pd != data.end(); pd++)
            dVar += sqr(dMean - *pd);

        return dVar / data.size();
    }

    double sd() const { return sqrt(variance()); }

    double median()
    {
        if(data.size()==0)
            return 0;

        std::sort(data.begin(), data.end());

        if(data.size() %2 == 1)
            return data[data.size()/2];
        else
            return 0.5*(data[(data.size()/2) - 1] + data[data.size()/2]);
    }

    enum eConfBounds { eSymConf95, eSymConf99, eSymConf9995 };
    double confBound(const eConfBounds confBound) const
    {
        double t_val = HUGE;
        bool bFewPoints = (data.size()<100);//Use values for 20
        switch(confBound)
        {
        case eSymConf95:
            t_val = bFewPoints ? 1.724718 : 1.644854;
            break;
        case eSymConf99:
            t_val = bFewPoints ? 2.52798 : 2.32635;
            break;
        case eSymConf9995:
            t_val = bFewPoints ? 3.8495 : 3.2905;
            break;
        }

        double dXPercentSymConfidenceBound = t_val*sd()/sqrt((double)(data.size()-1));

        return dXPercentSymConfidenceBound;
    }

    static void writeTSVheader(const char * dataName, std::ostream & file)
    {
        file << dataName << "Name \tCount\tmean \t99% conf\ts.d.\t" << dataName << " median\t";
    }



    void writeTSVdata(std::ostream & file, bool bCount= false)
    {
	if(szName.size()) file << szName << '\t';
        if(bCount) file << data.size() << '\t';
        file << mean() << '\t'
             << confBound(eSymConf99) << '\t'
             << sd() << '\t'
             << median() << '\t';
    }

    void printData() const
    {
        for(CDynArray<double>::const_iterator pd = data.begin(); pd != data.end(); pd++)
        {
            std::cout << *pd << ", ";
        }
        std::cout << '\n';
    }        
        
    void pp(const char * seperator = ", ", std::ostream & s = std::cout) const
    {
        data.pp(seperator, s);
    }
    
    void cropTop(double dProportion)
    {
        if(data.size()==0)
            return;

        std::sort(data.begin(), data.end());
        
        data.resize((int)((1-dProportion)*(double)data.size()));
    }
};

class CRTErrorStats {
    CStats statsR, statsT;
    int nFail, nLargeErrors, nTotalRuns;
    std::ostream & results;
public:
    const std::string strName;

    CRTErrorStats(std::ostream & results, const char * szName) : nFail(0), nLargeErrors(0), nTotalRuns(0), results(results), strName(szName) {
    }

    void addFail() {
        nFail++;
        nTotalRuns++;
    }

    void addRT(double dErrInR, double dErrInT) {
        nTotalRuns++;
        double dSuccessThresh = 0.25;
        if (dErrInR > dSuccessThresh)
            nLargeErrors++;
        //else {
        statsR.add(dErrInR);
        statsT.add(dErrInT);
        //}
    }

    ~CRTErrorStats() {
        pp(results);
    }

    void pp(std::ostream & os) {
        if (nTotalRuns <= 1)
            return;

        double dSuccessRate = 1 - ((double) nFail / (double) nTotalRuns);
        double dSolutionsWithErrors = 1 - ((double) (nLargeErrors + nFail) / (double) nTotalRuns);
        std::ostringstream ss;

        ss << "\"" << strName << "\"\tRuns:" << nTotalRuns << '\t' << statsR.mean() << '\t' << statsR.median() << '\t' << statsT.mean() << '\t' << statsT.median() << "\tSuccesses" << dSuccessRate << "\tUsefulSuccesses" << dSolutionsWithErrors << endl;
        cout << ss.str();
        os << ss.str();
    }
};

#endif
