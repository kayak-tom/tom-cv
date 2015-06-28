#ifndef PR_H
#define PR_H

/**
 * @class CEvalPR
 * @brief Log some data, then print Precision/recall on destructor.
 */
#include <util/vectorUtil.h>
#include <util/pp.h>
//#include "optimiserStats.h"

enum eSuccessStatus { eFail, eSuccess }; //Repeated in lines.h

//Numbers that initialise to 0, useful for concise counting code in places
class CInt
{
    int i;
public:
    CInt() : i(0) {}
    operator int & () { return i; }
    operator const int () const { return i; }
};

class CDouble
{
    double d;
public:
    CDouble() : d(0) {}
    operator double & () { return d; }
    operator const double () const { return d; }
};


template<typename TNum>
class CEvalPR 
{
    TNum nTotalFound, nTotalGT, nTotalCorrect, nTotalCandidates;
    std::string label;
    boost::mutex mxPREvalVines;
public:
    CEvalPR(const std::string label) : nTotalFound(0), nTotalGT(0), nTotalCorrect(0), nTotalCandidates(0), label(label) {}

    void logMeasurements(const TNum nTotalFound_in, const TNum nTotalGT_in, const TNum nTotalCorrect_in, const TNum nTotalCandidates_in) {
        boost::mutex::scoped_lock lock(mxPREvalVines);
        nTotalCorrect += nTotalCorrect_in;
        nTotalGT += nTotalGT_in;
        nTotalFound += nTotalFound_in;
        nTotalCandidates += nTotalCandidates_in;

        const bool bVerbose = nTotalGT_in>0;

        const double dPrecision = nTotalCorrect_in/(double)nTotalFound_in;
        const double dRecall = nTotalCorrect_in/(double)nTotalGT_in;
        if(bVerbose)
        { 
            cout << label << " one frame proportion of 2D canes with ground truth: "; 
            if(nTotalFound_in>0) 
            {
                COUT(dPrecision);
            }
            else 
            {
                cout << "-" << endl;
            }
            
            COUT(dRecall);
        }
    }

    ~CEvalPR() {
        const bool bVerbose = nTotalGT>0;

        const double dPrecision = nTotalCorrect/(double)nTotalFound;
        const double dRecall = nTotalCorrect/(double)nTotalGT;
        if(bVerbose)
        {
            cout << label << ": all frames ground truth results: Precision=" << dPrecision << " Recall=" << dRecall << endl; 
        
            const double dNumIncorrectlyClassified = (double)((nTotalFound-nTotalCorrect) /*returned incorrect*/ + (nTotalGT-nTotalCorrect) /* not returned correct */);
            const double dEmpiricalClassificationError = dNumIncorrectlyClassified / nTotalCandidates;
            cout << label << ": " << TO_STRING(dEmpiricalClassificationError) << TO_STRING(nTotalCorrect) << TO_STRING(nTotalFound) << TO_STRING(nTotalGT) << TO_STRING(nTotalCandidates) << endl;
            
        }
    }
};

class CEvalNum
{
    std::string label;
    boost::mutex mxPREvalVines;
    std::vector<double> adNums;
public:
    CEvalNum(const std::string & label) : label(label) {}
    
    void log(const double dNum)
    {
        boost::mutex::scoped_lock lock(mxPREvalVines);
        CHECKBADNUM(dNum);
        adNums.push_back(dNum);
    }
    
    ~CEvalNum() 
    {
        if(adNums.size()==0)
            return;
        
        double dMean=0, dSD=0;
        meanSD(adNums, dMean, dSD);
        
        cout << label << ": mean=" << dMean << " SD=" << dSD << " Count=" << adNums.size() << endl;
    }
};

class CEvalFailure;

/*
 * Count success/fail at a point, print a summary at program end
 * */
class CEvalResults
{
    friend class CEvalFailure;
    
    static CEvalResults s_evalResults;
    
    class CSuccessCount
    {
        int anSuccessFail[2];
    public:
        CSuccessCount() 
        {
            anSuccessFail[0] = 0;
            anSuccessFail[1] = 0;
        }
        void log(const eSuccessStatus result)
        {
            anSuccessFail[result]++;
        }
        
        int total() const { return anSuccessFail[eSuccess] + anSuccessFail[eFail]; }
        
        std::string toString() const
        {
            return ::toString(anSuccessFail[eSuccess]) + "/" + ::toString(total()) + "=" + ::toString(100.0*anSuccessFail[eSuccess]/(double)total()) + "%";
        }
    };
    
    typedef std::map<std::string, CSuccessCount > TCounts;
    TCounts aSuccessFailCounts;
    void log_int(const std::string & label, const eSuccessStatus result)
    {
        aSuccessFailCounts[label].log(result);
    }
    
public:
    static void log(const std::string & label, const eSuccessStatus result)
    {
        s_evalResults.log_int(label, result);
    }
    
    static void log(const std::string & label, const bool bSuccessful)
    {
        s_evalResults.log_int(label, bSuccessful ? eSuccess : eFail);
    }
    
    ~CEvalResults() 
    {
        for(const TCounts::value_type & val : aSuccessFailCounts)
        {
            cout << val.first << ":\t" << val.second.toString() << endl;
        }
    }
};

/*
 * Count success/fail at various points in a function, print a summary at program end
 * */
class CEvalFailure
{
    static CEvalFailure s_evalFailures;
    
    class CFailureSummary
    {
        typedef std::map<std::string, CInt> TFailureSummary;
        TFailureSummary aFailureCounts;
        CEvalResults::CSuccessCount successCount;
    public:
        void log(const eSuccessStatus result, const std::string & reason)
        {
            successCount.log(result);
            aFailureCounts[reason]++;
        }
        
        std::string toString() const
        {
            std::ostringstream ss;
            
            ss << successCount.toString() << endl;
            
            const double dTotal = (double)successCount.total();
            
            for(const TFailureSummary::value_type & val : aFailureCounts)
            {
                ss << '\t' << val.first << ": " << val.second << " = " << std::setprecision(3) << (100.0*(int)val.second/dTotal) << "%" << endl;
            }
            return ss.str();
        }
    };
    
    typedef std::map<std::string, CFailureSummary > TFailureSummaries;
    TFailureSummaries aSuccessFailCounts;
    
    boost::mutex mxLogFailures;
    
    void log_int(const std::string & label, const eSuccessStatus result, const std::string & reason)
    {
        boost::mutex::scoped_lock lock(mxLogFailures);
        aSuccessFailCounts[label].log(result, reason);
    }
    
public:
    static void log(const std::string & label, const eSuccessStatus result, const std::string & reason)
    {
        s_evalFailures.log_int(label, result, reason);
    }
    
    ~CEvalFailure() 
    {
        for(const TFailureSummaries::value_type & val : aSuccessFailCounts)
        {
            cout << val.first << ":\n" << val.second.toString() << endl;
        }
    }
};

#endif // PR_H