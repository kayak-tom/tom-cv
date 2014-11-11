#include <util/exception.h>
#include <util/vectorUtil.h>
#include <algorithm>
#include <util/convert.h>
#include <boost/foreach.hpp>

void printVector(std::vector<double> & adNumbers)
{
    cout << adNumbers.size() << " elements: ";
    BOOST_FOREACH(const double d, adNumbers)
    {
        cout << d << ' ';
    }
    cout << endl;
}

double median(std::vector<double> & adNumbers)
{
    const bool bVerbose = false;
    
    if(bVerbose)
        printVector(adNumbers);
    
    if(IS_DEBUG) CHECK(adNumbers.size()==0, "Median not defined for 0 element vec");
    const int nMax = (int)adNumbers.size() / 2;
    std::vector<double>::iterator pMid = adNumbers.begin()+nMax;
    //std::nth_element(adNumbers.begin(), pMid, adNumbers.end());
    std::partial_sort(adNumbers.begin(), pMid, adNumbers.end());
    double dMedian = *(adNumbers.begin()+nMax);
    if(bVerbose) cout << "Element " << nMax << " of " << adNumbers.size() << " is " << dMedian << endl;
    if(adNumbers.size() % 2 == 0)
    {
        //std::vector<double>::iterator pMid2 = adNumbers.begin()+(nMax-1);
        //std::nth_element(adNumbers.begin(), pMid2, adNumbers.end());
        const double dMedian2 = *(adNumbers.begin()+(nMax-1));

        if(IS_DEBUG) CHECK(dMedian2 > dMedian, "Sort has failed somewhere");
    
        dMedian = 0.5 * (dMedian + dMedian2);
        if(bVerbose) cout << "Element " << (nMax-1) << " of " << adNumbers.size() << " is " << dMedian2 << endl;
    }
    //CHECKBADNUM(dMedian); // TODO: restore me
    if(bVerbose) cout << "Median " << dMedian << endl;

    if(bVerbose)
        printVector(adNumbers);

    return dMedian;
}

double percentile(const int nPercentile, std::vector<double> & aNumbers)
{
    const int nSize = (int)aNumbers.size();

    if(IS_DEBUG) CHECK(nSize == 0, "percentile: 0 elements supplied");

    const int medianIdx = (nSize*nPercentile)/100;

    std::vector<double>::iterator pPercentile = aNumbers.begin() + medianIdx;
    std::nth_element(aNumbers.begin(), pPercentile, aNumbers.end());
    double c = *pPercentile;

    return c;
}

void meanSD(std::vector<double> & adNumbers, double & dMean, double & dSD)
{
    dMean=dSD=0;
    if(adNumbers.size()==0)
        return;
    if(adNumbers.size()==1)
    {
        dMean=adNumbers[0];
        return;
    }
    double dSum=0,dSumSq=0;
    BOOST_FOREACH(const double d, adNumbers)
    {
        dSum += d;
        dSumSq += sqr(d);
    }
    
    dMean = dSum/adNumbers.size();
    dSD=sqrt( dSumSq/adNumbers.size() - sqr(dMean) );
}
