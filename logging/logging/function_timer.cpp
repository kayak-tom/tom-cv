#include <util/exception.h>

#include "function_timer.h"
#include "imageViewUI.h"
#include <boost/scoped_ptr.hpp>
#include <fstream>
#include <vector>
#include <map>
#include <boost/thread/mutex.hpp>
#include <util/convert.h>

/*
 Measures clock time (seconds) taken by various slow functions, for testing and for checking for time variability (e.g. cost increasing as the number of frames grows)
* Logs to its own logfile/spreadsheet
* Doesn't save times from autotests, or if waiting in waitKey
* 
*/

class CTimedEvent
{
    double dStartTime, dEndTime;
public:
    CTimedEvent() : dStartTime(-1), dEndTime(-1) {}
    CTimedEvent(const double dStartTime, const double dEndTime) : dStartTime(dStartTime), dEndTime(dEndTime)
    {
        
    }
    
    double getElapsedTime() const { return dEndTime-dStartTime; }
    
    double getStartTime() const { return dStartTime; }
    double getEndTime() const { return dEndTime; }
    //int getTime() const { return timeStep; }
};

/*class CTimerLog_base
{
public:
    virtual void startTimingEvents(const std::string & pathToOutputDir) = 0;

    static CTimerLog_base * getTimerLog();
};*/

class CTimerLog //: public CTimerLog_base
{
    friend class CEventTimer;
    static boost::scoped_ptr<CTimerLog> pTimerLog;
    const double dTimeOrigin;
    double dWaitkeyTime; //don't log events if a waitkey falls in them
    const std::string pathToOutputDir;
    std::ofstream timingEventLog;
    typedef std::map<std::string, std::vector<CTimedEvent> > TTimeEvents;
    TTimeEvents aTimes;    
    boost::mutex mxTimeFunctions;
    
    void addDebugFlag(std::ostream & os) const 
    {
        if(IS_DEBUG)
            os << "DEBUG" << endl;
        else
            os << "RELEASE" << endl;
    }
    
    
public:
    CTimerLog(const std::string & pathToOutputDir) : dTimeOrigin(UI().getCurrentTimeSecs()), dWaitkeyTime(-1), pathToOutputDir(pathToOutputDir), timingEventLog((pathToOutputDir + "/timingEventLog.tsv").c_str())
    {
        addDebugFlag(timingEventLog);
        timingEventLog << "Timestep\tFnName\tDuration\tStartTime\tEndTime" << endl;
    }
    
    ~CTimerLog() 
    {
        if(aTimes.size()==0)
            return;
            
        // Make timing summary
        std::string strTimingSummary = pathToOutputDir + "/timing.tsv";
        std::ofstream timingSummary(strTimingSummary.c_str());
        addDebugFlag(timingSummary);
        
        timingSummary << "Event\tTotal\t#Calls\tMean\tRelSD\tMin\tMax\t" << endl;
        
        typedef std::map<double, std::string, std::greater<double> > TSortedTotals;
        TSortedTotals aTotals; //will wipe duplicates but doesn't really matter
        
        //BOOST_FOREACH(const TTimeEvents::value_type & timedEvent, aTimes)
        for(auto timedEvent : aTimes)
        {
            double dTotal=0, dSS=0, dMin=HUGE, dMax=0;
            for(const auto & t: timedEvent.second)
            {
                dTotal += t.getElapsedTime();
                dSS += sqr(t.getElapsedTime());
                if(t.getElapsedTime() < dMin)
                    dMin = t.getElapsedTime();
                if(t.getElapsedTime() > dMax)
                    dMax = t.getElapsedTime();
            }
            
            const double dN = (double)timedEvent.second.size();
            const double dMean = dTotal/dN;
            const double dSD = sqrt(dSS/dN - sqr(dMean));
            
            timingSummary << timedEvent.first << '\t' << dTotal << '\t' << dN << '\t' << dMean << '\t' << (dSD/dMean) << '\t' << dMin << '\t' << dMax << '\t' << endl;
            
            aTotals[dTotal] = timedEvent.first;
        }
        
        // Make list of total times
        std::string strTotals = pathToOutputDir + "/totalTimes.tsv";
        std::ofstream totals(strTotals.c_str());
        addDebugFlag(totals);
        totals << "Event\tTotalTime\t" << endl;
        for(const TSortedTotals::value_type & timedEvent: aTotals)
        {
            totals << timedEvent.second << '\t' << timedEvent.first << endl;
        }
    }
    
    static void log(const std::string & eventName, const double dStartTime, const double dEndTime)
    {
        if(!pTimerLog)
            return;
        
        pTimerLog->logEvent(eventName, dStartTime, dEndTime);
    }
    
    void pauseLogging()
    {
        dWaitkeyTime = UI().getCurrentTimeSecs();
    }
private:
    void logEvent(const std::string & eventName, const double dStartTime, const double dEndTime)
    {
        if(dStartTime <= dWaitkeyTime && dEndTime >= dWaitkeyTime) 
            return;
			
        boost::mutex::scoped_lock timerLock(mxTimeFunctions);
        
        CTimedEvent timedEvent(dStartTime-dTimeOrigin, dEndTime-dTimeOrigin);
        aTimes[eventName].push_back(timedEvent);
        
        timingEventLog << eventName << '\t' << timedEvent.getElapsedTime() << '\t' << timedEvent.getStartTime() << '\t' << timedEvent.getEndTime() << endl;
    }
    
};


boost::scoped_ptr<CTimerLog> CTimerLog::pTimerLog;


bool CEventTimer::isTiming()
{
    return (bool)CTimerLog::pTimerLog;
}

void CEventTimer::startTimingEvents()
{
    CHECK(CTimerLog::pTimerLog, "CTimerLog::pTimerLog already init");
    CTimerLog::pTimerLog.reset(new CTimerLog(UI().getLogPath()));
}

void CEventTimer::pauseLogging()
{
    if(CTimerLog::pTimerLog)
        CTimerLog::pTimerLog->pauseLogging();
}


CEventTimer::CEventTimer(const std::string & eventName) : eventName(eventName), dStartTime(UI().getCurrentTimeSecs()) {
    //s.startTimer();
}

double CEventTimer::getElapsedTime() const
{
    const double dEndTime = UI().getCurrentTimeSecs();
    if(dEndTime > dStartTime)
        return dEndTime - dStartTime;
    else
        return 0;
}

CEventTimer::~CEventTimer() {
    //s.stopTimer();
    const double dEndTime = UI().getCurrentTimeSecs();
    
    CTimerLog::log(eventName, dStartTime, dEndTime);
    /*if(adTimes.size()>0) {
        double dMin=HUGE, dMax = 0, dMean = 0;
        BOOST_FOREACH(const double dTime, adTimes) {
            if(dTime < dMin)
                dMin = dTime;
            if(dTime > dMax)
                dMax = dTime;

            dMean += dTime;
        }
        dMean /= adTimes.size();

        ALWAYS_VERBOSE;
        COUT2("Timing results: ", eventName);
        COUT(adTimes.size());
        COUT(dMean);
        COUT(dMin);
        COUT(dMax);
    }*/    
}
