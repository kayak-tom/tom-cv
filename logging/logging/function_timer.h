#ifndef FUNCTIONTIMER_H
#define FUNCTIONTIMER_H

#include <string>

class CEventTimer
{
    const std::string eventName;
    const double dStartTime;
public:
    CEventTimer(const std::string & eventName);
    ~CEventTimer();
    
    double getElapsedTime() const;

    static void startTimingEvents();
    static void pauseLogging();
    
    //Only true after autotests have finished
    static bool isTiming();
};


#endif // FUNCTIONTIMER_H
