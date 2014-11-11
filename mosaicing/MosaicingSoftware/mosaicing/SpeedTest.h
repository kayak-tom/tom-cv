/* \file SpeedTest.h \brief May be used for timing sections of code. */
#pragma once

#ifndef __GNUC__
#include <windows.h>
#include <cstdlib>
#include <iostream>

using namespace std;

typedef struct {
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} stopWatch;

//! Unused, but may be used for timing sections of code.
class CStopWatch {

private:
    stopWatch timer;
    LARGE_INTEGER frequency;
    double LIToSecs( LARGE_INTEGER & L) ;
public:
    CStopWatch() ;
    void startTimer( ) ;
    void stopTimer( ) ;
    double getElapsedTime() ;
};
#else
class CStopWatch {

private:
    int start, stop;
    double LIToSecs( int & L) ;
public:
    CStopWatch() ;
    void startTimer( ) ;
    void stopTimer( ) ;
    double getElapsedTime() ;
};

#endif
