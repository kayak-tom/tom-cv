/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */
#pragma once

#ifndef _SPEED_TEST
#define _SPEED_TEST


#ifndef __GNUC__
#include <windows.h>
#include <cstdlib>
#include <iostream>

typedef struct {
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} stopWatch;

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

unsigned long long rdtsc(void);

class CStopWatch {

private:
    int start, stop;
    long long startAsm,stopAsm;
    double LIToSecs( int & L) ;
public:
    CStopWatch() ;
    void startTimer( ) ;
    void stopTimer( ) ;
    double getElapsedTime() ;
};

#endif

//#define int_TIME2(label, line, ...) CStopWatch timer_ ## line; timer_ ## line.startTimer();\
//                                    const double dStartTimeSeconds_ ## line = UI().getElapsedTime();\
//                                    __VA_ARGS__; \
//                                    timer_ ## line.stopTimer(); \
//                                    cout << "TIME_FUNCTION: " << label << " took " << timer_ ## line.getElapsedTime() << " CPU seconds, " << (UI().getElapsedTime() - dStartTimeSeconds_ ## line) << " actual seconds" << endl;
//
//#define int_TIME1(label, line, ...) int_TIME2(label, line, __VA_ARGS__)
//
//#define TIME(label, ...) int_TIME1(label, __LINE__, __VA_ARGS__)

void testTiming();
int testFastNorms();

#endif //_SPEED_TEST