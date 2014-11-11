#include "SpeedTest.h"
#include <time.h>

#ifndef __GNUC__
using namespace std;

double CStopWatch::LIToSecs( LARGE_INTEGER & L) {
    return ((double)L.QuadPart /(double)frequency.QuadPart) ;
}

CStopWatch::CStopWatch(){
    timer.start.QuadPart=0;
    timer.stop.QuadPart=0;
    QueryPerformanceFrequency( &frequency ) ;
}

void CStopWatch::startTimer( ) {
    QueryPerformanceCounter(&timer.start) ;
}

void CStopWatch::stopTimer( ) {
    QueryPerformanceCounter(&timer.stop) ;
}

double CStopWatch::getElapsedTime() {
    LARGE_INTEGER time;
    time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
    return LIToSecs( time) ;
}
#else
double CStopWatch::LIToSecs( int & L) {
    return ((double)L /(double)CLOCKS_PER_SEC) ;
}

CStopWatch::CStopWatch(){
    start=0;
    stop=0;
    //QueryPerformanceFrequency( &frequency ) ;
}

void CStopWatch::startTimer( ) {
    start = clock() ;
}

void CStopWatch::stopTimer( ) {
    stop = clock() ;
}

double CStopWatch::getElapsedTime() {
    int time;
    time = stop - start;
    return LIToSecs( time) ;
}
#endif
