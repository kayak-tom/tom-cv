/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "SpeedTest.h"
#include <time.h>
#include "util/exception.h"
//#include "util/random.h"
#include "util/fastnorms.h"
#include "util/convert.h"
#include <iostream>

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

#include <stdlib.h> //causes errors in some Windows setups.

/*unsigned long long rdtsc(void)
{
    unsigned result[2];
    asm("rdtsc" : "=a" (result[0]), "=d" (result[1]));
    unsigned long long res = ((unsigned long long)result[1]) << 32;
    res += (unsigned long long)result[0];
    return res;
}*/

static inline unsigned long long read_rdtsc(void)
{
unsigned long long d;
__asm__ __volatile__ ("rdtsc" : "=A" (d) );
return d;
}

double CStopWatch::LIToSecs( int & L) {
    return ((double)L /(double)CLOCKS_PER_SEC) ;
}

CStopWatch::CStopWatch(){
    start=0;
    stop=0;
    startAsm = 0;
    stopAsm = 0;
    //QueryPerformanceFrequency( &frequency ) ;
}

void CStopWatch::startTimer( ) {
    start = clock() ;
    startAsm = read_rdtsc();
}

void CStopWatch::stopTimer( ) {
    stop = clock() ;
    stopAsm = read_rdtsc();
}

double CStopWatch::getElapsedTime() {
    int nTimeClock;
    nTimeClock = stop - start;
    double dTimeAsm = ((0.5/0.7) * 5.0e-10) * (double)(stopAsm - startAsm); // DEBUGONLY( * 0.5);
    if(nTimeClock<0) return 0;
    if(nTimeClock == 0)
    {
        if(dTimeAsm<0) return 0;
        return dTimeAsm;
    }
    double dTimeClock = LIToSecs( nTimeClock);

    //if(IS_DEBUG) CHECK((dtime>0.3) && (dTime2/dtime > 1.2 || dTime2 / dtime < 0.8), "Invald cpu clock time conversion");

//#ifndef _DEBUG //Fails in debug mode
    if((dTimeAsm>0.3) && (dTimeClock/dTimeAsm > 1.2 || dTimeClock / dTimeAsm < 0.8)) 
    {
        REPEAT(1, std::cout << "Invald cpu clock time conversion time=" << dTimeClock << "dTimeAsm=" << dTimeAsm <<  std::endl;);
    }
    //cout << dTimeAsm << " " << dTimeClock << endl;

    return dTimeClock < 3.0/CLOCKS_PER_SEC ? dTimeAsm : dTimeClock;
//#else
//    return dTimeClock;
//#endif
}
#endif

void testTiming()
{
    CStopWatch s, t;
    cout << " Testing timing :)\n";

    for(int nSpeed = 10; nSpeed < 20000; nSpeed += 10000)
    {
        t.startTimer();
        int nTotalTime = 0;
        double dTotalTime = 0;
        do
        {
            s.startTimer();
            for(int i=0; i<nSpeed; i++)
            {
                double dRandom = 1;  // CRandom::Uniform(0.02);
                while(dRandom > 0)
                    dRandom += rand() / (100.0*MAX_INT);
            }
            s.stopTimer();
            dTotalTime += s.getElapsedTime();
            if(dTotalTime > nTotalTime)
            {
                cout << nTotalTime << " seconds" << endl;
                nTotalTime = (int)dTotalTime + 1;
            }
        } while(nTotalTime < 20);

        t.stopTimer();
        cout << t.getElapsedTime() << " seconds by longer timer" << endl;
    }
}

template<typename CHAR, bool bNegVals>
int testFastNorms_int() {
    cout << "Start" << endl;

    const int LENGTH=128;

    CHAR ac1[LENGTH], ac2[LENGTH];
    int an1[LENGTH], an2[LENGTH];

    unsigned char * pc1 = reinterpret_cast<unsigned char *>(ac1);
    unsigned char * pc2 = reinterpret_cast<unsigned char *>(ac2);

    for(int i=0;i<LENGTH;i++)
    {
        int val1 = ((i+250) % 128);
        int val2 = ((i * 17) % 127);
        if(bNegVals) 
        {
            val1 -= 64;
            val2 -= 64;
        }
        ac1[i] = val1;
        ac2[i] = val2;
        an1[i] = val1;
        an2[i] = val2;
    }

    CStopWatch s; int nLength=0;

    s.startTimer();
    nLength = L1distNew<unsigned int, eL1>(pc1,pc2,LENGTH);
    s.stopTimer();
    cout << "Time New=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1dist<CHAR, eL1>(ac1,ac2,LENGTH);
    s.stopTimer();
    cout << "Time=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1dist<int, eL1>(an1,an2,LENGTH);
    s.stopTimer();
    cout << "Time unsigned int=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1dist<int, eL1>((int*)(void*)an1, (int*)(void*)an2, LENGTH);
    s.stopTimer();
    cout << "Time int=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1distNew<unsigned long long, eL1>(pc1, pc2,LENGTH);
    s.stopTimer();
    cout << "Time New 64=" << s.getElapsedTime() << " Length=" << nLength << endl;

    /*s.startTimer();
    nLength = L1distNewApprox<unsigned long long>(ac1,ac2,LENGTH);
    s.stopTimer();
    cout << "Time New Approx 64=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1distNewApprox<int>(ac1,ac2,LENGTH);
    s.stopTimer();
    cout << "Time New Approx=" << s.getElapsedTime() << " Length=" << nLength << endl;*/

    s.startTimer();
    nLength = L1distNew<unsigned int, eSSD>(pc1, pc2,LENGTH);
    s.stopTimer();
    cout << "SSD Time New=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1distNew<unsigned long long, eSSD>(pc1, pc2,LENGTH);
    s.stopTimer();
    cout << "SSD Time New 64=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1dist<CHAR, eSSD>(ac1,ac2,LENGTH);
    s.stopTimer();
    cout << "SSD Time=" << s.getElapsedTime() << " Length=" << nLength << endl;

    s.startTimer();
    nLength = L1dist<int, eSSD>(an1,an2,LENGTH);
    s.stopTimer();
    cout << "SSD Time unsigned int=" << s.getElapsedTime() << " Length=" << nLength << endl;

    int ss1=sumSquares(ac1, LENGTH);
    int ss2=sumSquares(ac2, LENGTH);
    //L2distCos returns a.^2 + 2a.b + b.^2
    //We want a.^2 - 2a.b + b.^2 = 2a.^2 + 2b.^2 - L2distCos

    s.startTimer();
    nLength = L2distParallel<unsigned int>(pc1,pc2,LENGTH);
    s.stopTimer();
    cout << "SSD Time L2Cos=" << s.getElapsedTime() << " Length=" <<  2*ss1 + 2*ss2 - nLength << endl;

    s.startTimer();
    nLength = L2distParallel<unsigned long long>(pc1,pc2,LENGTH);
    s.stopTimer();
    cout << "SSD Time L2Cos 64=" << s.getElapsedTime() << " Length=" << 2*ss1 + 2*ss2 - nLength  << endl;

    return 0;
}
int testFastNorms() {
    testFastNorms_int<unsigned char, false>();
    testFastNorms_int<char, false>();
    testFastNorms_int<unsigned char, true>();
    return testFastNorms_int<char, true>();
}