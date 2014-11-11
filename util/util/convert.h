/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/* maths / utility functions */

#pragma once

#define MAX_INT 0x7FFFFFFF
#define MIN_INT (-MAX_INT)
//0x80000000

#include<string.h>
#include "exception.h"
#include "optimisation_attributes.h"

#ifndef USE_MATH_DEFINES
#define USE_MATH_DEFINES
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include<math.h>
#include<cmath>

#ifndef HUGE
#define HUGE 10e+50
#endif

#ifndef __GNUC__
#include <float.h>
namespace std
{
    inline bool isnan(double d)
    {
        return _isnan(d) != 0;
    }
    inline bool isinf(double d)
    {
        return _finite(d) == 0;
    }
}
#endif

#include "util/floor.h"

//Array copying (memcpy) TODO: Speedups--pointers + wmemcpy

/*inline void arrayCopy(  char const * aSrc, double * aDst, unsigned int nCount )
{
    for(unsigned int i=nCount; i>0; i--)
    {
        *aDst = (double)(*aSrc);
        aDst++; aSrc++;
    }
}*/

template<typename DOUBLE, typename INT>
inline void arrayCopy( DOUBLE const * aSrc, INT * aDst, unsigned int nCount )
{
    for(unsigned int i=nCount; i>0; i--)
    {
        *aDst = (INT)doubleToInt((DOUBLE)*aSrc);
        aDst++; aSrc++;
    }
}

/*inline void arrayCopy(  float const * aSrc, int * aDst, unsigned int nCount )
{
    for(unsigned int i=nCount; i>0; i--)
    {
        *aDst = (int)doubleToInt(*aSrc);
        aDst++; aSrc++;
    }
}

inline void arrayCopy(  double const * aSrc, unsigned char * aDst, unsigned int nCount )
{
    for(unsigned int i=nCount; i>0; i--)
    {
        *aDst = (unsigned char)doubleToInt(*aSrc);
        aDst++; aSrc++;
    }
}*/

template<class DOUBLE, class CHAR>
inline void arrayCopyScale(  DOUBLE const * aSrc, CHAR * aDst, DOUBLE dScale, unsigned int nCount )
{
    if(dScale==1) return arrayCopy(  aSrc, aDst, nCount );

    for(unsigned int i=nCount; i>0; i--)
    {
        volatile int val = doubleToInt(dScale * *aSrc);
        if(IS_DEBUG) CHECK((int)((CHAR)val) != val, "Conversion to char overflowed");
        *aDst = (CHAR)val;
        aDst++; aSrc++;
    }
}

//template<class CHAR>
//inline void arrayCopyScale(  float const * aSrc, CHAR * aDst, float dScale, unsigned int nCount )
//{
//    if(dScale==1) return arrayCopy(  aSrc, aDst, nCount );
//
//    for(unsigned int i=nCount; i>0; i--)
//    {
//        *aDst = (CHAR)floatToInt(dScale * *aSrc);
//        aDst++; aSrc++;
//    }
//}

inline void arrayCopy(  double const * aSrc, double * aDst, unsigned int nCount)
{
    memcpy(aDst, aSrc, nCount*sizeof(double));
}

inline void arrayCopy(  char const * aSrc, char * aDst, unsigned int nCount)
{
    memcpy(aDst, aSrc, nCount*sizeof(char));
}

//Define dist functions for each type here:
inline int cosDist(const char * ac1, const char * ac2, int nLen)
{
    int d = 128*128;
    for(int i=0;i<nLen;i++)
        d -= (int)ac1[i] * (int)ac2[i]; //todo: pointer arithmatic

    return d > 0 ? d : 0;
}

inline double cosDist(const double * ad1, const double * ad2, int nLen)
{
    double d=1;
    for(int i=0;i<nLen;i++)
        d -= ad1[i] * ad2[i]; //todo: pointer arithmatic

    return d;

}

inline double euclidDist(const double * ad1, const double * ad2, int nLen)
{
    double d=0;
    for(int i=0;i<nLen;i++)
    {
        volatile double temp = ad1[i] - ad2[i];
        d += temp*temp;
    }

    return d;
}

template<class CHAR>
inline int euclidDist(const CHAR * ac1, const CHAR * ac2, int nLen)
{
    int d=0;
    for(int i=nLen; EXPECT(i>0, 1); i--)
    {
        volatile int temp = (int)(*ac1) - (int)(*ac2);
        d += temp*temp;
        ac1++;
        ac2++;
    }

    return d;
}
template<class CHAR, int INNERLOOP, int LENGTH>
inline int euclidDistUnroll(const CHAR * ac1, const CHAR * ac2, int);

//Go fast by unrolling inner loop. Manually unrolling also doesn't help
template<class CHAR, int INNERLOOP, int LENGTH>
inline int euclidDistUnroll(const CHAR * ac1, const CHAR * ac2, int)
{
    int d=0;
    int i=LENGTH/INNERLOOP;
    int numLeft = LENGTH - i*INNERLOOP;
    for(; i>0; i--) //__builtin_expect doesn't help
    {
        for(int j=INNERLOOP; j>0; j--)
        {
            volatile int temp = (int)(*ac1) - (int)(*ac2); //Indexing into arrays is fractionally worse
            d += temp*temp;
            ac1++;
            ac2++;
        }
    }
    //There's a few left over...
    for(; numLeft>0; numLeft--)
    {
        volatile int temp = (int)(*ac1) - (int)(*ac2);
        d += temp*temp;
        ac1++;
        ac2++;
    }

    return d;
}

template<class CHAR>
inline int L1DistSlower(const CHAR * ac1, const CHAR * ac2, int nLen)
{
    int d = 0;
    for(int i=nLen;i>0;i--)
    {
        volatile int temp = (int)*ac1 - (int)*ac2;
        d += (temp>0) ? temp : -temp;
        ac1++;
        ac2++;
    }

    return d;
}

template<class CHAR>
inline int L1Dist(const CHAR * ac1, const CHAR * ac2, int nLen)
{
    int d = 0;
    for(int i=nLen;i>0;i--)
    {
        volatile int temp = (int)*ac1 - (int)*ac2;
        if(temp>0)
            d += temp;
        else
            d -= temp;

        ac1++;
        ac2++;
    }

    return d;
}

template<class CHAR>
inline int MaxDistSlower(const CHAR * ac1, const CHAR * ac2, int nLen)
{
    int d = 0;
    for(int i=nLen; i>0; i--)
    {
        volatile int dist = (int)*ac1 - (int)*ac2;
        d = dist > 0 ? (dist>d ? dist : d) : (-dist>d ? -dist : d);
        ac1++;
        ac2++;
    }

    return d;
}

template<class CHAR>
inline int MaxDist(const CHAR * ac1, const CHAR * ac2, int nLen)
{
    int d = 0;
    for(int i=nLen; i>0; i--)
    {
        volatile int dist = (int)*ac1 - (int)*ac2;
        if(dist<0) dist = -dist;
        if(d<dist) d=dist;
        ac1++;
        ac2++;
    }

    return d;
}

#define COUNTHITS(x)

class intLookup
{
    static const int nLookupSize = 256*256*8; //*128 is highest possible
    static int anSQRT[nLookupSize];
    static double adUcharToDouble[256];
    static float afUcharToFloat[256];

    static bool bSetup;
public:
    COUNTHITS(static int nHit; static int nMiss);
    inline static int Sqrt(int n)
    {
        if(IS_DEBUG) CHECK((n<0 || !bSetup), "intLookup::Sqrt OOB or not setup");
        if(n >= nLookupSize) //never actually happens
        {
            COUNTHITS(nMiss++);
            return doubleToInt(sqrt((double)n));
        }
        COUNTHITS(nHit++);
        return anSQRT[n];
    };

    static void Setup()
    {
        for(int i=0; i<nLookupSize; i++)
            anSQRT[i] = doubleToInt(sqrt((double)i));

        for(int i=0;i<256;i++)
        {
            adUcharToDouble[i] = (double)i;
            afUcharToFloat[i] = (float)i;
        }

        bSetup = true;
    };

    static inline float uchar2float(unsigned char c) { return afUcharToFloat[c]; }
    static inline double uchar2double(unsigned char c) { return adUcharToDouble[c]; }
};

inline double sqr(const double x) { return x*x; }
inline double cube(const double x) { return x*x*x; }
inline double fourthpow(const double x) { return sqr(sqr(x)); }

inline float sqr(const float x) { return x*x; }
inline float cube(const float x) { return x*x*x; }
inline float fourthpow(const float x) { return sqr(sqr(x)); }

inline int sqr(const int x) { return x*x; }

double nCr(double n,double r);

inline bool zero(double d)
{
    //Decide if determinants, etc. are too close to 0 to bother with
    const double EPSILON = 1e-3;
    return (d<EPSILON) && (d>-EPSILON);
}

inline bool is64bit() HARD_INLINE;

inline bool is64bit()
{
    return (sizeof(void*) == 8);
}

inline int idFromPtr(const void * p)
{

    union { const void * p; int an[sizeof(void*)/sizeof(int)]; } uConvert;
    uConvert.p = p;

/*#ifndef __GNUC__
#  pragma warning(push)
#  pragma warning(disable:4127) //Conditional expression is constant
#endif*/

    if(is64bit())
        return uConvert.an[0] ^ uConvert.an[1];
    else
        return uConvert.an[0];

/*#ifndef __GNUC__
#  pragma warning(pop)
#endif*/
}

template<typename T> inline void setConstant(T * array, T val, int nCount)
{
    for(; nCount>0; nCount--)
    {
        *array = val;
        array++;
    }
}

template<typename T> inline void setZero(T * array, int nCount)
{
    setConstant<T>(array, 0, nCount);
}

/**
 * @brief val2 is within a fraction thresh of val1
 * @param val1
 * @param val2
 * @param thresh
 * @return 
 */
inline bool within(double val1, double val2, double thresh)
{
	if(val1<0)
	    return val1*(1+thresh)<=val2 && val1*(1-thresh)>=val2;
	else
		return val1*(1+thresh)>=val2 && val1*(1-thresh)<=val2;
}

inline void checkProb(double p)
{
    if(IS_DEBUG) CHECK(p<0 || p>1, "Probability outside 0...1");
}

inline void checkProb(float p)
{
    if(IS_DEBUG) CHECK(p<0 || p>1, "Probability outside 0...1");
}

inline void checkTotalProb(double dP1, double dP2=0, double dP3=0, double dP4=0, double dP5=0)
{
    DEBUGONLY(
    checkProb(dP1);
    checkProb(dP2);
    checkProb(dP3);
    checkProb(dP4);
    checkProb(dP5);

    double sum = dP1 + dP2 + dP3 + dP4 + dP5;
    if(!zero(1-sum))
    {
        cout << sum << ": ";
        THROW("Total probabilities don't sum to 1");
    });
}

inline double diffMod(const double d1, const double d2, const double mod/* = 2*M_PI*/)
{
    const double diff = fabs(d1-d2);
    return (diff < mod*0.5) ? diff : (mod-diff);
}

template<typename T> inline T max3(const T a, const T b, const T c)
{
    return std::max<T>(a, std::max<T>(b, c));
}
template<typename T> inline T max4(const T a, const T b, const T c, const T d)
{
    return std::max<T>(std::max<T>(a, b), std::max<T>(c, d));
}
template<typename T> inline T min3(const T a, const T b, const T c)
{
    return std::min<T>(a, std::min<T>(b, c));
}
template<typename T> inline T min4(const T a, const T b, const T c, const T d)
{
    return std::min<T>(std::min<T>(a, b), std::min<T>(c, d));
}
template<typename T> inline T median3(const T a, const T b, const T c)
{
    const T nMax=max3(a,b,c);
    if(a==nMax)
        return std::max<T>(b,c);
    else if (b==nMax)
        return std::max<T>(a,c);
    //else c is max
    return std::max<T>(a,b);
}

//Clips to a range so that x->closest in [a,b]
template<typename T> inline T clip(const T x, const T a, const T b)
{
    if(IS_DEBUG) 
    {
        CHECK(a > b, "Clip range end error");
        CHECKBADNUM(x);
    }
    if(x<a)
        return a;
    if(x>b)
        return b;

    return x;
}

inline double log_b(const double d, const double b) { return log(d)/log(b); }
inline double log_2(const double d) { return log_b(d,2.0); }

inline double exp_b(const double d, const double b) { return pow(b, d); }

#define ARRAYZ(T, atName, length) ARRAY(T, atName, length); setZero(PTR(atName), length)
