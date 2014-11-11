#pragma once

#include<string.h>
#include "GRCException.h"
#include<math.h>

//! For initialising integers.
#define MAX_INT 0x7FFFFFFF
#define MIN_INT (-MAX_INT) //0x80000000 doesn't work under linux

//! For initialising doubles.
#ifndef HUGE
#define HUGE 10e+50
#endif

// Conversion functions
//! Fast (not nan/inf etc. safe) typecast
inline int doubleToInt(double d)
{
#ifndef __GNUC__
  int i;

  _asm fld d
  _asm fistp i

  return i;
#else
  return (int)d;
#endif
}

inline int doubleToInt(float d)
{
#ifndef __GNUC__
  int i;

  _asm fld d
  _asm fistp i

  return i;
#else
  return (int)d;
#endif
}

template<class CHAR>
inline int euclidDist(const CHAR * ac1, const CHAR * ac2, int nLen)
{
    int d=0;
    for(int i=nLen; i>0; i--)
    {
        int temp = (int)(*ac1) - (int)(*ac2);
        d += temp*temp;
        ac1++;
        ac2++;
        //std::cout << temp << ' ';
    }

    return d;
}
