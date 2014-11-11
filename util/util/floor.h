/* 
 * File:   floor.h
 * Author: tom
 *
 * Created on 27 January 2011, 16:40
 */

#pragma once
#ifndef FLOOR_H
#define	FLOOR_H

/* This one is fast but doesn't work

 typedef double ff_lreal;
typedef float  ff_real;
typedef unsigned long ff_uint32;
typedef long ff_int32;

const ff_lreal _double2fixmagic = 68719476736.0*1.5;     //2^36 * 1.5,  (52-_shiftamt=36) uses limited precisicion to floor
const ff_int32 _shiftamt        = 16;                    //16.16 fixed point representation,

#if BigEndian_
	#define iexp_				0
	#define iman_				1
#else
	#define iexp_				1
	#define iman_				0
#endif //BigEndian_

// ================================================================================================
// Real2Int
// ================================================================================================
inline ff_int32 Real2Int(ff_lreal val)
{
#if DEFAULT_CONVERSION
	return val;
#else
	val		= val + _double2fixmagic;
	return ((ff_int32*)&val)[iman_] >> _shiftamt;
#endif
}

// ================================================================================================
// Real2Int
// ================================================================================================
inline ff_int32 Real2Int(ff_real val)
{
#if DEFAULT_CONVERSION
	return val;
#else
	return Real2Int ((ff_lreal)val);
#endif
}*/

/*OLD conversions:*/

 // Conversion functions
#define FLOAT_TO_INT(in,out)  \
            __asm__ __volatile__ ("fistpl %0" : "=m" (out) : "t" (in) : "st") ;

inline int doubleToIntOld(double d) PURE_FN HARD_INLINE;
inline int doubleToIntOld(double d)
{
#ifndef __GNUC__ 
	#ifndef _WIN64
	  int i;

	  _asm fld d
	  _asm fistp i

	  return i;
	#else
	  return (int)d;
	#endif
#else

  volatile int n;
  FLOAT_TO_INT(d, n)
  return n;
#endif
}
inline int doubleToIntOld(float d) PURE_FN HARD_INLINE;
inline int doubleToIntOld(float d)
{
#ifndef __GNUC__ // Windows
    #ifndef _WIN64
		int i;

	  _asm fld d
	  _asm fistp i

	  return i;
	#else //Win64
	  return (int)d;
	#endif
#else
  volatile int n;
  FLOAT_TO_INT(d, n)
  return n;
#endif
}
/*inline int intFloor(double d) HARD_INLINE;
inline int intFloor(double d)
{
    int f = doubleToInt(d);

    if((double)f>d)
        f--;

    if(IS_DEBUG) CHECK(f != (int)floor(d), "intFloor failed");
    return f;
}*/

inline int doubleToInt(float d) PURE_FN HARD_INLINE;
inline int doubleToInt(double d) PURE_FN HARD_INLINE;
inline int intFloor(double val) PURE_FN HARD_INLINE;
inline int intFloor(float val) PURE_FN HARD_INLINE;

inline int doubleToInt(float d) { return doubleToIntOld(d); }
inline int doubleToInt(double d) { return doubleToIntOld(d); }

inline int intFloor(double val)
{
    int n = doubleToInt(val);
    return((double)n > val) ? (n-1) : n;
}

inline int intFloor(float val)
{
    int n = doubleToInt(val);
    return((double)n > val) ? (n-1) : n;
}


#endif	/* FLOOR_H */

