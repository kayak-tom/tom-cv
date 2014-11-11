/*
 * optimisation_attributes.h
 *
 *  Created on: 16 Apr 2010
 *      Author: tom
 */
#pragma once

#ifndef OPTIMISATION_ATTRIBUTES_H_
#define OPTIMISATION_ATTRIBUTES_H_

#ifdef __GNUC__
#  define HAVE_GCC_44 (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
#  define PURE_FN __attribute__((pure)) //Function has no consequences
#
#  ifdef __OPTIMIZE__
#    define HAVE_OPT(...) __VA_ARGS__
#  else
#    define HAVE_OPT(...)
#  endif
#
#  define HARD_INLINE HAVE_OPT(__attribute__((always_inline)))
#
#else
#  define HAVE_GCC_44 0
#  define HARD_INLINE
#  define PURE_FN
#endif

#if HAVE_GCC_44
#  define HOT __attribute__((hot)) HAVE_OPT(__attribute__((optimize("-O3"))) __attribute__((optimize("-ffast-math")))) //optimise for speed, even in debug
#  define COLD  __attribute__((cold)) //optimise for size
#else
#  define HOT
#  define COLD
#endif

//#define EXPECT(var, expectedVal) __builtin_expect(var, expectedVal)
#define EXPECT(var, expectedVal) var

#endif /* OPTIMISATION_ATTRIBUTES_H_ */
