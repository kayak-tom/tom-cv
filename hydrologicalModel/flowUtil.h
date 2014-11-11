/* Drainage functions, etc.
 */
 
#ifndef _FLOWUTIL_H
#define _FLOWUTIL_H

inline double sqr(const double t) { return t*t; }

inline double pseudoHuber(const double x, const double t=0.01)
{
	const double dHuber = t * (sqrt( 1 + sqr(x/t) ) - 1);
	return dHuber;
}

#endif
