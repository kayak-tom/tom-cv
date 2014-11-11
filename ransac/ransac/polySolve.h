/* 
 * File:   polySolve.h
 * Author: tom
 *
 * Created on 8 April 2011, 4:12 PM
 */

#ifndef POLYSOLVE_H
#define	POLYSOLVE_H

#define MAXDEGREE	10
#define MDP1	 MAXDEGREE+1
void rpoly_ak1(const double op[MDP1], int* Degree, double zeror[MAXDEGREE], double zeroi[MAXDEGREE]);
  
#endif	/* POLYSOLVE_H */

