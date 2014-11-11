/* 
 * File:   houghForE.h
 * Author: tom
 *
 * Created on 13 September 2011, 3:28 PM
 */

#ifndef HOUGHFORE_H
#define	HOUGHFORE_H

#include "util/Simple2dPoint.h"
#include <Eigen/Core>

class C3x3MatModel;

class CHoughForE
{
public:
    virtual int findRTHough(const T2dPoints & points0, const T2dPoints & points1, const Eigen::Matrix3d & E_exact, const CMask & mask_exact, C3x3MatModel & model, CMask & mask) = 0;
    virtual ~CHoughForE() {};
    
    static CHoughForE * makeHough();
};

#endif	/* HOUGHFORE_H */

