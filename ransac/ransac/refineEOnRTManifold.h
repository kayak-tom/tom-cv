/* 
 * File:   refineEOnRTManifold.h
 * Author: tom
 *
 * Created on 28 June 2011, 4:07 PM
 */

#ifndef REFINEEONRTMANIFOLD_H
#define	REFINEEONRTMANIFOLD_H

#include "geom/geom.h"
#include "refiner.h"

class CRefineEOnRTManifold// : public CImCorrModelRefiner {
{
public:
    //CRefineEOnRTManifold() : CImCorrModelRefiner(5, p1, p2) {}
    
    //virtual ~CRefineEOnRTManifold();
    
private:
    
    //virtual bool fitModel_int(CMask & mask, CModel & pModel, bool bVerbose);

public:

    static double refineRobustOnAll(const T2dPoints & p1, const T2dPoints & p2, C3dRotation & R, C3dPoint & t);
    static double refineRobustOnMask(const T2dPoints & p1, const T2dPoints & p2, C3dRotation & R, C3dPoint & t, CMask & mask);

    static double refineLSOnMask(const T2dPoints & p1, const T2dPoints & p2, C3dRotation & R, C3dPoint & t_in, CMask & mask);

};

double refineWeightedLinearRobustOnAll(const T2dPoints & p1, const T2dPoints & p2, Eigen::Matrix3d & E, const bool bRefineF);
double refineWeightedLinearRobustOnMask(const T2dPoints & p1, const T2dPoints & p2, Eigen::Matrix3d & E, CMask & mask, const bool bRefineF);
double refineWeightedLinearLSOnMask(const T2dPoints & p1, const T2dPoints & p2, Eigen::Matrix3d & E, CMask & mask, const bool bRefineF);

#endif	/* REFINEEONRTMANIFOLD_H */

