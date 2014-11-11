/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#ifndef GEOG_EIGEN_H
#define GEOG_EIGEN_H

#include "geom.h"
#include <Eigen/Core>

/* Decompose a homography matrix into a rotation, translation, and plane normal, by Levenberg least-squares optimisation.

   H = (R - t * n^t)s

   |n|==1

   Note SIGN AMBIGUITY--same solution if you flip signs of n and t.

   Note DEGENERATE if t==0. I think the rotation is still useful.

   If d is distance to plane, t*d is camera translation vector.

   Returns true on success.

   Uses SVD to find s (2nd s.v.). n and R^-1 * t are perpendicular to 2nd col of V, this is *not currently used*

   Closed form solution--TODO!

   See similar LM refinement in CEssentialMatGoldStandard in libransac
*/
bool decomposeHomography(const Eigen::Matrix3d & H, C3dRotation & rotation, C3dPoint & planeNormal, C3dPoint & camMotion) HOT;

//Create a homography from rotation, translation, normal.
void makeH(const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, Eigen::Matrix3d & HAtParams);

void makeE(const C3dRotation & rotation, const C3dPoint & cameraDir, Eigen::Matrix3d & EAtParams);
void makeE(const C3dRotation & rotation, const Eigen::Vector3d & cameraDir, Eigen::Matrix3d & EAtParams);
void makeE(const Eigen::Matrix3d & R, const Eigen::Vector3d & cameraDir, Eigen::Matrix3d & EAtParams);

//Refine E by minimising Sampson error (dist to corresponding epiline) by Levenberg-Marquardt refinement of R and T. Returns condition number (1.0/min diagonal of JTJ inverse)
double refineRT_LM_SampsonError(const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, const bool bAllResiduals, bool bVerbose) HOT;

//As above with robust cost function
class CRobustLMforEParams
{
public:
    const bool bRobustErr;
    const double dRobustifyThresh;
    const double dRobustifyScale;
    const double dHighResidConditionScale;
    const double dHighResidCutoff;

    CRobustLMforEParams() : bRobustErr(false), dRobustifyThresh(-1), dRobustifyScale(-1), dHighResidConditionScale(1), dHighResidCutoff(HUGE)
    {

    }

    CRobustLMforEParams(const bool bRobustErr, const double dRobustifyThresh, const double dRobustifyScale, const double dHighResidConditionScale, const double dHighResidCutoff)
     : bRobustErr(bRobustErr), dRobustifyThresh(dRobustifyThresh), dRobustifyScale(dRobustifyScale), dHighResidConditionScale(dHighResidConditionScale), dHighResidCutoff(dHighResidCutoff)
    {
        if(bRobustErr)
        {
            if(IS_DEBUG) CHECK(dRobustifyThresh <= 0 || dRobustifyScale < 0, "CRobustLMforEParams: Bad values");
        }

        REPEAT(1,
        if(dRobustifyThresh > 1)
            std::cout << "Warning: Unrealistically high val for dRobustifyThresh\n";
        if(dRobustifyScale > 1)
            std::cout << "Warning: dRobustifyScale should normally be less than 1, unless we are increasing penalty for outliers\n";

        if(IS_DEBUG) CHECK(dHighResidConditionScale < 1, "Not increasing penalty for bad scales");
        if(dHighResidCutoff > 0.1)
            std::cout << "Warning: Possible unrealistically high val for dHighResidCutoff\n");
    }
};

double refineRT_LM_SampsonError_robust(const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, const CRobustLMforEParams & robustRefinementParams, bool bVerbose) HOT;

//Normalises matrices first
double EssentialMat_SSD(const Eigen::Matrix3d & E1, const Eigen::Matrix3d & E2, bool bVerbose = true);

//Refining R and T to minimising reprojection err. *without* adjusting 3d points (reconstructed each time with a new cam instead)
void refineRT_ReprojErr_NotPoints(const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, bool bVerbose);

//full bundle adjust to refine points and camera
typedef CDynArray<C3dPoint> T3dPoints;
void refine3d(const CPointVec2d & p1, const CPointVec2d & p2, T3dPoints & points3d, C3dRotation & rotation, C3dPoint & cameraDir, bool bVerbose);
#endif //GEOG_EIGEN_H
