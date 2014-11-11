/* 
 * File:   guidedFeatureMatching.h
 * Author: tom
 *
 * Created on 9 March 2012, 10:24 AM
 */

#ifndef GUIDEDFEATUREMATCHING_H
#define	GUIDEDFEATUREMATCHING_H

#include "description/descriptor.h"
#include <Eigen/Core>
#include "geom/geom.h"

class CLikelihoodFn;

class CGuidedFeatureMatching {
public:
    typedef Eigen::Matrix<double,6,6> TCovMat66;
    typedef Eigen::Matrix<double,9,9> TCovMat99;

    /*
     * Match features, with a prior from the IMU
     * 
     * pDS1, pDS2: descriptor sets
     * MS : matching settings (currently not used)
     * t : 3d translation vector
     * thetaPhiPsi : 3d vector with Euler angles in radians (TODO: confirm which order these are applied)
     * fullCovariance : 6x6 covariance matrix for ( t | fullCovariance )
     */
    static const CBoWCorrespondences * guidedFeatureMatch(const CMatchableDescriptors * pDS1, const CMatchableDescriptors * pDS2, const CCamCalibMatrix & K, const CMatchableDescriptors::CMatchSettings & MS, const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const TCovMat66 & fullCovariance);

    //As for guidedFeatureMatch, but points are contrained to lie on plane with normal n
    static const CBoWCorrespondences * guidedFeatureMatch_Homography(const CMatchableDescriptors * pDS1, const CMatchableDescriptors * pDS2, const CCamCalibMatrix & K, const CMatchableDescriptors::CMatchSettings & MS, const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const Eigen::Vector3d & n, const TCovMat99 & fullCovariance);
    
private:
    static const CBoWCorrespondences * guidedFeatureMatch_int(CLikelihoodFn * pFn, const CMatchableDescriptors * pDS1, const CMatchableDescriptors * pDS2, const CCamCalibMatrix & K, const CMatchableDescriptors::CMatchSettings & MS);

};

C3dRotation phiThetaPsiToQuat(const Eigen::Vector3d & phiThetaPsi);
Eigen::Matrix3d makeH(const Eigen::Vector3d & translationSample, const Eigen::Vector3d & phiThetaPsiSample, const Eigen::Vector3d & normalSample);
inline Eigen::Vector2d homog(const Eigen::Vector3d & v3) { return v3.head<2>() / v3(2); }

#endif	/* GUIDEDFEATUREMATCHING_H */

