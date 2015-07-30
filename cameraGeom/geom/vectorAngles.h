
#ifndef GEOMETRY_H
#define GEOMETRY_H
/*
 * Points, lines, planes, distances to/from etc.
 */

#include "newCamera.h"

class C3dBoundedLine;
class C2dBoundedLine;
class C2dLine;
class C3dLine;

/**
 * @brief Shorten the length of a residual vector by dRadius (at most). Used when a residual is computed to a 3D line, then truncated to the boundary of the cylinder.
 * @param resid
 * @param dRadius
 */
void shortenResidual(TEigen3dPoint & resid, const double dRadius);

//boost::optional<const C2dBoundedLine> projectBoundedLine(const CWorldCamera & P, const C3dBoundedLine & line);

/**
 * @brief Project ray from p out the the depth of the bounded line then find the closest point. 
 * @param p
 * @param P
 * @param line
 * @return Point ON line (hence on the 3D component) which is closest to the measured point's ray
 */
const C3dWorldPoint threeDPointNearLineAndMeasurement(const C2dImagePointPx & p, const CWorldCamera & P, const C3dBoundedLine & line);
const C3dWorldPoint threeDPointNearLineAndMeasurement(const C2dImagePointPx & p, const CWorldCamera & P, const C3dLine & line);

const C3dBoundedLine threeDLineBetweenLineAndMeasurement(const C2dImagePointPx & p, const CWorldCamera & P, const C3dBoundedLine & line);

/**
 * @brief Ensure that new wire is broken into chunks (about) dSegmentLength long
 * @param line
 * @return 
 */
void splitLine(const C3dBoundedLine & line, T3dPointVector & line_split, const double dSegmentLength);

//Return true if incompatible
bool incompatLines(const C2dBoundedLine & line1, const C2dBoundedLine & line2);

template<class TPoint>
double angleBetween3Points(const TPoint & p1, const TPoint & p2, const TPoint & p3);

template<class TPoint>
double angleBetweenVectors(const TPoint & direction1, const TPoint & direction2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/ );

template<class TPoint>
double angleBetweenUnitVectors(const TPoint & dir1, const TPoint & dir2, const bool bDirected);


/**
 * @brief Returns angle in range -Pi...Pi
 * @param direction1
 * @param direction1
 * @return 
 */
double signedAngleBetweenVectors(const TEigen2dPoint & direction1, const TEigen2dPoint & direction2);
double signedAngleBetweenUnitVectors(const TEigen2dPoint & direction1, const TEigen2dPoint & direction2);

double safe_acos(const double dCosAngle);

template<class TVecType>
double tanAngleBetweenVectors(const TVecType & seg1vec, const TVecType & seg2vec);

//Return a vector perpendicular to v, of the same magnitude
TEigen2dPoint perpendicular(const TEigen2dPoint & v);

optional<C3dWorldPoint> reconstructPointToLine(const CWorldCamera & P_point, const C2dImagePointPx & x, const CWorldCamera & P_line, const C2dLine & line);

#endif // GEOMETRY_H
