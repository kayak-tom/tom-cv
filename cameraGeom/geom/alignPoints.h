/*
 * alignPoints.h
 *
 *  Created on: 10/12/2009
 *      Author: tom
 */

#ifndef ALIGNPOINTS_H_
#define ALIGNPOINTS_H_

#include "util/dynArray.h"
#include "geom.h"

class C3dPointMatch
{
    C3dPoint point1, point2;
public:
    C3dPointMatch(const C3dPoint & point1, const C3dPoint & point2) :
        point1(point1), point2(point2)    {};
    C3dPointMatch() {};

    const C3dPoint & p1() const { return point1; };
    const C3dPoint & p2() const { return point2; };
    ~C3dPointMatch() {  };

    bool zero() const { return point1.getX() + point1.getY() + point2.getX() + point2.getY()==0; }
    void setZero() { point1=C3dPoint(); point2=C3dPoint(); }
    bool operator==(const C3dPointMatch & m) const { return point1 == m.p1() && point2 == m.p2(); }
};

typedef CDynArray<C3dPointMatch> T3dPointMatchVector;

class CAlignPoints
{
public:
    CAlignPoints();
    virtual ~CAlignPoints();

    virtual bool alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s) = 0;
};

class CAlignPointsUmeyama : public CAlignPoints
{
public:
    CAlignPointsUmeyama(){};
    virtual ~CAlignPointsUmeyama(){};

    virtual bool alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s);
};

class CAlignPointsNew : public CAlignPoints
{
public:
    CAlignPointsNew() {};
    virtual ~CAlignPointsNew() {};

    virtual bool alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s);
};

class CAlignPointsUmeyama_Weighted : public CAlignPoints
{
public:
    CAlignPointsUmeyama_Weighted() {};
    virtual ~CAlignPointsUmeyama_Weighted() {};

    virtual bool alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s);
};

double median(std::vector<double> & adNumbers);

#endif /* ALIGNPOINTS_H_ */
