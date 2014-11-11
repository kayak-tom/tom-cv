/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * Simple2dPoint.h -- defines simple types often used in RANSAC
 *
 * Location in an image (pair of doubles). For E estimation use a calibrated point!
 *
 *
 *  Created on: 6/10/2009
 *      Author: tom
 */

#ifndef SIMPLE2DPOINT_H_
#define SIMPLE2DPOINT_H_

#include "dynArray.h"
#include "exception.h"
#include <iostream>

#define NN_MAX 12 //Global upper-limit on N and M in N-M correspondences (so we can statically-allocate memory)

class CCamCalibMatrix;
class CLocation;

class CSimple2dPoint
{
    friend void testK(const CCamCalibMatrix&, CLocation, bool);
protected:
    double x, y;
    void correctRD(const CCamCalibMatrix & K);//pixel coords->calibrated coords
    void uncorrectRD(const CCamCalibMatrix & K); //calibrated coords->pixel coords
public:
    CSimple2dPoint(const CLocation & l);

    inline CSimple2dPoint(double x, double y) : x(x), y(y) {};
    inline CSimple2dPoint() : x(0), y(0) {};
    inline double getX() const { return x; };
    inline double getY() const { return y; };

    inline bool operator==(const CSimple2dPoint & p) const { return x==p.x && y==p.y; };

    //Predicate for data structures--used for duplicate detection
    inline bool operator<(const CSimple2dPoint & p) const { return x<p.x || (x==p.x && y<p.y); };

    void calibrate(const CCamCalibMatrix & K, const bool CORRECT_RD);//pixel coords->calibrated coords
    void uncalibrate(const CCamCalibMatrix & K);//calibrated coords->pixel coords
};

typedef CDynArray<int> TSubSet;
typedef CDynArray<CSimple2dPoint> T2dPoints;

class CPointIds
{
    int nId1, nId2;
public:
    CPointIds(int n1, int n2) : nId1(n1), nId2(n2) {}
    CPointIds() : nId1(0), nId2(0) {}
    int id1() const { return nId1; }
    int id2() const { return nId2; }
};

typedef CDynArray<CPointIds> TPointIdentifiers; //uniquely identify each image point as an integer, used to handle N-N matches
class CPointIdentifiers : public TPointIdentifiers
{
public:
    CPointIdentifiers(int n) : TPointIdentifiers(n) {}
    CPointIdentifiers() {}

    inline bool incompatible(int i, int j) const
    {
        if(IS_DEBUG) CHECK(i==j, "Asking if a point is incompat with itself")
        if(IS_DEBUG) CHECK(i<0 || j<0, "Asking if a point is incompat with itself")

        const CPointIds & pi = get(i);
        const CPointIds & pj = get(j);
        return (pi.id1() == pj.id1()) || (pi.id2() == pj.id2());
    }
};

class CMask : public CDynArray<unsigned char> //Todo: what's faster??
{
public:
    CMask(int count) : CDynArray<unsigned char>(count, 0) {}
    CMask(const CMask & copyFrom) : CDynArray<unsigned char>(copyFrom.size(), 0)
    {
         copyFrom.copyInto(*this);
    }

    int countInliers() const
    {
        int nInlierCount = 0;

        for (int i = 0; i < (int)size(); i++)
            if ((*this)[i]) nInlierCount++;

        return nInlierCount;
    }
    void pp() const
    {
        std::cout << countInliers() << "/" << size()<< ": ";
        for (int i = 0; i < (int)size(); i++)
            std::cout << (((*this)[i]) ? '1' : '0');
        std::cout << std::endl;
    }
    int checksum() const
    {
        int nChecksum = 0;
        for (int i = 0; i < (int)size(); i++)
            nChecksum ^= (((*this)[i]) ? i+1 : 0);
        return nChecksum;
    }
    void setZero()
    {
        for (iterator p=begin(); p != end(); p++)
            *p = 0;
    }
};

typedef CDynArray<double> TInlierProbs;
class CInlierProbs : public TInlierProbs
{
public:
    CInlierProbs(int n) : TInlierProbs(n) {}
    CInlierProbs() {}
};


#endif /* SIMPLE2DPOINT_H_ */
