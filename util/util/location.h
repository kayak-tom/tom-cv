/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once
#ifndef LOCATION_H
#define LOCATION_H

#include "util/convert.h"
#include "util/dynArray.h"
#include "util/Simple2dPoint.h"

class CLocation //a descriptor location in an image. Equiv to 1 integer pre-point, but has subpixel precision (if SUBPIX_RES>1).
{
    static inline int PX_SCALE_INT()
    {
#ifndef SUBPIX_RES
#define SUBPIX_RES 30
#endif
        return SUBPIX_RES;
    } //Will work with images sized up to 1000x1000
    static inline double PX_SCALE_DB()
    {
        return (double)PX_SCALE_INT();
    } //Will work with images sized up to 1000x1000
    static inline double PX_SCALE_DB_INV()
    {
        return 1.0 / PX_SCALE_DB();
    } //Will work with images sized up to 1000x1000

    union {
        unsigned short asx[sizeof(int)/sizeof(unsigned short)];
        int nId;
    } data;
    int & thisAsInt() {
        return data.nId;
        //return *reinterpret_cast<int*>(reinterpret_cast<void*>(data.asx));
        /*union { int i; unsigned short aus[sizeof(int)/sizeof(unsigned short)]; } asInt;
        asInt.i = 0;
        asInt.aus[0] = data.asx[0];
        asInt.aus[1] = data.asx[1];
        return asInt.i;
        union { int * pi; unsigned short * as; } asInt;
        asInt.pi = 0;
        asInt.as = (unsigned short*)data.asx;
        return *(asInt.pi);*/
    }
    const int & thisAsInt() const {
        return data.nId;
        //return *reinterpret_cast<const int*>(data.asx);
    }

    //CHECK_SIZES_EQUAL(unsigned short[2], CLocation*, CLocation_1);
    //CHECK_SIZES_EQUAL(int, CLocation*, CLocation_2);
public:
    CHECK_SIZES_EQUAL(int, unsigned short[2], CLocation_2);

    int id() const { return thisAsInt(); }

    CLocation(double x, double y)
    {
        data.nId = 0;
        data.asx[0] = (unsigned short) (doubleToInt(PX_SCALE_DB() * x));
        data.asx[1] = (unsigned short) (doubleToInt(PX_SCALE_DB() * y));
    }
    CLocation(int x, int y)
    {
        data.nId = 0;
        data.asx[0] = (unsigned short) (PX_SCALE_INT() * x);
        data.asx[1] = (unsigned short) (PX_SCALE_INT() * y);
    }

    CLocation()
    {
        thisAsInt() = 0;
    }

    CLocation(int nId)
    {
        thisAsInt() = nId;
    }

    void operator=(const CLocation & loc)
    {
        thisAsInt() = loc.id();
    }

    bool operator==(const CLocation & loc) const
    {
        return loc.id() == id();
    }

    int x() const
    {
        return ((int) data.asx[0]) / PX_SCALE_INT();
    }
    int y() const
    {
        return ((int) data.asx[1]) / PX_SCALE_INT();
    }

    int xAsIntFast() const
    {
        return ((int) data.asx[0]);
    }
    int yAsIntFast() const
    {
        return ((int) data.asx[1]);
    }

    double dx() const
    {
        return ((double) data.asx[0]) * PX_SCALE_DB_INV();
    }
    double dy() const
    {
        return ((double) data.asx[1]) * PX_SCALE_DB_INV();
    }

    bool zero() const
    {
        return id() == 0;
    }

    //So we can use set to find duplicates etc.
    bool operator<(CLocation l) const
    {
        return id() < l.id();
    }

    static bool SUBPIX_SUPPORT()
    {
        return PX_SCALE_INT() > 1;
    }
};

template<int BINS>
class locationHash
{
public:
    int operator()(CLocation loc) const
    {
        unsigned int locAsInt = (unsigned int)loc.id();
        return locAsInt % BINS;
    }
};

//One feature correspondence
class CCorrespondence
{
    CLocation Loc1, Loc2;
    double dPP;
public:
    CCorrespondence(const CLocation &L1, const CLocation &L2, double dPP) : Loc1(L1), Loc2(L2), dPP(dPP) {}
    CCorrespondence()
    {
        //Used in mosaic work where correspondences are copied around THROW( "CCorrespondence::CCorrespondence Needed but should never actually be called");
        DEBUGONLY(REPEAT(1, std::cout << "Error--CCorrespondence::CCorrespondence needed but should never actually be called (except when mosaicing)\n"));
    }

    void zero() { dPP = 0; } //Flag to delete later
    double priorProb() const { return dPP; }

    const CLocation & Location1() const { return Loc1; }
    const CLocation & Location2() const { return Loc2; }
};

class CCamCalibMatrix;
class CPointIdentifiers;

class CBoWCorrespondences : public CDynArray<CCorrespondence>
{
public:
    CBoWCorrespondences() { CDynArray<CCorrespondence>::reserve(600); }

    //For debugging:
    bool contains(const CCorrespondence * pC) const
    {
        for(const_iterator pCorr = begin(); pCorr != end(); pCorr++)
            if(pCorr->Location1() == pC->Location1() && pCorr->Location2() == pC->Location2())
                return true;
        return false;
    }

    //Converts points to calibrated points, and also sets up pointIds (which tells RANSAC which correspondences are incompatible with each other)
    void calibrate(const CCamCalibMatrix & K, T2dPoints & aCalibratedPoints1, T2dPoints & aCalibratedPoints2, CInlierProbs & adArrLikelihood, CPointIdentifiers & pointIds, const bool CORRECT_RD) const;
};
#endif