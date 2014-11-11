/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "structSizeParams.h"
#include "image/convert_OpenCV.h"
#include <map>

class CSSObserver {
    static double s_dTypicalSS;
    static int s_nSSsObserved;
public:

    static void observeSS(double dSS) {
        if (dSS > 0) {
            std::cout << "Observed " << dSS << std::endl;
            s_dTypicalSS = (s_nSSsObserved * s_dTypicalSS + dSS) / (1 + s_nSSsObserved);
            s_nSSsObserved++;
        }
    }

    static double typicalSS() {
        return s_dTypicalSS;
    }
};

class C3dPointLocLoc {
public:
    C3dPoint point;
    CLocation loc1, loc2;

    C3dPointLocLoc(const C3dPoint & point, const CLocation loc1, const CLocation loc2) : point(point), loc1(loc1), loc2(loc2) {
        if(IS_DEBUG) CHECK(loc1.zero() || loc2.zero(), "C3dPointLocLoc: bad loc");
    }
};

typedef CDynArray<C3dPointLocLoc> T3dLocations;

class C3dPointCollection {
    typedef map2_NF<CLocation, const C3dPoint *> TPointMap;
    CDynArray<C3dPoint> aPoints; // points live here
    TPointMap locMap1, locMap2;
    double dStructureSize;

    int originalSize() const {
        return locMap1.size();
    }
public:
    //Must know 3d point count as use pointers into dynarray!! Todo fix

    C3dPointCollection(int n3dPoints) : dStructureSize(-1) {
        aPoints.reserve(n3dPoints);
    }

    double structureSize() const {
        return dStructureSize;
    }

    typedef TPointMap::const_iterator const_iterator;

    const_iterator begin1() const {
        return locMap1.begin();
    }

    const_iterator end1() const {
        return locMap1.end();
    }

    const_iterator begin2() const {
        return locMap2.begin();
    }

    const_iterator end2() const {
        return locMap2.end();
    }

    const C3dPoint * find1(CLocation loc) const {
        TPointMap::const_iterator pLoc = locMap1.find(loc);
        if (pLoc == locMap1.end())
            return 0;
        return pLoc->second;
    }

    const C3dPoint * findRandom(CLocation loc) const//This is to break SCORE by making a measurement from a differnt 3d point
    {
        TPointMap::const_iterator pLoc = locMap1.find(loc);
        if (pLoc == locMap1.end())
            return 0;
        return locMap1.begin()->second; //return WRONG point
    }

    const C3dPoint * find2(CLocation loc) const {
        //if (locMap2.count(loc)) return locMap2[loc];
        //return 0;
        TPointMap::const_iterator pLoc = locMap2.find(loc);
        if (pLoc == locMap2.end())
            return 0;
        return pLoc->second;
    }

    void insert(const C3dPoint & point, CLocation l1, CLocation l2) {
        if(IS_DEBUG) CHECK(l1.zero() || l2.zero(), "C3dPointCollection::insert: Bad params");
        if(IS_DEBUG) CHECK(/*!bSyntheticData &&*/(locMap1.find(l1) != locMap1.end() || locMap2.find(l2) != locMap2.end()), "C3dPointCollection::insert: Collision! This happens with S.D. as 2 points may be different as doubles but cast to the same place.");

        //aPoints MUST be pre-allocated
        aPoints.push_back(point);
        locMap1.init(l1, aPoints.end() - 1);
        locMap2.init(l2, aPoints.end() - 1);
    }

    void get3dLocations(const CDynArray<CLocation> & locIm1, T3dLocations & locationsOfFeature1) const {
        for (CDynArray<CLocation>::const_iterator pLoc = locIm1.begin(); pLoc != locIm1.end(); pLoc++) {
            const C3dPoint * pPoint = 0;
            if (true) {
                pPoint = find1(*pLoc);
            } else {
                REPEAT(10, cout << "WARNING: Breaking SCORE by returning random point rather than actual point\n");
                pPoint = findRandom(*pLoc);
            }

            if (pPoint) {
                //Find its locations in both images
                const_iterator pLocIm2 = begin2();
                for (; pLocIm2 != end2(); pLocIm2++)
                    if (pLocIm2->second == pPoint) break;
                if (pLocIm2 == end2()) {
                    std::cout << "ERROR: Missing point\n";
                    continue;
                }

                locationsOfFeature1.push_back(C3dPointLocLoc(*pPoint, *pLoc, pLocIm2->first));
                //cout << "Point found in im 1";
            }
            /*else if(pPoint2)
                    cout << "Point found in im 2";*/
        }
    }
};
