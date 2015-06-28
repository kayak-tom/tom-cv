/* 
 * File:   closestPoint.h
 * Author: hilandtom
 *
 * Created on January 18, 2012, 7:02 AM
 */

#ifndef CLOSESTPOINT_H
#define    CLOSESTPOINT_H

#include <Eigen/Core>
#include <util/dynArray.h>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>

template<typename TFloat, int DIMS, int NCLOSEST>
class CExactClosestPoint
{
public:
    typedef Eigen::Matrix<TFloat, DIMS, 1> TVec;
protected:
    CDynArray<TVec> aAllPoints;
public:
    void addUniquePoint(const TVec & vec)
    {
        aAllPoints.push_back(vec);
    }
    void addPoint(const TVec & vec)
    {
        BOOST_FOREACH(const TVec & existingVec, aAllPoints)
        {
            if(existingVec == vec)
                return;                
        }
        aAllPoints.push_back(vec);
    }
    
    void findNClosest(const TVec & vec, TVec * aNClosest, TFloat * adDistancesSq ) const
    {
        if(IS_DEBUG) CHECK(aAllPoints.size() < NCLOSEST, "Too few points to find N closest");
        
        setConstant(adDistancesSq, HUGE, NCLOSEST);
        
        BOOST_FOREACH(const TVec & p, aAllPoints)
        {
            TFloat dDistSq = (p-vec).squaredNorm();
            
            if(dDistSq >= adDistancesSq[NCLOSEST-1])
                continue;
            
            int nInsertPoint = NCLOSEST-1;
            while(nInsertPoint>0 && dDistSq < adDistancesSq[nInsertPoint-1])
                nInsertPoint--;
                
            for(int i=NCLOSEST-1; i>nInsertPoint; i--)
            {
                adDistancesSq[i] = adDistancesSq[i-1];
                aNClosest[i] = aNClosest[i-1];
            }
            adDistancesSq[nInsertPoint] = dDistSq;
            aNClosest[nInsertPoint] = p;            
        }
    }   
    
    virtual ~CExactClosestPoint() {}
};

template<typename TFloat, int DIMS, int NCLOSEST>
class CExactClosestPointFast : public CExactClosestPoint<TFloat, DIMS, NCLOSEST>
{
    typedef CExactClosestPoint<TFloat, DIMS, NCLOSEST> TExactClosest;
    boost::multi_array<TVec, DIMS> aPointMap;    
    TVec origin, scale_inv;
    int anBins[DIMS], MAX_BINS;
    static const int MAX_DIMS = 4;
    boost::multi_array< TExactClosest, DIMS > aSubBins;
    
    const TExactClosest & closestBin(int * anBins) const
    {
        switch(DIMS)
        {
            case 1:
                return aSubBins[anBins[0]];
            case 2:
                return aSubBins[anBins[0]];
            case 3:
                return aSubBins[anBins[0]][anBins[1]][anBins[2]];
            case 4:
                return aSubBins[anBins[0]][anBins[1]][anBins[2]][anBins[3]];
        }        
        THROW("Bad num of dims");
    }
    
    TExactClosest & closestBin(int * anBins)
    {
        switch(DIMS)
        {
            case 1:
                return aSubBins[anBins[0]];
            case 2:
                return aSubBins[anBins[0]];
            case 3:
                return aSubBins[anBins[0]][anBins[1]][anBins[2]];
            case 4:
                return aSubBins[anBins[0]][anBins[1]][anBins[2]][anBins[3]];
        }        
        THROW("Bad num of dims");
    }
    
public:
    CExactClosestPointFast(int nMaxBins) : MAX_BINS(nMaxBins) {}
    
    void train()
    {
        //Find ranges
        TVec min_point, max_point;
        min_point.setConstant(HUGE);
        max_point.setConstant(-HUGE);
        
        BOOST_FOREACH(const TVec & p, aAllPoints)
        {
            min_point = min_point.min(p);
            max_point = max_point.max(p);
        }
        
        TVec range = max_point - min_point;
        max_point += range * 0.1;//Expand a bit
        min_point -= range * 0.1;
        range = max_point - min_point;//Reset expanded range
        
        TFloat volume = 1;
        for(int i=0;i<DIMS;i++)
            volume *= range(i);
        
        TFloat dVolPerBin = volume/MAX_BINS;
        TFloat dSideLength = pow(dVolPerBin, 1.0/DIMS);
        
        for(int i=0;i<DIMS;i++)
        {
            anBins[i] = ceil(range/dSideLength);
            scale_inv(i) = anBins[i]/range(i);
        }
        
        switch(DIMS)
        {
            case 1:
                aSubBins.resize(boost::extents[anBins[0]]);
                break;
            case 2:
                aSubBins.resize(boost::extents[anBins[0]][anBins[1]]);
                break;
            case 3:
                aSubBins.resize(boost::extents[anBins[0]][anBins[1]][anBins[2]]);
                break;
            case 4:
                aSubBins.resize(boost::extents[anBins[0]][anBins[1]][anBins[2]][anBins[3]]);
                break;
        }
        
        //For each vertex, add N closest
        for(int i0 = 0; i0 <= anBins[0]; i0++)
        for(int i1 = 0; i1 <= ((DIMS>=2) ? anBins[1] : 0); i1++)
        for(int i2 = 0; i2 <= ((DIMS>=3) ? anBins[2] : 0); i2++)
        for(int i3 = 0; i3 <= ((DIMS>=4) ? anBins[3] : 0); i3++)
        {
            int anBin[MAX_DIMS] = {i0,i1,i2,i3};
            TVec point = indexToPoint(anBin);
            TVec aNClosest[NCLOSEST];
            double adDistances[NCLOSEST];
            findNClosest(point, aNClosest, adDistances);
            //Add to all 2^DIMS surrounding bins
            //todo
        }
        THROW("Todo: not implemented");
        
        
        //For each point, add to corresponding bin [almost always unnecessary]
        //Would adding the point to its neighbours' bins help?
        int anBin[DIMS];

        BOOST_FOREACH(const TVec & p, aAllPoints)
        {
            if(toIndex(vec, anBin))
            {
                const TExactClosest & closestPoints = closestBin(anBin);
                closestPoints.findNClosest(vec, aNClosest, adDistancesSq); //will include vec

                //add all to this bin
            }
            else
            {
                THROW("Should have found bins encompassing all points")
            }
        
        }
    }
    
    //Return false if OOB
    bool toIndex(const TVec & vec, int * anBin) const
    {
        TVec indices = (vec - origin).array() * scale_inv.array();
        
        for(int i=0;i<DIMS;i++)
        {
            anBin[i] = intFloor(indices(i));
            if(anBin[i] < 0 || anBin[i] > anBins[i])
                return false; //OOB
        }
        return true;
    }
        //Return false if OOB
    TVec indexToPoint(int * anBin) const
    {
        TVec indices;
        for(int i=0;i<DIMS;i++)
            indices(i) = anBin[i];
        
        return (indices.array() / scale_inv.array()) + origin ;
    }
    
    void findNClosest_fast(const TVec & vec, TVec * aNClosest, TFloat * adDistancesSq ) const
    {
        int anBin[DIMS];
        if(toIndex(vec, anBin))
        {
            const TExactClosest & closestPoints = closestBin(anBin);
            closestPoints.findNClosest(vec, aNClosest, adDistancesSq);
        }
        else
        {
            findNClosest(vec, aNClosest, adDistancesSq);
        }
    }
    
};



#endif    /* CLOSESTPOINT_H */

