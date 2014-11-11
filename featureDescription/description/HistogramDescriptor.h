/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "Histogram.h"

extern int HIST_FACTOR;
extern int HIST_NORM_LB;

#define TEMPLATE_HV template<int HIST_BINS_H, int HIST_BINS_S, int HIST_BINS_I, int HIST_RADIUS, eHistCompMethod HIST_COMP_METHOD, bool HIST_CUMULATIVE,\
	     class descriptorType>

#define CTHistAndVectorDescriptor CHistAndVectorDescriptor<HIST_BINS_H, HIST_BINS_S, HIST_BINS_I, HIST_RADIUS, HIST_COMP_METHOD, HIST_CUMULATIVE, descriptorType>
#define CTHistDescriptor CHistDescriptor<HIST_BINS_H, HIST_BINS_S, HIST_BINS_I, HIST_RADIUS, HIST_COMP_METHOD, HIST_CUMULATIVE>

TEMPLATE_HV
class CHistAndVectorDescriptor : public descriptorType, public CTHist
{
public:
    CHistAndVectorDescriptor(const typename descriptorType::elType * aDescriptor, const ImageRGB * pIm, CvPoint point) : descriptorType(aDescriptor), CTHist(pIm, point) {};
	CHistAndVectorDescriptor(const double * aDescriptor, double dScale, const ImageRGB * pIm, CvPoint point) : descriptorType(aDescriptor, dScale), CTHist(pIm, point)
    {
        //double dist=abs(Compare(this));
        //double dist=abs(CTVectorSpaceDescriptor::distance(this));
        //DEBUGONLY(double dist=abs((double)distance(this)));
        //if(IS_DEBUG) CHECK(dist > 512, "CHistAndVectorDescriptor:CHistAndVectorDescriptor: distance non-zero");
    };
    CHistAndVectorDescriptor(const float * aDescriptor, float dScale, const ImageRGB * pIm, CvPoint point) : descriptorType(aDescriptor, dScale), CTHist(pIm, point) {};

    CDescriptor::TDist distance(const CDescriptor * pd) const
	{
        const CTHist * pHist = CAST<const CHistAndVectorDescriptor *>(pd);
        if(IS_DEBUG) CHECK(!pHist, "CHistAndVectorDescriptor:distance: Descriptor not a CHistAndVectorDescriptor");
        CDescriptor::TDist dSURFdist = descriptorType::distance(pd) ;
        CDescriptor::TDist dHDist = (CDescriptor::TDist)((HIST_FACTOR * CTHist::Compare(pHist))/10) - HIST_NORM_LB;
        if( dHDist < 0 ) dHDist = 0;


        CDescriptor::TDist dDist = dSURFdist + dHDist; //note 0*, need to turn down distinc.
        if(IS_DEBUG) CHECK(dDist < 0, "CHistAndVectorDescriptor:distance: distance negative");
        return dDist;
	};

#ifndef __GNUC__
    static int DescriptorLength() { return descriptorType::DescriptorLength() + HIST_BINS; };
#else
    static int DescriptorLength() { return descriptorType::DescriptorLength() + HIST_BINS_H+HIST_BINS_S+HIST_BINS_I; };
#endif
    static int margin() { return (HIST_RADIUS+1); };

    //Hacky k-means constructor that unconcatenates vector back into 2 char arrays for hist and SURF
    //Avoid using as we are converting signed/unsigned chars
    CHistAndVectorDescriptor(const double * aDescriptor) : descriptorType(aDescriptor, 1), CTHist(aDescriptor + descriptorType::DescriptorLength()) {};
    CHistAndVectorDescriptor(const char * aDescriptor) : descriptorType(aDescriptor), CTHist(aDescriptor + descriptorType::DescriptorLength()) {};

	virtual int size() const { return sizeof(this); }
};

H_TEMPLATE
class CHistDescriptor : private CTHist, public CDescriptor
{
public:
    //First 2 params dont do anything
    CHistDescriptor(const double * aDescriptor, double dScale, const ImageRGB * pIm, CvPoint point) : CTHist(pIm, point) { /*dScale; aDescriptor;*/ };

    CDescriptor::TDist distance(const CDescriptor * pd) const
	{
        const CTHist * pHist = CAST<const CHistDescriptor *>(pd);
        if(IS_DEBUG) CHECK(!pHist, "CHistDescriptor:distance: Descriptor not a CHistDescriptor");
        CDescriptor::TDist dHDist = (CDescriptor::TDist)((HIST_FACTOR * CTHist::Compare(pHist))/10);
        return dHDist;
    };

    static inline int DescriptorLength() { return HIST_BINS_H+HIST_BINS_S+HIST_BINS_I; };
    static inline int SURFDescriptorLength() { return 0; };
    typedef char elType;
    inline char * DescriptorVector() const { return (char *)(void *)this; };

	virtual int size() const { return sizeof(this); }
	virtual int length() const { return HIST_BINS; }
};

H_TEMPLATE
class CHistAndLocationDescriptor : public CTHistDescriptor
{
    CLocation Location;
public:
    CHistAndLocationDescriptor(const double * aDescriptor, double dScale, const ImageRGB * pIm, CvPoint point) : CTHistDescriptor(aDescriptor, dScale, pIm, point), Location(point.x, point.y) {};
    virtual CLocation location() const { return Location; };

	virtual int size() const { return sizeof(this); }
};

TEMPLATE_HV
class CHistVectorLocationDescriptor : public CTHistAndVectorDescriptor
{
    CLocation Location;
public:
    CHistVectorLocationDescriptor(const double * aDescriptor, double dScale, const ImageRGB * pIm, CvPoint point) : CTHistAndVectorDescriptor(aDescriptor, dScale, pIm, point), Location(point.x, point.y) {};

    CHistVectorLocationDescriptor(const double * aDescriptor, int x, int y) : CTHistAndVectorDescriptor(aDescriptor), Location(x, y)  {};
    CHistVectorLocationDescriptor(const char * aDescriptor, int x, int y) : CTHistAndVectorDescriptor(aDescriptor), Location(x, y)  {};
    CHistVectorLocationDescriptor(const char * aDescriptor) : CTHistAndVectorDescriptor(aDescriptor), Location(0, 0)  {};

    virtual CLocation location() const { return Location; };

	virtual int size() const { return sizeof(this); }
};

