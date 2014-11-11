/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#ifndef PATCH_DESCRIPTOR_H
#define PATCH_DESCRIPTOR_H

#include "descriptor.h"
#include "patch.h"

#ifndef PATCHFACTORY_H_
#pragma message("Don't #include this, its slow. You probably want featureExtract/patchFeatureExtractor.h instead")
#endif

#define CTPatchDescriptor CPatchDescriptor<dataType, PATCH_RADIUS, PATCH_COMP_METHOD, ImType, imDataType>

PN_TEMPLATE
class CPatchDescriptor : public CTPatchWithNorm, public CDescriptor
{
public:
	//typedef CTPatch::elType elType;

	CPatchDescriptor(double dScale, const IplImage * pIm, CLocation point, const CPatchParams & PATCH_PARAMS) : CTPatchWithNorm(pIm, CLocation(dScale*point.dx(), dScale*point.dy()), PATCH_PARAMS) { /*dScale; aDescriptor;*/ };

    CDescriptor::TDist distance(const CDescriptor * pd) const
	{
        const CTPatchWithNorm * pPatch = CAST<const CPatchDescriptor *>(pd);
        if(IS_DEBUG) CHECK(!pPatch, "CPatchDescriptor:distance: Descriptor not a CPatchDescriptor");
        CDescriptor::TDist dPatchDist = (CDescriptor::TDist)(CTPatchWithNorm::distance(pPatch));
        return dPatchDist;
    };

    static inline int DescriptorLength() { return CTPatch::SIZE; };
    static inline int SURFDescriptorLength() { return 0; };
    inline char * DescriptorVector() const { return (char *)(void *)this; };

	virtual bool drawable() const { return true; }
    virtual uchar val(int x, int y, int nChannel) const { return CTPatchWithNorm::val(x, y, nChannel); };
	virtual int diameter() const { return CTPatchWithNorm::DIAMETER; };

	virtual int size() const { return sizeof(this); }
};

PN_TEMPLATE
class CPatchAndLocationDescriptor : public CTPatchDescriptor
{
    CLocation Location;
public:
    CPatchAndLocationDescriptor(double dScale, const IplImage * pIm, CLocation point, const CPatchParams & PATCH_PARAMS) : CTPatchDescriptor(dScale, pIm, point, PATCH_PARAMS), Location(point) {};
    virtual CLocation location() const { return Location; };

	virtual int size() const { return sizeof(this); }
	virtual int length() const { return CTPatchDescriptor::SIZE; }
};
#endif
