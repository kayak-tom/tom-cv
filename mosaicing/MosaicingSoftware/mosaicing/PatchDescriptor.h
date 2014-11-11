#pragma once
#ifndef _PATCHDESCRIPTOR_H
#define _PATCHDESCRIPTOR_H

#include "Descriptor.h"
#include "Patch.h"

namespace grc
{
#define cTPatchDescriptor cPatchDescriptor<PATCH_RADIUS, PATCH_COMP_METHOD, ImType, imDataType, PATCH_ORIENT, NORMALISE_PATCH, PATCH_SCALE>

//! Descriptor representing a patch of an image, and its location.
P_TEMPLATE
class cPatchDescriptor : private cTPatch, public cDescriptor
{
public:
    //! Distance from another descriptor.
    cDescriptor::TDist distance(const cDescriptor * pd) const
	{
        const cTPatch * pPatch = CAST<const cPatchDescriptor *>(pd);
        if(IS_DEBUG) CHECK(!pPatch, "cPatchDescriptor:distance: Descriptor not a cPatchDescriptor");
        cDescriptor::TDist dPatchDist = (cDescriptor::TDist)(cTPatch::distance(pPatch));
        return dPatchDist;
    };

    //! Fast approximation to distance from another descriptor.
    cDescriptor::TDist fastDistance(const cDescriptor * pd) const
	{
        const cTPatch * pPatch = CAST<const cPatchDescriptor *>(pd);
        if(IS_DEBUG) CHECK(!pPatch, "cPatchDescriptor:fastDistance: Descriptor not a cPatchDescriptor");
        cDescriptor::TDist dPatchDist = (cDescriptor::TDist)(cTPatch::fastDistance(pPatch));
        return dPatchDist;
    };

    //! Return how close we may be to the edge of the image.
    int margin() { return PATCH_SCALE*((PATCH_RADIUS*17/10)); }; //15/10 ~= root 2

    //! Compute the patch descriptor at this point
    cPatchDescriptor(const IplImage * pIm, CvPoint2D32f point) : cTPatch(pIm, cvPointFrom32f(point)), cDescriptor(point) {};
};

}

#endif //_PATCHDESCRIPTOR_H
