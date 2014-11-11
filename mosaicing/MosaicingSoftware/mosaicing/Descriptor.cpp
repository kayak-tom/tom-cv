#include "Descriptor.h"

//#define TOrientedNormalisedPatchDescriptor cPatchDescriptor<5, grc::ePatchEuclidSquared, grc::RGB, unsigned char, false /*Oriented*/, true /*Normalised*/, 2 /*subsample*/>

namespace grc
{

Enums::eDescriptorInvarianceType cDescriptor::eDescriptorType_s = (Enums::eDescriptorInvarianceType)-1;
Enums::eDescriptorSize cDescriptor::eDescriptorSize_s = (Enums::eDescriptorSize)-1;

#define RETURN_NEW_DESCRIPTOR2(size, oriented, normalised) \
return new cPatchDescriptor<size, grc::ePatchEuclidSquared, grc::RGB, unsigned char, oriented, normalised, 2 /*subsample*/>(pRGBImage, Corner);

#define RETURN_NEW_DESCRIPTOR(oriented, normalised) \
    switch(eDescriptorSize_s){\
    case Enums::e9x9:\
        RETURN_NEW_DESCRIPTOR2(4, oriented, normalised)\
    case Enums::e11x11:\
        RETURN_NEW_DESCRIPTOR2(5, oriented, normalised)\
    case Enums::e13x13:\
        RETURN_NEW_DESCRIPTOR2(6, oriented, normalised)\
    case Enums::e15x15:\
        RETURN_NEW_DESCRIPTOR2(7, oriented, normalised)}

//! Descriptor factory method
/*! Complex because instantiating one of many template classes--so need to have build them all */
cDescriptor * cDescriptor::newDescriptor(const IplImage * pRGBImage, const CvPoint2D32f &Corner)
{
    //const int aPatchRadii[4] = {4, 5, 6, 7};

    switch(eDescriptorType_s)
    {
    case Enums::eSimplePatch:
        RETURN_NEW_DESCRIPTOR(false, false);
    case Enums::eOriented:
        RETURN_NEW_DESCRIPTOR(true, false);
    case Enums::eNormalised:
        RETURN_NEW_DESCRIPTOR(false, true);
    case Enums::eOrientedNormalised:
        RETURN_NEW_DESCRIPTOR(true, true);
    };
    throw new GRCException("cDescriptor::newDescriptor: Unknown descriptor type");
}

//!Will be called once only to init margin
int cDescriptor::getMargin_int()
{
   static int callCount = 0;
   if(callCount>0)
        throw new GRCException("cDescriptor::getMargin_int: This should only ever be called once");

   cDescriptor * desc = newDescriptor(0, cvPoint2D32f(0,0));

   int margin = desc->margin();

   delete desc;

   callCount++;

   return margin;
}

}
