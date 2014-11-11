#pragma once
#ifndef _DESCRIPTOR
#define _DESCRIPTOR

#include "util/opencv.h"
#include <vector>
#include "GRCException.h"
#include "Enums.h"

namespace grc {

#ifdef _DEBUG
#define CAST dynamic_cast //!<safe polymorphic typecast
#else
#define CAST static_cast //!<fast polymorphic typecast
#endif

//! Represents a feature in an image and a descriptor of the patch around it
class cDescriptor
{
    CvPoint2D32f loc_;//!< Location of feature

    //! Will be called once only to initialise margin
    static int getMargin_int();
public:
    typedef int TDist;

    virtual TDist distance(const cDescriptor  * pd) const = 0;
    virtual TDist fastDistance(const cDescriptor  * pd) const { return distance(pd); };
    
    cDescriptor(CvPoint2D32f loc) : loc_(loc) {};

    //! Point location
    CvPoint2D32f location() const { return loc_; };

    virtual ~cDescriptor() {};

    //! This is descriptor-implementation specific and is needed by the corner detector
    static int getMargin() { static int margin_s = getMargin_int(); return margin_s; }; 

    //! Return minimum distance to image edge for this descriptor (needed for corner detection).
    virtual int margin() = 0;

    static Enums::eDescriptorInvarianceType eDescriptorType_s; //!< Descriptor type
    static Enums::eDescriptorSize eDescriptorSize_s; //!< Descriptor patch size

    //! Descriptor factory method -- makes descriptors of type eDescriptorType_s and size eDescriptorSize_s.
    static cDescriptor * newDescriptor(const IplImage * pRGBImage, const CvPoint2D32f &Corner);
};

typedef std::vector<cDescriptor *> TDescriptorVector;

//! Set of descriptors (e.g. in an image), owns memory.
class DescriptorSet : public TDescriptorVector
{
public:

    ~DescriptorSet()
    {
        for(iterator ppDesc = begin(); ppDesc != end(); ppDesc++)
            delete *ppDesc;
    };

};

} //namespace grc

#include "PatchDescriptor.h"

#endif // _DESCRIPTOR
