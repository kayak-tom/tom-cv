#pragma once
#ifndef _TRANSFORM_INFO
#define _TRANSFORM_INFO

#include "util/location.h"
#include "Transform.h"
#include "CorrespondenceSet.h"

//class CBoWCorrespondences;

namespace grc {

//! Contains transform and enough information to refine with Bundle Adjustment
/*! Contains transform and enough information (correspondences) to refine with Bundle Adjustment. This
 *  info is also used for internal consistency checks.
 */
class TransformInfo 
{
    Transform * trans_; //!<The transform we're wrapping
    size_t id1_, id2_; //!< Ids of the frames this transform is between.
    CorrespondenceSet correspondences_; //!< Holds a copy of the correspondences used to compute this transform that can be changed by BA
public:
    size_t id1() const { return id1_; };
    size_t id2() const { return id2_; };

    TransformInfo(Transform * trans, const CorrespondenceSet * correspondencesIn, size_t id1, size_t id2);
    TransformInfo(Transform * trans, size_t id1, size_t id2);

    //! Copy constructor--the engine and renderer work only on copies of the original.
    TransformInfo(const TransformInfo & transInfo);

    ~TransformInfo() { delete trans_; };

    const Transform * transform() const { return trans_; };
    Transform * transform() { return trans_; };

    const CorrespondenceSet & correspondences() const { return correspondences_; };

    //! Use to accumulate (multiply) transforms together
    /*! T:a->b * T:b->c = T:a->c. Internal consistency checks used to verify order is correct */
    TransformInfo * accumulate(const TransformInfo * T2) const;

    //! Use to accumulate (multiply) transforms together
    /*! T:a->b * (T:c->b)^-1 = T:a->c. Internal consistency checks used to verify order is correct */
    TransformInfo * accumulateInverse(const TransformInfo * T2) const;
};

class TransformInfo2
{
    Transform * trans_; //!<The transform we're wrapping
    size_t id1_, id2_; //!< Ids of the frames this transform is between.
    CBoWCorrespondences correspondences_; //!< Holds a copy of the correspondences used to compute this transform that can be changed by BA
public:
    size_t id1() const { return id1_; };
    size_t id2() const { return id2_; };

    TransformInfo2(Transform * trans, const CBoWCorrespondences * correspondencesIn, size_t id1, size_t id2);
    TransformInfo2(Transform * trans, size_t id1, size_t id2);

    //! Copy constructor--the engine and renderer work only on copies of the original.
    TransformInfo2(const TransformInfo2 & transInfo);
    TransformInfo2() : trans_(0),id1_(-1), id2_(-1) {}

    ~TransformInfo2() { delete trans_; };

    const Transform * transform() const { return trans_; };
    Transform * transform() { return trans_; };

    const CBoWCorrespondences & correspondences() const { return correspondences_; };

    //! Use to accumulate (multiply) transforms together
    /*! T:a->b * T:b->c = T:a->c. Internal consistency checks used to verify order is correct */
    TransformInfo2 * accumulate(const TransformInfo2 * T2) const;

    //! Use to accumulate (multiply) transforms together
    /*! T:a->b * (T:c->b)^-1 = T:a->c. Internal consistency checks used to verify order is correct */
    TransformInfo2 * accumulateInverse(const TransformInfo2 * T2) const;
};

}


#endif // _TRANSFORM_INFO
