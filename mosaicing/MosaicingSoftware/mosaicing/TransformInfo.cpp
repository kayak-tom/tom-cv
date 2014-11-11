//! 
#include "TransformInfo.h"
#include "GRCException.h"
#include "util/location.h"

namespace grc {

//! 
TransformInfo::TransformInfo(Transform * trans, const CorrespondenceSet * correspondencesIn, size_t id1, size_t id2) : trans_(trans), id1_(id1), id2_(id2), correspondences_(*correspondencesIn)
{
    // correspondences are copied
}

//! copy constructor
TransformInfo::TransformInfo(const TransformInfo & transInfo)
: correspondences_(transInfo.correspondences()),
    id1_(transInfo.id1()),
    id2_(transInfo.id2()),  
    trans_(0)
{
    const Transform * trans = transInfo.transform();
    
    if(!trans)
        throw new GRCException("TransformInfo::TransformInfo: Transform doesn't exist");
    
    if(const PerspectiveTransform * pt = dynamic_cast<const PerspectiveTransform *>(trans))
        trans_ = new PerspectiveTransform(*pt);
    else if(const AffineTransform * at = dynamic_cast<const AffineTransform *>(trans))
        trans_ = new AffineTransform(*at);
    else if(const SimilarityTransform * st = dynamic_cast<const SimilarityTransform *>(trans))
        trans_ = new SimilarityTransform(*st);
    else
        throw new GRCException("TransformInfo::TransformInfo: Transform copy failed");
}

TransformInfo::TransformInfo(Transform * trans, size_t id1, size_t id2) : trans_(trans), id1_(id1), id2_(id2)
{
    // no correspondences, no translation (this is either the identity taking last im to itself, or is only for rendering)
    /*if(id1 != id2)
        throw new GRCException( "TransformInfo::TransformInfo: No correspondences but different ids (can't refine transform)" );*/
}

TransformInfo * TransformInfo::accumulate(const TransformInfo * T2) const
{
    Transform * accumulatedTrans = trans_->accumulate(T2->transform());
    
    //We can only accumulate a->b * b->c = a->c
    if(id2() != T2->id1())
        throw new GRCException( "TransformInfo::accumulate: Id's do not match--invalid transformation accumulation" );

    return new TransformInfo(accumulatedTrans, id1(), T2->id2()); 
};

TransformInfo * TransformInfo::accumulateInverse(const TransformInfo * T2) const
{
    Transform * accumulatedTrans =trans_->accumulateInverse(T2->transform());
    
    //We can only accumulate (b->a)^-1 * b->c = a->c
    if(id1() != T2->id1())
        throw new GRCException( "TransformInfo::accumulateInverse: Id's do not match--invalid transformation accumulation" );

    return new TransformInfo(accumulatedTrans, T2->id2(), id2()); 
};
	

//!
TransformInfo2::TransformInfo2(Transform * trans, const CBoWCorrespondences * correspondencesIn, size_t id1, size_t id2) : trans_(trans), id1_(id1), id2_(id2)
{
	correspondencesIn->copyInto(correspondences_);
}

//! copy constructor
TransformInfo2::TransformInfo2(const TransformInfo2 & transInfo)
:   id1_(transInfo.id1()),
    id2_(transInfo.id2()),
    trans_(0)
{
	transInfo.correspondences().copyInto(correspondences_);

    const Transform * trans = transInfo.transform();

    if(!trans)
        throw new GRCException("TransformInfo2::TransformInfo2: Transform doesn't exist");

    if(const PerspectiveTransform * pt = dynamic_cast<const PerspectiveTransform *>(trans))
        trans_ = new PerspectiveTransform(*pt);
    else if(const AffineTransform * at = dynamic_cast<const AffineTransform *>(trans))
        trans_ = new AffineTransform(*at);
    else if(const SimilarityTransform * st = dynamic_cast<const SimilarityTransform *>(trans))
        trans_ = new SimilarityTransform(*st);
    else
        throw new GRCException("TransformInfo2::TransformInfo2: Transform copy failed");
}

TransformInfo2::TransformInfo2(Transform * trans, size_t id1, size_t id2) : trans_(trans), id1_(id1), id2_(id2)
{
    // no correspondences, no translation (this is either the identity taking last im to itself, or is only for rendering)
    /*if(id1 != id2)
        throw new GRCException( "TransformInfo2::TransformInfo2: No correspondences but different ids (can't refine transform)" );*/
}

TransformInfo2 * TransformInfo2::accumulate(const TransformInfo2 * T2) const
{
    Transform * accumulatedTrans = trans_->accumulate(T2->transform());

    //We can only accumulate a->b * b->c = a->c
    if(id2() != T2->id1())
        throw new GRCException( "TransformInfo2::accumulate: Id's do not match--invalid transformation accumulation" );

    return new TransformInfo2(accumulatedTrans, id1(), T2->id2());
};

TransformInfo2 * TransformInfo2::accumulateInverse(const TransformInfo2 * T2) const
{
    Transform * accumulatedTrans =trans_->accumulateInverse(T2->transform());

    //We can only accumulate (b->a)^-1 * b->c = a->c
    if(id1() != T2->id1())
        throw new GRCException( "TransformInfo2::accumulateInverse: Id's do not match--invalid transformation accumulation" );

    return new TransformInfo2(accumulatedTrans, T2->id2(), id2());
};
}
