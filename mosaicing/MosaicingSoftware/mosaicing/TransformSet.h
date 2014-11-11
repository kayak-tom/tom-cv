#pragma once
#ifndef _TRANSFORM_SET
#define _TRANSFORM_SET

#include "util/set2.h"
#include "TransformInfo.h"
#include <iostream>

namespace grc {

    //! Pair of ids--identify transforms, etc. in STL containers.
    typedef std::pair<size_t, size_t> IdPair_int;
    
    //! Wrapper class around a pair of ids, used to identify a pairwise transformation.
    class IdPair : protected IdPair_int 
    {
        friend bool std::less<size_t>::operator ()(const size_t &, const size_t &) const;
        friend bool std::less<IdPair>::operator ()(const IdPair &, const IdPair &) const;
    public:
        IdPair(size_t im1, size_t im2) : IdPair_int(im1, im2) {};
        //! Nicer more meaningful names than .first and .second, especially when a member of a pair itself. 
        size_t im1Id() const { return first; };

        //! Nicer more meaningful names than .first and .second, especially when a member of a pair itself. 
        size_t im2Id() const { return second; };
    };

    //! Parent of TransformSet
    typedef map2_NF<IdPair, TransformInfo2 *> TTransformSet;

    //! Set of transforms, wrapped in TransformInfo wrappers. Owns all its memory (copies are made) to be dealt with by the renderer.
    class TransformSet : public TTransformSet
    {
        //Sometimes we only need to render the latest frame:
        bool onlyRenderLatestTransform_;
        IdPair latestTransform_;
        Transform * mosaicTransform_;
    public:
        //! Create empty set of transforms
        TransformSet();
        
        //! Destructor deletes all members
        ~TransformSet();

        //! Erase member (if exists), and delete
        void erase(const IdPair & id);
        
        bool exists(IdPair ids) const;

        //! Used to test whether to render incrementally.
        bool onlyRenderLatestTransform() const;

        //! Returns latest transform (the only one to be rendered when we render incrementally.)
        const TransformInfo2 * latestTransform() const;

        //! Returns transform to apply to the entire old mosaic (to avoid the new image sliding off the screen) when we render incrementally.
        const Transform * mosaicTransform() const;

        //! Setup incremental rendering. The engine calls this when it decides the frame can be incrementally rendered.
		 /*! \param mosaicTrans The transform to apply to the entire old mosaic
		  * \param latestTrans The index of the transform to apply to the most recent frame
          */
        void setFastRender(Transform * mosaicTrans, IdPair latestTrans);
        void setLatestTransform(IdPair latestTrans);
    };

}


#endif // _TRANSFORM_SET
