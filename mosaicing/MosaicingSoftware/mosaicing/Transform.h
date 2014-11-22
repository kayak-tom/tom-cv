#pragma once
#ifndef _TRANSFORM
#define _TRANSFORM

#include "util/opencv.h"
#include <vector>
#include "Enums.h"

namespace grc {

    //! Represents a transform
    /*! Derived classes implement all transform-specific code--that is: estimation from 
     *  correspondences, accumulation, application to points or for warping one image 
     *  into another. Must be compatible with LM parameter refinement.
     */
    class Transform : public CvMat
    {
    public:
        //! Safer cast from polymorphic type to non-polymorphic struct
        operator CvMat *() const { return (CvMat *)this; };

        //! Easy access to matrix elements.
        virtual double & operator()(int x, int y) = 0;
        double operator()(int x, int y) const;

        //! A parameter is a pointer to a number (double) in a transform.
        /*! Refinements are scaled depending on the number, to make the sizes of derivatives relatively similar (e.g. of the matrix elements representing translation vs. those representing skew). */
        class param
        {
            friend class TransformErrFunction;
            double *paramRef_, scale_, scale_inv_;
        public:
            param(double *paramRef, double scale) : paramRef_(paramRef), scale_(scale), scale_inv_(1. / scale) {};
        };
        typedef std::vector<param> TParamLocations;

        //! Add all the tunable parameters for this transform to the parameter vector
        virtual void addParams(TParamLocations * paramVector) = 0;

        static Enums::eTransformType eTransType_s; //!<Controls the type of transform used throughout the application.
        static int warpMethod_s;//!<Controls whether to use NN or Bilinear interpolation when warping images.
        
        //! Factory methods for new transforms
        static Transform * newTransform();
        static Transform * copyTransform(const Transform * T);
        
        //! Multiply together transforms. Allocs new memory
        virtual Transform * accumulate(const Transform * T2) const = 0;
        
        //! Multiply by an inverted transform. Allocs new memory
        virtual Transform * accumulateInverse(const Transform * T2) const = 0;

        //! Returns the component of the transform representing translation. Used for deciding how overlapping transforms are.
        virtual CvPoint2D64f translation() const = 0;
        
        //! Apply a shift (translate) this transform (for mosaic positioning)
        virtual void shift(double x, double y);

        //! Use this transform to warp images
        virtual void applyToImage(const IplImage * sourceIm, IplImage * destIm) const = 0;
        //! Use this transform to transform a point
        virtual CvPoint2D64f applyToPoint(CvPoint2D64f p) const = 0;

        //Estimation:
        
        //! Apply to a matrix of points
        virtual void applyToPoints(const CvMat * positions, CvMat * newPositions) const = 0;
        
        //! Return a suitable threshhold (in pixels) for inliers for this transform type.
        virtual double getRANSACThreshhold() const = 0;

        //! Return the minimum number of points required to estimateFromPoints.
        virtual size_t getRANSACSetSize() const { return 4; };
        
        //! Estimate transform from a set of points (use a least-squares estimate)
        virtual void estimateFromPoints(const CvMat * points1, const CvMat * points2) = 0;
        
        virtual ~Transform() {}
    };

    //!Implements a perspective transform
    /*! Implemented as a 3x3 homegeous matrix. */
    class PerspectiveTransform : public Transform
    {
        static const int ROWS=3;
        static const int COLS=3;
        double data_[ROWS*COLS]; //!< Transform matrix data
    public:
        PerspectiveTransform();
        PerspectiveTransform(const PerspectiveTransform & T);

        virtual double & operator()(int x, int y);
        virtual CvPoint2D64f applyToPoint(CvPoint2D64f p) const;
        virtual PerspectiveTransform * accumulate(const Transform * T2) const;
        virtual PerspectiveTransform * accumulateInverse(const Transform * T2) const;

        //! Add all the tunable parameters for this transform to the parameter vector
        virtual void addParams(TParamLocations * paramVector);

        virtual void applyToImage(const IplImage * sourceIm, IplImage * destIm) const;

        //!Apply to a matrix of points
        virtual void applyToPoints(const CvMat * positions, CvMat * newPositions) const;
        
        //!Estimate transform from a set of points
        virtual void estimateFromPoints(const CvMat * points1, const CvMat * points2);
        virtual double getRANSACThreshhold() const { return 3.5; };

        virtual CvPoint2D64f translation() const;
    };

    class SimilarityTransform;

    //! Implements Affine Transform
    /*! Affine transform: approximation to perspective transform with
     * [0 0 1] as bottom row of matrix (so no homegeous division is necessary). */
    class AffineTransform : public Transform
    {
        static const int ROWS=2;
        static const int COLS=3;
        double data_[ROWS*COLS];//!< Transform matrix data

        AffineTransform * accumulateInt(const Transform * T2, bool bInvert) const;
    public:
        AffineTransform();
        AffineTransform(const AffineTransform & T);

        virtual double & operator()(int x, int y);
        virtual CvPoint2D64f applyToPoint(CvPoint2D64f p) const;
        virtual AffineTransform * accumulate(const Transform * T2) const;
        virtual AffineTransform * accumulateInverse(const Transform * T2) const;

        //! Add all the tunable parameters for this transform to the parameter vector
        virtual void addParams(TParamLocations * paramVector);

        SimilarityTransform * getSimilarityTransform() const;

        virtual void applyToImage(const IplImage * sourceIm, IplImage * destIm) const;

        //!Apply to a matrix of points
        virtual void applyToPoints(const CvMat * positions, CvMat * newPositions) const;
        
        //!Estimate transform from a set of points
        virtual void estimateFromPoints(const CvMat * points1, const CvMat * points2);
        virtual double getRANSACThreshhold() const { return 5; };

        virtual CvPoint2D64f translation() const;
    };

    //! Similarity (Rotation/Scale/translation) transform
    /*!  Implemented using getSimilarityTransform() and getAffineTransform to convert too/from affine transform
     */
    class SimilarityTransform : public Transform
    {
        double theta_, //!< Similarity transform angle
            scale_;//!< Similarity transform scale factor
        CvPoint2D64f translation_;//!< Similarity transform translation

        //! Allow easy rendering using cvWarpAffine
        AffineTransform * getAffineTransform() const; 
    public:
        SimilarityTransform(double theta_, double scale, const CvPoint2D64f & translation);
        SimilarityTransform(const SimilarityTransform & T);
        SimilarityTransform();

        virtual double & operator()(int x, int y);
        virtual CvPoint2D64f applyToPoint(CvPoint2D64f p) const;
        virtual SimilarityTransform * accumulate(const Transform * T2) const;
        virtual SimilarityTransform * accumulateInverse(const Transform * T2) const;

        //! Add all the tunable parameters for this transform to the parameter vector
        virtual void addParams(TParamLocations * paramVector);

        virtual void applyToImage(const IplImage * sourceIm, IplImage * destIm) const;
        virtual void shift(double x, double y);

        //!Apply to a matrix of points
        virtual void applyToPoints(const CvMat * positions, CvMat * newPositions) const;
        
        //!Estimate transform from a set of points
        virtual void estimateFromPoints(const CvMat * points1, const CvMat * points2);
        virtual double getRANSACThreshhold() const { return 3; };
        virtual size_t getRANSACSetSize() const { return 3; };

        virtual CvPoint2D64f translation() const;
    };

    //Max branching factor for finding matches
    //static const int MAX_BF=6;
    typedef std::vector<int> TCloseIdSet;
}

#endif // _TRANSFORM
