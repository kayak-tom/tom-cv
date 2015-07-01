#ifndef POINTWITHTHICKNESS_H
#define POINTWITHTHICKNESS_H

#include "lines.h"

#define THROW_NO_THICKNESS THROW("Functions using polyline thicknesses should only be called for polylines with thicknesses (and functions without thicknesses should not call functions with)")

class CConstructFromLMVector
{

};

template<class TLineType_in>
class CPolylineControlPoint
{
public:
    typedef TLineType_in TLineType;
    typedef typename TLineType::TVecType TVecType;
    enum ePolylineParamIndices { ePolyX=0, ePolyY=1, ePolyZ=(TVecType::RowsAtCompileTime==3 ? 2 : -1), PARAMS_PER_POLY_CONTROL_POINT=TVecType::RowsAtCompileTime, ePolyWidth=-1 }; //Used in LM optimisations
private:
    TVecType x; //World or image pixel coordinates, depending on TVecType
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    /*friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        const char * names[] = {"X", "Y", "Z", "\0"};
        CHECK(x.size() > 3, "Point dim error");

        for (int i = 0; i < x.size(); i += 1) {
            ar & boost::serialization::make_nvp(names[i], x[i]);
        }
    }*/

    CPolylineControlPoint<TLineType> scale(const double s) const { return CPolylineControlPoint<TLineType>(x*s); }

    CPolylineControlPoint() {}
    CPolylineControlPoint(const TVecType & x/*, const CConstructFromPoint &*/) : x(x) {
        checkBadNum();
    }
    CPolylineControlPoint(const TVecType & x, const CConstructFromLMVector &) : x(x) {
        checkBadNum();
    }

    CPolylineControlPoint(const TVecType &, const double) {
        THROW_NO_THICKNESS;
    }

    static TVecType uninit() {
        return TVecType::uninit();
    }
    static CPolylineControlPoint uninitCP() {
        return CPolylineControlPoint(uninit());
    }

    TVecType & getPoint() {//'Ref' is important as g2o optimiser refers back to polyline points.
        return x;
    }
    const TVecType & getPoint() const {//'Ref' is important as g2o optimiser refers back to polyline points.
        return x;
    }

    TVecType asVector() const {
        return x;
    }

    static const bool HAS_THICKNESS = false;

    const double & getWidth() const {
        THROW_NO_THICKNESS;
    }
    static const double minWidth() {
        THROW_NO_THICKNESS;
    }
    const double getRad() const {
        THROW_NO_THICKNESS;
    }

    void addLocalisationErrors(const double dGaussianNoiseSD);

    double getWidth_default(const double dDefault) const {
        return (dDefault <= 0) ? 1 : dDefault;
    }
    
    double distanceToPoint(const TVecType & vec) const;
    
    std::ostream & pp(std::ostream & s) const;

    void checkBadNum() const;
};

template<class TLineType_in>
class CPolylineControlPointWithThickness //Don't inherit from CPolylineControlPoint to avoid any poor alignment issues (may not matter)
{
public:
    typedef TLineType_in TLineType;
    typedef typename TLineType::TVecType TVecType;
    enum ePolylineParamIndices { ePolyX=0, ePolyY=1, ePolyZ=(TVecType::RowsAtCompileTime==3 ? 2 : -1), ePolyWidth=TVecType::RowsAtCompileTime, PARAMS_PER_POLY_CONTROL_POINT=TVecType::RowsAtCompileTime+1 }; //Used in LM optimisations
    typedef Eigen::Matrix<double, PARAMS_PER_POLY_CONTROL_POINT, 1> TLMVector;
private:
    TVecType x; //World or image pixel coordinates, depending on TVecType
    double dWidth; //Metres or pixels, depending on TVecType
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    /*friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        const char * names[] = {"X", "Y", "Z", "\0"};
        CHECK(x.size() > 3, "Point dim error");

        for (int i = 0; i < x.size(); i += 1) {
            ar & boost::serialization::make_nvp(names[i], x(i));
        }
        ar & boost::serialization::make_nvp("W", dWidth);
    }*/

    CPolylineControlPointWithThickness() {}
    CPolylineControlPointWithThickness(const TVecType & x/*, const CConstructFromPoint &*/) {
        THROW_NO_THICKNESS;
    }
    CPolylineControlPointWithThickness(const TLMVector & lmVec, const CConstructFromLMVector &);
    CPolylineControlPointWithThickness(const TVecType & x, const double dWidth);

    static TVecType uninit() {
        return TVecType::uninit();
    }
    static CPolylineControlPointWithThickness uninitCP() {
        return CPolylineControlPointWithThickness(uninit(), 0.1);
    }

    CPolylineControlPointWithThickness<TLineType> scale(const double s) const { return CPolylineControlPointWithThickness<TLineType>(x*s, dWidth*s); }

    TVecType & getPoint() {//'Ref' is important as g2o optimiser refers back to polyline points.
        return x;
    }
    const TVecType & getPoint() const {//'Ref' is important as g2o optimiser refers back to polyline points.
        return x;
    }

    const double & getWidth() const {//'Ref' is important for 'serialisable'.
        if(IS_DEBUG) {
            CHECKBADNUM(dWidth);
            CHECK(dWidth<=0, "Negative width");
        }
        return dWidth;
    }
    double & getWidth() { //'Ref' is important as g2o optimiser refers back to polyline points.
        return dWidth;
    }
    const double getRad() const {
        return 0.5*getWidth();
    }

    TLMVector asVector() const {
        TLMVector lmVec;
        lmVec.template head<TVecType::RowsAtCompileTime>() = x;
        lmVec.template tail<1>()(0) = dWidth;
        return lmVec;
    }

    static const bool HAS_THICKNESS = true;

    void addLocalisationErrors(const double dGaussianNoiseSD);

    double getWidth_default(const double dDefault) const {
        return (dDefault<=0) ? getWidth() : dDefault;
    }

    double distanceToPoint(const TVecType & vec) const;

    std::ostream & pp(std::ostream & s) const;

    void checkBadNum() const;
};

typedef CPolylineControlPoint<C2dBoundedLine> T2dPolylineControlPoint;
typedef CPolylineControlPointWithThickness<C2dBoundedLine> T2dPolylineControlPointWithThickness;
typedef CPolylineControlPoint<C3dBoundedLine> T3dPolylineControlPoint;
typedef CPolylineControlPointWithThickness<C3dBoundedLine> T3dPolylineControlPointWithThickness;

class C2dPolylineControlPoint : public T2dPolylineControlPoint
{
    typedef T2dPolylineControlPoint Base;
public:
    C2dPolylineControlPoint() {}
    C2dPolylineControlPoint(const Base::TVecType & x) : T2dPolylineControlPoint(x) {
    }
    C2dPolylineControlPoint(const Base::TVecType & x, const double dWidth) : T2dPolylineControlPoint(x, dWidth) {
    }
    C2dPolylineControlPoint(const Base::TVecType & x, const CConstructFromLMVector & lm_unused) : T2dPolylineControlPoint(x, lm_unused) {
    }

    C2dPolylineControlPoint(const Base & cp) : Base(cp) {  }

    //using T2dPolylineControlPoint::uninit;

    void draw(cv::Mat & image, const cv::Scalar & col, const int nRad=4, const bool bDrawX=false) const;
};

class C2dPolylineControlPointWithThickness : public T2dPolylineControlPointWithThickness
{
    typedef T2dPolylineControlPointWithThickness Base;
public:
    C2dPolylineControlPointWithThickness() {}
    C2dPolylineControlPointWithThickness(const Base::TVecType & x) : T2dPolylineControlPointWithThickness(x) {
    }
    C2dPolylineControlPointWithThickness(const Base::TVecType & x, const double dThickness) : T2dPolylineControlPointWithThickness(x, dThickness) {
    }
    C2dPolylineControlPointWithThickness(const Base::TLMVector & x, const CConstructFromLMVector & lm_unused) : T2dPolylineControlPointWithThickness(x, lm_unused) {
    }

    C2dPolylineControlPointWithThickness(const Base & cp) : Base(cp) {}

    //using T2dPolylineControlPointWithThickness::uninit;

    void draw(cv::Mat & image, const cv::Scalar & col, const bool bDrawX=false) const;
};


class C3dPolylineControlPoint : public T3dPolylineControlPoint
{
    typedef T3dPolylineControlPoint Base;
public:
    typedef C2dPolylineControlPoint T2dPoint;

    C3dPolylineControlPoint() {}
    C3dPolylineControlPoint(const Base::TVecType & x) : T3dPolylineControlPoint(x) {
    }
    C3dPolylineControlPoint(const Base::TVecType & x, const double dWidth) : T3dPolylineControlPoint(x, dWidth) {
    }
    C3dPolylineControlPoint(const Base::TVecType & x, const CConstructFromLMVector & lm_unused) : T3dPolylineControlPoint(x, lm_unused) {
    }

    C3dPolylineControlPoint(const Base & cp) : Base(cp) {  }

    //using T3dPolylineControlPoint::uninit;
    //static C3dPolylineControlPoint uninit() { return T3dPolylineControlPoint::uninit(); }

    void draw(cv::Mat & image, const CWorldCamera & P, const cv::Scalar & col, const double dRad=0.01, const bool bDrawX=false) const;
};

class C3dPolylineControlPointWithThickness : public T3dPolylineControlPointWithThickness
{
    typedef T3dPolylineControlPointWithThickness Base;
public:
    typedef C2dPolylineControlPointWithThickness T2dPoint;

    C3dPolylineControlPointWithThickness() {}
    C3dPolylineControlPointWithThickness(const Base::TVecType & x) : T3dPolylineControlPointWithThickness(x) {
    }
    C3dPolylineControlPointWithThickness(const Base::TVecType & x, const double dThickness) : T3dPolylineControlPointWithThickness(x, dThickness) {
    }
    C3dPolylineControlPointWithThickness(const Base::TLMVector & x, const CConstructFromLMVector & lm_unused) : T3dPolylineControlPointWithThickness(x, lm_unused) {
    }
    C3dPolylineControlPointWithThickness(const Base & cp) : Base(cp) { }
    //C3dPolylineControlPointWithThickness(const Base & cp) { *this = cp; }

    //using T3dPolylineControlPointWithThickness::uninit;

    void draw(cv::Mat & image, const CWorldCamera & P, const cv::Scalar & col, const bool bDrawX=false) const;

    double volume() const;

    C3dWorldPoint closestPointOnSurface(const C3dWorldPoint & p) const;
};

template<class TControlPoint>
double distanceBetween(const TControlPoint & cp1, const TControlPoint & cp2)
{
    double dDist = (cp1.getPoint()-cp2.getPoint()).norm();
    if(!TControlPoint::HAS_THICKNESS)
        return dDist;
        
    dDist -= cp1.getRad();
    dDist -= cp2.getRad();
    
    return std::max<double>(0, dDist);
}


template<class TLineType_in>
void CPolylineControlPointWithThickness<TLineType_in>::checkBadNum() const
{
    for(int i=0; i<x.size(); i++)
        CHECKBADNUM(x(i));
    CHECKBADNUM(dWidth);
}


template<class TLineType_in>
void CPolylineControlPoint<TLineType_in>::checkBadNum() const
{
    for(int i=0; i<x.size(); i++)
        CHECKBADNUM(x(i));
}

#endif //POINTWITHTHICKNESS_H
