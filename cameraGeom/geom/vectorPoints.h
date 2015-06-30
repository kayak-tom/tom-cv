#ifndef VECTORPOINTS_H
#define VECTORPOINTS_H

//#include <boost/serialization/serialization.hpp>

typedef Eigen::Matrix<double, 3, 1> TEigen3dPoint;
typedef Eigen::Matrix<double, 4, 1> TEigen3dPointHomog;
typedef Eigen::Matrix<double, 3, 3> TEigen3dRotation;
typedef Eigen::Matrix<double, 3, 1> TEigen2dPointHomog;
typedef Eigen::Matrix<double, 2, 1> TEigen2dPoint;
typedef Eigen::Matrix<double, 3, 4, Eigen::RowMajor> TCamMat;
typedef Eigen::Matrix<double, 3, 3> TFMat;

/*define constructors etc. required when inheriting from Eigen*/
#define INHERIT_FROM_EIGEN_VECTOR(TClassName,TVecType) \
    public: \
    typedef TVecType Base; \
    void checkInit() const { if(IS_DEBUG) CHECK(squaredNorm() > HUGE, "vector is uninitialised"); checkForNan(); } \
    void checkForNan() const { CHECKNAN(sum()); } \
    \
    TClassName() : Base() { \
        setConstant(0); \
    } \
    \
    TClassName(const Base & M) : Base(M) { \
        checkForNan(); \
    }\
    \
    TClassName(const TVecType ## Homog & vhomog) { \
        *this = vhomog.segment(0, Base::RowsAtCompileTime) / vhomog(Base::RowsAtCompileTime); \
        checkForNan(); \
    }\
    \
    template<typename OtherDerived> \
    TClassName & operator=(const Eigen::MatrixBase <OtherDerived>& other) { \
        this->Base::operator=(other); \
        checkForNan(); \
        return *this; \
    } \
    \
    TVecType ## Homog homog() const { \
        TVecType ## Homog vhomog; \
        vhomog << ((Base &)*this), 1; \
        checkForNan(); \
        return vhomog; \
    }\
    template<class OtherDerived>\
    TClassName(const Eigen::EigenBase<OtherDerived > & other) : TVecType(other)  {} \
    static TClassName uninit() { return Constant(-10000); }
    
class C2dImagePointPx;

class C3dWorldPoint : public TEigen3dPoint
{
    INHERIT_FROM_EIGEN_VECTOR(C3dWorldPoint, TEigen3dPoint);
public:
    typedef C2dImagePointPx T2dPoint;//3D points project to C2dImagePointPx

    C3dWorldPoint(const double x, const double y, const double z);

    //MAKE_SERIALISABLE;
};

//bool depthIsReasonable(const double dZ, const double dMargin=0);

typedef std::vector<C3dWorldPoint,Eigen::aligned_allocator<C3dWorldPoint> > T3dPointVector;
typedef std::vector<TEigen3dPoint,Eigen::aligned_allocator<TEigen3dPoint> > T3dVecVector;


//2d point in normalised homogeneous image coordinates
class C2dImagePoint : public TEigen2dPoint
{
    INHERIT_FROM_EIGEN_VECTOR(C2dImagePoint, TEigen2dPoint);
private:
    //void operator=(const C2dImagePointPx &); //Do not allow copying between normalised and unnormalised image points
    C2dImagePoint(const C2dImagePointPx &);
public:
    C2dImagePoint(const double x, const double y);
};

inline void testAssignmentsCompilingOk()
{
    C2dImagePoint x,y;
    const C2dImagePoint z = x+y;
    //C2dImagePoint x2=2.0*z;
    cout << z;
}

typedef std::vector<C2dImagePoint,Eigen::aligned_allocator<C2dImagePoint> > T2dImPointVector;
#define TGenericPointVector std::vector<T,Eigen::aligned_allocator<T> >

//2d point in pixel coordinates.
class C2dImagePointPx : public TEigen2dPoint
{
    INHERIT_FROM_EIGEN_VECTOR(C2dImagePointPx, TEigen2dPoint);
private:
    void operator=(const C2dImagePoint &); //Do not allow copying between normalised and unnormalised image points
    C2dImagePointPx(const C2dImagePoint &);
public:
    C2dImagePointPx(const cv::Point2f & p) : TEigen2dPoint((double)p.x,(double)p.y) 
    {
    }

    C2dImagePointPx(const cv::Point & p) : TEigen2dPoint((double)p.x,(double)p.y) 
    {
    }
        
    C2dImagePointPx(const double x, const double y) : TEigen2dPoint(x,y) 
    {
    }

    operator cv::Point () const;

    void operator= (const cv::Point & p) {
        (*this)(0)=p.x;
        (*this)(1)=p.y;
    }
};

typedef std::vector<C2dImagePointPx, Eigen::aligned_allocator<C2dImagePointPx> > T2dImPointPxVector;


#endif // CALIBRATEROBOT_H
