/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "geom.h"
#include <iomanip>

#ifndef USE_MATH_DEFINES
#define USE_MATH_DEFINES
#endif
#include <math.h>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include "util/exception.h"
#include "util/dynArray.h"
#include "util/random.h"
#include <boost/smart_ptr.hpp>

using namespace std;

#define MAT_FROM_QUAT(Mfn) \
Mfn(1,1) = 1 - 2*sqr(Q[1]) - 2*sqr(Q[2]);      Mfn(1,2) = 2*Q[0]*Q[1] - 2*Q[2]*Q[3];      Mfn(1,3) = 2*Q[0]*Q[2] + 2*Q[1]*Q[3]; \
Mfn(2,1) = 2*Q[0]*Q[1] + 2*Q[2]*Q[3];         Mfn(2,2) = 1 - 2*sqr(Q[0]) - 2*sqr(Q[2]);    Mfn(2,3) = 2*Q[1]*Q[2] - 2*Q[0]*Q[3]; \
Mfn(3,1) = 2*Q[0]*Q[2] - 2*Q[1]*Q[3];         Mfn(3,2) = 2*Q[1]*Q[2] + 2*Q[0]*Q[3];     Mfn(3,3) = 1 - 2*sqr(Q[0]) - 2*sqr(Q[1]);

#define ACCESS_CAM_AS_MAT(r, c) cam.adCam[(r-1)*4 + (c-1)]
#define ACCESS_EIGENMAT_AS_MAT(r, c) rot(r-1, c-1)

void C3dRotationQuat::asMat(Eigen::Matrix3d & rot) const //Windows requires full template param spec
{
    MAT_FROM_QUAT(ACCESS_EIGENMAT_AS_MAT);
}

#ifdef NEW_VERSION_BELOW
//Homography decomposition my Levenberg-Marquardt
//H=s(R - t n')
//Todo: find closed-form soln.
const int NUM_PARAMS = 11, NUM_MODEL_VARS = 9;
typedef Eigen::Matrix<double, NUM_PARAMS, 1 > TParamVector;
typedef Eigen::Matrix<double, NUM_MODEL_VARS, 1 > TResidVector;
typedef Eigen::Matrix<double, NUM_MODEL_VARS, NUM_PARAMS> TJMatrix;
typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;

void f(const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, const double s, Eigen::Matrix3d HAtParams) {
    Eigen::Vector3d t, n;
    t << camMotion.getX(), camMotion.getY(), camMotion.getZ();
    n << planeNormal.getX(), planeNormal.getY(), planeNormal.getZ();
    Eigen::Matrix3d rot;
    rotation.asMat(rot);
    HAtParams = s * (rot - t * n.transpose());
}

void getParamVector(const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, const double s, TParamVector & params) {
    for (int i = 0; i < 4; i++)
        params[i] = rotation[i];

    params[4] = planeNormal.getX();
    params[5] = planeNormal.getY();
    params[6] = planeNormal.getZ();

    params[7] = camMotion.getX();
    params[8] = camMotion.getY();
    params[9] = camMotion.getZ();

    params[10] = s;
}

void setParams(const TParamVector & params, C3dRotation & rotation, C3dPoint & planeNormal, C3dPoint & camMotion, double & s) {
    rotation = C3dRotation(params[0], params[1], params[2], params[3]); //Normalisation normally redundent
    planeNormal = C3dPoint(params[4], params[5], params[6]);
    planeNormal.normalise(); //Normalisation normally redundent

    camMotion = C3dPoint(params[7], params[8], params[9]);
    s = params[10];
}

void getJAndResidAtParams(const Eigen::Matrix3d & H, const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, const double s, const double epsilon, TJMatrix & J, TResidVector & residual) {
    TParamVector paramVec, tempParamVec;
    getParamVector(rotation, planeNormal, camMotion, s, paramVec);
    Eigen::Matrix3d HAtParams, HAtParamsPlusEps;
    f(rotation, planeNormal, camMotion, s, HAtParams);

    cout << HAtParams << endl << "=HAtParams\n";

    //Now vary each param and find deriv.
    for (int nParam = 0; nParam < NUM_MODEL_VARS; nParam++) {
        TParamVector tempParamVec = paramVec;
        tempParamVec(nParam) += epsilon;
        C3dRotation tempRotation;
        C3dPoint tempPlaneNormal;
        C3dPoint tempCamMotion;
        double temps;
        setParams(tempParamVec, tempRotation, tempPlaneNormal, tempCamMotion, temps);
        f(tempRotation, tempPlaneNormal, tempCamMotion, temps, HAtParamsPlusEps);

        const Eigen::Matrix3d derivs = (HAtParamsPlusEps - HAtParams) / epsilon;
        for (int nVar = 0; nVar < NUM_MODEL_VARS; nVar++) {
            J(nVar, nParam) = derivs(nVar / 3, nVar % 3);
        }
    }

    const Eigen::Matrix3d resids = (H - HAtParams);
    for (int nVar = 0; nVar < NUM_MODEL_VARS; nVar++) {
        residual(nVar) = resids(nVar / 3, nVar % 3);
    }
}

//Adjust 10 parameters to refine 9 Eigen::Matrix elements

void decomposeHomography(const Eigen::Matrix3d & H, C3dRotation & rotation, C3dPoint & planeNormal, C3dPoint & camMotion) {
    rotation = C3dRotation(); //Identity
    planeNormal = C3dPoint(0, 0, -1);
    camMotion = C3dPoint();
    double s = 1;

    const int MAX_ITERS = 25;

    TParamVector params;

    for (int nIter = 0; nIter < MAX_ITERS; nIter++) {
        const double lambda = 0.25;
        const double epsilon = 0.0001;

        getParamVector(rotation, planeNormal, camMotion, s, params);

        TJMatrix J;
        TResidVector residual;
        getJAndResidAtParams(H, rotation, planeNormal, camMotion, s, epsilon, J, residual);

        TJTJMatrix JTJ = J.transpose() * J;
        JTJ.diagonal() *= (1 + lambda);

        const TJTJMatrix inv = JTJ.inverse();

        const TParamVector paramUpdateVec = inv * J.transpose() * residual;

        params += paramUpdateVec;

        setParams(params, rotation, planeNormal, camMotion, s);

        cout << "Error = " << residual.squaredNorm() << endl;
        cout << "s = " << s << endl;
        cout << "n = " << planeNormal << endl;
        cout << "t = " << camMotion << endl;
        cout << "R = " << rotation << endl;
    }
}
#else
//Homography decomposition by Levenberg-Marquardt
//H=s(R - t n')
//Todo: find closed-form soln.
const int NUM_PARAMS = 10, NUM_MODEL_VARS = 9;
typedef Eigen::Matrix<double, NUM_PARAMS, 1 > TParamVector;
typedef Eigen::Matrix<double, NUM_MODEL_VARS, 1 > TResidVector;
typedef Eigen::Matrix<double, NUM_MODEL_VARS, NUM_PARAMS> TJMatrix;
typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;

void makeH(const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, Eigen::Matrix3d & HAtParams) {
    Eigen::Vector3d t, n;
    t << camMotion.getX(), camMotion.getY(), camMotion.getZ();
    n << planeNormal.getX(), planeNormal.getY(), planeNormal.getZ();
    Eigen::Matrix3d rot;
    rotation.asMat(rot);
    HAtParams = (rot - t * n.transpose());
}

void getParamVector(const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, TParamVector & params) {
    for (int i = 0; i < 4; i++)
        params[i] = rotation[i];

    params[4] = planeNormal.getX();
    params[5] = planeNormal.getY();
    params[6] = planeNormal.getZ();

    params[7] = camMotion.getX();
    params[8] = camMotion.getY();
    params[9] = camMotion.getZ();
}

void setParams(const TParamVector & params, C3dRotation & rotation, C3dPoint & planeNormal, C3dPoint & camMotion) {
    rotation = C3dRotation(params[0], params[1], params[2], params[3]); //Normalisation normally redundent, which is ok
    planeNormal = C3dPoint(params[4], params[5], params[6]);
    planeNormal.normalise(); //Normalisation normally redundent

    camMotion = C3dPoint(params[7], params[8], params[9]);
}

void getJAndResidAtParams(const Eigen::Matrix3d & H, const C3dRotation & rotation, const C3dPoint & planeNormal, const C3dPoint & camMotion, const double epsilon, TJMatrix & J, TResidVector & residual) {
    TParamVector paramVec, tempParamVec;
    getParamVector(rotation, planeNormal, camMotion, paramVec);
    Eigen::Matrix3d HAtParams, HAtParamsPlusEps;
    makeH(rotation, planeNormal, camMotion, HAtParams);

    //cout << HAtParams << endl << "=HAtParams\n";

    //Now vary each param and find deriv.
    for (int nParam = 0; nParam < NUM_PARAMS; nParam++) {
        TParamVector tempParamVec = paramVec;
        tempParamVec(nParam) += epsilon;
        C3dRotation tempRotation;
        C3dPoint tempPlaneNormal;
        C3dPoint tempCamMotion;
        setParams(tempParamVec, tempRotation, tempPlaneNormal, tempCamMotion);
        makeH(tempRotation, tempPlaneNormal, tempCamMotion, HAtParamsPlusEps);

        //cout << HAtParamsPlusEps-HAtParams << endl << "=HAtParamsPlusEps-HAtParams" << endl;
        const Eigen::Matrix3d derivs = (HAtParamsPlusEps - HAtParams) / epsilon;
        for (int nVar = 0; nVar < NUM_MODEL_VARS; nVar++) {
            J(nVar, nParam) = derivs(nVar / 3, nVar % 3);
        }
    }

    const Eigen::Matrix3d resids = (H - HAtParams);
    for (int nVar = 0; nVar < NUM_MODEL_VARS; nVar++) {
        residual(nVar) = resids(nVar / 3, nVar % 3);
    }
}

/* Decompose a homography Eigen::Matrix into a rotation, translation, and plane normal, by Levenberg least-squares optimisation.

   H = (R - t * n^t)s

   |n|==1

   Note SIGN AMBIGUITY--same solution if you flip signs of n and t.

   Note DEGENERATE if t==0. I think the rotation is still useful.

   If d is distance to plane, t*d is camera translation vector.

   Returns true on success.

   Uses SVD to find s (2nd s.v.). n and R^-1 * t are perpendicular to 2nd col of V, this is *not currently used*

   Closed form solution--TODO!
 */
//Adjust 10 parameters to refine 9 Eigen::Matrix elements. Need to use Levenberg algorithm because Lev-Mar assumes derivatives nonzero somewhere (todo: I'm sure there's a proper way of coping with this)

bool decomposeHomography(const Eigen::Matrix3d & H_in, C3dRotation & rotation, C3dPoint & planeNormal, C3dPoint & camMotion) {
    rotation = C3dRotation(); //Identity
    planeNormal = C3dPoint(0, 0, -1);
    camMotion = C3dPoint();

#if EIGEN_VERSION_AT_LEAST(2,90,0)
    Eigen::JacobiSVD< Eigen::Matrix3d > svdH(H_in);
#else
    Eigen::SVD< Eigen::Matrix3d > svdH(H_in);
#endif

    const double s = svdH.singularValues()(1);
    const Eigen::Matrix3d H = H_in / s;

    //cout << H;
    cout << "s = " << svdH.singularValues()(1) << endl;

    const int MAX_ITERS = 200;

    TParamVector params;

    double lambda = 0.025, dErr = HUGE;
    const double eps = sqr(0.000001);
    int nIter = 0;
    for (; nIter < MAX_ITERS && dErr > eps && dErr <= HUGE /*breakout on overflow. Needed??*/; nIter++) {
        const double epsilon = 0.0001;

        getParamVector(rotation, planeNormal, camMotion, params);

        TJMatrix J;
        TResidVector residual;
        getJAndResidAtParams(H, rotation, planeNormal, camMotion, epsilon, J, residual);

        TJTJMatrix JTJ = J.transpose() * J;
        //JTJ.diagonal() *= (1+lambda);
        JTJ.diagonal().array() += lambda;

        const TJTJMatrix inv = JTJ.inverse();

        TParamVector paramUpdateVec = inv * J.transpose() * residual;

        params += paramUpdateVec;

        setParams(params, rotation, planeNormal, camMotion);

        dErr = residual.squaredNorm();

        lambda *= 0.6; //Unsophisticated, TODO
    }

    if (dErr > eps) {
        cout << "Homography convergence failed after " << nIter << " iterations" << endl;
        cout << "Error = " << dErr << endl;
        cout << "n = " << planeNormal << endl;
        cout << "t = " << camMotion << endl;
        cout << "R = " << rotation << endl;
        return false;
    }
    cout << "Homography convergence after " << nIter << " iterations, ";
    cout << "error = " << dErr << endl;
    cout << "n = " << planeNormal << endl;
    cout << "t = " << camMotion << endl;
    cout << "R = " << rotation << endl;

    return true;
}

#endif

inline void doubleMemCopy(int n, const double * pdSrc, double * pdDest) {
    for (int ii = n; ii > 0; ii--) {
        *pdDest = *pdSrc;
        pdSrc++;
        pdDest++;
    }
}

template<int R, int C>
inline void doubleMemCopy(const Eigen::Matrix<double, R, C> & mat, double * pdDest) {
    //    int n = mat.rows()*mat.cols();
    //    for(int ii=n; ii>0; ii--)
    //    {
    //        *pdDest = *pdSrc; pdSrc++; pdDest++;
    //    }
    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            *pdDest = mat(r, c);
            pdDest++;
        }
}

template<int R, int C>
inline void doubleMemCopyTrans(const Eigen::Matrix<double, C, R> & mat, double * pdDest) {
    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            *pdDest = mat(c, r);
            pdDest++;
        }
}

int dist(CLocation l1, CLocation l2) {
    double dx = l1.dx() - l2.dx();
    double dy = /*3*/(l1.dy() - l2.dy());
    return (int) ((dx * dx + dy * dy));
}

void CCamera::calibrate(const CCamCalibMatrix & calib) {
    matMul < 3, 3, 4 > (calib.adK, adCam, adCam);
    /*Eigen::Matrix3d K_t(calib.adK);
    Eigen::Matrix<double, 4, 3> Cam_t(adCam);
    cout << K << "=K\n";
    cout << Cam << "=cam\n";
    //Cam = K * Cam;
    Cam *= K;
    cout << Cam << "=calib cam\n";

    doubleMemCopy<3,4>(Cam, adCam);*/
}

//Convert image coordinates to world
C3dPoint CCamera::imageToWorld(const C2dPoint & imPoint, const double dDepth) const
{
    const double x=imPoint.getX();
    const double y=imPoint.getY();
    const double det = (-at1(1, 1)*at1(3, 2)*y+at1(1, 1)*at1(2, 2)-at1(1, 2)*at1(2, 1)+at1(1, 2)*at1(3, 1)*y+x*at1(2, 1)*at1(3, 2)-x*at1(3, 1)*at1(2, 2));
    const double d=dDepth;
    const double X = (at1(1, 2)*at1(2, 3)*d-at1(1, 3)*d*at1(2, 2)+at1(1, 2)*at1(2, 4)-at1(1, 4)*at1(2, 2)+at1(3, 3)*d*x*at1(2, 2)-at1(2, 3)*d*x*at1(3, 2)+y*at1(1, 3)*d*at1(3, 2)-at1(3, 3)*d*at1(1, 2)*y+at1(3, 4)*x*at1(2, 2)-at1(2, 4)*x*at1(3, 2)+y*at1(1, 4)*at1(3, 2)-at1(3, 4)*at1(1, 2)*y)  /det;
    const double Y = -(-at1(1, 4)*at1(2, 1)+at1(2, 1)*at1(3, 4)*x+at1(2, 1)*at1(3, 3)*d*x-at1(1, 3)*d*at1(2, 1)-at1(2, 3)*d*at1(3, 1)*x+at1(3, 1)*y*at1(1, 3)*d+at1(1, 1)*at1(2, 4)-at1(3, 4)*at1(1, 1)*y-at1(2, 4)*at1(3, 1)*x-at1(3, 3)*d*at1(1, 1)*y+at1(1, 1)*at1(2, 3)*d+at1(3, 1)*y*at1(1, 4))/det;
    
    C3dPoint X3d(X,Y,dDepth);
    
    if(IS_DEBUG)
    {
        C2dPoint check = X3d.photo(*this);
        if(IS_DEBUG) CHECK((check-imPoint).sum_square() > 0.0001, "imageToWorld failed");
    }
    
    return X3d;
}
//Points must be CALIBRATED

/*void correctRD(Matrix & points, const double * adRD, int nRDcoeffs)
{
    if(IS_DEBUG) CHECK(points.ncols() == 2,"correctRD(): Need 3d homo centred points atm");

    int numPoints = points.ncols();
    for(int n=1; n <= numPoints; n++)
    {
//        double dRadius = sqrt(points.column(n).sum_square());
        double x=points(1, n), y=points(2, n);
        double dRadius_sq = (SQR(x)+SQR(y));
        double dRadCorrectionFactor = radCorrectionFactor(dRadius_sq, adRD, nRDcoeffs);
        //points.column(n) *= dRadCorrectionFactor;
        points(1, n) = x * dRadCorrectionFactor;
        points(2, n) = y * dRadCorrectionFactor;
    }
}*/
C2dPoint::C2dPoint(const Eigen::Vector2d & v) : CSimple2dPoint(v(0), v(1)) {
    if(IS_DEBUG) CHECK(std::isnan(x + y), "C2dPoint: NaN");
}

C2dPoint::C2dPoint(const Eigen::Vector3d & v) : CSimple2dPoint(v(0), v(1)) {
    double inv = 1.0 / v(2);
    x *= inv;
    y *= inv;
    if(IS_DEBUG) CHECK(std::isnan(x + y) || std::isinf(x + y), "C2dPoint: NaN/inf");
}

C3dPoint::C3dPoint(const Eigen::Vector3d & v) : x(v(0)), y(v(1)), z(v(2)) {
    if(IS_DEBUG) CHECK(std::isnan(x + y + z), "C3dPoint: NaN");
}

C3dPoint::C3dPoint(const Eigen::Vector3f & v) : x(v(0)), y(v(1)), z(v(2)) {
    if(IS_DEBUG) CHECK(std::isnan(x + y + z), "C3dPoint: NaN");
}

void C2dPoint::asVector(Eigen::Vector2d & vec) const {
    vec(0) = x;
    vec(1) = y;
}

void C2dPoint::asVector(Eigen::Vector3d & vec) const {
    vec(0) = x;
    vec(1) = y;
    vec(2) = 1;
}

void C3dPoint::asVector(Eigen::Vector3d & vec) const {
    vec(0) = x;
    vec(1) = y;
    vec(2) = z;
}

/*int roundToInt(double d)
{
    int i=(int)d;
    if(d-(double)i>0.5) i++;
    return i;
}
CLocation colVecToLoc(const C3dPoint & v)
{
    double dScale = 1.0/v.getZ();
#if EXACT_HACK > 1
    return CLocation(doubleToInt(EXACT_HACK*v.getX()*dScale), doubleToInt(EXACT_HACK*v.getY()*dScale));
#else
    return CLocation(doubleToInt(v.getX()*dScale), doubleToInt((v.getY()*dScale)));
#endif
}*/

CLocation colVecToLoc(const C2dPoint & v) {
#if EXACT_HACK > 1
    return CLocation(doubleToInt(EXACT_HACK * v.getX()), doubleToInt(EXACT_HACK * v.getY()));
#else
    return CLocation(doubleToInt(v.getX()), doubleToInt(v.getY()));
#endif
}

//Reconstruct 3d pos of a point

C3dPoint reconstruct(const CCamera &P, const CCamera &Pp, const C2dPoint &p, const C2dPoint &pp) {
    //typedef Eigen::RowVector4d Vec4d;
    typedef Eigen::Matrix<double, 1, 4, Eigen::RowMajor + Eigen::AutoAlign> Vec4d;

    const Vec4d Prow1(P.rowData(0));
    const Vec4d Prow2(P.rowData(1));
    const Vec4d Prow3(P.rowData(2));
    const Vec4d Pprow1(Pp.rowData(0));
    const Vec4d Pprow2(Pp.rowData(1));
    const Vec4d Pprow3(Pp.rowData(2));

    Eigen::Matrix4d A;
    A.row(0) = p.getX() * Prow3 - Prow1;
    A.row(1) = p.getY() * Prow3 - Prow2;
    A.row(2) = pp.getX() * Pprow3 - Pprow1;
    A.row(3) = pp.getY() * Pprow3 - Pprow2;

#if EIGEN_VERSION_AT_LEAST(2,90,0)
    Eigen::JacobiSVD< Eigen::Matrix4d > svdA(A, Eigen::ComputeFullV);// speed-up by not computing U
#else
    Eigen::SVD< Eigen::Matrix4d > svdA(A);
#endif

    const Eigen::Matrix4d & V = svdA.matrixV();

    
    C3dPoint p3d(V(0, 3), V(1, 3), V(2, 3));
    p3d /= V(3, 3);

    //cout << p3d.depth(P) << "-";
    //The SVs tell us very little cout << svdA.singularValues().transpose() << " ";

    //C3dPoint p3d_deeper = p3d * 1.1;
    //C2dPoint newx = p3d_deeper.photo(P); this would be the same
    //C2dPoint newpp = p3d_deeper.photo(Pp);
    //cout << (pp-newpp).sum_square() << endl;
    
    return p3d;
}

C2dPoint operator*(const CCamera & P, const C3dPoint & p) {
    const double * R = P.adCam;
    double xNew = R[0] * p.getX() + R[1] * p.getY() + R[2] * p.getZ() + R[3];
    double yNew = R[4] * p.getX() + R[5] * p.getY() + R[6] * p.getZ() + R[7];
    double zNew = R[8] * p.getX() + R[9] * p.getY() + R[10] * p.getZ() + R[11];
    return C2dPoint(C3dPoint(xNew, yNew, zNew));
};

C2dPoint C3dPoint::photo(const CCamera & P) const {
    return P * * this;
}

bool C3dPoint::testInFront(const CCamera & P) const {
    const double * R = P.adCam;
    double zNew = R[8] * x + R[9] * y + R[10] * z + R[11];

    bool bInFront = zNew > 0;

    if(IS_DEBUG) CHECK((depth(P) > 0) != bInFront, "Project: Depth inconsistency");
    return bInFront;
}

double det3d(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return a * e * i - a * f * h - b * d * i + b * f * g + c * d * h - c * e*g;
}

double C3dPoint::depth(const CCamera & P) const {
    const double * R = P.adCam;

    double m3_ss_inv = 1.0/sqrt(sqr(R[8]) + sqr(R[9]) + sqr(R[10])); //(top of p161) Todo: normally sqrting 1 here
    if(IS_DEBUG) CHECK(!zero(sqr(R[8]) + sqr(R[9]) + sqr(R[10]) - 1), "Depth computation will fail")

            double sdM =
#ifdef _DEBUG
            sign(det3d(R[0], R[1], R[2], R[4], R[5], R[6], R[8], R[9], R[10])); //=sign(M.determinant());
    if (sdM == -1)
        std::cout << P;
    if(IS_DEBUG) CHECK(sdM == -1, "Not a problem, but if this is never -1 we don't need to compute it");
#else
            1;
#endif

    double zProj = R[8] * x + R[9] * y + R[10] * z + R[11];
    return (sdM * zProj * m3_ss_inv);
}

/*
 * UNRELIABLE
bool testDepthOk(const Eigen::Matrix &P, const ColumnVector &Q)
{
    if(IS_DEBUG) CHECK(!zero(Q(4)-1), "testDepthOk: Coord not homog.");
    ColumnVector m = P.submatrix(3,3,1,3).t(); //p.row(3) UNUSED
    ColumnVector p = P*Q;
    return (p(3) >= 0);
    //P.submatrix(1,3,1,3)
}*/

bool testPair(const CCamera &P, const CCamera &Pp, const C2dPoint &p1, const C2dPoint &p2) {
    C3dPoint Q = reconstruct(P, Pp, p1, p2);
    return Q.testInFront(P, Pp);
}

/*bool testPointInFront(const Eigen::Matrix &P, const Eigen::Matrix &Pp, const C3dPoint &Q)
{
    return getDepth(P, Q.asCV())>0 && getDepth(Pp, Q.asCV())>0;
}*/

//Inefficient! Do all points at once.
/*double getDepth(const Eigen::Matrix &P, const C3dPoint &Q)
{
    return getDepth(P, Q.asCV());
    / * /if(IS_DEBUG) CHECK(!zero(Q(4)-1), "getDepth: Coord not homog.");
//    RowVector m3 = P.submatrix(3,3,1,3); //todo: need the homo bit here..?
    Matrix M = P.submatrix(1,3,1,3);
    RowVector m3 = M.row(3); //(top of p161)
    double sdM=sign(M.determinant());
    C3dPoint Qhomo(4);
    Qhomo(1) = Q(1); Qhomo(2) = Q(2); Qhomo(3) = Q(3); Qhomo(4) = 1;
    C3dPoint p = P*Qhomo;
    return (sdM*p(3)/sqrt(m3.sum_square()));* /
}*/

/*inline ColumnVector locToMat(const CLocation loc)
{
    ColumnVector v(3); v << loc.x() << loc.y() << 1;
    return v;
};*/

/*Matrix getCanonicalP()
{
    return R0 | Origin;
}*/

bool clip(CLocation loc, const int IM_WIDTH, const int IM_HEIGHT) {
    int x = loc.x();
    int y = loc.y();
    return x < 0 || x >= IM_WIDTH || y < 0 || y >= IM_HEIGHT;
}

ostream& operator<<(ostream& s, const C3dPoint& X) {
    s << "(" << X.getX() << ", " << X.getY() << ", " << X.getZ() << ")" << flush;
    return s;
}

ostream& operator<<(ostream& s, const C2dPoint& X) {
    //s << setprecision(23);
    s << "(" << X.getX() << ", " << X.getY() << ")" << flush;
    //s << setprecision(3);
    return s;
}

ostream& operator<<(ostream& s, const CCamCalibMatrix & X) {
    s << "(" << X.adK[0] << ", " << X.adK[1] << ", " << X.adK[2] << ")\n";
    s << "(" << X.adK[3] << ", " << X.adK[4] << ", " << X.adK[5] << ")\n";
    s << "(" << X.adK[6] << ", " << X.adK[7] << ", " << X.adK[8] << ")\n" << endl;
    s << "Inverse: \n";
    s << "(" << X.K_inv.adK_inv[0] << ", " << X.K_inv.adK_inv[1] << ", " << X.K_inv.adK_inv[2] << ")\n";
    s << "(" << X.K_inv.adK_inv[3] << ", " << X.K_inv.adK_inv[4] << ", " << X.K_inv.adK_inv[5] << ")\n";
    s << "(" << X.K_inv.adK_inv[6] << ", " << X.K_inv.adK_inv[7] << ", " << X.K_inv.adK_inv[8] << ")\n" << endl; //Todo: RD coeffs
    return s;
}

ostream& operator<<(ostream& s, const C3dRotationMat& X) {
    s << "(" << X[0] << ", " << X[1] << ", " << X[2] << ")\n";
    s << "(" << X[3] << ", " << X[4] << ", " << X[5] << ")\n";
    s << "(" << X[6] << ", " << X[7] << ", " << X[8] << ")\n" << flush;
    return s;
}

ostream& operator<<(ostream& s, const CCamera& X) {
    s << "(" << X.adCam[0] << ", " << X.adCam[1] << ", " << X.adCam[2] << ", " << X.adCam[3] << ")\n";
    s << "(" << X.adCam[4] << ", " << X.adCam[5] << ", " << X.adCam[6] << ", " << X.adCam[7] << ")\n";
    s << "(" << X.adCam[8] << ", " << X.adCam[9] << "," << X.adCam[10] << ", " << X.adCam[11] << ")\n";
    return s;
}

ostream& operator<<(ostream& s, const C3dRotationQuat& X) {
    s << "(" << X[0] << ", " << X[1] << ", " << X[2] << ", " << X[3] << " (" << X.angle() << "))" << flush;
    return s;
}

ostream& operator<<(ostream& s, const CLocation& X) {
    double x = X.dx(), y = X.dy();

    if (x == 0 && y == 0)
        cout << "Zero being printed";

    s << "(" << x << ", " << y << ")";
    return s;
}

C3dPoint C3dRotationMat::operator*(const C3dPoint & p) const {
    C3dPoint rotated(p);
    rotated.rotate(*this);
    return rotated;
};

C3dPoint crossproduct(const C3dPoint & p1, const C3dPoint & p2) {
    return C3dPoint(p1.getY() * p2.getZ() - p1.getZ() * p2.getY(),
            p1.getZ() * p2.getX() - p1.getX() * p2.getZ(),
            p1.getX() * p2.getY() - p1.getY() * p2.getX());
}

double dotproduct(const C3dPoint & p1, const C3dPoint & p2) {
    return p1.getX() * p2.getX() + p1.getY() * p2.getY() + p1.getZ() * p2.getZ();
}

double dotproduct(const C2dPoint & p1, const C2dPoint & p2) {
    return p1.getX() * p2.getX() + p1.getY() * p2.getY();
}

bool C3dRotationMat::operator==(const C3dRotationMat & p) const {
    double ssd = 0;
    for (int i = 0; i < 9; i++)
        ssd += sqr(p.R[i] - R[i]);
    return zero(ssd);
};

C3dPoint C3dRotationMat::headingVector() const {
    return C3dPoint(R[2], R[5], R[8]);
};

CCamera C3dRotationMat::operator|(const C3dPoint & p) const {
    CCamera P;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            P.adCam[i * 4 + j] = R[3 * i + j];

    P.adCam[0 + 3] = p.getX();
    P.adCam[4 + 3] = p.getY();
    P.adCam[8 + 3] = p.getZ();

    return P;
};

bool C3dPoint::operator==(const C3dPoint & p) const {
    return zero(sqr(p.x - x) + sqr(p.y - y) + sqr(p.z - z));
};

C3dPoint operator*(double p, const C3dPoint & T) {
    return T*p;
};

void printMat2(const C3dPoint &G, const char * caption) {
    cout << caption << G << endl;
}

void printMat2(const C3dRotation &G, const char * caption) {
    cout << caption << endl << G << endl;
}

///////////////// Quaternion: ///////////////////////////
//From http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/ethan.htm
#define MAT_TO_QUAT(Rot) \
cout << "This appears broken\n";\
if( trace > 0.0000001 ) { \
double s = 0.5 / sqrt(trace);\
Q[3] = 0.25 / s;\
Q[0] = ( Rot(3,2) - Rot(2,3) ) * s;\
Q[1] = ( Rot(1,3) - Rot(3,1) ) * s;\
Q[2] = ( Rot(2,1) - Rot(1,2) ) * s;\
} else {\
if ( Rot(1,1) > Rot(2,2) && Rot(1,1) > Rot(3,3) )\
{\
    double s = 2.0 * sqrt( 1.0 + Rot(1,1) - Rot(2,2) - Rot(3,3));\
    Q[3] = (Rot(3,2) - Rot(2,3) ) / s;\
    Q[0] = 0.25 * s;\
    Q[1] = (Rot(1,2) + Rot(2,1) ) / s;\
    Q[2] = (Rot(1,3) + Rot(3,1) ) / s;\
} else if (Rot(2,2) > Rot(3,3)) {\
    double s = 2.0 * sqrt( 1.0 + Rot(2,2) - Rot(1,1) - Rot(3,3));\
    Q[3] = (Rot(1,3) - Rot(3,1) ) / s;\
    Q[0] = (Rot(1,2) + Rot(2,1) ) / s;\
    Q[1] = 0.25 * s;\
    Q[2] = (Rot(2,3) + Rot(3,2) ) / s;\
} else {\
    double s = 2.0 * sqrt( 1.0 + Rot(3,3) - Rot(1,1) - Rot(2,2) );\
    Q[3] = (Rot(2,1) - Rot(1,2) ) / s;\
    Q[0] = (Rot(1,3) + Rot(3,1) ) / s;\
    Q[1] = (Rot(2,3) + Rot(3,2) ) / s;\
    Q[2] = 0.25 * s;\
}\
}

C3dRotationQuat::C3dRotationQuat(const CCamera & cam) {
    // This appears broken    double trace = cam.adCam[0] + cam.adCam[5] + cam.adCam[10] + 1.0;
    //MAT_TO_QUAT(ACCESS_CAM_AS_MAT)

    Eigen::Matrix3d rot;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            rot(r, c) = cam.rowData(r)[c];

    fromMat(rot);

    //normalise();
    lengthOk();
}

C3dRotationQuat::C3dRotationQuat(const Eigen::Matrix3d & rot) {
    fromMat(rot);
}

void C3dRotationQuat::fromMat(const Eigen::Matrix3d & rot) {
    /*    double trace = rot.trace();

        MAT_TO_QUAT(ACCESS_EIGENMAT_AS_MAT)

        normalise();
        lengthOk();

        BROKEN -- try this from boost instead
     */

    double fTrace = rot.trace();
    double fRoot;

    //From http://www.geometrictools.com/LibFoundation/Mathematics/Wm4Quaternion.inl
    double m_afTuple[4];
    if (fTrace > (double) 0.0) //0 is w
    {
        // |w| > 1/2, may as well choose w > 1/2
        fRoot = sqrt(fTrace + (double) 1.0); // 2w
        m_afTuple[0] = ((double) 0.5) * fRoot;
        fRoot = ((double) 0.5) / fRoot; // 1/(4w)
        m_afTuple[1] = (rot(2, 1) - rot(1, 2)) * fRoot;
        m_afTuple[2] = (rot(0, 2) - rot(2, 0)) * fRoot;
        m_afTuple[3] = (rot(1, 0) - rot(0, 1)) * fRoot;
    } else {
        // |w| <= 1/2
        int i = 0;
        if (rot(1, 1) > rot(0, 0)) {
            i = 1;
        }
        if (rot(2, 2) > rot(i, i)) {
            i = 2;
        }
        //        int j = ms_iNext[i];
        //        int k = ms_iNext[j];
        int j = (i + 1);
        j %= 3;
        int k = (j + 1);
        k %= 3;

        fRoot = sqrt(rot(i, i) - rot(j, j) - rot(k, k)+(double) 1.0);
        //double* apfQuat[3] = { &m_afTuple[1], &m_afTuple[2], &m_afTuple[3] };
        m_afTuple[i + 1] = ((double) 0.5) * fRoot;
        fRoot = ((double) 0.5) / fRoot;
        m_afTuple[0] = (rot(k, j) - rot(j, k)) * fRoot;
        m_afTuple[j + 1] = (rot(j, i) + rot(i, j)) * fRoot;
        m_afTuple[k + 1] = (rot(k, i) + rot(i, k)) * fRoot;
    }

    Q[0] = m_afTuple[1];
    Q[1] = m_afTuple[2];
    Q[2] = m_afTuple[3];
    Q[3] = m_afTuple[0];

    lengthOk();
}

void C3dRotationQuat::normalise() {
    double length = ((sqr(Q[0]) + sqr(Q[1]) + sqr(Q[2]) + sqr(Q[3])));
    double check = length - 1;
    if (check > 0.0000001 || check < -0.0000001) {
        double scale = 1.0 / sqrt(length);
        Q[0] *= scale;
        Q[1] *= scale;
        Q[2] *= scale;
        Q[3] *= scale;
    }
}

void C3dRotationQuat::lengthOk() const {
#ifdef _DEBUG
    double err = ((sqr(Q[0]) + sqr(Q[1]) + sqr(Q[2]) + sqr(Q[3]) - 1)*100);
    if (!zero(err)) {
        cout << *this << endl << err << endl;
        if(IS_DEBUG) CHECK(1, "New quaternion not normalised");
    }
#endif
}

void C3dRotationQuat::quatMult(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dRotationQuat & Qres) {
    //v(4)=dotproduct(a,quatConj(b));
    //v.rows(1,3)=crossproduct(a_vec,b_vec)   +   a(4)*b_vec    +   b(4)*a_vec;

    Qres.Q[0] = Q1[1] * Q2[2] - Q1[2] * Q2[1] + Q1[3] * Q2[0] + Q2[3] * Q1[0];
    Qres.Q[1] = Q1[2] * Q2[0] - Q1[0] * Q2[2] + Q1[3] * Q2[1] + Q2[3] * Q1[1];
    Qres.Q[2] = Q1[0] * Q2[1] - Q1[1] * Q2[0] + Q1[3] * Q2[2] + Q2[3] * Q1[2];

    Qres.Q[3] = Q1[3] * Q2[3] - (Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2]);
}

void C3dRotationQuat::quatMultIntoVec(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dPoint & Qres) {
    //v(4)=dotproduct(a,quatConj(b));
    //v.rows(1,3)=crossproduct(a_vec,b_vec)   +   a(4)*b_vec    +   b(4)*a_vec;

    Qres.x = Q1[1] * Q2[2] - Q1[2] * Q2[1] + Q1[3] * Q2[0] + Q2[3] * Q1[0];
    Qres.y = Q1[2] * Q2[0] - Q1[0] * Q2[2] + Q1[3] * Q2[1] + Q2[3] * Q1[1];
    Qres.z = Q1[0] * Q2[1] - Q1[1] * Q2[0] + Q1[3] * Q2[2] + Q2[3] * Q1[2];

    if(IS_DEBUG) CHECK(!zero(Q1[3] * Q2[3] - (Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2])), "Bad rotation");
}

//Quaternion * vector

void C3dRotationQuat::quatMultByVec(const C3dRotationQuat & Q1, const C3dPoint & vec, C3dRotationQuat & Qres) {
    //v(4)=dotproduct(a,quatConj(b));
    //v.rows(1,3)=crossproduct(a_vec,b_vec)   +   a(4)*b_vec    +   b(4)*a_vec;

    Qres.Q[0] = Q1[1] * vec.getZ() - Q1[2] * vec.getY() + Q1[3] * vec.getX();
    Qres.Q[1] = Q1[2] * vec.getX() - Q1[0] * vec.getZ() + Q1[3] * vec.getY();
    Qres.Q[2] = Q1[0] * vec.getY() - Q1[1] * vec.getX() + Q1[3] * vec.getZ();

    Qres.Q[3] = -(Q1[0] * vec.getX() + Q1[1] * vec.getY() + Q1[2] * vec.getZ());
}

//Quaternion conj. * vector

void C3dRotationQuat::conjQuatMultVec(const C3dRotationQuat & Q1, const C3dPoint & vec, C3dRotationQuat & Qres) {
    //v(4)=dotproduct(a,quatConj(b));
    //v.rows(1,3)=crossproduct(a_vec,b_vec)   +   a(4)*b_vec    +   b(4)*a_vec;

    Qres.Q[0] = Q1[2] * vec.getY() - Q1[1] * vec.getZ() + Q1[3] * vec.getX();
    Qres.Q[1] = Q1[0] * vec.getZ() - Q1[2] * vec.getX() + Q1[3] * vec.getY();
    Qres.Q[2] = Q1[1] * vec.getX() - Q1[0] * vec.getY() + Q1[3] * vec.getZ();

    Qres.Q[3] = (Q1[0] * vec.getX() + Q1[1] * vec.getY() + Q1[2] * vec.getZ());
}

//Quaternion * conj. quaternion

void C3dRotationQuat::quatMultConj(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dRotationQuat & Qres) {
    //v(4)=dotproduct(a,quatConj(b));
    //v.rows(1,3)=crossproduct(a_vec,b_vec)   +   a(4)*b_vec    +   b(4)*a_vec;

    Qres.Q[0] = ((Q1[2] * Q2[1] - Q1[1] * Q2[2]) - Q1[3] * Q2[0]) + Q2[3] * Q1[0];
    Qres.Q[1] = ((Q1[0] * Q2[2] - Q1[2] * Q2[0]) - Q1[3] * Q2[1]) + Q2[3] * Q1[1];
    Qres.Q[2] = ((Q1[1] * Q2[0] - Q1[0] * Q2[1]) - Q1[3] * Q2[2]) + Q2[3] * Q1[2];

    Qres.Q[3] = Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2] + Q1[3] * Q2[3]; //just dot prod
}

//Quaternion * conj. quaternion

void C3dRotationQuat::quatMultConjIntoVec(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dPoint & Qres) {
    //v(4)=dotproduct(a,quatConj(b));
    //v.rows(1,3)=crossproduct(a_vec,b_vec)   +   a(4)*b_vec    +   b(4)*a_vec;

    Qres.x = ((Q1[2] * Q2[1] - Q1[1] * Q2[2]) - Q1[3] * Q2[0]) + Q2[3] * Q1[0];
    Qres.y = ((Q1[0] * Q2[2] - Q1[2] * Q2[0]) - Q1[3] * Q2[1]) + Q2[3] * Q1[1];
    Qres.z = ((Q1[1] * Q2[0] - Q1[0] * Q2[1]) - Q1[3] * Q2[2]) + Q2[3] * Q1[2];

    if(IS_DEBUG) CHECK(!zero(Q1[0] * Q2[0] + Q1[1] * Q2[1] + Q1[2] * Q2[2] + Q1[3] * Q2[3]), "Bad rotation (probably scale overflow creating a massive vector)"); //just dot prod
}

C3dRotationQuat::C3dRotationQuat(C3dPoint axis, double angle) {

    Q[3] = cos(0.5 * angle);

    double dLength = axis.length();
    if (Q[3] < 1 && Q[3] > -1) {
        axis /= (dLength);
        axis *= sqrt(1 - sqr(Q[3]));
        //cout << axis.sum_square() + sqr(Q[3]) << endl;
    } else
        axis *= 0;

    if (angle < 0) {
        angle = -angle;
        axis *= -1;
    }
    double dRotations = angle * 0.5 * M_1_PI;
    if ((int) (floor(dRotations)) % 2 == 1) axis *= -1;

    Q[0] = axis.getX();
    Q[1] = axis.getY();
    Q[2] = axis.getZ();

    lengthOk();
}

void C3dRotationQuat::setRandom(double dAngle) {
    C3dPoint axis;
    axis.setRandomNormal();
    *this = C3dRotationQuat(axis, dAngle);
}

void C3dRotationQuat::setRandom() {
    for (int i = 0; i < 4; i++)
        Q[i] = CRandom::Normal();
    normalise();
}

void C3dPoint::setRandom(double dDepth) {
    x = CRandom::Uniform(-0.5, 0.5);
    y = CRandom::Uniform(-0.5, 0.5);
    z = dDepth + CRandom::Uniform(-0.5, 0.5);
}

void C3dPoint::setRandomPlanar(double dDepth) {
    setRandom(0);
    z = dDepth;
}

C3dRotationQuat C3dRotationQuat::operator*(const C3dRotationQuat & R2) const {
    C3dRotationQuat multRot;

    quatMult(*this, R2, multRot);

    return multRot;
}

C3dPoint C3dRotationQuat::operator*(const C3dPoint & p) const {
    C3dPoint multPoint;

    //v=q*v*q_conjugate

    C3dRotationQuat temp;
    quatMultByVec(*this, p, temp);
    quatMultConjIntoVec(temp, *this, multPoint);

    /*if(!zero(multPoint.sum_square() -p.sum_square()))
    {
        cout << "multPoint=" << multPoint << endl;
        cout << "p=" << p << endl;
        cout << "this=" << *this << endl;
        quatMultByVec(*this, p, temp);
        quatMultConjIntoVec(temp, *this, multPoint);
    }*/

    return multPoint;
}

CCamera C3dRotationQuat::operator|(const C3dPoint & p) const {
    CCamera cam;
    //PRINTMAT(asMat());
    //PRINTMAT(p.asCV());

    MAT_FROM_QUAT(ACCESS_CAM_AS_MAT);

    cam.adCam[0 + 3] = p.getX();
    cam.adCam[4 + 3] = p.getY();
    cam.adCam[8 + 3] = p.getZ();

    //cout << cam << "=CAM\n";
    return cam;
}

bool C3dRotationQuat::operator==(const C3dRotationQuat & p) const {
    //return p.Q[0] == Q[0] && p.Q[1] == Q[1] && p.Q[2] == Q[2] && p.Q[3] == Q[3];
    return zero(sqr(p.Q[0] - Q[0]) + sqr(p.Q[1] - Q[1]) + sqr(p.Q[2] - Q[2]) + sqr(p.Q[3] - Q[3]))
            || zero(sqr(p.Q[0] + Q[0]) + sqr(p.Q[1] + Q[1]) + sqr(p.Q[2] + Q[2]) + sqr(p.Q[3] + Q[3]));
}

void C3dRotationQuat::toPhiThetaPsi2(double & phi, double & theta, double & psi) const {
    const double qx = Q[0];
    const double qy = Q[1];
    const double qz = Q[2];
    const double qw = Q[3];

    double n = 1.0; //norm();
    double s = n > 0 ? 2. / (n * n) : 0.;

    double m00, /*m01, m02,*/ m10, /*m11, m12,*/ m20, m21, m22;

    double xs = qx*s;
    double ys = qy*s;
    double zs = qz*s;

    double wx = qw*xs;
    double wy = qw*ys;
    double wz = qw*zs;

    double xx = qx*xs;
    double xy = qx*ys;
    double xz = qx*zs;

    double yy = qy*ys;
    double yz = qy*zs;

    double zz = qz*zs;

    m00 = 1.0 - (yy + zz);
    //m11 = 1.0 - (xx + zz);
    m22 = 1.0 - (xx + yy);

    m10 = xy + wz;
    //m01 = xy - wz;

    m20 = xz - wy;
    //m02 = xz + wy;
    m21 = yz + wx;
    //m12 = yz - wx;

    phi = atan2(m21, m22);
    theta = atan2(-m20, sqrt(m21 * m21 + m22 * m22));
    psi = atan2(m10, m00);

    double sphi = sin(phi);
    double stheta = sin(theta);
    double spsi = sin(psi);
    double cphi = cos(phi);
    double ctheta = cos(theta);
    double cpsi = cos(psi);

    double _r[3][3] = {//create rotational Matrix
        {cpsi*ctheta, cpsi * stheta * sphi - spsi*cphi, cpsi * stheta * cphi + spsi * sphi},
        {spsi*ctheta, spsi * stheta * sphi + cpsi*cphi, spsi * stheta * cphi - cpsi * sphi},
        { -stheta, ctheta*sphi, ctheta * cphi}
    };

    double _w = sqrt(max<double>(0, 1 + _r[0][0] + _r[1][1] + _r[2][2])) / 2.0;
    double _x = sqrt(max<double>(0, 1 + _r[0][0] - _r[1][1] - _r[2][2])) / 2.0;
    double _y = sqrt(max<double>(0, 1 - _r[0][0] + _r[1][1] - _r[2][2])) / 2.0;
    double _z = sqrt(max<double>(0, 1 - _r[0][0] - _r[1][1] + _r[2][2])) / 2.0;
    C3dRotationQuat qBack(
            (_r[2][1] - _r[1][2]) >= 0 ? fabs(_x) : -fabs(_x),
            (_r[0][2] - _r[2][0]) >= 0 ? fabs(_y) : -fabs(_y),
            (_r[1][0] - _r[0][1]) >= 0 ? fabs(_z) : -fabs(_z),
            _w);

    if (qBack != *this) {
        cout << "Quat: " << *this << endl;
        cout << "Converted back to: " << qBack << endl;
        //THROW( "Quaternion conversion failed")
        cout << "Warning: Quaternion conversion from toro failed\n";
    }
}

void C3dRotationQuat::toPhiThetaPsi(double & phi, double & theta, double & psi) const {
    //From http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/index.htm
    //Todo: check that inverse is the same
    const double qx = Q[2];
    const double qy = Q[1];
    const double qz = Q[0];
    const double qw = Q[3];
    const double eps = 1e-7;

    const double xy_plus_zw = 2 * (qx * qy + qz * qw);
    if (fabs(xy_plus_zw - 1.0) < eps)// (north pole)
    {
        theta = 2 * atan2(qx, qw);
        phi = M_PI / 2;
        psi = 0;
    } else if (fabs(xy_plus_zw + 1.0) < eps)// (south pole)
    {
        theta = -2 * atan2(qx, qw);
        phi = -M_PI / 2;
        psi = 0;
    } else {
        theta = atan2(2 * qy * qw - 2 * qx*qz, 1 - 2 * sqr(qy) - 2 * sqr(qz));
        phi = asin(xy_plus_zw);
        psi = atan2(2 * qx * qw - 2 * qy*qz, 1 - 2 * sqr(qx) - 2 * sqr(qz));
    }


    double sphi = sin(phi);
    double stheta = sin(theta);
    double spsi = sin(psi);
    double cphi = cos(phi);
    double ctheta = cos(theta);
    double cpsi = cos(psi);

    double _r[3][3] = {//create rotational Matrix
        {cpsi*ctheta, cpsi * stheta * sphi - spsi*cphi, cpsi * stheta * cphi + spsi * sphi},
        {spsi*ctheta, spsi * stheta * sphi + cpsi*cphi, spsi * stheta * cphi - cpsi * sphi},
        { -stheta, ctheta*sphi, ctheta * cphi}
    };

    double _w = sqrt(max<double>(0, 1 + _r[0][0] + _r[1][1] + _r[2][2])) / 2.0;
    double _x = sqrt(max<double>(0, 1 + _r[0][0] - _r[1][1] - _r[2][2])) / 2.0;
    double _y = sqrt(max<double>(0, 1 - _r[0][0] + _r[1][1] - _r[2][2])) / 2.0;
    double _z = sqrt(max<double>(0, 1 - _r[0][0] - _r[1][1] + _r[2][2])) / 2.0;
    C3dRotationQuat qBack(
            (_r[2][1] - _r[1][2]) >= 0 ? fabs(_x) : -fabs(_x),
            (_r[0][2] - _r[2][0]) >= 0 ? fabs(_y) : -fabs(_y),
            (_r[1][0] - _r[0][1]) >= 0 ? fabs(_z) : -fabs(_z),
            _w);

    if (qBack != *this) {
        cout << "Quat: " << *this << endl;
        cout << "Converted back to: " << qBack << endl;
        //THROW( "Quaternion conversion failed")
        cout << "Warning: Quaternion conversion failed\n";
    }
}

C3dPoint C3dRotationQuat::headingVector() const {
    C3dPoint p(0, 0, 1);
    //    DEBUGONLY(C3dPoint q = p;)
    p.rotate(*this);
    //DEBUGONLY(double cosAng=dotproduct(p, q);
    //double ang = fabs(acos(cosAng));
    //)
    //if(IS_DEBUG) CHECK(!zero(ang - fabs(angle())), "Quat rotation failed");

    return p;
}

void C3dPoint::rotate(const C3dRotationQuat & R) {
    //v=q*v*q_conjugate

    C3dRotationQuat temp;
    C3dRotationQuat::quatMultByVec(R, *this, temp);
    C3dRotationQuat::quatMultConjIntoVec(temp, R, *this);
}

void C3dPoint::rotateInv(const C3dRotationQuat & R) {
    //v=q_conjugate*v*q

    C3dRotationQuat temp;
    C3dRotationQuat::conjQuatMultVec(R, *this, temp);
    C3dRotationQuat::quatMultIntoVec(temp, R, *this);
}

bool isRotMat(const Eigen::Matrix3d & R) {
    Eigen::Matrix3d RRt = R * R.transpose();

    return zero(RRt.trace() - 3) && zero(R.determinant() - 1);
}

void Xmat(Eigen::Vector3d const & translation, Eigen::Matrix3d & transMat) {
    transMat << 0, -translation(2), translation(1),
            translation(2), 0, -translation(0),
            -translation(1), translation(0), 0;
}

void makePp(CCamera & Pp, const int nPossibility, const C3dRotation & R_poss, const C3dRotation & R_poss2, const C3dPoint & t_norm) {
    switch (nPossibility) {
        case 0:
            Pp = R_poss | t_norm;
            return;
        case 1:
            Pp = R_poss | -t_norm;
            return;
        case 2:
            Pp = R_poss2 | t_norm;
            return;
        case 3:
            Pp = R_poss2 | -t_norm;
            return;
        default:
            THROW("Bad camera possibility id");
    }
}

void getCamsFromE(const Eigen::Matrix3d & E, CCamera aPp[4]) {
    Eigen::Matrix3d D_getRfromE;
    D_getRfromE << 0, -1, 0, 1, 0, 0, 0, 0, 1;

#if EIGEN_VERSION_AT_LEAST(2,90,0)
    Eigen::JacobiSVD<Eigen::Matrix3d > svdE(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
#else
    Eigen::SVD<Eigen::Matrix3d > svdE(E);
#endif

    Eigen::Matrix3d U = svdE.matrixU();
    Eigen::Matrix3d V = svdE.matrixV();

    const Eigen::Matrix3d & UV = U * V;
    if (!isRotMat(UV)) {
        /* cout << (U * V).determinant() << " = det(UV) ";
        cout << "Using -E to force Rot Mat\n"; */

        if (!isRotMat(U))
            U *= -1;
        else
            V *= -1;

        if (!isRotMat(U * V)) {
            cout << (U * V).determinant() << " = det(UV)\n";
            THROW("CAM_SELECT_ERROR: Rot Mat still not a rot mat\n");
        }
    }

    if(IS_DEBUG) CHECK((svdE.singularValues()(0) / svdE.singularValues()(1)) > 1.2, "E is not an essential matrix");
    if(IS_DEBUG) CHECK(!zero(0.01*svdE.singularValues()(2)/svdE.singularValues()(1)), "E is not an essential matrix");

    const Eigen::Matrix3d & V_t = V.transpose();

    C3dPoint t_norm(U(0, 2), U(1, 2), U(2, 2)); // = U.column(3); // not sure if this is correct up to scale

    const Eigen::Matrix3d & R1 = U * D_getRfromE * V_t;
    C3dRotation R_poss(R1);
    const Eigen::Matrix3d & R2 = U * (D_getRfromE.transpose()) * V_t;
    C3dRotation R_poss_t(R2);
    //if(IS_DEBUG) CHECK(!isRotMat(R_poss) || !isRotMat(R_poss_t), "CSLAMLocMatch::getStructure: Calculated rotation aren't rotation matrices");
    //isRotMat
    /*#define IRM(R) if(!isRotMat(R)) { cout << #R " is not a rotation Eigen::Matrix:"; PRINTMAT(R); }
     IRM(R_poss);
     IRM(R_poss_t);
     IRM(U);
     IRM(V);
     SVD(-E, D, U, V);
     IRM(U);
     IRM(V);*/

    for (int i = 0; i < 4; i++) //iterate over 4 possibilities
    {
        makePp(aPp[i], i, R_poss, R_poss_t, t_norm);
    }

    ////////////////////////// Check E is good (can we compute it from 2 cameras??)
    //Actually we're now using getFundam... to get an F, so this won't work
    /*
     Eigen::Matrix Tx = Xmat(t_norm);
     Eigen::Matrix Erecon = Tx*R_poss;
     DiagonalMatrix D_eIdeal(3); D_eIdeal=0; D_eIdeal(1)=1; D_eIdeal(2)=1;
     E=U*D_eIdeal*V.t();
     if(!matrixEqUpToSign(Erecon, E))
     {
     PRINTMAT(E);
     PRINTMAT(Erecon);
     }
     if(IS_DEBUG) CHECK(!matrixEqUpToSign(Erecon, E), "CSLAMLocMatch::getStructure: E from T, R not equal to original E");
     Todo: uncomment this block once estimating E again */
    //E = 0;
    //PRINTMAT(Erecon);
    /*            DiagonalMatrix Dideal(3);
     Dideal=1;
     Dideal(3,3) = 0;
     Eigen::Matrix Erecon2 = U*Dideal*V.t();
     if(IS_DEBUG) CHECK(matfpeq(Erecon2, E), "CSLAMLocMatch::getStructure: E from T, R not equal to original E")*/
    //PRINTMAT(Erecon2);

    ////////////////////////// Check rot angle is good

    /*double dAngle1 = R_poss.angle();
    cout << dAngle1 << "=angle 1, ";
    double dAngle2 = R_poss_t.angle();
    cout << dAngle2 << "=angle 2\n";*/

    /*if(dAngle1 > 3)
     {
     Pp[0] = 0;
     Pp[1] = 0;
     }
     else if(dAngle2 > 3)
     {
     Pp[2] = 0;
     Pp[3] = 0;
     }
     ////////////////////////*/
}

void getCamsFromE(const double * pdE, CCamera * aPp) {
    Eigen::Matrix3d E(pdE);
    E.transposeInPlace();
    getCamsFromE(E, aPp);
}

const int PURE_ROT_1 = 10, PURE_ROT_2 = 20;

void removePointsBehind(CPointVec2d & ap1, CPointVec2d & ap2, const CCamera & Pp, const CCamera * pPp2 = 0) {
    //Now remove all points behind cam
    CPointVec2d::iterator pGoodPoint1 = ap1.begin();
    CPointVec2d::iterator pGoodPoint2 = ap2.begin();
    CPointVec2d::const_iterator pPoint1 = ap1.begin();
    CPointVec2d::const_iterator pPoint2 = ap2.begin();
    int nNewSize = 0;
    static const CCamera Pagain;
    for (; pPoint1 != ap1.end(); pPoint1++, pPoint2++) {
        if (testPair(Pagain, Pp, *pPoint1, *pPoint2)
                || (pPp2 && testPair(Pagain, *pPp2, *pPoint1, *pPoint2))) {
            *pGoodPoint1 = *pPoint1;
            pGoodPoint1++;
            *pGoodPoint2 = *pPoint2;
            pGoodPoint2++;
            nNewSize++;
        }
    }

    ap1.resize(nNewSize);
    ap2.resize(nNewSize);
}

int chooseCamMat(const CCamera * aPp, const int N_CAM_MATS, CPointVec2d & ap1, CPointVec2d & ap2) {
    static const int CAM_NOT_SET = -1;
    static const CCamera P;

    int nNumPoints = ap1.size();

    //CDynArray<int> anChosenCamMat(nNumPoints, CAM_NOT_SET),
    ARRAY(int, anCamMatTimesChosen, N_CAM_MATS);
    for (int nCamMatIdx = 0; nCamMatIdx < N_CAM_MATS; nCamMatIdx++)
        anCamMatTimesChosen[nCamMatIdx] = 0;

    //todo: To speed up: once have found cam, try that 1st for other points.
    int nPoint = 0;
    for (; nPoint < nNumPoints; nPoint++) {
        int nChosenCamMat = CAM_NOT_SET;
        for (int nCamMatIdx = 0; nCamMatIdx < N_CAM_MATS; nCamMatIdx++) {
            if (testPair(P, aPp[nCamMatIdx], ap1[nPoint], ap2[nPoint])) {
                if (nChosenCamMat == CAM_NOT_SET || nChosenCamMat == nCamMatIdx) {
                    nChosenCamMat = nCamMatIdx;
#ifndef _DEBUG
                    break;
#endif
                } else {
                    cout << "One of 2 compatible cameras: ";
                    cout << (aPp[nCamMatIdx]);
                    cout << "Two of 2 compatible cameras: ";
                    cout << (aPp[nChosenCamMat]);
                    cout << nChosenCamMat << " and " << nCamMatIdx << " are both compatible\n";
                    nChosenCamMat = CAM_SELECT_ERROR;
                }
            }
        }

        if (nChosenCamMat == CAM_NOT_SET)
            cout << "WARNING: getRelPos:: Failed to choose a camera\n";
        else if (nChosenCamMat == CAM_SELECT_ERROR)
            cout << "WARNING: getRelPos:: Multiple compatible camera matrices\n";
        else //if (nChosenCamMat != CAM_NOT_SET && nChosenCamMat != CAM_SELECT_ERROR)
        {
            //cout << "Cam " << nChosenCamMat << " chosen\n";
            anCamMatTimesChosen[nChosenCamMat]++;

            const int LEADER_MARGIN = 5 + nPoint / 2;
            if (nPoint > LEADER_MARGIN) {
                int nCamMatIdx = 0;
                for (; nCamMatIdx < N_CAM_MATS; nCamMatIdx++) {
                    if (nChosenCamMat != nCamMatIdx) {
                        if (anCamMatTimesChosen[nChosenCamMat] < LEADER_MARGIN + anCamMatTimesChosen[nCamMatIdx])
                            break;
                    }
                }
                if (nCamMatIdx == N_CAM_MATS)
                    break; //We've a clear winner
            }
        }
    }

    int nBestCamMat = CAM_SELECT_ERROR, nSecondBestCamMat = CAM_SELECT_ERROR;
    int nNumTimesChosen = 0, nNumTimesChosenSecond = 0, nTotalTimesChosen = 0;
    for (int nCam = 0; nCam < N_CAM_MATS; nCam++) {
        nTotalTimesChosen += anCamMatTimesChosen[nCam];
        if (anCamMatTimesChosen[nCam] > nNumTimesChosen) {
            nSecondBestCamMat = nBestCamMat;
            nNumTimesChosenSecond = nNumTimesChosen;

            nNumTimesChosen = anCamMatTimesChosen[nCam];
            nBestCamMat = nCam;
        } else if (anCamMatTimesChosen[nCam] > nNumTimesChosenSecond) {
            nNumTimesChosenSecond = anCamMatTimesChosen[nCam];
            nSecondBestCamMat = nCam;
        }
    }

    //if(nNumTimesChosenSecond * 3 > nNumTimesChosen)
    {
        for (int nCamMatIdx = 0; nCamMatIdx < N_CAM_MATS; nCamMatIdx++)
            if (anCamMatTimesChosen[nCamMatIdx] > 1)
                cout << "Cam " << nCamMatIdx << " chosen " << anCamMatTimesChosen[nCamMatIdx] << "times\n";
    }

    //if(IS_DEBUG) CHECK(nCamMat == CAM_SELECT_ERROR || anCamMatTimesChosen[nCamMat]==0, "getRelPos:: Failed to find camera pair");
    if (nBestCamMat == CAM_SELECT_ERROR || anCamMatTimesChosen[nBestCamMat] == 0) {
        cout << "WARNING: No camera mat found\n";
        return CAM_SELECT_ERROR;
    } else if (nNumTimesChosen < nNumTimesChosenSecond * 2)
        //We've selected 0,1 or 2,3 (same rotation, different T dir.)
    {
        cout << "Possible pure rotation detected. Warning: Scale cannot be tracked through pure rotation\n";
        removePointsBehind(ap1, ap2, aPp[nBestCamMat], aPp + nSecondBestCamMat);

        if (max<int>(nBestCamMat, nSecondBestCamMat) < 2)
            return PURE_ROT_1;
        else if (min<int>(nBestCamMat, nSecondBestCamMat) >= 2)
            return PURE_ROT_2;
        else {
            //No reconstructions available
            ap1.resize(0);
            ap2.resize(0);
            cout << "Cannot choose camera, 2 candidates with different rotations.\n";
            return CAM_SELECT_ERROR;
        }
    }

    if (nPoint == nNumPoints) {
        if (anCamMatTimesChosen[nBestCamMat] < nPoint / 2) cout << "Best cam mat got only " << anCamMatTimesChosen[nBestCamMat] << " of " << nNumPoints << " votes\n";

        if (anCamMatTimesChosen[nBestCamMat] < 2 + nNumPoints / 3) {
            ap1.resize(0);
            ap2.resize(0);
            cout << "This is insufficient\n";
            return CAM_SELECT_ERROR;
        }
    }
    removePointsBehind(ap1, ap2, aPp[nBestCamMat]);

    if(IS_DEBUG) CHECK(ap1.size() < anCamMatTimesChosen[nBestCamMat], "Error selecting points in front");

    return nBestCamMat;
}

bool chooseCamFromE(const double * pdE, CPointVec2d & ap1, CPointVec2d & ap2, CCamera & Pp, const bool bAllowPureRotation) {
    Eigen::Matrix3d E(pdE);
    E.transposeInPlace();
    return chooseCamFromE(E, ap1, ap2, Pp, bAllowPureRotation);
}

//Can now return a pure rotation.

bool chooseCamFromE(const Eigen::Matrix3d & E, CPointVec2d & ap1, CPointVec2d & ap2, CCamera & Pp, const bool bAllowPureRotation) {
    const int N_CAM_MATS = 4;
    CCamera aPp[N_CAM_MATS];
    getCamsFromE(E, aPp);
    int nCamMat = chooseCamMat(aPp, N_CAM_MATS, ap1, ap2);
    if (nCamMat == CAM_SELECT_ERROR) {
        cout << "Error selecting camera\n";
        return false;
    } else if (nCamMat == PURE_ROT_1) {
        Pp = aPp[0];
        if (bAllowPureRotation)
            Pp.makePureRot(); //Otherwise try to let refinement decide... Should really refine both and choose best
    } else if (nCamMat == PURE_ROT_2) {
        Pp = aPp[2];
        if (bAllowPureRotation)
            Pp.makePureRot();
    } else
        Pp = aPp[nCamMat];

    return true;
}

bool chooseCamFromRT(const C3dRotation & R, const C3dPoint & t, CPointVec2d & ap1, CPointVec2d & ap2, CCamera & Pp, const bool bAllowPureRotation) {
    const int N_CAM_MATS = 2;
    CCamera aPp[N_CAM_MATS];
    aPp[0] = R | t;
    aPp[1] = R | -t;

    int nCamMat = chooseCamMat(aPp, N_CAM_MATS, ap1, ap2);
    if (nCamMat == CAM_SELECT_ERROR) {
        cout << "Error selecting camera\n";
        return false;
    } else if (nCamMat == PURE_ROT_1) {
        Pp = aPp[0];
        if (bAllowPureRotation)
            Pp.makePureRot();
        else
            return false;
    } else if (nCamMat == PURE_ROT_2) {
        THROW("Shouldn't be here")
    } else
        Pp = aPp[nCamMat];

    return true;

}

C3dPoint::C3dPoint(const C2dPoint & p, const double depth) : x(p.getX()), y(p.getY()), z(depth) {
};

double camReprojectionErr(const CCamera & P, const CCamera & Pp, const CPointVec2d & p1, const CPointVec2d & p2, CMask * pMask) {
    Eigen::VectorXd vResiduals(p1.size());
    camReprojectionErr(P, Pp, p1, p2, vResiduals, pMask);
    return vResiduals.sum();
}

void camReprojectionErr(const CCamera & P, const CCamera & Pp, const CPointVec2d & p1, const CPointVec2d & p2, Eigen::VectorXd & vResiduals, CMask * pMask) {
    int nBehind = 0;

    const int NUM_POINTS = (int) p1.size();
    const int RESID_SIZE = (int)vResiduals.rows();
    for (int i = 0; i < NUM_POINTS; i++) {
        if (pMask && !(*pMask)[i])
            continue;

        C2dPoint point1 = p1[i];
        C2dPoint point2 = p2[i];
        C3dPoint Q = reconstruct(P, Pp, point1, point2);
        DEBUGONLY(if (!Q.testInFront(P, Pp))
                nBehind++);
        C2dPoint im1PointResid = P * Q - point1;
        C2dPoint im2PointResid = Pp * Q - point2;
        //cout << im1PointResid << endl << im2PointResid << endl;
        if (RESID_SIZE == NUM_POINTS)
            vResiduals(i) = im1PointResid.sum_square() + im2PointResid.sum_square();
        else if (RESID_SIZE == 1)
            vResiduals(0) += im1PointResid.sum_square() + im2PointResid.sum_square();
        else {
            vResiduals(i * 4 + 0) = im1PointResid.getX();
            vResiduals(i * 4 + 1) = im1PointResid.getY();
            vResiduals(i * 4 + 2) = im2PointResid.getX();
            vResiduals(i * 4 + 3) = im2PointResid.getY();
        }
    }

    if (nBehind > 0)
        cout << "Warning: " << nBehind << " reconstructed 3d points no longer in front (refinement still usually improves accuracy)\n";
}

void testK(const CCamCalibMatrix & K, CLocation l, const bool CORRECT_RD) {
    CCamera P;
    P.calibrate(K);

    C2dPoint p(l);
    C2dPoint p_original(p);
    p.calibrate(K, false); //This can correct for RD
    C3dPoint p3_calib(p, 1);
    C2dPoint p_backInIm = p3_calib.photo(P);
    if(IS_DEBUG) CHECK(!zero((p_backInIm - p_original).sum_square()), "Photo and reproj back into image failed");
    if (K.canCorrectRD() && CORRECT_RD) {
        p.correctRD(K);
        C3dPoint p3_calib_correctRD(p, 1);
        C2dPoint p_correctRD_backInIm = p3_calib_correctRD.photo(P);
        double dDistMoved = (p_correctRD_backInIm - p_original).length();
        cout << "RD correction moved " << l << " " << dDistMoved << " pixels to " << p_correctRD_backInIm << endl;
        if(IS_DEBUG) CHECK(dDistMoved > 60, "Radial distortion coeff. moved point more than 60 pixels, probably wrong!!");
    } else
        cout << "No RD correction\n";
}

void calibrateBoWCorr(const CCamCalibMatrix & K, const CBoWCorrespondences * pCorr, CPointVec2d & aCalibratedPoints1, CPointVec2d & aCalibratedPoints2, CInlierProbs & adArrLikelihood, CPointIdentifiers & pointIds, const bool CORRECT_RD) {
    //cout << "Deprecated: Use CBoWCorrespondences::calibrate\n";

    static bool bFirstTime = true;
    if (bFirstTime) {
        bFirstTime = false;
        testK(K, CLocation(1, 1), CORRECT_RD);
        testK(K, CLocation(100, 0), CORRECT_RD);
        testK(K, CLocation(100, 100), CORRECT_RD);
        testK(K, CLocation(0, 1000), false);
    }

    pCorr->calibrate(K, aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, CORRECT_RD);
}

//Estimate SD of camera motion from

CRelPoseSD::CRelPoseSD(const TPointVec2d & m0, const TPointVec2d & m1, const CCamCalibMatrix & K, const double dFeatureLocalisationAccuracy, double dPointDepthMean, double dPointDepthSD) : dBaselineFromPointDepthMean(dPointDepthMean), dBaselineFromPointDepthVar(dPointDepthSD) {
    const int NUM_POINTS = m0.size();
    if (NUM_POINTS <= 5) {
        makeUninformative();
        return;
    }

    double dAvMotion = 0;
    for (int i = 0; i < NUM_POINTS; i++) {
        dAvMotion += sqrt((m1[i] - m0[i]).sum_square());
    }
    dAvMotion /= NUM_POINTS;
    double dAvMotionPx = dAvMotion * K.focalLength();

    const int NUM_DOF = 5; // Min num of points to fit model
    const double dExpectedMotionError = M_SQRT2 * dFeatureLocalisationAccuracy; //Because 2 measurements each have LE sqrt(2)
    double dPropToErr = dExpectedMotionError / (dAvMotionPx * sqrt((double) (NUM_POINTS - NUM_DOF)));

    dRelOrientationSD = 9.87 * dPropToErr; //9.87 comes from fitRotationParam.xls
    dCameraMotionAngleSD = 23.47 * dPropToErr; //23.47 from fitTranslationParam.xls
}

void C2dPoint::addNoise(double dSD) {
    x += CRandom::Normal(0, dSD);
    y += CRandom::Normal(0, dSD);
}

void C3dPoint::addNoise(double dSD) {
    x += CRandom::Normal(0, dSD);
    y += CRandom::Normal(0, dSD);
    z += CRandom::Normal(0, dSD);
}

void C3dPoint::setRandomNormal() {
    x = CRandom::Normal();
    y = CRandom::Normal();
    z = CRandom::Normal();
}

std::ostream& operator<<(std::ostream& s, const C3dPose & X) {
    s << X.t << " " << X.R;
    return s;
}
