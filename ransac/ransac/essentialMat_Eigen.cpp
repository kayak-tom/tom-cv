/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */
/* Based on MATLAB code from http://vis.uky.edu/~stewe/FIVEPOINT
 * Uses fivepointSetupMatricesHelper function so is for academic use only.

Half of the cost is the eigenvector computation, and most of that is the cost of solving
a 10th order polynomial. Increasing the epsilon convergence parameter for the Eigen
eigenvector computation speeds things up. Also add EIGEN_ALWAYS_INLINE_ATTRIB to operator()'s in src/Core/Coeffs.h

The algorithm described in "An efficient solution to the five-point
relative pose problem", D. Nister, 2004 may be faster; it avoids the eigenvector
computation but adds another 10th order poly. to solve instead.

 */
#define OVERRIDE_EPS
extern double SCALE_EPS;

#include "util/exception.h"
#include "polySolve.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

#ifdef EIGEN_POLY //is relatively slow anyway
#include <unsupported/Eigen/Polynomials>
#endif

#include "calibrated_fivepoint_helper.h"

#include "essentialMat_Eigen.h"
#include "util/convert.h"

#include <iostream>
#include "util/opencv.h"
#include <boost/scoped_array.hpp>
#include "makeBasis_GramSchmidt.h"

using namespace std;

//USING_PART_OF_NAMESPACE_EIGEN;
using namespace Eigen;

typedef Matrix<double, 10, 10 > Matrix10d;

ostream & operator<<(ostream& s, const CSimple2dPoint& X) {
    s << "(" << X.getX() << ", " << X.getY() << ")" << flush;
    return s;
}

double getResids(const Matrix3d & E, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1)
{
    double dErrSq = 0;
    for (int i = 0; i < 5; i++) {
        Eigen::Vector3d x(m0[anHypSet[i]].getX(), m0[anHypSet[i]].getY(), 1);
        Eigen::RowVector3d xp(m1[anHypSet[i]].getX(), m1[anHypSet[i]].getY(), 1);
        dErrSq += xp * E * x;
    }
    return dErrSq;
}

class CSortReal
{
public:
    bool operator()(const std::complex<double> & r1, const std::complex<double> & r2)
    {
        return r1.real() < r2.real();
    }
};

//FASTER MATRIC INVERSION for sizes e.g. 10x10: http://en.wikipedia.org/wiki/Matrix_inverse#Blockwise_inversion
template<int N>
class CInv;

//Specialisations for when will not be called

template<> class CInv < 4 > {
public:

    static void invert(const Eigen::Matrix<double, 4, 4 > & M, Eigen::Matrix<double, 4, 4 > & M_inv) {
        M_inv = M.inverse();
    };
};

template<> class CInv < 3 > {
public:

    static void invert(const Eigen::Matrix<double, 3, 3 > & M, Eigen::Matrix<double, 3, 3 > & M_inv) {
        M_inv = M.inverse();
    };
};

template<> class CInv < 2 > {
public:

    static void invert(const Eigen::Matrix<double, 2, 2 > & M, Eigen::Matrix<double, 2, 2 > & M_inv) {
        M_inv = M.inverse();
    };
};

template<int N> class CInv {

    template<const int SPLIT>
    class CInvInt {
    public:

        static void invertSplit(const Eigen::Matrix<double, N, N> & M, Eigen::Matrix<double, N, N> & M_inv) {
            const Matrix<double, SPLIT, SPLIT> & A = M.template block<SPLIT, SPLIT > (0, 0);
            const Matrix<double, SPLIT, N - SPLIT> & B = M.template block<SPLIT, N - SPLIT > (0, SPLIT);
            const Matrix<double, N - SPLIT, SPLIT> & C = M.template block < N - SPLIT, SPLIT > (SPLIT, 0);
            const Matrix<double, N - SPLIT, N - SPLIT> & D = M.template block < N - SPLIT, N - SPLIT > (SPLIT, SPLIT);

            const Matrix<double, SPLIT, SPLIT> & A_inv = A.inverse();

            const Matrix<double, N - SPLIT, SPLIT> & C_A_inv = -C*A_inv; //minus simplifies later
            const Matrix<double, SPLIT, N - SPLIT> & A_inv_B = -A_inv*B; //minus simplifies later
            const Matrix<double, N - SPLIT, N - SPLIT> & DCAB_int = D + C_A_inv * B; //+ cancels - above

            //This is SLOWER slightly
            //Matrix<double, N-SPLIT, SPLIT> C_A_inv = (C*A_inv).lazy(); C_A_inv *= -1; //minus simplifies later
            //Matrix<double, SPLIT, N-SPLIT> A_inv_B = (A_inv*B).lazy(); A_inv_B *= -1; //minus simplifies later
            //Matrix<double, N-SPLIT, N-SPLIT> DCAB_int = (C_A_inv * B).lazy(); DCAB_int += D;  //+ cancels - above

            Matrix<double, N - SPLIT, N - SPLIT> DCAB_inv;
            CInv < N - SPLIT>::invert(DCAB_int, DCAB_inv);

            M_inv.template block<SPLIT, SPLIT > (0, 0) = A_inv + A_inv_B * DCAB_inv * C_A_inv;
            M_inv.template block<SPLIT, N - SPLIT > (0, SPLIT) = A_inv_B * DCAB_inv; //minus added earlier
            M_inv.template block < N - SPLIT, SPLIT > (SPLIT, 0) = DCAB_inv * C_A_inv; //minus added earlier
            M_inv.template block < N - SPLIT, N - SPLIT > (SPLIT, SPLIT) = DCAB_inv;
        }
    };
public:

    static void invert(const Eigen::Matrix<double, N, N> & M, Eigen::Matrix<double, N, N> & M_inv) {
        if (N == 9 || N <= 6)
            //Best to split into 3,3 or 3,3,3 or 3,2. Also--don't want to split 5 into 4.
            CInvInt < 3 > ::invertSplit(M, M_inv);
        else
            CInvInt < 4 > ::invertSplit(M, M_inv);

        //std::cout << M_inv * M << endl;
    }
};

bool checkUpright(const Eigen::Vector3d & v, const double dUprightThresh) {
    return (fabs(v(1) - 1) < dUprightThresh);
}

bool isRotMat2(const Eigen::Matrix3d & R) {
    return zero(R.determinant() - 1);
}

bool upright(const Eigen::Matrix3d & E, const double dUprightThresh) //Todo: code duplication in geom.cpp...
{
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
    if (!isRotMat2(UV)) {
        /*DEBUGONLY(
        cout << (U * V).determinant() << " = det(UV) ";
        cout << "upright: Using -E to force Rot Mat\n");*/

        if (!isRotMat2(U))
            U *= -1;
        else
            V *= -1;
    }

    //if(IS_DEBUG) CHECK ((svdE.singularValues()(0) / svdE.singularValues()(1)) > 1.2, "E is not an essential matrix");
    if (!zero(svdE.singularValues()(2)) || (svdE.singularValues()(0) / svdE.singularValues()(1)) > 1.2) {
        cout << "E is not an essential matrix\n";
        return false;
    }

    const Eigen::Matrix3d & V_t = V.transpose();

    //C3dPoint t_norm(U(0, 2), U(1, 2), U(2, 2)); // = U.column(3); // not sure if this is correct up to scale

    if (fabs(U(1, 2)) > dUprightThresh)
        return false;

    const Eigen::Matrix3d & R1 = U * D_getRfromE * V_t;
    const Eigen::Matrix3d & R2 = U * (D_getRfromE.transpose()) * V_t;

    Vector3d vUpright;
    vUpright << 0, 1, 0;
    //cout << R1*vUpright << endl;
    //cout << R2*vUpright << endl;
    return checkUpright(R1*vUpright, dUprightThresh) || checkUpright(R2*vUpright, dUprightThresh);
}

//Implement http://en.wikipedia.org/wiki/Convolution#Discrete_convolution

template <int N, int M>
class CConv {
public:
    static const int LENGTH = N + M - 1;

    Eigen::Matrix<double, 1, LENGTH> operator()(const Eigen::Matrix<double, 1, N> & p1, const Eigen::Matrix<double, 1, M> & p2) {
        Eigen::Matrix<double, 1, LENGTH> res;
        res.setZero();
        for (int n = 0; n < LENGTH; n++)
            for (int m = std::max<int>(n - M + 1, 0); m <= std::min<int>(n, N - 1); m++) {
                //cout << n << ':' << m << '+' << n-m << endl;
                res(n) += p1[m] * p2[n - m];
            }

        return res;
    }
};

void makeBasisForE(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, Matrix<double, 9, 4 > & EE) {
    const int nPoints = 5;

    Matrix<double, 5, 1> Q1_col1;
    Matrix<double, 5, 1> Q2_col1;
    Matrix<double, 5, 1> Q1_col2;
    Matrix<double, 5, 1> Q2_col2;

    Matrix< double, nPoints, 9 > Q;

    for (int i = 0; i < nPoints; i++) {
        Q1_col1(i) = m0[anHypSet[i]].getX();
        Q1_col2(i) = m0[anHypSet[i]].getY();
        Q2_col1(i) = m1[anHypSet[i]].getX();
        Q2_col2(i) = m1[anHypSet[i]].getY();

        //And put 1s in the last col
        Q(i, 8) = 1;
    }

    Q.col(0) = Q1_col1.array() * Q2_col1.array();
    Q.col(1) = Q1_col2.array() * Q2_col1.array();
    Q.col(2) = Q2_col1;
    Q.col(3) = Q1_col1.array() * Q2_col2.array();
    Q.col(4) = Q1_col2.array() * Q2_col2.array();
    Q.col(5) = Q2_col2;
    Q.col(6) = Q1_col1;
    Q.col(7) = Q1_col2;
    //already Q.col(8) = 1;

    Eigen::JacobiSVD< Matrix< double, nPoints, 9 > > svdQ(Q, Eigen::ComputeFullV);

    /*OpenCV SVD: *
    CvMat Q_t_cv = cvMat(9, nPoints, CV_64F, Q.data());
    double Q_cv_data[9 * nPoints];
    CvMat Q_cv = cvMat(nPoints, 9, CV_64F, Q_cv_data);
    cvTranspose(&Q_t_cv, &Q_cv);

    if(IS_DEBUG) CHECK(cvmGet(&Q_cv, 2,3) != Q(2,3), "Matrix wrap failed");

    CvMat V_t_cv = cvMat(9, 9, CV_64F, V.data());
    double Diag_data[9];
    CvMat Diag = cvMat(5, 1, CV_64F, Diag_data);
    cvSVD( &Q_cv, &Diag, NULL, &V_t_cv, CV_SVD_MODIFY_A | CV_SVD_U_T | CV_SVD_V_T);//
     */

    EE = svdQ.matrixV().block < 9, 4 > (0, 5); //NOT & because access EE.data() as an array
}

/*template<typename N> 
inline void makeOrth(N & m1, const N & m2)*/

using namespace cv;

int
cvSolveCubic_noMalloc( const double* c, double* r )
{
    int n = 0;

    double a0 = 1., a1, a2, a3;
    double x0 = 0., x1 = 0., x2 = 0.;
    //int coeff_count;

    /*if( !CV_IS_MAT(coeffs) )
        CV_Error( !coeffs ? CV_StsNullPtr : CV_StsBadArg, "Input parameter is not a valid matrix" );

    if( !CV_IS_MAT(roots) )
        CV_Error( !roots ? CV_StsNullPtr : CV_StsBadArg, "Output parameter is not a valid matrix" );

    if( (CV_MAT_TYPE(coeffs->type) != CV_32FC1 && CV_MAT_TYPE(coeffs->type) != CV_64FC1) ||
        (CV_MAT_TYPE(roots->type) != CV_32FC1 && CV_MAT_TYPE(roots->type) != CV_64FC1) )
        CV_Error( CV_StsUnsupportedFormat,
        "Both matrices should be floating-point (single or double precision)" );

    coeff_count = coeffs->rows + coeffs->cols - 1;

    if( (coeffs->rows != 1 && coeffs->cols != 1) || (coeff_count != 3 && coeff_count != 4) )
        CV_Error( CV_StsBadSize,
        "The matrix of coefficients must be 1-dimensional vector of 3 or 4 elements" );

    if( (roots->rows != 1 && roots->cols != 1) ||
        roots->rows + roots->cols - 1 != 3 )
        CV_Error( CV_StsBadSize,
        "The matrix of roots must be 1-dimensional vector of 3 elements" );

    if( CV_MAT_TYPE(coeffs->type) == CV_32FC1 )
    {
        const float* c = coeffs->data.fl;
        if( coeffs->rows > 1 )
            step = coeffs->step/sizeof(c[0]);
        if( coeff_count == 4 )
            a0 = c[0], c += step;
        a1 = c[0];
        a2 = c[step];
        a3 = c[step*2];
    }
    else
    {*/
        // = coeffs->data.db;
        a0 = c[0];
        a1 = c[1];
        a2 = c[2];
        a3 = c[3];
//    }

    if( a0 == 0 )
    {
        if( a1 == 0 )
        {
            if( a2 == 0 )
                n = a3 == 0 ? -1 : 0;
            else
            {
                // linear equation
                x0 = -a3/a2;
                n = 1;
            }
        }
        else
        {
            // quadratic equation
            double d = a2*a2 - 4*a1*a3;
            if( d >= 0 )
            {
                d = sqrt(d);
                double q1 = (-a2 + d) * 0.5;
                double q2 = (a2 + d) * -0.5;
                if( fabs(q1) > fabs(q2) )
                {
                    x0 = q1 / a1;
                    x1 = a3 / q1;
                }
                else
                {
                    x0 = q2 / a1;
                    x1 = a3 / q2;
                }
                n = d > 0 ? 2 : 1;
            }
        }
    }
    else
    {
        a0 = 1./a0;
        a1 *= a0;
        a2 *= a0;
        a3 *= a0;

        double Q = (a1 * a1 - 3 * a2) * (1./9);
        double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) * (1./54);
        double Qcubed = Q * Q * Q;
        double d = Qcubed - R * R;

        if( d >= 0 )
        {
            double theta = acos(R / sqrt(Qcubed));
            double sqrtQ = sqrt(Q);
            double t0 = -2 * sqrtQ;
            double t1 = theta * (1./3);
            double t2 = a1 * (1./3);
            x0 = t0 * cos(t1) - t2;
            x1 = t0 * cos(t1 + (2.*CV_PI/3)) - t2;
            x2 = t0 * cos(t1 + (4.*CV_PI/3)) - t2;
            n = 3;
        }
        else
        {
            double e;
            d = sqrt(-d);
            e = pow(d + fabs(R), 0.333333333333);
            if( R > 0 )
                e = -e;
            x0 = (e + Q / e) - a1 * (1./3);
            n = 1;
        }
    }


    r[0] = x0;
    r[1] = x1;
    r[2] = x2;

    return n;
}


int calcEssentialMat_7point(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & models, int & FAIL_COUNT) 
{
    if(IS_DEBUG) CHECK(m0.size() != m1.size() || models.numModels() > 0, "calcEssentialMat_5point_Eigen: Bad parameters");
    if(IS_DEBUG) CHECK((int) anHypSet.size() != 7, "calcEssentialMat_7point: Insufficient points");
    
    Matrix<double, 9, 2 > EE;
    makeBasisForE_GramSchmidt_Vectorise<7>(anHypSet, m1, m0, EE);
    
    //For OpenCV:
        // f1, f2 is a basis => lambda*f1 + mu*f2 is an arbitrary f. matrix.
    // as it is determined up to a scale, normalize lambda & mu (lambda + mu = 1),
    // so f ~ lambda*f1 + (1 - lambda)*f2.
    // use the additional constraint det(f) = det(lambda*f1 + (1-lambda)*f2) to find lambda.
    // it will be a cubic equation.
    // find c - polynomial coefficients.
    Matrix<double, 9, 1 > f1 = EE.col(0);
    Matrix<double, 9, 1 > f2 = EE.col(1);
    Vector4d c;
    Vector3d roots;
    
    f1 -= f2;
    
    double t0 = f2(4)*f2(8) - f2(5)*f2(7);
    double t1 = f2(3)*f2(8) - f2(5)*f2(6);
    double t2 = f2(3)*f2(7) - f2(4)*f2(6);

    c(3) = f2(0)*t0 - f2(1)*t1 + f2(2)*t2;

    c(2) = f1(0)*t0 - f1(1)*t1 + f1(2)*t2 -
           f1(3)*(f2(1)*f2(8) - f2(2)*f2(7)) +
           f1(4)*(f2(0)*f2(8) - f2(2)*f2(6)) -
           f1(5)*(f2(0)*f2(7) - f2(1)*f2(6)) +
           f1(6)*(f2(1)*f2(5) - f2(2)*f2(4)) -
           f1(7)*(f2(0)*f2(5) - f2(2)*f2(3)) +
           f1(8)*(f2(0)*f2(4) - f2(1)*f2(3));

    t0 = f1(4)*f1(8) - f1(5)*f1(7);
    t1 = f1(3)*f1(8) - f1(5)*f1(6);
    t2 = f1(3)*f1(7) - f1(4)*f1(6);

    c(1) = f2(0)*t0 - f2(1)*t1 + f2(2)*t2 -
           f2(3)*(f1(1)*f1(8) - f1(2)*f1(7)) +
           f2(4)*(f1(0)*f1(8) - f1(2)*f1(6)) -
           f2(5)*(f1(0)*f1(7) - f1(1)*f1(6)) +
           f2(6)*(f1(1)*f1(5) - f1(2)*f1(4)) -
           f2(7)*(f1(0)*f1(5) - f1(2)*f1(3)) +
           f2(8)*(f1(0)*f1(4) - f1(1)*f1(3));

    c(0) = f1(0)*t0 - f1(1)*t1 + f1(2)*t2;

    // solve the cubic equation; there can be 1 to 3 roots ...
    //int n = cvSolveCubic( c.data(), roots.data() );
    
    //cv::Mat cv_c(4, 1, CV_64FC1, c.data());
    //cv::Mat cv_roots(4, 1, CV_64FC1, roots.data());
    //int n = cv::solveCubic(cv_c, cv_roots);
    int n = cvSolveCubic_noMalloc(c.data(), roots.data());

    if( n > 3 )
        return 0;

    for( int k = 0; k < n; k++ )
    {
        Matrix3d E;
        // for each root form the fundamental matrix
        double lambda = roots(k), mu = 1.;
        //double s = f1(8)*roots(k) + f2(8);
        //Want lambda^2+mu^2 = 1
        double dScale = 1.0/sqrt(sqr(lambda) + 1);
        lambda *= dScale; mu *= dScale;

        // normalize each matrix, so that F(3,3) (~fmatrix[8]) == 1
        /*if( fabs(s) > DBL_EPSILON )
        {
            mu = 1./s;
            lambda *= mu;
            E(2,2) = 1.;
        }
        else
            E(2,2) = 0.;*/

        for(int i = 0; i < 9; i++ )
            E.data()[i] = f1(i)*lambda + f2(i)*mu;
        
        Matrix3d EET = E * E.transpose();
        const double dTrace = EET.trace();
        EET.diagonal().array() -= 0.5*dTrace;
        EET *= E;
        if(EET.squaredNorm() < 0.001)
        {
            //const double dScale = sqrt(2.0/dTrace);
            makeClosestE(E);
            C3x3MatModel & m = static_cast<C3x3MatModel &>(models.addModel());
            for(int r=0;r<3;r++)
                for(int c=0;c<3;c++)
                    m(r,c) = E(r,c);
        }   
        else FAIL_COUNT++;
        
        //Eigen::JacobiSVD<Matrix3d> svd(E);
        //cout << " " << svd.singularValues().transpose() << endl;

    }

    return models.numModels();
}


/* finds complex roots of a polynomial using Durand-Kerner method:
   http://en.wikipedia.org/wiki/Durand%E2%80%93Kerner_method */
#ifndef USE_OLD_OPENCV
double solvePoly2( const Mat& coeffs0, Mat& roots0, int maxIters )
{
    typedef Complex<double> C;

    double maxDiff = 0;
    int iter, i, j, n;

    CV_Assert( coeffs0.dims <= 2 &&
               (coeffs0.cols == 1 || coeffs0.rows == 1) &&
               (coeffs0.depth() == CV_32F || coeffs0.depth() == CV_64F) &&
               coeffs0.channels() <= 2 );
    n = coeffs0.cols + coeffs0.rows - 2;

    if( ((roots0.rows != 1 || roots0.cols != n) &&
        (roots0.rows != n || roots0.cols != 1)) ||
        (roots0.type() != CV_32FC2 && roots0.type() != CV_64FC2) )
        roots0.create( n, 1, CV_64FC2 );

    AutoBuffer<C> buf(n*2+2);
    C *coeffs = buf, *roots = coeffs + n + 1;
    Mat coeffs1(coeffs0.size(), CV_MAKETYPE(CV_64F, coeffs0.channels()), coeffs0.channels() == 2 ? coeffs : roots);
    coeffs0.convertTo(coeffs1, coeffs1.type());
    if( coeffs0.channels() == 1 )
    {
        const double* rcoeffs = (const double*)roots;
        for( i = 0; i <= n; i++ )
            coeffs[i] = C(rcoeffs[i], 0);
    }

    C p(1, 0), r(0.02, 0.04);

    for( i = 0; i < n; i++ )
    {
        roots[i] = p;
        p = p * r;
    }

    maxIters = maxIters <= 0 ? 1000 : maxIters;
    for( iter = 0; iter < maxIters; iter++ )
    {
        maxDiff = 0;
        for( i = 0; i < n; i++ )
        {
            p = roots[i];
            C num = coeffs[n], denom = 1;
            for( j = 0; j < n; j++ )
            {
                num = num*p + coeffs[n-j-1];
                if( j != i ) denom = denom * (p - roots[j]);
            }
            num /= denom;
            roots[i] = p - num;
            maxDiff = max(maxDiff, abs(num));
        }
        if( maxDiff <= 0 )
            break;
    }

    if( coeffs0.channels() == 1 )
    {
        const double verySmallEps = 1e-100;
        for( i = 0; i < n; i++ )
            if( fabs(roots[i].im) < verySmallEps )
                roots[i].im = 0;
    }

    Mat(roots0.size(), CV_64FC2, roots).convertTo(roots0, roots0.type());
    return maxDiff;
}
#endif

template<class T>
void sortComplex(T & v)
{
    std::vector<std::complex<double> > vRoots(v.rows());

    for (int i = 0; i < v.rows(); i++)
        vRoots[i] = v(i);

    std::sort(vRoots.begin(), vRoots.end(), CSortReal());

    for (int i = 0; i < v.rows(); i++)
        v(i) = vRoots[i];
}

int calcEssentialMat_5point_Eigen(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & pdEssentialMat, const double dUprightThresh, int & FAIL_COUNT, const bool bUseNister04) 
{
    if(IS_DEBUG) CHECK(m0.size() != m1.size() || pdEssentialMat.numModels() > 0, "calcEssentialMat_5point_Eigen: Bad parameters");
    if(IS_DEBUG) CHECK((int) m0.size() < 5, "calcEssentialMat_5point: Insufficient points");
    int nMats = 0;

    Matrix<double, 9, 4 > EE;
    makeBasisForE_GramSchmidt_Vectorise<5>(anHypSet, m0, m1, EE); //Makes NO DIFFERENCE to success rate
    //makeBasisForE(anHypSet, m0, m1, EE);

    Matrix<double, 10, 20, 0 > A;

    /* Test with E from MATLAB*/
    /*EE(0,0) =     -0.3772 ;EE(0,1) =  -0.4926 ; EE(0,2) =   0.2394 ;    EE(0,3) =  -0.4407;
      EE(1,0) =   0.5807   ; EE(1,1) =  -0.3477  ; EE(1,2) =  -0.5481   ; EE(1,3) =   0.0750;
      EE(2,0) =  -0.0893 ;   EE(2,1) =   0.5395  ; EE(2,2) = -0.0460 ;    EE(2,3) =   0.0895;
      EE(3,0) =-0.5015  ;    EE(3,1) =   0.1742  ; EE(3,2) =  -0.6112  ;  EE(3,3) =  -0.1781 ;
      EE(4,0) = -0.2724   ;  EE(4,1) = -0.2674   ; EE(4,2) =  0.1841  ;   EE(4,3) =   0.8230 ;
      EE(5,0) =   0.3913  ;  EE(5,1) =   0.0017  ; EE(5,2) =   0.3133 ;   EE(5,3) =   -0.2162 ;
      EE(6,0) =  0.1600  ;   EE(6,1) =   0.2617 ;  EE(6,2) =   -0.0235  ; EE(6,3) =   0.1592;
      EE(7,0) =   0.0726   ; EE(7,1) =   0.2634  ; EE(7,2) =    0.3618  ; EE(7,3) =   -0.1045;
      EE(8,0) =   -0.0530  ; EE(8,1) =    -0.3251; EE(8,2) =   -0.0560  ; EE(8,3) =     0.0061   ;*/

    fivepointSetupMatricesHelper(EE, A);
    DEBUGONLY(double dSum = A.sum();
    if (isnan(dSum) || isinf(dSum) || dSum > HUGE || dSum < -HUGE) {
        cout << "EIGEN: Warning: A not setup fully (Sum=" << dSum << ")\n"; //This happens with synthetic data because the 2 images are actually identical.
    }
    dSum = EE.sum();
    if (isnan(dSum) || isinf(dSum) || dSum > HUGE || dSum < -HUGE) {
        cout << "EIGEN: Warning: EE not setup fully (Sum=" << dSum << ")\n"; //This happens with synthetic data because the 2 images are actually identical.
    })
    //cout << 'A' << A << endl;


    const int anPermuteA[20] = {0, 3, 1, 2, 4, 10, 6, 12, 5, 11, 7, 13, 16, 8, 14, 17, 9, 15, 18, 19};
    Matrix<double, 10, 20 > Aperm;
    for (int j = 0; j < 20; j++) {
        for (int i = 0; i < 10; i++) {
            Aperm(i, j) = A(i, anPermuteA[j]);
        }
    }

    //cout << "Aperm: " << Aperm;

    const Matrix<double, 10, 10 > & A_lhs = Aperm.block < 10, 10 > (0, 0);
    //Matrix<double, 10, 10 > A_rhs = A_lhs.inverse() * Aperm.block < 10, 10 > (0, 10);
    Matrix<double, 10, 10 > A_rhs = A_lhs.partialPivLu().solve(Aperm.block < 10, 10 > (0, 10));
    
    //Check its working: std::cout << A_rhs - A_lhs.inverse() * Aperm.block < 10, 10 > (0, 10) << endl;
    
    
    Matrix<double, 4, 10 > SOLS;
    Eigen::EigenSolver<Matrix10d>::EigenvalueType eigvals;

    if (bUseNister04) {

        /*/LaRREF(Aperm);
                //cout << "Aperm RREF: " << Aperm;

                LaVectorLongInt pivots( 10);
                LUFactorizeIP(A_lhs, pivots);
                LaLUInverseIP(A_lhs, pivots);
         //cout << A_lhs * A_rhs;

                Aperm(LaIndex(0,9),LaIndex(10,19)).inject(A_lhs * A_rhs);
         //cout << "Aperm RREF: " << Aperm;*/

        Matrix<double, 3, 13 > B;
        B(0, 0) = -A_rhs(5, 0);
        B(0, 1) = A_rhs(4, 0) - A_rhs(5, 1);
        B(0, 2) = A_rhs(4, 1) - A_rhs(5, 2);
        B(0, 3) = A_rhs(4, 2);
        B(0, 4) = -A_rhs(5, 3);
        B(0, 5) = A_rhs(4, 3) - A_rhs(5, 4);
        B(0, 6) = A_rhs(4, 4) - A_rhs(5, 5);
        B(0, 7) = A_rhs(4, 5);
        B(0, 8) = -A_rhs(5, 6);
        B(0, 9) = A_rhs(4, 6) - A_rhs(5, 7);
        B(0, 10) = A_rhs(4, 7) - A_rhs(5, 8);
        B(0, 11) = A_rhs(4, 8) - A_rhs(5, 9);
        B(0, 12) = A_rhs(4, 9);

        B(1, 0) = -A_rhs(7, 0);
        B(1, 1) = A_rhs(6, 0) - A_rhs(7, 1);
        B(1, 2) = A_rhs(6, 1) - A_rhs(7, 2);
        B(1, 3) = A_rhs(6, 2);
        B(1, 4) = -A_rhs(7, 3);
        B(1, 5) = A_rhs(6, 3) - A_rhs(7, 4);
        B(1, 6) = A_rhs(6, 4) - A_rhs(7, 5);
        B(1, 7) = A_rhs(6, 5);
        B(1, 8) = -A_rhs(7, 6);
        B(1, 9) = A_rhs(6, 6) - A_rhs(7, 7);
        B(1, 10) = A_rhs(6, 7) - A_rhs(7, 8);
        B(1, 11) = A_rhs(6, 8) - A_rhs(7, 9);
        B(1, 12) = A_rhs(6, 9);

        B(2, 0) = -A_rhs(9, 0);
        B(2, 1) = A_rhs(8, 0) - A_rhs(9, 1);
        B(2, 2) = A_rhs(8, 1) - A_rhs(9, 2);
        B(2, 3) = A_rhs(8, 2);
        B(2, 4) = -A_rhs(9, 3);
        B(2, 5) = A_rhs(8, 3) - A_rhs(9, 4);
        B(2, 6) = A_rhs(8, 4) - A_rhs(9, 5);
        B(2, 7) = A_rhs(8, 5);
        B(2, 8) = -A_rhs(9, 6);
        B(2, 9) = A_rhs(8, 6) - A_rhs(9, 7);
        B(2, 10) = A_rhs(8, 7) - A_rhs(9, 8);
        B(2, 11) = A_rhs(8, 8) - A_rhs(9, 9);
        B(2, 12) = A_rhs(8, 9);

        typedef Matrix<double, 1, 5 > RVector5d;

        const RowVector4d & b11 = B.block < 1, 4 > (0, 0);
        const RowVector4d & b12 = B.block < 1, 4 > (0, 4);
        const RVector5d & b13 = B.block < 1, 5 > (0, 8);
        const RowVector4d & b21 = B.block < 1, 4 > (1, 0);
        const RowVector4d & b22 = B.block < 1, 4 > (1, 4);
        const RVector5d & b23 = B.block < 1, 5 > (1, 8);
        const RowVector4d & b31 = B.block < 1, 4 > (2, 0);
        const RowVector4d & b32 = B.block < 1, 4 > (2, 4);
        const RVector5d & b33 = B.block < 1, 5 > (2, 8);

        CConv < 4, 4 > conv44;
        CConv < 4, 5 > conv45;
        CConv < 7, 5 > conv75;
        CConv < 8, 4 > conv84;
        CConv < 5, 4 > conv54;

        const Matrix<double, 1, 11 > poly = (conv75(conv44(b11, b22), b33) - conv84(conv45(b11, b23), b32)) +
                                            (conv84(conv45(b12, b23), b31) - conv75(conv44(b12, b21), b33)) +
                                            (conv84(conv54(b13, b21), b32) - conv84(conv54(b13, b22), b31));

        Matrix<double, 3, 3 > bt;
        SOLS.setZero();
        
        const bool bUseOpencvPolysolve = false, bUseNETLIB_RPOLY=true;
        if(bUseOpencvPolysolve)
        {
#ifndef USE_OLD_OPENCV
            //Doesn't converge
            const int iters=30000;
            
            {
                cv::Mat coeffs3(3,1,CV_64FC1), roots3(2,1,CV_64FC2);
                coeffs3.at<double>(0) = 1;
                coeffs3.at<double>(1) = 0;
                coeffs3.at<double>(2) = 1;
                cv::solvePoly(coeffs3,roots3, iters);

                cout << roots3 << endl;
            }            
              
            cv::Mat coeffs(11,1,CV_64FC1), roots(10,1,CV_64FC2);
            for (int i = 0; i < 11; i++)
                coeffs.at<double>(i) = poly(i);
//                coeffs.at<cv::Vec2d>(i) = cv::Vec2d(poly(i), 0);
            
            cout << coeffs << endl;
            

            solvePoly2(coeffs,roots, iters);
            
            for (int i = 0; i < 10; i++)
            {
                cv::Vec2d root=roots.at<cv::Vec2d>(i);
                //if(root(1) < 0)
                  //  root(1) = 0; //Try one of each pair of imaginary roots too
                eigvals(i) = std::complex<double>(root(0), root(1));
            }
            //cout << roots << endl;
            std::cout << eigvals.transpose() << endl;
#else
                        THROW("Opencv Polysolve won't work in old opencv (doesn't work anyway)");
#endif
        }
        else if(bUseNETLIB_RPOLY)
        {
            int degree=10;
            Eigen::Matrix<double, 10, 1> real,imag;
            rpoly_ak1(poly.data(), &degree, real.data(), imag.data());
            
            std::vector<std::complex<double> > vRoots(10);
            
            for (int i = 0; i < 10; i++)
                eigvals(i) = std::complex<double>(real(i), imag(i));

            //sortComplex(eigvals);  std::cout << eigvals.transpose() << endl;

        } 
        else
        {
#ifdef EIGEN_POLY
            
            Matrix<double, 1, 11 > poly_reversed;
            for (int i = 0; i < 11; i++)
                poly_reversed(i) = poly(10 - i);

            PolynomialSolver<double, 10> psolve(poly_reversed);

            eigvals = psolve.roots();
            //sortComplex(eigvals); std::cout << eigvals.transpose() << endl << endl;
#else
            THROW("Need to #define EIGEN_POLY (not really necessary--is slower than netlib)");
#endif
        }

        for (int i = 0; i < 10; i++) {
            if (eigvals(i).imag() == 0) {
                double root = eigvals(i).real();
                double root_2 = sqr(root);
                double root_3 = root*root_2;
                double root_4 = sqr(root_2);

                bt(0, 0) = (b11(0) * root_3) + (b11(1) * root_2) + b11(2) * root + b11(3);
                bt(0, 1) = (b12(0) * root_3) + (b12(1) * root_2) + b12(2) * root + b12(3);
                bt(0, 2) = (b13(0) * root_4) + (b13(1) * root_3) + (b13(2) * root_2) + b13(3) * root + b13(4);
                bt(1, 0) = (b21(0) * root_3) + (b21(1) * root_2) + b21(2) * root + b21(3);
                bt(1, 1) = (b22(0) * root_3) + (b22(1) * root_2) + b22(2) * root + b22(3);
                bt(1, 2) = (b23(0) * root_4) + (b23(1) * root_3) + (b23(2) * root_2) + b23(3) * root + b23(4);
                bt(2, 0) = (b31(0) * root_3) + (b31(1) * root_2) + b31(2) * root + b31(3);
                bt(2, 1) = (b32(0) * root_3) + (b32(1) * root_2) + b32(2) * root + b32(3);
                bt(2, 2) = (b33(0) * root_4) + (b33(1) * root_3) + (b33(2) * root_2) + b33(3) * root + b33(4);

#if EIGEN_VERSION_AT_LEAST(2,90,0)
                Eigen::JacobiSVD<Eigen::Matrix3d > svdE(bt, Eigen::ComputeFullV);
#else
                Eigen::SVD<Eigen::Matrix3d > svdE(E);
#endif

                Vector3d xy1;
                for (int j = 0; j < 3; j++)
                    xy1(j) = svdE.matrixV()(j, 2); // V_trans(2,j); //if(IS_DEBUG) CHECK

                if (fabs(xy1(2)) > 0) {
                    double dScale = 1.0 / xy1(2);
                    SOLS(0, i) = xy1(0) * dScale;
                    SOLS(1, i) = xy1(1) * dScale;
                }
                SOLS(2, i) = root; //in matlab this is complex but the complex rows are never used
                SOLS(3, i) = 1;
            }
        }
        //cout << "SOLS\n" << SOLS ;
    } else /*GB stuff*/ {

        //Uses Grobner basis version

        const int anColIndices[6] = {0, 1, 2, 4, 5, 7};

        Matrix10d M; // init to 0
        M.setZero();

        for (int i = 0; i < 6; i++) {
            M.row(i) = A_rhs.row(anColIndices[i]);

            DEBUGONLY(if (std::isnan(M(i, 0)))
                    cout << A_rhs << endl;)
                if(IS_DEBUG) CHECK(std::isnan(M(i, 0)), "NaN detected");
        }
        M *= -1;
        M(6, 0) = 1;
        M(7, 1) = 1;
        M(8, 3) = 1;
        M(9, 6) = 1;
        //cout << "\nM=\n" << M;

        Eigen::EigenSolver<Matrix10d> eigenVV(M);

        eigvals = eigenVV.eigenvalues(); // both complex
#if EIGEN_VERSION_AT_LEAST(2,90,0)
        const Eigen::EigenSolver<Matrix10d>::EigenvectorsType & eigenVecs = eigenVV.eigenvectors();
#else
        const Eigen::EigenSolver<Matrix10d>::EigenvectorType & eigenVecs = eigenVV.eigenvectors();
#endif

        SOLS = eigenVecs.block < 4, 10 > (6, 0).real();
    }

    const Matrix<double, 9, 10 > & Evec /*(9, 10)*/ = EE*SOLS;

    for (int i = 0; i < 10; i++) {
        if (0 == eigvals(i).imag()) {
            bool bUprightnessIsOK = dUprightThresh < 0;

            Matrix3d E;
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++) {
                    E(r, c) = Evec(c + 3 * r, i);
                }

            // do normalisation:
            double dLen = sqrt(E.squaredNorm());
            double dLenInv = M_SQRT2 / dLen;
            E *= dLenInv;

            if (!bUprightnessIsOK) {
                bUprightnessIsOK = upright(E, dUprightThresh); //Todo: could modify E a bit to make upright? but then shouldn't be using full-on E estimation anyway
                //cout << (bUprightnessIsOK ? "E is upright" : "E is not upright") << endl;
            }

            if (bUprightnessIsOK) {
                
                //makeClosestE(E);
                //double dErrSq = getResids(E, anHypSet, m0, m1);
                double dErrSq = getETEResid(E);
                
                if (dErrSq > 0.00001) 
                    FAIL_COUNT++;
                else
                {
                    nMats++;
                    C3x3MatModel & model = dynamic_cast<C3x3MatModel &> (pdEssentialMat.addModel());
                    int j = 0;
                    for (int r = 0; r < 3; r++)
                        for (int c = 0; c < 3; c++) {
                            model[j] = E(r, c);
                            j++;
                        }
                }
            }
        }
    }

    //cout << nMats << " found\n";
    if(IS_DEBUG) CHECK(nMats != pdEssentialMat.numModels(), "Counting models failed");

    return nMats;
}

#include "geom/geom.h"
#include "util/random.h"
#include "geom/geom_eigen.h"


int C5ptEssentialMat::getModels_int(const TSubSet & anHypSet, CModels & aModels) {
    //if(IS_DEBUG) CHECK((int)p1.size() != nPoints || (int)p2.size() != nPoints, "C5ptEssentialMat::getModels_int: Bad number of points" );
    int temp;
    return calcEssentialMat_5point_Eigen(anHypSet, p1, p2, dynamic_cast<T3x3MatModels &> (aModels), dUprightThresh, temp, true);
}

int C7ptEssentialMat_GS::getModels_int(const TSubSet & anHypSet, CModels & aModels) {
    int temp;
    return calcEssentialMat_7point(anHypSet, p1, p2, dynamic_cast<T3x3MatModels &> (aModels), temp);
}

int CMCEssentialMat::getModels_int(const TSubSet & anHypSet, CModels & aModels) {
    for (int j = 0; j < nHypotheses; j++) {

        C3dRotation R(CRandom::Normal(), CRandom::Normal(), CRandom::Normal(), CRandom::Normal());
        C3dPoint t(CRandom::Normal(), CRandom::Normal(), CRandom::Normal());
        t.normalise();

        Eigen::Matrix3d E;
        makeE(R, t, E);
        C3x3MatModel & model = static_cast<C3x3MatModel &> (aModels.addModel());
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                model(r, c) = E(r, c);
    }

    return nHypotheses;
}

int C2ptMCEssentialMat::getModels_int(const TSubSet & anHypSet, CModels & aModels) {
    for (int j = 0; j < nHypotheses; j++) {

        C3dRotation R(CRandom::Normal(), CRandom::Normal(), CRandom::Normal(), CRandom::Normal());

        ARRAY(C3dPoint, x, nPoints);
        ARRAY(C3dPoint, xp, nPoints);
        
        for (int i = 0; i < nPoints; i++) {
            x[i] = C3dPoint(C2dPoint(p1[anHypSet[i]]), 1);
            xp[i] = C3dPoint(C2dPoint(p2[anHypSet[i]]), 1);

            x[i].rotate(R);
        }
        
        C3dPoint t = crossproduct(crossproduct(x[0], xp[0]), crossproduct(x[1], xp[1]));

        Eigen::Matrix3d E;
        makeE(R, t, E);
        C3x3MatModel & model = static_cast<C3x3MatModel &> (aModels.addModel());
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                model(r, c) = E(r, c);
    }


    return nHypotheses;
}