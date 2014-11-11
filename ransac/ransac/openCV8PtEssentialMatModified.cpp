/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCV8PtEssentialMatModified.cpp
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#include <Eigen/SVD>

#include "openCV8PtEssentialMatModified.h"
#include "util/opencv.h"
#include "ransac.h"

using namespace std;

COpenCV8PtEssentialMatModified::COpenCV8PtEssentialMatModified(const T2dPoints & p1, const T2dPoints & p2, double d8ptCutoff, bool bFindF) : CImCorrModelRefiner(8, p1, p2), d8ptCutoff(d8ptCutoff), bFindF(bFindF) {
    //if(IS_DEBUG) CHECK((d8ptCutoff > 0) == bFindF, "8pt cutoff irrelevant (set negative) when finding F")
}

COpenCV8PtEssentialMatModified::~COpenCV8PtEssentialMatModified() {

}

bool SVsEqual(double sv1, double sv2, double dEssential8ptCutoff) {
    if(IS_DEBUG) CHECK(sv1<sv2 /*|| zero(sv2)*/, "SVsEqual: Supply nonzero SVs in sorted order");
    if (!zero(sv2)) {
        bool bNearlyEqual = (sv1 / sv2 < dEssential8ptCutoff); //sv1 less than 4.5*sv2 sv2==0 => degenerate config
        if (!bNearlyEqual)
            cout << "SV ratio test failed\n";
        return bNearlyEqual;
    }
    return false;
}

bool BoW_EssentialMatrix_8Point(const T2dPoints & apprime, const T2dPoints & ap, CMask & mask, C3x3MatModel & fmatrix, bool bEssential, double dEssential8ptCutoff) HOT;

bool BoW_EssentialMatrix_8Point(const T2dPoints & apprime, const T2dPoints & ap, CMask & mask, C3x3MatModel & fmatrix, bool bEssential, double dEssential8ptCutoff) {
    const int count = (int) ap.size();

    const int nInlierCount = mask.countInliers();

    if(IS_DEBUG) CHECK(nInlierCount < 8, "Too few points for 8pt algorithm");

    int result = 0;
    CvMat* A = 0;

    double w[9], v[9 * 9];
    CvMat W = cvMat(1, nInlierCount == 8 ? 8 : 9, CV_64F, w);
    CvMat V = cvMat(9, 9, CV_64F, v);
    CvMat U, F0, TempF;

    int i, good_count = 0;
    CvPoint2D64f m0c = {0, 0}, m1c = {0, 0};
    double t, scale0 = 0, scale1 = 0;
    double* a;
    int a_step;

    CV_FUNCNAME("icvFMatrix_8Point");

    //__BEGIN__;

    //assert( m0 && m1 && fmatrix );

    try {
        // compute centers and average distances for each of the two point sets
        for (i = 0; i < count; i++)
            if (mask[i]) {
                double x = ap[i].getX(), y = ap[i].getY();
                m0c.x += x;
                m0c.y += y;

                x = apprime[i].getX(), y = apprime[i].getY();
                m1c.x += x;
                m1c.y += y;
                good_count++;
            }

        if (good_count < 8) goto todo__END__;

        // calculate the normalizing transformations for each of the point sets:
        // after the transformation each set will have the mass center at the coordinate origin
        // and the average distance from the origin will be ~sqrt(2).
        // Todo: if(bEssential) we probably dont need to do this. Is also slow!
        t = 1. / good_count;
        m0c.x *= t;
        m0c.y *= t;
        m1c.x *= t;
        m1c.y *= t;

        for (i = 0; i < count; i++)
            if (mask[i]) {
                double x = ap[i].getX() - m0c.x, y = ap[i].getY() - m0c.y;
                scale0 += sqrt(x * x + y * y);

                x = fabs(apprime[i].getX() - m1c.x), y = fabs(apprime[i].getY() - m1c.y);
                scale1 += sqrt(x * x + y * y);
            }

        scale0 *= t;
        scale1 *= t;

        if (scale0 < FLT_EPSILON || scale1 < FLT_EPSILON) goto todo__END__;

        scale0 = M_SQRT2 / scale0;
        scale1 = M_SQRT2 / scale1;

        CV_CALL(A = cvCreateMat(good_count, 9, CV_64F));
        a = A->data.db;
        a_step = A->step / sizeof (a[0]);

        // form a linear system: for each selected pair of points m0 & m1,
        // the row of A(=a) represents the equation: (m1, 1)'*F*(m0, 1) = 0
        for (i = 0; i < count; i++) {
            if (mask[i]) {
                double x0 = (ap[i].getX() - m0c.x) * scale0;
                double y0 = (ap[i].getY() - m0c.y) * scale0;
                double x1 = (apprime[i].getX() - m1c.x) * scale1;
                double y1 = (apprime[i].getY() - m1c.y) * scale1;

                a[0] = x1 * x0;
                a[1] = x1 * y0;
                a[2] = x1;
                a[3] = y1 * x0;
                a[4] = y1 * y0;
                a[5] = y1;
                a[6] = x0;
                a[7] = y0;
                a[8] = 1;
                a += a_step;
            }
        }

        cvSVD(A, &W, 0, &V, CV_SVD_MODIFY_A + CV_SVD_V_T);

        for (i = 0; i < 8; i++) {
            if (fabs(w[i]) < FLT_EPSILON) break;
        }

        if (i < 7) goto todo__END__;

        F0 = cvMat(3, 3, CV_64F, v + 9 * 8); // take the last column of v as a solution of Af = 0

        // make F0 singular (of rank 2) by decomposing it with SVD,
        // zeroing the last diagonal element of W and then composing the matrices back.

        // use v as a temporary storage for different 3x3 matrices
        W = U = V = TempF = F0;
        W.data.db = v;
        U.data.db = v + 9;
        V.data.db = v + 18;
        TempF.data.db = v + 27;

        // apply the transformation that is inverse
        // to what we used to normalize the point coordinates
        {
            double tt0[] = {scale0, 0, -scale0 * m0c.x, 0, scale0, -scale0 * m0c.y, 0, 0, 1};
            double tt1[] = {scale1, 0, -scale1 * m1c.x, 0, scale1, -scale1 * m1c.y, 0, 0, 1};
            CvMat T0, T1;
            T0 = T1 = F0;
            T0.data.db = tt0;
            T1.data.db = tt1;

            // F0 <- T1'*F0*T0
            cvGEMM(&T1, &F0, 1., 0, 0., &TempF, CV_GEMM_A_T);
            //F0.data.db = fmatrix.asDouble9();

            cvGEMM(&TempF, &T0, 1., 0, 0., &F0, 0);

            // make F(3,3) = 1
            /*if( fabs(F0.data.db[8]) > FLT_EPSILON )
             cvScale( &F0, &F0, 1./F0.data.db[8] );*/
        } //**** Moved from...

        Eigen::Matrix3d E(F0.data.db), Diag110;

        //cout << "E before svd\n" << E/E.norm() << endl;


        Diag110.setIdentity();
        Diag110(2, 2) = 0;

        Eigen::JacobiSVD<Eigen::Matrix3d> svd(E, Eigen::ComputeFullU + Eigen::ComputeFullV);
        //cout << "svds\n" << svd.singularValues().transpose() << endl;

        if (!bEssential) {
            Diag110(0, 0) = svd.singularValues()(0);
            Diag110(1, 1) = svd.singularValues()(1);
        }

        E = svd.matrixU() * Diag110 * svd.matrixV().transpose();
        
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                fmatrix(r, c) = E(r, c);

        //cout << "E after svd\n" << E/E.norm() << endl;
        /*

        cvSVD(&F0, &W, &U, &V, CV_SVD_MODIFY_A + CV_SVD_U_T + CV_SVD_V_T);

        //SVs not a good measure of ER uncertainty as depends on whether we have picked noisy inliers (removing inliers until we have 8 will improve certainty)
        cout << good_count << " inliers, sv's = " << W.data.db[8] << ", " << W.data.db[4] << ", " << W.data.db[0] << "\n ";

        W.data.db[8] = 0.;

        /////// Make Essential Matrix: (After undoing normalisation)
        if (bEssential) {
            if (!SVsEqual(W.data.db[0], W.data.db[4], dEssential8ptCutoff)) {
                cout << "WARNING: 8-pt alg failed because " << W.data.db[0] << "and " << W.data.db[4] << " are not close enough (? trying anyway)\n";
                result = 0;
                //EXIT;
            }

            W.data.db[4] = 1.;
            W.data.db[0] = 1.;
        }

        // F0 <- U*diag([1, 1, 0])*V'
        cvGEMM(&U, &W, 1., 0, 0., &TempF, CV_GEMM_A_T);
        cvGEMM(&TempF, &V, 1., 0, 0., &F0, CV_GEMM_B_T);

        cout << "E using old method\n" << Eigen::Matrix3d(F0.data.db) << endl;
        
        ****...here
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                fmatrix(r, c) = cvmGet(&F0, r, c);**/

        result = 1;
    } catch (CException pEx) {
        cvReleaseMat(&A);
        throw pEx;
    }

todo__END__:
    ;
exit:
    ;

    cvReleaseMat(&A);
    return result != 0;
}

bool COpenCV8PtEssentialMatModified::fitModel_int(CMask & mask, CModel & model, bool bVerbose) {
    return BoW_EssentialMatrix_8Point(p1, p2, mask, dynamic_cast<C3x3MatModel &> (model), !bFindF, d8ptCutoff);
}
