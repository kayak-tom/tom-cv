/* 
 * Nonlinear refinement of rotation and translation, on a set of correspondences
 */

#include "geom/geom_eigen.h"
#include "refineEOnRTManifold.h"
#include <Eigen/Dense>
#include "makeBasis_GramSchmidt.h"
#include "time/SpeedTest.h"

#define SEGMENT_T head<3>()
#define SEGMENT_Q tail<4>()

template<typename TM> void normalise(TM & M) {
    M /= M.norm();
}

typedef Eigen::Matrix<double, 3, 1 > T3Vec;

typedef std::vector<T3Vec> THomogPointVec;


static const double ROBUST_THRESH = 0.005, ROBUST_THRESH_INV = 1.0 / ROBUST_THRESH; //Threshhold for all robust norms

class CLeastSquares {
public:

    inline double operator()(const double) const {
        return 1;
    }
};

/*class CHuber_TorrMurray {
public:

    inline double operator()(const double d) const {
        const double b = 0.02;

        const double d_abs = fabs(d);
        if (d_abs < b)
            return 1;
        else if (d_abs < 3 * b)
            return b / d; //I think this is possibly wrong--should use 1/root(d)
        else
            return 0; //TODO, probably best to just return SIGMA/d;
    }
};*/

class CHuber {
public:

    inline double operator()(const double d) const {
        const double b_sq = sqr(ROBUST_THRESH);

        const double d_abs = fabs(d);
        if (d_abs < ROBUST_THRESH)
            return 1;
        else
            return sqrt(2 * ROBUST_THRESH * d_abs - b_sq) / d_abs;
    }
};

class CB_Z { //Blake-Zisserman Gaussian + uniform
public:

    inline double operator()(const double d) const {
        const double SD_INV = 1.0 / (0.5 * ROBUST_THRESH);
        const double eps = exp(-sqr(ROBUST_THRESH * SD_INV)); //Equally likely to be inlier or outlier at thresh
        const double d_abs = fabs(d) + 1e-12;
        const double zeroPoint = log(1 + eps); //Needed for LS...

        const double dCost_sq = zeroPoint - log(exp(-sqr(d * SD_INV)) + eps);
        if(IS_DEBUG) CHECK(dCost_sq < 0, "Cost computation failed");

        return sqrt(dCost_sq) / d_abs; //TODO: should this sqrt be here?! Cross-ref g2oEdges.h
    }
};

class CPseudoHuber {
public:

    inline double operator()(const double d) const {
        const double b_sq = sqr(ROBUST_THRESH);
        const double d_abs = fabs(d) + 1e-12;

        //C(delta) = 2*b^2*(sqrt(1+(delta/b)^2) - 1);
        return sqrt(2 * b_sq * (sqrt(1 + sqr(d * ROBUST_THRESH_INV)) - 1)) / d_abs;
    }
};

class CL1 {
public:

    inline double operator()(const double d) const {
        return 1 / (sqrt(fabs(d)) + 1e-12);
    }
};

template<class CRho>
class CRefineEOnManifold {
    static const int NUM_PARAMS = 7;
    static const int TANGENTS_DIMS_Q = 3;
    static const int TANGENTS_DIMS_T = 2;
    static const int N_PARAMS_TS = TANGENTS_DIMS_Q + TANGENTS_DIMS_T;
    static const int NUM_MODEL_VARS = Eigen::Dynamic;
    typedef Eigen::Matrix<double, NUM_PARAMS, 1 > TParamVector;
    typedef Eigen::Matrix<double, N_PARAMS_TS, 1 > TParamVectorTS;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS, 1 > TResidVector;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS, N_PARAMS_TS> TJMatrix;
    typedef Eigen::Matrix<double, N_PARAMS_TS, N_PARAMS_TS> TJTJMatrix;
    typedef Eigen::Matrix<double, 3, 3 > TE;

    typedef Eigen::Matrix<double, 3, 1 > TTangentAtT;
    typedef Eigen::Matrix<double, 4, 1 > TTangentAtQ;

    static void computeDerivativesNumerically(const TParamVector & params, const TResidVector & residuals,
            const THomogPointVec & p1, const THomogPointVec & p2, const TTangentAtQ * aTangentsAtQ, const TTangentAtT * aTangentsAtT, const CMask & mask, TJMatrix & J) {
        TResidVector resids_plus(residuals.size());

        double delta = 1e-7;
        for (int nParam = 0; nParam < TANGENTS_DIMS_T; nParam++) {
            TParamVector params_plus = params;
            params_plus.SEGMENT_T += delta * aTangentsAtT[nParam];
            computeResiduals(mask, p1, p2, params_plus, resids_plus);

            J.col(nParam) = (resids_plus - residuals).array() / delta;
        }

        for (int nParam = TANGENTS_DIMS_T; nParam < N_PARAMS_TS; nParam++) {
            TParamVector params_plus = params;
            params_plus.SEGMENT_Q += delta * aTangentsAtQ[nParam - TANGENTS_DIMS_T];
            computeResiduals(mask, p1, p2, params_plus, resids_plus);

            J.col(nParam) = (resids_plus - residuals).array() / delta;
        }

        //cout << J << endl << endl;
    }

    static void RTtoParams(const C3dRotation & q, const Eigen::Vector3d & t, TParamVector & params) {
        for (int i = 0; i < 3; i++)
            params(i) = t(i);
        for (int i = 0; i < 4; i++)
            params(i + 3) = q[i];
    }

    static void paramsToRT(const TParamVector & params, C3dRotation & q, Eigen::Vector3d & t) {
        t = params.SEGMENT_T;
        q = C3dRotation(params(3), params(4), params(5), params(6));
    }

    /*inline static void setVars(const TParamVector & params, TFloat & tt1, TFloat & tt2, TFloat & tt3, TFloat & s1, TFloat & s2, TFloat & s3, TFloat & s4) const {
        tt1 = params(0);
        tt2 = params(1);
        tt3 = params(2);
        s1 = params(3);
        s2 = params(4);
        s3 = params(5);
        s4 = params(6);
    }*/
    template<int SIZE>
    static void makeOrthonormalBasis(const Eigen::Matrix<double, SIZE, 1 > & q, Eigen::Matrix<double, SIZE, 1 > * aTangentsAtQ) {
        for (int i = 0; i < SIZE - 1; i++) {
            aTangentsAtQ[i] -= q.dot(aTangentsAtQ[i]) * q;
            for (int j = 0; j < i; j++) {
                aTangentsAtQ[i] -= aTangentsAtQ[j].dot(aTangentsAtQ[i]) * aTangentsAtQ[j];
            }
            normalise(aTangentsAtQ[i]);
        }
    }

    static void makeTangentVecs(const TParamVector & params, TTangentAtQ * aTangentsAtQ, TTangentAtT * aTangentsAtT) {
        TTangentAtT t = params.SEGMENT_T;
        makeOrthonormalBasis < 3 > (t, aTangentsAtT);
        TTangentAtQ q = params.SEGMENT_Q;
        makeOrthonormalBasis < 4 > (q, aTangentsAtQ);
    }

    template<int SIZE>
    static void makeRandom(Eigen::Matrix<double, SIZE, 1 > * aTangentsAtQ) {
        for (int i = 0; i < SIZE - 1; i++) {
            aTangentsAtQ[i].setRandom();
        }
    }

    static void makeRandom(TTangentAtQ * aTangentsAtQ, TTangentAtT * aTangentsAtT) {
        makeRandom < 3 > (aTangentsAtT);
        makeRandom < 4 > (aTangentsAtQ);
    }

    template<bool bSquared>
    static double sampsonsErr(const T3Vec & x, const T3Vec & xp, TE & E) {
        T3Vec xpE = xp.transpose() * E;
        T3Vec Ex = E * x;

        double numerator = xpE.dot(x);
        double denom = xpE.segment(0, 2).squaredNorm() + Ex.segment(0, 2).squaredNorm();

        if (bSquared)
            return sqr(numerator) / denom;
        else
            return numerator / sqrt(denom);
    }

    static void computeResiduals(const CMask & mask, const THomogPointVec & p1, const THomogPointVec & p2, const TParamVector & params, TResidVector & residuals) {
        CRho rho;

        TE E_temp;
        C3dRotation q_temp;
        Eigen::Vector3d t_temp;
        paramsToRT(params, q_temp, t_temp);
        makeE(q_temp, t_temp, E_temp); //Don't need to re-normalise here most of the time...

        const int NUM_RESIDS = (int) p1.size();
        for (int i = 0; i < NUM_RESIDS; i++) {
            if (mask[i]) {
                double dSE = sampsonsErr < false > (p1[i], p2[i], E_temp);
                residuals(i) = dSE * rho(dSE);
                CHECKNAN(residuals(i));
                /*if (IS_DEBUG && isnan(residuals(i))) {
                    cout << dSE << "=SE, rho=" << rho(dSE) << endl;
                    for (double d = -0.1; d < 0.1; d += 0.002)
                        cout << d << "\t" << rho(d) << endl;

                    THROW("residual is nan");
                }*/
            } else
                residuals(i) = 0;
        }
    }
public:

    static double refine(const THomogPointVec & p1, const THomogPointVec & p2, C3dRotation & q_in, Eigen::Vector3d & t_in, CMask & mask) {
        if(IS_DEBUG) CHECK(mask.countInliers()==0, "No points to use");
        CHECK(t_in.squaredNorm() == 0, "T is zero, should be set to something nonzero even if we have pure rotation");
        
        double dBestResid = HUGE;
        const int NUM_ITERS = 1;
        C3dRotation q_best = q_in;
        Eigen::Vector3d t_best = t_in;
        CStopWatch s;
        s.startTimer();
        double dRandStartParam = 0.05; //Don't want over-adjustement
        for (int i = 0; i < NUM_ITERS; i++) {
            C3dRotation q = q_in;
            Eigen::Vector3d t = t_in;

            if (i > 0) {
                C3dRotation qrand;
                qrand.setRandom(dRandStartParam);
                q = q*qrand;

                Eigen::Vector3d trand = dRandStartParam * Eigen::Vector3d::Random();
                t += trand;
                normalise(t);
            }

            double dErr = refine_int(p1, p2, q, t, mask);
            if (dErr < dBestResid) {
                dBestResid = dErr;
                q_best = q;
                t_best = t;
            }
        }
        s.stopTimer();
        double dTime = s.getElapsedTime() / NUM_ITERS;
        cout << "Time: " << dTime << endl;

        q_in = q_best;
        t_in = t_best;

        return dBestResid;
    }
    static double refine_int(const THomogPointVec & p1, const THomogPointVec & p2, C3dRotation & q, Eigen::Vector3d & t, CMask & mask) HOT;

private:

    static void linearWeight(const T3Vec & x, const T3Vec & xp, TE & E, double & w, double & r) {
        T3Vec xpE = xp.transpose() * E;
        T3Vec Ex = E * x;

        double numerator = xpE.dot(x);
        double denom = xpE.segment(0, 2).squaredNorm() + Ex.segment(0, 2).squaredNorm();

        w = 1 / (sqrt(denom) + 1e-8);
        r = numerator;

        CHECKNAN(w);
        CHECKNAN(r);
    }

    static double refineLinear1iter(const THomogPointVec & p1, const THomogPointVec & p2, Eigen::Matrix3d & E, CMask & mask, const bool bF) {
        const int NUM_POINTS = p1.size();
        const int NUM_INLIERS = mask.countInliers();
        if(IS_DEBUG) CHECK(NUM_INLIERS <= 7, "Degenerate--too few inliers");
        CRho rho;

        typedef Eigen::Matrix<double, Eigen::Dynamic, 9 > TAMat;
        TAMat A(NUM_INLIERS, 9);
        int nRow = 0;
        //Setup and weight constraint matrix:
        for (int i = 0; i < NUM_POINTS; i++) {
            if (mask[i]) {
                double dw, dr;
                linearWeight(p1[i], p2[i], E, dw, dr);
                const double dSampsonsErr = dw*dr;
                const double dgamma = rho(dSampsonsErr);
                CHECKNAN(dgamma);

                const double x0 = p1[i](0);
                const double y0 = p1[i](1);
                const double x1 = p2[i](0);
                const double y1 = p2[i](1);

                A(nRow, 0) = x1 * x0;
                A(nRow, 1) = x1 * y0;
                A(nRow, 2) = x1;
                A(nRow, 3) = y1 * x0;
                A(nRow, 4) = y1 * y0;
                A(nRow, 5) = y1;
                A(nRow, 6) = x0;
                A(nRow, 7) = y0;
                A(nRow, 8) = 1;

                A.row(nRow) *= dw*dgamma;
                nRow++;
            }
        }
        CHECKNAN(A.sum());

        Eigen::JacobiSVD<TAMat> svdA(A, Eigen::ComputeFullV);

        //Convert to F
        Eigen::Matrix<double, 9, 1 > lastCol = svdA.matrixV().col(8);
        Eigen::Matrix3d F(lastCol.data());
        F.transposeInPlace();

        CHECKNAN(lastCol.sum());

        //Project to either E or F
        E = F;
        makeClosestE(E, bF);

        return (A * lastCol).squaredNorm();

        /*double dTotalErrProjected = 0;
        const bool bVerbose = false;
        double dTotalErrUnprojected = 0;
        for (int i = 0; i < NUM_POINTS; i++) {
            dTotalErrProjected += sampsonsErr < true > (p1[i], p2[i], E);
            if (bVerbose) dTotalErrUnprojected += sampsonsErr < true > (p1[i], p2[i], F);
        }
        if (bVerbose) {
            cout << "Err projected: " << dTotalErrProjected << endl; //These are a bit high because include the full residuals from outliers
            if (!bF) cout << "Err not projected: " << dTotalErrUnprojected << endl;
        }
        return dTotalErrProjected;*/
    }
public:

    static void testNorms(const THomogPointVec & p1, const THomogPointVec & p2, Eigen::Matrix3d & E, CMask & mask) {
        CPseudoHuber pseudoHuber;
        CL1 l1;
        CHuber huber;
        CB_Z bz;

        for (int i = 0; i < (int) p1.size(); i++) {
            double dSE_sq = sampsonsErr < true > (p1[i], p2[i], E);
            double dSE = sampsonsErr < false > (p1[i], p2[i], E);
            int isUnderThresh = dSE_sq < sqr(ROBUST_THRESH);
            cout << (int) mask[i] << " " << isUnderThresh << " " << dSE_sq << " Huber weight=" << huber(dSE)
                    << " PHuber weight=" << pseudoHuber(dSE)
                    << " BZ weight=" << bz(dSE)
                    << " L1 weight=" << l1(dSE) << endl;
        }
    }

    //Iteratively-reweighted 8-point algorithm

    static double refineLinear(const THomogPointVec & p1, const THomogPointVec & p2, Eigen::Matrix3d & E, CMask & mask, bool bF) {

        //testNorms(p1,p2,E,mask);
        CStopWatch s;
        s.startTimer();

        double dErr = HUGE, dNewErr = HUGE;
        const double EPS = sqr(0.005);
        int nIter = 1;
        for (; nIter < 20; nIter++) {

            dNewErr = refineLinear1iter(p1, p2, E, mask, bF);
            //cout << dNewErr << " error" << endl;
            if (fabs(dErr - dNewErr) < EPS) {
                break;
            }

            dErr = dNewErr;
        }

        if (bF)
            makeClosestE(E);

        s.stopTimer();
        double dTime = s.getElapsedTime();
        cout << "Time: " << dTime << endl;

        cout << "Weighted linear refinement converged after " << nIter << " iterations, error=" << dNewErr << endl;

        return dNewErr;
    }

};

template<class CRho>
double CRefineEOnManifold<CRho>::refine_int(const THomogPointVec & p1, const THomogPointVec & p2, C3dRotation & q, Eigen::Vector3d & t, CMask & mask) {

    CHECK(mask.countInliers() < 5, "Too few inliers for refinement");


    const bool bCGD = false;

    const int NUM_RESIDS = (int) p1.size();
    TResidVector residuals(NUM_RESIDS);
    TParamVector params;

    RTtoParams(q, t, params);

    TTangentAtQ aTangentsAtQ[TANGENTS_DIMS_Q];
    TTangentAtT aTangentsAtT[TANGENTS_DIMS_T];
    makeRandom(aTangentsAtQ, aTangentsAtT);
    computeResiduals(mask, p1, p2, params, residuals);

    double lambda = 0.25, dErr = residuals.squaredNorm(), dStepLenSq = HUGE;

    const double eps = sqr(0.0001);
    int nIter = 0;
    const int MAX_ITERS = 150;
    for (; nIter < MAX_ITERS && dErr > eps; nIter++) { //Error will not actually go to 0
        TJMatrix J(NUM_RESIDS, 5);

        makeTangentVecs(params, aTangentsAtQ, aTangentsAtT);
        computeDerivativesNumerically(params, residuals, p1, p2, aTangentsAtQ, aTangentsAtT, mask, J);

        CHECKNAN(J.sum());
        TJTJMatrix JTJ;
        if (!bCGD)
            JTJ = J.transpose() * J;

        //cout << params.transpose() << endl;
        TParamVectorTS paramUpdateVec;
        for (;;) {
            if (!bCGD) {
                TJTJMatrix JTJ_weighted = JTJ;
                JTJ_weighted.diagonal().array() *= (1 + lambda);
                paramUpdateVec = JTJ_weighted.llt().solve(J.transpose() * residuals); // inv * J.transpose() * residuals;
            } else {
                THROW("Not implemented");
            }

            CHECKNAN(paramUpdateVec.sum());

            TParamVector paramsNew = params;
            for (int i = 0; i < TANGENTS_DIMS_T; i++)
                paramsNew.SEGMENT_T -= paramUpdateVec(i) * aTangentsAtT[i];
            for (int i = 0; i < TANGENTS_DIMS_Q; i++)
                paramsNew.SEGMENT_Q -= paramUpdateVec(i + TANGENTS_DIMS_T) * aTangentsAtQ[i];

            paramsNew.SEGMENT_T /= paramsNew.SEGMENT_T.norm();

            computeResiduals(mask, p1, p2, paramsNew, residuals);

            double dErrNew = residuals.squaredNorm();
            dStepLenSq = paramUpdateVec.squaredNorm();
            if (dErrNew < dErr) {
                lambda *= 0.25;
                dErr = dErrNew;
                params = paramsNew;
                break;
            } else {
                lambda *= 5;
                if (lambda > 100) {
                    if (IS_DEBUG && dStepLenSq > sqr(0.001)) {
                        cout << "ERROR lambda going huge, params: " << params.transpose() << endl << endl;
                        cout << "J: " << J << endl << endl;
                        cout << "dErr: " << dErr << " (most likely a false min?)" << endl << endl;
                    }
                    dStepLenSq = 0; //force breakout
                    break;
                }
            }
        }
        if (dStepLenSq < eps) //Step is tiny
            break;
    }

    cout << "Convergence after " << nIter << " iterations, ";
    cout << "error = " << dErr << ", ";
    cout << "RMS error = " << sqrt(dErr / residuals.size()) << endl;
    //cout << "params = " << params.transpose() << endl;

    paramsToRT(params, q, t);

    return dErr;
}

double CRefineEOnRTManifold::refineRobustOnAll(const T2dPoints & p1, const T2dPoints & p2, C3dRotation & R, C3dPoint & t) {
    CMask mask(p1.size());
    mask.setConst(true);
    return refineRobustOnMask(p1, p2, R, t, mask);
}

void pointsToEigen(const T2dPoints & p1, const T2dPoints & p2, THomogPointVec & ap1, THomogPointVec & ap2) {
    ap1.reserve(p1.size());
    ap2.reserve(p1.size());
    for (int i = 0; i < (int) p1.size(); i++) {
        //if(mask[i])
        {
            T3Vec x(p1[i].getX(), p1[i].getY(), 1);
            T3Vec xp(p2[i].getX(), p2[i].getY(), 1);

            ap1.push_back(x);
            ap2.push_back(xp);
        }
    }
}

typedef CB_Z TRobustNorm; //Performs poorly: 2%
//0.0614739typedef CHuber TRobustNorm; //Performs poorly: 14% Huber, 7% pseudoHuber
//typedef CLeastSquares TRobustNorm; //Obviously performs poorly (2%)
//typedef CL1 TRobustNorm; //From correct point: 35%
//typedef CPseudoHuber TRobustNorm;

double CRefineEOnRTManifold::refineRobustOnMask(const T2dPoints & p1, const T2dPoints & p2, C3dRotation & R, C3dPoint & t_in, CMask & mask) {
    THomogPointVec ap1, ap2;
    pointsToEigen(p1, p2, ap1, ap2);
    Eigen::Vector3d t = t_in.asVector();
    const double res = CRefineEOnManifold<TRobustNorm>::refine(ap1, ap2, R, t, mask);
    t_in = C3dPoint(t);
    return res;
}

double CRefineEOnRTManifold::refineLSOnMask(const T2dPoints & p1, const T2dPoints & p2, C3dRotation & R, C3dPoint & t_in, CMask & mask) {
    THomogPointVec ap1, ap2;
    pointsToEigen(p1, p2, ap1, ap2);
    Eigen::Vector3d t = t_in.asVector();
    const double res = CRefineEOnManifold<CLeastSquares>::refine(ap1, ap2, R, t, mask);
    t_in = C3dPoint(t);
    return res;
}

double refineWeightedLinearRobustOnAll(const T2dPoints & p1, const T2dPoints & p2, Eigen::Matrix3d & E, const bool bF) {
    CMask mask(p1.size());
    mask.setConst(true);
    return refineWeightedLinearRobustOnMask(p1, p2, E, mask, bF);
}

double refineWeightedLinearRobustOnMask(const T2dPoints & p1, const T2dPoints & p2, Eigen::Matrix3d & E, CMask & mask, const bool bF) {
    THomogPointVec ap1, ap2;
    pointsToEigen(p1, p2, ap1, ap2);
    return CRefineEOnManifold<TRobustNorm>::refineLinear(ap1, ap2, E, mask, bF);
}

double refineWeightedLinearLSOnMask(const T2dPoints & p1, const T2dPoints & p2, Eigen::Matrix3d & E, CMask & mask, const bool bF) {
    THomogPointVec ap1, ap2;
    pointsToEigen(p1, p2, ap1, ap2);
    return CRefineEOnManifold<CLeastSquares>::refineLinear(ap1, ap2, E, mask, bF);
}
/* 
 * DONE Use Marquardt
 * DONE Use Huber/robust Huber
 * TODO Use CGD
 * TODO Output stats
 */