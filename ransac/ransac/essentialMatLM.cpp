#include "essentialMatLM.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Cholesky> 
#include <Eigen/Geometry>
#include <boost/scoped_array.hpp>
#include <util/convert.h>
#include "geom/geom.h"
#include "geom/geom_eigen.h"
#include "util/random.h"
#include "time/SpeedTest.h"
#include <fstream>
#include "essentialMat_Eigen.h"

#include "LMParams.h"

//#include "opencv2/opencv.hpp"
#include "util/opencv.h"

template<typename TMat> inline void normalise(TMat & M) {
    const typename TMat::Scalar SS = M.squaredNorm();
    if (fabs(SS - 1) < 0.1)
        M *= 1.5 - 0.5 * SS; //normalise (binom approx)
    else
        M /= sqrt(SS);
}

template<typename TMat> inline void setRandomUnif01(TMat & M) {
    for (int r = 0; r < TMat::RowsAtCompileTime; r++)
        for (int c = 0; c < TMat::ColsAtCompileTime; c++)
            M(r, c) = CRandom::Uniform();
}

class CGDStatsFull : public CGDStats {
    int anResultTypeCounts[TERMINTION_TYPES];
    CStats aResultStats[TERMINTION_TYPES], aIterStats[TERMINTION_TYPES];
public:

    virtual void add(int nIters, eTerminationResult result) {
        anResultTypeCounts[result]++;
        aIterStats[result].add(nIters);
    }

    CGDStatsFull() {
        setZero(anResultTypeCounts, TERMINTION_TYPES);

        aResultStats[eSuccess].setName("Num successes/model");
        aResultStats[eFalseMin].setName("Num false minima/model");
        aResultStats[eRepeatedMin].setName("Num duplicates/model");
        aIterStats[eSuccess].setName("Iters/success    ");
        aIterStats[eFalseMin].setName("Iters/FalseMin   ");
        aIterStats[eRepeatedMin].setName("Iters/RepeatedMin");
    }

    virtual void reset() {
        for (int res = 0; res < TERMINTION_TYPES; res++)
            aResultStats[res].add(anResultTypeCounts[res]);

        setZero(anResultTypeCounts, TERMINTION_TYPES);
    }

    virtual void writeResults(CParam & alg, const char * filename) {
        std::ofstream file(filename, std::ios_base::app);
        alg.printParam(file);
        CStats::writeTSVheader("", file);
        file << endl;
        for (int res = 0; res < TERMINTION_TYPES; res++) {
            aResultStats[res].writeTSVdata(file, true);
            file << endl;
        }

        for (int res = 0; res < TERMINTION_TYPES; res++) {
            aIterStats[res].writeTSVdata(file, true);
            file << endl;
        }

        file.close();
    }
};

enum eParamMode {
    eMinimal, eMinimalDegen, eFull, eMinimalOnManifold
};

template<typename TFloat, int NUM_POINTS, eParamMode PARAM_MODE> class CLevMarForEHypothesis {
public:
    typedef Eigen::Matrix<TFloat, 3, 1 > THomoImPoint;
private:

    static const int NUM_PARAMS = (PARAM_MODE != eFull) ? NUM_POINTS : 7;
    typedef Eigen::Matrix<TFloat, NUM_PARAMS, 1 > TParamVector;
    typedef Eigen::Matrix<TFloat, NUM_POINTS, 1 > TResidVector;
    typedef Eigen::Matrix<TFloat, NUM_POINTS, NUM_PARAMS> TJMatrix;
    typedef Eigen::Matrix<TFloat, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;
    typedef Eigen::Matrix<TFloat, 3, 3 > TEssentialMat;

    inline static TFloat FastNormal() {
        return (TFloat) CRandom::FasterNormal2();
    }

    class CParamVector : public TParamVector {
    public:

        inline void setVars(TFloat & t1, TFloat & t2, TFloat & s1, TFloat & s2, TFloat & s3) const {
            if(IS_DEBUG) CHECK(NUM_PARAMS > 5, "Wrong setVars called");
            const CParamVector & params = *this;
            t1 = params(0);
            t2 = params(1);
            s1 = params(2);
            s2 = (NUM_POINTS > 3) ? params(3) : 0;
            s3 = (NUM_POINTS > 4) ? params(4) : 0;
        }

        inline void setVars(TFloat & tt1, TFloat & tt2, TFloat & tt3, TFloat & s1, TFloat & s2, TFloat & s3, TFloat & s4) const {
            if (PARAM_MODE == eFull) {
                const CParamVector & params = *this;
                if(IS_DEBUG) CHECK(!zero(params.squaredNorm() - 2), "Params not normalised");
                tt1 = params(0);
                tt2 = params(1);
                tt3 = params(2);
                s1 = params(3);
                s2 = params(4);
                s3 = params(5);
                s4 = params(6);
            } else {
                setVars(tt1, tt2, s1, s2, s3);
                s4 = 1;
                tt3 = 1;
            }
        }

        inline static const TFloat approxInvRoot(const TFloat f) {
            if (fabs(f - 1) < 0.1)
                return 1.5 - 0.5 * f;
            return 1.0f / std::sqrt((float)f);
        }

        void Normalise() {
            if (PARAM_MODE != eFull)
                return;

            /* Why doesn't any of this work in gcc?
            //Eigen::Matrix<TFloat,3,1> & t = block<3,1>(0,0);
            EigenNormalise(TParamVector::block<3,1>(0,0));
            //TParamVector::block<3,1>(0,0) /= (TParamVector::block<3,1>(0,0)).stableNorm();
			
            //Eigen::Matrix<TFloat,4,1> & q = block<4,1>(3,0);
            EigenNormalise(TParamVector::block<4,1>(3,0));
            //TParamVector::block<4,1>(3,0) /= (TParamVector::block<4,1>(3,0)).stableNorm();
             */
            TParamVector & params = *this;

            // = 1.0 / (TParamVector::block(0, 0, 3, 1)).stableNorm();
            // = 1.0 / (TParamVector::block(3, 0, 4, 1)).stableNorm();
            const TFloat length_t_sq = TParamVector::block(0, 0, 3, 1).squaredNorm();
            const TFloat length_s_sq = TParamVector::block(0, 0, 3, 1).squaredNorm();
            const TFloat length_t_inv = approxInvRoot(length_t_sq);
            const TFloat length_s_inv = approxInvRoot(length_s_sq);

            for (int i = 0; i < NUM_POINTS; i++)
                params(i) *= (i < 2) ? length_t_inv : length_s_inv;

            if(IS_DEBUG) CHECK(!zero(TParamVector::squaredNorm() - 2), "Params not Normalised");
        }

        void randomise() {
            TParamVector & params = *this;
            if (PARAM_MODE != eFull) {
                for (int i = 0; i < 2; i++)
                    params(i) = FastNormal();

                const bool bRandomDirection = false;

                if (bRandomDirection) {
                    TFloat t3 = FastNormal();
                    for (int i = 0; i < 2; i++)
                        params(i) /= (params(0) + t3); //So that vector is prop. to t1 t2 t3
                }

                for (int i = 2; i < NUM_POINTS; i++)
                    params(i) = FastNormal();

                if (bRandomDirection) {
                    TFloat s4 = FastNormal();
                    for (int i = 2; i < NUM_POINTS; i++)
                        params(i) /= params(2) + s4;
                }

                //params.block<NUM_POINTS-2,1>(2,0) /= (params(2)+s4); //So that vector is prop. to s1 s2 s3 s4
            } else {
                for (int i = 0; i < NUM_PARAMS; i++)
                    params(i) = FastNormal();

                Normalise();
            }
        }

        void operator=(const TParamVector & in) {
            ((TParamVector &) * this) = in;
        }
    };

    //This class represents the state needed (basis vectors, etc) for opt. on a manifold.

    class CEForOptOnManifold {
        Eigen::Matrix<TFloat, 3, 1 > t, tau[2];
        Eigen::Matrix<TFloat, 4, 1 > q, chi[3];
    public:

        CEForOptOnManifold(CParamVector & initialParams) {
            if (PARAM_MODE != eMinimalOnManifold)
                return;

            for (int i = 0; i < 3; i++)
                setRandomUnif01(chi[i]);
            //chi[i].setRandom();


            TFloat t1, t2, s1, s2, s3;
            initialParams.setVars(t1, t2, s1, s2, s3);

            t << t1, t2, 1 - t1;
            q << s1, s2, s3, 1 - s1; // (undo move from 7D to 5D)

            reset(initialParams, true);
        }

        void computeBasis() {
            Eigen::Matrix<TFloat, 3, 1 > b_use(1, 0, 0);
            if (fabs(t.transpose() * b_use) > 0.95) //Need a unit vector NOT PLL to t
                b_use = Eigen::Matrix<TFloat, 3, 1 > (0, 1, 0);

            tau[0] = b_use.cross(t);
            //tau[0] /= tau[0].stableNorm();
            normalise(tau[0]);
            tau[1] = tau[0].cross(t);

            for (int i = 0; i < 3; i++) {
                //chi[i].setRandom();
                chi[i] -= ((TFloat) (chi[i].transpose() * q)) * q;
                for (int j = 0; j < i; j++)
                    chi[i] -= ((TFloat) (chi[i].transpose() * chi[j])) * chi[j];

                //chi[i] /= chi[i].stableNorm();
                normalise(chi[i]);

                //cout << chi[i].transpose() * q << endl;
                if(IS_DEBUG) CHECK(!zero(((TFloat) (chi[i].transpose() * q))*0.001), "Not perp");
            }
            if(IS_DEBUG) CHECK(!zero(((TFloat) (chi[0].transpose() * chi[2]))*0.001), "Not perp");
        }

        void updateE(CParamVector & params) {
            for (int i = 0; i < 2; i++)
                t += params(i) * tau[i];
            for (int i = 0; i < 3; i++)
                q += params(2 + i) * chi[i];

            reset(params, true);
        }

        void reset(CParamVector & params, const bool bCB) {
            Normalise();
            if (bCB) computeBasis();
            params.setZero();
        }

        void updateEToVars(const CParamVector & params, TFloat & tt1, TFloat & tt2, TFloat & tt3, TFloat & s1, TFloat & s2, TFloat & s3, TFloat & s4) const {
            Eigen::Matrix<TFloat, 3, 1 > t_temp = t;
            Eigen::Matrix<TFloat, 4, 1 > q_temp = q;
            for (int i = 0; i < 2; i++)
                t_temp += params(i) * tau[i];
            for (int i = 0; i < 3; i++)
                q_temp += params(2 + i) * chi[i];

            tt1 = t_temp(0);
            tt2 = t_temp(1);
            tt3 = t_temp(2);
            s1 = q_temp(0);
            s2 = q_temp(1);
            s3 = q_temp(2);
            s4 = q_temp(3);
        }

        void setVars(TFloat & tt1, TFloat & tt2, TFloat & tt3, TFloat & s1, TFloat & s2, TFloat & s3, TFloat & s4) const {
            tt1 = t(0);
            tt2 = t(1);
            tt3 = t(2);
            s1 = q(0);
            s2 = q(1);
            s3 = q(2);
            s4 = q(3);
        }

        void setTau(int nBasisVec, TFloat & tau1, TFloat & tau2, TFloat & tau3) {
            tau1 = tau[nBasisVec](0);
            tau2 = tau[nBasisVec](1);
            tau3 = tau[nBasisVec](2);
        }

        void setChi(int nBasisVec, TFloat & chi1, TFloat & chi2, TFloat & chi3, TFloat & chi4) {
            chi1 = chi[nBasisVec](0);
            chi2 = chi[nBasisVec](1);
            chi3 = chi[nBasisVec](2);
            chi4 = chi[nBasisVec](3);
        }

        void Normalise() {
            normalise(t);
            normalise(q);
            //t /= t.stableNorm();
            //q /= q.stableNorm();
        }

        const Eigen::Matrix<TFloat, 4, 1 > & getq() const {
            return q;
        }

        const Eigen::Matrix<TFloat, 3, 1 > & gett() const {
            return t;
        }
    };

    //Essential matrix as a rotation and translation

    class CEssentialMatRT {
        C3dRotation q;
        C3dPoint t;

        void makeE_fixed(const C3dRotation & q_before, TEssentialMat & EAtParams) const {
            Eigen::Vector3d vt;
            t.asVector(vt);
            Eigen::Matrix3d rot, X_mat;
            Xmat(vt, X_mat);
            C3dRotation q_total = q_before*q;
            q_total.asMat(rot);
            EAtParams = (rot * X_mat).cast<TFloat > ();
            if (EAtParams(2, 2) < 0)
                EAtParams *= (TFloat) - 1; //So we can detect duplicates more easily. Probably now unneeded
        }

    public:

        void print() const {
            cout << q << " - " << t << endl;
        }

        CEssentialMatRT() {
        }

        CEssentialMatRT(const CEForOptOnManifold & E) : q(E.getq()), t(E.gett()) {

        }

        CEssentialMatRT(const CParamVector & params) {
            if (PARAM_MODE == eMinimal) {
                TFloat t1, t2, s1, s2, s3;
                params.setVars(t1, t2, s1, s2, s3);

                t = C3dPoint(t1, t2, 1 - t1);
                t.normalise();
                q = C3dRotation(s1, s2, s3, 1 - s1);
            } else if (PARAM_MODE == eMinimalDegen || PARAM_MODE == eFull) {
                TFloat tt1, tt2, tt3, s1, s2, s3, s4;
                params.setVars(tt1, tt2, tt3, s1, s2, s3, s4);

                t = C3dPoint(tt1, tt2, tt3);
                q = C3dRotation(s1, s2, s3, s4);
            } else
                THROW("Don't use this c'tor");
        }

        double distance(const CEssentialMatRT & existingSolution) const {
            double dist = 0;
            dist += std::min<double>((t - existingSolution.t).sum_square(), (t + existingSolution.t).sum_square()); //sign ambiguity

            /* DO NOT look at rotations: there is a rotation ambiguity
             * double qDist1 = 0, qDist2 = 0;
            for (int i = 0; i < 4; i++) {
                qDist1 += sqr(q[i] - existingSolution.q[i]);
                qDist2 += sqr(q[i] + existingSolution.q[i]);
            }
            dist += std::min<double>(qDist1, qDist2); //sign ambiguity*/
            return dist;
        }

        void makeEfromParams(const THomoImPoint * ax, const THomoImPoint * axp, const C3dRotation & q_before, TEssentialMat & E) const {
            if (IS_DEBUG) {
                //C3dRotation q_conj(-s1,-s2,-s3,1-s1);
                for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
                    C3dPoint x(ax[nPoint]); //Possibly where the float-cast is happening
                    C3dPoint xp(axp[nPoint]);
                    xp.rotate(q.conj());
                    double dResid = dotproduct(crossproduct(xp, t), x);
                    CHECK(!zero(dResid), "rot trans form of E not working");
                }

            }
            makeE_fixed(q_before, E);
        }
    };

    typedef std::vector<CEssentialMatRT> TExistingSolns;

    class CExistingSolns : public TExistingSolns {
    public:

        bool contains(const CParamVector & params, const CEForOptOnManifold & E_atManifold, const double DIST_TO_EXISTING_SOLN_THRESH, bool bVerbose) {
            if (TExistingSolns::size() == 0 || NUM_POINTS != 5)
                return false;

            CEssentialMatRT newHyp;
            if (PARAM_MODE == eMinimalOnManifold)
                newHyp = CEssentialMatRT(E_atManifold);
            else
                newHyp = CEssentialMatRT(params);

            if (bVerbose) {
                cout << "New: ";
                newHyp.print();
            }

            const bool bDebugDist = false;

            for (typename TExistingSolns::const_iterator pHyp = TExistingSolns::begin(); pHyp != TExistingSolns::end(); pHyp++) {
                if (bVerbose) pHyp->print();
                const double dDist = pHyp->distance(newHyp);
                if (bDebugDist)
                    cout << dDist << ' ';
                else {
                    if (dDist < DIST_TO_EXISTING_SOLN_THRESH)
                        return true;
                }
            }
            if (bDebugDist) cout << endl;

            return false;
        }

    };

    static void calcResiduals(const THomoImPoint * ax, const THomoImPoint * axp, const CParamVector & params, TResidVector & vResiduals, CEForOptOnManifold & E_atManifold) {
        TEssentialMat E;

        if (PARAM_MODE == eMinimal) {
            TFloat tt1, tt2, s1, s2, s3;
            params.setVars(tt1, tt2, s1, s2, s3);

            const TFloat t1 = (s2 * s1);
            const TFloat t2 = (t1 * tt1);
            const TFloat t3 = (s2 * tt2);
            const TFloat t4 = (s3 * tt1);
            const TFloat t5 = (s3 * s1);
            const TFloat t6 = (t1 * tt2);
            const TFloat t7 = (t5 * tt2);
            const TFloat t8 = (t5 * tt1);
            const TFloat t11 = 2 * s2 * tt1;
            const TFloat t12 = (s3 * s3);
            const TFloat t13 = (s1 * s1);
            const TFloat t14 = 2 * t13;
            const TFloat t15 = t12 * tt1;
            const TFloat t16 = s2 * s2;
            const TFloat t17 = t16 * tt1;
            const TFloat t18 = 2 * s1;
            const TFloat t19 = t13 * tt1;
            const TFloat t20 = 2 * t19;
            const TFloat t21 = s1 * tt1;
            const TFloat t22 = 2 * t21;
            const TFloat t23 = 2 * t2;
            const TFloat t24 = 2 * t8;
            const TFloat t25 = t11 - 1 + t12 + tt1 - t14 - t15 + t16 - t17 + t18 + t20 - t22 - t23 + t24;
            const TFloat t27 = t12 * tt2;
            const TFloat t28 = t16 * tt2;
            const TFloat t30 = 2 * s1 * tt2;
            const TFloat t32 = 2 * t13 * tt2;
            const TFloat t34 = s3 * s2;
            const TFloat t37 = 1 - t17 + t30 - t18 + t22 + t16 + t15 - t32 - t12 - tt1 - 2 * t34 * tt2;
            const TFloat t38 = t34 * tt1;
            E(0, 0) = 2 * (t1 - s3 - t2 - t3 + t4 + t5 + t6 - t7 - t8); //2 * t1 - 2 * s3 - 2 * t2 - 2 * t3 + 2 * t4 + 2 * t5 + 2 * t6 - 2 * t7 - 2 * t8;
            E(0, 1) = t25;
            E(0, 2) = 2 * t4 - t27 - t24 - t28 - t23 + tt2 - t30 + t32;
            E(1, 0) = t37;
            E(1, 1) = 2 * (-t1 + t19 - t8 + t4 + t2 - t21 - s3 + t5 + t38); //-2 * t1 + 2 * t19 - 2 * t8 + 2 * t4 + 2 * t2 - 2 * t21 - 2 * s3 + 2 * t5 + 2 * t38;
            E(1, 2) = 2 * s3 * tt2 - tt1 + t22 + t15 - t17 - 2 * t7 + 2 * t6;
            E(2, 0) = -t14 - t22 + t20 - tt2 + 2 * t34 + t18 + t28 - 2 * t38 - t27 + t30;
            E(2, 1) = t23 - t22 + 2 * s2 - 2 * t5 + t24 - 2 * t1 - t17 + tt1 + t15 - t11;
            E(2, 2) = 2 * (-t38 + t6 - t3 + t19 + t7 - t21); // -2 * t38 + 2 * t6 - 2 * t3 + 2 * t19 + 2 * t7 - 2 * t21;
        } else {
            TFloat tt1, tt2, tt3, s1, s2, s3, s4;
            if (PARAM_MODE == eFull || PARAM_MODE == eMinimalDegen)
                params.setVars(tt1, tt2, tt3, s1, s2, s3, s4);
            else if (PARAM_MODE == eMinimalOnManifold)
                E_atManifold.updateEToVars(params, tt1, tt2, tt3, s1, s2, s3, s4);
            else
                THROW("Unhandled set vars");

            const TFloat s1_sq = sqr(s1);
            const TFloat s2_sq = sqr(s2);
            const TFloat s3_sq = sqr(s3);
            const TFloat s4_sq = sqr(s4);

            E << -2 * s2 * s4 * tt2 - 2 * s3 * s1 * tt2 - 2 * s3 * s4 * tt3 + 2 * s1 * s2 * tt3, -s1_sq * tt3 - s4_sq * tt3 + s2_sq * tt3 + 2 * s3 * s1 * tt1 + 2 * s2 * s4 * tt1 + s3_sq * tt3, -2 * s1 * s2 * tt1 - s2_sq * tt2 + s1_sq * tt2 + 2 * s3 * s4 * tt1 - s3_sq * tt2 + s4_sq * tt2,
                    s2_sq * tt3 - 2 * s2 * s3 * tt2 + s4_sq * tt3 - s3_sq * tt3 - s1_sq * tt3 + 2 * s1 * s4 * tt2, -2 * s3 * s4 * tt3 + 2 * s3 * s2 * tt1 - 2 * s1 * s2 * tt3 - 2 * s1 * s4 * tt1, s3_sq * tt1 - s2_sq * tt1 + 2 * s3 * s4 * tt2 + 2 * s1 * s2 * tt2 - s4_sq * tt1 + s1_sq * tt1,
                    s2_sq * tt2 + 2 * s1 * s4 * tt3 + 2 * s3 * s2 * tt3 - s4_sq * tt2 + s1_sq * tt2 - s3_sq * tt2, 2 * s2 * s4 * tt3 + s3_sq * tt1 + s4_sq * tt1 - s1_sq * tt1 - s2_sq * tt1 - 2 * s3 * s1 * tt3, 2 * s3 * s1 * tt2 - 2 * s1 * s4 * tt1 - 2 * s3 * s2 * tt1 - 2 * s2 * s4 * tt2;
        }

        for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
            vResiduals(nPoint) = axp[nPoint].transpose() * E * ax[nPoint];
        }

        //cout << "Rnew: " << vResiduals.transpose() << endl;
        //calcResiduals_old(ax, axp, params, vResiduals);
    }

    static void calcResiduals_old(const THomoImPoint * ax, const THomoImPoint * axp, const CParamVector & params, TResidVector & vResiduals) {
        TFloat t1, t2, s1, s2, s3;
        params.setVars(t1, t2, s1, s2, s3);

        for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
            const TFloat x = ax[nPoint](0);
            const TFloat y = ax[nPoint](1);
            const TFloat xp = axp[nPoint](0);
            const TFloat yp = axp[nPoint](1);
            const TFloat zp = (NUM_POINTS == 5) ? 1 : axp[nPoint](2);
            vResiduals(nPoint) = -2 * s1 * s2 * xp * t1 + 2 * s2 * s1 * yp * t1 * y + 2 * s3 * xp * t1 + 2 * s1 * yp * t1 + 2 * s3 * yp * t1 * y + 2 * s2 * xp * t1 * y - yp * x * t1 + 2 * s1 * s2 * xp * x + 2 * s1 * s3 * xp * x - 2 * s1 * xp * t1 * y - 2 * s1 * s3 * xp * x * t1 - 2 * s2 * xp * t2 * x - 2 * s3 * s1 * yp * t1 * y + 2 * s2 * s3 * yp * t1 * y - 2 * s1 * s2 * xp * t1 * y - 2 * s1 * s1 * xp * y + 2 * s1 * s1 * xp * t2 + s3 * s3 * yp * x * t1 - s2 * s2 * yp * x * t1 - s3 * s3 * xp * t1 * y - s2 * s2 * xp * t1 * y - 2 * s1 * s3 * xp * t2 * x - s3 * s3 * yp * x - 2 * s1 * s1 * yp * t2 * x + 2 * s1 * s1 * yp * t1 * y + 2 * s1 * s1 * xp * t1 * y - 2 * s3 * yp * y + 2 * s3 * yp * t2 - 2 * s3 * s1 * yp * t2 + 2 * s1 * yp * x * t1 - 2 * s1 * s2 * xp * x * t1 + s2 * s2 * xp * y + s3 * s3 * yp * t1 - s2 * s2 * xp * t2 + s3 * s3 * xp * y + s2 * s2 * yp * x - s2 * s2 * yp * t1 - s3 * s3 * xp * t2 - 2 * s3 * s2 * yp * t2 * x - xp * y - 2 * s1 * s3 * xp * t1 - 2 * s1 * xp * t2 + 2 * s3 * s1 * zp * t2 - 2 * s3 * s1 * zp * y - 2 * s2 * zp * t1 * y + 2 * s2 * s1 * zp * t2 - 2 * s2 * s1 * zp * y - 2 * s1 * zp * x * t1 + 2 * s1 * s1 * zp * x * t1 + 2 * s3 * s2 * zp * x - 2 * s3 * s2 * zp * t1 - s2 * s2 * zp * t1 * y + s2 * s2 * zp * t2 * x + s3 * s3 * zp * t1 * y - s3 * s3 * zp * t2 * x - 2 * s1 * zp * t1 * y + 2 * s1 * zp * t2 * x + 2 * s3 * s1 * yp * y + xp * t2 + 2 * s1 * xp * y + 2 * s1 * yp * t2 * x - 2 * s2 * zp * t2 + 2 * s2 * zp * y + 2 * s1 * zp * x - 2 * s1 * zp * t1 - 2 * s1 * s1 * zp * x + 2 * s1 * s1 * zp * t1 + zp * t1 * y - zp * t2 * x - 2 * s1 * yp * x + 2 * s3 * xp * x * t1 + 2 * s2 * s1 * yp * t2 + yp * x - yp * t1 - 2 * s3 * xp * x + xp * t1 * y + 2 * s1 * s2 * xp * t2 * x - 2 * s1 * yp * t1 * y + 2 * s1 * s3 * xp * t1 * y - 2 * s2 * s1 * yp * y + 2 * s2 * s1 * zp * t1 * y - 2 * s3 * s2 * zp * x * t1 + 2 * s3 * s1 * zp * t1 * y;
        }

        cout << "R_old: " << vResiduals.transpose() << endl;
    }

    static void makeJ(const THomoImPoint * ax, const THomoImPoint * axp, const CParamVector & params, TJMatrix & J, CEForOptOnManifold & E_atManifold) {
        TEssentialMat D[NUM_PARAMS];
        if (PARAM_MODE == eMinimal) {
            TEssentialMat & Ds1 = D[2],
                    & Ds2 = D[3],
                    & Ds3 = D[4],
                    & Dtt1 = D[0],
                    & Dtt2 = D[1];

            TFloat tt1, tt2, s1, s2, s3;
            params.setVars(tt1, tt2, s1, s2, s3);

            Ds1(0, 0) = s2 - (s2 * tt1) + s3 + (s2 * tt2) - (s3 * tt2) - (s3 * tt1);
            Ds1(0, 1) = -2 * s1 + 1 + 2 * s1 * tt1 - tt1 - s2 * tt1 + s3 * tt1;
            Ds1(0, 2) = -s3 * tt1 - s2 * tt1 - tt2 + 2 * s1 * tt2;
            Ds1(1, 0) = tt2 - 1 + tt1 - 2 * s1 * tt2;
            Ds1(1, 1) = -s2 + 2 * s1 * tt1 - s3 * tt1 + s2 * tt1 - tt1 + s3;
            Ds1(1, 2) = tt1 - s3 * tt2 + s2 * tt2;
            Ds1(2, 0) = -2 * s1 - tt1 + 2 * s1 * tt1 + 1 + tt2;
            Ds1(2, 1) = s2 * tt1 - tt1 - s3 + s3 * tt1 - s2;
            Ds1(2, 2) = s2 * tt2 + 2 * s1 * tt1 + s3 * tt2 - tt1;
            Ds1 *= 2;

            if (NUM_POINTS > 3) {
                Ds2(0, 0) = s1 - s1 * tt1 - tt2 + s1 * tt2;
                Ds2(0, 1) = tt1 + s2 - s2 * tt1 - s1 * tt1;
                Ds2(0, 2) = -s2 * tt2 - s1 * tt1;
                Ds2(1, 0) = -s2 * tt1 + s2 - s3 * tt2;
                Ds2(1, 1) = -s1 + s1 * tt1 + s3 * tt1;
                Ds2(1, 2) = -s2 * tt1 + s1 * tt2;
                Ds2(2, 0) = s3 + s2 * tt2 - s3 * tt1;
                Ds2(2, 1) = s1 * tt1 + 1 - s1 - s2 * tt1 - tt1;
                Ds2(2, 2) = -s3 * tt1 + s1 * tt2 - tt2;
                Ds2 *= 2;
            }
            if (NUM_POINTS > 4) {
                Ds3(0, 0) = -1 + tt1 + s1 - s1 * tt2 - s1 * tt1;
                Ds3(0, 1) = s3 - s3 * tt1 + s1 * tt1;
                Ds3(0, 2) = tt1 - s3 * tt2 - s1 * tt1;
                Ds3(1, 0) = s3 * tt1 - s3 - s2 * tt2;
                Ds3(1, 1) = -s1 * tt1 + tt1 - 1 + s1 + s2 * tt1;
                Ds3(1, 2) = tt2 + s3 * tt1 - s1 * tt2;
                Ds3(2, 0) = -s2 * tt1 + s2 - s3 * tt2;
                Ds3(2, 1) = -s1 + s1 * tt1 + s3 * tt1;
                Ds3(2, 2) = -s2 * tt1 + s1 * tt2;
                Ds3 *= 2;
            }

            Dtt1(0, 0) = -2 * s2 * s1 + 2 * s3 - 2 * s3 * s1;
            Dtt1(0, 1) = 2 * s2 + 1 - s3 * s3 - s2 * s2 + 2 * s1 * s1 - 2 * s1 - 2 * s2 * s1 + 2 * s3 * s1;
            Dtt1(0, 2) = -2 * s2 * s1 + 2 * s3 - 2 * s3 * s1;
            Dtt1(1, 0) = -s2 * s2 + 2 * s1 + s3 * s3 - 1;
            Dtt1(1, 1) = 2 * s1 * s1 - 2 * s3 * s1 + 2 * s3 + 2 * s2 * s1 - 2 * s1 + 2 * s3 * s2;
            Dtt1(1, 2) = -s2 * s2 + 2 * s1 + s3 * s3 - 1;
            Dtt1(2, 0) = -2 * s1 + 2 * s1 * s1 - 2 * s3 * s2;
            Dtt1(2, 1) = 2 * s2 * s1 - 2 * s1 + 2 * s3 * s1 - s2 * s2 + 1 + s3 * s3 - 2 * s2;
            Dtt1(2, 2) = -2 * s1 + 2 * s1 * s1 - 2 * s3 * s2;
            Dtt2(0, 0) = -2 * s2 + 2 * s2 * s1 - 2 * s3 * s1;
            Dtt2(0, 1) = 0;
            Dtt2(0, 2) = -s3 * s3 - s2 * s2 + 1 - 2 * s1 + 2 * s1 * s1;
            Dtt2(1, 0) = 2 * s1 - 2 * s1 * s1 - 2 * s3 * s2;
            Dtt2(1, 1) = 0;
            Dtt2(1, 2) = 2 * s3 - 2 * s3 * s1 + 2 * s2 * s1;
            Dtt2(2, 0) = -1 + s2 * s2 - s3 * s3 + 2 * s1;
            Dtt2(2, 1) = 0;
            Dtt2(2, 2) = 2 * s2 * s1 - 2 * s2 + 2 * s3 * s1;
        } else if (PARAM_MODE == eFull || PARAM_MODE == eMinimalDegen) {
            TFloat tt1, tt2, tt3, s1, s2, s3, s4;
            params.setVars(tt1, tt2, tt3, s1, s2, s3, s4);

            const TFloat s1_sq = sqr(s1);
            const TFloat s2_sq = sqr(s2);
            const TFloat s3_sq = sqr(s3);
            const TFloat s4_sq = sqr(s4);

            int ntt1 = -1, ntt2 = -1, ntt3 = -1, ns1 = -1, ns2 = -1, ns3 = -1, ns4 = -1;
            if (PARAM_MODE == eFull) {
                ntt1 = 0, ntt2 = 1, ntt3 = 2, ns1 = 3, ns2 = 4, ns3 = 5, ns4 = 6;
            } else {
                ntt1 = 0, ntt2 = 1, ns1 = 2;
                if (NUM_POINTS > 3) ns2 = 3;
                if (NUM_POINTS > 3) ns3 = 4;
            }

            //tt1
            D[ntt1] << 0, 2 * s3 * s1 + 2 * s2*s4, -2 * s2 * s1 + 2 * s3*s4,
                    0, 2 * s3 * s2 - 2 * s1*s4, s3_sq - s2_sq - s4_sq + s1_sq,
                    0, s3_sq + s4_sq - s1_sq - s2_sq, -2 * s1 * s4 - 2 * s3*s2;

            //tt2
            D[ntt2] << -2 * s2 * s4 - 2 * s3*s1, 0, -s2_sq + s1_sq - s3_sq + s4_sq,
                    -2 * s3 * s2 + 2 * s1*s4, 0, 2 * s3 * s4 + 2 * s2*s1,
                    s2_sq - s4_sq + s1_sq - s3_sq, 0, 2 * s3 * s1 - 2 * s2*s4;

            //tt3
            if (ntt3 >= 0) D[ntt3] << -2 * s3 * s4 + 2 * s2 * s1, -s1_sq - s4_sq + s2_sq + s3_sq, 0,
                    s2_sq + s4_sq - s3_sq - s1_sq, -2 * s3 * s4 - 2 * s2 * s1, 0,
                    2 * s1 * s4 + 2 * s3 * s2, 2 * s2 * s4 - 2 * s3 * s1, 0;

            //s1
            D[ns1] << -2 * s3 * tt2 + 2 * s2*tt3, -2 * s1 * tt3 + 2 * s3*tt1, -2 * s2 * tt1 + 2 * s1*tt2,
                    -2 * s1 * tt3 + 2 * s4*tt2, -2 * s2 * tt3 - 2 * s4*tt1, 2 * s2 * tt2 + 2 * s1*tt1,
                    2 * s4 * tt3 + 2 * s1*tt2, -2 * s1 * tt1 - 2 * s3*tt3, 2 * s3 * tt2 - 2 * s4*tt1;

            //s2
            if (ns2 >= 0) D[ns2] << 2 * s1 * tt3 - 2 * s4 * tt2, 2 * s2 * tt3 + 2 * s4 * tt1, -2 * s1 * tt1 - 2 * s2 * tt2,
                    -2 * s3 * tt2 + 2 * s2 * tt3, -2 * s1 * tt3 + 2 * s3 * tt1, -2 * s2 * tt1 + 2 * s1 * tt2,
                    2 * s2 * tt2 + 2 * s3 * tt3, 2 * s4 * tt3 - 2 * s2 * tt1, -2 * s3 * tt1 - 2 * s4 * tt2;

            //s3    ****GCC's warnings here can be ignored****
            if (ns3 >= 0) D[ns3] << -2 * s1 * tt2 - 2 * s4 * tt3, 2 * s1 * tt1 + 2 * s3 * tt3, 2 * s4 * tt1 - 2 * s3 * tt2,
                    -2 * s2 * tt2 - 2 * s3 * tt3, -2 * s4 * tt3 + 2 * s2 * tt1, 2 * s3 * tt1 + 2 * s4 * tt2,
                    -2 * s3 * tt2 + 2 * s2 * tt3, -2 * s1 * tt3 + 2 * s3 * tt1, -2 * s2 * tt1 + 2 * s1 * tt2;

            //s4
            if (ns4 >= 0) D[ns4] << -2 * s2 * tt2 - 2 * s3 * tt3, -2 * s4 * tt3 + 2 * s2 * tt1, 2 * s3 * tt1 + 2 * s4 * tt2,
                    2 * s4 * tt3 + 2 * s1 * tt2, -2 * s1 * tt1 - 2 * s3 * tt3, 2 * s3 * tt2 - 2 * s4 * tt1,
                    2 * s1 * tt3 - 2 * s4 * tt2, 2 * s2 * tt3 + 2 * s4 * tt1, -2 * s1 * tt1 - 2 * s2 * tt2;
        } else if (PARAM_MODE == eMinimalOnManifold) {
            TFloat tt1, tt2, tt3, s1, s2, s3, s4;
            E_atManifold.setVars(tt1, tt2, tt3, s1, s2, s3, s4); //independent of params
            const TFloat s1_sq = sqr(s1);
            const TFloat s2_sq = sqr(s2);
            const TFloat s3_sq = sqr(s3);
            const TFloat s4_sq = sqr(s4);

            //1st 2 params are for tangent space of t
            for (int nt_basis_vec = 0; nt_basis_vec < 2; nt_basis_vec++) {
                TFloat tau1, tau2, tau3;
                E_atManifold.setTau(nt_basis_vec, tau1, tau2, tau3);

                D[nt_basis_vec] << 2 * s1 * s2 * tau3 - 2 * s2 * s4 * tau2 - 2 * s3 * s1 * tau2 - 2 * s3 * s4*tau3, 2 * s3 * s1 * tau1 + 2 * s2 * s4 * tau1 + s2_sq * tau3 + s3_sq * tau3 - s1_sq * tau3 - s4_sq*tau3, -s2_sq * tau2 + s1_sq * tau2 + s4_sq * tau2 + 2 * s3 * s4 * tau1 - 2 * s1 * s2 * tau1 - s3_sq*tau2,
                        2 * s1 * s4 * tau2 - 2 * s2 * s3 * tau2 - s3_sq * tau3 - s1_sq * tau3 + s2_sq * tau3 + s4_sq*tau3, -2 * s1 * s4 * tau1 - 2 * s3 * s4 * tau3 - 2 * s1 * s2 * tau3 + 2 * s2 * s3*tau1, -s4_sq * tau1 + s3_sq * tau1 - s2_sq * tau1 + 2 * s1 * s2 * tau2 + s1_sq * tau1 + 2 * s3 * s4*tau2,
                        2 * s1 * s4 * tau3 + s2_sq * tau2 - s4_sq * tau2 + s1_sq * tau2 - s3_sq * tau2 + 2 * s2 * s3*tau3, 2 * s2 * s4 * tau3 + s3_sq * tau1 - 2 * s3 * s1 * tau3 + s4_sq * tau1 - s2_sq * tau1 - s1_sq*tau1, 2 * s3 * s1 * tau2 - 2 * s2 * s4 * tau2 - 2 * s2 * s3 * tau1 - 2 * s1 * s4*tau1;
            }

            //Next 3 params are for tangent space of s
            for (int nq_basis_vec = 0; nq_basis_vec < 3; nq_basis_vec++) {
                TFloat chi1, chi2, chi3, chi4;
                E_atManifold.setChi(nq_basis_vec, chi1, chi2, chi3, chi4);
                D[2 + nq_basis_vec] << 2 * tt3 * chi1 * s2 - 2 * tt2 * s2 * chi4 - 2 * tt2 * chi2 * s4 + 2 * tt3 * s1 * chi2 - 2 * tt2 * chi3 * s1 - 2 * tt3 * s3 * chi4 - 2 * tt3 * chi3 * s4 - 2 * tt2 * s3*chi1, 2 * tt3 * s2 * chi2 - 2 * tt3 * s1 * chi1 - 2 * tt3 * s4 * chi4 + 2 * tt1 * chi2 * s4 + 2 * tt3 * s3 * chi3 + 2 * tt1 * s3 * chi1 + 2 * tt1 * chi3 * s1 + 2 * tt1 * s2*chi4, 2 * tt2 * s1 * chi1 + 2 * tt2 * s4 * chi4 + 2 * tt1 * s3 * chi4 - 2 * tt2 * s2 * chi2 - 2 * tt1 * s1 * chi2 - 2 * tt1 * chi1 * s2 - 2 * tt2 * s3 * chi3 + 2 * tt1 * chi3*s4,
                        2 * tt3 * s2 * chi2 + 2 * tt3 * s4 * chi4 - 2 * tt3 * s1 * chi1 + 2 * tt2 * s1 * chi4 + 2 * tt2 * chi1 * s4 - 2 * tt2 * s2 * chi3 - 2 * tt2 * chi2 * s3 - 2 * tt3 * s3*chi3, 2 * tt1 * s2 * chi3 - 2 * tt3 * s3 * chi4 - 2 * tt3 * chi3 * s4 - 2 * tt3 * s1 * chi2 - 2 * tt1 * chi1 * s4 + 2 * tt1 * chi2 * s3 - 2 * tt3 * chi1 * s2 - 2 * tt1 * s1*chi4, -2 * tt1 * s2 * chi2 + 2 * tt2 * s1 * chi2 + 2 * tt2 * chi1 * s2 - 2 * tt1 * s4 * chi4 + 2 * tt1 * s3 * chi3 + 2 * tt2 * s3 * chi4 + 2 * tt2 * chi3 * s4 + 2 * tt1 * s1*chi1,
                        2 * tt3 * s1 * chi4 + 2 * tt3 * chi1 * s4 + 2 * tt2 * s2 * chi2 - 2 * tt2 * s4 * chi4 + 2 * tt2 * s1 * chi1 - 2 * tt2 * s3 * chi3 + 2 * tt3 * s2 * chi3 + 2 * tt3 * chi2*s3, 2 * tt3 * s2 * chi4 + 2 * tt3 * chi2 * s4 - 2 * tt1 * s2 * chi2 - 2 * tt3 * s3 * chi1 - 2 * tt3 * chi3 * s1 + 2 * tt1 * s4 * chi4 - 2 * tt1 * s1 * chi1 + 2 * tt1 * s3*chi3, 2 * tt2 * s3 * chi1 + 2 * tt2 * chi3 * s1 - 2 * tt2 * s2 * chi4 - 2 * tt2 * chi2 * s4 - 2 * tt1 * s2 * chi3 - 2 * tt1 * chi2 * s3 - 2 * tt1 * s1 * chi4 - 2 * tt1 * chi1*s4;

                //D[2+nq_basis_vec] << 2*chi1*(s2+b*chi2)*tt3+(2*(s1+b*chi1))*chi2*tt3-2*chi2*(s4+b*chi4)*tt2-(2*(s2+b*chi2))*chi4*tt2-2*chi3*(s1+b*chi1)*tt2-(2*(s3+b*chi3))*chi1*tt2-2*chi3*(s4+b*chi4)*tt3-(2*(s3+b*chi3))*chi4*tt3,     2*chi3*(s1+b*chi1)*tt1+(2*(s3+b*chi3))*chi1*tt1+2*chi2*(s4+b*chi4)*tt1+(2*(s2+b*chi2))*chi4*tt1+(2*(s2+b*chi2))*tt3*chi2+(2*(s3+b*chi3))*tt3*chi3-(2*(s1+b*chi1))*tt3*chi1-(2*(s4+b*chi4))*tt3*chi4, -(2*(s2+b*chi2))*tt2*chi2+(2*(s1+b*chi1))*tt2*chi1+(2*(s4+b*chi4))*tt2*chi4+2*chi3*(s4+b*chi4)*tt1+(2*(s3+b*chi3))*chi4*tt1-2*chi1*(s2+b*chi2)*tt1-(2*(s1+b*chi1))*chi2*tt1-(2*(s3+b*chi3))*tt2*chi3,
                //					 2*chi1*(s4+b*chi4)*tt2+(2*(s1+b*chi1))*chi4*tt2-2*chi2*(s3+b*chi3)*tt2-(2*(s2+b*chi2))*chi3*tt2-(2*(s3+b*chi3))*tt3*chi3-(2*(s1+b*chi1))*tt3*chi1+(2*(s2+b*chi2))*tt3*chi2+(2*(s4+b*chi4))*tt3*chi4, -2*chi1*(s4+b*chi4)*tt1-(2*(s1+b*chi1))*chi4*tt1-2*chi3*(s4+b*chi4)*tt3-(2*(s3+b*chi3))*chi4*tt3-2*chi1*(s2+b*chi2)*tt3-(2*(s1+b*chi1))*chi2*tt3+2*chi2*(s3+b*chi3)*tt1+(2*(s2+b*chi2))*chi3*tt1,    -(2*(s4+b*chi4))*tt1*chi4+(2*(s3+b*chi3))*tt1*chi3-(2*(s2+b*chi2))*tt1*chi2+2*chi1*(s2+b*chi2)*tt2+(2*(s1+b*chi1))*chi2*tt2+(2*(s1+b*chi1))*tt1*chi1+2*chi3*(s4+b*chi4)*tt2+(2*(s3+b*chi3))*chi4*tt2,
                //					 2*chi1*(s4+b*chi4)*tt3+(2*(s1+b*chi1))*chi4*tt3+(2*(s2+b*chi2))*tt2*chi2-(2*(s4+b*chi4))*tt2*chi4+(2*(s1+b*chi1))*tt2*chi1-(2*(s3+b*chi3))*tt2*chi3+2*chi2*(s3+b*chi3)*tt3+(2*(s2+b*chi2))*chi3*tt3, 2*chi2*(s4+b*chi4)*tt3+(2*(s2+b*chi2))*chi4*tt3+(2*(s3+b*chi3))*tt1*chi3-2*chi3*(s1+b*chi1)*tt3-(2*(s3+b*chi3))*chi1*tt3+(2*(s4+b*chi4))*tt1*chi4-(2*(s2+b*chi2))*tt1*chi2-(2*(s1+b*chi1))*tt1*chi1, 2*chi3*(s1+b*chi1)*tt2+(2*(s3+b*chi3))*chi1*tt2-2*chi2*(s4+b*chi4)*tt2-(2*(s2+b*chi2))*chi4*tt2-2*chi2*(s3+b*chi3)*tt1-(2*(s2+b*chi2))*chi3*tt1-2*chi1*(s4+b*chi4)*tt1-(2*(s1+b*chi1))*chi4*tt1;
            }
        } else
            THROW("Unhandled PARAM_MODE");

        for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
            for (int nParam = 0; nParam < NUM_PARAMS; nParam++)
                J(nPoint, nParam) = axp[nPoint].transpose() * D[nParam] * ax[nPoint];
        }
        //cout << J << "=J_new" << endl;
        //makeJ_old(ax, axp, params, J);
    }

    static void makeJ_old(const THomoImPoint * ax, const THomoImPoint * axp, const CParamVector & params, TJMatrix & J) {
        TFloat t1, t2, s1, s2, s3;
        params.setVars(t1, t2, s1, s2, s3);

        for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
            const TFloat x = ax[nPoint](0);
            const TFloat y = ax[nPoint](1);
            const TFloat xp = axp[nPoint](0);
            const TFloat yp = axp[nPoint](1);
            const TFloat zp = (NUM_POINTS == 5) ? 1 : axp[nPoint](2);

            J(nPoint, 0) = -2 * s2 * xp * s1 + 2 * s2 * s3 * yp * y - 2 * s1 * s2 * xp * y + 2 * s1 * s3 * xp * y - 2 * s1 * s2 * xp * x - 2 * s1 * s3 * xp * x + 2 * s1 * s1 * xp * y + s3 * s3 * yp * x + 2 * s3 * yp * y - yp + 2 * s3 * xp + 2 * s1 * yp - s2 * s2 * xp * y - s3 * s3 * xp * y - s2 * s2 * yp * x + xp * y + 2 * s3 * s1 * zp * y + 2 * s2 * s1 * zp * y - 2 * s3 * s2 * zp * x - 2 * s3 * s2 * zp - s2 * s2 * zp * y + s3 * s3 * zp * y - 2 * s1 * zp * y - 2 * s3 * s1 * yp * y - 2 * s1 * xp * y - 2 * s1 * yp * y + 2 * s2 * xp * y + 2 * s1 * s1 * yp * y - 2 * s2 * zp * y - 2 * s1 * zp * x + 2 * s1 * s1 * zp * x + 2 * s1 * yp * x + 2 * s1 * s1 * zp + zp * y - yp * x + 2 * s3 * xp * x - 2 * s1 * s3 * xp + s3 * s3 * yp - yp * s2 * s2 - 2 * s1 * zp + 2 * s2 * s1 * yp * y;
            J(nPoint, 1) = -2 * s2 * xp * x + 2 * xp * s1 * s1 - 2 * s1 * s3 * xp * x - 2 * s1 * s1 * yp * x + 2 * s3 * yp - 2 * s1 * s3 * yp - s2 * s2 * xp - s3 * s3 * xp - 2 * s3 * s2 * yp * x - 2 * xp * s1 + 2 * s3 * s1 * zp + 2 * s2 * s1 * zp + s2 * s2 * zp * x - s3 * s3 * zp * x + 2 * s1 * zp * x + xp + 2 * s1 * yp * x - 2 * s2 * zp - zp * x + 2 * s1 * yp * s2 + 2 * s1 * s2 * xp * x;
            J(nPoint, 2) = -2 * s3 * xp * t1 - 2 * s2 * xp * t1 - 2 * s3 * yp * t1 * y - 2 * s2 * xp * t1 * y + 2 * s2 * xp * x + 2 * yp * x * t1 - 2 * s2 * yp * y + 4 * s1 * xp * t1 * y + 2 * s2 * xp * t2 * x + 2 * s2 * yp * t1 * y - 2 * s3 * xp * t2 * x - 2 * s2 * xp * x * t1 + 2 * s3 * xp * t1 * y + 2 * yp * t2 * x + 2 * s2 * yp * t2 - 2 * yp * t1 * y + 2 * s3 * yp * y - 2 * s3 * yp * t2 + 2 * s3 * zp * t2 - 2 * s3 * zp * y - 2 * zp * x * t1 + 2 * xp * y + 4 * s1 * xp * t2 + 2 * s2 * zp * t1 * y + 4 * s1 * zp * x * t1 - 2 * xp * t2 - 4 * s1 * xp * y - 4 * s1 * yp * t2 * x + 2 * s2 * zp * t2 - 2 * s2 * zp * y - 4 * s1 * zp * x + 4 * s1 * zp * t1 - 2 * zp * t1 * y + 2 * zp * t2 * x + 2 * s3 * zp * t1 * y - 2 * s3 * xp * x * t1 + 2 * zp * x - 2 * zp * t1 - 2 * yp * x + 2 * yp * t1 + 2 * s3 * xp * x - 2 * xp * t1 * y + 4 * s1 * yp * t1 * y;
            if (NUM_POINTS > 3) J(nPoint, 3) = -2 * s1 * xp * t1 + 2 * s3 * yp * t1 * y - 2 * s2 * xp * t1 * y - 2 * s1 * xp * t1 * y + 2 * s2 * yp * x - 2 * s2 * yp * t1 - 2 * s1 * xp * x * t1 - 2 * s3 * yp * t2 * x + 2 * s1 * xp * t2 * x - 2 * s2 * xp * t2 + 2 * s1 * xp * x - 2 * xp * t2 * x - 2 * s2 * yp * x * t1 + 2 * s1 * zp * t2 + 2 * s3 * zp * x - 2 * s3 * zp * t1 + 2 * s2 * zp * t2 * x + 2 * s1 * yp * t2 - 2 * s2 * zp * t1 * y + 2 * s1 * zp * t1 * y - 2 * s1 * zp * y - 2 * s1 * yp * y + 2 * s2 * xp * y - 2 * zp * t1 * y - 2 * s3 * zp * x * t1 + 2 * zp * y - 2 * zp * t2 + 2 * xp * t1 * y + 2 * s1 * yp * t1 * y;
            if (NUM_POINTS > 4) J(nPoint, 4) = -2 * s1 * xp * t1 + 2 * s3 * yp * t1 + 2 * s1 * xp * t1 * y - 2 * s3 * yp * x - 2 * yp * y + 2 * yp * t2 + 2 * s3 * xp * y - 2 * s3 * xp * t2 - 2 * xp * x + 2 * xp * t1 + 2 * s2 * yp * t1 * y - 2 * s3 * xp * t1 * y - 2 * s1 * xp * x * t1 - 2 * s1 * xp * t2 * x - 2 * s2 * yp * t2 * x + 2 * yp * t1 * y + 2 * s1 * xp * x + 2 * s3 * yp * x * t1 + 2 * xp * x * t1 + 2 * s1 * zp * t2 + 2 * s2 * zp * x - 2 * s2 * zp * t1 - 2 * s3 * zp * t2 * x - 2 * s1 * yp * t2 + 2 * s1 * zp * t1 * y - 2 * s1 * zp * y + 2 * s1 * yp * y + 2 * s3 * zp * t1 * y - 2 * s2 * zp * x * t1 - 2 * s1 * yp * t1 * y;
        }
        cout << J << "=J_old" << endl;
    }

    //Hessian of the residuals, scaled by r, not of the sum-of-squared residuals
    //See derivation in http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm

    static void make_rH(const THomoImPoint * ax, const THomoImPoint * axp, const CParamVector & params, const TResidVector & vResiduals, TJTJMatrix & H) {
        TFloat t1, t2, s1, s2, s3;
        params.setVars(t1, t2, s1, s2, s3);

        H.setZero();

        for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
            const TFloat x = ax[nPoint](0);
            const TFloat y = ax[nPoint](1);
            const TFloat xp = axp[nPoint](0);
            const TFloat yp = axp[nPoint](1);
            const TFloat zp = (NUM_POINTS == 5) ? 1 : axp[nPoint](2);
            const TFloat ri = vResiduals(nPoint);

            H(0, 2) += ri * (-2 * s2 * xp - 2 * s2 * xp * y + 2 * s3 * xp * y + 2 * s2 * zp * y - 2 * s2 * xp * x - 2 * s3 * xp * x + 2 * s3 * zp * y - 2 * zp * y + 4 * s1 * xp * y - 2 * zp * x + 4 * s1 * zp * x + 4 * s1 * zp + 2 * yp - 2 * s3 * yp * y - 2 * xp * y - 2 * yp * y + 4 * s1 * yp * y + 2 * yp * x - 2 * s3 * xp + 2 * s2 * yp * y - 2 * zp);
            H(1, 2) += ri * (2 * s2 * zp + 2 * s3 * zp + 2 * zp * x + 4 * xp * s1 - 2 * s3 * xp * x - 4 * s1 * yp * x - 2 * s3 * yp - 2 * xp + 2 * yp * x + 2 * yp * s2 + 2 * s2 * xp * x);
            H(2, 2) += ri * (4 * zp * x * t1 + 4 * xp * t1 * y - 4 * zp * x + 4 * zp * t1 + 4 * xp * t2 - 4 * xp * y - 4 * yp * t2 * x + 4 * yp * t1 * y);
            if (NUM_POINTS > 3) {
                H(0, 3) += ri * (-2 * xp * s1 + 2 * s3 * yp * y - 2 * s1 * xp * y + 2 * s1 * zp * y - 2 * s1 * xp * x - 2 * s3 * zp * x - 2 * s3 * zp - 2 * s2 * zp * y - 2 * zp * y - 2 * s2 * xp * y - 2 * s2 * yp * x + 2 * xp * y - 2 * yp * s2 + 2 * s1 * yp * y);
                H(1, 3) += ri * (2 * s1 * zp + 2 * s2 * zp * x - 2 * xp * x - 2 * zp - 2 * s2 * xp - 2 * s3 * yp * x + 2 * s1 * yp + 2 * s1 * xp * x);
                H(2, 3) += ri * (-2 * xp * t1 - 2 * xp * t1 * y + 2 * xp * x - 2 * yp * y + 2 * zp * t1 * y + 2 * xp * t2 * x + 2 * yp * t1 * y - 2 * xp * x * t1 + 2 * yp * t2 + 2 * zp * t2 - 2 * zp * y);
                H(3, 3) += ri * (-2 * xp * t1 * y - 2 * zp * t1 * y + 2 * yp * x - 2 * yp * t1 - 2 * xp * t2 - 2 * yp * x * t1 + 2 * zp * t2 * x + 2 * xp * y);
            }
            if (NUM_POINTS > 4) {
                H(0, 4) += ri * (2 * s2 * yp * y + 2 * s1 * xp * y - 2 * s1 * xp * x + 2 * s1 * zp * y - 2 * s2 * zp * x - 2 * s2 * zp + 2 * s3 * zp * y + 2 * s3 * yp * x + 2 * yp * y + 2 * xp - 2 * s3 * xp * y - 2 * s1 * yp * y + 2 * xp * x - 2 * xp * s1 + 2 * s3 * yp);
                H(1, 4) += ri * (2 * s1 * zp - 2 * s3 * zp * x - 2 * s1 * xp * x + 2 * yp - 2 * s1 * yp - 2 * s3 * xp - 2 * s2 * yp * x);
                H(2, 4) += ri * (-2 * xp * t1 + 2 * zp * t1 * y - 2 * yp * t1 * y - 2 * xp * t2 * x + 2 * xp * t1 * y + 2 * yp * y - 2 * yp * t2 + 2 * zp * t2 - 2 * zp * y - 2 * xp * x * t1 + 2 * xp * x);
                H(3, 4) += ri * (-2 * zp * x * t1 + 2 * yp * t1 * y - 2 * yp * t2 * x + 2 * zp * x - 2 * zp * t1);
                H(4, 4) += ri * (2 * yp * t1 + 2 * zp * t1 * y - 2 * yp * x + 2 * xp * y - 2 * xp * t2 - 2 * xp * t1 * y + 2 * yp * x * t1 - 2 * zp * t2 * x);
            }
            //H(3, 0) += ri * (-2 * xp * s1 + 2 * s3 * yp * y - 2 * s1 * xp * y + 2 * s1 * zp * y - 2 * s1 * xp * x - 2 * s3 * zp * x - 2 * s3 * zp - 2 * s2 * zp * y - 2 * zp * y - 2 * s2 * xp * y - 2 * s2 * yp * x + 2 * xp * y - 2 * yp * s2 + 2 * s1 * yp * y);
            //H(3, 1) += ri * (2 * s1 * zp + 2 * s2 * zp * x - 2 * xp * x - 2 * zp - 2 * s2 * xp - 2 * s3 * yp * x + 2 * s1 * yp + 2 * s1 * xp * x);
            //H(3, 2) += ri * (-2 * xp * t1 - 2 * xp * t1 * y + 2 * xp * x - 2 * yp * y + 2 * zp * t1 * y + 2 * xp * t2 * x + 2 * yp * t1 * y - 2 * xp * x * t1 + 2 * yp * t2 + 2 * zp * t2 - 2 * zp * y);
            //H(2, 0) += ri * (-2 * s2 * xp - 2 * s2 * xp * y + 2 * s3 * xp * y + 2 * s2 * zp * y - 2 * s2 * xp * x - 2 * s3 * xp * x + 2 * s3 * zp * y - 2 * zp * y + 4 * s1 * xp * y - 2 * zp * x + 4 * s1 * zp * x + 4 * s1 * zp + 2 * yp - 2 * s3 * yp * y - 2 * xp * y - 2 * yp * y + 4 * s1 * yp * y + 2 * yp * x - 2 * s3 * xp + 2 * s2 * yp * y - 2 * zp);
            //H(2, 1) += ri * (2 * s2 * zp + 2 * s3 * zp + 2 * zp * x + 4 * xp * s1 - 2 * s3 * xp * x - 4 * s1 * yp * x - 2 * s3 * yp - 2 * xp + 2 * yp * x + 2 * yp * s2 + 2 * s2 * xp * x);
            //H(4, 0) += ri * (2 * s2 * yp * y + 2 * s1 * xp * y - 2 * s1 * xp * x + 2 * s1 * zp * y - 2 * s2 * zp * x - 2 * s2 * zp + 2 * s3 * zp * y + 2 * s3 * yp * x + 2 * yp * y + 2 * xp - 2 * s3 * xp * y - 2 * s1 * yp * y + 2 * xp * x - 2 * xp * s1 + 2 * s3 * yp);
            //H(4, 1) += ri * (2 * s1 * zp - 2 * s3 * zp * x - 2 * s1 * xp * x + 2 * yp - 2 * s1 * yp - 2 * s3 * xp - 2 * s2 * yp * x);
            //H(4, 2) += ri * (-2 * xp * t1 + 2 * zp * t1 * y - 2 * yp * t1 * y - 2 * xp * t2 * x + 2 * xp * t1 * y + 2 * yp * y - 2 * yp * t2 + 2 * zp * t2 - 2 * zp * y - 2 * xp * x * t1 + 2 * xp * x);
            //H(4, 3) += ri * (-2 * zp * x * t1 + 2 * yp * t1 * y - 2 * yp * t2 * x + 2 * zp * x - 2 * zp * t1);
        }
        for (int r = 1; r < NUM_POINTS; r++)
            for (int c = 0; c < r; c++)
                H(r, c) = H(c, r);
    }

    static bool getHypothesis(const CEssentialMatGradientDescParams & PARAMS, const THomoImPoint * ax, const THomoImPoint * axp, TEssentialMat & E, CExistingSolns & aExistingSolns, CGDStats & stats) {
        CParamVector params;
        C3dRotation q_before;

        ARRAY(THomoImPoint, axp_localCopy, NUM_POINTS);
        for (int i = 0; i < NUM_POINTS; i++)
            axp_localCopy[i] = axp[i];

        if (NUM_POINTS < 5) //Apply a random rotation to xp;
        {
            q_before = C3dRotation(FastNormal(), FastNormal(), FastNormal(), FastNormal());
            C3dRotation q_conj(-q_before[0], -q_before[1], -q_before[2], q_before[3]);

            for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
                C3dPoint xp(axp[nPoint]);
                xp.rotate(q_conj);
                axp_localCopy[nPoint] = THomoImPoint((TFloat) xp.getX(), (TFloat) xp.getY(), (TFloat) xp.getZ());
            }
        }

        params.randomise();
        CEForOptOnManifold E_atManifold(params);

        if (getHypothesisGradientDesc(PARAMS, ax, PTR(axp_localCopy), params, aExistingSolns, stats, E_atManifold)) {
            if (PARAM_MODE != eMinimalOnManifold) {
                CEssentialMatRT E_rt(params);
                E_rt.makeEfromParams(ax, PTR(axp_localCopy), q_before, E);
                aExistingSolns.push_back(E_rt);
            } else {
                CEssentialMatRT E_rt(E_atManifold);
                E_rt.makeEfromParams(ax, PTR(axp_localCopy), q_before, E);
                aExistingSolns.push_back(E_rt);
            }
            return true;
        } else
            return false;

    }

    class C3ptBracket {
        //CParamClass aParams[3];
        //TResidVector aResids[3];
        TFloat adAlpha[3], adErr[3];
        int nActive, nSet;

        //DO NOT NEED for quad. fitting anyway
        //CParamClass & activeParam() { return aParams[nActive]; }

        TFloat & activealpha() {
            return adAlpha[nActive];
        }
        //DO NOT NEED for quad. fitting anyway
        //TResidVector & activeResid() { return aResids[nActive]; }

        TFloat & activeErr() {
            return adErr[nActive];
        }

        void copy(const int nFrom, const int nTo) {
            adAlpha[nTo] = adAlpha[nFrom];
            adErr[nTo] = adErr[nFrom];
        }
    public:

        C3ptBracket() : nActive(0), nSet(0) {
        }

        //When we know the residuals at the start point:

        C3ptBracket(const CParamVector & params, TResidVector & vResiduals) : nActive(0), nSet(1) {
            //activeParam() = params;
            //activeResid() = vResiduals;
            activeErr() = vResiduals.squaredNorm();
            activealpha() = 0;
            nActive = 1;
        }

        enum eDirection {
            eAdvance, eRetreat, eDone
        };

        eDirection addResiduals(const CParamVector & params, const TFloat alpha, const THomoImPoint * ax, const THomoImPoint * axp, CEForOptOnManifold & E_atManifold) {
            nSet++;
            //activeParam() = params;
            TResidVector residuals;
            calcResiduals(ax, axp, params, residuals, E_atManifold);
            activeErr() = residuals.squaredNorm();
            activealpha() = alpha;

            if (nSet == 1) {
                nActive = 1;
                return eAdvance;
            } else if (nSet == 2) {
                if (activeErr() < adErr[0]) {
                    //Advance until an upturn. Already set 0...1
                    nActive = 2;
                    return eAdvance;
                } else {
                    //Retreat until <= initial error
                    //Copy 1..2
                    copy(1, 2);
                    //nActive == 1 now
                    return eRetreat;
                }
            } else //We either have a 3-bracket, are strictly increasing, or are strictly decreasing.
            {
                if (nActive == 1) //we're retreating until < r(0)
                {
                    if (activeErr() < adErr[0])
                        return eDone;
                    else {
                        copy(1, 2);
                        return eRetreat;
                    }
                } else //We're increasing alpha until error starts increasing again
                {
                    if (activeErr() > adErr[1])
                        return eDone;
                    else {
                        copy(1, 0);
                        copy(2, 1);
                        return eAdvance;
                    }
                }
            }
        }

        TFloat findMinAlphaByQuadratic() {
            //Method by http://www-personal.umich.edu/~murty/611/611slides9.pdf
            if(IS_DEBUG) CHECK(adErr[0] < adErr[1] || adErr[2] < adErr[1], "Not a 3-bracket");

            const TFloat alpha1 = adAlpha[0]; //notation to match paper
            const TFloat alpha2 = adAlpha[1];
            const TFloat alpha3 = adAlpha[2];
            const TFloat f_alpha1 = adErr[0];
            const TFloat f_alpha2 = adErr[1];
            const TFloat f_alpha3 = adErr[2];
            const TFloat alpha1_sq = sqr(alpha1);
            const TFloat alpha2_sq = sqr(alpha2);
            const TFloat alpha3_sq = sqr(alpha3);

            TFloat alphaMin = (alpha2_sq - alpha3_sq) * f_alpha1 + (alpha3_sq - alpha1_sq) * f_alpha2 + (alpha1_sq - alpha2_sq) * f_alpha3;
            alphaMin /= 2 * ((alpha2 - alpha3) * f_alpha1 + (alpha3 - alpha1) * f_alpha2 + (alpha1 - alpha2) * f_alpha3);

            if(IS_DEBUG) CHECK(alphaMin < alpha1 || alphaMin > alpha3, "alpha min failed");
            return alphaMin;
        }
    };

    //Method by http://www-personal.umich.edu/~murty/611/611slides9.pdf
    //residuals and params should be full on way in and out

    static TFloat find3ptBracket(const CEssentialMatGradientDescParams::CSDParams & PARAMS, const THomoImPoint * ax, const THomoImPoint * axp, const TJMatrix & J, const TFloat dErr, const CParamVector & delX, CParamVector & params, TResidVector & vResiduals, CEForOptOnManifold & E_atManifold) {
        C3ptBracket bracket(params, vResiduals);

        TFloat alphaScale;
        //TFloat alpha = PARAMS.ALPHA
        TFloat alpha = 2 * vResiduals.squaredNorm() / sqrt(delX.squaredNorm());

        CParamVector newParams;
        newParams = params + alpha*delX;
        newParams.Normalise();

        if (bracket.addResiduals(newParams, alpha, ax, axp, E_atManifold) == C3ptBracket::eAdvance)
            alphaScale = 2;
        else
            alphaScale = 0.5;

        do {
            alpha *= alphaScale;
            CParamVector alpha_delX;
            alpha_delX = alpha*delX;
            newParams = params + alpha_delX;
            newParams.Normalise();
            //if(alpha_delX.squaredNorm()<PARAMS.SMALL_STEP_DROPOUT) //2e-14, 11.4 secs 
            if (alpha < PARAMS.SMALL_STEP_DROPOUT) //2e-8, 10 secs 
                return vResiduals.squaredNorm(); //Force false minima
        } while (bracket.addResiduals(newParams, alpha, ax, axp, E_atManifold) != C3ptBracket::eDone);

        TFloat alphaMin = bracket.findMinAlphaByQuadratic();
        params = params + alphaMin*delX;
        params.Normalise();

        //cout << (delX/sqrt(delX.squaredNorm())).transpose() << endl;

        calcResiduals(ax, axp, params, vResiduals, E_atManifold);

        return vResiduals.squaredNorm();
    }

    //Set better params and return new error. residuals should be full

    static TFloat minSteepestDescent(const CEssentialMatGradientDescParams::CSDParams & PARAMS, const THomoImPoint * ax, const THomoImPoint * axp, const TJMatrix & J, const TFloat dErr, CParamVector & params, TResidVector & vResiduals, CEForOptOnManifold & E_atManifold) {
        CParamVector delX;
        delX = -J.transpose() * vResiduals;
        return find3ptBracket(PARAMS, ax, axp, J, dErr, delX, params, vResiduals, E_atManifold);
    }

    //Set better params and return new error. residuals should be full

    static TFloat minCGD(const CEssentialMatGradientDescParams::CSDParams & PARAMS, const THomoImPoint * ax, const THomoImPoint * axp, const TJMatrix & J, const TFloat dErr, CParamVector & params, TResidVector & vResiduals, CParamVector & LAMBDAx_lastIter, CParamVector & oldDelX, CEForOptOnManifold & E_atManifold) {
        CParamVector delX;
        delX = -J.transpose() * vResiduals;

        TFloat beta = 0;
        if (PARAMS.BETA_FN == CEssentialMatGradientDescParams::CSDParams::eFletcherReeves)
            beta = delX.squaredNorm() / oldDelX.squaredNorm();
        else if (PARAMS.BETA_FN == CEssentialMatGradientDescParams::CSDParams::ePolakRibiere) {
            beta = (delX.transpose()*(delX - oldDelX));
            beta /= oldDelX.squaredNorm();
        } else if (PARAMS.BETA_FN == CEssentialMatGradientDescParams::CSDParams::eHestenesStiefel) {
            beta = (delX.transpose()*(delX - oldDelX));
            beta /= (LAMBDAx_lastIter.transpose()*(delX - oldDelX));
        } else THROW("BETA_FN unhandled");

        if (beta < 0) beta = 0;

        LAMBDAx_lastIter = delX + beta*LAMBDAx_lastIter;

        oldDelX = delX;
        return find3ptBracket(PARAMS, ax, axp, J, dErr, LAMBDAx_lastIter, params, vResiduals, E_atManifold);
    }

    //Set better params and return new error. residuals should be full

    static TFloat minLevMar(const CEssentialMatGradientDescParams::CLMParams & PARAMS, const THomoImPoint * ax, const THomoImPoint * axp, const TJMatrix & J, const TFloat dErr, CParamVector & params, TResidVector & vResiduals, TFloat & lambda, CEForOptOnManifold & E_atManifold) {
        const TFloat BOOST = (TFloat) PARAMS.BOOST, DROP = (TFloat) PARAMS.DROP;
        const bool USE_HESSIAN = PARAMS.USE_HESSIAN;
        const bool USE_MARQUARDT = PARAMS.USE_MARQUARDT;

        TJTJMatrix JTJ = J.transpose() * J;
        TResidVector vNewResids;

        if (USE_HESSIAN) {
            TJTJMatrix H;
            make_rH(ax, axp, params, vResiduals, H);
            JTJ += H;
            //cout << JTJ << endl;
            if (IS_DEBUG) {
                const TFloat delta = 0.001f;
                const int k = 2;
                const TFloat x = vResiduals.squaredNorm();
                CParamVector newParams = params;
                newParams(k) -= delta;
                calcResiduals(ax, axp, newParams, vNewResids, E_atManifold);
                const TFloat x_minus = vNewResids.squaredNorm();
                newParams(k) += 2 * delta;
                calcResiduals(ax, axp, newParams, vNewResids, E_atManifold);
                const TFloat x_plus = vNewResids.squaredNorm();

                //const TFloat d = (x_plus - x_minus) / (2 * delta);
                const TFloat d2 = (x_plus + x_minus - 2 * x) / (delta * delta);
                //cout << "d=" << d <<endl;
                //cout << "d2=" << d2 <<endl;
                if (USE_HESSIAN && fabs(d2) > 0.001) {
                    CHECK(!zero(fabs(1 - 2 * JTJ(k, k) / d2)), "2nd derivative failed");
                }
            }
        }

        TFloat lambda_old = 0;
        const bool DEBUG_LAMBDA = false;

        for (;;) //for(int nInnerLoopCount = 0; nInnerLoopCount<10; nInnerLoopCount++)
        {
            if (USE_MARQUARDT)
                JTJ.diagonal() *= (1 + lambda) / (1 + lambda_old);
            else
                JTJ.diagonal().array() += lambda - lambda_old;

            TParamVector paramUpdateVec;
            if (NUM_POINTS < 5) {
                TJTJMatrix inv = JTJ.inverse(); //50% of cost
                paramUpdateVec = inv * J.transpose() * vResiduals;
                //const TJTJMatrix inv2 = JTJ.transpose().inverse();
            } else {
                const TParamVector & temp = J.transpose() * vResiduals;
                //paramUpdateVec = JTJ.partialPivLu().solve(temp);
                paramUpdateVec = JTJ.llt().solve(temp);
                //paramUpdateVec = JTJ.ldlt().solve(temp); No difference from llt
                //cout << paramUpdateVec << endl;
                if(IS_DEBUG) CHECK(!zero((JTJ * paramUpdateVec - temp).squaredNorm()), "LLT failed");
            }

            CParamVector newParams;
            newParams = params - paramUpdateVec;
            newParams.Normalise();

            calcResiduals(ax, axp, newParams, vNewResids, E_atManifold);

            const TFloat dNewErr = vNewResids.squaredNorm();
            //cout << "New err = " << dNewErr << endl;
            //if(DEBUG_LAMBDA)
            //cout << nIter << "- << nInnerLoopCount" << ": " << dNewErr << "=err, lambda=" << lambda << endl;

            if (dNewErr <= dErr && !DEBUG_LAMBDA) {
                params = newParams;
                vResiduals = vNewResids;
                lambda *= DROP;
                return dNewErr;
            }
            if (lambda > 1e+5)
                return dErr;
            lambda_old = lambda;
            lambda *= BOOST;
            //Does not help: TFloat dUpdateSize = paramUpdateVec.squaredNorm();
            //if(dUpdateSize < sqr(DROPOUT_THRESH) )
            //return false; //Fail (We're close enough to the minimum. Could breakout without using new (worse) params)
        }
        THROW("Broke from infinite loop");
    }

    static bool getHypothesisGradientDesc(const CEssentialMatGradientDescParams & PARAMS, const THomoImPoint * ax, const THomoImPoint * axp, CParamVector & params, CExistingSolns & aExistingSolns, CGDStats & stats, CEForOptOnManifold & E_atManifold) {
        const int MAX_ITERS = PARAMS.MAX_ITERS;

        TFloat lambda = (TFloat) PARAMS.LM.LAMBDA;
        const bool VERBOSE = PARAMS.VERBOSE;

        TResidVector vResiduals;
        TJMatrix J;

        calcResiduals(ax, axp, params, vResiduals, E_atManifold);
        TFloat dErr = vResiduals.squaredNorm();

        CParamVector LAMBDAx_lastIter, oldDelX;
        LAMBDAx_lastIter.setZero();
        oldDelX.setOnes();

        for (int nIter = 1; nIter <= MAX_ITERS; nIter++) {
            //Need to map E=(q,t) -> chi[3], tau[2]
            //params is the alphas and betas. params.setZero()

            makeJ(ax, axp, params, J, E_atManifold);

            TFloat dNewErr = 1e+10;
            switch (PARAMS.GDAlg) {
                case CEssentialMatGradientDescParams::eLevenbergMarquardt:
                    dNewErr = minLevMar(PARAMS.LM, ax, axp, J, dErr, params, vResiduals, lambda, E_atManifold); //vResiduals, lambda should probably be the state of a class
                    break;
                case CEssentialMatGradientDescParams::eSteepestDescent:
                    dNewErr = minSteepestDescent(PARAMS.SD, ax, axp, J, dErr, params, vResiduals, E_atManifold);
                    break;
                case CEssentialMatGradientDescParams::eCGD:
                    dNewErr = minCGD(PARAMS.SD, ax, axp, J, dErr, params, vResiduals, LAMBDAx_lastIter, oldDelX, E_atManifold);
                    break;

                default:
                    THROW("GD alg not implemented yet");
            }

            if (PARAM_MODE == eMinimalOnManifold)
                E_atManifold.updateE(params);

            //Err is proportional to number of points (RESIDS). Should ALWAYS go to 0 for a solution (are there false minima?)
            if (dNewErr < NUM_POINTS * 1e-20) // Must be < 1e-8 so that these points actually end up as inliers, changing to e-6 actually makes almost no difference to speed
            {
                if (VERBOSE) cout << "Converged in " << nIter << " iterations\n";
                stats.add(nIter, CGDStats::eSuccess);
                //aExistingSolns.contains(params, E_atManifold, PARAMS.DIST_TO_EXISTING_SOLN_THRESH, false);
                return true;
            }

            //If it doesn't look like it will hit 0 then stop: dNewErr / (dErr-dNewErr) is num of steps required at this rate of descent
            if (dNewErr > (dErr - dNewErr) * PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA/* this is about 1*/ * (MAX_ITERS - nIter)) {
                if (VERBOSE) cout << "Converged to false minimum\n";
                stats.add(nIter, CGDStats::eFalseMin);
                return false;
            }

            if (aExistingSolns.contains(params, E_atManifold, PARAMS.DIST_TO_EXISTING_SOLN_THRESH, false)) {
                if (VERBOSE) cout << "Already found this minima\n";
                stats.add(nIter, CGDStats::eRepeatedMin);
                return false;
            }

            dErr = dNewErr;
        }
        if (VERBOSE) cout << "Failed to converge\n";
        return false;
    }

public:

    static int makeHypotheses(const CEssentialMatGradientDescParams & PARAMS, const int nHypotheses, const THomoImPoint * ax, const THomoImPoint * axp, CModels & pModels, CGDStats & stats) {
        ARRAY(TEssentialMat, aE, nHypotheses);

        int nModelFound = 0;
        CExistingSolns aExistingSolns;
        aExistingSolns.reserve((NUM_POINTS == 5) ? nHypotheses : 0);
        for (int nHyp = 0; nHyp < nHypotheses; nHyp++) {
            const int MAX_ATTEMPTS = PARAMS.MAX_FALSE_MINIMA; //fails at false minima
            int nAttempts = 0;
            for (; nAttempts < MAX_ATTEMPTS; nAttempts++) {
                if (getHypothesis(PARAMS, ax, axp, aE[nModelFound], aExistingSolns, stats))
                    break;
            }
            if (nAttempts == MAX_ATTEMPTS)
                break;

            bool bAlreadyFound = false;
#ifdef _DEBUG
            for (int i = 0; i < NUM_POINTS; i++) {
                TFloat dResid = axp[i].transpose() * aE[nModelFound] * ax[i];
                //cout << dResid << ' ';
                if(IS_DEBUG) CHECK(!zero(dResid), "xp E x not 0");
            }
            //cout << endl;

            if (NUM_POINTS == 5) //Might find same minimum multiple times
            {
                TFloat dMinDist = 100;
                for (int nOldHyp = 0; nOldHyp < nModelFound; nOldHyp++) {
                    const TFloat dDist = (aE[nOldHyp] - aE[nModelFound]).squaredNorm();

                    if (dMinDist > dDist)
                        dMinDist = dDist;

                    if (zero(dDist)) {
                        THROW("Should already have detected duplicate solution\n"); // << aE[nOldHyp]  << endl << aE[nModelFound] << endl;
                        bAlreadyFound = true;
                        break;
                    }
                }
                //cout << nModelFound << ": " << dMinDist << endl; Actually pretty obvious which are identical
            }
#endif
            if (!bAlreadyFound) {
                C3x3MatModel & E = static_cast<C3x3MatModel &> (pModels.addModel());
                for (int r = 0; r < 3; r++)
                    for (int c = 0; c < 3; c++)
                        E(r, c) = aE[nModelFound](r, c);
                nModelFound++;
            }
        }

        return nModelFound;
    }
};

//For comparison:
//int calcEssentialMat_5point_Eigen(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & pdEssentialMat, const double dUprightThresh, int & FAIL_COUNT, const bool bNister);

typedef double TFloat;
typedef CLevMarForEHypothesis<TFloat, 5, eMinimal>::THomoImPoint THomoImPoint;

int findE_NRforBasis(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, TEModels & models, CGDStats & stats);

double SCALE_EPS = 1;

void testELM() {
    int FAIL_COUNT = 0;

    SCALE_EPS = 1e+6; //seems to work well
    /*if(SCALE_EPS < 1e+12)
    SCALE_EPS *= 10;
    cout << SCALE_EPS << endl;
    double eps = Eigen::NumTraits<double>::epsilon() * SCALE_EPS;
    cout << eps << endl;*/


    const int NUM_POINTS = 7;
    CEssentialMatGradientDescParams PARAMS(0, 0);

    ARRAY(THomoImPoint, ax, NUM_POINTS);
    ARRAY(THomoImPoint, axp, NUM_POINTS);

    TSubSet anHypSet(NUM_POINTS);
    T2dPoints m0(NUM_POINTS), m1(NUM_POINTS);
    for (int i = 0; i < NUM_POINTS; i++)
        anHypSet[i] = i;

    //PARAMS.USE_MARQUARDT = false;
    //PARAMS.USE_HESSIAN = false; //Much slower to use H

    //PARAMS.LAMBDA = 0.5;
    PARAMS.VERBOSE = IS_DEBUG;

    PARAMS.GDAlg = CEssentialMatGradientDescParams::eLevenbergMarquardt;

    enum eAlgs {
        eMoM, // (17.4 with 10 starts)
        eMin, //14.6 / 10 starts
        eMinDegen, // 12.13 / 5 starts
        eNister, //7.15
        eStewenius, //18.69
        eGDBasis //4.24 or 6.39 on R3
    };

    eAlgs alg = eGDBasis;


    if (PARAMS.GDAlg != CEssentialMatGradientDescParams::eLevenbergMarquardt) {
        PARAMS.MAX_ITERS = 1000;
        PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA = 6;
        PARAMS.MAX_ITERS = 80;
        PARAMS.SD.BETA_FN = CEssentialMatGradientDescParams::CSDParams::ePolakRibiere;
    }


    //For New LM on R3
    bool bNewBasisR3 = false;
    if (alg == eGDBasis) {
        if (bNewBasisR3) //best 7.07
        {
            PARAMS.LM.LAMBDA = 0.005;
            PARAMS.DIST_TO_EXISTING_SOLN_THRESH = 0.05;
            PARAMS.LM.DROP = 0.35;
            PARAMS.LM.BOOST = 5;
            PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 8;
        } else { //For New LM
            PARAMS.LM.LAMBDA = 0.0044;
            PARAMS.LM.DROP = 0.25;
            PARAMS.LM.BOOST = 10;

            PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 10;
            PARAMS.DIST_TO_EXISTING_SOLN_THRESH = 0.97;
        }
    }

    PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 10;

    CGDStats stats;
    T3x3MatModels models10;
    TEModels models;

    cv::Mat points0(7, 2, CV_32FC1), points1(7, 2, CV_32FC1);

    //for(PARAMS.DIST_TO_EXISTING_SOLN_THRESH = 0.95; PARAMS.DIST_TO_EXISTING_SOLN_THRESH <= 1; PARAMS.DIST_TO_EXISTING_SOLN_THRESH = PARAMS.DIST_TO_EXISTING_SOLN_THRESH + 0.005)
    //for(PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA = 0.3; PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA <= 15; PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA = PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA + 0.1)
    for (PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 5; PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS <= 20; PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS + 1)
        //for(PARAMS.MAX_FALSE_MINIMA = 1; PARAMS.MAX_FALSE_MINIMA < 20; PARAMS.MAX_FALSE_MINIMA = PARAMS.MAX_FALSE_MINIMA + 1)
        //for(PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA = 0.01; PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA < 20; PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA = PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA * 2)
        //for(PARAMS.MAX_ITERS = 4; PARAMS.MAX_ITERS < 20; PARAMS.MAX_ITERS = PARAMS.MAX_ITERS + 1)
        //for(PARAMS.LM.DROP = 0.15; PARAMS.LM.DROP < 0.5; PARAMS.LM.DROP = PARAMS.LM.DROP + 0.025)
        //for(PARAMS.LM.BOOST = 8; PARAMS.LM.BOOST < 12; PARAMS.LM.BOOST = PARAMS.LM.BOOST + 0.5)
        //for(PARAMS.LM.LAMBDA = 0.003; PARAMS.LM.LAMBDA < 400; PARAMS.LM.LAMBDA = PARAMS.LM.LAMBDA *1.1)
        //for(PARAMS.SD.SMALL_STEP_DROPOUT = 1e-16; PARAMS.SD.SMALL_STEP_DROPOUT < 1; PARAMS.SD.SMALL_STEP_DROPOUT = PARAMS.SD.SMALL_STEP_DROPOUT * 2)
        //PARAMS.LM.LAMBDA = 0;
        //for(PARAMS.FAILURE_PROB_THRESH = 0.1; PARAMS.FAILURE_PROB_THRESH < 1.0; PARAMS.FAILURE_PROB_THRESH = PARAMS.FAILURE_PROB_THRESH+0.1)
        //for(PARAMS.MAX_SOLUTIONS = 1; PARAMS.MAX_SOLUTIONS < 10; PARAMS.MAX_SOLUTIONS = PARAMS.MAX_SOLUTIONS + 1)
    {
        int numCalls = 0;
        CStopWatch s;
        s.startTimer();
        int totalModels = 0;
        while (totalModels < 1000000) {
            numCalls++;

            for (int i = 0; i < NUM_POINTS; i++) {
                ax[i] = THomoImPoint((TFloat) CRandom::FasterNormal2(), (TFloat) CRandom::FasterNormal2(), 1.0);
                axp[i] = THomoImPoint((TFloat) CRandom::FasterNormal2(), (TFloat) CRandom::FasterNormal2(), 1.0);
                m0[i] = CSimple2dPoint(ax[i](0), ax[i](1));
                m1[i] = CSimple2dPoint(axp[i](0), axp[i](1));
                if (NUM_POINTS == 7) {
                    points0.at<float>(i, 0) = (float)ax[i](0);
                    points0.at<float>(i, 1) = (float)ax[i](1);
                    points1.at<float>(i, 0) = (float)axp[i](0);
                    points1.at<float>(i, 1) = (float)axp[i](1);
                }
            }
            models10.reset();
            models.reset();
            int numModels = 0;

            if (NUM_POINTS == 7) {
                //Test timing
                //cv::findFundamentalMat(points1, points1, CV_FM_7POINT);
                //cvFindFundamentalMat()

                calcEssentialMat_7point(anHypSet, m0, m1, models10, FAIL_COUNT);
                numModels++;
            } else {
                if (alg == eMoM)
                    numModels = CLevMarForEHypothesis<TFloat, NUM_POINTS, eMinimalOnManifold>::makeHypotheses(PARAMS, PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS, PTR(ax), PTR(axp), models, stats);
                else if (alg == eMin)
                    numModels = CLevMarForEHypothesis<TFloat, NUM_POINTS, eMinimal>::makeHypotheses(PARAMS, PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS, PTR(ax), PTR(axp), models, stats);
                else if (alg == eMinDegen)
                    numModels = CLevMarForEHypothesis<TFloat, NUM_POINTS, eMinimalDegen>::makeHypotheses(PARAMS, PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS, PTR(ax), PTR(axp), models, stats);
                else if (alg == eNister || alg == eStewenius)
                    numModels = calcEssentialMat_5point_Eigen(anHypSet, m0, m1, models10, -1, FAIL_COUNT, alg == eNister);
                else if (alg == eGDBasis) {
                    if (bNewBasisR3)
                        numModels = findE_LMforBasisR3(PARAMS, anHypSet, m0, m1, models, stats);
                    else {
                        numModels = findE_LMforBasis(PARAMS, anHypSet, m0, m1, models, stats);
                        //numModels = findE_NRforBasis(PARAMS, anHypSet, m0, m1, models, stats);
                    }
                }
            }

            totalModels += numModels;
            stats.reset();
        };
        s.stopTimer();
        cout << "Time: " << s.getElapsedTime() << ", " << FAIL_COUNT << "=fails, num calls=" << numCalls << endl;
    }

    stats.writeResults(PARAMS.GDAlg, "E_GradientDescResults.tsv");
};

int CEssentialMatLM::getModels_int(const TSubSet & anHypSet, CModels & pModels) {
    //testELM();

    ARRAY(THomoImPoint, ax, numPoints());
    ARRAY(THomoImPoint, axp, numPoints());
    for (int i = 0; i < numPoints(); i++) {
        const int idx = anHypSet[i];
        ax[i] = THomoImPoint(p1[idx].getX(), p1[idx].getY(), 1.0);
        axp[i] = THomoImPoint(p2[idx].getX(), p2[idx].getY(), 1.0);
    }

    CEssentialMatGradientDescParams PARAMS(0, 0);
    PARAMS.VERBOSE = IS_DEBUG;

    CGDStats noStats;

    switch (numPoints()) {
        case 3:
            return CLevMarForEHypothesis<TFloat, 3, eMinimalDegen>::makeHypotheses(PARAMS, nHypotheses, PTR(ax), PTR(axp), pModels, noStats);
        case 4:
            return CLevMarForEHypothesis<TFloat, 4, eMinimalDegen>::makeHypotheses(PARAMS, nHypotheses, PTR(ax), PTR(axp), pModels, noStats);
        case 5:
            return CLevMarForEHypothesis<TFloat, 5, eMinimalOnManifold>::makeHypotheses(PARAMS, nHypotheses, PTR(ax), PTR(axp), pModels, noStats);
    }

    THROW("Unhandled numPoints()");
}

