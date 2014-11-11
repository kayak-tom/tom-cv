#define OVERRIDE_EPS
extern double SCALE_EPS;

//#define EIGEN_DONT_VECTORIZE

#include "util/exception.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>

#include "essentialMat_Eigen.h"
#include "util/convert.h"
#include "util/random.h"

#include "essentialMatLM.h" 
#include "LMParams.h" 
#include "util/dynArray.h" 
#include "makeBasis_GramSchmidt.h"
#include <fstream>

using namespace Eigen;
static const double ROOT_2 = sqrt(2.0);

void addModel(const Eigen::Matrix3d & Emodel, TEModels & models) {
    double dScale = ROOT_2;
    if (Emodel(0, 0) < 0) dScale = -ROOT_2;

    C3x3MatModel & M = static_cast<C3x3MatModel &> (models.addModel());

    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            M(r, c) = dScale * Emodel(r, c);
}


static const int NUM_PARAMS = 4, E_ROWS = 3, E_COLS = 3; //4 rows allows vectorisation, but not worthwhile

template<bool bAllResids, typename TFloat>
class CLevMarETE {
    typedef Eigen::Matrix<TFloat, NUM_PARAMS, 1 > TParamVec;
public:
    typedef CFixedArray<TParamVec, 100 > TParamsAlreadyFound; //TODO resize me

    static const int NUM_TS_PARAMS = NUM_PARAMS - 1,
            NUM_RESIDUALS = bAllResids ? 9 : 3;

    typedef Eigen::Matrix<double, 9, 1 > TEFullvec;
    typedef Eigen::Matrix<TFloat, NUM_RESIDUALS, 1 > TEvec;
    typedef Eigen::Matrix<TFloat, E_ROWS, E_COLS, Eigen::RowMajor> TE;
    typedef Eigen::Matrix<TFloat, E_ROWS, E_ROWS, Eigen::RowMajor> TEET;
    //typedef Eigen::Matrix<TFloat, E_ROWS, RESID> TEETE;
    typedef Eigen::Matrix<TFloat, 2, 3 > TE23; //Rows 2 and 3 from a 3x3 mat
    typedef Eigen::Matrix<double, 3, 3 > TE33;

    typedef Eigen::Matrix<TFloat, NUM_RESIDUALS, NUM_TS_PARAMS> TJ;
    typedef Eigen::Matrix<TFloat, NUM_TS_PARAMS, NUM_TS_PARAMS> TJTJ;

    typedef Eigen::Matrix<TFloat, NUM_TS_PARAMS, 1 > TTSParamVec;

    static const int RESID_MAT_ROWS = bAllResids ? E_ROWS : 1;
    typedef Eigen::Matrix<TFloat, RESID_MAT_ROWS, E_ROWS, Eigen::RowMajor> TResidTemp; //All, or top row, of EET Sometimes smaller if only using 3d resids
    typedef Eigen::Matrix<TFloat, RESID_MAT_ROWS, E_COLS> TResids; //All, or top row of full expression EETE + ... 

    static inline void toMat(const TEFullvec & vec, TE & mat) {
        int i = 0;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) {
                mat(r, c) = (TFloat) vec(i);
                i++;
            }

        if (E_ROWS == 4)
            mat.row(3).setZero();
        if (E_COLS == 4)
            mat.col(3).setZero();
    }
private:

    static inline void toVec(const TResids & mat, TEvec & vec) {
        int i = 0;
        for (int r = 0; r < ((RESID_MAT_ROWS >= 3) ? 3 : 1); r++)
            for (int c = 0; c < 3; c++) {
                vec(i) = mat(r, c);
                i++;
            }
    }

    class CEfunction {
        TE Ehere, aEbasis[NUM_PARAMS], aE_TSbasis[NUM_TS_PARAMS];
        TParamVec params, aParamTangentVecs[NUM_TS_PARAMS];

        void updateBasis(const bool bUseBinomialApprox) {
            //Make 3 vectors perp. to params.
            for (int i = 0; i < NUM_TS_PARAMS; i++) {
                aParamTangentVecs[i] -= aParamTangentVecs[i].dot(params) * params;
                for (int j = 0; j < i; j++)
                    aParamTangentVecs[i] -= aParamTangentVecs[i].dot(aParamTangentVecs[j]) * aParamTangentVecs[j];

                if (bUseBinomialApprox)
                    aParamTangentVecs[i] *= 1.5 - 0.5 * aParamTangentVecs[i].squaredNorm();
                else
                    aParamTangentVecs[i] /= aParamTangentVecs[i].stableNorm();

                //Make linear combination of basis vecs for E in this direction.
                aE_TSbasis[i].setZero();
                for (int j = 0; j < NUM_PARAMS; j++)
                    aE_TSbasis[i] += aEbasis[j] * aParamTangentVecs[i](j);
            }
        }
        TResidTemp EET;
        TE33 modelSoCanReturnRef;
    public:

        CEfunction(const Matrix<double, 9, 4 > & EE) {
            for (int i = 0; i < NUM_PARAMS; i++) {
                toMat(EE.col(i), aEbasis[i]);
            }

            for (int i = 0; i < NUM_TS_PARAMS; i++) {
                //aParamTangentVecs[i].setRandom();
                for (int r = 0; r < TParamVec::RowsAtCompileTime; r++)
                    for (int c = 0; c < TParamVec::ColsAtCompileTime; c++)
                        aParamTangentVecs[i](r, c) = CRandom::Uniform();
            }

            //Set a starting point:
            for (int i = 0; i < NUM_PARAMS; i++)
                params(i) = CRandom::FasterNormalCLT();

            TTSParamVec paramUpdate;
            paramUpdate.setZero();
            updateE_init(paramUpdate, false); //Initialise E
        }

        inline const TE33 & model() {
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                    modelSoCanReturnRef(r, c) = Ehere(r, c);
            return modelSoCanReturnRef;
        }

        const TParamVec & getParams() const {
            return params;
        }

        static TFloat traceEET(const TE & Ein) {
            return Ein.squaredNorm();
        }

        static TFloat traceE1E2T(const TE & E1, const TE & E2) {
            return (E1.array() * E2.array()).sum();
        }

        //Full EET

        inline static void makeEET(const TE & Ein, TEET & EET) {
            EET = Ein * Ein.transpose();
            //cout << EET << endl;
            const TFloat halfTrace = 0.5 * EET.trace();
            for (int i = 0; i < 3; i++)
                EET(i, i) -= halfTrace;
            //cout << EET << endl;
        }

        //For reduced residuals:

        inline static void makeEET(const TE & Ein, Eigen::Matrix<TFloat, 1, TResidTemp::ColsAtCompileTime> & EET) {
            EET = Ein.row(0) * Ein.transpose();
            EET(0) -= 0.5 * (EET(0) + Ein.row(1).squaredNorm() + Ein.row(2).squaredNorm());
        }

        inline static void makeEET23(const TE & Ein, TE23 & EET23) {
            EET23 = Ein.block < 2, 3 > (1, 0) * Ein.transpose();
            TFloat halfTrace = 0.5 * (EET23(0, 1) + EET23(1, 2) + Ein.row(0).squaredNorm());
            EET23(0, 1) -= halfTrace;
            EET23(1, 2) -= halfTrace;
        }

        void getResiduals(const TE & Ein, TResids & residuals) {
            //residuals = Ein * Ein.transpose() * Ein - 0.5 * traceEET(Ein) * Ein;

            //Faster? We remember and re-use EET
            makeEET(Ein, EET);
            residuals = EET * Ein;
            //cout << residuals << endl;
        }

        //Check that other residuals are also going to zero

        bool checkModelIsGood(const TFloat dResids) const {
            TE23 EET23;
            makeEET23(Ehere, EET23);
            TFloat dSS = (EET23 * Ehere).squaredNorm();
            cout << dSS << " " << dResids << endl;
            return dSS < 10000 * dResids;
        }

        void getResiduals(TResids & residuals) {
            getResiduals(Ehere, residuals);
        }

        void newE(const TParamVec & newParams, TE & E) const {
            E.setZero();
            for (int i = 0; i < NUM_PARAMS; i++)
                E += newParams(i) * aEbasis[i];
        }

        void updateE_init(const TTSParamVec & paramUpdate, const bool bUseBinomialApprox) {
            for (int i = 0; i < NUM_TS_PARAMS; i++)
                params += paramUpdate(i) * aParamTangentVecs[i];

            if (bUseBinomialApprox) {
                params *= 1.5 - 0.5 * params.squaredNorm();
                //cout << params.squaredNorm() << endl;
            } else {
                params /= params.stableNorm();
                //cout << params.transpose() << endl;
            }

            newE(params, Ehere);
            updateBasis(bUseBinomialApprox);
        }

        TParamVec newParams;
        TE Enew;

        void getNewResiduals(const TTSParamVec & paramUpdate, TResids & residuals) {

            newParams = params;
            for (int i = 0; i < NUM_TS_PARAMS; i++)
                newParams += paramUpdate(i) * aParamTangentVecs[i];

            TFloat dSN = newParams.squaredNorm();
            if (fabs(dSN - 1) < 0.1)
                newParams *= 1.5 - 0.5 * dSN; //normalise (binom approx)
            else
                newParams /= sqrt(dSN);

            if (IS_DEBUG && newParams.squaredNorm() < 0.9)
                cout << "Error, normalisation failed" << endl;

            newE(newParams, Enew);

            getResiduals(Enew, residuals);
        }

        //Once error is lower, accept most recent params and Ehere

        void acceptNew() {
            //TFloat dDiff = (params-newParams).squaredNorm();
            params = newParams;
            Ehere = Enew;

            //params /= params.stableNorm(); //Check working correctly

            //if(dDiff > 0.2)
            //updateBasis(false);
        }

        static inline void computeAplusAT(const TEET & A, TEET & AplusAT) {
            AplusAT = A + A.transpose();
            const TFloat trace = A(0, 0) + A(1, 1) + A(2, 2);
            for (int i = 0; i < 3; i++)
                AplusAT(i, i) -= trace;
        }

        static inline void computeAplusAT(const TEET & A, Eigen::Matrix<TFloat, 1, E_ROWS> & AplusAT) {
            AplusAT = A.row(0) + A.col(0).transpose();
            AplusAT(0) -= A.trace();
        }

        //Only works at E

        void getJ(TJ & J) const {
            const TE & E0 = Ehere;

            /*if (IS_DEBUG) {
                TE Etest = E0 * E0.transpose();
                Etest.diagonal().array() -= 0.5 * traceEET(E0);
                CHECK((Etest - EET).squaredNorm() > 1e-10, "Miss-remembered EET-.5 tr");
            }*/

            for (int i = 0; i < NUM_TS_PARAMS; i++) {
                const TE & E1 = aE_TSbasis[i];

                //cout << E1 * E0.transpose() << endl;
                //cout << E0 * E1.transpose() << endl;

                //TE deriv_old = ((E1 * E0.transpose() + E0 * E1.transpose()) * E0 + E0 * E0.transpose() * E1
                //        - (traceE1E2T(E0, E1) * E0 + 0.5 * traceEET(E0) * E1));

                TEET A = E0 * E1.transpose();
                TResidTemp AplusAT;
                computeAplusAT(A, AplusAT);

                TResids deriv = AplusAT * E0 + EET*E1;

                TEvec vecTemp;
                toVec(deriv, vecTemp);
                J.col(i) = vecTemp;

                //and test:
                /*if(false){
                    cout << J.col(i) << endl;
                    TFloat delta = 0.001;
                    TE r_plusDelta, r;
                    getResiduals(E0, r);
                    TE E0_delta = (E0 + delta * E1);
                    getResiduals(E0_delta, r_plusDelta);
                    r_plusDelta -= r;
                    r_plusDelta.array() /= delta;
                    cout << r_plusDelta << endl;
                }*/
            }
        }


    };

    static TFloat minLevMar(const CEssentialMatGradientDescParams::CLMParams & PARAMS, const TJ & J, CEfunction & E, TFloat & lambda, TResids & vResiduals, const TFloat dErr) {
        const TFloat BOOST = (double) PARAMS.BOOST, DROP = (double) PARAMS.DROP;
        //const bool USE_HESSIAN=PARAMS.USE_HESSIAN;
        const bool USE_MARQUARDT = PARAMS.USE_MARQUARDT;
        TTSParamVec paramsTS;
        paramsTS.setZero();

        //TE vResiduals;
        //E.getResiduals(vResiduals); //following normalisation
        //const double dErr = vResiduals.squaredNorm();

        TJTJ JTJ = J.transpose() * J;
        TResids vNewResids;

        TFloat lambda_old = 0;

        for (;;) {
            if (USE_MARQUARDT)
                JTJ.diagonal() *= (1 + lambda) / (1 + lambda_old);
            else
                JTJ.diagonal().array() += lambda - lambda_old;

            TTSParamVec paramUpdateVec;
            TEvec residualVec;
            toVec(vResiduals, residualVec);

            if (sizeof (TFloat) == sizeof (double)) {
                TJTJ inv = JTJ.inverse(); //50% of cost
                paramsTS = -inv * J.transpose() * residualVec;
            } else {
                const TTSParamVec & temp = J.transpose() * residualVec;
                //paramUpdateVec = JTJ.partialPivLu().solve(temp);
                paramsTS = -JTJ.llt().solve(temp);
                //cout << paramUpdateVec << endl;
                if(IS_DEBUG) CHECK(!zero((JTJ * paramsTS + temp).squaredNorm()), "LLT failed");
            }

            E.getNewResiduals(paramsTS, vNewResids);

            const TFloat dNewErr = vNewResids.squaredNorm();
            //cout << "New err = " << dNewErr << endl;
            //if(DEBUG_LAMBDA)
            //cout << nIter << "- << nInnerLoopCount" << ": " << dNewErr << "=err, lambda=" << lambda << endl;

            if (dNewErr <= dErr || lambda == 0) {
                vResiduals = vNewResids;
                lambda *= DROP;
                E.acceptNew();
                //E.updateE(paramsTS, true);
                return dNewErr;
            }
            lambda_old = lambda;
            lambda *= BOOST;

            if (lambda > 1e+5)
                return dErr;
            //Does not help: TFloat dUpdateSize = paramUpdateVec.squaredNorm();
            //if(dUpdateSize < sqr(DROPOUT_THRESH) )
            //return false; //Fail (We're close enough to the minimum. Could breakout without using new (worse) params)
        }
        THROW("Broke from infinite loop");

    }

    /*static bool contains(const TEModels & models, TE Enew, const TFloat THRESH) {
        if(Enew(0) < 0)
            Enew *= -1;

        for (int nModel = 0; nModel < models.numModels(); nModel++) {
            const C3x3MatModel & m = static_cast<const C3x3MatModel &> (models.getData(nModel));
            if(IS_DEBUG) CHECK(sizeof (TE) != 9 * sizeof (TFloat), "Dodgy cast will fail");
            TE E(m.asDouble9());
            E.transposeInPlace();
            if ((Enew - E).squaredNorm() < THRESH)
                return true;
        }
        return false;
    }*/
public:

    //bool calcEssentialMat_LM_Basis(const CEssentialMatGradientDescParams & PARAMS, const Matrix<double, 9, 4 > & EE, TEModels & models, CGDStats & stats) HOT;

    static bool paramsAlreadyFound(const TParamVec & params, const TParamsAlreadyFound & aParamsAlreadyFound, const TFloat THRESH) {
        //cout << params.stableNorm() << endl;
        for (typename TParamsAlreadyFound::const_iterator pParamsFound = aParamsAlreadyFound.begin(); pParamsFound != aParamsAlreadyFound.end(); pParamsFound++) {
            TFloat dDist = params.dot(*pParamsFound);
            if (fabs(dDist) > THRESH)
                return true;
        }
        return false;
    }

    static void rememberParams(const TParamVec & params, TParamsAlreadyFound & aParamsAlreadyFound) //Don't try here again
    {
        aParamsAlreadyFound.push_back(params); 
        //REPEAT(1, cout << "MUST PUT THIS BACK\n";);

        if(IS_DEBUG && (params.squaredNorm() < 0.95 || params.squaredNorm() > 1.05))
            cout << "Param vector isn't normalised " << params.norm() << endl;
        //if(IS_DEBUG) CHECK(params.squaredNorm() < 0.95 || params.squaredNorm() > 1.05, "Param vector isn't normalised")
    }

    static bool calcEssentialMat_LM_Basis(const CEssentialMatGradientDescParams & PARAMS, const Matrix<double, 9, 4 > & EE, TEModels & models, CGDStats & stats, const Eigen::Vector3d * am0, const Eigen::RowVector3d * am1, TParamsAlreadyFound & aParamsAlreadyFound) {
        CEfunction E(EE);

        TFloat lambda = (double) PARAMS.LM.LAMBDA;

        TResids resids;
        E.getResiduals(resids);
        TFloat dErr = resids.squaredNorm();

        bool /*bAlreadyChecked = false,*/ bBad = false;

        const int MAX_ITERS = PARAMS.MAX_ITERS;
        for (int nIter = 0; nIter < MAX_ITERS; nIter++) {
            TJ J;
            E.getJ(J);

            TFloat dNewErr = 1e+10;
            //switch(PARAMS.GDAlg)
            //{
            //case CEssentialMatGradientDescParams::eLevenbergMarquardt:
            dNewErr = minLevMar(PARAMS.LM, J, E, lambda, resids, dErr);
            /*	break;
            default:
                    THROW("GD alg not implemented yet");
            }*/

            bool bFalseMin = false;
            //cout << nIter << ": ";
            //if(!E.checkModelIsGood(dNewErr)) 
            //  bFalseMin=true;
            //else

            //Err is proportional to number of points (RESIDS). Should ALWAYS go to 0 for a solution (are there false minima?)
            if (dNewErr < NUM_RESIDUALS * 1e-20) // Must be < 1e-8 so that these points actually end up as inliers, changing to e-6 actually makes almost no difference to speed
            {
                TE33 Emodel = E.model();

                //cout << "B:" << Emodel << endl << endl;
                //makeClosestE(Emodel);
                //cout << "A:" << Emodel << endl << endl;
                TFloat dResid = getETEResid(Emodel);

                //TFloat dResid = getResids<double>(Emodel.transpose(), am0, am1);
                if (dResid > 0.00001) {
                    //cout << "Bad model\n";
                    bFalseMin = true;
                } else {
                    if (IS_DEBUG) {
                        if (bBad)
                            cout << "Good model flagged as bad\n";
                        Eigen::JacobiSVD<TE33> svd(Emodel);
                        if (svd.singularValues()(0) > 1.2 * svd.singularValues()(1) || svd.singularValues()(2) > 0.1)
                            cout << "BAD MODEL FOUND: " << svd.singularValues().transpose() << endl;
                    }

                    addModel(Emodel, models);

                    if (PARAMS.VERBOSE) cout << "Converged in " << nIter << " iterations\n";

                    //cout << Emodel << endl << endl;

                    stats.add(nIter, CGDStats::eSuccess);
                    rememberParams(E.getParams(), aParamsAlreadyFound);

                    return true;
                }
            }

            //Detect when we are heading for bad solutions -- Not worthwhile
            /*if(dNewErr < 1e-5 && !bAlreadyChecked)
            {
                bAlreadyChecked = true;
                const TFloat dResid = getETEResid(E.model());
                if(IS_DEBUG) cout << dNewErr << ": " << dResid / dNewErr << endl;
                if (dResid > 20*dNewErr )
                {
                    bFalseMin = true;
                    if(IS_DEBUG) cout << "Bad model\n";
                }
                else
                {
                    if(IS_DEBUG) cout << "Good model\n";
                }
            }*/


            //If it doesn't look like it will hit 0 then stop: dNewErr / (dErr-dNewErr) is num of steps required at this rate of descent
            if (bFalseMin || dNewErr > (dErr - dNewErr) * PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA/* this is about 1*/ * (MAX_ITERS - nIter)
                    || E.model().row(0).squaredNorm() < 1e-4) {
                if (PARAMS.VERBOSE) cout << "Converged to false minimum\n";
                stats.add(nIter, CGDStats::eFalseMin);

                rememberParams(E.getParams(), aParamsAlreadyFound); //Stay away from this FM in future
                return false;
            }

            /*if (contains(models, E.model(), PARAMS.DIST_TO_EXISTING_SOLN_THRESH)) {
                if (PARAMS.VERBOSE) cout << "Already found this minima\n";
                stats.add(nIter, CGDStats::eRepeatedMin);
                return false;
            }*/
            if (paramsAlreadyFound(E.getParams(), aParamsAlreadyFound, PARAMS.DIST_TO_EXISTING_SOLN_THRESH)) {
                if (PARAMS.VERBOSE) cout << "Already found this minima\n";
                stats.add(nIter, CGDStats::eRepeatedMin);
                return false;
            }
            dErr = dNewErr;
        }
        if (PARAMS.VERBOSE) cout << "Failed to converge\n";
        return false;
    }
};


typedef Matrix<double, 3, 3, Eigen::RowMajor> TE;

double getETEResid(const TE & E) {
    TE EET = E * E.transpose();
    EET.diagonal().array() -= 0.5 * EET.trace();
    EET *= E;
    return EET.squaredNorm();
}

void makeSC(const double a, const double b, double & s, double & c) {
    const double r = sqrt(sqr(a) + sqr(b));
    double r_inv = 1.0 / r;
    s = -b*r_inv;
    c = a*r_inv;
}

void makeSC2(const double a, const double b, double & s, double & c) {
    const double r = sqrt(sqr(a) + sqr(b));
    double r_inv = 1.0 / r;
    c = b*r_inv;
    s = a*r_inv;
}

void makeJacobiSC(const double akk, const double all, const double akl, double & s, double & c) {
    //from http://en.wikipedia.org/wiki/Jacobi_rotation
    const double beta = (all - akk) / (2 * akl);
    double t = 1.0 / (fabs(beta) + sqrt(sqr(beta) + 1));
    if (beta < 0) t *= -1;

    c = 1.0 / sqrt(sqr(t) + 1);
    s = c*t;
}

bool calcEssentialMat_Iterative_Basis(const Matrix<double, 9, 4 > & EE, TEModels & models, CGDStats & stats, const Eigen::Vector3d * am0, const Eigen::RowVector3d * am1, CLevMarETE < false, double>::TParamsAlreadyFound & aParamsAlreadyFound) {

    const int NUM_BASIS_VECS = 4;
    TE E, aEbasis[NUM_BASIS_VECS];
    Vector4d alpha;
    alpha.setRandom();

    for (int i = 0; i < NUM_BASIS_VECS; i++) {
        CLevMarETE < false, double>::toMat(EE.col(i), aEbasis[i]);
        alpha(i) = CRandom::FasterNormalCLT();
    }

    std::ofstream residFile;

    const bool bLog = false;
    const bool bVerbose = false;

    const double dLastResid = HUGE;

    if (bLog)
        residFile.open("resids.tsv", std::ios_base::app);

    TE Diag110;
    Diag110.setIdentity();
    Diag110(2, 2) = 0;

    const int MAX_ITERS = 125;

    for (int nIter = 0; nIter < MAX_ITERS; nIter++) {
        EigenNormalise(alpha); 

        if (CLevMarETE < false, double>::paramsAlreadyFound(alpha, aParamsAlreadyFound, 0.975))
            return false;

        E.setZero();
        for (int i = 0; i < NUM_BASIS_VECS; i++) {
            E += alpha(i) * aEbasis[i];
        }


        if (bVerbose) {
            cout << alpha.transpose() << endl << endl;
            cout << E << endl << endl;

            cout << getETEResid(E) << endl << endl;
        }
        if (bLog)
            residFile << getETEResid(E) << '\t';

        enum eUpdateMethods {
            eSVD, eAddDelta, eJacobiDiagonalise, eGivens, eGivensDiag, eSimple
        };

        const eUpdateMethods eUpdateMethod = eJacobiDiagonalise;
        if (eUpdateMethod == eSVD) {
            JacobiSVD<TE> svd(E, Eigen::ComputeFullU + Eigen::ComputeFullV);
            E = svd.matrixU() * Diag110 * svd.matrixV().transpose();
        } else if (eUpdateMethod == eSimple) {
            
            for(int i=0;i<10;i++)
            {
                TE EET = E * E.transpose();
                TE Res = EET;
                
                Res.diagonal().array() -= 0.5;
                Res *= E;

                EET.diagonal().array() += 0.5;
                E = 2*E - EET * E;
                
                cout << "Res: \n" << Res << endl;
            }
            
        } else if (eUpdateMethod == eAddDelta || eUpdateMethod == eJacobiDiagonalise) {
            TE EET = E * E.transpose();
            TE Res = EET;
            Res.diagonal().array() -= 0.5;
            Res *= E;

            if (bVerbose)
            {
                cout << "R: \n" << Res << endl << endl;
                cout << "EET: \n" << EET << endl << endl;
            }

            double dResid = Res.squaredNorm();
            if (dResid < 1e-16)
                break; //converged;

            if (dResid > (dLastResid - dResid)*(MAX_ITERS - nIter)) //Unlikely to make it to 0 (false minima)
            {
                CLevMarETE < false, double>::rememberParams(alpha, aParamsAlreadyFound);
                return false;
            }

            if (eUpdateMethod == eJacobiDiagonalise) {
                //E = -2*E*E.transpose()*E;
                //E -= R;

                //Form the svd of EET by 3 Jacobi rotations

                TE & T = EET, RV, RU;
                RV.setIdentity();
                RU.setIdentity();

                for (int i = 0; i < 3; i++) {
                    TE J;
                    double s, c;
                    double t10 = fabs(T(1, 0));
                    double t20 = fabs(T(2, 0));
                    double t12 = fabs(T(1, 2));

                    if (t10 > t20 && t10 > t12) {
                        makeJacobiSC(T(0, 0), T(1, 1), T(1, 0), s, c);
                        J << c, -s, 0, s, c, 0, 0, 0, 1;
                    } else if (t20 > t12) {
                        makeJacobiSC(T(0, 0), T(2, 2), T(2, 0), s, c);
                        J << c, 0, -s, 0, 1, 0, s, 0, c;
                    } else {
                        makeJacobiSC(T(1, 1), T(2, 2), T(2, 1), s, c);
                        J << 1, 0, 0, 0, c, -s, 0, s, c;
                    }

                    T = J * T * J.transpose();
                    //RV = J * RV;
                    RU *= J.transpose();

                    //cout << "T: \n" << T << endl;
                    //cout << "RU: \n" << RU << endl; //This is -U
                    //cout << "RV: \n" << RV << endl;
                }

                TE DiagOfE; DiagOfE.setIdentity();
                
                int idx = -1;
                if(T.diagonal().minCoeff(&idx) < 1e-14) //we're there
                    break;
                
                DiagOfE(idx, idx) = 0;
                                
                TE D_inv; D_inv.setZero(); D_inv.diagonal() = T.diagonal().array().sqrt().cwiseInverse();
                
                RV = D_inv * RU.transpose() * E;
                CHECKNAN(RV.sum());

                E = RU * DiagOfE * RV;
                
                if(false){
                    cout << "Eproj: " << RU * DiagOfE * RV << endl;

                    JacobiSVD<TE> svd(E, Eigen::ComputeFullU + Eigen::ComputeFullV);

                    cout << "Usvd: \n" << svd.matrixU() << endl;
                    cout << "Vsvd: \n" << svd.matrixV() << endl;
                    cout << "SVs: \n" << svd.singularValues() << endl;

                    E = svd.matrixU() * Diag110 * svd.matrixV().transpose();
                    cout << "E: \n" << E << endl;
                }
            } else {
                TE Q = EET + E.transpose() * E + E*E;
                Q.diagonal().array() -= 0.5;

                TE Delta = -Res * Q.inverse();

                if (bVerbose) {
                    cout << "Q: \n" << Q << endl << endl;
                    cout << "Delta: \n" << Delta << endl << endl;

                    cout << "EET: \n" << EET << endl << endl;
                    cout << "ETE \n" << E.transpose() * E << endl << endl;
                }

                E += Delta;
            }
        } else if (eUpdateMethod == eGivens) {
            TE P, G1, G2, G3;
            P << 0, 0, -1, 1, 0, 0, 0, 1, 0;

            TE T = P * E;

            cout << "T:\n" << T << endl;

            double s, c;

            makeSC(T(0, 0), T(1, 0), s, c);
            G1 << c, -s, 0, s, c, 0, 0, 0, 1;
            T = G1*T;

            cout << "G1:\n" << G1 << endl;
            cout << "T:\n" << T << endl;

            makeSC(T(1, 1), T(2, 1), s, c);
            G2 << 1, 0, 0, 0, c, -s, 0, s, c;
            T = G2*T;

            cout << "G2:\n" << G2 << endl;
            cout << "T:\n" << T << endl;

            makeSC(T(2, 2), T(0, 2), s, c);
            G3 << c, 0, -s, 0, 1, 0, s, 0, c;
            T = G3*T;

            cout << "G3:\n" << G3 << endl;
            cout << "T:\n" << T << endl;
        } else if (eUpdateMethod == eGivensDiag) {
            TE G1, G2, G3;

            TE T = E;

            cout << "T0:\n" << T << endl;

            for (int i = 0; i < 10; i++) {
                double t00 = fabs(T(0, 0));
                double t11 = fabs(T(1, 1));
                double t22 = fabs(T(2, 2));

                double s, c;

                if (t22 < t11 && t22 < t00) {
                    makeSC2(T(0, 0), T(1, 0), s, c);
                    G1 << c, -s, 0, s, c, 0, 0, 0, 1;
                    T = G1*T;
                } else if (t11 > t00) {
                    makeSC2(T(1, 1), T(2, 1), s, c);
                    G2 << 1, 0, 0, 0, c, -s, 0, s, c;
                    T = G2*T;
                } else {
                    makeSC2(T(2, 2), T(0, 2), s, c);
                    G3 << c, 0, s, 0, 1, 0, -s, 0, c;
                    T = G3*T;
                }
                cout << "T:\n" << T << endl;
            }

        }

        for (int i = 0; i < NUM_BASIS_VECS; i++) {
            alpha(i) = (aEbasis[i].array() * E.array()).sum();
        }
    }

    if (bLog) {
        residFile << endl;
        residFile.close();
    }

    EigenNormalise(alpha);
    CLevMarETE < false, double>::rememberParams(alpha, aParamsAlreadyFound);

    addModel(E, models);
    return true;
}

int findE_LMforBasis(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, TEModels & models, CGDStats & stats) {
    Matrix<double, 9, 4 > EE;
    makeBasisForE_GramSchmidt_Vectorise < 5 > (anHypSet, m0, m1, EE);
    CLevMarETE < false, double>::TParamsAlreadyFound aParamsAlreadyFound;

    //cout << "###############\n";

    Eigen::Vector3d am0[5];
    Eigen::RowVector3d am1[5];
    for (int i = 0; i < 5; i++) {
        am0[i] = Eigen::Vector3d(m1[anHypSet[i]].getX(), m1[anHypSet[i]].getY(), 1);
        am1[i] = Eigen::RowVector3d(m0[anHypSet[i]].getX(), m0[anHypSet[i]].getY(), 1);
    }

    for (int i = 0; i < PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS; i++) {
        const bool bUseIterative = false;

        if (bUseIterative)
            calcEssentialMat_Iterative_Basis(EE, models, stats, am0, am1, aParamsAlreadyFound);
        else
            CLevMarETE < false, double>::calcEssentialMat_LM_Basis(PARAMS, EE, models, stats, am0, am1, aParamsAlreadyFound);

        if (models.numModels() >= PARAMS.MAX_SOLUTIONS)
            break;
        /*if(i>2) Doesn't help
        {
                double dProbSuccessNextTimeEst = models.numModels()/i;
                if(PARAMS.FAILURE_PROB_THRESH > dProbSuccessNextTimeEst)
                        break;
        }*/
    }

    return models.numModels();
}

int CEssentialMatLMbasis::getModels_int(const TSubSet & anHypSet, CModels & models) {

    static CEssentialMatGradientDescParams PARAMS(0, 0);
    //REPEAT(1, { PARAMS.VERBOSE = IS_DEBUG;    PARAMS.LM.LAMBDA = 0.0002;    PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 8;    PARAMS.DIST_TO_EXISTING_SOLN_THRESH = 0.1;});
    //New params
   

    REPEAT(1,
    { 
        PARAMS.LM.LAMBDA = 0.0044; PARAMS.LM.DROP = 0.25; PARAMS.LM.BOOST = 10; PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 10; PARAMS.DIST_TO_EXISTING_SOLN_THRESH = 0.975;
        //PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 100;
        //PARAMS.MAX_ITERS=100;
    });

    CGDStats noStats;
    TEModels & modelsCast = static_cast<TEModels &> (models);
    return findE_LMforBasis(PARAMS, anHypSet, p1, p2, modelsCast, noStats);
}

