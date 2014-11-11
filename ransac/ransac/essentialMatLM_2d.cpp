#define OVERRIDE_EPS
extern double SCALE_EPS;

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

using namespace Eigen;

static const int NUM_PARAMS = 3;
typedef Eigen::Matrix<double, NUM_PARAMS, 1 > TParamVec;
typedef CFixedArray<TParamVec, 20> TParamsAlreadyFound;

template<bool bAllResids>
class CLevMarETE_R3 {
    static const int NUM_RESIDUALS = bAllResids ? 9 : 3;

    typedef Eigen::Matrix<double, 9, 1 > TEFullvec;
    typedef Eigen::Matrix<double, NUM_RESIDUALS, 1 > TEvec;
    typedef Eigen::Matrix3d TE;
    //typedef Eigen::Matrix<double, 2, 3 > TE23; //Rows 2 and 3 from a 3x3 mat

    typedef Eigen::Matrix<double, NUM_RESIDUALS, NUM_PARAMS> TJ;
    typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJ;

    static const int RESID_MAT_ROWS =  bAllResids ? 3 : 1;
    typedef Eigen::Matrix<double, RESID_MAT_ROWS, 3> TE_var; //Sometimes smaller if only using 3d resids
    typedef Eigen::Matrix<double, 3, RESID_MAT_ROWS> TResids;
    
    static const TE & toMat(const TEFullvec & vec) {
        return reinterpret_cast<const TE &> (vec);
    }

    static const TEvec & toVec(const TE_var & mat) {
        return reinterpret_cast<const TEvec &> (mat);
    }

    class CEfunction {
        TE Ehere, aEbasis[4];
        TParamVec params;
        
        TE_var EET;
    public:

        CEfunction(const Matrix<double, 9, 4 > & EE) {
            for (int i = 0; i < 4; i++) {
                aEbasis[i] = toMat(EE.col(i));
            }

            //Set a starting point:
            for (int i = 0; i < NUM_PARAMS; i++)
                params(i) = CRandom::FasterNormal2();
            //params /= CRandom::FasterNormal2(); //To ensure directions are uniform
            
            //aEbasis[3] *= CRandom::FasterNormal2();

            newE(params, Ehere);
        }

        const TE & model() const {
            return Ehere;
        }

        const TParamVec & getParams() const { return params; }

        static double traceEET(const TE & Ein) {
            return Ein.squaredNorm();
        }

        static double traceE1E2T(const TE & E1, const TE & E2) {
            return (E1.array() * E2.array()).sum();
        }

        inline static void makeEET(const TE & Ein, TE & EET)
        {
            EET = Ein * Ein.transpose();
            EET.diagonal().array() -= 0.5 * EET.trace();
        }
        inline static void makeEET(const TE & Ein, Eigen::Matrix<double, 1, 3> & EET)
        {
            EET = Ein.row(0) * Ein.transpose();
            EET(0) -= 0.5 * (EET(0) + Ein.block<2,3>(1,0).squaredNorm());
        }
        /*inline static void makeEET23(const TE & Ein, TE23 & EET23)
        {
            EET23 = Ein.block<2,3>(1,0) * Ein.transpose();
            double halfTrace = 0.5 * (EET23(0,1) + EET23(1,2) + Ein.row(0).squaredNorm());
            EET23(0,1) -= halfTrace;
            EET23(1,2) -= halfTrace;
        }*/
        
        void getResiduals(const TE & Ein, TE_var & residuals) {
            //residuals = Ein * Ein.transpose() * Ein - 0.5 * traceEET(Ein) * Ein;

            //Faster? We remember and re-use EET
            makeEET(Ein, EET);
            residuals = EET * Ein;
        }
        
        /*//Check that other residuals are also going to zero
        bool checkModelIsGood(const double dResids) const
        {
            TE23 EET23;
            makeEET23(Ehere, EET23);
            double dSS = (EET23 * Ehere).squaredNorm();
            cout << dSS << " " << dResids << endl;
            return dSS < 10000*dResids;
        }*/

        void getResiduals(TE_var & residuals) {
            getResiduals(Ehere, residuals);
        }

        void newE(const TParamVec & newParams, TE & E) const {
            E = aEbasis[3];
            for (int i = 0; i < NUM_PARAMS; i++)
                E += newParams(i) * aEbasis[i];
        }


        TParamVec newParams;
        TE Enew;

        void getNewResiduals(const TParamVec & paramUpdate, TE_var & residuals) {

            newParams = params + paramUpdate;
            
            newE(newParams, Enew);

            getResiduals(Enew, residuals);
        }

        //Once error is lower, accept most recent params and Ehere

        void acceptNew() {
            //double dDiff = (params-newParams).squaredNorm();
            params = newParams;
            Ehere = Enew;

            //params /= params.stableNorm(); //Check working correctly

            //if(dDiff > 0.2)
                //updateBasis(false);
        }

        static inline void computeAplusAT(const TE & A, TE & AplusAT)
        {
            AplusAT = A + A.transpose();
            AplusAT.diagonal().array() -= A.trace();
        }
        
        static inline void computeAplusAT(const TE & A, Eigen::Matrix<double, 1, 3> & AplusAT)
        {
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

            for (int i = 0; i < NUM_PARAMS; i++) {
                const TE & E1 = aEbasis[i];

                //cout << E1 * E0.transpose() << endl;
                //cout << E0 * E1.transpose() << endl;

                //TE deriv_old = ((E1 * E0.transpose() + E0 * E1.transpose()) * E0 + E0 * E0.transpose() * E1
                //        - (traceE1E2T(E0, E1) * E0 + 0.5 * traceEET(E0) * E1));

                TE A = E0 * E1.transpose();
                TE_var AplusAT;
                computeAplusAT(A, AplusAT);
                
                TE_var deriv = AplusAT * E0 + EET*E1;

                J.col(i) = toVec(deriv);

                //and test:
                /*if(false){
                    cout << J.col(i) << endl;
                    double delta = 0.001;
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

    static double minLevMar(const CEssentialMatGradientDescParams::CLMParams & PARAMS, const TJ & J, CEfunction & E, double & lambda, TE_var & vResiduals, const double dErr) {
        const double BOOST = PARAMS.BOOST, DROP = PARAMS.DROP;
        //const bool USE_HESSIAN=PARAMS.USE_HESSIAN;
        const bool USE_MARQUARDT = PARAMS.USE_MARQUARDT;
        TParamVec paramsTS;
        paramsTS.setZero();

        //TE vResiduals;
        //E.getResiduals(vResiduals); //following normalisation
        //const double dErr = vResiduals.squaredNorm();

        TJTJ JTJ = J.transpose() * J;
        TE_var vNewResids;

        double lambda_old = 0;

        for (;;) {
            if (USE_MARQUARDT)
                JTJ.diagonal() *= (1 + lambda) / (1 + lambda_old);
            else
                JTJ.diagonal().array() += lambda - lambda_old;

            //TParamVec paramUpdateVec;
            const TEvec residualVec(toVec(vResiduals));
            if (true) {
                TJTJ inv = JTJ.inverse(); 
                paramsTS = -inv * J.transpose() * residualVec;
            } else {
                const TParamVec & temp = J.transpose() * residualVec;
                //paramUpdateVec = JTJ.partialPivLu().solve(temp);
                paramsTS = -JTJ.llt().solve(temp);
                //cout << paramUpdateVec << endl;
                if(IS_DEBUG) CHECK(!zero((JTJ * paramsTS + temp).squaredNorm()), "LLT failed");
            }

            E.getNewResiduals(paramsTS, vNewResids);

            const double dNewErr = vNewResids.squaredNorm();
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

    static bool contains(const TEModels & models, TE Enew, const double THRESH) {
        if(Enew(0) < 0)
            Enew *= -1;

        for (int nModel = 0; nModel < models.numModels(); nModel++) {
            const C3x3MatModel & m = static_cast<const C3x3MatModel &> (models.getData(nModel));
            if(IS_DEBUG) CHECK(sizeof (TE) != 9 * sizeof (double), "Dodgy cast will fail");
            TE E(m.asDouble9());
            E.transposeInPlace();
            if ((Enew - E).squaredNorm() < THRESH)
                return true;
        }
        return false;
    }

    //bool calcEssentialMat_LM_Basis(const CEssentialMatGradientDescParams & PARAMS, const Matrix<double, 9, 4 > & EE, TEModels & models, CGDStats & stats) HOT;
    static bool paramsAlreadyFound(const TParamVec & params, const TParamsAlreadyFound & aParamsAlreadyFound, const double THRESH)
    {
        for(TParamsAlreadyFound::const_iterator pParamsFound = aParamsAlreadyFound.begin(); pParamsFound != aParamsAlreadyFound.end(); pParamsFound++)
        {
            double dDist = (params-*pParamsFound).squaredNorm();
            if(dDist < THRESH)
                return true;
        }
        return false;
    }

    static void rememberParams(const TParamVec & params, TParamsAlreadyFound & aParamsAlreadyFound) //Don't try here again
    {
        aParamsAlreadyFound.push_back(params);
        if(IS_DEBUG)
            cout << params.squaredNorm() << endl;
    }

public:
    static bool calcEssentialMat_LM_Basis(const CEssentialMatGradientDescParams & PARAMS, const Matrix<double, 9, 4 > & EE, TEModels & models, CGDStats & stats, const Eigen::Vector3d * am0, const Eigen::RowVector3d * am1, TParamsAlreadyFound & aParamsAlreadyFound) {
        CEfunction E(EE);
        //double dErr = 1e+10;
        double lambda = PARAMS.LM.LAMBDA;

        TE_var resids;
        E.getResiduals(resids);
        double dErr = resids.squaredNorm();

        const int MAX_ITERS = PARAMS.MAX_ITERS;
        for (int nIter = 0; nIter < MAX_ITERS; nIter++) {
            TJ J;
            E.getJ(J);

            double dNewErr = 1e+10;
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
            if (dNewErr < NUM_RESIDUALS * 1e-16) // Must be < 1e-8 so that these points actually end up as inliers, changing to e-6 actually makes almost no difference to speed
            {
                C3x3MatModel & M = static_cast<C3x3MatModel &> (models.addModel());
                double dScale = 1; //Keep signs the same for faster comparison
                TE Emodel = E.model();
                
                makeClosestE(Emodel);
                double dResid = getResids(Emodel, am0, am1);
                if(dResid > 0.01)
                    bFalseMin = true;
                else
                {
                    if (IS_DEBUG) {
                        Eigen::JacobiSVD<TE> svd(Emodel);
                        if( svd.singularValues()(0) > 1.2*svd.singularValues()(1) || svd.singularValues()(2) > 0.1)
                            cout << "BAD MODEL FOUND: " << svd.singularValues().transpose() << endl;
                    }

                    //cout << dResid << endl;

                    if (Emodel(0, 0) < 0) dScale = -1;
                    for (int r = 0; r < 3; r++)
                        for (int c = 0; c < 3; c++)
                            M(r, c) = dScale * Emodel(r, c);

                    if (PARAMS.VERBOSE) cout << "Converged in " << nIter << " iterations\n";
                    //cout << Emodel << endl;
                    //cout << "p" << E.getParams().transpose() << endl << endl;
                    
                    stats.add(nIter, CGDStats::eSuccess);
                    rememberParams(E.getParams(), aParamsAlreadyFound);
                
                    return true;
                }
            }
            
            //If it doesn't look like it will hit 0 then stop: dNewErr / (dErr-dNewErr) is num of steps required at this rate of descent
            if (bFalseMin || dNewErr > (dErr - dNewErr) * PARAMS.DROPOUT_EARLY_ON_FALSE_MINIMA/* this is about 1*/ * (MAX_ITERS - nIter)) {
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

int findE_LMforBasisR3(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, TEModels & models, CGDStats & stats) {
    Matrix<double, 9, 4 > EE;
    makeBasisForE(anHypSet, m0, m1, EE);
    TParamsAlreadyFound aParamsAlreadyFound;
    
    //cout << "#########################\n";

    Eigen::Vector3d am0[5];
    Eigen::RowVector3d am1[5];
    for(int i=0; i<5;i++)
    {
        am0[i] = Eigen::Vector3d(m1[anHypSet[i]].getX(), m1[anHypSet[i]].getY(), 1);
        am1[i] = Eigen::RowVector3d(m0[anHypSet[i]].getX(), m0[anHypSet[i]].getY(), 1);
    }

    for (int i = 0; i < PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS; i++) {
        CLevMarETE_R3<false>::calcEssentialMat_LM_Basis(PARAMS, EE, models, stats, am0, am1, aParamsAlreadyFound);

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


/*int CEssentialMatLMbasis::getModels_int(const TSubSet & anHypSet, CModels & models) {

    static CEssentialMatGradientDescParams PARAMS(0, 0);
    REPEAT(1, { PARAMS.VERBOSE = IS_DEBUG;    PARAMS.LM.LAMBDA = 0.0002;    PARAMS.MAX_NUM_SUCCESSFUL_ATTEMPTS = 8;    PARAMS.DIST_TO_EXISTING_SOLN_THRESH = 0.1;});
        
    CGDStats noStats;
    TEModels & modelsCast = static_cast<TEModels &>(models);
    return findE_LMforBasis(PARAMS, anHypSet, p2, p1, modelsCast, noStats);
}*/
