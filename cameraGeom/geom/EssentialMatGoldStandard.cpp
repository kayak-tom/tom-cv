/*
 * EssentialMatGoldStandard.cpp
 *
 *  Created on: 21/11/2009
 *      Author: hilandtom
 */

#include "geom.h"
#include "geom_eigen.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/LU>

#include "util/dynArray.h"
#include <algorithm>

using namespace std;
//extern double edSVDCondition, edMin, edNormalisedMin, edRatio, edSVDMin, edJTJMin, edJTJNormalisedMin, edJTJRatio, edJDet;

void makeE(const Eigen::Matrix3d & R, const Eigen::Vector3d & cameraDir, Eigen::Matrix3d & EAtParams)
{
    Eigen::Matrix3d X_mat;
    Xmat(cameraDir, X_mat);
    EAtParams = X_mat * R;
}

void makeE(const C3dRotation & rotation, const C3dPoint & cameraDir, Eigen::Matrix3d & EAtParams)
{
    Eigen::Vector3d t;
    t << cameraDir.getX(), cameraDir.getY(), cameraDir.getZ();
    makeE(rotation, t, EAtParams);
}

void makeE(const C3dRotation & rotation, const Eigen::Vector3d & cameraDir, Eigen::Matrix3d & EAtParams)
{
    Eigen::Matrix3d rot;
    rotation.asMat(rot);
    makeE(rot, cameraDir, EAtParams);
}
static const int NUM_PARAMS = 7; //New GCC gives an error if this is defined in the class?!

template<bool bAllResids> // if bAllResids get full N-d residual vector
class CLevMarForE
{
    static const int NUM_MODEL_VARS = bAllResids ? Eigen::Dynamic : 1;
    typedef Eigen::Matrix<double, NUM_PARAMS, 1> TParamVector;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS,1> TResidVector;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS, NUM_PARAMS> TJMatrix;
    typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;

    static void getParamVector(const C3dRotation & rotation, const C3dPoint & cameraDir, TParamVector & params)
    {
        for(int i=0;i<4;i++)
            params[i] = rotation[i];

        params[4] = cameraDir.getX();
        params[5] = cameraDir.getY();
        params[6] = cameraDir.getZ();
    }

    static void setParams(const TParamVector & params, C3dRotation & rotation, C3dPoint & cameraDir)
    {
        rotation = C3dRotation(params[0], params[1], params[2], params[3]); //Normalisation normally redundent, which is ok
        cameraDir = C3dPoint(params[4], params[5], params[6]);
        cameraDir.normalise(); //Normalisation normally redundent, which is ok
    }

    /*static void cost(const Eigen::Matrix3d & E, const CPointVec2d & p1, const CPointVec2d & p2, TResidVector & vResiduals)
    {
        const int NUM_POINTS = (int)p1.size();

        vResiduals(0) = 0;

        for(int i=0; i<NUM_POINTS; i++)
        {
            Eigen::Vector3d v1; p1[i].asVector(v1);
            Eigen::Vector3d v2; p2[i].asVector(v2);
//            const Eigen::Matrix<double, 1, 1> & cost = v2.transpose() * E * v1;
            Eigen::RowVector3d v2_E_normalised = v2.transpose() * E; / *cout << v2_E_normalised.squaredNorm() << endl;* / v2_E_normalised.normalize();
            const Eigen::Matrix<double, 1, 1> & cost = v2_E_normalised * v1;

            if(bAllResids)
                vResiduals(i) = cost(0, 0);
            else
                vResiduals(0) += sqr(cost(0, 0)); //Todo: What norm??
        }
    }*/
    static inline double robustCost(const double d, const CRobustLMforEParams & robustRefinementParams)
    {
        const double fd = fabs(d);
        const double THRESH = robustRefinementParams.dRobustifyThresh;
        if(fd < THRESH)
            return d;

        const double absRobustRes = (fd - THRESH) * robustRefinementParams.dRobustifyScale + THRESH;

        return (d>0) ? absRobustRes : -absRobustRes;
    }

    template<bool bRobustCost>
    static void cost(const Eigen::Matrix3d & E, const CPointVec2d & p1, const CPointVec2d & p2, const CRobustLMforEParams & robustRefinementParams, TResidVector & vResiduals)
    {
        const int NUM_POINTS = (int)p1.size();

        vResiduals(0) = 0;

        for(int i=0; i<NUM_POINTS; i++) //use the same normalised score as OpenCV
        {
            Eigen::Vector3d v1; p1[i].asVector(v1);
            Eigen::Vector3d v2; p2[i].asVector(v2);

            Eigen::Vector3d abc = E*v1;
            double s1 = sqr(abc(0)) + sqr(abc(1));
            const Eigen::Matrix<double, 1, 1> d1 = v2.transpose() * abc;
            //double d1_squared = d1(0) * d1(0);

//            double dResid = d1_squared / s1;
            double dResid = d1(0) / sqrt(s1);

            abc = v2.transpose()*E;

            double s0 = sqr(abc(0)) + sqr(abc(1));
            const Eigen::Matrix<double, 1, 1> d0 = v1.transpose() * abc;
            //double d0_squared = d0(0) * d0(0);

//            dResid += d0_squared / s0;
            dResid += d0(0) / sqrt(s0);

            if(bAllResids)
                vResiduals(i) = bRobustCost ? robustCost(dResid, robustRefinementParams) : dResid;
            else
                vResiduals(0) += sqr(dResid); //Todo: What norm??
        }
    }

    static void calcResiduals(const bool bRefineE, const CPointVec2d & p1, const CPointVec2d & p2, const C3dRotation & rotation, const C3dPoint & cameraDir, const CRobustLMforEParams & robustRefinementParams, TResidVector & vResiduals)
    {
        if(bRefineE)
        {
            Eigen::Matrix3d EAtParams;
            makeE(rotation, cameraDir, EAtParams);
            //cout << EAtParams << "\n=E at params\n";
            if(robustRefinementParams.bRobustErr)
                cost<true>(EAtParams, p1, p2, robustRefinementParams, vResiduals); //Set residuals
            else
                cost<false>(EAtParams, p1, p2, robustRefinementParams, vResiduals); //Set residuals

            /* yes these are the same and are been minimised!
            CMask m(p1.size());
            for(int i=0; i<m.size(); i++)
                m[i] = 1;

            double dSSR = E_sumSquaredResiduals(EAtParams, p1, p2, m);
            cout << dSSR  << ", " << vResiduals.squaredNorm() << "=cost, E_sumSquaredResiduals\n";
            */
        }
        else
        {
            const CCamera P, Pp = rotation | cameraDir; 
            camReprojectionErr(P, Pp, p1, p2, (Eigen::VectorXd &)vResiduals, 0);
        }
    }

    /*static void calcResiduals(const CPointVec2d & p1, const CPointVec2d & p2, const C3dRotation & rotation, const C3dPoint & cameraDir, TResidVector & vResiduals)
    {
        CCamera Pp = rotation | cameraDir;
        makeE(rotation, cameraDir, EAtParams);
        //cout << EAtParams << "\n=E at params\n";
        cost(EAtParams, p1, p2, vResiduals); //Set residuals

        / * yes these are the same and are been minimised!
        CMask m(p1.size());
        for(int i=0; i<m.size(); i++)
            m[i] = 1;

        double dSSR = E_sumSquaredResiduals(EAtParams, p1, p2, m);
        cout << dSSR  << ", " << vResiduals.squaredNorm() << "=cost, E_sumSquaredResiduals\n";
        * /
    }*/

    static void getJAtParams(const bool bRefineE, const C3dRotation & rotation, const C3dPoint & cameraDir, const CPointVec2d & p1, const CPointVec2d & p2, const double epsilon, const CRobustLMforEParams & robustRefinementParams, TJMatrix & J, TResidVector & vResiduals)
    {
        const int RESIDS = (int)vResiduals.rows();
        TParamVector paramVec, tempParamVec;
        getParamVector(rotation, cameraDir, paramVec);

        //Now vary each param and find deriv.
        TResidVector vResidualsAtEps(RESIDS);
        for(int nParam = 0; nParam < NUM_PARAMS; nParam++)
        {
            TParamVector tempParamVec = paramVec; tempParamVec(nParam) += epsilon;
            C3dRotation tempRotation; C3dPoint tempcameraDir;
            setParams(tempParamVec, tempRotation, tempcameraDir);

            calcResiduals(bRefineE, p1, p2, tempRotation, tempcameraDir, robustRefinementParams, vResidualsAtEps);

            J.col(nParam) = (vResidualsAtEps - vResiduals)/epsilon;
        }
    }

public:
    /* Use Levenberg-1944 as Marquardt-1963 a little less stable when matrix is near singular. Probably makes little difference (Lampton-1997)
     * Damping strategy from \cite{Lampton-1997}. BOOST and DROP factors appear to make little difference in performance, neither does Marquardt;s method
    */    //void fromRotMat(const Matrix & Rot);

    //Adjust 10 parameters to refine 9 matrix elements. Need to use Levenberg algorithm because Lev-Mar assumes derivatives nonzero somewhere (todo: I'm sure there's a proper way of coping with this)
    static double refineE(const bool bRefineE, const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, const CRobustLMforEParams & robustRefinementParams, bool bVerbose)
    {
        if(IS_DEBUG) CHECK(p1.size() < 5, "Too few points for R and T refinement, five needed (more better)");
        if(IS_DEBUG) CHECK(p1.size() != p2.size(), "Too few points for R and T refinement, five needed (more better)");

        const int MAX_ITERS = 300, NUM_POINTS=p1.size(), 
            RESIDS = bAllResids ? (bRefineE ? NUM_POINTS: NUM_POINTS*4 ) : 1;

        TParamVector params;

        double lambda = 0.1;
        const double eps=sqr(0.000001), BOOST=5, DROP=0.3;
        int nIter=0;

        TResidVector vResiduals(RESIDS), vNewResids(RESIDS);
        TJMatrix J(RESIDS, NUM_PARAMS);

        calcResiduals(bRefineE, p1, p2, rotation, cameraDir, robustRefinementParams, vResiduals);
        double dErr = vResiduals.squaredNorm();
        TJTJMatrix JTJ;

        for(; nIter < MAX_ITERS && dErr > eps && dErr <= HUGE /*breakout on overflow. Needed??*/; nIter++)
        {
            const double epsilon = 0.00005;

            getJAtParams(bRefineE, rotation, cameraDir, p1, p2, epsilon, robustRefinementParams, J, vResiduals); //vResiduals MUST ALREADY BE FILLED

            //Repeated from getJAndResidAtParams but doesn't matter...
            getParamVector(rotation, cameraDir, params);

            double dNewErr = HUGE;
            TParamVector newParams; newParams.setZero();


            for(;;)
            {
                JTJ = J.transpose()*J;

                //cout << JTJ << endl;
                //JTJ.diagonal() *= (1+lambda);
                JTJ.diagonal().array() += lambda;

                const TJTJMatrix inv = JTJ.inverse();

                TParamVector paramUpdateVec = inv*J.transpose()*vResiduals;

                newParams = params - paramUpdateVec;

                C3dRotation newRotation;
                C3dPoint newCameraDir;
                setParams(newParams, newRotation, newCameraDir);

                calcResiduals(bRefineE, p1, p2, newRotation, newCameraDir, robustRefinementParams, vNewResids);

                dNewErr = vNewResids.squaredNorm(); //*DOES NOT REALLY MATTER that this has been squared twice for 1d cost version */

                //cout << "New Residual = " << dNewErr << endl;
                if(dNewErr <= dErr)
                {
                    rotation = newRotation;
                    cameraDir = newCameraDir;
                    break;
                }
                else
                {
                    lambda *= BOOST;
                    double dUpdateSize = paramUpdateVec.squaredNorm();
                    if(dUpdateSize < cube(epsilon) )
                        break; //We're close enough to the minimum. breakout WITHOUT adopting newRotation, newCameraDir
                }
            }

            //Err is proportional to number of points (RESIDS)
            if(dErr-dNewErr < RESIDS * 1e-9) //Includes dNewErr>dErr, i.e. breakout from previous loop
                break;

            vResiduals = vNewResids;
            lambda *= DROP;
            dErr = dNewErr;
            params = newParams;
            //setParams(params, rotation, cameraDir);
        }
        return 1;

        /*const TParamVector & jtjDiag = JTJ.diagonal();
        const TParamVector & curvature = JTJ.inverse().diagonal();
        double dMinCurvature = fabs(jtjDiag.minCoeff());// /sqrt(jtjDiag.squaredNorm())); // fabs(curvature.minCoeff()) / sqrt(curvature.squaredNorm());
        const bool bIllConditioned = dMinCurvature < 1e-5;

        if(nIter==0)
        {
            cout << "No refinement needed\n";
            return 1;
        }

        if(bVerbose)
        {
            /// *if(!bIllConditioned)
            {
                cout << "JTJ qx qy qz qw x y z (" << jtjDiag.minCoeff() << "): " ;
                cout << jtjDiag.transpose() << endl;
                cout << "Inv JTJ qx qy qz qw x y z (" << curvature.minCoeff() << "): " ;
                cout << curvature.transpose() << endl;
            }// * /

            if(bRefineE)
                cout << "Essential matrix";
            else
                cout << "Camera";

            cout << " refined after " << nIter << " iterations" << endl;
            cout << "Error = " << dErr << endl;
            cout << "t = " << cameraDir << endl;
            cout << "R = " << rotation << endl;

            CDynArray<double> aResids;
            for(int i=0;i<RESIDS;i++)
                aResids.push_back(vResiduals(i));
            std::sort(aResids.begin(), aResids.end(), sortResids);
            aResids.pp(",");
        }

        if(!bIllConditioned) //Todo: Remove outlier and continue
        {
            double dRelErr = dErr/NUM_POINTS;
            if(bVerbose)
                cout << "Relative err = " << dRelErr << "...";
            if(dRelErr > robustRefinementParams.dHighResidCutoff / *0.00002* /) //Todo: Work out better threshhold
            {
                double dCondition = robustRefinementParams.dHighResidConditionScale * (dRelErr / robustRefinementParams.dHighResidCutoff);
                if(bVerbose)
                    cout << "Residuals too high, ill-conditioned (" << dCondition << ")\n";
                return dCondition;
            }
            else
            {
                //cout << "err less than thresh = " << robustRefinementParams.dHighResidCutoff << endl;
            }

        }

        / *if(bIllConditioned)
        {
            cout << "Refinement ill-conditioned\n";
            cout << "JTJ qx qy qz qw x y z (" << jtjDiag.minCoeff() << "): " ;
            cout << jtjDiag.transpose() << endl;
            cout << "Inv JTJ qx qy qz qw x y z (" << curvature.minCoeff() << "): " ;
            cout << curvature.transpose() << endl;

//            return false;
        }
        else* / if (!bIllConditioned && bVerbose)
        {
#if EIGEN_VERSION_AT_LEAST(2,90,0)
            Eigen::JacobiSVD<TJTJMatrix> svd( JTJ );
#else
            Eigen::SVD<TJTJMatrix> svd( JTJ );
#endif
            cout << svd.singularValues()(6) << "=min SV\n";
            cout << svd.singularValues()(0) / svd.singularValues()(6) << "=SV condition\n";
        }

        //Very approx model of condition: if worse than 1e-5 probably ill conditioned. Not necessarily true though
        double dHeuristicScaleVar = 1;
        
        if(dMinCurvature < 1e-12)
            return 100;

        for(double dCurvatureThresh = 1e-5; dMinCurvature <= dCurvatureThresh; dCurvatureThresh*=0.1, dHeuristicScaleVar*=2);

        return dHeuristicScaleVar;



        //Set edSVDCondition, edMin, edNormalisedMin, edRatio, edJTJMin, edJTJNormalisedMin, edJTJRatio, edSVDMin;
/ *
        Eigen::SVD<TJTJMatrix> svd( JTJ );
        edSVDMin = svd.singularValues()(6);
        edSVDCondition = svd.singularValues()(0) / edSVDMin;


        edMin = curvature.minCoeff();
        edRatio = edMin / curvature.maxCoeff();
        edNormalisedMin = edMin / sqrt(curvature.squaredNorm());

        edJTJMin = jtjDiag.minCoeff();
        edJTJRatio = edJTJMin / jtjDiag.maxCoeff();
        edJTJNormalisedMin = edJTJMin / sqrt(jtjDiag.squaredNorm());
        edJDet = JTJ.determinant();

        return 1.0/dMinCurvature;*/
    }
private:
    inline static bool sortResids(const double d1, const double d2)
    {
        return fabs(d1)<fabs(d2);
    };
};

double refineRT_LM_SampsonError(const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, const bool bAllResiduals, bool bVerbose)
{
    const CRobustLMforEParams robustRefinementParams;
    if(bAllResiduals)
        return CLevMarForE<true>::refineE(true, p1,p2,rotation, cameraDir, robustRefinementParams, bVerbose);
    else
        return CLevMarForE<false>::refineE(true, p1,p2,rotation,cameraDir, robustRefinementParams, /*0.00004, 5,*/ bVerbose);
}

double refineRT_LM_SampsonError_robust(const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, const CRobustLMforEParams & robustRefinementParams, bool bVerbose)
{
    return CLevMarForE<true>::refineE(true, p1,p2,rotation, cameraDir, robustRefinementParams, bVerbose);
}

double EssentialMat_SSD(const Eigen::Matrix3d & E1, const Eigen::Matrix3d & E2, bool bVerbose)
{
    Eigen::Matrix3d E1norm = E1/(E1.norm());
    Eigen::Matrix3d E2norm = E2/(E2.norm());

    double dDiff = std::min<double>((E1norm-E2norm).squaredNorm(), (E1norm+E2norm).squaredNorm());

    if(bVerbose && dDiff > 0.01)
    {
        cout << "Matrices not equal: SSD=" << dDiff << endl;
        cout << "E1=" << endl << E1 << endl;
        cout << "E2 normalised=" << endl << E2 << endl;
    }

    return dDiff;
}

void refineRT_ReprojErr_NotPoints(const CPointVec2d & p1, const CPointVec2d & p2, C3dRotation & rotation, C3dPoint & cameraDir, bool bVerbose)
{
    const CRobustLMforEParams robustRefinementParams;
    CLevMarForE<true>::refineE(false, p1,p2,rotation,cameraDir, robustRefinementParams, bVerbose);
}

//template<bool bAllResids> const int NUM_PARAMS = 7;
//template<bool bAllResids> const int CLevMarForE<true>::NUM_PARAMS = 7;
