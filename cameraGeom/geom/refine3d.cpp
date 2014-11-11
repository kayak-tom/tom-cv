/*
 * refine3d
 *
 *  Created on: 21/11/2009
 *      Author: hilandtom
 */

#include "geom.h"
#include "geom_eigen.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>
//#include <Eigen/Sparse>

using namespace std;

class CLevMarRefine3d
{
    static const int NUM_PARAMS = Eigen::Dynamic, NUM_MODEL_VARS = Eigen::Dynamic;
    typedef Eigen::Matrix<double, NUM_PARAMS, 1> TParamVector;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS,1> TResidVector;

//    typedef Eigen::SparseMatrix<double/*, NUM_MODEL_VARS, NUM_PARAMS*/> TJMatrix;
//    typedef Eigen::SparseMatrix<double/*, NUM_PARAMS, NUM_PARAMS*/> TJTJMatrix;

	typedef Eigen::Matrix<double, NUM_MODEL_VARS, NUM_PARAMS> TJMatrix;
    typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;

    static void getParamVector(const C3dRotation & rotation, const C3dPoint & cameraDir, const T3dPoints & points3d, TParamVector & params)
    {
        for(int i=0;i<4;i++)
            params[i] = rotation[i];

        params[4] = cameraDir.getX();
        params[5] = cameraDir.getY();
        params[6] = cameraDir.getZ();

        int i=7;
        for(T3dPoints::const_iterator pPoint = points3d.begin(); pPoint != points3d.end(); pPoint++, i+=3)
        {
            params[i] = pPoint->getX();
            params[i+1] = pPoint->getY();
            params[i+2] = pPoint->getZ();
        }
    }

    static void setParams(const TParamVector & params, C3dRotation & rotation, C3dPoint & cameraDir, T3dPoints & points3d)
    {
        rotation = C3dRotation(params[0], params[1], params[2], params[3]); //Normalisation normally redundent, which is ok
        cameraDir = C3dPoint(params[4], params[5], params[6]);
        cameraDir.normalise(); //Normalisation normally redundent, which is ok

        int i=7;
        for(T3dPoints::iterator pPoint = points3d.begin(); pPoint != points3d.end(); pPoint++, i+=3)
        {
            *pPoint = C3dPoint(params[i], params[i+1], params[i+2]);
        }
    }

    static void cost(const CPointVec2d & p1, const CPointVec2d & p2, const C3dRotation & rotation, const C3dPoint & cameraDir, const T3dPoints & points3d, TResidVector & vResiduals)
    {
        const int NUM_POINTS = (int)p1.size();

        vResiduals(0) = 0;

        CCamera P, Pp = rotation | cameraDir;

        for(int i=0; i<NUM_POINTS; i++) //use the same normalised score as OpenCV
        {
            C2dPoint resid1 = p1[i] - P*points3d[i];
            C2dPoint resid2 = p2[i] - Pp*points3d[i];
            vResiduals[4*i+0] = resid1.getX();
            vResiduals[4*i+1] = resid1.getY();
            vResiduals[4*i+2] = resid2.getX();
            vResiduals[4*i+3] = resid2.getY();
        }
    }

    static void calcResiduals(const CPointVec2d & p1, const CPointVec2d & p2, const C3dRotation & rotation, const C3dPoint & cameraDir, const T3dPoints & points3d, TResidVector & vResiduals)
    {
        cost( p1, p2, rotation, cameraDir, points3d, vResiduals); //Set residuals
    }

    static void getJAtParams(const C3dRotation & rotation, const C3dPoint & cameraDir, const CPointVec2d & p1, const CPointVec2d & p2, const T3dPoints & points3d, const double epsilon, TJMatrix & J, TResidVector & vResiduals)
    {
        const int NUM_POINTS=(int)p1.size(), RESIDS = (int)vResiduals.rows(), NUM_PARAMS = 7 + NUM_POINTS*3;

        TParamVector paramVec(NUM_PARAMS), tempParamVec(NUM_PARAMS);
        getParamVector(rotation, cameraDir, points3d, paramVec);

        //Now vary each param and find deriv.
        TResidVector vResidualsAtEps(RESIDS);

        T3dPoints tempPoints3d(NUM_POINTS);

        for(int nParam = 0; nParam < NUM_PARAMS; nParam++)
        {
            tempParamVec = paramVec; tempParamVec(nParam) += epsilon;
            C3dRotation tempRotation;
            C3dPoint tempcameraDir;
            setParams(tempParamVec, tempRotation, tempcameraDir, tempPoints3d);

            calcResiduals(p1, p2, tempRotation, tempcameraDir, points3d, vResidualsAtEps);

            J.col(nParam) = (vResidualsAtEps - vResiduals)/epsilon;
        }
    }

public:
    /* Use Levenberg-1944 as Marquardt-1963 a little less stable when matrix is near singular. Probably makes little difference (Lampton-1997)
     * Damping strategy from \cite{Lampton-1997}. BOOST and DROP factors appear to make little difference in performance, neither does Marquardt;s method
    */
    //Adjust 10 parameters to refine 9 matrix elements. Need to use Levenberg algorithm because Lev-Mar assumes derivatives nonzero somewhere (todo: I'm sure there's a proper way of coping with this)
    static void refine3d(const CPointVec2d & p1, const CPointVec2d & p2, T3dPoints & points3d, C3dRotation & rotation, C3dPoint & cameraDir, bool bVerbose)
    {
        const int MAX_ITERS = 300, NUM_POINTS=p1.size(), 
            RESIDS = NUM_POINTS*4,
            PARAMS = 7 + NUM_POINTS*3;

        TParamVector params(PARAMS);
        points3d.clear(); points3d.resize(NUM_POINTS);

        double lambda = 0.1;
        const double eps=sqr(0.000001), BOOST=5, DROP=0.3;
        int nIter=0;

        TResidVector vResiduals(RESIDS), vNewResids(RESIDS);
        TJMatrix J(RESIDS, PARAMS);

        CCamera P, Pp = rotation | cameraDir;
        for(int i = 0; i < NUM_POINTS; i++)
        {
            points3d[i] = reconstruct(P, Pp, p1[i], p2[i]);
        }

        calcResiduals(p1, p2, rotation, cameraDir, points3d, vResiduals);
        double dErr = vResiduals.squaredNorm();

        for(; nIter < MAX_ITERS && dErr > eps && dErr <= HUGE /*breakout on overflow. Needed??*/; nIter++)
        {
            const double epsilon = 0.00005;

            getJAtParams(rotation, cameraDir, p1, p2, points3d, epsilon, J, vResiduals); //vResiduals MUST ALREADY BE FILLED

            //Repeated from getJAndResidAtParams but doesn't matter...
            getParamVector(rotation, cameraDir, points3d, params);

            double dNewErr = HUGE;
            TParamVector newParams(PARAMS); newParams.setZero();

            for(;;)
            {
                TJTJMatrix JTJ = J.transpose()*J;
                //cout << JTJ << endl;
                //JTJ.diagonal() *= (1+lambda);
                JTJ.diagonal().array() += lambda;
                /*for(int i=0; i<PARAMS; i++)
                    JTJ(i,i) += lambda;

                Eigen::SparseLU<TJTJMatrix> JTJ_LU;
                const TJTJMatrix inv = JTJ_LU*/

                const TJTJMatrix inv = JTJ.inverse();

                TParamVector paramUpdateVec = inv*J.transpose()*vResiduals;

                newParams = params - paramUpdateVec;

                C3dRotation newRotation;
                C3dPoint newCameraDir;
                T3dPoints newPoints3d(NUM_POINTS);

                setParams(newParams, newRotation, newCameraDir, newPoints3d);

                calcResiduals(p1, p2, newRotation, newCameraDir, newPoints3d, vNewResids);

                dNewErr = vNewResids.squaredNorm();

                //cout << "New Residual = " << dNewErr << endl;
                if(dNewErr <= dErr)
                {
                    rotation = newRotation;
                    cameraDir = newCameraDir;

                    for(int i = 0; i < NUM_POINTS; i++)
                        points3d[i] = newPoints3d[i];

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

        if(bVerbose)
        {
            cout << "Camera and 3d points";

            cout << " refined after " << nIter << " iterations" << endl;
            cout << "Error = " << dErr << endl;
            cout << "t = " << cameraDir << endl;
            cout << "R = " << rotation << endl;
        }
    }
};

void refine3d(const CPointVec2d & p1, const CPointVec2d & p2, T3dPoints & points3d, C3dRotation & rotation, C3dPoint & cameraDir, bool bVerbose)
{
    CLevMarRefine3d::refine3d(p1, p2, points3d, rotation, cameraDir, bVerbose);
}
