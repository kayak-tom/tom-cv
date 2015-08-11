#include <util/exception.h>
#include <util/pp.h>
#include "reconstruct3d.h"
#include <geom/levMarNumerical.h>
#include <Eigen/Dense>

template<typename TFloat, int N, int M>
void pseudoInv(const Eigen::Matrix<TFloat, N, M> & A, Eigen::Matrix<TFloat, M, N> & A_dag)
{
    Eigen::JacobiSVD<Eigen::Matrix<TFloat, N, M> > svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    Eigen::Matrix<TFloat, M, N> D_inv; D_inv.setZero();
    D_inv.diagonal() = svd.singularValues().array().inverse();
    A_dag = svd.matrixV() * D_inv * svd.matrixU().transpose();
    
    cout << A*A_dag << endl;
}

optional<C3dWorldPoint> reconstruct(const CWorldCamera &P, const CWorldCamera &Pp, const C2dImagePointPx &p, const C2dImagePointPx &pp, const double dReprojectionErrorPx, int bVerbose) {

    typedef Eigen::Matrix<double, 1, 4, Eigen::RowMajor + Eigen::AutoAlign> Vec4d;
    
    if(bVerbose>1)
    {
        COUT(p);
        COUT(pp);
    }    

    const Vec4d Prow1 = P.getPMat().row(0);
    const Vec4d Prow2 = P.getPMat().row(1);
    const Vec4d Prow3 = P.getPMat().row(2);
    const Vec4d Pprow1 = Pp.getPMat().row(0);
    const Vec4d Pprow2 = Pp.getPMat().row(1);
    const Vec4d Pprow3 = Pp.getPMat().row(2);
    
    Eigen::Matrix4d A;
    A.row(0) = p.x() * Prow3 - Prow1;
    A.row(1) = p.y() * Prow3 - Prow2;
    A.row(2) = pp.x() * Pprow3 - Pprow1;
    A.row(3) = pp.y() * Pprow3 - Pprow2;
    
    if(bVerbose>1)
    {
        COUT(Pp.getPMat().row(2));
        COUT(Pprow3);
        COUT(A);
    }

    Eigen::JacobiSVD< Eigen::Matrix4d > svdA(A, Eigen::ComputeFullV);// speed-up by not computing U

    const TEigen3dPointHomog V = svdA.matrixV().col(3);
    
    if(bVerbose>1)
    {
        COUT(svdA.matrixV());
        COUT(V.transpose());
    }
    
    //C3dWorldPoint p3d=(V.head<3>()/V(3)).eval();
    C3dWorldPoint p3d(V);
    
    COUT(p3d.transpose());
    COUT(C3dWorldPoint(V).transpose());
    
    optional<const C2dImagePointPx> x_reproj = P.projectToPx(p3d);
    if(!x_reproj)
    {
        COUT2("Failed to project back to 2D", p3d);
        return optional<C3dWorldPoint>();
    }
    optional<const C2dImagePointPx> xp_reproj = Pp.projectToPx(p3d);
    if(!xp_reproj)
    {
        COUT2("Failed to project back to 2D", p3d);
        return optional<C3dWorldPoint>();
    }
    C2dImagePointPx err1 = *x_reproj - p;
    //Equal to err1 C2dImagePointPx err2 = *xp_reproj - pp;
    const double dErrThresh = sqr(dReprojectionErrorPx);
    if(err1.squaredNorm() > dErrThresh/* || err2.squaredNorm() > dErrThresh*/)
    {
        COUT2("RE too high", err1);
        //COUT2("RE too high", err2);
        return optional<C3dWorldPoint>();
    }
    
    if(IS_DEBUG)
    {
        C2dImagePointPx err2 = *xp_reproj - pp;
        CHECK((err1.squaredNorm() > 1) && !within(err1.squaredNorm(), err2.squaredNorm(), 0.25), "Point isn't central (numerical errors are enough to make it a bit out...)");
    }
    
    //REPEAT(5, test3DRecon(P,Pp));
    
    return optional<C3dWorldPoint>(p3d);
}

class CLMForX : public CLMFunction
{
    const CWorldCamera & P;
    const CWorldCamera & Pp;
    const C2dImagePointPx & p;
    const C2dImagePointPx & pp;

public:
    CLMForX(const CWorldCamera &P, const CWorldCamera &Pp, const C2dImagePointPx &p, const C2dImagePointPx &pp) 
    : P(P) , Pp(Pp), p(p), pp(pp)
    {
        
    }
    
    virtual Eigen::VectorXd init()
    {
        Eigen::VectorXd Xinit = P.pxToWorld_z(p, 1.0);// + Pp.pxToWorld(pp, 0.5);
        return Xinit;
    }
    
    virtual int inputs() const { return 3; }
    virtual int values() const { return 4; }
    virtual eLMSuccessStatus function(const Eigen::VectorXd &x, Eigen::VectorXd &resids, bool bVerbose = false, const int nParamChanged = -1)
    {
        C3dWorldPoint X = x;
        
        const optional<const C2dImagePointPx> im_x = P.projectToPx(X); 
        const optional<const C2dImagePointPx> im_xp = Pp.projectToPx(X);
        
        if(!im_x || !im_xp)
            return eLMFail;
            
        C2dImagePointPx err1 = *im_x - p;
        C2dImagePointPx err2 = *im_xp - pp;
        
        resids.segment<2>(0) = err1*0.1; //0.1 keeps values sensible
        resids.segment<2>(2) = err2*0.1;
        return eLMSuccess;
    }
};

optional<const C3dWorldPoint> reconstructLM(const CWorldCamera &P, const CWorldCamera &Pp, const C2dImagePointPx &p, const C2dImagePointPx &pp) {
    CLMForX lmForX(P, Pp, p, pp);
    CLevMar LM(lmForX);
    Eigen::VectorXd Xinit = lmForX.init();
    LM.minimise(Xinit);
    
    if(LM.residuals().norm() > 12)
        return optional<const C3dWorldPoint>();
    
    TEigen3dPoint Xtemp = Xinit;
    C3dWorldPoint X = Xtemp;
    return optional<const C3dWorldPoint>(X);
}

/* More tests in autotest */
void test3DRecon(const CWorldCamera & P1, const CWorldCamera & P2)
{
    const bool bVerbose = false;
    C2dImagePointPx x(250,250);
    C3dWorldPoint X = P1.pxToWorld_depth(x, 0.5);
    C2dImagePointPx xp = P2.projectToPx_checked(X);
    //C2dImagePointPx x_recon = P1.projectToPx(X);
    
    optional<C3dWorldPoint> X_recon1 = reconstruct(P1,P2,x,xp, 10, bVerbose);
    
    if(bVerbose) cout << X.transpose() << " = " << X_recon1->transpose() << endl;

    optional<const C3dWorldPoint> X_recon2 = reconstructLM(P1,P2,x,xp);
    if(bVerbose) cout << X.transpose() << " = " << X_recon2->transpose() << endl;
    
    CHECK( !zero((X-*X_recon1).squaredNorm()), "2 point 3D reconstruction failed");
    CHECK( !zero((X-*X_recon2).squaredNorm()), "2 point 3D reconstruction LM failed");
}

