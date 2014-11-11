#include "mvNormalSampler.h"

#include "util/random.h"
#include "util/convert.h"
#include <Eigen/QR>
#include <Eigen/Dense>

template<class MatrixType>
void checkIsPosDef(const MatrixType & information)
{
    typedef Eigen::Matrix<double, MatrixType::RowsAtCompileTime, 1> VecType;
    
    //cout << information.sum() << "=info sum" << endl;
    //std::cerr << information << "=info/covariance\n" << endl;
    CVALGRIND_CHECK_INIT::checkInit(information.sum());

    if(information.squaredNorm() > 0)//TB June 2013: turned check always on
    {
        if((information - information.transpose()).squaredNorm() > 1e-6)
        {
            cout << information << "=info/covariance\n" << endl;
            THROW("Covariance matrix not symmetric");
        }
        
        //Eigen::EigenSolver<MatrixType> eigenVV(information);
        Eigen::SelfAdjointEigenSolver <MatrixType> eigenVV(information);
        const VecType & eigvals = eigenVV.eigenvalues(); 
        if(eigvals.minCoeff() < 0)
        {
            cout << "Eigenvalues: " << eigvals.transpose() << endl;
            CHECK(eigvals.minCoeff() <= 0, "Covariance matrix not +ve definite");
        }
        //const VecType & eigvals_imag = eigenVV.eigenvalues().imag(); // both complex
        //CHECK(!zero(eigvals_imag.maxCoeff()) || !zero(eigvals_imag.minCoeff()), "Covariance matrix not +ve definite (imag EVals)");
    }
}


template<class MatrixType>
MatrixType matrixSqrt(const MatrixType & cov)
{
    typedef Eigen::Matrix<double, MatrixType::RowsAtCompileTime, 1> VecType;

    checkIsPosDef(cov);

    //First find sqrt of covariance
    Eigen::SelfAdjointEigenSolver<MatrixType> eigenVV(cov);
    const VecType & eigvals = eigenVV.eigenvalues(); 
    const MatrixType & eigVectors = eigenVV.eigenvectors(); 

    MatrixType Diag; Diag.resize(cov.rows(), cov.cols()); Diag.setZero();
    MatrixType Diag_test = Diag;

    for(int i=0;i<cov.rows();i++)
    {
        Diag(i,i) = sqrt(eigvals(i));
        Diag_test(i,i) = eigvals(i);
    }

    MatrixType covarianceSqrt = eigVectors * Diag * eigVectors.transpose();

#ifdef _DEBUG
    MatrixType covTest = eigVectors * Diag_test * eigVectors.transpose();

    try {
        //TODO: check ok...
        CHECK(!zero(0.01*(covTest - cov).squaredNorm()), "Matrix eigendecomposition failed");
        //TODO: restore CHECK(!zero(0.01*(covSqrt*covSqrt - cov).squaredNorm()), "Matrix sqrt failed");
    }
    catch(...)
    {
        cout << "eigvals: " << eigvals.transpose() << endl;
        cout << "eigVectors: " << eigVectors << endl;
        cout << "Diag: " << Diag << endl;
        cout << "Diag_test: " << Diag_test << endl;
        cout << "covSqrt: " << covarianceSqrt << endl;
        cout << "covTest: " << covTest << endl;
        cout << "cov: " << cov << endl;
        throw;
    }
#endif
    if(IS_DEBUG) CHECK(std::isnan(covarianceSqrt.sum()), "covSqrt is nan");
    return covarianceSqrt;
}


template<class MatrixType>
CMVNormalSampler<MatrixType>::CMVNormalSampler(const VecType & mean, const MatrixType & cov) : mean(mean), mvNormal(mean.size())
{
    const bool bVerbose = false;
    if(bVerbose)
        cout << "About to sample from MV distn with mean: " << mean.transpose() << " cov: \n" << cov << endl;
    
    CVALGRIND_CHECK_INIT::checkInit(cov.sum());
    CVALGRIND_CHECK_INIT::checkInit(mean.sum());

    if(IS_DEBUG) CHECK(cov.squaredNorm() < 1e-32, "cov is 0"); //May be tiny when projecting accurate 3d points into image
    if(IS_DEBUG) CHECK(std::isnan(mean.sum()), "mean is nan");
    if(IS_DEBUG) CHECK(std::isnan(cov.sum()), "cov is nan");

    covSqrt = matrixSqrt<MatrixType>(cov);
}

template<class MatrixType>
typename CMVNormalSampler<MatrixType>::VecType CMVNormalSampler<MatrixType>::sample()
{
    for(int i=0; i<mean.size(); i++)
        mvNormal(i) = CRandom::Normal();
    
    //cout << "mvNormal: " << (mvNormal).transpose() << endl;
    //cout << "Resid: " << (covSqrt*mvNormal).transpose() << endl;
    if(IS_DEBUG && mean.size() > 10)
    {
        double dVar = mvNormal.squaredNorm()/mean.size();
        CHECK(dVar < 0.1 || dVar > 10, "Error computing normal RVs")
    }        

    VecType sample = mean + covSqrt*mvNormal;
    if(sample.squaredNorm() == 0)
    {
        THROW("Sample is 0 (extremely unlikely)");
    }
    
    return sample;
}

/*
void makeClosestPosDef(Eigen::MatrixXd & P)
{
    Eigen::MatrixXd Pt = 0.5*(P + P.transpose());
    //TODO should use Higham 2002
    / *
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Pt);
    
    //By spectral theorem, want Pt = U D U^T
    Eigen::MatrixXd D(P.rows(), P.rows());
    D.setZero();
    D.diagonal() = svd.singularValues();
    
    P = svd.matrixU() * D * svd.matrixU().transpose(); //TODO use V as well* /
    
    Eigen::RealSchur<Eigen::MatrixXd> schur(Pt);
    Eigen::MatrixXd T = schur.matrixT();
    //cout << T << endl;
    for(int r=0;r<T.rows();r++)
        for(int c=0;c<T.rows();c++)
            if(r!=c) 
                T(r,c) = 0;
            else
                T(r,c) = max<double>(T(r,c), 1e-6);
        
    P = schur.matrixU() * T * schur.matrixU().transpose();
    
    cout << P - Pt << endl;
    
    checkIsPosDef(P);
}*/


template class CMVNormalSampler<Eigen::Matrix2d>;
template class CMVNormalSampler<Eigen::Matrix3d>;
template class CMVNormalSampler<Eigen::MatrixXd>;

template void checkIsPosDef<>(const Eigen::Matrix2d & information);
template void checkIsPosDef<>(const Eigen::Matrix3d & information);

//TB: export matrix sqrt properly
template Eigen::Matrix2d matrixSqrt<>(const Eigen::Matrix2d & cov);
template Eigen::Matrix3d matrixSqrt<>(const Eigen::Matrix3d & cov);
template Eigen::MatrixXd matrixSqrt<>(const Eigen::MatrixXd & cov);

