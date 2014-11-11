#include "drawEllipse.h"
#include "util/convert.h"


//From #include "geom/mvNormalSampler.h"
template<class MatrixType>
MatrixType matrixSqrt(const MatrixType & cov);

void drawEllipse(cv::Mat & im, const Eigen::Matrix2d & A, const Eigen::Vector2d & b, const double sigma, const cv::Scalar col, const int thickness, const bool bMarkCentre)
{
    if(bMarkCentre)
    {
        const cv::Point centre(cvRound(b.x()), cvRound(b.y()));
        cv::circle(im, centre, 2, col, -1);
    }
    
    const int nNumPoints = 100;
    
    const Eigen::Matrix2d rootA = matrixSqrt(A);
    
    cv::Point p_last;
    for(int i=0;i<nNumPoints;i++)
    {
        const double dTheta = i*(2*M_PI/nNumPoints);
        const Eigen::Vector2d p1 (sin(dTheta), cos(dTheta));
        
        const Eigen::Vector2d p_sigma = p1 * sigma;
        
        const Eigen::Vector2d p_trans = rootA*p_sigma + b;
        const cv::Point p_im(cvRound(p_trans.x()), cvRound(p_trans.y()));
        
        if(i>0)
            cv::line(im, p_last, p_im, col, thickness);
        
        p_last = p_im;
    }
}