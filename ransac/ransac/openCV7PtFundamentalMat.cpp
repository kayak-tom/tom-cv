/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCV8PtEssentialMatModified.cpp
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#include "openCV7PtFundamentalMat.h"
#include "util/opencv.h"
#include "ransac.h"

#ifndef USE_OLD_OPENCV
#include "opencv2/calib3d/calib3d.hpp"
#endif

using namespace std;

int COpenCV7PtFundamentalMat::getModels_int(const TSubSet & anHypSet, CModels & pModels)
{

	//CvMat H = cvMat(3, 3, CV_64FC1, const_cast<double *>(model.asDouble9()));

	double adPoints1_data[2*7];
	double adPoints2_data[2*7];

	cv::Mat CvMat1(nPoints, 2, CV_64FC1, adPoints1_data);
	cv::Mat CvMat2(nPoints, 2, CV_64FC1, adPoints2_data);

	for(int i=0; i<nPoints; i++)
	{
		CvMat1.at<double>(i, 0) = p1[anHypSet[i]].getX();
		CvMat1.at<double>(i, 1) = p1[anHypSet[i]].getY();
		CvMat2.at<double>(i, 0) = p2[anHypSet[i]].getX();
		CvMat2.at<double>(i, 1) = p2[anHypSet[i]].getY();
	}
	//vector<uchar> status;

	cv::Mat allSolutionsMat = cv::findFundamentalMat( CvMat1, CvMat2, cv::FM_7POINT);

	for(int offset = 0; offset < allSolutionsMat.rows/3; offset++)
	{
		double dL2Norm = 0;
		C3x3MatModel & model = dynamic_cast<C3x3MatModel &>(pModels.addModel());
		for(int r=0; r<3; r++)
			for(int c=0; c<3; c++)
			{
				model(r, c) = allSolutionsMat.at<double>(r+offset, c);
				dL2Norm += sqr(model(r, c));
			}
		dL2Norm = sqrt(dL2Norm);
		double dLenInv = M_SQRT2/dL2Norm;

		for(int r=0; r<3; r++)
			for(int c=0; c<3; c++)
				model(r, c) = dLenInv * model(r, c);
	}

	return pModels.numModels();
}

/*
bool COpenCVHomography::fitModel_int(CMask & mask, CModel & pModel, bool bVerbose)
{
    int nInliers = mask.countInliers();
	if(nInliers < 4) return false;

    const int nPoints = p1.size();

	CvMat H = cvMat(3, 3, CV_64FC1, const_cast<double *>(pModel.asDouble9()));

	CvMat * pCvMat1 = cvCreateMat(nInliers, 2, CV_64FC1);
	CvMat * pCvMat2 = cvCreateMat(nInliers, 2, CV_64FC1);

    int nInlier = 0;
    for(int i=0; i<nPoints; i++)
	{
        if(mask[i])
        {
		    cvmSet(pCvMat1, nInlier, 0, p1[i].getX());
		    cvmSet(pCvMat1, nInlier, 1, p1[i].getY());
		    cvmSet(pCvMat2, nInlier, 0, p2[i].getX());
		    cvmSet(pCvMat2, nInlier, 1, p2[i].getY());
            nInlier++;
        }
	}
    if(IS_DEBUG) CHECK(nInlier != nInliers, "Counting error");

	cvFindHomography( pCvMat1, pCvMat2,
	                       &H, / *int method=* / 0,
	                       / * ransacReprojThreshold=* / 0, 0);

    cvReleaseMat(&pCvMat1);
	cvReleaseMat(&pCvMat2);

	//cvReleaseMat(&pCvMask);

	return true;
}
*/
