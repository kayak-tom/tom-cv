/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCV8PtEssentialMatModified.cpp
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#include "openCV4PtHomography.h"
#include "util/opencv.h"
#include "ransac.h"

#ifndef USE_OLD_OPENCV
#include "opencv2/calib3d/calib3d.hpp"
#endif
//#include "opencv/cv.h"//cvFindHomography has been moved here in the latest OpenCV

using namespace std;

int COpenCV4ptHomography::getModels_int(const TSubSet & anHypSet, CModels & pModels)
{
	C3x3MatModel & model = dynamic_cast<C3x3MatModel &>(pModels.addModel());

	//CvMat H = cvMat(3, 3, CV_64FC1, const_cast<double *>(model.asDouble9()));

	//double adPoints1_data[2*4];
	//double adPoints2_data[2*4];

	//CvMat CvMat1 = cvMat(nPoints, 2, CV_64FC1, adPoints1_data);
	//CvMat CvMat2 = cvMat(nPoints, 2, CV_64FC1, adPoints2_data);
	cv::Mat CvMat1(nPoints, 2, CV_64FC1);
	cv::Mat CvMat2(nPoints, 2, CV_64FC1);

	for(int i=0; i<nPoints; i++)
	{
        CvMat1.at<double>(i,0) = p1[anHypSet[i]].getX();
        CvMat1.at<double>(i,1) = p1[anHypSet[i]].getY();
        CvMat2.at<double>(i,0) = p2[anHypSet[i]].getX();
        CvMat2.at<double>(i,1) = p2[anHypSet[i]].getY();
        
        
		/*cvmSet(&CvMat1, i, 0, p1[anHypSet[i]].getX());
		cvmSet(&CvMat1, i, 1, p1[anHypSet[i]].getY());
		cvmSet(&CvMat2, i, 0, p2[anHypSet[i]].getX());
		cvmSet(&CvMat2, i, 1, p2[anHypSet[i]].getY());*/
	}

	cv::Mat temp = cv::findHomography(CvMat1, CvMat2);
    
    for(int r=0;r<3;r++)
        for(int c=0;c<3;c++)
            model(r,c) = temp.at<double>(r,c);
    

	return pModels.numModels();
}

bool COpenCVHomography::fitModel_int(CMask & mask, CModel & pModel, bool bVerbose)
{
	C3x3MatModel & model = dynamic_cast<C3x3MatModel &>(pModel);
    int nInliers = mask.countInliers();
	if(nInliers < 4) return false;

    const int nPoints = p1.size();

	//cv::Mat H(3, 3, CV_64FC1, const_cast<double *>(model.asDouble9()));

	cv::Mat pCvMat1(nInliers, 2, CV_64FC1);
	cv::Mat pCvMat2(nInliers, 2, CV_64FC1);

    int nInlier = 0;
    for(int i=0; i<nPoints; i++)
	{
        if(mask[i])
        {
            pCvMat1.at<double>(nInlier, 0) = p1[i].getX();
            pCvMat1.at<double>(nInlier, 1) = p1[i].getY();
            pCvMat2.at<double>(nInlier, 0) = p2[i].getX();
            pCvMat2.at<double>(nInlier, 1) = p2[i].getY();
		    
            /*cvmSet(pCvMat1, nInlier, 0, p1[i].getX());
		    cvmSet(pCvMat1, nInlier, 1, p1[i].getY());
		    cvmSet(pCvMat2, nInlier, 0, p2[i].getX());
		    cvmSet(pCvMat2, nInlier, 1, p2[i].getY());*/
            nInlier++;
        }
	}
    if(IS_DEBUG) CHECK(nInlier != nInliers, "Counting error");

	cv::Mat temp = cv::findHomography( pCvMat1, pCvMat2);
                 
    for(int r=0;r<3;r++)
        for(int c=0;c<3;c++)
            model(r,c) = temp.at<double>(r,c);

    //cvReleaseMat(&pCvMat1);
	//cvReleaseMat(&pCvMat2);

	return true;
}
