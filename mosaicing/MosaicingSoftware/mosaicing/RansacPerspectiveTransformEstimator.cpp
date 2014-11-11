#include "util/exception.h"
pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "RansacPerspectiveTransformEstimator.h"
#include "util/opencv_highgui.h"

namespace grc {

Transform * RansacTransformEstimator::getTransform(const CorrespondenceSet * corresp) const
{
        Transform * result = Transform::newTransform();
		Transform * testTransform = Transform::newTransform();

        const double RANSAC_INLIER_MAX_DIST_PX = testTransform->getRANSACThreshhold();
        const size_t RANSAC_SAMPLE_SET = testTransform->getRANSACSetSize();

        if (corresp->size() < RANSAC_SAMPLE_SET) {
			cvSetIdentity(result);
			return result;
		}

		int numInliers = 0;
		int bestNumInliers = 0;
		bool *inliers = new bool[corresp->size()];
		bool *bestInliers = new bool[corresp->size()];

		CvMat *sample1 = cvCreateMat(3,RANSAC_SAMPLE_SET,CV_64FC1);
		CvMat *sample2 = cvCreateMat(3,RANSAC_SAMPLE_SET,CV_64FC1);

		cvSet(sample1, cvScalarAll(1));
		cvSet(sample2, cvScalarAll(1));


		CvMat *all1 = cvCreateMat(3, corresp->size(), CV_64FC1);
		CvMat *all2 = cvCreateMat(3, corresp->size(), CV_64FC1);
		CvMat *pred = cvCreateMat(3, corresp->size(), CV_64FC1);

		for (size_t p = 0; p < corresp->size(); p++) {
			cvmSet(all1, 0, p, (*corresp)[p].pointInIm1().x);		
			cvmSet(all1, 1, p, (*corresp)[p].pointInIm1().y);		
			cvmSet(all1, 2, p, 1);		
			cvmSet(all2, 0, p, (*corresp)[p].pointInIm2().x);		
			cvmSet(all2, 1, p, (*corresp)[p].pointInIm2().y);		
			cvmSet(all2, 2, p, 1);	
		}

		int * indecies = new int[RANSAC_SAMPLE_SET];
		int numIters = maxRansacIters_;

		for (int i = 0; i < numIters; i++) {
			// Pick four points
			for (size_t j = 0; j < RANSAC_SAMPLE_SET; j++) {
				int ix = rand() % corresp->size();
				bool newIx = true;
				for (size_t k = 0; k < j; k++) {
					if (ix == indecies[k]) newIx = false;
				}
				if (newIx) {
					indecies[j] = ix;
					cvmSet(sample1, 0, j, (*corresp)[ix].pointInIm1().x );
					cvmSet(sample1, 1, j, (*corresp)[ix].pointInIm1().y );
					cvmSet(sample2, 0, j, (*corresp)[ix].pointInIm2().x );
					cvmSet(sample2, 1, j, (*corresp)[ix].pointInIm2().y );
				} else {
					j--;
				}
			}
//			cvFindHomography(sample1, sample2, testTransform);
            testTransform->estimateFromPoints(sample1, sample2);

/*
			if (fabs(cvmGet(testTransform, 0, 2)) > 20) {
				cout << "Transform: " << endl;
				for (int r = 0; r < 3; r++) {
					cout << " ";
					for (int c = 0; c < 3; c++) {
						cout << cvmGet(testTransform, r, c) << " ";
					}
					cout << endl;
				}
				cout << "Computed from " << endl;


			}
*/
			numInliers = 0;
            
            testTransform->applyToPoints(all1, pred);

			for (size_t f = 0; f < corresp->size(); f++) {
                //Lots of potential speedups here: write cvmFastGet, div!
				double dx = cvmGet(pred, 0, f) - cvmGet(all2, 0, f);
				double dy = cvmGet(pred, 1, f) - cvmGet(all2, 1, f);
				inliers[f] = (dx*dx + dy*dy) < (RANSAC_INLIER_MAX_DIST_PX*RANSAC_INLIER_MAX_DIST_PX);
				if (inliers[f]) {
					numInliers++;
				}
			}

			if (numInliers > bestNumInliers) {
				bestNumInliers = numInliers;
				for (size_t f = 0; f < corresp->size(); f++) {
					bestInliers[f] = inliers[f];
				}
				double ratio = double(numInliers) / double(corresp->size());
				ratio = ratio*ratio; // Squared
				ratio = ratio*ratio; // To the fourth, which is what we need
				numIters = int(log(0.01) / log(1-ratio));
			}
		}

		/*if (bestNumInliers < 5) {
			IplImage *img = cvCreateImage(cvSize(320,240), IPL_DEPTH_8U, 3);
			cvSetZero(img);
			for (size_t p = 0; p < corresp->size(); p++) {
				CvPoint p1 = cvPoint(int((*corresp)[p].first.x+0.5), int((*corresp)[p].first.y+0.5));
				CvPoint p2 = cvPoint(int((*corresp)[p].second.x+0.5), int((*corresp)[p].second.y+0.5));
				cvCircle(img, p1, 3, CV_RGB(255,0,0), -1);
				cvCircle(img, p2, 3, CV_RGB(0,255,0), -1);
				cvLine(img, p1, p2, CV_RGB(0, 0, 255));
			}
			cvSaveImage("debug.jpg", img);
			cvReleaseImage(&img);
			result = NULL;
		} */

        //Numbers (approx) + formula come from Brown and Lowe
        const int MIN_INLIERS = 6; //was 8
        const double MIN_PROP_INLIERS = 0.3;

        const int minInliers = MIN_INLIERS + (int)( MIN_PROP_INLIERS * (double)corresp->size() );
        
    	//cout << "RANSAC found " << bestNumInliers << " (of " << corresp->size() << ") inliers in " << numIters << " iterations: " << endl;
        if(bestNumInliers < minInliers)
        {
            result = NULL;
        } else {

			CvMat *p1 = cvCreateMat(3, bestNumInliers, CV_64FC1);
			CvMat *p2 = cvCreateMat(3, bestNumInliers, CV_64FC1);
			
			cvSet(p1, cvScalarAll(1));
			cvSet(p2, cvScalarAll(1));
			
			int ix = 0;
			for (size_t f = 0; f < corresp->size(); f++) {
				if (bestInliers[f]) {
					cvmSet(p1, 0, ix, cvmGet(all1, 0, f) );
					cvmSet(p1, 1, ix, cvmGet(all1, 1, f) );
					cvmSet(p2, 0, ix, cvmGet(all2, 0, f) );
					cvmSet(p2, 1, ix, cvmGet(all2, 1, f) );
					ix++;
				}
			}
			//cvFindHomography(p1, p2, result);
            result->estimateFromPoints(p1, p2);


	/*	DEBUGGING code commented out for now, remove once we're sure this is working
			if ( (fabs(cvmGet(result, 0, 0) - 1) > 0.1) ||
				 (fabs(cvmGet(result, 0, 1) - 0) > 0.1) ||
				 (fabs(cvmGet(result, 0, 2) - 0) > 20) ||
				 (fabs(cvmGet(result, 1, 0) - 0) > 0.1) ||
				 (fabs(cvmGet(result, 1, 1) - 1) > 0.1) ||
				 (fabs(cvmGet(result, 1, 2) - 0) > 20) ||
				 (fabs(cvmGet(result, 2, 0) - 0) > 0.1) ||
				 (fabs(cvmGet(result, 2, 1) - 0) > 0.1) ||
				 (fabs(cvmGet(result, 2, 2) - 1) > 0.1) ) {
				for (int r = 0; r < 3; r++) {
					for (int c = 0; c < 3; c++) {
						cout << cvmGet(result, r, c) << " ";
					}
					cout << endl;
				}	
				cout << endl;

				cout << "P1: " << endl;
				for (int r= 0; r < 2; r++) {
					for (int c = 0; c < p1->cols; c++) {
						cout << cvmGet(p1, r, c) << " ";
					}
					cout << endl;
				}
				cout << endl;
				
				cout << "P2: " << endl;
				for (int r= 0; r < 2; r++) {
					for (int c = 0; c < p2->cols; c++) {
						cout << cvmGet(p2, r, c) << " ";
					}
					cout << endl;
				}
				cout << endl;

				
				IplImage *img = cvCreateImage(cvSize(320,240), IPL_DEPTH_8U, 3);
				cvSetZero(img);
				for (size_t p = 0; p < corresp->size(); p++) {
					CvPoint p1 = cvPoint(int((*corresp)[p].first.x+0.5), int((*corresp)[p].first.y+0.5));
					CvPoint p2 = cvPoint(int((*corresp)[p].second.x+0.5), int((*corresp)[p].second.y+0.5));
					if (bestInliers[p]) {
						cvCircle(img, p1, 3, CV_RGB(255,0,0), -1);
						cvCircle(img, p2, 3, CV_RGB(0,255,0), -1);
					} else {
						cvCircle(img, p1, 3, CV_RGB(255,0,0));
						cvCircle(img, p2, 3, CV_RGB(0,255,0));
					}
					cvLine(img, p1, p2, CV_RGB(255, 255, 255));
				}
				cvSaveImage("debug.jpg", img);
				cvReleaseImage(&img);
			}
		*/
		}
		cvReleaseMat(&pred);
		cvReleaseMat(&all1);
		cvReleaseMat(&all2);
		cvReleaseMat(&sample1);
		cvReleaseMat(&sample2);
		delete testTransform;
		
		delete [] inliers;
		delete [] bestInliers;
        delete [] indecies;

		return result;
	}

}
pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
