#include "util/exception.h"

pragma_warning(push)
pragma_warning(disable : 4996) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "BaySACTransformEstimator.h"
#include "util/opencv_highgui.h"
#include "ransac/ransac.h"
#include "image/convert_OpenCV.h"
#include <boost/filesystem.hpp>
#include <sstream>
#include "util/location.h"

#include <boost/math/distributions/normal.hpp>
#include "cvGeom.h" //may conflict...

#include<boost/scoped_array.hpp>

//#ifndef USE_OLD_OPENCV
//#include "opencv2/opencv.hpp"
//#endif


namespace grc {

    //d^2 is uniformly distributed

    inline double unifSqPdf(double d, double dMaxSq_inv) {
        if(IS_DEBUG) CHECK(d < 0 || dMaxSq_inv <= 0, "Bad params");
        //return 2*d*dMaxSq_inv;
        return dMaxSq_inv; //independent of d
    }
    //d is uniformly distributed. NOT CORRECT PRIOR

    inline double unifPdf(double d, double dMax_inv) {
        if(IS_DEBUG) CHECK(d < 0 || dMax_inv <= 0, "Bad params");
        return dMax_inv;
    }

    inline double diam_unifAreaPdf(double d, double dMax_inv) {
        if(IS_DEBUG) CHECK(d < 0 || dMax_inv <= 0, "Bad params");
        return 2 * d*dMax_inv;
    }

    inline double halfNormalPdf(double d, double s) {
        if(IS_DEBUG) CHECK(d < 0 || s <= 0, "Bad params");
        boost::math::normal_distribution<> N01(0, s);
        /*const double dUnifFrac = 0.2;
        return (1-dUnifFrac) * 2*boost::math::pdf(N01, d)
                        + dUnifFrac  * unifSqPdf(d, 1.0/sqr(500));*/
        return 2 * boost::math::pdf(N01, d);
    }

    inline double normalPdf2d(double dx, double dy, double s) {
        if(IS_DEBUG) CHECK(s <= 0, "Bad params");
        boost::math::normal_distribution<> N01(0, s);
        /*const double dUnifFrac = 0.2;
        return (1-dUnifFrac) * 2*boost::math::pdf(N01, d)
                        + dUnifFrac  * unifSqPdf(d, 1.0/sqr(500));*/
        return boost::math::pdf(N01, dx) * boost::math::pdf(N01, dy); //can be simplified...
    }

    Transform * BaySACTransformEstimator::getTransform(const CBoWCorrespondences * corresp, size_t imageId1, size_t imageId2, const Transform * lastTrans, double dLastTransScale) const {
        C3x3MatModel H;

        T2dPoints calibratedPoints1, calibratedPoints2;
        CInlierProbs adArrLikelihood, adMMLikelihood;
        CPointIdentifiers pointIds;

        /* Test verifies that updates are independent of the order they are applied in, and appear not to violate law of total prob.
        if(lastTrans)
        {
                CBoWCorrespondences * pCorr = const_cast<CBoWCorrespondences *>(corresp);
                pCorr->clear();
                for(int i=0; i<3; i++)
                        pCorr->push_back(CCorrespondence(CLocation(1,1), CLocation(i,1), 0.2));
                for(int i=0; i<2; i++)
                        pCorr->push_back(CCorrespondence(CLocation(2,i), CLocation(0,1), 0.2));
        }*/

        corresp->calibrate(K, calibratedPoints1, calibratedPoints2, adArrLikelihood, pointIds, true);

        const int N = corresp->size();

        CMask inliers(N);

        if (lastTrans) {
            //Bayes update on inlier likelihood
            const double D_MAX = 500; //TODO: image width
            const double dMax_inv = 1.0 / (D_MAX);
            const double dMaxSq_inv = sqr(dMax_inv);

            typedef CFixedArray<int, NN_MAX> TFastIntVec;
            typedef CFixedArray<int, 2 * (NN_MAX - 1) > TLongerFastIntVec; // a N-M match can be incompatible with (N-1) + (M-1) others
            typedef map2<int, TFastIntVec, std::less<int> /*, CIndividualPool_NoFree_Allocator<std::pair<int, TFastIntVec>, 512 >*/ > TPointCollisionMaps;
            TPointCollisionMaps samePointLeft, samePointRight;

            std::vector<TLongerFastIntVec> aaIncompatable(N);

            //Todo: copied from BaySAC init
            ARRAY(TFastIntVec *, avVectorsLeft, N); //Remember locations to prevent a set lookup
            ARRAY(TFastIntVec *, avVectorsRight, N);

            for (int i = 0; i < N; i++) {
                //add i to a list for its left and right points
                TFastIntVec & vLeft = samePointLeft.initOrGet(pointIds[i].id1());
                vLeft.push_back(i);
                TFastIntVec & vRight = samePointRight.initOrGet(pointIds[i].id2());
                vRight.push_back(i);

                //Convert references to pointers...
                //Prob not std c++...
                avVectorsLeft[i] = &vLeft;
                avVectorsRight[i] = &vRight;
            }
            for (int i = 0; i < N; i++) {
                TLongerFastIntVec & aIncompat = aaIncompatable[i];
                TFastIntVec & aIncompatLeft = *(avVectorsLeft[i]),
                        & aIncompatRight = *(avVectorsRight[i]);

                for (TFastIntVec::const_iterator pLeft = aIncompatLeft.begin(); pLeft != aIncompatLeft.end(); pLeft++) {
                    int nLeft = *pLeft;
                    if (nLeft != i)
                        aIncompat.push_back(nLeft);
                }
                for (TFastIntVec::const_iterator pRight = aIncompatRight.begin(); pRight
                        != aIncompatRight.end(); pRight++) {
                    int nRight = *pRight;
                    if (nRight != i)
                        aIncompat.push_back(nRight);
                }
            }

            //dLastTransScale = 2 if this time gap is twice as far, etc.
            const double s = 150 * sqrt(dLastTransScale); //, s_inv = 1.0/s;
            const bool bVerbose = false;

            for (int i = 0; i < N; i++) {
                const CCorrespondence & corr = corresp->get(i);
                CLocation l1 = corr.Location1();
                CvPoint2D64f p1 = cvPoint2D64f(l1.dx(), l1.dy());

                CLocation l2 = corr.Location2();
                CvPoint2D64f p2 = cvPoint2D64f(l2.dx(), l2.dy());

                CvPoint2D64f pp = lastTrans->applyToPoint(p1);

                double d = sqrt(SSD(p2, pp)) / dLastTransScale;

                const double P_inlier = adArrLikelihood[i];
                //const double dP_inlier_num = halfNormalPdf(d, s) * P_inlier;
                const double dP_inlier_num = normalPdf2d(p2.x - pp.x, p2.y - pp.y, s) * P_inlier;
                //const double dP_inlier_num_test = halfNormalPdf(d, s) * P_inlier;
                const double dP_inlier_denom = dP_inlier_num + (1 - P_inlier) * unifSqPdf(d, dMaxSq_inv);
                adArrLikelihood[i] = (dP_inlier_num / dP_inlier_denom);
                //adArrLikelihood[i] *= (adArrLikelihood[i] > P_inlier) ? 0.99 : 1.01; //0.99 prevents probs summing to just over 1 emerging

                if (bVerbose) {
                    cout << "Prior prob of d=" << P_inlier << endl;
                    cout << "d=" << d << endl;
                    cout << "Posterior prob of d=" << adArrLikelihood[i] << endl;
                }
                if(IS_DEBUG) CHECK(adArrLikelihood[i] >= 1 || adArrLikelihood[i] <= 0, "Point prob update failed");

                TLongerFastIntVec & aIncompat = aaIncompatable[i];
                if (aIncompat.size() > 0) {
                    const double dUpdate = (1 - adArrLikelihood[i]) / (1 - P_inlier);
                    //const double dUpdate = 1/(1+adArrLikelihood[i]-P_inlier); //Todo : derive... Venn diagram works
                    //const double dUpdate = unifSqPdf(d, dMaxSq_inv)/dP_inlier_denom;
                    //const double dUpdate = unifSqPdf(d, dMaxSq_inv)/(adArrLikelihood[i]*halfNormalPdf(d, s)   +   (1 - adArrLikelihood[i])*unifPdf(d, dMax_inv));

                    for (TLongerFastIntVec::const_iterator pIncompat = aIncompat.begin(); pIncompat != aIncompat.end(); pIncompat++) {
                        const int j = *pIncompat;
                        const double P_j_inlier = adArrLikelihood[j];
                        if (bVerbose)
                            cout << "Prior prob of incompat. j: " << adArrLikelihood[j] << "...";

                        adArrLikelihood[j] *= dUpdate;

                        if (bVerbose)
                            cout << "posterior prob of j:" << adArrLikelihood[j] << endl;

                        if(IS_DEBUG) CHECK(adArrLikelihood[j] >= 1 || adArrLikelihood[j] <= 0, "Point prob update failed (buggy--set BOW.BOWCorrespondenceProb.MAX_PRIOR=0.4 to fix)");

                        if (false && dUpdate > 1) //And update (drop) points incompatible with j (only bother if we're boosing probs, necessary so that probs never sum to more than 1)
                        {
                            const double dUpdate_j = (1 - adArrLikelihood[j]) / (1 - P_j_inlier);
                            TLongerFastIntVec & aIncompat_j = aaIncompatable[j];
                            for (TLongerFastIntVec::const_iterator pIncompat_j = aIncompat_j.begin(); pIncompat_j != aIncompat_j.end(); pIncompat_j++) {
                                const int k = *pIncompat_j;

                                if (bVerbose)
                                    cout << "Prior prob of incompat. k: " << adArrLikelihood[k] << "...";

                                adArrLikelihood[k] *= dUpdate_j;
                                if(IS_DEBUG) CHECK(dUpdate_j > 0 && dUpdate_j < 1 /* because only dropping probs here */ && adArrLikelihood[k] >= 1 || adArrLikelihood[k] <= 0, "Point prob update failed");

                                if (bVerbose)
                                    cout << "posterior prob of k:" << adArrLikelihood[k] << endl;
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < N; i++) {
                adArrLikelihood[i] *= 0.99; //Stop overly likely points
                if (bVerbose) {
                    const CCorrespondence & corr = corresp->get(i);
                    CLocation l1 = corr.Location1();
                    CvPoint2D64f p1 = cvPoint2D64f(l1.dx(), l1.dy());

                    CLocation l2 = corr.Location2();
                    CvPoint2D64f p2 = cvPoint2D64f(l2.dx(), l2.dy());

                    CvPoint2D64f pp = lastTrans->applyToPoint(p1);

                    double d = sqrt(SSD(p2, pp)) / dLastTransScale;
                    if(adArrLikelihood[i] >= 1 || adArrLikelihood[i] <= 0)
                    {
                        cout << i << ": d: " << d << " PP: " << adArrLikelihood[i] << " incompat count: " << aaIncompatable[i].size() << endl;
                        THROW("Point prob update failed");
                    }
                }
            }
        }

        int nInliers = getH(calibratedPoints1, calibratedPoints2, adArrLikelihood, pointIds, PARAMS, H, inliers, K.focalLength(), 1);

        if (imSource) {
            CvPtr<IplImage> im1(imSource->createImage());
            CvPtr<IplImage> im2(imSource->createImage());
            imSource->loadImage(imageId1, im1);
            imSource->loadImage(imageId2, im2);
            IplImage * pImOut = 0;
            markCorrespondences(corresp, inliers, im1, im2, &pImOut);
            ostringstream filename;
            boost::filesystem::create_directories("corresp");
            filename << "corresp/corresp" << imageId1 << ',' << imageId2 << ".jpg";
            cvSaveImage(filename.str().c_str(), pImOut);
            cvReleaseImage(&pImOut);
        }

        int nMinInliers = getMinInliers(adArrLikelihood, 0.2 /*20% overlap*/, 4 /*points for homography*/);
        if (nInliers > nMinInliers) {
            const_cast<BaySACTransformEstimator *> (this)->nRansacSuccess++;
            const_cast<BaySACTransformEstimator *> (this)->nTotalInliers += nInliers;

            CvMatFixed < 3, 3 > cvK, cvKinv, Hp, temp;
            Transform * result = new PerspectiveTransform(); // Transform::newTransform();

            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++) {
                    cvmSet(cvK, r, c, K(r, c));
                    cvmSet(cvKinv, r, c, K.inv(r, c));
                    cvmSet(Hp, r, c, H(r, c));
                }
            cvMatMul(cvK, Hp, temp);
            cvMatMul(temp, cvKinv, (CvMat *) (*result));

            //H <- K*H*Kinv;
            /*for(int r=0;r<3;r++)
                    for(int c=0;c<3;c++)
                    {
                            //Compute H(r,c)

                            double dEl = 0;
                            for(int i=0;i<3;i++)
                            {
                                    dEl += K(r, i) * H(i, c);
                            }

                            cvmSet((CvMat *)(*result), r, c, dEl);
                    }

            / *for(int r=0;r<3;r++)
                    for(int c=0;c<3;c++)
                            cvmSet((CvMat *)(*result), r, c, H(r, c));*/

            return result;
        } else {
            const_cast<BaySACTransformEstimator *> (this)->nRansacFail++;
            cout << "Too few inliers (" << nInliers << ") found, min " << nMinInliers << endl;
            return 0;
        }
    }

    BaySACTransformEstimator::~BaySACTransformEstimator() {
        cout << "RANSAC Successes: " << nRansacSuccess << " Fail: " << nRansacFail << " av inliers on success: " << nTotalInliers / (double) nRansacSuccess << endl;
    }

}

pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
