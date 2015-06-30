/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */
#define BOOST_NO_CXX11_SCOPED_ENUMS

#define DEPRECATED(x) THROW( "Calling deprecated code")

//#define SHOW_CORR
//#define TUNE "SlamConfig2"
#define LESS_THREADS IS_DEBUG //For easier debugging in gdb
//#define RUN100

static int TOTAL_RANSAC_INLIERS = 0; //For tuning
static double CLOSE_TO_CV = 0;
#ifndef _DEBUG
//#define _DEBUG
#endif

#include "Board.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "util/exception.h"
#include "util/opencv.h"

#include <algorithm>
#include <queue>
#include "time/SpeedTest.h"
#include "util/random.h"
#include "util/set2.h"
#include "imageSource/imageSourceFromDir.h"

#include <string.h>
#include <fstream>
#include <time.h>

#include "logging/redirectCout.h"

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <boost/smart_ptr.hpp>

#include "geom/geom.h"
#include "geom/geom_eigen.h"

#include "plot.h"

#include "autotest.h"
#include "geom/pointAlignment.h"
#include "cImageSimulator.h"

#include "threeDPointCollection.h"
#include "scale/BoWSpeedo.h"
#include "scale/EdgeScaleObserver.h"

#include "bow/bagOfWords.h"

#include "featureExtract/featureExtractor.h"

#include "ransac/ransac.h"

#include "bowslamParams.h"
#include "params/config.h"

#include "slamMap.h"
#include "svdScaleOptimiser.h"
#include <boost/math/distributions/normal.hpp>

#include "newgui.h"
#include "util/cout_TS.h"

#include "ransac/refineEOnRTManifold.h"
//#include "image/convert_OpenCV.h"

pragma_warning(push)
pragma_warning(disable : 4127) //const conditional expression

using namespace std;

typedef int TTime;

static const TTime FRAMERATE_MS = 1;

const int STACK_LIM = 10000; //safety check for infinite loops
#ifdef STACKCH
#define STACKCHECK(x) x
#else
#define STACKCHECK(x);
#endif

CvFont font, font_small;

enum eSCTypes {
    eGridSlide, eCubePhotos, ePlanarRot, eGridZoom, eSimCube
};
const eSCTypes eSCType = eSimCube;

const int nSyntheticTimesteps = 16; //8 or 16

CvMat * pointVectorToCvMat(const CPointVec2d & aPoints) {
    int n = (int) aPoints.size();
    CvMat * pMat = cvCreateMat(2, n, CV_64FC1);

    for (int c = 0; c < n; c++) {
        cvmSet(pMat, 0, c, aPoints[c].getX());
        cvmSet(pMat, 1, c, aPoints[c].getY());
    }

    return pMat;
}

void testRDCorrection(const CCamCalibMatrix & K, const IplImage * pFrame, const char * szFilename) {
    const int x = pFrame->width;
    const int y = pFrame->height;

    CvPtr<IplImage> greyIm(cvCreateImage(cvSize(x, y), IPL_DEPTH_8U, 1));
    CvPtr<IplImage> correctedIm(cvCreateImage(cvSize(x, y), IPL_DEPTH_8U, 1));

    if (pFrame->nChannels == 3)
        cvCvtColor(pFrame, greyIm, CV_RGB2GRAY);
    else
        cvCopy(pFrame, greyIm);

    /*for(int lx=20; lx < y; lx += 20)
            cvLine(greyIm, cvPoint( 0, lx ), cvPoint( x, lx ), CV_RGB(0,0,0));*/

    cvSetZero(correctedIm);

    CCamera P;
    P.calibrate(K);

    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++) {
            C2dPoint imPoint(i, j);
            imPoint.calibrate(K, true);

            C3dPoint p3_calib(imPoint, 1);
            C2dPoint p_backInIm = p3_calib.photo(P);

            int nNewX = cvRound(p_backInIm.getX());
            int nNewY = cvRound(p_backInIm.getY());
            if (nNewX >= 0 && nNewX < x && nNewY >= 0 && nNewY < y) {
                const uchar val = CIplPx<uchar>::getGrey(greyIm, i, j);
                CIplPx<uchar>::setGrey(correctedIm, nNewX, nNewY, val);
            }
        }
    cvSaveImage(szFilename, correctedIm);
}

CBoWCorrespondences * getSyntheticCorr(int T1, int T2, CRunGuiAp &gui, const CImParams & IM_PARAMS) {
    const int IM_WIDTH = IM_PARAMS.IM_WIDTH;
    const int IM_HEIGHT = IM_PARAMS.IM_HEIGHT;

    //N_FEATURES_RAND = 6; //hack for speedup
    CBoWCorrespondences * pCorr = new CBoWCorrespondences();
    if (T1 + 1 != T2 && T1 + 2 != T2) return pCorr;
    //	if(!(T1+1 == T2 && 4 < T2) && !(T1==0 && T2==4)) return pCorr; //test can cope without structure for a few frames

    const int nPerturb = 0;

    switch (eSCType) {
        case eGridSlide:
            for (int x = 100; x < 300; x += 40)
                for (int y = 100; y < 300; y += 40) {
                    CLocation L1(x + 25 * T1, y);
                    CLocation L2(x + 25 * T2, y + CRandom::perturb(nPerturb));
                    pCorr->push_back(CCorrespondence(L1, L2, 0.5));
                }
            break;
        case eGridZoom:
            for (int x = 1; x <= 16; x++)
                for (int y = 1; y <= 16; y++) {
                    CLocation L1(x * (20 + T1), y * (20 + T1));
                    CLocation L2(x * (20 + T2), y * (20 + T2) + CRandom::perturb(nPerturb));
                    pCorr->push_back(CCorrespondence(L1, L2, 0.5));
                }
            break;
        case eSimCube: //Todo: copy to AutoTest
            /*Matrix P(3,4), Pp(3,4);
             P=0; P(1,1)=P(2,2)=P(3,3) = 1;*/
            const double X_ROT_SCALE = -0.04, Y_ROT_SCALE = -0.5, Z_ROT_SCALE = 0.05;
            //const double X_ROT_SCALE = 0, Y_ROT_SCALE = 0, Z_ROT_SCALE= 0;

            int x[] = {0, 1, 2, 3, 3, 3, 3, 2, 0, 0, 0, 0, 1, 2, 3, 3};
            int z[] = {0, 0, 0, 0, 1, 2, 3, 3, 3, 2, 1, 0, 0, 0, 0, 1}; /*, bow, bAlignedStructure, *this, gui*/
            //int x[] = {0,0,0,0,0,0,0};
            //int z[] = {0,1,3,6,10,14,18};
            C3dPoint t1(x[T1] * 0.01, 0, z[T1] * -0.01);
            C3dPoint t2(x[T2] * 0.01, 0, z[T2] * -0.01);

            C3dPoint axis(X_ROT_SCALE, Y_ROT_SCALE, Z_ROT_SCALE);
            C3dRotation R1(axis, X_ROT_SCALE * T1); //= getRotMat(X_ROT_SCALE * T1, Y_ROT_SCALE * (T1 > 3 && T1 < 12), Z_ROT_SCALE * T1);
            C3dRotation R2(axis, X_ROT_SCALE * T2); //= getRotMat(X_ROT_SCALE * T2, Y_ROT_SCALE * (T2 > 3 && T2 < 12), Z_ROT_SCALE * T2);
            CCamera P = R1 | -(R1 * t1); //Todo: should Xmat(t*) be normalised?? I don't think it makes any difference
            CCamera Pp = R2 | -(R2 * t2);

            CCamCalibMatrix K;
            //getCalibrationMat(K, &adRDcoeffs, nRDcoeffs);
            P.calibrate(K);
            Pp.calibrate(K);

            //PRINTMAT(P);
            //PRINTMAT(Pp);
            const bool bClip = false;

            const double CUBE_RAD = bClip ? 0.2 : 0.1;
            const double CUBE_STEP = 0.04; //0.05
            const double CUBE_DEPTH = bClip ? 0.3 : 0.2;

            for (double z = CUBE_DEPTH - CUBE_RAD; z <= CUBE_DEPTH + CUBE_RAD; z += CUBE_STEP)
                for (double x = -CUBE_RAD; x <= CUBE_RAD; x += CUBE_STEP)
                    for (double y = -CUBE_RAD; y <= CUBE_RAD; y += CUBE_STEP) {
                        C3dPoint Q(x, y, z);
                        C2dPoint p = Q.photo(P);
                        C2dPoint pp = Q.photo(Pp);

                        CLocation L1(colVecToLoc(p));
                        CLocation L2(colVecToLoc(pp));
                        ///////////////// Check reconstruction is working:
                        /*PRINTMAT(p);
                         PRINTMAT(pp);
                         printLoc(L1);
                         printLoc(L2);

                         //Now calibrate image points and reconstruct:
                         p(2) = int(p(2));
                         p(1) = int(p(1));
                         pp(2) = int(pp(2));
                         pp(1) = int(pp(1));
                         p = K.i()*p;
                         pp = K.i()*pp;
                         Q=reconstruct(K.i()*P, K.i()*Pp, p, pp);
                         PRINTMAT(Q);*/
                        ////////////////////
                        if (bClip) if (CLIP(L1) || CLIP(L2)) continue; //Apply clipping frame

                        if (!Q.testInFront(P, Pp)) {
                            cout << "Synthetic point behind cam...";
                            continue;
                        }
                        //Equal for x-y planar motion
                        //cout << getDepth(P,Q) << ' ';T_ab_0
                        //cout << getDepth(Pp,Q) << "=depth from P Pp\n";

                        if (!testPair(P, Pp, p, pp)) {
                            cout << "Error: Point in front of cam, but reconstructed point behind cam\n";
                            //cout << (P);
                            //cout << (Pp);
                            cout << (reconstruct(P, Pp, p, pp));
                            continue;
                        }
                        /*if(!testPair(K.i()*P, K.i()*Pp, K.i()*p, K.i()*pp))
                         {
                         cout << "Error: Point in front of cam, but reconstructed calibrated point behind cam\n";
                         PRINTMAT(K.i()*P);
                         PRINTMAT(K.i()*Pp);
                         PRINTMAT(reconstruct(K.i()*P, K.i()*Pp, K.i()*p, K.i()*pp));
                         continue;
                         }
                         C3dPoint Qcalib=reconstruct(K.i()*P, K.i()*Pp, K.i()*p, K.i()*pp);*/

                        //Equal for x-y planar motion
                        //cout << getDepth(K.i()*P,Qcalib) << ' ';
                        //cout << getDepth(K.i()*Pp,Qcalib) << "=depth from  K.i()*P  K.i()*Pp\n";

                        ////////////////////
                        L1 = CLocation(L1.x() + CRandom::perturb(nPerturb), L1.y() + CRandom::perturb(nPerturb));

                        pCorr->push_back(CCorrespondence(L1, L2, 0.5));
                    }
            break;
    }
    if (gui.showImages()) {
        CvPtr<IplImage> img(cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, 3));
        cvSetZero(img);
        for (CBoWCorrespondences::iterator pCorrespondence = pCorr->begin(); pCorrespondence != pCorr->end(); pCorrespondence++) {
            //CCorrespondence * pCorr = *ppCorr;
            cvCircle(img, locToCvPoint(pCorrespondence->Location1()), 1, CV_RGB(0, 255, 0), 2);
            cvCircle(img, locToCvPoint(pCorrespondence->Location2()), 1, CV_RGB(0, 0, 255), 2);
        }
        gui.showIm(img, CRunGuiAp::eWindowSimIm);
    }

    return pCorr;
}

/*Could also: dynamically don't remove too many corr.
 remove corr with displacement vectors far from median
 */
void removeZeroCorrs(CBoWCorrespondences *& pCorr) {
    CBoWCorrespondences::iterator pGoodCorrespondence = pCorr->begin();
    for (CBoWCorrespondences::iterator pCorrespondence = pCorr->begin(); pCorrespondence < pCorr->end(); pCorrespondence++) {
        if (pCorrespondence->priorProb() > 0) {
            *pGoodCorrespondence = *pCorrespondence;
            pGoodCorrespondence++;
        }
    }
    int nCorrLeft = pGoodCorrespondence - pCorr->begin();
    pCorr->resize(nCorrLeft);

    DEBUGONLY(
    for (CBoWCorrespondences::iterator pCorrespondence = pCorr->begin(); pCorrespondence < pCorr->end(); pCorrespondence++) {
        if(IS_DEBUG) CHECK(pCorrespondence->priorProb() == 0, "All corrs should be init. now");
    })
}

void limitTrackLength(CBoWCorrespondences * pCorr, const CBoWSLAMParams & BOWSLAMPARAMS) {
    //if(!BOWSLAMPARAMS.StructureSize.REMOVE_DISTANT_POINTS) return;
    DEPRECATED((
            const int nMaxTrackLen = BOWSLAMPARAMS.MAX_TRACK_LEN;
            const int nTargetCorresps = BOWSLAMPARAMS.TARGET_CORRESPONDENCES;

    if ((int) (pCorr->size()) <= nTargetCorresps) return;

            CBoWCorrespondences::iterator ppEnd = pCorr->end();
            CDynArray<int> vTrackLengths; vTrackLengths.reserve(pCorr->size());
        for (CBoWCorrespondences::iterator pCorrespondence = pCorr->begin(); pCorrespondence < ppEnd; pCorrespondence++) {
            vTrackLengths.push_back(dist(pCorrespondence->Location1(), pCorrespondence->Location2()));
        }
    sort(vTrackLengths.begin(), vTrackLengths.end(), less<int> ());
            double nMaxTrackLenSq = max<double> (nMaxTrackLen * nMaxTrackLen, vTrackLengths[nTargetCorresps]);
            cout << (int) (sqrt(nMaxTrackLenSq)) << "=Max track length\n";
    for (CBoWCorrespondences::iterator pCorrespondence = pCorr->begin(); pCorrespondence < ppEnd; pCorrespondence++) {
        if (dist(pCorrespondence->Location1(), pCorrespondence->Location2()) > nMaxTrackLenSq) {
            pCorrespondence->zero();
        }
    }

    removeZeroCorrs(pCorr);
            cout << pCorr->size() << " correspondences remaining\n";
}
))
}

class CMxLockLater : private boost::noncopyable {
    boost::mutex & mx;
    boost::unique_lock<boost::mutex> * pLock;
public:

    CMxLockLater(boost::mutex & mx) : mx(mx), pLock(0) {
    };

    void lock() {
        if (!pLock)
            pLock = new boost::unique_lock<boost::mutex > (mx);
    }

    bool locked() const {
        return pLock != 0;
    }

    ~CMxLockLater() {
        delete pLock;
    }
};

double CSSObserver::s_dTypicalSS = 0;
int CSSObserver::s_nSSsObserved = 0;

class CSLAMLocMatch;
class CSLAMLocation;
class CBoWMap;

/*class CMapPos
{
public:
        enum ePositionAlg
        {
                ePosBA, ePosCameraMat, ePos3dAlignment, ePos1dAlignment, ePosExact, ePosExtrapolated
        }; //How did we recontruct this point?
        enum ePosType
        {
                eRelative, eAbsolute, eOrigin
        };

        class CMPScale
        {
                double dScale;
                CScaleDistn scaleDistn; //<= 0, sum to give total badness. Big is good
                int nInlierCount;
                const CMapPos * pRelMapPosForThisScale;
                double dStructureSize; //The reconstructed 'size' of this link--scale should not grow outside a range where scale*dStructureSize is too big or too small
                //C3dPoints * pStructureScaledFrom;
        public:
                CMPScale(const double dScale,const CScaleDistn scaleDistn,const int nInlierCount, const CMapPos * pRelMapPosForThisScale, const double dStructureSize = -1) : dScale(dScale), scaleDistn(scaleDistn), nInlierCount(nInlierCount), pRelMapPosForThisScale(pRelMapPosForThisScale), dStructureSize(dStructureSize) {}
                CMPScale() : dScale(1.0), scaleDistn(0), nInlierCount(0), pRelMapPosForThisScale(0) {}
                double scale() const { return dScale; }
                CScaleDistn scaleDistn() const { return scaleDistn; }
                double inlierCount() const { return nInlierCount; }
                const CMapPos * relMapPosForThisScale() const { return pRelMapPosForThisScale; }
                double structureSize() const { return dStructureSize; }
        };
private:
        class CMPScales
        {
                CMPScales(const CMPScales &) {}

                CDynArray<CMPScale> aMPScales;
                int nBestScale;
                CMPScales(const CMPScale &mpScale) : aMPScales(1), nBestScale(0) { aMPScales[0] = mpScale; }

                static double MAX_STRUCT_SIZE, MIN_STRUCT_SIZE;
                static bool CONSTRAIN_STRUCTURE_SIZE;
        public:
                static void initStructSize(const CStructureSizeParams & BOWSLAMPARAMS)
                {
                        MAX_STRUCT_SIZE = BOWSLAMPARAMS.STRUCT_SIZE_UB, MIN_STRUCT_SIZE = BOWSLAMPARAMS.STRUCT_SIZE_LB; CONSTRAIN_STRUCTURE_SIZE = BOWSLAMPARAMS.CONSTRAIN_STRUCTURE_SIZE;
                }
                CMPScales() : nBestScale(-1) { aMPScales.reserve(10 / *2+2*CBoWSLAM_Config::CLinkSelection::NUM_TO_LINK* /); }

                void addScale(const CMPScale &mpScale)
                {
                        DEBUGONLY(for(int i=0; i<(int)aMPScales.size(); i++))
                                if(IS_DEBUG) CHECK(aMPScales[i].relMapPosForThisScale() == mpScale.relMapPosForThisScale(), "Adding a duplicate map scale ?!");

                        aMPScales.push_back(mpScale);
                }

                bool updateCorrectIdx(const CMapPos * pPosRelTo)
                {
                        if(nBestScale != -1 && aMPScales[nBestScale].relMapPosForThisScale() == pPosRelTo) return true;
                        if(!pPosRelTo)
                        {
                                nBestScale = -1;
                                return false;
                        }

                        nBestScale = -1;
                        for(int i=0; i<(int)aMPScales.size(); i++)
                                if(aMPScales[i].relMapPosForThisScale() == pPosRelTo)
                                {
                                        nBestScale = i;
                                        break;
                                }

                        if(nBestScale == -1)
                        {
                                if(pPosRelTo->isOrigin())
                                {
                                        nBestScale = aMPScales.size();
                                        CMPScale originScale(1.0, 0.0, 0, pPosRelTo, -1);
                                        addScale(originScale);
                                }

                                else
                                {
                                        cout << "Could not find this map pos (" << (pPosRelTo ? pPosRelTo->timeTrigFrom() : 0) << " to "  << (pPosRelTo ? pPosRelTo->time() : 0) << ") Options: ";
                                        for(int i=0; i<(int)aMPScales.size(); i++)
                                        {
                                                const CMapPos * pOption = aMPScales[i].relMapPosForThisScale();
                                                if(pOption)
                                                        cout <<  pOption->timeTrigFrom() << " to "  << pOption->time() << ", ";
                                                else
                                                        cout << "NULL, ";
                                        }
                                        //if(IS_DEBUG) CHECK(1, "This map pos is unknown");
                                        cout << "Prob Error: This map pos is unknown--but we will just return a huge badness + never use it";
                                }
                        }
                        return (nBestScale != -1);
                }

                void printWhatImScale(const CMapPos * pPosRelTo)
                {
                        nBestScale = -1;
                        for(int i=0; i<(int)aMPScales.size(); i++)
                                if(aMPScales[i].relMapPosForThisScale() == pPosRelTo)
                                {
                                        nBestScale = i;
                                        break;
                                }
                }

                void scale(double & dPrevScale, const CMapPos * pBestPos, bool bQuiet = false) const
                {
                        const_cast<CMPScales *>(this)->updateCorrectIdx(pBestPos);
                        scale(dPrevScale, bQuiet);
                }

                double getMyRelScale() const { return nBestScale == -1 ? -1 : aMPScales[nBestScale].scale(); }

                void scale(double & dPrevScale, bool bQuiet = false) const
                {
                        //if(IS_DEBUG) CHECK(nBestScale == -1, "This map position does not know the scale of the link it is relative to");
                        if(nBestScale == -1)
                        {
                                cout << "This map position does not know the scale of the link it is relative to" << endl;
                                dPrevScale = 1;
                        }
                        else
                        {
                                double dMyScale = aMPScales[nBestScale].scale();
                                dPrevScale *= dMyScale;

                                if(CONSTRAIN_STRUCTURE_SIZE && dMyScale > 0)
                                {
                                        double dStructureSize = aMPScales[nBestScale].structureSize();

                                        if(dStructureSize>0)
                                        {
                                                if(dPrevScale * dStructureSize > MAX_STRUCT_SIZE)
                                                        dPrevScale = MAX_STRUCT_SIZE/dStructureSize;
                                                else if(dPrevScale * dStructureSize < MIN_STRUCT_SIZE)
                                                        dPrevScale = MIN_STRUCT_SIZE/dStructureSize;
                                        }
                                        else
                                        {
                                                dStructureSize = CSSObserver::typicalSS();
                                                if(dStructureSize>0)
                                                {
                                                        if(dPrevScale * dStructureSize > 2*MAX_STRUCT_SIZE)
                                                                dPrevScale = MAX_STRUCT_SIZE/dStructureSize;
                                                        else if(dPrevScale * dStructureSize < 0.5*MIN_STRUCT_SIZE)
                                                                dPrevScale = MIN_STRUCT_SIZE/dStructureSize;
                                                }
                                        }
                                }
                        }
                }

                CScaleDistn scaleDistn(const CMapPos * pBestPos) const {
                        const_cast<CMPScales *>(this)->updateCorrectIdx(pBestPos);
                        return scaleDistn();
                }
                CScaleDistn scaleDistn() const {
                        //if(IS_DEBUG) CHECK(nBestScale == -1, "This map position does not know the scale of the link it is relative to");
                        if(nBestScale == -1)
                        {
                                cout << "This map position does not know the scale of the link it is relative to" << endl;
                                return CScaleDistn::TOO_BAD_LOGBADNESS;
                        }

                        return aMPScales[nBestScale].scaleDistn();
                }

                double inlierCount() const {
                        if(IS_DEBUG) CHECK(nBestScale == -1, "This map position does not know the scale of the link it is relative to");
                        return aMPScales[nBestScale].inlierCount();
                }
        };

        C3dPoint x;
        C3dRotation Ori;
        //const TTime nTimeTrigFrom, nTimeTrigTo, nTime;
        const ePositionAlg ePosAlg;
        ePosType eIsRelative;
        const CSLAMLocation * pLocRelTo;
        CMPScales scales;
        const CSLAMLocMatch * pParent; //only needed for time() method
        bool linkedInLoop_int(int nTime) const;
public:
        static void initStructSize(const CStructureSizeParams & BOWSLAMPARAMS) { CMPScales::initStructSize(BOWSLAMPARAMS); }

        bool isOrigin() const { return eIsRelative == eOrigin; }

        CMapPos(const C3dRotation & Ori, const C3dPoint &x, const CSLAMLocMatch * pLinkTrigAccross, ePositionAlg ePosAlg, ePosType ePosCalcMethod, const CSLAMLocation * pLocRelTo = 0);

        const CMapPos * const * address() const;
        const CSLAMLocMatch * parent() const { return pParent; }

        void notifyPositionChanged();

        void addScale(const CMapPos::CMPScale & mpScale);

        void printPath() const;

        ePositionAlg getPosAlg() const
        {
                return ePosAlg;
        }

        int inlierCount() const
        {
                return scales.inlierCount();
        }
        C3dPoint pos(bool bQuiet = false) const;
        C3dRotation orientation() const;
        double scale(bool bQuiet = false) const;

        int pathLength() const;
        CScaleDistn scaleDistn(STACKCHECK(int stackDepth = 0)) const;

        TTime timeTrigFrom() const;
        TTime time() const;
        TTime relativeToTime() const;
        TTime alignedTo(int depth = 0) const;

        CvPoint cvPos(double dMapScale, C2dPoint & exactPoint) const //in x,z plane
        {
                const C3dPoint & myPos = pos(true);
                exactPoint = C2dPoint(myPos.getX(), myPos.getZ());
                CvPoint p = cvPoint(myPos.getX() * dMapScale, myPos.getZ() * dMapScale) + ORIGIN;
                return p;
        }

        //This is angle around y (vertical) axis:
        double heading() const
        {
                const C3dPoint & myOri = orientation().headingVector();
                return atan2(myOri.getX(), myOri.getZ());
        }

        //3 points giving a little arrow
        void heading(double dPxLen, CvPoint * aPoints) const
        {
                const C3dPoint & myOri = orientation().t().headingVector() * dPxLen;
                aPoints[0] = cvPoint(doubleToInt(myOri.getX()), doubleToInt(myOri.getZ()));
                C3dPoint axis(0,1,0);
                C3dRotation rLeft(axis, -0.3);
                C3dRotation rRight(axis, 0.3);
                const C3dPoint & myOriLeft = 0.8*(rLeft*myOri);
                const C3dPoint & myOriRight = 0.8*(rRight*myOri);
                aPoints[1] = cvPoint(doubleToInt(myOriLeft.getX()), doubleToInt(myOriLeft.getZ()));
                aPoints[2] = cvPoint(doubleToInt(myOriRight.getX()), doubleToInt(myOriRight.getZ()));
        }

        //3 points giving a little arrow
        void heading(double dPxLen, C2dPoint * aPoints) const
        {
                const C3dPoint & myOri = orientation().t().headingVector() * dPxLen;
                aPoints[0] = C2dPoint((myOri.getX()), (myOri.getZ()));
                C3dPoint axis(0,1,0);
                C3dRotation rLeft(axis, -0.3);
                C3dRotation rRight(axis, 0.3);
                const C3dPoint & myOriLeft = 0.8*(rLeft*myOri);
                const C3dPoint & myOriRight = 0.8*(rRight*myOri);
                aPoints[1] = C2dPoint((myOriLeft.getX()), (myOriLeft.getZ()));
                aPoints[2] = C2dPoint((myOriRight.getX()), (myOriRight.getZ()));

                if((aPoints[1] - aPoints[2]).sum_square() > sqr(2*dPxLen) || (aPoints[0] - aPoints[2]).sum_square() > sqr(2*dPxLen))
                        cout << "Error: " << aPoints[0] << aPoints[1] << aPoints[2] << rLeft << rRight << endl;
        }

        void heading(double dPxLen, double & dX, double & dY) const
        {
                const C3dPoint & myOri = orientation().t().headingVector() * dPxLen;
                dX = myOri.getX(), dY = myOri.getZ();
        }

        void pp() const
        {
                cout << "T=" << time() << ' ';
                cout << (pos());
                cout << (orientation());
        }
        bool linkedInLoop(int nTime) const;
};

double CMapPos::CMPScales::MAX_STRUCT_SIZE = HUGE; //Todo: get rid of hacks and statics
double CMapPos::CMPScales::MIN_STRUCT_SIZE = 0;
bool CMapPos::CMPScales::CONSTRAIN_STRUCTURE_SIZE = false;*/

class CScaleRange;

class C3dPoints : public CEdge {
public:
    //virtual CScale badness() const;
    virtual const CScale SLAMscale() const;
    //Return length^2, NOT including estimated scale
    virtual double measureObject(CBoWSpeedo::CBoWObjectOccurance * pObjOccurance);
    virtual void findReconstructedPairs(CBoWSpeedo::TObjectsAndLocations & localObjectFinder) const;

    NSLAMMap::CSLAMMap & slamMap;

    virtual void notifyUpdatedORScale() {
        slamMap.adjustScale(time1(), time2(), SCOREscale);
    }

    TTime time1() const;
    TTime time2() const;
private:
    const C3dPointCollection * p3dPoints;
    const CSLAMLocMatch * pParentLink;
    //const CSLAMLocation * pLocAlignedTo; //could get from parent link
    C3dRotation Rab; //Rotation from E matrix
    C3dPoint T_ba_b_dir; //Direction from E matrix
    CRelPoseSD relPoseSD; //Error in Rab and T_ba_b_dir

    void ScaleAndAlign3d(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, C3dPointMatchVector & vPointMatches, const CSLAMLocMatch * pLastLink, int nMatchingEndThis, int nMatchingEndOther, const double dConditionNum, /*CScaleRange & scales,*/ CBoWMap & map);
    void ScaleAndAlign3d_Umeyama(const C3dPointMatchVector & vPointMatches, const CSLAMLocMatch * pOtherLink);

    void get3dcorr(const C3dPoints * p3dPoint2, int nMatchingEnd1, int nMatchingEnd2, C3dPointMatchVector & vPointMatches) const;
    //bool resolveScale1d(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, const C3dPointMatch & pMatch, const C3dPoints * p_0aStruct, double &c, int nMatchingEndThis, int nMatchingEndOther, const bool bNewMethod);
    //int resolveScale1d_median(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, const C3dPointMatchVector & vPointMatches, const C3dPoints * p_0aStruct, double & dScale, bool * abInliers, int nMatchingEndThis, int nMatchingEndOther, CScaleRange & scales);
    void alignPointMatchVector(C3dPointMatchVector & vPointMatches, const C3dPoints * p_0aStruct, const int nMatchingEndThis, const int nMatchingEndOther) const;
    int resolveScale1d(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, C3dPointMatchVector & vMatches, const C3dPoints * p_0aStruct, double &c, double & dVariance, int nMatchingEndThis, int nMatchingEndOther) const;

    void addPosEstimatesToMap(const C3dPoints * pLastStruct, int nMatchingEndThis, int nMatchingEndOther, CRelScale & scaleDistn, CBoWMap & map) const;
    void getMatchingEnds(const CSLAMLocMatch * pOtherLink, int & nMatchingEndThis, int &nMatchingEndOther) const;
public:
    C3dPoints(const C3dPointCollection ** p3dPoints, const CSLAMLocMatch * pLink, const C3dRotation &Rab, const C3dPoint &T_ba_b_dir, const CRelPoseSD & relPoseSD, const double dConditionNum, CBoWMap & map, CMxLockLater & mxLockLater);
    C3dPoints(const C3dPoints * pStructToExtrapFrom, const CSLAMLocMatch * pLink, CBoWMap & map); //extrapolate

    virtual ~C3dPoints();

    double structureSize() const {
        return p3dPoints ? p3dPoints->structureSize() : -1;
    }

    /*const CMapPos * const * getAddress(const CMapPos * pMP) const
    {
            if(pPosT1 == pMP) return &pPosT1;
            else if(pPosT2 == pMP) return &pPosT2;
            THROW( "Map position not recognised")
    }
    void notifyPositionChanged()
    {
            if(pPosT1)
                    pPosT1->notifyPositionChanged();
            if(pPosT2)
                    pPosT2->notifyPositionChanged();
    }*/

    const C3dRotation & rot_ab() const {
        return Rab;
    }

    const C3dPoint & dir_ba_b() const {
        return T_ba_b_dir;
    }

    inline const C3dNormalisedPoseWithSD pose12() const {
        return C3dNormalisedPoseWithSD(Rab, -(Rab.t() * T_ba_b_dir), relPoseSD);
    }

    inline const C3dNormalisedPoseWithSD pose21() const {
        return C3dNormalisedPoseWithSD(Rab.t(), T_ba_b_dir, relPoseSD);
    }

    void ScaleAndAlign(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, const CSLAMLocMatch * pOtherLink, const double dConditionNum, CBoWMap & map);

    /*const CSLAMLocation * alignedTo() const
    {
            if(IS_DEBUG) CHECK(bAlignedAndScaled && (!pLocAlignedTo), "alignedTo: bAlignedAndScaled appears wrong");
            cout << "maybe get rid of this\n";
            return pLocAlignedTo;
    }

    bool aligned() const
    {
            if(IS_DEBUG) CHECK(bAlignedAndScaled && (!pLocAlignedTo), "aligned: pLocAlignedTo missing");
            return bAlignedAndScaled;
    }*/

    const C3dPoint * findMatch2(CLocation l) const {
        return p3dPoints->find2(l);
    }

    const C3dPoint * findMatch1(CLocation l) const {
        return p3dPoints->find1(l);
    }

    /*void AddT2Scale(const CMapPos::CMPScale & mpScale);
    void AddT1Scale(const CMapPos::CMPScale & mpScale);

    const CMapPos * positionT1() const
    {
            return pPosT1;
    }

    const CMapPos * positionT2() const
    {
            return pPosT2;
    }*/
    const CSLAMLocMatch * parentLink() const {
        return pParentLink;
    }
};

const bool bKeepAllLinks = true; //Do we keep unaligned links or not?
//A link between 2 locations

class CSLAMLocMatch {
public:

    enum eStructFoundType {
        eSamePlace, e3dStructFound, eNoMatch, ePureRotation
    };
private:
    CSLAMLocation * pLoc1, *pLoc2;
    C3dPoints * pStructureAndRelPos;
    CBoWCorrespondences const * pCorr;
    C3dPoints * getStructure(CPointVec2d & pointsCam1, CPointVec2d & pointsCam2, CInlierProbs &, CPointIdentifiers & pointIds, bool bNearby, CBoWMap & map, CRunGuiAp &gui, const int nSpareCores, CMxLockLater & mxLockLater, CSLAMLocMatch::eStructFoundType & eMotionType);
    C3dPoints * getStructure2(CPointVec2d & pointsCam1, CPointVec2d & pointsCam2, CInlierProbs &, CPointIdentifiers & pointIds, bool bNearby, CBoWMap & map, CRunGuiAp &gui, const int nSpareCores, CMxLockLater & mxLockLater, CSLAMLocMatch::eStructFoundType & eMotionType);
    C3dPoints * getStructureFromRANSACInliers(CBoWMap & map, CPointVec2d & pointsCam1, CPointVec2d & pointsCam2, C3x3MatModel & E, CMask & mask, CSLAMLocMatch::eStructFoundType & eMotionType, CMxLockLater & mxLockLater);
public:
    CSLAMLocMatch(CSLAMLocation * pLoc1, CSLAMLocation * pLoc2);

    CSLAMLocMatch(CSLAMLocation * pLoc1, CSLAMLocation * pLoc2, CBoWMap &); //Extrapolate

    void AddCorrespondences(const CBoWSLAMParams & BOWSLAMPARAMS, CBoWSpeedo & bow, CRunGuiAp &gui);
    eStructFoundType BoWCorrToStructure(CBoWMap &, CRunGuiAp & gui, const int nSpareCores, CMxLockLater & mxLockLater);

    ~CSLAMLocMatch() {
        delete pStructureAndRelPos;
        delete pCorr; //should be 0 unless ex. thrown
    }

    /*const CMapPos * const * getAddress(const CMapPos * pMP) const
    {
            return pStructureAndRelPos->getAddress(pMP);
    }*/

    CSLAMLocation * Loc1() const {
        return pLoc1;
    }

    CSLAMLocation * Loc2() const {
        return pLoc2;
    }

    /*void notifyPositionChanged()
    {
            if(pStructureAndRelPos)
                    pStructureAndRelPos->notifyPositionChanged();
    }*/

    CSLAMLocation * otherLoc(const CSLAMLocation * pOneLoc) const {
        if(IS_DEBUG) CHECK(!pOneLoc || !this || !pLoc1 || !pLoc2, "otherLoc: Uninitialised parameter/variable");
        return pLoc2 == pOneLoc ? pLoc1 : pLoc2;
    }

    C3dPoints * structureLink() const {
        return pStructureAndRelPos;
    }

    int timeDiff() const;
};

/*class CDijkstraQueue
{
public:
        typedef pair<const CSLAMLocation *, CScaleDistn > TLocBadnessPair;
private:
        typedef multimap<CScaleDistn, const CSLAMLocation *, std::greater<CScaleDistn> > TDijkstraMap;
        typedef map<const CSLAMLocation *, CScaleDistn > TDijkstraLocSet;
        TDijkstraMap mapQueue;
        TDijkstraLocSet mapLocBadnesses;
public:
        const TLocBadnessPair popTop()
        {
                if(mapQueue.size()==0) return TLocBadnessPair(0,0);
                TDijkstraMap::iterator top = mapQueue.begin();
                mapQueue.erase(top);

                TDijkstraLocSet::iterator ltop = mapLocBadnesses.find(top->second);
                mapLocBadnesses.erase(ltop);

                return *ltop;
        }

        void add(const CSLAMLocation * pLoc, CScaleDistn dBadness)
        {
                mapLocBadnesses[pLoc] = dBadness;
                mapQueue.insert(pair<CScaleDistn, const CSLAMLocation *>(dBadness, pLoc));
        }

        void update(const CSLAMLocation * pLoc, CScaleDistn dBadness)
        {
                CScaleDistn dOldBadness = mapLocBadnesses[pLoc];
                //pLoc, dBadness
                TDijkstraMap::iterator found = mapQueue.lower_bound(dOldBadness);
                while(found->second != pLoc && found->first == dOldBadness)
                        found++;

                CHECK(found == mapQueue.end() || found->first != dOldBadness, "Updating a non-existant loc");
                mapQueue.erase(found);

                add(pLoc, dBadness);
        }

        bool exists(const CSLAMLocation * pLoc) const
        {
                return mapLocBadnesses.find(pLoc) != mapLocBadnesses.end();
        }

};*/

class CSLAMLocation {
    const TTime START_FRAME;
    TTime nTime;
    typedef CDynArray<CSLAMLocMatch *> TLinkVec;
    TLinkVec vLinks; //The position of this frame is defined by the relative positions accross each link
    //typedef CDynArray<const CMapPos *> vPositions;
    /*const CMapPos * pCachedBestPos; //Cache a pointer to best position
    int nUpdateId, nUpdateIdCount;
    static int s_nUpdateId;*/
public:
    //void clearPosition() { pCachedBestPos = 0; }

    CSLAMLocation(const TTime START_FRAME, TTime nTime) : START_FRAME(START_FRAME),
    nTime(nTime)//, nUpdateId(-1), nUpdateIdCount(0)
    {
    }

    ~CSLAMLocation() {
        disconnectLinks(0);
    }

    void disconnectLinks(CBoWSpeedo::CScaleObserver * pSO) {
        //Delete all links now not referenced
        for (TLinkVec::iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++) {
            CSLAMLocMatch * pLink = *ppLink;
            if (pLink) {
                pLink->otherLoc(this)->unlink(pLink);
                if (pSO)
                    pSO->removeEdge(pLink->structureLink());
                delete pLink;
                *ppLink = 0; //should be unnecessary
            }
        }
        vLinks.clear();
    }

    void addNearby(const int nRecurseDepth, set2<TTime, std::greater<TTime> > & nearbyTimes) const {
        nearbyTimes.insert(nTime);
        if (nRecurseDepth > 0)
            for (TLinkVec::const_iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++) {
                CSLAMLocMatch * pLink = *ppLink;
                if (pLink)
                    pLink->otherLoc(this)->addNearby(nRecurseDepth - 1, nearbyTimes);
            }
    }

    void pp() const {
        /*cout << "Location at time " << nTime << ": Aligned links:\n";

        //Print my links and their quality
        for (TLinkVec::const_iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++)
        {
                CSLAMLocMatch * pLink = *ppLink;
                if (pLink && pLink->structureLink() && pLink->structureLink()->aligned())
                {
                        cout << "link to " << pLink->otherLoc(this)->time() << " aligned to " << pLink->structureLink()->alignedTo()->time() << endl;
                }
        }

        cout << "Unaligned links:\n";

        //Print my links and their quality
        for (TLinkVec::const_iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++)
        {
                CSLAMLocMatch * pLink = *ppLink;
                if (pLink && (!pLink->structureLink() || !pLink->structureLink()->aligned()))
                {
                        if (pLink->otherLoc(this)->time() < nTime) cout << "link to " << pLink->otherLoc(this)->time() << " with no aligned structure " << endl;
                }
        }*/

        //Print my positions
        /*const CMapPos * pPos = bestPosition();
        if (pPos) pPos->pp();

        cout << "All my positions:\n";
        vPositions positions;
        getAllPositions(positions, true);
        for (vPositions::const_iterator ppPos = positions.begin(); ppPos != positions.end(); ppPos++)
        {
                (*ppPos)->pp();
                //cout << "(from " << (*ppPos)->getPosAlg() << ")\n";
                cout << "(from " << (*ppPos)->relativeToTime() << ", badness = " << (*ppPos)->scaleDistn() << ")\n";
        }*/
    }

    void link(CSLAMLocMatch * pLink, const CBoWMap & map) {
        vLinks.push_back(pLink);
    }

    void unlink(CSLAMLocMatch * pLink) {
        for (TLinkVec::iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++) {
            if (pLink == *ppLink) {
                //*ppLink = 0;
                vLinks.erase(ppLink); //as may have removed one location without another
                return;
            }
        }
    }

    /*void cleanupLinks()
    {
            vLinks.cleanup();
    }*/

    const CSLAMLocMatch * getLink(TTime id) const {
        for (TLinkVec::const_iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++) {
            if ((*ppLink)->otherLoc(this)->time() == id)
                return *ppLink;
        }
        cout << "ERROR: Link doesn't exist in list. Seems to occur occasionally when extrapolating links" << endl;
        return 0;//THROW("Link doesn't exist in list");
    }

    //Return all links with positions.

    void goodLinks(/*CSLAMMap & slamMap,*/ CDynArray<CSLAMLocMatch *> & aLinks) const {
        for (TLinkVec::const_iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++) {
            CSLAMLocMatch * pLink = *ppLink;
            //cout << "found link from " << pLink->Loc1()->time() << " to " <<  pLink->Loc2()->time();

            const C3dPoints * pStrucLink = pLink->structureLink();
            if (pStrucLink) {
                /*int nLinkTime = pStrucLink->alignedTo()->time();
                if (nLinkTime < nEarliestTime)
                {
                        nEarliestTime = nLinkTime;
                        aLinks.clear();
                }
                if (nLinkTime == nEarliestTime)*/
                aLinks.push_back(pLink); //Good that we align all links--will pull in other parts of map, etc.
            }
        }
    }

public:

    TTime time() const {
        return nTime;
    }
};

class CBoundingBox {
    double dMinX, dMinY, dMaxX, dMaxY;
    double dShiftX, dShiftY, dScale;
    const int nWidth, nHeight, nMargin;
    bool bDone;
public:
    static const bool FLIP = true; //Depends on .origin.

    CBoundingBox(int nWidth, int nHeight, int nMargin) : dMinX(-1), dMinY(-1), dMaxX(1), dMaxY(1), dShiftX(0), dShiftY(0), dScale(1), nWidth(nWidth - 2 * nMargin), nHeight(nHeight - 2 * nMargin), nMargin(nMargin), bDone(false) {
    }

    void includePose(const C3dPose & pose) {
        if(IS_DEBUG) CHECK(bDone, "Can't add new pose now");

        double dXrange = dMaxX - dMinX;
        double dYrange = dMaxY - dMinY;
        double dRange = 2 * max3(dXrange, dYrange, 10.0); //stop it growing massive when there's a couple of outliers

        const C3dPoint & x = pose.t;

        if (x.getX() < dMinX) {
            if (x.getX() > dMinX - dRange)
                dMinX = x.getX();
        } else if (x.getX() > dMaxX) {
            if (x.getX() < dMaxX + dRange)
                dMaxX = x.getX();
        }

        double y = x.getZ();

        if (FLIP)
            y = -y;

        if (y < dMinY) {
            if (x.getY() > dMinY - dRange)
                dMinY = y;
        } else if (y > dMaxY) {
            if (x.getY() < dMaxY + dRange)
                dMaxY = y;
        }

        //cout << "Observed " << x << endl;
    }

    void doneObserving() {
        if(IS_DEBUG) CHECK(bDone, "Already done observing");
        bDone = true;
        double dXrange = dMaxX - dMinX;
        double dYrange = dMaxY - dMinY;
        double dScaleX = nWidth / dXrange;
        double dScaleY = nHeight / dYrange;
        //dScale = std::min<double>(dScaleX, dScaleY);
        dShiftX = dMinX; //Shift first
        dShiftY = dMinY; //Shift first
        /*cout << dMaxX << endl;
        cout << dMaxY << endl;
        cout << dScale << endl;
        cout << dShiftX << endl;
        cout << dShiftY << endl;*/
        if (dScaleX > dScaleY) {
            dScale = dScaleY;
            double dXPxSpan = dXrange * dScale;
            dShiftX = dShiftX + ((dXPxSpan - nWidth) / 2) / dScale;
        } else {
            dScale = dScaleX;
            double dYPxSpan = dYrange * dScale;
            dShiftY = dShiftY + ((dYPxSpan - nWidth) / 2) / dScale;
        }
        cout << dMaxX << endl;
        cout << dMaxY << endl;
        cout << dScale << endl;
        cout << dShiftX << endl;
        cout << dShiftY << endl;
    }

    inline CvPoint cvPos(C2dPoint & exactPoint) const //in x,z plane
    {
        return cvPos(exactPoint.getX(), exactPoint.getY());
    }

    inline CvPoint cvPos(double x, double y) const //in x,z plane
    {
        if (FLIP)
            y = -y;
        if(IS_DEBUG) CHECK(!bDone, "Not done observing");
        CvPoint p = cvPoint(nMargin + (x - dShiftX) * dScale, nMargin + (y - dShiftY) * dScale);
        return p;
    }
};

class CPoseWrapper {
public:
    const C3dPose & pose;
    const CBoundingBox & bb;

    CPoseWrapper(const C3dPose & pose, const CBoundingBox & bb) : pose(pose), bb(bb) {
    }

    CvPoint cvPos(C2dPoint & exactPoint) const //in x,z plane
    {
        const C3dPoint & myPos = pose.t;
        exactPoint = C2dPoint(myPos.getX(), myPos.getZ());

        CvPoint p = bb.cvPos(myPos.getX(), myPos.getZ());
        return p;
    }

    //This is angle around y (vertical) axis:

    double heading() const {
        const C3dPoint & myOri = pose.R.headingVector();
        return atan2(myOri.getX(), myOri.getZ());
    }

    //3 points giving a little arrow

    void heading(double dPxLen, CvPoint * aPoints) const {
        const C3dPoint & myOri = pose.R.t().headingVector() * dPxLen;
        aPoints[0] = cvPoint(doubleToInt(myOri.getX()), doubleToInt(myOri.getZ()));
        C3dPoint axis(0, 1, 0);
        C3dRotation rLeft(axis, -0.3);
        C3dRotation rRight(axis, 0.3);
        const C3dPoint & myOriLeft = 0.8 * (rLeft * myOri);
        const C3dPoint & myOriRight = 0.8 * (rRight * myOri);
        aPoints[1] = cvPoint(doubleToInt(myOriLeft.getX()), doubleToInt(myOriLeft.getZ()));
        aPoints[2] = cvPoint(doubleToInt(myOriRight.getX()), doubleToInt(myOriRight.getZ()));
        if (bb.FLIP) {
            for (int i = 0; i < 3; i++)
                aPoints[i].y = -aPoints[i].y;
        }
    }

    //3 points giving a little arrow

    void heading(double dPxLen, C2dPoint * aPoints) const {
        const C3dPoint & myOri = pose.R.t().headingVector() * dPxLen;
        aPoints[0] = C2dPoint((myOri.getX()), (myOri.getZ()));
        C3dPoint axis(0, 1, 0);
        C3dRotation rLeft(axis, -0.3);
        C3dRotation rRight(axis, 0.3);
        const C3dPoint & myOriLeft = 0.8 * (rLeft * myOri);
        const C3dPoint & myOriRight = 0.8 * (rRight * myOri);
        aPoints[1] = C2dPoint((myOriLeft.getX()), (myOriLeft.getZ()));
        aPoints[2] = C2dPoint((myOriRight.getX()), (myOriRight.getZ()));

        if ((aPoints[1] - aPoints[2]).sum_square() > sqr(2 * dPxLen) || (aPoints[0] - aPoints[2]).sum_square() > sqr(2 * dPxLen))
            cout << "Error: " << aPoints[0] << aPoints[1] << aPoints[2] << rLeft << rRight << endl;
    }

    void heading(double dPxLen, double & dX, double & dY) const {
        const C3dPoint & myOri = pose.R.t().headingVector() * dPxLen;
        dX = myOri.getX(), dY = myOri.getZ();
    }

    void pp() const {
        cout << (pose.t);
        cout << (pose.R);
    }

    /*static CvScalar getColour(CMapPos::ePositionAlg ePosAlg, bool bDim = false)
    {
            int nVal = bDim ? 180 : 255;
            switch (ePosAlg)
            {
            case CMapPos::ePosBA:
                    return CV_RGB(nVal,0,0);
            case CMapPos::ePosCameraMat:
                    return CV_RGB(0,nVal,0);
            case CMapPos::ePos3dAlignment:
                    return CV_RGB(nVal,0,nVal);
            case CMapPos::ePos1dAlignment:
                    return CV_RGB(0,0,nVal);
            case CMapPos::ePosExact:
            default:
                    return CV_RGB(0,0,0);
            }
    }*/

    static CvScalar getColour(const CScale & LB) {
        //int nVal = 255;
        int nLB = (int) 256. * LB.relBadness();

        return CV_RGB(nLB, 0, 255 - nLB);
    }

    void mark(const CScale & badness, IplImage * img, LibBoard::Board & board, bool bAlignedTo0Only, const TTime nTime, int nFrameNumbers, const C3dPose * pSourcePose) const {
        if (bAlignedTo0Only && !badness.hasScale()) return; //NB atm we only get positions aligned to 0...

        CvScalar col = getColour(badness);

        C2dPoint exactPos;
        CvPoint p = cvPos(exactPos);

        //Iterate thru+mark links:
        /* Mark all links:
         for(TLinkVec::const_iterator ppLink = vLinks.begin(); ppLink != vLinks.end(); ppLink++)
         {
         const CSLAMLocMatch * pLink = *ppLink;
         if(pLink)
         {
         const CMapPos * pOtherMapPos =  pLink->otherLoc(this)->position();
         CvPoint otherEnd = pOtherMapPos ? pOtherMapPos->cvPos(dMapScale) : cvORIGIN;
         cvLine(img, p, otherEnd, CV_RGB(0,255,0));
         }
         }*/
        //Mark 1 (earliest link) //NB we should be marking FOREWARDS if we want the links we got the structure from
        /*const CSLAMLocMatch * pLink = earliestLink();// *ppLink;
         if(pLink)
         {
         const CMapPos * pOtherMapPos =  pLink->otherLoc(this)->bestPosition();
         if(pOtherMapPos || pLink->otherLoc(this)->time() == 0) //Even so, this isn't necessarily the location we're getting a position from
         {
         CvPoint otherEnd = pOtherMapPos ? pOtherMapPos->cvPos(dMapScale) : cvORIGIN;
         cvLine(img, p, otherEnd, CV_RGB(0,255,0));
         }
         }*/
        //CMapPos::ePositionAlg ePosAlg = pPosition->getPosAlg();
        //...and put a blob+a num on the point itself
        //cvCircle(img, p, 1, col, 2);

        //And mark the heading
        const int ARROW_LEN = 10;
        CvPoint aArrowEnds[3];
        heading(ARROW_LEN, aArrowEnds);
        CvPoint cvpheading = p + aArrowEnds[0];
        CvPoint headingL = p + aArrowEnds[1];
        CvPoint headingR = p + aArrowEnds[2];
        cvLine(img, p, cvpheading, col, 2);
        cvLine(img, headingL, cvpheading, col, 2);
        cvLine(img, headingR, cvpheading, col, 2);

        board.setLineWidth(0.5 * EPS_SCALE);
        board.setPenColorRGBi(col.val[2], col.val[1], col.val[0]);
        exactPos *= EPS_SCALE;
        double x = exactPos.getX();
        double y = exactPos.getY();

        C2dPoint aArrEnds[3];
        heading(1.2 * EPS_SCALE, aArrEnds);
        /*if(((aArrEnds[0]).sum_square() > 1) || ((aArrEnds[1]).sum_square() > 1) || ((aArrEnds[2]).sum_square() > 1))
        {
                cout << "Error getting heading: " << aArrEnds[0] << aArrEnds[1] << aArrEnds[2] << exactPos << endl;
                pPosition->heading(0.75*EPS_SCALE, aArrEnds);
        }*/

        board.drawLine(x, y, x + aArrEnds[0].getX(), y + aArrEnds[0].getY(), 3);
        board.drawLine(x + aArrEnds[1].getX(), y + aArrEnds[1].getY(), x + aArrEnds[0].getX(), y + aArrEnds[0].getY(), 3);
        board.drawLine(x + aArrEnds[2].getX(), y + aArrEnds[2].getY(), x + aArrEnds[0].getX(), y + aArrEnds[0].getY(), 3);
        /*static CDynArray<LibBoard::Point> points(5);
        points[0] = LibBoard::Point(x, y);
        points[1] = LibBoard::Point(x+aArrEnds[0].getX(), y+aArrEnds[0].getY());
        points[2] = LibBoard::Point(x+aArrEnds[1].getX(), y+aArrEnds[1].getY());
        points[3] = LibBoard::Point(x+aArrEnds[0].getX(), y+aArrEnds[0].getY());
        points[4] = LibBoard::Point(x+aArrEnds[2].getX(), y+aArrEnds[2].getY());
        board.drawPolyline(points, 0);*/
        //pPosition->heading(0.75, dX, dY);
        //board.drawArrow(x, y, x+dX, y+dY, true, 0);

        if (nFrameNumbers && nTime % nFrameNumbers == 0)// && BOWSLAMPARAMS.Im.IM_SOURCE != eSyntheticImages)
        {
            double dX = 0, dY = 0;
            char text[20];
            sprintf_s(text, 20, "%d", nTime);
            CvPoint p_text = p;
            if (aArrowEnds[0].x > abs(aArrowEnds[0].y)) {
                p_text.y += 15; //cout << "Shifting time=" << nTime;
                dY = EPS_SCALE;

            } else if (aArrowEnds[0].x < -abs(aArrowEnds[0].y)) {
                p_text.y -= 15; //cout << "Shifting time=" << nTime;
                dY = -EPS_SCALE;
            } else {
                p_text.x += 10;
                dX = EPS_SCALE;
            }
            static CvPoint p_text_lastTime = cvPoint(-1000, -1000);

            if (SSD<int>(p_text, p_text_lastTime) > 400) {
                cvPutText(img, text, p_text, &font, CV_RGB(0, 150, 0));

                board.setPenColorRGBi(0, 150, 0);
                board.drawText(x + dX, y + dY, text, -1);
                p_text_lastTime = p_text;
            }
        }

        if (pSourcePose) {
            CPoseWrapper poseSource(*pSourcePose, bb);
            C2dPoint exactPosSource;
            CvPoint p_source = poseSource.cvPos(exactPosSource);
            exactPosSource *= EPS_SCALE;
            cvLine(img, p, p_source, CV_RGB(255, 255, 0));
            board.setPenColorRGBi(255, 255, 0);
            board.drawLine(x, y, exactPosSource.getX(), exactPosSource.getY(), 6);
        }
    }

    void markRobot(IplImage * pIm, LibBoard::Board & board, const TTime nTime, const int nFrameNumbers) const {
        C2dPoint exactPos;
        CvPoint p = cvPos(exactPos);

        ///////////
        const int ARROW_LEN = 10;
        CvPoint aArrowEnds[3];
        heading(ARROW_LEN, aArrowEnds);

        CvPoint cvpheading = p + aArrowEnds[0];
        CvPoint headingL = p + aArrowEnds[1];
        CvPoint headingR = p + aArrowEnds[2];
        CvScalar RED = CV_RGB(255, 0, 0);
        cvLine(pIm, p, cvpheading, RED, 4);
        cvLine(pIm, headingL, cvpheading, RED, 4);
        cvLine(pIm, headingR, cvpheading, RED, 4);
        //////////


        board.setLineWidth(0.4 * EPS_SCALE);
        board.setPenColorRGBi(255, 0, 0);
        double x = exactPos.getX() * EPS_SCALE;
        double y = exactPos.getY() * EPS_SCALE;

        double dX1 = 0, dY1 = 0;
        heading(1.75 * EPS_SCALE, dX1, dY1);
        board.drawArrow(x, y, x + dX1, y + dY1, true, 2);
        if (!nFrameNumbers) {
            double dX = 0, dY = 0;
            char text[20];
            sprintf_s(text, 20, "%d", nTime);
            if (dX1 > abs(dY1))
                dY = EPS_SCALE;
            else if (dX1 < -abs(dY1))
                dY = -EPS_SCALE;
            else
                dX = EPS_SCALE;

            board.setPenColorRGBi(0, 150, 0);
            board.drawText(x + dX, y + dY, text, 1);
        }
    }

    void markOrigin(LibBoard::Board & board) const {
        const C3dPose origin;
        CPoseWrapper pw(origin, bb);

        C2dPoint exactPos;
        pw.cvPos(exactPos);

        board.setLineWidth(0.8 * EPS_SCALE);
        board.setPenColorRGBi(255, 0, 0);
        double x = exactPos.getX() * EPS_SCALE;
        double y = exactPos.getY() * EPS_SCALE;

        board.drawLine(x - EPS_SCALE, y - EPS_SCALE, x + EPS_SCALE, y + EPS_SCALE, 2);
        board.drawLine(x - EPS_SCALE, y + EPS_SCALE, x + EPS_SCALE, y - EPS_SCALE, 2);
    }
};

void copyPoint(CvMat * cvM, int nOldCol, int nNewCol) {
    cvmSet(cvM, 0, nNewCol, cvmGet(cvM, 0, nOldCol));
    cvmSet(cvM, 1, nNewCol, cvmGet(cvM, 1, nOldCol));
}

void downweightBadInlierSet(const CMask & abInliers, CInlierProbs & adArrLikelihood) {
    int nCorrespondences = (int) abInliers.size();

    for (int nOldCorr = 0; nOldCorr < nCorrespondences; nOldCorr++) {
        if (!abInliers[nOldCorr]) {
            adArrLikelihood[nOldCorr] *= 0.5; //Todo--could just be able to do a Bayes update here? But the actual update would be tiny.
        }
    }
}

int removeBadInliers(CMask & abInliers, CPointVec2d & cvM1, CPointVec2d & cvM2, CBoWCorrespondences * pCorr, CPointIdentifiers & pointIds, CInlierProbs & adArrLikelihood) {
    int nCorrespondences = (int) cvM1.size(), nNewCorr = 0;

    for (int nOldCorr = 0; nOldCorr < nCorrespondences; nOldCorr++) {
        if (!abInliers[nOldCorr]) {
            //Was an outlier, so keep it
            cvM1[nNewCorr] = cvM1[nOldCorr];
            cvM2[nNewCorr] = cvM2[nOldCorr];
            pointIds[nNewCorr] = pointIds[nOldCorr];
            adArrLikelihood[nNewCorr] = adArrLikelihood[nOldCorr];
            nNewCorr++;
        } else
            (*pCorr)[nOldCorr].zero();
    }

    cvM1.resize(nNewCorr);
    cvM2.resize(nNewCorr);
    abInliers.resize(nNewCorr);
    pointIds.resize(nNewCorr);
    adArrLikelihood.resize(nNewCorr);

    removeZeroCorrs(pCorr);

    return nNewCorr;
}

/*bool canReconstructAccurately(const CLocation L1, const CLocation L2, const int MIN_TRACK_LEN_RECONSTRUCT)
{
        int nLengthSq = sqr(L1.dx() - L2.dx()) + sqr(L1.dy() - L2.dy());
        if (nLengthSq < sqr(MIN_TRACK_LEN_RECONSTRUCT))
        {
                //cout << "Track too short for accurate reconstruction--length^2 = " << nLengthSq << endl;
                return false;
        }
        return true;
}*/

class CBoWMap : public NSLAMMap::CSLAMMap {
public:
    const CBoWSLAMParams & BOWSLAMPARAMS;
    CImageSource * pImSource;
private:
    typedef map2<TTime, CSLAMLocation *> TMap;
    TMap locations;

    CvPtr<IplImage> img;
    LibBoard::Board board;

    CRunGuiAp &gui;

    bool locExists(TTime nLoc) {
        TMap::const_iterator ppLoc = locations.find(nLoc);
        return ppLoc != locations.end() && ppLoc->second != 0;
    }
    TTime nLastLocation, // = location robot was at at time nLastLocationTime
    nLastLocationTime, nLastTime;

    CBoWSpeedo::CScaleObserver * pSO;
public:

    CBoWMap(const CBoWSLAMParams & BOWSLAMPARAMS_in, CRunGuiAp &gui, CBoWSpeedo::CScaleObserver * pSO, CImageSource * pImSource)
    : NSLAMMap::CSLAMMap(BOWSLAMPARAMS_in.Mapping.SET_ORIGIN == CBoWSLAMParams::CMappingParams::eAtRobot, BOWSLAMPARAMS_in.Mapping.VERBOSE),
    BOWSLAMPARAMS(BOWSLAMPARAMS_in), pImSource(pImSource),
    img(cvCreateImage(cvSize(BOWSLAMPARAMS.Im.IM_HEIGHT + BOWSLAMPARAMS.Im.IM_WIDTH, BOWSLAMPARAMS.Im.IM_HEIGHT), 8, 3)), gui(gui), nLastLocation(-1), nLastLocationTime(-1), nLastTime(-1), pSO(pSO) {
        cout << "Initialising map...\n";
        img->origin = 0;
        board.setFont(LibBoard::Fonts::Helvetica, 16 * EPS_SCALE);

        //Check RANSAC thresh is sensible. Should be now its set in terms of pixel distance
        const double RANSAC_THRESH_PX = BOWSLAMPARAMS.RANSAC.E_INLIER_THRESH_PX;
        const double dThreshRelErr = RANSAC_THRESH_PX / BOWSLAMPARAMS.Corner.CORNER_LOCALISATION_SD();
        cout << "RANSAC threshhold == " << dThreshRelErr << " times av corner localisation error\n";
        CHECK(dThreshRelErr < 1, "RANSAC threshhold too small (?)");
        CHECK(dThreshRelErr > 20, "RANSAC threshhold too big (?)");
    }

    virtual ~CBoWMap() {
        cout << "Deleting map\n";
        TMap::iterator ppLocEnd = locations.end();
        for (TMap::iterator ppLoc = locations.begin(); ppLoc != ppLocEnd; ppLoc++) {
            CSLAMLocation * pos = ppLoc->second;
            delete pos;
        }
    }

    const CSLAMLocMatch * getBestParent(const CSLAMLocation * pLoc1) {
        //What link is pLoc1 positioned from?
        TTime id = positionSource(pLoc1->time());
        return pLoc1->getLink(id);
    }

    void drawMap(const IplImage * pFrame = 0) {
        if (!haveMap()) return;

        const int nComponent = currentComponent();
        aComponentData[nComponent].setDirty(); //Ensure refresh before updating scales

        if (BOWSLAMPARAMS.Output.PRINT_SPEEDS) {
            cout << "Speeds\n";
            printSpeeds(nComponent);
        }
        if (BOWSLAMPARAMS.Output.PRINT_POS_SOURCES) //Very slow
        {
            pp(nComponent);
        }

        if (BOWSLAMPARAMS.Output.OPTIMISE_SCALES) {
            //CSVDScaleOptimiser scaleOptimiser;
            CCollapsedSVDScaleOptimiser scaleOptimiser2;
            optimiseScaleAroundLoops(scaleOptimiser2, nComponent, BOWSLAMPARAMS.Optimise.SPANNER_T, BOWSLAMPARAMS.Optimise.VERBOSE);

            if (BOWSLAMPARAMS.Output.PRINT_SPEEDS) {
                cout << "Optimised speeds\n";
                printSpeeds(nComponent);
            }

            if (BOWSLAMPARAMS.TORO.SAVE_TORO_MAP && BOWSLAMPARAMS.TORO.TORO_CONNECTIVITY != BOWSLAMPARAMS.Optimise.SPANNER_T)
                cout << "Warning: BOWSLAMPARAMS.Optimise.SPANNER_T != BOWSLAMPARAMS.TORO.TORO_CONNECTIVITY; scales are being optimised around larger/smaller cycles than Toro will optimise\n" << BOWSLAMPARAMS.Optimise.SPANNER_T << BOWSLAMPARAMS.TORO.TORO_CONNECTIVITY << endl;
        }

        std::map<TTime, double> aAngleFromVertical, aRelativeDistFromGroundPlane, aDistFromGroundPlane;

        CBoundingBox bb(img->height, img->height, 20);
        TMap::iterator ppLocEnd = locations.end();
        for (TMap::iterator ppLoc = locations.begin(); ppLoc != ppLocEnd; ppLoc++) {
            const CSLAMLocation * pos = ppLoc->second;
            if(IS_DEBUG) CHECK(!pos, "All locations in map should be initialised");
            const TTime nTime = ppLoc->first;
            const C3dPose * pPose = position(nTime, nComponent);
            if (pPose)
                bb.includePose(*pPose);
        }

        CvRect boundingBox = cvRect(BOWSLAMPARAMS.Im.IM_WIDTH, 0, BOWSLAMPARAMS.Im.IM_HEIGHT, BOWSLAMPARAMS.Im.IM_HEIGHT);
        CvRect frameBoundingBox = cvRect(0, 0, BOWSLAMPARAMS.Im.IM_WIDTH, BOWSLAMPARAMS.Im.IM_HEIGHT);

        IplImage subImageMap = *img;
        CIplPx<uchar>::cropImageToRect(subImageMap, boundingBox);
        IplImage subImageFrame = *img;
        CIplPx<uchar>::cropImageToRect(subImageFrame, frameBoundingBox);

        //CPlot::setWhite(img); //won't work any more
        cvSet(&subImageMap, CV_RGB(255, 255, 255));

        int nCurrentRobotPos = currentPos();
        if (pFrame) {
            cout << subImageFrame.height << ' ' << subImageFrame.width << endl;
            cout << pFrame->height << ' ' << pFrame->width << endl;
            if (BOWSLAMPARAMS.Im.IM_CHANNELS == 1)
                cvCvtColor(pFrame, &subImageFrame, CV_GRAY2RGB);
            else
                cvCopy(pFrame, &subImageFrame);

            //For each edge pointing to current edge, mark objects
            int nPosSource = positionSource(nCurrentRobotPos);

            const CEdgeScaleObserver * pSpeedoScaleObserver = dynamic_cast<CEdgeScaleObserver *> (pSO);
            if (nPosSource >= 0 && pSpeedoScaleObserver) {
                pSpeedoScaleObserver->markObjects(nPosSource, nCurrentRobotPos, &subImageFrame);
            }

        }

        //cvCircle(&subImageMap, cvPoint(10, 10), 5, CV_RGB(0,255,255)); //Origin is top-left

        board.clear(255, 255, 255);

        bb.doneObserving();

        const bool bMarkIfAlignedTo0 = false;

        ostringstream numpy_arr_height;
        numpy_arr_height << "heights = array([";
        ostringstream numpy_arr_relheight;
        numpy_arr_relheight << "relheights = array([";
        ostringstream numpy_arr_angle;
        numpy_arr_angle << "angles = array([";

        int nFirstPos = firstPosThisComponent();
        const int nFrameNumbers = BOWSLAMPARAMS.Output.PRINT_FRAME_NUMS;

        for (TTime nTime = nFirstPos; nTime <= nCurrentRobotPos; nTime++) {
            const C3dPose * pPose = position(nTime, nComponent);

            if (pPose) {
                //Todo--ATM many have no proper map pos... pos->markAll(img, dMapScale, bMarkIfAlignedTo0, true);
                CPoseWrapper pose(*pPose, bb);
                const CScale & qual = positionQuality(nTime);
                //if(IS_DEBUG) CHECK(qual != positionQuality(nTime), "Copy CScaleDistn failing");

                const C3dPose * pSourcePose = 0;
                if (BOWSLAMPARAMS.Output.PRINT_POSITION_SOURCES) {
                    int nPosSource = positionSource(nTime);
                    if (nPosSource >= 0)
                        pSourcePose = position(nPosSource, nComponent);
                }

                pose.mark(qual, &subImageMap, board, bMarkIfAlignedTo0, nTime, nFrameNumbers, pSourcePose);

                if (nCurrentRobotPos == nTime) {
                    pose.markOrigin(board); //once...
                    pose.markRobot(&subImageMap, board, nTime, nFrameNumbers);
                }

                if (BOWSLAMPARAMS.Output.PRINT_HEURISTIC_ERRORS) {
                    double dDist = sqrt(sqr(pose.pose.t.getX()) + sqr(pose.pose.t.getZ()));
                    double dHeight = pose.pose.t.getY();
                    numpy_arr_height << "[" << nTime << "," << dHeight << "], ";
                    if (dDist > 0) {
                        numpy_arr_relheight << "[" << nTime << "," << dHeight / dDist << "], ";
                        aDistFromGroundPlane[nTime] = dHeight / dDist;
                    }
                    C3dPoint vertical(0, 1, 0);
                    vertical.rotate(pose.pose.R);
                    double dAngle = acos(vertical.getY());
                    aAngleFromVertical[nTime] = dAngle;
                    numpy_arr_angle << "[" << nTime << "," << dAngle << "], ";
                }
            }
        }

        if (BOWSLAMPARAMS.Output.PRINT_HEURISTIC_ERRORS) {
            numpy_arr_height << "]\n" << flush;
            numpy_arr_relheight << "]\n" << flush;
            numpy_arr_angle << "]\n" << flush;
            cout << numpy_arr_height.str() << numpy_arr_relheight.str() << numpy_arr_angle.str();
        }

        if (gui.showImages())
            gui.showIm(img, CRunGuiAp::eWindowMap, nCurrentRobotPos);

        if (BOWSLAMPARAMS.Output.SAVE_EPS)
            gui.saveEps(board, nCurrentRobotPos);

        if (false && aDistFromGroundPlane.size() > 0) {
            //Now plot angle errors
            const int SIZE = 400;
            const CvSize cvSIZE = cvSize(SIZE, SIZE);
            CvPtr<IplImage> angleErr(cvCreateImage(cvSIZE, IPL_DEPTH_8U, 3));
            CvPtr<IplImage> heightErr(cvCreateImage(cvSIZE, IPL_DEPTH_8U, 3));

            CPlot::plot(heightErr, aDistFromGroundPlane, true);
            CPlot::plot(angleErr, aAngleFromVertical, true);

            gui.showIm(heightErr, CRunGuiAp::eHeightPlot, nCurrentRobotPos);
            gui.showIm(angleErr, CRunGuiAp::eAnglePlot, nCurrentRobotPos);
        }

        if (BOWSLAMPARAMS.TORO.SAVE_TORO_MAP) {
            char szFilename[100];
            gui.getToroMapFilename(nCurrentRobotPos, BOWSLAMPARAMS.TORO.TORO_CONNECTIVITY, szFilename);
            saveToroMap(szFilename, BOWSLAMPARAMS.TORO.TORO_CONNECTIVITY, nComponent, BOWSLAMPARAMS.TORO.EXTRA_UNINF_EDGES, BOWSLAMPARAMS.TORO.TWO_D);

            /*gui.getToroMapFilename(nCurrentRobotPos, 0, szFilename);
            saveToroMap(szFilename, 0, nComponent, BOWSLAMPARAMS.TORO.EXTRA_UNINF_EDGES, BOWSLAMPARAMS.TORO.TWO_D); //all edges
            gui.getToroMapFilename(nCurrentRobotPos, 25, szFilename);
            saveToroMap(szFilename, 25, nComponent, BOWSLAMPARAMS.TORO.EXTRA_UNINF_EDGES, BOWSLAMPARAMS.TORO.TWO_D);*/

            //gui.getToroMapFilename2d(nCurrentRobotPos, 25, szFilename);
            //saveToroMap(szFilename, 25, nComponent, BOWSLAMPARAMS.TORO.EXTRA_UNINF_EDGES, true);
        }
    }

    double upright(const C3dRotation & R) const {
        C3dPoint vertical(0, 1, 0);
        vertical.rotate(R);
        return vertical.getY();
    }

    double forwards(const C3dRotation & R) const {
        C3dPoint forwards(0, 0, 1);
        forwards.rotate(R);
        return forwards.getZ();
    }

    double sidewards(const C3dRotation & R) const {
        C3dPoint forwards(0, 0, 1);
        forwards.rotate(R);
        return fabs(forwards.getX());
    }

    double score() {
        //Camera moves forwards, then turns 90 degrees. Maximise num of RANSAC inliers *while maintaining accuracy*
        if (!haveMap()) return 0;
        const int nComponent = currentComponent();

        double dScore = 0;
        //int nLocations = locations.size();

        TMap::const_iterator ppLocEnd = locations.end();

        int numPositions = 0;
        for (TMap::const_iterator ppLoc = locations.begin(); ppLoc != ppLocEnd; ppLoc++) {
            if (positionQuality(ppLoc->first).hasScale() && getComponentId(ppLoc->first) == nComponent)
                numPositions++;
        }
        if (numPositions < 140)
            return 0.01 * numPositions;

        int nPos1 = 1; //For NZi3
        int nPos2 = 47; //
        int nPos3 = 69; //
        int nPos4 = 126;
        /*int nPos1 = 71; //For UA
        int nPos2 = 241; //
        int nPos3 = 331; //
        int nPos4 = 496; //*/

        const C3dPose * pPose1 = 0;
        const C3dPose * pPose2 = 0;
        const C3dPose * pPose3 = 0;
        const C3dPose * pPose4 = 0;

        numPositions = 0;
        C3dPoint poseLastTime;
        int nPrevPoseTime = -1;
        CDynArray<C3dPoint> aSteps;

        for (TMap::const_iterator ppLoc = locations.begin(); ppLoc != ppLocEnd; ppLoc++) {
            int nTime = ppLoc->first;
            if (positionQuality(nTime).hasScale() && getComponentId(nTime) == nComponent) {
                const C3dPose * pPose = position(nTime, nComponent);
                if (numPositions == nPos1)
                    pPose1 = pPose;
                else if (numPositions == nPos2)
                    pPose2 = pPose;
                else if (numPositions == nPos3)
                    pPose3 = pPose;
                else if (numPositions == nPos4)
                    pPose4 = pPose;

                if (pPose) {
                    if (nPrevPoseTime >= 0) {
                        C3dPoint step = pPose->t - poseLastTime;
                        step /= (nTime - nPrevPoseTime);
                        cout << step << "=step before rot\n";
                        step.rotate(pPose->R);
                        cout << step << "=step after rot\n";
                        aSteps.push_back(step);
                    }

                    poseLastTime = pPose->t;
                    nPrevPoseTime = nTime;
                }
                numPositions++;
            }
        }

        //They're in a different component
        if (!pPose1)
            return 1;
        if (!pPose2)
            return 2;
        if (!pPose3)
            return 3;
        if (!pPose4)
            return 4;

        C3dPoint dist_1side = pPose2->t - pPose1->t;
        C3dPoint dist_2side = pPose4->t - pPose3->t;
        double dLen1 = dist_1side.length();
        double dLen2 = dist_2side.length();
        double dMinLen = min(dLen1, dLen2);
        double dMaxLen = max(dLen1, dLen2);
        if (dMinLen / dMaxLen > 0.1) {
            ofstream SDlog("ScaleDrift.tsv", ios::app);
            SDlog << dMaxLen / dMinLen << endl;
            SDlog.close();
        }

        //now look at limiting SD too :)

        dScore += upright(pPose1->R);
        dScore += upright(pPose2->R);
        dScore += upright(pPose3->R);

        dScore += forwards(pPose1->R);
        dScore += forwards(pPose2->R);
        dScore += 2 * sidewards(pPose3->R);

        double dist01 = fabs(pPose1->t.getZ());
        double dist02 = pPose2->t.getZ();
        double dist03 = pPose3->t.getZ() + fabs(pPose3->t.getX());

        double predict02 = dist01 * 2;
        double predict03 = dist01 * 4;

        if (pPose1->t.getZ() > 0) {
            dScore += 1 - (fabs(dist02 - predict02) / max<double>(predict02, dist02));
            dScore += 1 - (fabs(dist03 - predict03) / max<double>(predict03, dist03));
        }

        //Measure scale drift
        double dist34 = fabs(pPose3->t.getX() - pPose4->t.getX());
        double sd_weight = 1 - (fabs(dist01 - dist34) / max<double>(dist01, dist34));
        cout << "sd_weight: " << sd_weight << endl;
        dScore += 10 * sd_weight;
        dScore += 0.1 * sqrt(numPositions); //stop this dropping off...

        dScore += 0.01 * CLOSE_TO_CV;
        CLOSE_TO_CV = 0;

        dScore += log(max<int>(TOTAL_RANSAC_INLIERS, 1));

        TOTAL_RANSAC_INLIERS = 0;

        if (aSteps.size() < 5) {
            dScore *= 0.01;
        } else {
            C3dPoint mean, var;
            aSteps.getMeanVar(mean, var);
            var /= mean.length();
            mean.normalise();
            dScore += 10 * max<double>(mean.getZ(), 0);
            dScore += 10 / (0.2 + var.sum_square());
            cout << "Mean " << mean << endl;
            cout << "Var " << var << endl;
        }

        return dScore;
    }

    double score_old() {
        if (!haveMap()) return 0;
        const int nComponent = currentComponent();

        double dScore = 0;
        int nLocations = locations.size();

        TMap::const_iterator ppLocEnd = locations.end();
        dScore -= locations.size(); //Stop it just registering everything
        for (TMap::const_iterator ppLoc = locations.begin(); ppLoc != ppLocEnd; ppLoc++) {
            if (positionQuality(ppLoc->first).hasScale())
                dScore += 0.6;
        }

        int prevTime = -1;
        C3dPoint lastDiff; //amount it moved in 1 time step last time
        for (int i = prevTime + 1; i < nLocations; i++) {
            const CScale & LB = positionQuality(i);
            if (LB.hasScale()) {
                const C3dPose * pPose = position(i, nComponent);
                if (!pPose)
                    continue;

                const C3dPose & pose = *pPose;
                C3dRotation rot = pose.R;
                C3dPoint pos = pose.t;

                /*C3dPoint rotAxis(0,1,0);
                const double circleLength = 300.123;
                double theta = i * 2*M_PI / circleLength; //2*M_PI / 250=0.0251
                C3dRotation simRot(rotAxis, theta);

                double angleErr = abs((simRot * rot.t()).angle());
                dScore -= angleErr;

                C3dPoint pos = pose.t;

                double circlePos = (i / circleLength);
                circlePos -= floor(circlePos);
                if(circlePos < 0.8 || circlePos > 2)
                {
                        if(pos.sum_square() > 200 && pos.sum_square() < sqr(circleLength/2))
                                dScore += 3;
                }*/

                //1) Check rotation axis is near-vertical
                const bool bUpright = true;
                if (bUpright) {
                    C3dPoint vertical(0, 1, 0);
                    vertical.rotate(rot);
                    dScore += vertical.getY() > 0.6 ? 1 : 0;
                    dScore += vertical.getY() > 0.7 ? 1 : 0;
                }
                if (prevTime >= 0) {
                    //C3dRotation lastrot = position(prevTime, nComponent)->R;
                    C3dPoint lastpos = position(prevTime, nComponent)->t;

                    C3dPoint diff = pos - lastpos;
                    //if (diff.getZ() > 2 * fabs(diff.getX() + diff.getY())) dScore++;

                    if (bUpright) {
                        C3dPoint forewards(0, 0, 1);
                        forewards.rotate(rot);

                        if (diff.length() > 0) {
                            C3dPoint diffNorm = diff;
                            diffNorm /= diff.length();
                            dScore += dotproduct(diffNorm, forewards); //forwards-moving
                        }
                    }
                    C3dPoint dnorm = diff;

                    /*if(dnorm.sum_square()>0)
                    {
                            dnorm.normalise();

                            if((i>0 && i<38) || i>380)
                            {
                                    dScore += dnorm.getZ() > 0.7 ? 2 : 0;
                            }
                            if(i>60 && i<100)
                            {
                                    dScore += dnorm.getX() > 0.7 ? 2 : 0;
                            }
                            if(i>120 && i<290)
                            {
                                    dScore += dnorm.getZ() < -0.7 ? 2 : 0;
                            }
                            if(i>300 && i<360)
                            {
                                    dScore += dnorm.getX() < -0.7 ? 2 : 0;
                            }
                    }*/

                    double dTime = i * (int)BOWSLAMPARAMS.Im.VideoFile.LOAD_SUBSET / 30.0;
                    if (dTime > 0 && dTime < 11) {
                        dScore += dnorm.getZ() > 0.7 ? 2 : 0;
                    } else if (dTime > 13)
                        dScore += dnorm.getZ() < -0.4 ? 2 : 0;




                    if (prevTime > 0) //constant speed
                    {
                        double dp = dotproduct(diff, lastDiff) / lastDiff.length();
                        if (dp / (i - prevTime) > 0.7 && dp / (i - prevTime) < 1.3) dScore += 1.5;
                    }
                    lastDiff = diff * (1.0 / (i - prevTime));

                    double scale = diff.length() / (i - prevTime);
                    if (scale < 0.1 || scale > 10)
                        dScore -= 2;
                    if (scale > 0.25 && scale < 4)
                        dScore++;
                    if (scale > .75 && scale < 1.25)
                        dScore++;


                    double dBadness = LB.getG_sq();
                    /*if(dBadness > dScore*0.5)
                            dScore *= 0.5;
                    else*/
                    dScore -= dBadness;
                }

                prevTime = i;
            }
        }

        return dScore;
    }

    /*	void saveMap()
            {
                    char szFilename[200];
                    sprintf1(szFilename, 200, "maps/finalMap%s.eps", getTimeAndDate());
                    board.saveEPS(szFilename);
                    sprintf1(szFilename, 200, "maps/finalMap%s.png", getTimeAndDate());
                    cvSaveImage(szFilename, img);

    #ifdef _DEBUG
                    char szCommand[200];
                    sprintf(szCommand, "ps x | egrep 'gqX?view' > temp; if [ -s temp ]; then gqview -r \"%s\"; fi; rm temp", szFilename);
                    if(IS_DEBUG) CHECK(-1==system(szCommand), "Error displaying image");
    #endif
            }*/

    boost::mutex mxLink;

    void linkOneMT(CBoWSpeedo &bow, int nT1, int nT2, const int nSpareCores, CSLAMLocMatch::eStructFoundType & eSuccess) {
        CMxLockLater lockLater(mxLink);
        CSLAMLocMatch * pLocMatch = 0;
        eSuccess = CSLAMLocMatch::eNoMatch;
        try {
            try {
                if(IS_DEBUG) CHECK(nT1 >= nT2, "link: Can only link forward in time atm.");
                CSLAMLocation * pLoc1 = locations[nT1]; //these are reads so are ok in parallel
                CSLAMLocation * pLoc2 = locations[nT2];
                if(IS_DEBUG) CHECK(!pLoc1 || !pLoc2, "Linking non-existant location");
                pLocMatch = new CSLAMLocMatch(pLoc1, pLoc2); //link this frame to the previous with this transform
                pLocMatch->AddCorrespondences(BOWSLAMPARAMS, bow, gui);

                //This will lock the mutex, which shouldn't be unlocked until *this* fn returns
                eSuccess = pLocMatch->BoWCorrToStructure(*this, gui, nSpareCores, lockLater);
                if (eSuccess == CSLAMLocMatch::e3dStructFound || (eSuccess == CSLAMLocMatch::ePureRotation && BOWSLAMPARAMS.LinkSelection.ALLOW_ZERO_VELOCITY_LINKS)) {
                    CHECK(!lockLater.locked(), "Mutex hasn't locked");
                    //This is where we updateBestPosition:
                    pLoc1->link(pLocMatch, *this); //memory now owned by a combination of these nodes
                    pLoc2->link(pLocMatch, *this);

                    //Now add scales to OR: Need to UBP before we have an estimated scale to do OR with
                    //if(IS_DEBUG) CHECK(!pLoc1->bestPosition() || !pLoc2->bestPosition(), "Missing position for this link");
                    if (eSuccess == CSLAMLocMatch::e3dStructFound)
                        bow.addSpeedoEdge(pLocMatch->structureLink());
                } else {
                    lockLater.lock();
                    delete pLocMatch;
                    pLocMatch = 0;
                }
            } catch (CException pEx) {
                cout << "ERROR: Linking FAILED after exception caught: " << pEx.GetErrorMessage() << endl;
                throw;
            }
        } catch (...) {
            cout << "ERROR: Linking FAILED after exception caught, cleanup" << endl;
            lockLater.lock();
            delete pLocMatch;
        }
    }

    void linkExtrap(int nT1, int nT2) {
        //if(IS_DEBUG) CHECK(nT1 >= nT2, "link: Can only link forward in time atm.");
        //bool bAlignedStructure = false;
        new CSLAMLocMatch(locations[nT1], locations[nT2], *this); //link this frame to the previous with this transform
        //return bAlignedStructure;
    }

    /*template<bool bNEARBY>
    class CLinkCandidateSort;*/

    class CLinkCandidate {
        /*friend class CLinkCandidateSort<true> ;
        friend class CLinkCandidateSort<false> ;*/
        int nTemporallyClose, nMatchStrength, nGoodStructure, nId;
        CScale dPosScaleDistn;

        int temporalQualScore() const {
            if (nTemporallyClose == 0) return 0;
            return qualScore() + 2 * nTemporallyClose;
        }

        int qualScore() const {
            return 4 * nMatchStrength + 3 * nGoodStructure;
        }

    public:

        CLinkCandidate(int nTemporallyClose, int nMatchStrength, int nGoodStructure, int nId) :
        nTemporallyClose(nTemporallyClose), nMatchStrength(nMatchStrength), nGoodStructure(nGoodStructure), nId(nId) {
        }

        CLinkCandidate(int nId, int nMatchStrength) :
        nTemporallyClose(0), nMatchStrength(nMatchStrength), nGoodStructure(0), nId(nId) {
        }

        void setLB(const CScale & dLB) {
            dPosScaleDistn = dLB;
        }

        void setGoodStruct(int nGS) {
            nGoodStructure = nGS; //derived from log badness ranking
        }

        const CScale & posScaleDistn() const {
            return dPosScaleDistn;
        }

        int goodStructAndBoW() const {
            return nGoodStructure + nMatchStrength;
        }

        int id() const {
            return nId;
        }

    };

    /*template<bool bNEARBY>
    class CLinkCandidateSort: std::binary_function<const CLinkCandidate&, const CLinkCandidate&, bool>
    {
    public:
            bool operator()(const CLinkCandidate& lhs, const CLinkCandidate& rhs) const
            {
                    if (bNEARBY)
                            return lhs.temporalQualScore() > rhs.temporalQualScore();
                    else
                    //return lhs.qualScore() > rhs.qualScore();
                    {
                            //nGoodStructure is aligned-to id. Favour early ids, then similar-looking matches
                            return (lhs.nGoodStructure == BOWSLAMPARAMS.START_FRAME && lhs.nGoodStructure < rhs.nGoodStructure) || (lhs.nGoodStructure == rhs.nGoodStructure && lhs.nMatchStrength > rhs.nMatchStrength);
                    }
            }
    };*/
    class CLinkCandidateSortByBadness {
    public:

        bool operator()(const CLinkCandidate& lhs, const CLinkCandidate& rhs) const {
            return lhs.posScaleDistn().isMoreAccurateThan(rhs.posScaleDistn());
        }
    };

    class CLinkCandidateSortByBoWAndBadness {
    public:

        bool operator()(const CLinkCandidate& lhs, const CLinkCandidate& rhs) const {
            return lhs.goodStructAndBoW() > rhs.goodStructAndBoW();
        }
    };

    typedef CDynArray<CLinkCandidate> vLinkCandidates;

    class CLCManager //: private CBoWSLAM_Config::CLinkSelection
    {
        const CBoWSLAMParams::CLinkSelectionParams & LINK_PARAMS;
        CDynArray<int> anNearbyLinkCandidatesInOOP, anDistantLinkCandidatesInOOP;
        int nNumLinked, nNumLinkedLC;
        int nNumTriedLinking, nNumTriedLinkingLC;
        int nSamePlace;
        static const int DIFFERENT_PLACE = -1;
    public:

        CLCManager(const CBoWSLAMParams::CLinkSelectionParams & LINK_PARAMS) : LINK_PARAMS(LINK_PARAMS), nNumLinked(0), nNumLinkedLC(0), nNumTriedLinking(0), nNumTriedLinkingLC(0), nSamePlace(DIFFERENT_PLACE) {

        }

        int numLinked() const {
            return nNumLinked + nNumLinkedLC;
        }

        //Should this place be ignored because its in the same place as another?

        bool mergedWithExistingPlace(int & nMergePlace) const {
            nMergePlace = nSamePlace;

            if (LINK_PARAMS.KEEP_ALL_LINKS && numLinked())
                return false;
            else
                return (nSamePlace != DIFFERENT_PLACE); //Place will still be dropped if no links found
        }

        void addLCs(const vLinkCandidates & vLCs, bool bDistant = false) {
            //cout << "Added candidates:";
            for (vLinkCandidates::const_iterator pLC = vLCs.begin(); pLC < vLCs.end(); pLC++) {
                //cout << pLC->id() << ',';
                if(IS_DEBUG) CHECK(anNearbyLinkCandidatesInOOP.contains(pLC->id()), "Duplicate LC");
                if(IS_DEBUG) CHECK(anDistantLinkCandidatesInOOP.contains(pLC->id()), "Duplicate LC");

                if (bDistant)
                    anDistantLinkCandidatesInOOP.push_back(pLC->id());
                else
                    anNearbyLinkCandidatesInOOP.push_back(pLC->id());
            }
            //cout << endl;
        }

        //Returns several candidates that we can try to link concurrently
        //Each candidate is returned once, then flagged as USED
        //In TOTAL we return at most MAX_NUM_TO_TRY_LINKING + MAX_NUM_TO_TRY_LINKING_LC
        //At each time we return (NUM_TO_LINK - nNumLinked) + (NUM_TO_LINK_LC - nNumLinkedLC)

        void getLinkCandidates(CDynArray<TTime> & candidates) {
            if (nSamePlace != DIFFERENT_PLACE)
                return; //Stop linking now

            //Nearby
            {
                int nNumCandidatesNearby = min<int>(anNearbyLinkCandidatesInOOP.size() - nNumTriedLinking, LINK_PARAMS.MAX_NUM_TO_TRY_LINKING - nNumTriedLinking);
                if (nNumCandidatesNearby > LINK_PARAMS.NUM_TO_LINK - nNumLinked) nNumCandidatesNearby = LINK_PARAMS.NUM_TO_LINK - nNumLinked;

                for (int i = 0; i < nNumCandidatesNearby; i++)
                    candidates.push_back(anNearbyLinkCandidatesInOOP[nNumTriedLinking + i]);

                nNumTriedLinking += nNumCandidatesNearby;
            }

            //Distant
            {
                int nNumCandidatesDistant = min<int>(anDistantLinkCandidatesInOOP.size() - nNumTriedLinkingLC, LINK_PARAMS.MAX_NUM_TO_TRY_LINKING_LC - nNumTriedLinkingLC);
                if (nNumCandidatesDistant > LINK_PARAMS.NUM_TO_LINK_LC - nNumLinkedLC) nNumCandidatesDistant = LINK_PARAMS.NUM_TO_LINK_LC - nNumLinkedLC;

                for (int i = 0; i < nNumCandidatesDistant; i++)
                    candidates.push_back(anDistantLinkCandidatesInOOP[nNumTriedLinkingLC + i]);

                nNumTriedLinkingLC += nNumCandidatesDistant;
            }

            /*cout << "Returning candidates:";
            for(CDynArray<TTime>::iterator pn = candidates.begin(); pn != candidates.end(); pn++)
                    cout << *pn << ',';
            cout << endl;*/
        }

        void observeLinkingSuccess(TTime nSuccessfulCandidate) {
            if (anNearbyLinkCandidatesInOOP.contains(nSuccessfulCandidate))
                nNumLinked++;
            else
                nNumLinkedLC++;
        }

        void observeSamePlace(TTime nSuccessfulCandidate) {
            //observeLinkingSuccess(nSuccessfulCandidate);
            cout << nSuccessfulCandidate << " is the same place (zero dist)\n";
            if (nSamePlace != DIFFERENT_PLACE)
                cout << "TODO: Camera is in the same place as two previous locations\n";
            nSamePlace = nSuccessfulCandidate; //Todo--link to this time's parent
        }

        void pp() const {
            cout << "LCs: ";
            for (int i = 0; i < min<int>(anNearbyLinkCandidatesInOOP.size(), LINK_PARAMS.MAX_NUM_TO_TRY_LINKING); i++) {
                if (i > LINK_PARAMS.NUM_TO_LINK)
                    cout << '(';
                cout << anNearbyLinkCandidatesInOOP[i];
                if (i > LINK_PARAMS.NUM_TO_LINK)
                    cout << ')';
                cout << ',';
            }
            for (int i = 0; i < min<int>(anDistantLinkCandidatesInOOP.size(), LINK_PARAMS.MAX_NUM_TO_TRY_LINKING_LC); i++) {
                if (i > LINK_PARAMS.NUM_TO_LINK_LC)
                    cout << '(';
                cout << anDistantLinkCandidatesInOOP[i];
                if (i > LINK_PARAMS.NUM_TO_LINK_LC)
                    cout << ')';
                cout << ',';
            }
            cout << endl;
        }

    };

    void resetTime() //Tell the map the time has jumped forwards, so we don't assume we're in the same place
    {
        //nLastLocationTime = -1;
        nLastLocation = -1;
    }

    //Todo: don't add elements with only extrap. links...

    void getNearbyTimes(const int nCurrentTime, set2<TTime, std::greater<TTime> > & nearbyTimes) {
        int nRecurseDepth = BOWSLAMPARAMS.RECURSE_DEPTH;
        if (nLastLocation >= 0) {
            const CSLAMLocation * pLocLast = locations[nLastLocation];
            pLocLast->addNearby(nRecurseDepth - 1, nearbyTimes);
        }
        for (int nLastTime = nCurrentTime - 1; nLastTime >= nCurrentTime - BOWSLAMPARAMS.NEARBY_TIME; nLastTime--)
            if (nLastTime != nLastLocation && locations.exists(nLastTime)) {
                const CSLAMLocation * pLocNearby = locations[nLastTime];
                pLocNearby->addNearby(nRecurseDepth - 1, nearbyTimes);
            }

        if(IS_DEBUG) CHECK(nearbyTimes.exists(nCurrentTime), "Current time shouldn't yet be linked to nearby time");

        //cout << nCurrentTime << " nearby to: ";
        //nearbyTimes.pp();
    }

    const TBoWMatchVector * getBoWMatches(CBoWSpeedo & bow, const int nTime, const int nNumReturn) const {

        for (int breakout = 0; breakout < 20; breakout++) {
            const TBoWMatchVector * pMatches = bow.getMatches(nTime, nNumReturn); //returns ids of matching images, hopefully including last few
            if (pMatches->size())
                return pMatches;
            delete pMatches;
            pMatches = 0;
            cout << "Waiting for first round of clustering to finish...\n"; //This is happening when the WB is empty--possibly because ID lookup is wrong?
            sleep(1);
        }
        THROW("getBoWMatches: Clustering has probably crashed");
    }

    //Choose appropriate earlier nodes to link to:

    void linkToGoodMatches(CBoWSpeedo & bow, TTime nTime) {
        if (locations.size() <= 1 && !(BOWSLAMPARAMS.LinkSelection.DEBUG_DISABLE_LINKING))
            return;

        /*if(nTime % 4 == 0) //Test erasing...
        {
                locations.erase(nTime);
                bow.remove(nTime);
                nLastLocation = -1; // don't know where we're near any more...
                return;
        }*/

        vLinkCandidates /*vLCs_Nearby_Aligned0, vLCs_Nearby_Aligned0_unreliable,*/ vLCs_Nearby_goodPositions, vLCs_Nearby_noGoodPosition, vLCs_LoopClosure_goodPositions, vLCs_LoopClosure_Aligned;

        int nTryLinking = 60;
        boost::scoped_ptr<const TBoWMatchVector> pMatches(getBoWMatches(bow, nTime, nTryLinking));

        //const int nTimeNearby = FRAMERATE_MS * BOWSLAMPARAMS.NEARBY_TIME; //Last N_BOWSLAMPARAMS.NEARBY_FRAMES frames are temporally nearby
        set2<TTime, std::greater<TTime> > nearbyTimes, LCnearbyTimes;
        getNearbyTimes(nTime, nearbyTimes);

        //int nTryLinking = BOWSLAMPARAMS.LinkSelection.MAX_NUM_TO_TRY_LINKING*2; //todo--hack, debug this

        if (BOWSLAMPARAMS.LinkSelection.VERBOSE || BOWSLAMPARAMS.LinkSelection.DEBUG_DISABLE_LINKING)
            cout << pMatches->size() << " BoW matches to t=" << nTime << ": ";

        const TTime UNSET = -1;
        CLinkCandidate firstLClinkCandidate(UNSET, nTryLinking);
        bool bFoundAllNearbyMatches = false;

        int nBestLCCandidateScore = MAX_INT, nWorstNearbyToLink = MAX_INT, nBestNearbyRank = MAX_INT, nRank = 0;

        //First pass, find matches to nearby frames
        for (TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end(); pMatch++) {
            TTime lcTime = (TTime) pMatch->id();

            if (BOWSLAMPARAMS.LinkSelection.VERBOSE || BOWSLAMPARAMS.LinkSelection.DEBUG_DISABLE_LINKING)
                cout << lcTime << " (badness " << positionQuality(lcTime) << ") BoW score " << pMatch->MatchStrength() << '\n';

            if (lcTime >= nTime || !locExists(lcTime)) continue; //Do not link to 0 or future or locations that don't exist

            if(IS_DEBUG) CHECK(nTryLinking == 0, "Wrong number of matches returned")
            nTryLinking--;
            nRank++;

            if ((vLCs_Nearby_goodPositions.size() + vLCs_Nearby_noGoodPosition.size()) >= BOWSLAMPARAMS.LinkSelection.MAX_NUM_TO_TRY_LINKING)
                bFoundAllNearbyMatches = true;

            if (bFoundAllNearbyMatches && firstLClinkCandidate.id() != UNSET && !(BOWSLAMPARAMS.LinkSelection.DEBUG_DISABLE_LINKING))
                break;

            CLinkCandidate linkCandidate(lcTime, nTryLinking);

            const CScale & candidatePosBadness = positionQuality(lcTime);

            linkCandidate.setLB(candidatePosBadness);

            bool bNearby = nearbyTimes.exists(lcTime);

            if (bNearby) {
                if (nRank < nBestNearbyRank)
                    nBestNearbyRank = nRank;

                if (!bFoundAllNearbyMatches) {
                    if (candidatePosBadness.hasScale())
                        vLCs_Nearby_goodPositions.push_back(linkCandidate);
                    else
                        vLCs_Nearby_noGoodPosition.push_back(linkCandidate); //link to bad matches despite there being good matches so that these bad matches gain a position estimate

                    nWorstNearbyToLink = pMatch->MatchStrength();
                }
            } else {
                if (firstLClinkCandidate.id() == UNSET) {
                    firstLClinkCandidate = linkCandidate;
                    nBestLCCandidateScore = pMatch->MatchStrength();
                }
            }
        }

        if (BOWSLAMPARAMS.LinkSelection.VERBOSE)
            cout << endl;

        if (BOWSLAMPARAMS.LinkSelection.DEBUG_DISABLE_LINKING) {
            cout << "Return early\n";
            cout << nBestNearbyRank << " = best nearby rank\n";
            return;
        }
        //cout << "Continue\n";

        if (firstLClinkCandidate.id() != UNSET && (double) nBestLCCandidateScore > BOWSLAMPARAMS.LinkSelection.SIMILARITY_THRESH * nWorstNearbyToLink) {
            //getNearbyTimes(firstLClinkCandidate.id(), LCnearbyTimes);
            locations[firstLClinkCandidate.id()]->addNearby(BOWSLAMPARAMS.RECURSE_DEPTH, LCnearbyTimes);

            bool bGoodLCCandidates = false; //Are the matches that are not nearby all nearby each other?

            if (BOWSLAMPARAMS.LinkSelection.VERBOSE) {
                cout << "Best LC candidate is " << firstLClinkCandidate.id() << ", score = " << nBestLCCandidateScore << ", worst nearby score = " << nWorstNearbyToLink << ", nearby to ";
                LCnearbyTimes.pp();
                cout << endl;
            }

            TTime nScoreOfGoodMatches = nBestLCCandidateScore; //min<TTime>(nBestLCCandidateScore, nWorstNearbyToLink);
            TTime nNotSimilarAnymore = (int) ((double) nScoreOfGoodMatches * (double)BOWSLAMPARAMS.LinkSelection.SIMILARITY_THRESH);

            //Candidates must score at least X% of one of the nearby matches, and must be significantly more similar to the current frame than any other frames *not* in their neighbourhood

            //Second pass, can we try closing the loop?
            for (TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end(); pMatch++) {
                if (pMatch->MatchStrength() < nNotSimilarAnymore)
                    break;

                TTime lcTime = (TTime) pMatch->id();
                if (lcTime >= nTime || !locExists(lcTime) || nearbyTimes.exists(lcTime)) continue; //Do not link to 0 or future or locations that don't exist

                if (LCnearbyTimes.exists(lcTime)) //Edit here to admit more LC candidates
                {
                    //Nearby to LC candidate, so ok
                    bGoodLCCandidates = true;
                    continue;
                } else {
                    bGoodLCCandidates = false;
                    if (BOWSLAMPARAMS.LinkSelection.VERBOSE)
                        cout << "Similar but non-nearby LC candidate is " << lcTime << ", score " << pMatch->MatchStrength() << endl;
                    break;
                }
            }
            if (bGoodLCCandidates) {
                for (TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end(); pMatch++) {
                    if (pMatch->MatchStrength() < nNotSimilarAnymore)
                        break;

                    if (vLCs_LoopClosure_goodPositions.size() + vLCs_LoopClosure_Aligned.size() >= BOWSLAMPARAMS.LinkSelection.MAX_NUM_TO_TRY_LINKING_LC)
                        break;

                    TTime lcTime = (TTime) pMatch->id();
                    if (lcTime >= nTime || !locExists(lcTime) || nearbyTimes.exists(lcTime)) continue; //Do not link to 0 or future or locations that don't exist

                    CLinkCandidate linkCandidate(lcTime, nTryLinking--);
                    const CScale & candidatePosBadness = positionQuality(lcTime);
                    linkCandidate.setLB(candidatePosBadness);

                    if (BOWSLAMPARAMS.LinkSelection.VERBOSE)
                        cout << "Adding LC candidate " << lcTime << " with scale " << candidatePosBadness << endl;

                    if (candidatePosBadness.hasScale()) {
                        vLCs_LoopClosure_goodPositions.push_back(linkCandidate);
                    } else //don't bother positioning???? TODO: Need to cope with multiple maps here
                        vLCs_LoopClosure_Aligned.push_back(linkCandidate);
                }
            } else {
                if (BOWSLAMPARAMS.LinkSelection.VERBOSE)
                    cout << "No LC candidates added--not sufficiently distinctive\n";
            }
        } else {
            if (BOWSLAMPARAMS.LinkSelection.VERBOSE)
                cout << "No LC candidates added--none found with similarity greater than threshhold\n";
        }


        //Nearby hack
        /*if(vLCs_Nearby_Aligned0.size() == 0)
        {
                for(TTime lcTime = nTime-1; lcTime >= nTime-2; lcTime--)
                {
                        const CMapPos * pPos = locations[lcTime]->bestPosition();
                        if(pPos && pPos->alignedTo() == BOWSLAMPARAMS.START_FRAME)
                        {
                                vLCs_Nearby_Aligned0.push_back(linkCandidate);
                        }
                }
        }*/

        //Todo: make function:
        std::stable_sort(vLCs_Nearby_goodPositions.begin(), vLCs_Nearby_goodPositions.end(), CLinkCandidateSortByBadness());

        if (vLCs_Nearby_goodPositions.size() >= 2) {
            if(IS_DEBUG) CHECK(vLCs_Nearby_goodPositions[1].posScaleDistn().isMoreAccurateThan(vLCs_Nearby_goodPositions[0].posScaleDistn()), "Sorting by badness has failed");
        }

        int structScore = vLCs_Nearby_goodPositions.size();
        for (vLinkCandidates::iterator pLC = vLCs_Nearby_goodPositions.begin(); pLC != vLCs_Nearby_goodPositions.end(); pLC++)
            pLC->setGoodStruct(structScore--);
        std::stable_sort(vLCs_Nearby_goodPositions.begin(), vLCs_Nearby_goodPositions.end(), CLinkCandidateSortByBoWAndBadness());

        //#ifdef TUNE

        if (BOWSLAMPARAMS.TUNE_PARAMS && vLCs_Nearby_goodPositions.size() == 0 && nTime > 25)
            throw nTime;
        //#endif

#define PRINTVEC(v) cout << "Linking " << nTime << " to " # v ": " << v.size() << ':'; for(vLinkCandidates::iterator pLC = v.begin(); pLC != v.end(); pLC++) cout << pLC->id() << ','; cout << endl;

        /*PRINTVEC(vLCs_Nearby_Aligned0);*/

        CHECK(!locExists(nTime), "We are trying to link to a loc that doesn't exists");

        //First try linking to best nearby:
        CLCManager linkCandidateManager(BOWSLAMPARAMS.LinkSelection);
        linkCandidateManager.addLCs(vLCs_Nearby_goodPositions);
        linkCandidateManager.addLCs(vLCs_Nearby_noGoodPosition);

        //if(BOWSLAMPARAMS.TUNE_PARAMS){
        //vLCs_LoopClosure_goodPositions.clear();
        //vLCs_LoopClosure_Aligned.clear();}

        linkCandidateManager.addLCs(vLCs_LoopClosure_goodPositions, true);
        linkCandidateManager.addLCs(vLCs_LoopClosure_Aligned, true); //will be empty
        linkCandidateManager.pp();

        tryLink(linkCandidateManager, nTime, bow);

        nLastLocation = nLastLocationTime = nTime;
        int nExistingPlace = -1;
        if (linkCandidateManager.mergedWithExistingPlace(nExistingPlace)) {
            //We found that this location is in approx. the same place as another
            nLastLocation = nExistingPlace;
            nLastLocationTime = nTime;

            CSLAMLocation * pLoc = locations[nTime];
            if (!BOWSLAMPARAMS.LinkSelection.KEEP_ALL_LINKS)
                deleteId(nTime); //slow

            pLoc->disconnectLinks(pSO);
            delete pLoc;

            locations.erase(nTime);
            bow.remove(nTime);
        } else if (BOWSLAMPARAMS.EXTRAP && linkCandidateManager.numLinked() == 0) {
            //if(nBestNearbyRank < 5)
            {
                bool bSuccessfulExtrap = false;
                for (set2<TTime, std::greater<TTime> >::const_iterator pNearby = nearbyTimes.begin(); !bSuccessfulExtrap && pNearby != nearbyTimes.end(); pNearby++) {
                    const int nNearbyTime = *pNearby;
                    if (locations.exists(nNearbyTime)) {
                        if (positionQuality(nNearbyTime).hasScale()) {
                            cout << "Extrapolating link from time " << nTime << " to nearby time " << nNearbyTime << "...\n";
                            linkExtrap(nNearbyTime, nTime);
                            bSuccessfulExtrap = true;
                        } else
                            cout << "Nearby time has no position to extrapolate from\n";
                    }
                }
                if (!bSuccessfulExtrap) {
                    cout << "Nowhere to extrapolate to, but DO NOT ERASE, might want to start new component...\n";

                    /*locations.erase(nTime);
                    bow.remove(nTime);*/
                    nLastLocation = -1; // don't know where we're near any more...
                }
            }
            //else
            //cout << "Didn't extrapolate as best nearby rank was " << nBestNearbyRank << endl;
        }
    }

    // nNumToTryLinking is the max number of attempts at linking
    // nNumToTryLinkingLC is the max number of attempts at linking for LC

    void tryLink(CLCManager & linkCandidateManager, TTime nTime, CBoWSpeedo &bow) {
        const int nMaxToLink = BOWSLAMPARAMS.LinkSelection.MAX_NUM_TO_TRY_LINKING + BOWSLAMPARAMS.LinkSelection.MAX_NUM_TO_TRY_LINKING_LC;
        ARRAY(CSLAMLocMatch::eStructFoundType, aeSuccess, nMaxToLink);
        CDynArrayOwner<boost::thread> aThreads(nMaxToLink);

        const bool SINGLE_THREAD = true || IS_DEBUG || LESS_THREADS; //makes output a lot more readable!
        REPEAT(1, cout << "Single thread linking" << endl); //see line above

        for (;;) {
            CDynArray<TTime> anLinkCandidates;
            linkCandidateManager.getLinkCandidates(anLinkCandidates);

            const int nThreads = anLinkCandidates.size();

            if (nThreads == 0) break;

            DEBUGONLY(if (!SINGLE_THREAD && nThreads > 1) cout << "Starting " << nThreads << " linking threads\n";)

                const int nSpareCores = LESS_THREADS ? 0 : max<int>(BOWSLAMPARAMS.TOTAL_CORES - nThreads, 0);

            for (int i = 0; i < nThreads; i++) {
                aeSuccess[i] = CSLAMLocMatch::eNoMatch;

                DEBUGONLY(cout << "Linking " << anLinkCandidates[i] << " to " << nTime << endl;)

                if (i < nThreads - 1 && !SINGLE_THREAD) //debug mode--single thread
                    aThreads[i] = new boost::thread(boost::bind(&CBoWMap::linkOneMT, this, boost::ref(bow), anLinkCandidates[i], nTime, nSpareCores, boost::ref(aeSuccess[i])));
                else
                    linkOneMT(bow, anLinkCandidates[i], nTime, nSpareCores, aeSuccess[i]);
            }

            for (int i = 0; i < nThreads; i++) {
                if (i < nThreads - 1 && !SINGLE_THREAD) {
                    aThreads[i]->join();
                }

                if (aeSuccess[i] == CSLAMLocMatch::e3dStructFound || (aeSuccess[i] == CSLAMLocMatch::ePureRotation && BOWSLAMPARAMS.LinkSelection.ALLOW_ZERO_VELOCITY_LINKS))
                    linkCandidateManager.observeLinkingSuccess(anLinkCandidates[i]);
                else if (aeSuccess[i] == CSLAMLocMatch::eSamePlace || aeSuccess[i] == CSLAMLocMatch::ePureRotation) //Could use pure rotation to select link candidate that frame was registered to
                    linkCandidateManager.observeSamePlace(anLinkCandidates[i]); //For pure rotation DO NOT mark positions as the same or will be deleted
            }
        }
    }

    bool exists(TTime id) const {
        return locations.find(id) != locations.end();
    }

    CSLAMLocation * getLoc(TTime id) const {
        TMap::const_iterator pLoc = locations.find(id);
        if(IS_DEBUG) CHECK(pLoc == locations.end(), "Loc doesn't exist at this time");
        return pLoc->second;
    }

    //	void findNearbyMatches(const TBoWMatchVector * pMatches, vLinkCandidates & vLCs, TTime nTime, const eLinkMethods eLinkMeth)
    //	{
    //		//Find some nodes to link to: must be temporally nearby (not too close though), strong match, and have good structure itself
    //		int nMatches = pMatches->size();
    //		const int nTimeNearby = FRAMERATE_MS * BOWSLAMPARAMS.NEARBY_TIME; //Last N_BOWSLAMPARAMS.NEARBY_FRAMES frames are temporally nearby
    //		int nMatchStrength = nMatches;
    //
    //		for (TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end() && (int) vLCs.size() < BOWSLAMPARAMS.LinkSelection.MAX_NUM_TO_TRY_LINKING; pMatch++)
    //		{
    //			if ((TTime) pMatch->id() >= nTime || (TTime) pMatch->id() < tFirstStructureTime) continue; //Do not link to 0
    //
    //			int nStructureScore = 0;
    //			const CSLAMLocation * pLocation = locations[pMatch->id()];
    //			if (eLinkMeth != eLinkBoW)
    //			{
    //				const CMapPos * pMP = pLocation->bestPosition();
    //				if (pMP)
    //				{
    //					int nAlignedToTime = pMP->alignedTo();
    //					nStructureScore = ((nTime - nAlignedToTime) * pMatches->size()) / nTime;
    //				}
    //			}
    //
    //			int nTimeScore = 0; //higher is better, 0 is bad
    //			if (eLinkMeth == eLinkTemporally)
    //			{
    //				int nTimeBefore = nTime - pMatch->id();
    //				if (nTimeBefore < nTimeNearby)
    //				{
    //					const double dScaleTimeScore = (nMatches * 2.0) / (double) nTimeNearby;
    //					nTimeScore = (int) dScaleTimeScore * (nTimeNearby / 2.0 - abs(nTimeNearby / 2.0 - nTimeBefore));
    //				}
    //			}
    //			else if (eLinkMeth == eLinkML) //link with probabilities
    //			{
    //				THROW( "Deprecated")
    //				//nTimeScore = (int) (TIME_SCORE_WEIGHT * pLocation->prob()) + pLocation->prob() > 0.05 ? 1 : 0;
    //			}
    //			else if (eLinkMeth == eLinkBoW) //link ranked by what it is aligned to first, then BoW rank second.
    //			{
    //				//nGoodStructure is aligned-to id. Favour early ids, then similar-looking matches
    //				const CMapPos * pPos = pLocation->bestPosition();
    //				nStructureScore = nTime; //bad score
    //				if (pPos) nStructureScore = pPos->alignedTo();
    //			}
    //
    //			CLinkCandidate linkCandidate(nTimeScore, nMatchStrength--, nStructureScore, pMatch->id());
    //
    //			vLCs.push_back(linkCandidate);
    //		}
    //	}


    //Add a location, must be done in order of increasing time

    void mapLocation(CBoWSpeedo &bow, int nTime) {
        char pcOut[50];
        sprintf_s(pcOut, 50, "Mapping image, time %d...\n", nTime);
        cout << pcOut;

        if (nTime - nLastTime > BOWSLAMPARAMS.FRAMERATE_MS) {
            cout << "JUMP detected, reseting map...\n";
            resetTime();
        }

        nLastTime = nTime;

        //Add a location node
        CSLAMLocation * pNewLoc = new CSLAMLocation(BOWSLAMPARAMS.START_FRAME, nTime);
        locations.init(nTime, pNewLoc);

        //Now add edges (geometric matches)
        linkToGoodMatches(bow, nTime);

        sprintf_s(pcOut, 50, "Done mapping image, time %d\n", nTime);
        cout << pcOut;
    }

};

class CClusterDraw : public CClusterDrawer {
    const CRunGuiAp & gui;
    const int nWide;
    int nTime;
public:

    CClusterDraw(const CRunGuiAp & gui, const int nWide = 10) : gui(gui), nWide(nWide), nTime(0) {
    }

    virtual void drawCS(const CClusterSet * pCS, int nLevel) {
        if (pCS->size() == 0 || !pCS->get(0)->Centre()->drawable())
            return;

        char szOutDir[1000];
        gui.getClusterFilename(nTime++, nLevel, szOutDir);

        const int spacing = 2, bigspacing = 15;

        CvSize size = cvSize(-1, -bigspacing);
        int nDiameter = 0;
        for (CClusterSet::const_iterator ppCluster = pCS->begin(); ppCluster != pCS->end(); ppCluster++) {
            CCluster * pCluster = *ppCluster;
            nDiameter = pCluster->Centre()->diameter();
            size.height += bigspacing;
            size.height += (pCluster->Count() / nWide + ((pCluster->Count() % nWide > 0) ? 1 : 0)) * (nDiameter + spacing);
            cout << pCluster->Count() << " = cluster count\n";
        }
        size.width = nWide * (nDiameter + spacing) - spacing;

        //Now draw them:
        CvPtr<IplImage> pDescriptorImage(cvCreateImage(size, IPL_DEPTH_8U, 3));
        cvSet(pDescriptorImage, CV_RGB(255, 255, 255));

        int top = 0;
        for (CClusterSet::const_iterator ppCluster = pCS->begin(); ppCluster != pCS->end(); ppCluster++) {
            CCluster * pCluster = *ppCluster;
            int left = 0;
            const CDescriptorSet & DS = *(pCluster->Members());
            for (int i = 0; i < DS.Count(); i++) {
                if (i > 0 && i % nWide == 0) {
                    top += nDiameter + spacing;
                    left = 0;
                }
                const CDescriptor * pD = DS.get_const(i);
                for (int x = 0; x < nDiameter; x++)
                    for (int y = 0; y < nDiameter; y++) {
                        uchar r, g, b;
                        r = g = b = pD->val(x, y, 0);
                        CIplPx<uchar>::setRGB(pDescriptorImage, x + left, y + top, r, g, b);
                    }

                left += nDiameter + spacing;
            }
            top += nDiameter + spacing + bigspacing;
        }
        cvSaveImage(szOutDir, pDescriptorImage);

    }
};

class CBoWSLAM {
    static const int CAPTURE = -1;
    static const TTime FINISHED = (TTime) (-1);

public:
    const CBoWSLAMParams & BOWSLAMPARAMS; //not const cos set image size on-the-fly
private:

    int GetImageByIdx(CImageSource * pImageLoader, IplImage * pImgRGB, int imageIdx = CAPTURE) {
        int nId = -1;
        if (imageIdx == CAPTURE) {
            nId = clock();
            //Todo: from cam
            if(IS_DEBUG) CHECK(nId, "not implemented");
        } else {
            //nId = imageIdx * FRAMERATE_MS;
            if (!pImageLoader->loadNextImage(nId, pImgRGB)) //Id will be set (NB may jump)
                return FINISHED;

            if (nId == BOWSLAMPARAMS.Output.KILL_FRAME) {
                cout << "Triggering segfault to find mem use in Valgrind..." << endl;
                int * pn = 0;
                return *pn;
            }
        }
        /*if (BOWSLAMPARAMS.Im.IM_SOURCE == eSyntheticMatches)
                return nId < nSyntheticTimesteps ? nId : FINISHED;
        else*/
        return nId;

    }

    class locationIdsTS {
        boost::mutex mlocationIds;
        queue<TTime> locations;
        bool bFinished;
        boost::interprocess::interprocess_semaphore imReadySem;
        boost::interprocess::interprocess_semaphore imWaitForConsumerSem;
        const int READ_AHEAD_LIM;
        bool bWaitingForExplicitPost;
    public:

        locationIdsTS(const int READ_AHEAD_LIM_in) :
        bFinished(false), imReadySem(0), imWaitForConsumerSem(READ_AHEAD_LIM_in), READ_AHEAD_LIM(READ_AHEAD_LIM_in), bWaitingForExplicitPost(false) {
        }

        void push(TTime nLocationId) {
            if(IS_DEBUG) CHECK(bWaitingForExplicitPost, "locationIdsTS::push: bWaitingForExplicitPost indicates shouldn't be loading an image right now");
            //if(IS_DEBUG) CHECK(bFinished, "locationIdsTS::push: Finished");
            //cout << "Push " << nLocationId << endl;
            if (bFinished) return;

            {
                boost::mutex::scoped_lock scoped_lock(mlocationIds);
                locations.push(nLocationId);
                imReadySem.post();
            }
            imWaitForConsumerSem.wait();
        }

        void releaseFrame() {

            if (bWaitingForExplicitPost) {
                imWaitForConsumerSem.post();
                bWaitingForExplicitPost = false;
            }

        }

        TTime pop(bool bDelayPost) //bDelayPost stops the next image being grabbed, so we can release it after we've drawn the map
        {
            imReadySem.wait();

            if (!bDelayPost)
                imWaitForConsumerSem.post();
            else
                bWaitingForExplicitPost = true;

            if (finished()) {
                if(IS_DEBUG) CHECK((int) locations.size() > READ_AHEAD_LIM, "locationIdsTS::pop: Image reader is too far ahead");
                return FINISHED;
            }
            if(IS_DEBUG) CHECK(locations.size() == 0, "locationIdsTS::pop: Empty");
            boost::mutex::scoped_lock scoped_lock(mlocationIds);
            TTime nId = locations.front();
            locations.pop();
            //cout << "Pop " << nId << endl;
            return nId;
        }

        void doneFromConsumer() {
            boost::mutex::scoped_lock scoped_lock(mlocationIds);
            bFinished = true;
            imWaitForConsumerSem.post();
            imWaitForConsumerSem.post();
        }

        void doneFromProducer() {
            boost::mutex::scoped_lock scoped_lock(mlocationIds);
            bFinished = true;
            imReadySem.post();
            imReadySem.post();
        }

        void quit() {
            boost::mutex::scoped_lock scoped_lock(mlocationIds);
            locations = queue<TTime > ();
            bFinished = true;
            imReadySem.post();
        }

        bool finished() {
            boost::mutex::scoped_lock scoped_lock(mlocationIds);
            return locations.size() == 0 && bFinished;
        }

        void pp() const {
            cout << "Size: " << locations.size() << endl;
        }

    };
    CBoWSpeedo::CScaleObserver * pSpeedoScaleObserver;

    locationIdsTS unmappedLocations;

    CBoWSpeedo bow; //Delete first because might have clustering threads running that will want access to map for OR

    CvPtr<IplImage> pFrame;
    boost::scoped_ptr<CImageSource> pImageLoader;
    CBoWMap Map;

    boost::thread imLoadThread;
public:

    CBoWSLAM(CBoWSLAMParams & BOWSLAMPARAMS, CRunGuiAp & gui) : BOWSLAMPARAMS(BOWSLAMPARAMS),
    pSpeedoScaleObserver(new CEdgeScaleObserver(BOWSLAMPARAMS.Im)),
    //pSpeedoScaleObserver(new CBoWSpeedo::CNullScaleObserver),
    unmappedLocations(BOWSLAMPARAMS.Im.SAVE_FRAMES == CImParams::eDontSave ? BOWSLAMPARAMS.READ_AHEAD_LIM : 1000000),
    bow(BOWSLAMPARAMS.BOW, BOWSLAMPARAMS.BOWSpeedo, pSpeedoScaleObserver),
    pImageLoader(CImageSource::makeImageSource(BOWSLAMPARAMS.Im)), Map(BOWSLAMPARAMS, gui, pSpeedoScaleObserver, pImageLoader.get()),
    imLoadThread(boost::bind(&CBoWSLAM::loadImages, this, boost::ref(gui), pImageLoader.get())) {
        if (BOWSLAMPARAMS.Im.SAVE_FRAMES != CImParams::eDontSave) {
            imLoadThread.join();
            return;
        }

        if (BOWSLAMPARAMS.Im.IM_SOURCE == CImParams::eSimMap) {
            //Testing rel pose map alone
            C3dPoint straightAhead(0, 0, 1);
            C3dRotation rot0, rot45(C3dPoint(0, 1, 0), M_PI * 0.25), rot90(rot45 * rot45);
            C3dNormalisedPose pose(rot0, straightAhead);
            C3dNormalisedPose pose45(rot45, straightAhead);
            C3dNormalisedPose pose90(rot90, straightAhead);
            CRelPoseSD relPose;
            relPose.makeUninformative();
            C3dNormalisedPoseWithSD npose0(pose, relPose);
            C3dNormalisedPoseWithSD npose45(pose45, relPose);
            C3dNormalisedPoseWithSD npose90(pose90, relPose);

            CRelScale relScale123(0, 0.1), relScaleDouble(log(2), 0.1);

            Map.addBiDiScaleLink(0, 1, 2, npose0, npose45, relScale123); //1
            Map.addBiDiScaleLink(1, 2, 3, npose45, npose0, relScale123); //1
            Map.addBiDiScaleLink(2, 3, 4, npose0, npose0, relScaleDouble); //2
            Map.addBiDiScaleLink(3, 4, 5, npose0, npose0, relScaleDouble); //4
            Map.addBiDiScaleLink(4, 5, 6, npose0, npose0, relScaleDouble); //8
            Map.addBiDiScaleLink(5, 6, 7, npose0, npose0, relScaleDouble); //16
            Map.addBiDiScaleLink(6, 7, 8, npose0, npose45, relScale123);
            Map.addBiDiScaleLink(7, 8, 9, npose45, npose45, relScale123);
            Map.addBiDiScaleLink(8, 9, 10, npose45, npose45, relScale123);
            Map.addBiDiScaleLink(9, 10, 11, npose45, npose45, relScale123);
            Map.addBiDiScaleLink(10, 11, 12, npose45, npose0, relScale123);
            Map.addBiDiScaleLink(11, 12, 13, npose0, npose0, relScaleDouble.inverse()); //8
            Map.addBiDiScaleLink(12, 13, 14, npose0, npose0, relScaleDouble.inverse()); //4
            Map.addBiDiScaleLink(13, 14, 15, npose0, npose0, relScaleDouble.inverse()); //2
            Map.addBiDiScaleLink(14, 15, 16, npose0, npose0, relScaleDouble.inverse()); //1
            Map.addBiDiScaleLink(15, 16, 17, npose0, npose90, relScale123); //1
            Map.addBiDiScaleLink(16, 17, 18, npose90, npose0, relScaleDouble); //2
            Map.addBiDiScaleLink(17, 18, 19, npose0, npose0, relScaleDouble); //4
            Map.addBiDiScaleLink(18, 19, 20, npose0, npose0, relScaleDouble); //8
            Map.addBiDiScaleLink(19, 20, 21, npose0, npose0, relScale123); //8
            Map.addBiDiScaleLink(20, 21, 22, npose0, npose0, relScale123);
            Map.addBiDiScaleLink(21, 22, 23, npose0, npose0, relScale123);
            Map.addBiDiScaleLink(21, 22, 23, npose0, npose0, relScale123);

            CRelScale relScaleEighth(-log(8), 0.1), relScale16(log(16), 0.1);

            Map.addBiDiScaleLink(22, 23, 6, npose0, npose90, relScaleEighth); //1

            Map.addBiDiScaleLink(23, 6, 7, npose90, npose0, relScale16); //16

            //for(int i=0;i<4;i++)
            //Map.

            Map.drawMap();
        } else {
            CClusterDraw d(gui, 10);
            if (IS_DEBUG && BOWSLAMPARAMS.Im.IM_SOURCE != CImParams::eImageSim)
                bow.setClusterDebugger(&d);

            CStopWatch test;
            try {
                for (;;) {
                    TTime nId = getNextImId(BOWSLAMPARAMS.READ_AHEAD_LIM == 0); //if BOWSLAMPARAMS.READ_AHEAD_LIM==0 then we're making a video, so don't overwrite the frame that we'll use
                    int nTimeStep = 20;
                    if (nId % nTimeStep == 0) {
                        if (nId > 0) {
                            test.stopTimer();
                            double time = test.getElapsedTime() / nTimeStep;
                            cout << "Framerate=" << 1.0 / time << "Hz, time/frame=" << time << " id=" << nId << endl;
                            unmappedLocations.pp();
                        }

                        test.startTimer();
                    }

                    if (nId == FINISHED) return;

                    Map.mapLocation(bow, nId);

                    //if (gui.showImages())
                    if (nId > 0 && nId % BOWSLAMPARAMS.Output.PLOT_INTERVAL == 0)
                        if (!BOWSLAMPARAMS.MULTI_RUNS || nId >= 500)
                            Map.drawMap(pFrame);

                    unmappedLocations.releaseFrame(); //Can delay release of frame when using it to make video
                }
            } catch (int n) {
                unmappedLocations.doneFromConsumer();
                imLoadThread.join();
                cout << "Caught integer exception, must be tuning\n";
                throw score();
            } catch (CException pEx) {
                unmappedLocations.doneFromConsumer();
                imLoadThread.join();
                throw pEx;
            }

            sleep(3);
        }
    }

    double score() {
        return Map.score();
    }

    TTime getNextImId(bool bDelayPost) {
        return unmappedLocations.pop(bDelayPost);
    }

    ~CBoWSLAM() {
        cout << "Save all and delete CBoWSLAM...";
        saveAll();
        cout << "Finished\n";
    }

    void saveAll() {
        CStopWatch s;
        s.startTimer();

        s.stopTimer();
        double dLCTime = s.getElapsedTime();
        cout << "Full refresh took " << dLCTime << "secs\n";

        Map.drawMap(pFrame);
        //Map.saveMap();
    }

    void loadImages(CRunGuiAp &gui, CImageSource * pImageLoader) {
        if (BOWSLAMPARAMS.Im.IM_SOURCE == CImParams::eSimMap)
            return;

        CEdgeScaleObserver * pSO = dynamic_cast<CEdgeScaleObserver *> (pSpeedoScaleObserver);
        if (pSO)
            pSO->initImLoader(pImageLoader);

        boost::scoped_ptr<CFeatureExtractor> pFeatureExtractor(CFeatureExtractor::makeFeatureExtractor(BOWSLAMPARAMS.Im, BOWSLAMPARAMS.Corner, BOWSLAMPARAMS.PatchDescriptor, BOWSLAMPARAMS.DescriptorSetClustering));

        pFrame = cvCreateImage(BOWSLAMPARAMS.Im.SIZE(), IPL_DEPTH_8U, BOWSLAMPARAMS.Im.IM_CHANNELS);
        pFrame->origin = 0;

        boost::scoped_ptr<CImageSimulator> pImageSimulator;

        const CCamCalibMatrix & camCalibForSynthData = BOWSLAMPARAMS.Im.getCamCalibrationMat();
        if (BOWSLAMPARAMS.Im.IM_SOURCE == CImParams::eImageSim)
            pImageSimulator.reset(new CImageSimulator(camCalibForSynthData, BOWSLAMPARAMS.DescriptorSetClustering));

        if (/*IS_DEBUG &&*/ BOWSLAMPARAMS.Output.TEST_RD_CORRECTION && BOWSLAMPARAMS.CORRECT_RD) {
            const CCamCalibMatrix & K = BOWSLAMPARAMS.Im.getCamCalibrationMat();
            if (K.canCorrectRD() && BOWSLAMPARAMS.Im.IM_SOURCE != CImParams::eImageSim) {
                if (pImageLoader->loadImage(BOWSLAMPARAMS.Output.TEST_RD_CORRECTION, pFrame)) {
                    char szFilename[1000];
                    gui.getRDFilename(szFilename);
                    testRDCorrection(K, pFrame, szFilename);
                }
            }
        }

        try {
            try {
                DEBUGONLY(cvSetZero(pFrame));

                TTime nId = 0;
                for (int nIdx = 0; nId != FINISHED; nIdx++) {
                    if (boost::filesystem::exists("quit")) {
                        remove("quit");
                        cout << "Quit detected...";
                        unmappedLocations.quit();
                        nId = FINISHED;
                    } else {
                        CDescriptorSet * pDS = 0;
                        if (BOWSLAMPARAMS.Im.IM_SOURCE == CImParams::eImageSim) {
                            nId = nIdx;
                            pDS = pImageSimulator->getSimDescriptors(nId);
                        } else {
                            nId = GetImageByIdx(pImageLoader, pFrame, nIdx);
                            if (nId != FINISHED) {
                                try {
                                    pDS = pFeatureExtractor->getDescriptors(pFrame);
                                } catch (...) {
                                    cout << "Image loader thread caught exception, attempting to continue...\n";
                                    cvSaveImage("bad.png", pFrame);
                                    if (pDS) {
                                        delete pDS;
                                        pDS = 0;
                                    }
                                    continue;
                                }
                            }
                        }

                        if (pDS == 0) {
                            nId = FINISHED;
                            cout << "No more images\n";
                            unmappedLocations.doneFromProducer();
                        } else {
                            if (pDS->Count() < 100)
                                cout << pDS->Count() << " descriptors found\n";
                            CHECK(!pDS || pDS->Count() == 0, "Trying to add image with 0 descriptors");

                            if (!BOWSLAMPARAMS.MULTI_RUNS && gui.showImages()) {
                                if (BOWSLAMPARAMS.Im.IM_SOURCE == CImParams::eImageSim)
                                    cvSetZero(pFrame);

                                for (int nPt = 0; nPt < pDS->Count(); nPt++) {
                                    CvPoint centre = locToCvPoint((*pDS)[nPt]->location());
                                    cvCircle(pFrame, centre, 3, CV_RGB(255, 0, 0), -1);
                                    /* TODO: Add param:
                                     * CvPoint shift = cvPoint(15, 15);
                                    cvRectangle(pFrame, centre+shift, centre-shift, CV_RGB(255,0,0), 2);*/
                                }
                                //gui.showIm(pFrame, CRunGuiAp::eWindowIm, (BOWSLAMPARAMS.Im.IM_SOURCE == CImParams::eImageSim || BOWSLAMPARAMS.Output.PLOT_INTERVAL > 1) ? 0 : nIdx);
                            }
                            CStopWatch s;
                            s.startTimer();

                            bow.addImage(&pDS, nId); //now owns pDS

                            s.stopTimer();
                            cout << "Add image " << nId << " took " << s.getElapsedTime() << " seconds\n";

                            if (unmappedLocations.finished()) {
                                break;
                            }

                            if (nId >= BOWSLAMPARAMS.START_FRAME)
                                unmappedLocations.push(nId);
                        }
                    }
                }
            } catch (CException pEx) {
                cout << "Image loader thread caught: " << pEx.GetErrorMessage() << endl;
                unmappedLocations.quit();
            } catch (...) {
                cout << "Image loader thread caught unknown exception\n";
                unmappedLocations.quit();
            }
        } catch (...) {
            cout << "Error thrown from error catcher\n";
        }
    }
};

//Extrapolates a match between these locations

CSLAMLocMatch::CSLAMLocMatch(CSLAMLocation * pLoc1, CSLAMLocation * pLoc2, CBoWMap & map) :
pLoc1(pLoc1), pLoc2(pLoc2), pStructureAndRelPos(0), pCorr(0) {
    if(IS_DEBUG) CHECK(!pLoc1 || !pLoc2, "Loc 1 should have a position");

    if(IS_DEBUG) CHECK(!map.positionQuality(pLoc1->time()).hasScale(), "Loc 1 should now have a position");
    if(IS_DEBUG) CHECK(map.positionQuality(pLoc2->time()).hasScale(), "Loc 2 should not have a position");

    const CSLAMLocMatch * pLocToExtrapFrom = map.getBestParent(pLoc1);
    if(pLocToExtrapFrom)
    {
        C3dPoints * pStructToExtrapFrom = pLocToExtrapFrom->structureLink();
        pStructureAndRelPos = new C3dPoints(pStructToExtrapFrom, this, map);

        pLoc1->link(this, map);
        pLoc2->link(this, map);
    }
    else
    {
        cout << "Aborting extrapolated locMatch construction before link is created" << endl;
    }
    //TODO restore...if(IS_DEBUG) CHECK(!map.positionQuality(pLoc2->time()).hasScale(), "Loc 2 should now have a position");
    //if(IS_DEBUG) CHECK(!map.positionQuality(pLoc1->time()).hasScale(), "Loc 1 should now have a position");
}

//Sets up the match between 2 locations and the unscaled unaligned 3d points
//'new' me, DO NOT alloc on stack

CSLAMLocMatch::CSLAMLocMatch(CSLAMLocation * pLoc1, CSLAMLocation * pLoc2) :
pLoc1(pLoc1), pLoc2(pLoc2), pStructureAndRelPos(0), pCorr(0) {
    if(IS_DEBUG) CHECK(!pLoc1 || !pLoc2 /*|| pLoc1->time() >= pLoc2->time()*/, "CSLAMLocMatch::CSLAMLocMatch: Bad params/*--Loc1 must be earlier than Loc2*/");
}

void CSLAMLocMatch::AddCorrespondences(const CBoWSLAMParams & BOWSLAMPARAMS, CBoWSpeedo & bow, CRunGuiAp &gui) {
    if(IS_DEBUG) CHECK(pCorr, "Correspondences shouldn't exist yet");
    if(IS_DEBUG) CHECK(pLoc1->time() >= pLoc2->time(), "link: Can only link forward in time atm.");

    //Find correspondences, try to get relative position, link
    /*if (BOWSLAMPARAMS.Im.IM_SOURCE == eSyntheticMatches)
    {
            pCorr = getSyntheticCorr(pLoc1->time(), pLoc2->time(), gui, BOWSLAMPARAMS.Im);
    }
    else*/
    {
        //for (int N = 4; N <= 4 && (!pCorr || ((int) pCorr->size() < BOWSLAMPARAMS.TARGET_CORRESPONDENCES + RANSAC_HYP_SET_SIZE)); N++)
        {
            //if(N>1) cout << "Only found " << pCorr->size() << " correspondences with N=" << N-1 << endl;
            delete pCorr;
            pCorr = 0;
            for (;;) {
                CStopWatch s;
                s.startTimer();
                pCorr = bow.getCorrespondences(pLoc2->time(), pLoc1->time(), BOWSLAMPARAMS.BOWMatching);

                s.stopTimer();
                REPEAT(10, cout << "Correspondence extraction took " << s.getElapsedTime() << " seconds\n");

                if (pCorr) {
                    if (IS_DEBUG || pCorr->size() < 50)
                        cout << "Found " << pCorr->size() << " BF_NN=" << (int) BOWSLAMPARAMS.BOWMatching.MATCH_NN << " correspondences\n";
                    break;
                }
                cout << "Waiting for first round of clustering to finish...\n";
                sleep(1);
            }
        }
    }

}

//Should be threadsafe--getStructure locks while doing 3d align

CSLAMLocMatch::eStructFoundType CSLAMLocMatch::BoWCorrToStructure(CBoWMap & map, CRunGuiAp & gui, const int nSpareCores, CMxLockLater & mxLockLater) {
    if(IS_DEBUG) CHECK(!pCorr || pStructureAndRelPos /*|| pCorr->size() == 0*/, "CSLAMLocMatch::BoWCorrToStructure: No Correspondences found");
    if (map.BOWSLAMPARAMS.MAX_TRACK_LEN.isInit())
        THROW("Cannot limit track length here any more, could re-implement in BoW matching easily");
    //limitTrackLength(pCorr, map.BOWSLAMPARAMS);

    if (pLoc2->time() - pLoc1->time() > 50) cout << "Closing a loop..." << pLoc2->time() << ',' << pLoc1->time() << "\n";

    CPointVec2d calibratedPoints1, calibratedPoints2;
    CInlierProbs adArrLikelihood;
    CPointIdentifiers pointIds;

    pCorr->calibrate(map.BOWSLAMPARAMS.Im.getCamCalibrationMat(), calibratedPoints1, calibratedPoints2, adArrLikelihood, pointIds, (map.BOWSLAMPARAMS.Im.IM_SOURCE != CImParams::eImageSim) && map.BOWSLAMPARAMS.CORRECT_RD);

    CSLAMLocMatch::eStructFoundType eMotionType = eNoMatch;
    bool bNearby = (Loc2()->time() - Loc1()->time() < map.BOWSLAMPARAMS.NEARBY_TIME * map.BOWSLAMPARAMS.FRAMERATE_MS);
    pStructureAndRelPos = getStructure2(calibratedPoints1, calibratedPoints2, adArrLikelihood, pointIds, bNearby, map, gui, nSpareCores, mxLockLater, eMotionType); //should be scaled/aligned if poss
    delete pCorr;
    pCorr = 0;

    if(IS_DEBUG) CHECK((eMotionType == e3dStructFound || (eMotionType == ePureRotation && map.BOWSLAMPARAMS.LinkSelection.ALLOW_ZERO_VELOCITY_LINKS)) == (0 == pStructureAndRelPos), "Return not compatible with motion found");

    return eMotionType;
}

int CSLAMLocMatch::timeDiff() const {
    return abs(pLoc2->time() - pLoc1->time());
}

void C3dPoints::get3dcorr(const C3dPoints * p3dPoint2, int nMatchingEnd1, int nMatchingEnd2, C3dPointMatchVector & vPointMatches) const {
    //Iterate thru p3dpoints, looking up matches with first loc in this struct with second loc in previous structure:
    C3dPointCollection::const_iterator ppPointsEnd = nMatchingEnd1 == 1 ? p3dPoints->end1() : p3dPoints->end2();
    for (C3dPointCollection::const_iterator ppPoint = nMatchingEnd1 == 1 ? p3dPoints->begin1() : p3dPoints->begin2(); ppPoint != ppPointsEnd; ppPoint++) {
        const C3dPoint * pPointHere = ppPoint->second;
        if (pPointHere) {
            CLocation locOfPointHere = ppPoint->first;
            const C3dPoint * pPointThereAligned = 0;
            if (nMatchingEnd2 == 1)
                pPointThereAligned = p3dPoint2->findMatch1(locOfPointHere);
            else
                pPointThereAligned = p3dPoint2->findMatch2(locOfPointHere);

            if (pPointThereAligned) {
                vPointMatches.push_back(C3dPointMatch(*pPointThereAligned, *pPointHere));
            }
        }
    }
}

class CScaleRange {
    double dMinScale, dMaxScale, dBestGuessScale;
public:

    CScaleRange() : dMinScale(0), dMaxScale(HUGE), dBestGuessScale(1) {
    }

    bool inRange(double & dScale) const {
        return dScale >= dMinScale && dScale <= dMaxScale;
    }

    double bestGuessScale() const {
        return dBestGuessScale;
    }
};

void C3dPoints::getMatchingEnds(const CSLAMLocMatch * pOtherLink, int & nMatchingEndThis, int &nMatchingEndOther) const {
    //There are 4 possible tracks:
    int nEarlyTime = -1, nMidTime = -1, nLateTime = -1;
    if (pOtherLink->Loc1()->time() == pParentLink->Loc1()->time()) {
        nMatchingEndThis = 1;
        nMatchingEndOther = 1;
        nEarlyTime = pOtherLink->Loc2()->time();
        nMidTime = pOtherLink->Loc1()->time();
        nLateTime = pParentLink->Loc2()->time();
    } else if (pOtherLink->Loc1()->time() == pParentLink->Loc2()->time()) {
        nMatchingEndThis = 2;
        nMatchingEndOther = 1;
        nEarlyTime = pOtherLink->Loc2()->time();
        nMidTime = pOtherLink->Loc1()->time();
        nLateTime = pParentLink->Loc1()->time();
    } else if (pOtherLink->Loc2()->time() == pParentLink->Loc1()->time()) {
        nMatchingEndThis = 1;
        nMatchingEndOther = 2;
        nEarlyTime = pOtherLink->Loc1()->time();
        nMidTime = pOtherLink->Loc2()->time();
        nLateTime = pParentLink->Loc2()->time();
    } else if (pOtherLink->Loc2()->time() == pParentLink->Loc2()->time()) {
        nMatchingEndThis = 2;
        nMatchingEndOther = 2;
        nEarlyTime = pOtherLink->Loc1()->time();
        nMidTime = pOtherLink->Loc2()->time();
        nLateTime = pParentLink->Loc1()->time();
    } else
        THROW("These links don't join");
}

void C3dPoints::ScaleAndAlign(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, const CSLAMLocMatch * pOtherLink, const double dConditionNum, CBoWMap & map) {
    if(IS_DEBUG) CHECK(!pOtherLink, "C3dPoints::ScaleAndAlign: Bad param (no link to align to, or already aligned)");

    const C3dPoints * pLastStruct = pOtherLink->structureLink();

    if(IS_DEBUG) CHECK(!pLastStruct, "C3dPoints::ScaleAndAlign: Previous link has no structure to align this link to\n"); //if(!pLastStruct)	{		cout << "WARNING: Previous link has no structure to align this link to\n";		return;	}

    C3dPointMatchVector vPointMatches;

    int nMatchingEndThis = -1, nMatchingEndOther = -1;
    getMatchingEnds(pOtherLink, nMatchingEndThis, nMatchingEndOther);

    if (pLastStruct->p3dPoints && p3dPoints) //otherwise we are extrapolating or rotating
    {
        get3dcorr(pLastStruct, nMatchingEndThis, nMatchingEndOther, vPointMatches);
    }

    ScaleAndAlign3d(RS_PARAMS, vPointMatches, pOtherLink, nMatchingEndThis, nMatchingEndOther, dConditionNum, map);
}

void C3dPoints::alignPointMatchVector(C3dPointMatchVector & vPointMatches, const C3dPoints * p_0aStruct, const int nMatchingEndThis, const int nMatchingEndOther) const {
    if(IS_DEBUG) CHECK(!p_0aStruct, "getTandR_1d: no point");

    for (C3dPointMatchVector::iterator pMatch = vPointMatches.begin(); pMatch != vPointMatches.end(); pMatch++) {
        //This struct is from a to b. The other used to be 0a but may be bc, a0, or xb
        const C3dPoint & pOther = pMatch->p1(); // point at time a (other struct)
        const C3dPoint & pThis = pMatch->p2(); // point at time b (this struct)

        C3dPoint pOther_match, pThis_match; //todo--can simplify a bit

        if (nMatchingEndOther == 1)
            pOther_match = pOther;
        else
            pOther_match = p_0aStruct->rot_ab() * pOther + p_0aStruct->dir_ba_b();

        if (nMatchingEndThis == 1)
            pThis_match = pThis;
        else
            pThis_match = rot_ab() * pThis + dir_ba_b();

        *pMatch = C3dPointMatch(pOther_match, pThis_match);
    }
}

int C3dPoints::resolveScale1d(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, C3dPointMatchVector & vMatches, const C3dPoints * p_0aStruct, double &c, double & dVariance, int nMatchingEndThis, int nMatchingEndOther) const {
    alignPointMatchVector(vMatches, p_0aStruct, nMatchingEndThis, nMatchingEndOther);

    /*	switch(RS_PARAMS.ALIGNMENT_METHOD)
            {
            case CBoWSLAMParams::CResolveScaleParams::e1dRANSAC:
                    getTandR_1d_Ransac(vMatches, RS_PARAMS.ALIGNMENT_THETA_MAX_ERR, RS_PARAMS.RANSAC_1D_ERR, RS_PARAMS.MIN_INLIERS_1DALIGN, c, abInliers, dVariance);
                    break;
            case CBoWSLAMParams::CResolveScaleParams::e1dMedian:
                    getTandR_1d_Mean(vMatches, RS_PARAMS.ALIGNMENT_THETA_MAX_ERR, c, abInliers);
                    break;
            case CBoWSLAMParams::CResolveScaleParams::e1dMean:
                    getTandR_1d_Median(vMatches, RS_PARAMS.ALIGNMENT_THETA_MAX_ERR, c, abInliers);
                    break;

            default:
                    THROW("Unhandled 3d point align method")
            }
            int nInliers = abInliers.countInliers();
            if(dVariance <= 0)
                    dVariance = 1.0 / nInliers;*/
    //return getTandR_1d_RANSACLogs(vMatches, RS_PARAMS.ALIGNMENT_THETA_MAX_ERR, RS_PARAMS.RANSAC_1D_ERR, RS_PARAMS.MIN_INLIERS_1DALIGN, c, dVariance);
    return getDG(vMatches, RS_PARAMS.ALIGNMENT_THETA_MAX_ERR, c, dVariance, RS_PARAMS.VERBOSE);

    //return abInliers.countInliers();
}

//Returns sum of squared track lengths from motion of 0.1*scale

double getCondition(const CCamera & cam, double dScale, const C3dPointMatchVector & vPointMatches, bool * abInliers) {
    double dMyConditionMin = HUGE;
    for (int nAxis = 3/* or 1 */; nAxis <= 3; nAxis++) //We're only looking at camera optical axis.
    {
        double dMyCondition = 0;
        CCamera cam_offset = cam;
        const C3dPoint offset(0.1 * dScale * C3dPoint(nAxis == 1 ? 1 : 0, nAxis == 2 ? 1 : 0, nAxis == 3 ? 1 : 0));
        cam_offset.addOffset(offset);

        for (int i = 0; i < (int) vPointMatches.size(); i++) {
            if (abInliers[i]) {
                const C3dPoint & p2 = vPointMatches[i].p2(); //,  *p2 = vPointMatches[i]->p2();
                C2dPoint err = p2.photo(cam_offset) - p2.photo(cam);
                dMyCondition += err.sum_square();
            }
        }
        //cout << dMyCondition << ' ';
        if (dMyCondition < dMyConditionMin) dMyConditionMin = dMyCondition;
    }

    return dMyConditionMin;
}

TTime C3dPoints::time1() const {
    return pParentLink->Loc1()->time();
}

TTime C3dPoints::time2() const {
    return pParentLink->Loc2()->time();
}

void C3dPoints::ScaleAndAlign3d(const CBoWSLAMParams::CResolveScaleParams & RS_PARAMS, C3dPointMatchVector & vPointMatches, const CSLAMLocMatch * pLastLink, int nMatchingEndThis, int nMatchingEndOther, const double dConditionNum, /*CScaleRange & scales,*/ CBoWMap & map) {
    C3dPoints * pLastStruct = pLastLink->structureLink();

    double dExtrapolatedCVScale = (double) pParentLink->timeDiff() / (double) pLastLink->timeDiff(); //Constant velocity model
    const int nMaxTimeDiff = max<int>(pLastLink->timeDiff(), pParentLink->timeDiff());
    const bool NEARBY = nMaxTimeDiff < RS_PARAMS.MOTION_MODEL_MAX_TIME;
    if (!NEARBY)
        dExtrapolatedCVScale = 1; //Probably closing a loop

    const double dD_MM = log(dExtrapolatedCVScale); //0 if dExtrapolatedCVScale=1 etc.
    const double dG_sq_MM = sqr(RS_PARAMS.MOTION_MODEL_SD) * nMaxTimeDiff; //max<int>(nMaxTimeDiff, RS_PARAMS.MOTION_MODEL_MAX_TIME);

    CRelScale scaleDistn(CRelScale::eUninformativeScale, dExtrapolatedCVScale);

    int nInliers = 0;

    if ((int) vPointMatches.size() >= RS_PARAMS.MIN_INLIERS_1DALIGN) {
        double dMapPosLogScale_inv = 0, dLogScaleVariance = -1; //scales.bestGuessScale(); //Scale of 0 for rot-only links--then we reinit to 1.0. Trans is often dodgy.
        nInliers = resolveScale1d(RS_PARAMS, vPointMatches, pLastStruct, dMapPosLogScale_inv, dLogScaleVariance, nMatchingEndThis, nMatchingEndOther);
        if ((nInliers <= RS_PARAMS.MIN_INLIERS_1DALIGN)) {
            cout << "Not enough 3d points (" << nInliers << " of " << vPointMatches.size() << " remain, need " << (int) RS_PARAMS.MIN_INLIERS_1DALIGN << ")\n";
        } else {
            double dMapPosLogScale = dMapPosLogScale_inv;
            dMapPosLogScale_inv *= -1; //Swap signs

            if (dConditionNum != 1)
                cout << "Condition num " << dConditionNum << " applied...\n";

            double dLikelihood = 1;

            DEBUGONLY(cout << "Found scales " << exp(dMapPosLogScale) << ", adding to existing link " << endl;)

                    double dLogScaleVarianceMod = (dConditionNum / dLikelihood) * dLogScaleVariance;

            if (RS_PARAMS.MOTION_MODEL && NEARBY) {
                double dD_filtered = 0, dG_sq_filtered = 0;

                {
                    static int s_nCombinedScales = 0, s_nCombinedScalesBiggerMM = 0;
                    s_nCombinedScales++;
                    if (dLogScaleVarianceMod > dG_sq_MM) s_nCombinedScalesBiggerMM++;

                    PERIODIC(1000, cout << "MM more significant: " << s_nCombinedScalesBiggerMM << '/' << s_nCombinedScales << " = " << s_nCombinedScalesBiggerMM / (double) s_nCombinedScales << endl);
                }

                CVariablesFromMultiDistnsObserved::combineTwoObservations(dMapPosLogScale, dLogScaleVarianceMod, dD_MM, dG_sq_MM, dD_filtered, dG_sq_filtered);
                scaleDistn = CRelScale(dD_filtered, dG_sq_filtered);
            } else
                scaleDistn = CRelScale(dMapPosLogScale, dLogScaleVarianceMod);

            //Heuristic check that we're going at a constant speed
            /*double dCalcScale = exp(dMapPosLogScale);
            CLOSE_TO_CV += cube(1 - (fabs(dCalcScale - dExtrapolatedCVScale)/std::max<double>(dCalcScale, dExtrapolatedCVScale)));
            if(NEARBY && (true || dCalcScale/dExtrapolatedCVScale > 3 || dCalcScale/dExtrapolatedCVScale < .3333 ))
            {
                    cout << dCalcScale << "=calculated scale, " <<  1.0/dExtrapolatedCVScale << "=constant velocity scale. Heuristic scale check failed ";
                    if(dCalcScale > dExtrapolatedCVScale)
                            cout << "hi\n";
                    else
                            cout << "lo\n";
            }
            else
            {
                    scaleDistn = CRelScale(dMapPosLogScale, dLogScaleVarianceMod );
                    DEBUGONLY(cout << scaleDistn << endl);
            }*/
        }
    } else {
        if (p3dPoints) {
            cout << "Not enough 3d point matches to try to get a scale from (" << vPointMatches.size() << ", need " << (int) RS_PARAMS.MIN_INLIERS_1DALIGN << ")\n";
            cout << "Scale distn: " << scaleDistn << endl;
        } else {
            scaleDistn = CRelScale(CRelScale::eZeroScale, 0);
            cout << "Pure rotation: adding a scale 0\n";
        }
    }

    addPosEstimatesToMap(pLastStruct, nMatchingEndThis, nMatchingEndOther, scaleDistn, map);
}

//May be extrapolated or real

void C3dPoints::addPosEstimatesToMap(const C3dPoints * pLastStruct, int nMatchingEndThis, int nMatchingEndOther, CRelScale & scaleDistn, CBoWMap & map) const {
    
    //cout << "My pose12: " << pose12().SD. .getSLAMscale() << endl;
    //cout << "Other pose12: " << pLastStruct->pose12().getSLAMscale() << endl;
            
    if (nMatchingEndThis == 1 && nMatchingEndOther == 2) {
        map.addBiDiScaleLink(pLastStruct->time1(), pLastStruct->time2(), time2(), pLastStruct->pose12(), pose12(), scaleDistn);
    } else if (nMatchingEndThis == 1 && nMatchingEndOther == 1) {
        map.addBiDiScaleLink(pLastStruct->time2(), pLastStruct->time1(), time2(), pLastStruct->pose21(), pose12(), scaleDistn);
    } else if (nMatchingEndThis == 2 && nMatchingEndOther == 1) {
        map.addBiDiScaleLink(pLastStruct->time2(), pLastStruct->time1(), time1(), pLastStruct->pose21(), pose21(), scaleDistn);
    } else if (nMatchingEndThis == 2 && nMatchingEndOther == 2) {
        map.addBiDiScaleLink(pLastStruct->time1(), pLastStruct->time2(), time1(), pLastStruct->pose12(), pose21(), scaleDistn);
    }
}

void C3dPoints::ScaleAndAlign3d_Umeyama(const C3dPointMatchVector & vPointMatches, const CSLAMLocMatch * pOtherLink) {
    DEPRECATED(1);
}


//Extrapolating--need to restore elsewhere

C3dPoints::C3dPoints(const C3dPoints * pStructToExtrapFrom, const CSLAMLocMatch * pLink, CBoWMap & map) //extrapolate
: CEdge(pLink->Loc1()->time(), pLink->Loc2()->time()), slamMap(map), p3dPoints(0), pParentLink(pLink) {
    if(IS_DEBUG) CHECK(!pParentLink || !pStructToExtrapFrom, "C3dPoints::C3dPoints: bad params");
    //Rab is identity
    T_ba_b_dir = pStructToExtrapFrom->dir_ba_b(); //Maybe this should be id too--or maybe needs rotating. Test me by preventing any correspondences being returned for one frame
    relPoseSD.makeUninformative();
    double dMyTimeDiff = pParentLink->timeDiff();
    double dOldTimeDiff = abs((int) pStructToExtrapFrom->id1() - (int) pStructToExtrapFrom->id2());
    double dScale_ConstSpeed = max<double>(dMyTimeDiff, dOldTimeDiff) > 5 ? 1 : dMyTimeDiff / dOldTimeDiff;
    cout << "About to extrap rel scale with speed " << dScale_ConstSpeed << " from " << dMyTimeDiff << '/' << dOldTimeDiff << endl;
    CRelScale relScale(CRelScale::eUninformativeScale, dScale_ConstSpeed);

    int nMatchingEndThis = -1, nMatchingEndOther = -1;
    getMatchingEnds(pStructToExtrapFrom->parentLink(), nMatchingEndThis, nMatchingEndOther);

    addPosEstimatesToMap(pStructToExtrapFrom, nMatchingEndThis, nMatchingEndOther, relScale, map);
}


//Always added away from some earlier time (loc1) to the latest time.

C3dPoints::C3dPoints(const C3dPointCollection ** pp3dPoints, const CSLAMLocMatch * pLink, const C3dRotation &Rab, const C3dPoint &T_ba_b_dir, const CRelPoseSD & relPoseSD, const double dConditionNum, CBoWMap & map, CMxLockLater & mxLockLater) :
CEdge(pLink->Loc1()->time(), pLink->Loc2()->time()), slamMap(map), p3dPoints(*pp3dPoints), pParentLink(pLink), /*pLocAlignedTo(pLink->Loc2()),*/ Rab(Rab), T_ba_b_dir(T_ba_b_dir), relPoseSD(relPoseSD) {
    *pp3dPoints = 0;
    mxLockLater.lock();

    cout << "Started linking " << pParentLink->Loc1()->time() << " to " << pParentLink->Loc2()->time() << endl;
    /*cout << Rab << "=my rotation ab\n";
    cout << T_ba_b_dir << "=T_ba_b\n";
    cout << Rab.t()*T_ba_b_dir << "=T_ba_a\n";
    cout << (Rab.t()*Rab.t())*T_ba_b_dir << "=T_ba_a_rotated again\n";*/

    if(IS_DEBUG) CHECK(!pParentLink /*|| !zero(1-T_ba_b_dir.sum_square())*/, "C3dPoints::C3dPoints: T_ba_b_dir.sum_square() is not 1");

    CDynArray<CSLAMLocMatch *> aLinks;
    pParentLink->Loc1()->goodLinks(aLinks);
    pParentLink->Loc2()->goodLinks(aLinks); //catches links we have just added to this loc2

    if (aLinks.size() == 0)// && /*!aLinks[0]->Loc2()->bestPosition()*/ !pLink->Loc1()->bestPosition())
    {
        if (pLink->Loc1()->time() > 0)
            cout << "WARNING: No links to position at time " << pLink->Loc1()->time() << endl;
    } else {
        //Forward position and reverse position should never be w.r.t. each other, but may check this then doing updateBestPosition...
        for (CDynArray<CSLAMLocMatch *>::iterator ppLastLink = aLinks.begin(); ppLastLink < aLinks.end(); ppLastLink++) {
            DEBUGONLY(cout << "Linking (" << pParentLink->Loc1()->time() << " to " << pParentLink->Loc2()->time() << ") to " << (*ppLastLink)->Loc1()->time() << ", " << (*ppLastLink)->Loc2()->time() << endl;)
            ScaleAndAlign(map.BOWSLAMPARAMS.ResolveScale, *ppLastLink, dConditionNum, map); //won't select this link as hasn't been added yet
        }
    }
}

C3dPoints::~C3dPoints() {
    delete p3dPoints;
}

const CScale C3dPoints::SLAMscale() const {
    return slamMap.getSLAMscale(time1(), time2());
}

void C3dPoints::findReconstructedPairs(CBoWSpeedo::TObjectsAndLocations & localObjectFinder) const {
    //iterate over detected objects and mark num that were successfully reconstructed
    cout << "Edge " << nId1 << ',' << nId2 << ": localObjectFinder: Search for " << localObjectFinder.size() << " objects...\n";
    T3dLocations locationsOfFeature;

    for (CBoWSpeedo::TObjectsAndLocations::iterator ppOL = localObjectFinder.begin(); ppOL != localObjectFinder.end(); ppOL++) {
        CBoWSpeedo::CObjAndLocations * pOL = *ppOL;

        locationsOfFeature.clear();
        p3dPoints->get3dLocations(pOL->locations(), locationsOfFeature);
        int nWordCount = locationsOfFeature.size();
        //cout << "Edge " << nId1 << ',' << nId2 << ", " << nWordCount << " of this object found\n";

        if (nWordCount > 0) {
            pOL->setReconCounts(nWordCount);
        }
    }
    //cout << "Edge " << nId1 << ',' << nId2 << ": localObjectFinder: Found " << localObjectFinder.size() << " objects...\n";
}

double C3dPoints::measureObject(CBoWSpeedo::CBoWObjectOccurance * pObjOccurance) {
    T3dLocations locationsOfFeature1, locationsOfFeature2;

    //cout << pObjOccurance->getLocOfFeature1InIm1().size() << " occurances of feature 1" << "\n";
    //if(pObjOccurance->getLocOfFeature1InIm1().size()) cout << (pObjOccurance->getLocOfFeature1InIm1())[0] << "\n";
    //cout << pObjOccurance->getLocOfFeature2InIm1().size() << " occurances of feature 2\n";
    p3dPoints->get3dLocations(pObjOccurance->getLocOfFeature1InIm1(), locationsOfFeature1);
    p3dPoints->get3dLocations(pObjOccurance->getLocOfFeature2InIm1(), locationsOfFeature2);
    //cout << locationsOfFeature1.size() << " occurances of feature 1 reconstructed\n";
    //cout << locationsOfFeature2.size() << " occurances of feature 2 reconstructed\n";

    double dMinLengthSq = HUGE;
    for (T3dLocations::const_iterator p3dPoint1 = locationsOfFeature1.begin(); p3dPoint1 != locationsOfFeature1.end(); p3dPoint1++) {
        for (T3dLocations::const_iterator p3dPoint2 = locationsOfFeature2.begin(); p3dPoint2 != locationsOfFeature2.end(); p3dPoint2++) {
            C3dPoint diff = p3dPoint1->point;
            diff -= p3dPoint2->point;
            double dLengthSq = diff.sum_square();
            if (dLengthSq < dMinLengthSq) {
                dMinLengthSq = dLengthSq;
                pObjOccurance->setLocs(p3dPoint1->loc1, p3dPoint1->loc2, p3dPoint2->loc1, p3dPoint2->loc2);
            }
        }
    }
    if (dMinLengthSq > sqr(10)) return 0; //Much longer than baseline so not accurate

    return (dMinLengthSq == HUGE) ? 0 : dMinLengthSq;
}

class CRecon3dPoint {
    C3dPoint recon;
public:

    enum eReconStatus {
        eUninit, eInFront, eBehind, eAtInfinity
    };
private:
    eReconStatus reconStatus;
public:

    CRecon3dPoint() : reconStatus(eUninit) {

    }

    eReconStatus reconstructOne(const CCamera & Pp, const C2dPoint & p1, const C2dPoint & p2, const double DEPTH_THRESH) {
        const CCamera P;
        recon = reconstruct(P, Pp, p1, p2);

        double depth1 = recon.depth(P);
        if (depth1 < 0 && depth1 > -DEPTH_THRESH)
            reconStatus = eBehind;
        else {
            double depth2 = recon.depth(Pp);
            if (depth2 < 0 && depth2 > -DEPTH_THRESH)
                reconStatus = eBehind;
            else {
                if (depth1 < 0 && depth2 < 0)
                    reconStatus = eAtInfinity;
                else {
                    if (depth1 > DEPTH_THRESH || depth2 > DEPTH_THRESH)
                        reconStatus = eAtInfinity;
                    else
                        reconStatus = eInFront;
                }
            }

        }
        //		reconStatus = recon.testInFront(P, Pp) ? eInFront : eBehind;

        return reconStatus;
    }

    double depth() const {
        const CCamera P;
        return recon.depth(P);
    }

    eReconStatus status() const {
        return reconStatus;
    }

    const C3dPoint & point3d() const {
        return recon;
    }
};

static const int NUM_CAM_MATS = 4;

class C3dReconData //all reconstructions of a point
{
    CRecon3dPoint aRecon3dPoints[NUM_CAM_MATS];
    C2dPoint p1, p2;
    CLocation loc1, loc2;
public:

    C3dReconData(const C2dPoint & p1, const C2dPoint & p2, const CLocation loc1, const CLocation loc2) : p1(p1), p2(p2), loc1(loc1), loc2(loc2) {

    }

    C3dReconData() {

    }

    const C2dPoint & point2d1() const {
        return p1;
    }

    const C2dPoint & point2d2() const {
        return p2;
    }

    const CLocation l1() const {
        return loc1;
    }

    const CLocation l2() const {
        return loc2;
    }

    int tryReconstruct(const CCamera * aPp, int nStart, const double DEPTH_THRESH) {
        int nBestCamMat = -1;

        int aanOrder[4][4] = {
            {0, 1, 2, 3},
            {1, 0, 2, 3},
            {2, 3, 0, 1},
            {3, 2, 0, 1}}; //If the first fails try the other one in the pair next
        bool bAllAtInf = true;

        for (int i = 0; i < NUM_CAM_MATS; i++) {
            int nCamToRecon = aanOrder[nStart][i];
            CRecon3dPoint::eReconStatus reconStatus = reconstructOne(nCamToRecon, aPp[nCamToRecon], DEPTH_THRESH);

            if (reconStatus != CRecon3dPoint::eAtInfinity)
                bAllAtInf = false;
            else
                if (i == 1 && bAllAtInf) {
                if (nCamToRecon < 2)
                    return -1; //both at infinity
                else
                    return -2; //both at inf, trying second 2
            }


            if (reconStatus == CRecon3dPoint::eInFront) {
                if (nBestCamMat != -1)
                    cout << "Multiple compatible camera matrices\n";

                nBestCamMat = nCamToRecon;
                break;
            }

        }
        if (nBestCamMat >= 0)
            return nBestCamMat;
        else {
            //Return -1 if 1st 2 at inf, -1 if 2nd 2.

            //Will be 2nd 2 as haven't returned yet
            if (nStart < 2)
                return -2;
            else
                return -1;
        }
    }

    double depth(const int nCam) const {
        if(IS_DEBUG) CHECK(nCam < 0 || nCam >= 4, "Bad camera id");
        return aRecon3dPoints[nCam].depth();
    }

    CRecon3dPoint::eReconStatus reconstructOne(const int nCamToRecon, const CCamera & Pp, const double DEPTH_THRESH) {
        CRecon3dPoint::eReconStatus reconStatus = aRecon3dPoints[nCamToRecon].reconstructOne(Pp, p1, p2, DEPTH_THRESH);
        return reconStatus; // == CRecon3dPoint::eInFront || reconStatus == CRecon3dPoint::eAtInfinity;
    }

    bool inFrontOrInf(const int nCamToRecon) const {
        if(IS_DEBUG) CHECK(nCamToRecon < 0 || nCamToRecon >= 4, "Bad camera id");
        CRecon3dPoint::eReconStatus status = aRecon3dPoints[nCamToRecon].status();
        if(IS_DEBUG) CHECK(status == CRecon3dPoint::eUninit, "Haven't reconstructed point yet");
        return status == CRecon3dPoint::eInFront || status == CRecon3dPoint::eAtInfinity;
    }

    CRecon3dPoint::eReconStatus status(const int nCamToRecon) const {
        return aRecon3dPoints[nCamToRecon].status();
    }

    const C3dPoint & point3d(const int nCamToRecon) const {
        return aRecon3dPoints[nCamToRecon].point3d();
    }

    void pp() const {
        for (int i = 0; i < NUM_CAM_MATS; i++) {
            char state = '?';
            switch (aRecon3dPoints[i].status()) {
                case CRecon3dPoint::eInFront:
                    state = 'f';
                    break;
                case CRecon3dPoint::eBehind:
                    state = 'b';
                    break;
                case CRecon3dPoint::eAtInfinity:
                    state = 'i';
                    break;
                case CRecon3dPoint::eUninit:
                    state = '?';
                    break;
                default:
                    THROW("Unhandled enum option");
            }

            cout << state;
        }
        cout << '\n';
    }

    bool isRecon(int nCam) const {
        return aRecon3dPoints[nCam].status() != CRecon3dPoint::eUninit;
    }
};

typedef CDynArray<C3dReconData> TReconData;

class CReconData : public TReconData {
    CCamera aPp[NUM_CAM_MATS];
    static const int CAM_NOT_SET = -1;
    int nCam;
    const double DEPTH_THRESH;
    static double SCALE_DEPTH_THRESH;
public:

    CReconData(int nInliers, const C3x3MatModel & E, const double DEPTH_THRESH) : nCam(CAM_NOT_SET), DEPTH_THRESH(SCALE_DEPTH_THRESH*DEPTH_THRESH) {
        reserve(nInliers);
        getCamsFromE(E.asDouble9(), aPp);
    }

    bool chooseCam() {
        CFixedArray<int, NUM_CAM_MATS> aCounts(NUM_CAM_MATS, 0);

        //Try reconstructing points
        //const int nPointsToTest = min<int>(size(), 10);
        int nPointsAtInfinity = 0, n1st2AtInf = 0, n2nd2AtInf = 0;
        int nCamTryFirst = 0;
        for (int nPointTested = 0; nPointTested < size(); nPointTested++) {
            int nChosenCam = (*this)[nPointTested].tryReconstruct(aPp, nCamTryFirst, DEPTH_THRESH);
            if (nChosenCam >= 0) {
                aCounts[nChosenCam]++;

                if (nCamTryFirst != nChosenCam)
                    nCamTryFirst = aCounts.maxIdx();
                else if (nPointTested - nPointsAtInfinity > 15) {
                    int i = 0;
                    for (; i < NUM_CAM_MATS; i++) {
                        if (i != nCamTryFirst)
                            if (2 * aCounts[i] > aCounts[nCamTryFirst])
                                break;
                    }
                    if (i == NUM_CAM_MATS)
                        //Clear winner
                        break;
                }
            } else {
                nPointsAtInfinity++;
                if (nChosenCam == -1) {
                    //1st 2 at infinity
                    nCamTryFirst = 0;
                    n1st2AtInf++;
                } else {
                    //2nd 2 at inf
                    nCamTryFirst = 2;
                    n2nd2AtInf++;
                }
            }
        }

        nCam = aCounts.maxIdx();
        int nCamCount = aCounts[nCam];

        if (nPointsAtInfinity > 10 || nCamCount < 5)
            cout << nPointsAtInfinity << " points at infinity, best cam got " << nCamCount << " votes\n";

        if ((nCamCount <= 3 && nPointsAtInfinity > 10) || (nPointsAtInfinity > 6 * nCamCount)) {
            cout << "Possible pure rotation detected. Still handling as a good pose\n";
            cout << "Re-choosing cam based on points at inf: ";

            int nBestCamCount = 0;
            for (int i = 0; i < NUM_CAM_MATS; i++) {
                int nCamCountIncludingInf = aCounts[i] + (i < 2 ? n1st2AtInf : n2nd2AtInf);
                if (nCamCountIncludingInf > nBestCamCount) {
                    nBestCamCount = nCamCountIncludingInf;
                    nCam = i;
                }
                cout << i << endl;
            }
            if (nBestCamCount < 6) {
                cout << "Insufficient inliers in front or at inf found\n";
                return false;
            }
            //return CSLAMLocMatch::ePureRotation;
        }

        int i = 0;
        for (; i < NUM_CAM_MATS; i++) {
            if (i != nCamTryFirst)
                if (2 * aCounts[i] > nCamCount) {
                    cout << "Best cam got " << nCamCount << " votes, ";
                    cout << "Second best cam got " << aCounts[i] << " votes\n";
                    pp();
                    cout << "No cam is clearly better\n";

                    return false; //rubbish position
                }
        }

        //cout << "Chosen cam " << nCam << ", now reconstructing other points\n";

        for (iterator p = begin(); p < end(); p++) {
            if (!p->isRecon(nCam))
                p->reconstructOne(nCam, aPp[nCam], DEPTH_THRESH);
        }

        return true;
    }

    void choose2dPointsInFrontOrInf(CPointVec2d & aTestPoints1, CPointVec2d & aTestPoints2, CDynArray<double> & aDepths) const {
        //cout << "Finding points in front or at inf.\n";

        for (const_iterator p = begin(); p < end(); p++) {
            if (p->inFrontOrInf(nCam)) {
                aTestPoints1.push_back(p->point2d1());
                aTestPoints2.push_back(p->point2d2());
                double dDepth = p->depth(nCam);
                if (dDepth > 0)
                    aDepths.push_back(log(dDepth));
            }
        }
    }

    inline C3dRotation getR() const {
        return aPp[nCam].rotation();
    }

    inline C3dPoint getT() const {
        return aPp[nCam].translation();
    }

    void setNewCam(const CCamera & cam) {
        aPp[nCam] = cam;
        for (iterator p = begin(); p < end(); p++) {
            p->reconstructOne(nCam, cam, DEPTH_THRESH);
        }
    }

    //Also detects pure rotation and stationary cam.

    C3dPointCollection * get3dPointCollection(CSLAMLocMatch::eStructFoundType & eStatus) {
        int nCount = 0, nAtInf = 0;
        for (iterator p = begin(); p < end(); p++) {
            if (p->status(nCam) == CRecon3dPoint::eInFront)
                nCount++;
            else if (p->status(nCam) == CRecon3dPoint::eAtInfinity)
                nAtInf++;
        }

        if (nAtInf > 10)
            cout << nAtInf << " points at inf after refinement, " << nCount << " in front\n";

        {
            static double propPointsAtInf = 0.12;
            //These weights become important if we're ever stationary for long periods...
            propPointsAtInf = 0.95 * propPointsAtInf + (0.05 * nAtInf) / (nAtInf + nCount);

            if (propPointsAtInf > 0.15)
                SCALE_DEPTH_THRESH *= 1.015;
            else if (propPointsAtInf < 0.1)
                SCALE_DEPTH_THRESH *= 0.985;

            PERIODIC(100, cout << propPointsAtInf << " = prob points at inf, ";
                    cout << SCALE_DEPTH_THRESH << "=SCALE_DEPTH_THRESH " << DEPTH_THRESH << "=DEPTH_THRESH" << endl);
        }

        if (nCount < 3 && nAtInf < 10 * nCount) //can't do anything with less than this
        {
            cout << "Too few 3d points\n";
            eStatus = CSLAMLocMatch::eNoMatch;
            return 0;
        } else if (nAtInf >= 10 * nCount) {
            cout << "Pure rotation detected\n";
            eStatus = CSLAMLocMatch::ePureRotation;
            return 0;
        }

        C3dPointCollection * pv3dPoints = new C3dPointCollection(nCount);

        for (iterator p = begin(); p < end(); p++)
            if (p->status(nCam) == CRecon3dPoint::eInFront) {
                pv3dPoints->insert(p->point3d(nCam), p->l1(), p->l2());
            }

        eStatus = CSLAMLocMatch::e3dStructFound;
        return pv3dPoints;
    }

    void pp() const {
        int nCountInFront = 0;
        for (const_iterator p = begin(); p < end(); p++)
            if (p->status(nCam) == CRecon3dPoint::eInFront)
                nCountInFront++;
        if (nCountInFront * 2 < size()) {
            cout << nCountInFront << " of " << size() << "in front\n";
            for (const_iterator p = begin(); p < end(); p++)
                p->pp();
        }
    }
};

double CReconData::SCALE_DEPTH_THRESH = 1;

C3dPoints * CSLAMLocMatch::getStructureFromRANSACInliers(CBoWMap & map, CPointVec2d & pointsCam1, CPointVec2d & pointsCam2, C3x3MatModel & E, CMask & mask, CSLAMLocMatch::eStructFoundType & eMotionType, CMxLockLater & mxLockLater) {
    const int nPointsTotal = mask.size();
    if(IS_DEBUG) CHECK(nPointsTotal != pCorr->size(), "Size mismatch")
            const int nInliers = mask.countInliers();

    const CBoWSLAMParams & BOWSLAMPARAMS = map.BOWSLAMPARAMS;

    /* Depth ~ 1/angle so if angle < k*(feature localisation error (rads)) then depth uncertainty high so assume at infinity.
     * angle < k*(feature localisation error (rads)) === depth >  1/(k*(feature localisation error (rads)))
     */

    double DEPTH_THRESH = map.BOWSLAMPARAMS.Im.getCamCalibrationMat().focalLength() / (BOWSLAMPARAMS.DEPTH_THRESH * BOWSLAMPARAMS.Corner.CORNER_LOCALISATION_SD());
    static bool bLoggedDT = false;
    if (!bLoggedDT) {
        cout << "DEPTH_THRESH=" << DEPTH_THRESH << endl;
        bLoggedDT = true;
    }

    CReconData vReconData(nInliers, E, DEPTH_THRESH);

    //Put inliers into a structure
    for (int i = 0; i < nPointsTotal; i++) {
        if (mask[i]) {
            vReconData.push_back(C3dReconData(pointsCam1[i], pointsCam2[i], (*pCorr)[i].Location1(), (*pCorr)[i].Location2()));
        }
    }

    if (!vReconData.chooseCam()) {
        cout << "Failed to choose cam\n";
        return 0;
    } else {
        C3dPoint T = vReconData.getT();
        C3dRotation R = vReconData.getR();
        C3dPoint refinedT = T;
        C3dRotation refinedR = R;

        double dAngle = R.angle();
        if (dAngle > 0.2) {
            cout << R << "=R\n";
            cout << "T: " << T << endl;
            if (dAngle > 1.5) {
                cout << "Probable error choosing rotation\n";
                vReconData.pp();
            }
        }

        CPointVec2d aTestPoints1, aTestPoints2;
        CDynArray<double> aDepths;
        aTestPoints1.reserve(nInliers);
        aTestPoints2.reserve(nInliers);
        aDepths.reserve(nInliers);

        vReconData.choose2dPointsInFrontOrInf(aTestPoints1, aTestPoints2, aDepths);

        TOTAL_RANSAC_INLIERS += aTestPoints1.size(); //Maximise num of inliers that are actually in front
        TOTAL_RANSAC_INLIERS -= (nInliers - aTestPoints1.size())*3; //bad if outliers are getting through

        const CRobustLMforEParams robustRefinementParams(BOWSLAMPARAMS.RefineRT.ROBUST_COST,
                BOWSLAMPARAMS.RefineRT.ROBUST_COST_THRESH,
                BOWSLAMPARAMS.RefineRT.ROBUST_COST_SCALE,
                BOWSLAMPARAMS.RefineRT.ROBUST_COST_CONDITION_SCALE,
                BOWSLAMPARAMS.RefineRT.ROBUST_COST_CONDITION_THRESH);

        double dWellConditioned = refineRT_LM_SampsonError_robust(aTestPoints1, aTestPoints2, refinedR, refinedT,
                robustRefinementParams, BOWSLAMPARAMS.RefineRT.VERBOSE);

        double dTErrSq = (T - refinedT).sum_square();
        double dRErr = diff(R, refinedR);

        if (dRErr > 0.1 || dTErrSq > 0.1) {
            cout << "Substantial change in R or T with refinement\n";
            cout << R << " -> " << refinedR << endl;
            cout << T << " -> " << refinedT << endl;
        }

        if ((pLoc1->time() > pLoc2->time()) == (refinedT.getZ() < 0)) {
            cout << "Possible error: camera moving backwards??\nT:" << refinedT << "\nR:" << refinedR << endl;
            vReconData.pp();
        }

        if (dTErrSq != 0 && dRErr != 0) {
            vReconData.setNewCam(refinedR | refinedT);
        }

        double dPointDepthMean = 0, dPointDepthVar = 0;
        if (map.BOWSLAMPARAMS.ConstrainScale.CONSTRAIN_SCALE && map.BOWSLAMPARAMS.BOWSpeedo.NUM_OBJECTS == 0) {
            //aDepths.getMeanVar(dPointDepthMean, dPointDepthVar);
            grubbsInliers(aDepths, 0.01, dPointDepthMean, dPointDepthVar);
            //const double dPointDepthSD = sqrt(dPointDepthVar);
            cout << "Mean log depth: " << dPointDepthMean << " s.d. " << sqrt(dPointDepthVar) << " Actual: " << exp(dPointDepthMean) << endl;

            dPointDepthMean -= map.BOWSLAMPARAMS.ConstrainScale.BASELINE_TO_DEPTH_LNPARAM1;
            dPointDepthVar += map.BOWSLAMPARAMS.ConstrainScale.BASELINE_TO_DEPTH_LNPARAM2;
        } else {
            REPEAT(1, cout << "USING SCORE\n");
            dPointDepthVar = -2; // Use SCORE or something
        }

        CRelPoseSD relPoseSD(aTestPoints1, aTestPoints2, map.BOWSLAMPARAMS.Im.getCamCalibrationMat(), BOWSLAMPARAMS.Corner.CORNER_LOCALISATION_SD(), dPointDepthMean, dPointDepthVar);

        const C3dPointCollection * pv3dPoints = vReconData.get3dPointCollection(eMotionType);
        if(IS_DEBUG) CHECK((bool)pv3dPoints != (eMotionType == e3dStructFound), "Success does not match presence of structure");

        if (eMotionType == ePureRotation && refinedR.angle() < 0.2)
            eMotionType = eSamePlace;

        C3dPoints * p3dPoints = 0;

        bool bUseThisLink = (eMotionType == e3dStructFound || (eMotionType == ePureRotation && BOWSLAMPARAMS.LinkSelection.ALLOW_ZERO_VELOCITY_LINKS));

        if (bUseThisLink)
            p3dPoints = new C3dPoints(&pv3dPoints, this, refinedR, refinedT, relPoseSD, dWellConditioned, map, mxLockLater);

        return p3dPoints;
    }
}

C3dPoints * CSLAMLocMatch::getStructure2(CPointVec2d & pointsCam1, CPointVec2d & pointsCam2, CInlierProbs & adArrLikelihood, CPointIdentifiers & pointIds, bool bNearby, CBoWMap & map, CRunGuiAp &gui, const int nSpareCores, CMxLockLater & mxLockLater, CSLAMLocMatch::eStructFoundType & eMotionType) {
    //CTSOut cout;

    cout << "\nLinking " << Loc1()->time() << " to " << Loc2()->time() << "...";

    int nCorrespondences = pCorr->size();
    const CBoWSLAMParams & BOWSLAMPARAMS = map.BOWSLAMPARAMS;
    if(IS_DEBUG) CHECK((int) pointsCam1.size() != nCorrespondences || (int) pointsCam2.size() != nCorrespondences, "CSLAMLocMatch::getStructure: Different numbers of points?!");

    //const int MIN_INLIERS = bNearby ? BOWSLAMPARAMS.MIN_INLIERS_NEARBY : BOWSLAMPARAMS.MIN_INLIERS_DISTANT;
    const int MIN_INLIERS = bNearby ? getMinInliers(adArrLikelihood, 0.01 * (double)BOWSLAMPARAMS.MIN_INLIERS_NEARBY, 5) : getMinInliers(adArrLikelihood, 0.01 * (double)BOWSLAMPARAMS.MIN_INLIERS_DISTANT, 5);

    if (nCorrespondences < MIN_INLIERS)
        return 0; // cout << "Insufficient correspondences\n";

    C3dRotation R;
    C3dPoint T;

    CMask abInliers(nCorrespondences);

    int nInliers = 0;

    C3dPoints * p3dPoints = 0;

    eMotionType = CSLAMLocMatch::eNoMatch;

    try {
        C3x3MatModel Et;
        //cout << "Assuming upright\n";
        nInliers = getE((const T2dPoints &) pointsCam1, (const T2dPoints &) pointsCam2, adArrLikelihood, pointIds, BOWSLAMPARAMS.RANSAC, Et, abInliers, -1, BOWSLAMPARAMS.Im.getCamCalibrationMat().focalLength(), 1 + nSpareCores);
        
            
        CHECK(nInliers != abInliers.countInliers(), "Inlier count failed--RANSAC returned the wrong number of inliers?");
        
        C3x3MatModel E;
        if(nInliers > 5)
        {
            //Now do least-squares nonlinear refinement of solution:
            C3dPoint T;
            C3dRotation R;
            bool bSuccess = RTFromE(abInliers, Et, (const T2dPoints &) pointsCam1, (const T2dPoints &) pointsCam2, pointIds, R, T, BOWSLAMPARAMS.RANSAC.E_INLIER_THRESH_PX / BOWSLAMPARAMS.Im.getCamCalibrationMat().focalLength());
            if(!bSuccess)
            {
                cout << "Pure rotation detected, not setting T=(0,0,1)" << endl;
                T=C3dPoint(0,0,1);
            }

            abInliers.pp();

            CRefineEOnRTManifold::refineLSOnMask((const T2dPoints &) pointsCam1, (const T2dPoints &) pointsCam2, R, T, abInliers);

            cout << "Relative orientation: " << R << endl;
            cout << "Translation direction: " << T << endl; 
            abInliers.pp();

            Eigen::Matrix3d E_eigen;
            makeE(R,T,E_eigen);

            //Then go back to BoWSLAM's own resolution of which points are in front...

            //TODO get rid of this transpose hack...
            for(int r=0;r<3;r++)
                for(int c=0;c<3;c++)
                    E(r,c) = E_eigen(c,r);
        }
        else
        {
            //TODO get rid of this transpose hack...
            for(int r=0;r<3;r++)
                for(int c=0;c<3;c++)
                    E(r,c) = Et(c,r);
        }

        
        if (nInliers > MIN_INLIERS) {
            p3dPoints = getStructureFromRANSACInliers(map, pointsCam1, pointsCam2, E, abInliers, eMotionType, mxLockLater);
        } else {
            cout << "Insufficient inliers: " << nInliers << "/" << MIN_INLIERS << " (" << nCorrespondences << " total)" << endl;
        }
        
        if (BOWSLAMPARAMS.Output.OUTPUT_CORR) {
            static CvPtr<IplImage> pIm1(cvCreateImage(BOWSLAMPARAMS.Im.SIZE(), IPL_DEPTH_8U, BOWSLAMPARAMS.Im.IM_CHANNELS));
            static CvPtr<IplImage> pIm2(cvCreateImage(BOWSLAMPARAMS.Im.SIZE(), IPL_DEPTH_8U, BOWSLAMPARAMS.Im.IM_CHANNELS));

            map.pImSource->loadImage(Loc1()->time(), pIm1);
            map.pImSource->loadImage(Loc2()->time(), pIm2);

            static CvPtr<IplImage> pIm1RGB(cvCreateImage(BOWSLAMPARAMS.Im.SIZE(), IPL_DEPTH_8U, 3));
            static CvPtr<IplImage> pIm2RGB(cvCreateImage(BOWSLAMPARAMS.Im.SIZE(), IPL_DEPTH_8U, 3));

            cvCvtColor(pIm1, pIm1RGB, CV_GRAY2RGB);
            cvCvtColor(pIm2, pIm2RGB, CV_GRAY2RGB);

            static IplImage * pOut = 0;

            markCorrespondences(pCorr, abInliers, pIm1RGB, pIm2RGB, &pOut);

            const char * CORR_DIR = "corr";
            REPEAT(1,
            if (boost::filesystem::exists(CORR_DIR))
                    boost::filesystem::remove_all(CORR_DIR);
                    boost::filesystem::create_directory(CORR_DIR));

            char szFN[100];

            sprintf_s(szFN, 100, "%s/corr%d-%d.png", CORR_DIR, Loc1()->time(), Loc2()->time());

            cvSaveImage(szFN, pOut);
        }
    } catch (CException pEx) {
        cout << pEx.GetErrorMessage() << endl;
        //failed to find 8 inliers
    }

    if (eMotionType == eNoMatch) {
        //We've tried a few times to find a good F and have failed
        cout << "Failed to find a good R+Tdir\n";
    }

    if(IS_DEBUG) CHECK((eMotionType == e3dStructFound || (eMotionType == ePureRotation && BOWSLAMPARAMS.LinkSelection.ALLOW_ZERO_VELOCITY_LINKS)) == (0 == p3dPoints), "Success does not match presence of 3d structure");
    return p3dPoints;
}


double run(const bool bOutput, const char* szConfigFile, const char * szOtherConfig /*usually 0, used for tuning */) {
    CRunGuiAp gui(bOutput);

    char szLogfileName[1000];
    gui.getLogFilename(szLogfileName);

    redirectCout redir(szLogfileName, bOutput);

    cout << setprecision(4) << "Started--Build " __DATE__ " " __TIME__ "\n";
#ifdef _DEBUG
    cout << "Debug checks enabled (build with -O or undefine _DEBUG to disable)\n";
#else
    cout << "Debug checks disabled (build with -O0 or -D_DEBUG to enable)\n";
#endif

    //if (!test()) return 0;

    //LoadSettingsFromCfg(szConfigFile);
    CBoWSLAMParams BOWSLAMPARAMS(0, 0);
    config cfg(szConfigFile);
    BOWSLAMPARAMS.init(&cfg);

    if (szOtherConfig) {
        config cfg2(szOtherConfig);
        BOWSLAMPARAMS.init(&cfg2);

        if (BOWSLAMPARAMS.Im.RD2.isInit() && BOWSLAMPARAMS.Im.RD4.isInit()) {
            double dRD2 = BOWSLAMPARAMS.Im.RD2;
            double dRD4 = BOWSLAMPARAMS.Im.RD4;

            double adRD[2] = {dRD2, dRD4};
            const_cast<CCamCalibMatrix &> (BOWSLAMPARAMS.Im.getCamCalibrationMat()).setRDCoeffs(adRD, 2);
        }
    }

    BOWSLAMPARAMS.BOW.DescriptorBinning.RADIUS = BOWSLAMPARAMS.PatchDescriptor.radius();
    BOWSLAMPARAMS.MULTI_RUNS = !bOutput; //Todo: Same for tune
    BOWSLAMPARAMS.TUNE_PARAMS = (szOtherConfig != 0);

    cvInitFont(&font, CV_FONT_HERSHEY_PLAIN, 1.2, 1.2, 0, 2);
    cvInitFont(&font_small, CV_FONT_HERSHEY_PLAIN, 0.7, 0.7);

    double dScore = 0;

    try {
        CBoWSLAM slam(BOWSLAMPARAMS, gui);
        dScore = slam.score();
    } catch (double T) {
        dScore = T;
    } catch (CException pEx) {
        cout << pEx.GetErrorMessage() << endl;
        dScore = -1;
    }

    cout << "About to finish main, score=" << dScore << "...\n";

    cout << "Parameter use summary:" << endl;
    BOWSLAMPARAMS.printUseSummary();

    BOWSLAMPARAMS.printCfgFile();

    ofstream saveParams("../params/ParamDoc.html");
    BOWSLAMPARAMS.printHtmlDoc(saveParams);
    saveParams.close();

    cout << "cd ~/workspace/BoWSLAM/maps/" << gui.outputDir() << endl;

    return dScore;
}
 
/* void noBoostWarnings() {
   boost::system::generic_category.name();
    boost::system::system_category.name();
#ifndef _DEBUG
    boost::system::errno_ecat.name();
    boost::system::native_ecat.name();
    boost::system::posix_category.name();
#endif
 
} */

int tune(const char * szFolder, const char * szOtherParams);
void testRansac();
void testCorners();
void testCorners2();
void testEigenSpeed();
void testNorms();

int main(int argc, char* argv[]) {
    //cout << sizeof(CBoWSLAMParams) << endl;
    intLookup::Setup();

    if (argc < 2 || argc > 3) {
        cout << "\nNo config file specified\nUsage: " << argv[0] << " /path/to/config/file.cfg\n";
        cout << "Remember to specify an image source in the config file (Im.ImageDir.IMAGE_DIR=\"/path/to/directory/containing/images\")" << endl;
        return 1;
    }

    if (boost::filesystem::is_directory(argv[1])) {
        cout << "Starting tuning in directory " << argv[1] << endl;
        return tune(argv[1], "config/untunedParams.cfg") ? 0 : 1;
    }

    int nSuccess = 0, nRuns = (argc > 2) ? atoi(argv[2]) : 1;

    for (int nRun = 0; nRun < nRuns; nRun++)
        nSuccess += (run(nRuns == 1, argv[1], 0) ? 0 : 1);

    return nSuccess;
}

pragma_warning(pop)
