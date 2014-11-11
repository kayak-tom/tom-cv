#include "inlierCounter.h"
#include <fstream>

#include "util/stats.h"
#include "util/random.h"
#include "time/SpeedTest.h"

#include <Eigen/Core>

using namespace std;

inline bool CEssentialMatInlierCounter::isInlier_old(const C3x3MatModel & f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const double dThreshold_sq_use) const {
    double d0, d1, s0, s1;
    volatile const double m0x = p1.getX();
    volatile const double m0y = p1.getY();
    const double m1x = p2.getX();
    const double m1y = p2.getY();
    volatile double a = f[0] * m0x + f[1] * m0y + f[2]; //(a,b,c) = E*p1
    volatile double b = f[3] * m0x + f[4] * m0y + f[5];
    volatile double c = f[6] * m0x + f[7] * m0y + f[8];
    s1 = a * a + b * b;
    d1 = m1x * a + m1y * b + c; //p2 * (a,b,c)
    double d1_squared = d1 * d1;

    DEBUGONLY(if (pdResidual) *pdResidual = d1_squared / s1);

    if (EXPECT(d1_squared < dThreshold_sq_use * s1, 0)) {
        if (pdResidual) *pdResidual = d1_squared / s1;

        a = f[0] * m1x + f[3] * m1y + f[6]; //(a,b,c) = p2*E
        b = f[1] * m1x + f[4] * m1y + f[7];
        c = f[2] * m1x + f[5] * m1y + f[8];
        s0 = a * a + b * b;
        d0 = m0x * a + m0y * b + c;
        double d0_squared = d0 * d0;
        if (d0_squared < dThreshold_sq_use * s0) {
            if (pdResidual)
                *pdResidual += d0_squared / s0;

            return true;
        }
    }
    return false;
}

inline bool CEssentialMatInlierCounter::isInlier_vectoriseDouble(const Eigen::Matrix<double, 4, 2 > & f, const Eigen::Vector4d & fcol3, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const double dThreshold_sq_use) const {
    Eigen::Vector2d m0(p1.getX(), p1.getY());
    Eigen::Vector2d m1(p2.getX(), p2.getY());

    //volatile const double m0x = p1.getX();
    //volatile const double m0y = p1.getY();

    //const double m1x = p2.getX();
    //const double m1y = p2.getY();
    /*
            volatile double a = f[0] * m0x + f[1] * m0y + f[2]; //(a,b,c) = E*p1
            volatile double b = f[3] * m0x + f[4] * m0y + f[5];
            volatile double c = f[6] * m0x + f[7] * m0y + f[8];*/

    Eigen::Vector4d abc = f * m0 + fcol3;

    //s1 = a * a + b * b;
    const double s1 = abc.block < 2, 1 > (0, 0).transpose() * abc.block < 2, 1 > (0, 0);

    //d1 = m1x * a + m1y * b + c; //p2 * (a,b,c)
    const double d1 = m1.transpose() * abc.block < 2, 1 > (0, 0) + abc(2);
    const double d1_squared = sqr(d1);

    DEBUGONLY(if (pdResidual) *pdResidual = d1_squared / s1);

    if (EXPECT(d1_squared < dThreshold_sq_use * s1, 0)) {
        if (pdResidual) *pdResidual = d1_squared / s1;

        /*double a = f[0] * m1x + f[3] * m1y + f[6]; //(a,b,c) = p2*E
        double b = f[1] * m1x + f[4] * m1y + f[7];
        double c = f[2] * m1x + f[5] * m1y + f[8];*/
        abc = f * m1 + fcol3;

        double s0 = abc.block < 2, 1 > (0, 0).transpose() * abc.block < 2, 1 > (0, 0);
        double d0 = m0.transpose() * abc.block < 2, 1 > (0, 0) + abc(2);
        double d0_squared = sqr(d0);
        if (d0_squared < dThreshold_sq_use * s0) {
            if (pdResidual)
                *pdResidual += d0_squared / s0;

            return true;
        }
    }
    return false;
}

inline bool CEssentialMatInlierCounter::isInlier_vectorise(const Eigen::Matrix3f & f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const float dThreshold_sq_use) const {
    Eigen::Vector3f m0((float) p1.getX(), (float) p1.getY(), 1);
    Eigen::Vector3f m1((float) p2.getX(), (float) p2.getY(), 1);

    /*volatile const float m0x = (float)p1.getX();
    volatile const float m0y = (float)p1.getY();

    const float m1x = (float)p2.getX();
    const float m1y = (float)p2.getY();*/
    /*
            volatile double a = f[0] * m0x + f[1] * m0y + f[2]; //(a,b,c) = E*p1
            volatile double b = f[3] * m0x + f[4] * m0y + f[5];
            volatile double c = f[6] * m0x + f[7] * m0y + f[8];*/

    Eigen::Vector3f abc = f * m0;

    //s1 = a * a + b * b;
    const float s1 = abc.block < 2, 1 > (0, 0).transpose() * abc.block < 2, 1 > (0, 0);

    //d1 = m1x * a + m1y * b + c; //p2 * (a,b,c)
    const float d1 = m1.transpose() * abc;
    const float d1_squared = sqr(d1);

    DEBUGONLY(if (pdResidual) *pdResidual = d1_squared / s1);

    if (EXPECT(d1_squared < dThreshold_sq_use * s1, 0)) {
        if (pdResidual) *pdResidual = d1_squared / s1;

        /*float a = f[0] * m1x + f[3] * m1y + f[6]; //(a,b,c) = p2*E
        float b = f[1] * m1x + f[4] * m1y + f[7];
        float c = f[2] * m1x + f[5] * m1y + f[8];*/
        abc = f*m1;

        float s0 = abc.block < 2, 1 > (0, 0).transpose() * abc.block < 2, 1 > (0, 0);
        float d0 = m0.transpose() * abc;
        float d0_squared = d0 * d0;
        if (d0_squared < dThreshold_sq_use * s0) {
            if (pdResidual)
                *pdResidual += d0_squared / s0;

            return true;
        }
    }
    return false;
}

template<typename TFloat>
inline bool CEssentialMatInlierCounter::isInlier(const TFloat * f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const TFloat dThreshold_sq_use) {
    TFloat d1, s0, s1;
    volatile const TFloat m0x = p1.getX();
    volatile const TFloat m0y = p1.getY();
    const TFloat m1x = p2.getX();
    const TFloat m1y = p2.getY();
    volatile TFloat a = f[0] * m0x + f[1] * m0y + f[2]; //(a,b,c) = E*p1
    volatile TFloat b = f[3] * m0x + f[4] * m0y + f[5];
    volatile TFloat c = f[6] * m0x + f[7] * m0y + f[8];
    s1 = a * a + b * b;
    d1 = m1x * a + m1y * b + c; //p2 * (a,b,c)
    TFloat d1_squared = d1 * d1; //Sampson's err sq. numerator

    DEBUGONLY(if (pdResidual) *pdResidual = d1_squared / s1);

    if (EXPECT(d1_squared < dThreshold_sq_use * s1, 0)) {
        if (pdResidual) *pdResidual = d1_squared / s1;

        //TB Sep 2011: simplify: d1_squared == d0_squared
        a = f[0] * m1x + f[3] * m1y + f[6]; //(a,b,c) = p2*E^T
        b = f[1] * m1x + f[4] * m1y + f[7];
        //TB c = f[2] * m1x + f[5] * m1y + f[8];
        s0 = a * a + b * b;
        //TB d0 = m0x * a + m0y * b + c;
        //TB TFloat d0_squared = d0 * d0;
        if (d1_squared < dThreshold_sq_use * s0) {
            if (pdResidual && s0 < s1)
                *pdResidual = d1_squared / s0;

            return true;
        }
    }
    return false;
}
template<typename TFloat>
inline bool CEssentialMatInlierCounter::isInlier_SE(const TFloat * f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const TFloat dThreshold_sq_use) {
    TFloat d1, s0, s1;
    volatile const TFloat m0x = p1.getX();
    volatile const TFloat m0y = p1.getY();
    const TFloat m1x = p2.getX();
    const TFloat m1y = p2.getY();
    volatile TFloat a = f[0] * m0x + f[1] * m0y + f[2]; //(a,b,c) = E*p1
    volatile TFloat b = f[3] * m0x + f[4] * m0y + f[5];
    volatile TFloat c = f[6] * m0x + f[7] * m0y + f[8];
    s1 = a * a + b * b;
    d1 = m1x * a + m1y * b + c; //p2 * (a,b,c)
    TFloat d1_squared = d1 * d1; //Sampson's err sq. numerator

    a = f[0] * m1x + f[3] * m1y + f[6]; //(a,b,c) = p2*E^T
    b = f[1] * m1x + f[4] * m1y + f[7];
    s0 = a * a + b * b;
    
    if (pdResidual)
    {
        TFloat dErr = d1_squared / (s0+s1);
        *pdResidual = dErr;
        return dErr < dThreshold_sq_use;
    }
    else {
        return d1_squared < dThreshold_sq_use * (s0+s1);
    }
}

void CEssentialMatInlierCounter::timeInlierCounter() const {
    CStopWatch s;
    s.startTimer();
    int test = 0;
    C3x3MatModel m;
    Eigen::Matrix3f mf;
    mf.setRandom();
    Eigen::Matrix3d md;
    md.setRandom();
    /*for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                    mf(i,j) = m(i,j);*/
    float af[9];
    double ad[9];

    Eigen::Matrix<double, 4, 2 > f;
    Eigen::Vector4d fcol3;
    f.setRandom();
    fcol3.setRandom();

    for (int i = 0; i < 9; i++)
        af[i] = (float) (ad[i] = CRandom::Uniform());
    
    ad[0]=(double)af[0];//noop to supress warnings

    for (int nPoint = 0; nPoint < 1000000000; nPoint++) {
        CSimple2dPoint p1(CRandom::Uniform(), CRandom::Uniform());
        CSimple2dPoint p2(CRandom::Uniform(), CRandom::Uniform());
        //const bool bStatus = isInlier_vectorise(mf, p1, p2, 0, 0.123f); //0.1+microsecs
        //const bool bStatus = isInlier_vectoriseDouble(f, fcol3, p1, p2, 0, 0.123); //0.077 micros
        //const bool bStatus = isInlier<float>(af, p1, p2, 0, 0.123f); //0.051 micros
        const bool bStatus = isInlier_SE<double>(ad, p1, p2, 0, 0.123f); //0.044 micros
        test += bStatus;
    }
    s.stopTimer();
    cout << test << " time=" << s.getElapsedTime() << " nanoseconds" << endl;
}

void CEssentialMatInlierCounter::countInliers(const CModel & model_in, CRansacTerminator * pTerminator, const int nThread, const int nBGC, CMask & mask, double * adResiduals, int & nInlierCount, double dThresholdScale) const {
    //timeInlierCounter();
    const C3x3MatModel & model = dynamic_cast<const C3x3MatModel &> (model_in);
    pTerminator->reset();
    int nPoints = (int) ap.size();
    double dThreshold_sq_use = dThreshold_sq * dThresholdScale; //Allow threshhold to be adapted when doing topdown refinement
    if(IS_DEBUG) CHECK(mask.size() != nPoints, "Mask size must match number of points")
    if(IS_DEBUG) CHECK(mask.countInliers() != 0, "Mask should be zero'd before counting inliers")

    nInlierCount = 0;

    //For each point call isInlier
    T2dPoints::const_iterator pPoint1 = ap.begin();
    T2dPoints::const_iterator pPoint2 = ap_prime.begin();
    double * pdResiduals = adResiduals;
    CMask::iterator pMask = mask.begin();

    const bool TEST_TERMINATOR = false;

    ofstream * pInlierFile = 0;
    static CStats aSuccesses, aNumTrials;
    static int nBestInlierCount = 0;

    if (TEST_TERMINATOR) {
        pInlierFile = new ofstream("ransacTerminationLog.tsv", ios_base::app);
        if (pTerminator->bestGoodCount() == 0) {
            aSuccesses.add((nBestInlierCount > 20) ? 1 : 0);
            nBestInlierCount = 0;

            *pInlierFile << "RESET\n";
            *pInlierFile << "Av num trials: ";
            aNumTrials.writeTSVdata(*pInlierFile);
            *pInlierFile << endl << "Success rate: ";
            aSuccesses.writeTSVdata(*pInlierFile);
            *pInlierFile << endl;
        }
    }

    ofstream & inlierFile = *pInlierFile;
    bool bTerminated = false;

    for (int nPoint = nPoints; nPoint > 0; nPoint--, pPoint1++, pPoint2++, pMask++) {
        const bool bStatus = isInlier_SE<double>(model.asDouble9(), *pPoint1, *pPoint2, adResiduals ? pdResiduals++ : 0, dThreshold_sq_use);

        *pMask = bStatus;
        
        if (bStatus) {
            nInlierCount++;
            
            if(false){ // Just verify that everything is the correct way around...
                Eigen::Matrix3d E;
                for(int r=0;r<3;r++)
                    for(int c=0;c<3;c++)
                        E(r,c) = model(r,c);
                
                Eigen::Vector3d x(pPoint1->getX(), pPoint1->getY(), 1);
                Eigen::Vector3d xp(pPoint2->getX(), pPoint2->getY(), 1);
                cout << xp.transpose() * E * x << endl;
            }
        }

        if (!TEST_TERMINATOR) {
            if (pTerminator->observeStatus(nThread, bStatus)) {
                break;
            }
        } else {
            if (nPoint == 1) break; //at the end already

            if (pTerminator->observeStatus(nThread, bStatus)) {
                bTerminated = true;
                const int nTerminateAt = nPoints - (nPoint - 1);

                inlierFile << "T\t" << nInlierCount << '\t' << nTerminateAt << endl;
                aNumTrials.add(nTerminateAt);
                break;
            }
        }
    }

    if (TEST_TERMINATOR) {
        if (!bTerminated) {
            inlierFile << "E\t" << nInlierCount << '\t' << nPoints << endl;
            aNumTrials.add(nPoints);
        }

        if (nInlierCount > nBestInlierCount)
            nBestInlierCount = nInlierCount;

        inlierFile.close();
        delete pInlierFile;
    }
}

inline bool CHomographyInlierCounter::isInlier(const C3x3MatModel & H, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const double dThreshold_sq_use) const {
    volatile const double m0x = p1.getX();
    volatile const double m0y = p1.getY();
    const double m1x = p2.getX();
    const double m1y = p2.getY();
    volatile double predictX = H[0] * m0x + H[1] * m0y + H[2];
    volatile double predictY = H[3] * m0x + H[4] * m0y + H[5];
    volatile double predictW = H[6] * m0x + H[7] * m0y + H[8];
    volatile double predictW_inv = 1.0 / predictW;
    predictX *= predictW_inv;
    predictY *= predictW_inv;

    volatile double reproj_err = sqr(m1x - predictX) + sqr(m1y - predictY);

    if (pdResidual) *pdResidual = reproj_err;

    return dThreshold_sq_use > reproj_err;
    /*sum_i((x'i-(h11*xi + h12*yi + h13)/(h31*xi + h32*yi + h33))2+
              (y'i-(h21*xi + h22*yi + h23)/(h31*xi + h32*yi + h33))2) -> min*/
}

void CHomographyInlierCounter::countInliers(const CModel & model_in, CRansacTerminator * pTerminator, const int nThread, const int nBGC, CMask & mask, double * adResiduals, int & nInlierCount, double dThresholdScale) const {
    const C3x3MatModel & model = dynamic_cast<const C3x3MatModel &> (model_in);
    pTerminator->reset();
    int nPoints = (int) ap.size();
    if(IS_DEBUG) CHECK(mask.size() != nPoints, "Mask size must match number of points")
    if(IS_DEBUG) CHECK(mask.countInliers() != 0, "Mask should be zero'd before counting inliers")

    nInlierCount = 0;

    double dThreshold_sq_use = dThreshold_sq * dThresholdScale; //Allow threshhold to be adapted when doing topdown refinement

    //For each point call isInlier
    for (int nPoint = 0; nPoint < nPoints; nPoint++) {
        bool bStatus = isInlier(model, ap[nPoint], ap_prime[nPoint], adResiduals ? adResiduals + nPoint : 0, dThreshold_sq_use);

        if (bStatus) {
            mask[nPoint] = 1;
            nInlierCount++;
        }

        if (pTerminator->observeStatus(nThread, bStatus)) return;
    }
}
