/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "util/exception.h"
#include "util/dynArray.h"
#include "util/random.h"
#include <algorithm>

#include "alignPoints.h"

//enum eAlignmentMethods {eUmeyama, eUmeyamaRANSAC, e1dAlign, e1dAlignMedian, eAlign2dImTo3d};
// Only holds copies, no cleanup
class C3dPointMatchCopyVector : public T3dPointMatchVector
{
    //static CvRNG rng;
public:
    inline void randomSubset(T3dPointMatchVector & randSample) const
    {
        int sample_size = randSample.size(), count = (int)size();
        int i,j,k;
        const int max_random_iters = 1000;
        for(i = 0; i < sample_size; i++ )
        {
            for(k = max_random_iters; k>0; k-- )
            {
                int idxi = CRandom::Uniform(count);
                const C3dPointMatch & pRand3dCorr = *(begin() + idxi);
                for(j = 0; j < i; j++ )
                {
                    if(randSample[j] == pRand3dCorr ) //not quite ==, fp err allowed...TODO
                        break;
                }
                if( j == i ) //loop completed
                {
                    randSample[i] = pRand3dCorr;
                    break;
                }
            }
            if(IS_DEBUG) CHECK(k <= 0, "randomSubset: ERROR: Could not find sample_size points in max_random_iters iterations");
        }

        if(IS_DEBUG) CHECK( i < sample_size, "randomSubset: ERROR: Could not find sample_size points?!" ) ;
    };
};

//Owns point matches
class C3dPointMatchVector : public C3dPointMatchCopyVector
{
public:
    /*~C3dPointMatchVector()
    {
        for(iterator ppPM=begin(); ppPM<end(); ppPM++)
            delete *ppPM;
    };*/
};

void testForPlanarPoints(const T3dPointMatchVector & vPointMatches);

bool getTandR(const T3dPointMatchVector &vPointMatches, C3dRotation &R, C3dPoint &T, double &c);

//Returns # inliers or 0
/*int getTandR_RANSAC(const C3dPointMatchVector &vPointMatches, C3dRotation &R, C3dPoint &T, double &c);

bool getTandR_1d_Mean(const T3dPointMatchVector &vPointMatches, const C3dRotation &R0a, const C3dPoint &T_0a_0, const C3dPoint &T_ab_0_dir, double &c);

//Use median of good values as robust stat of translation
int getTandR_1d_Median(const T3dPointMatchVector &vPointMatches, const C3dRotation &R0a, const C3dPoint &T_0a_0, const C3dPoint &T_ab_0_dir, double &c);*/

int getTandR_1d_Mean(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c);
//Use median of good values as robust stat of translation
int getTandR_1d_Median(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c);
//Use median of good values as robust stat of translation
int getTandR_1d_Ransac(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, const double dRANSAC1dPropThresh, const int nRANSAC1dMinInliers, double &c, double & dVarconst, const bool bVerbose);
//Mean of inverse depths
int getTandR_1d_MeanID(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c);
//mean (and var) of logs
int getTandR_1d_MeanLogs(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c, double & dVar);

int getTandR_1d_RANSACLogs(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, const double dRANSAC1dPropThresh, const int nRANSAC1dMinInliers, double &c, double & dVar);

int saveScales(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, const char * fileName, const char * fileNameInliers, bool bLogs);

int getDG(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double & dMean, double & dVar, bool bVerbose);

double percentile(int nPercentile, CDynArray<double> & aNumbers);
double median(CDynArray<double> & aNumbers);
int RANSAC1d(const CDynArray<double> & aNumbers, const double p, const int nMinInliers, bool bQuiet, CMask &);

void grubbsInliers(CDynArray<double> & scales, double alpha, double & dMean, double & dVar);

