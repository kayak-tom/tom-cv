/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include <algorithm>
#include "geom.h"
#include "pointAlignment.h"
#include <fstream>

using namespace std;

//CvRNG C3dPointMatchCopyVector::rng = cvRNG(0);

double median(CDynArray<double> & aNumbers)
{
    return percentile(50, aNumbers);
}

double percentile(int nPercentile, CDynArray<double> & aNumbers)
{
    int nSize = aNumbers.size();

    if(IS_DEBUG) CHECK(nSize == 0, "percentile: 0 elements supplied");

    int medianIdx = (nSize*nPercentile)/100;

    CDynArray<double>::iterator pPercentile = aNumbers.begin() + medianIdx;
    nth_element(aNumbers.begin(), pPercentile, aNumbers.end());
    double c = *pPercentile;

    return c;
}

int RANSAC1d(CDynArray<double> & aNumbers, const double p, const int nMinInliers, bool bQuiet, double & dMean, double & dVar, const bool bLogs)
{
    int nMaxInliers = 0;//, nMinInliers = (int)(aNumbers.size()*dMinInlierRate);

    /*CDynArray<double> aNumbers;
    aNumbers_in.copyInto(aNumbers);*/

    std::sort(aNumbers.begin(), aNumbers.end());

    int middle = aNumbers.size()/2, nBestIdx = 0, nIdx = 0;

    if(!bQuiet)
    {
        for(CDynArray<double>::iterator pd = aNumbers.begin(); pd < aNumbers.end(); pd++)
            cout << *pd << ' ';
        cout << "\nmed=" << aNumbers[aNumbers.size()/2] << ' ';
    }

    CDynArray<double>::iterator pd_lo = aNumbers.begin();
    CDynArray<double>::iterator pd_hi = aNumbers.begin();
    CDynArray<double>::iterator pInliersBegin = aNumbers.end(), pInliersEnd = aNumbers.end();
    for(CDynArray<double>::iterator pd = aNumbers.begin(); pd < aNumbers.end(); pd++, nIdx++)
    {
        double d = *pd;
        double d_lo=d*(1-p);
        double d_hi=d*(1+p);
        if(d_lo>d_hi)
            std::swap(d_lo, d_hi);
        if(bLogs)
        {
            d_lo=d-p;
            d_hi=d+p;
        }
        while(*pd_lo < d_lo) pd_lo++; //this is in the range
        while(pd_hi < aNumbers.end() && *pd_hi < d_hi) pd_hi++; //this is not
        int inliers = (int)(pd_hi-pd_lo);
        if (inliers > nMaxInliers)
        {
            nBestIdx = nIdx;
            nMaxInliers = inliers;
            pInliersBegin = pd_lo;
            pInliersEnd = pd_hi;
        }
        else if (inliers == nMaxInliers)
        {
            //If closer to the middle choose this one
            if(abs(middle - nBestIdx) > abs(middle - nIdx))
            {
                nBestIdx = nIdx;
                nMaxInliers = inliers;
                pInliersBegin = pd_lo;
                pInliersEnd = pd_hi;
            }
        }
    }
    if(IS_DEBUG) CHECK(nMaxInliers == 0 || pInliersBegin == pInliersEnd, "1d RANSAC failed");

    bQuiet || cout << "\n1d RANSAC inliers=" << nMaxInliers << "/" << aNumbers.size();

    if(nMaxInliers < nMinInliers)
    {
        return 0;
    }

    dMean = 0;
    for(CDynArray<double>::const_iterator pInlier = pInliersBegin; pInlier < pInliersEnd; pInlier++)
        dMean += *pInlier;

    dMean /= nMaxInliers;
    bQuiet || cout << " ransac mean = " << dMean;

    if(nMaxInliers >= 4)
    {
        dVar = 0;
        for(CDynArray<double>::const_iterator pInlier = pInliersBegin; pInlier < pInliersEnd; pInlier++)
            dVar += sqr(dMean - (*pInlier));

        dVar /= (nMaxInliers*(nMaxInliers-3));
    }
    else
    {
        //The T-variable has infinite variance otherwise
        dVar = 1e+4;
    }
    bQuiet || cout << " ransac var = " << dVar << endl;

    //And set inlier/outlier statuses in mask:
    //mask.pp();

    /*CDynArray<double>::const_iterator pNumbersInOrder = aNumbers_in.begin(), pInliersBack = pInliersEnd-1;
    for(int i=0; i<(int)mask.size(); i++)
        if(mask[i])
        {
            mask[i] = (*pNumbersInOrder >= *pInliersBegin
                    && *pNumbersInOrder <= *pInliersBack );

            pNumbersInOrder++;
        }
    //mask.pp();

    DEBUGONLY(int nMaskCount =  mask.countInliers());
    if(IS_DEBUG) CHECK(nMaxInliers != nMaskCount, "Updating mask in 1d ransac failed");*/

    return nMaxInliers;
}

void testForPlanarPoints(const T3dPointMatchVector & vPointMatches)
{
    if((int)(vPointMatches.size()) > 4){
        const C3dPoint & p1first = vPointMatches[0].p1();
        const C3dPoint & p2first = vPointMatches[0].p2();
        const C3dPoint & p1second = vPointMatches[1].p1();
        const C3dPoint & p2second = vPointMatches[1].p2();
        C3dPoint norm1 = crossproduct(p1second - p1first, p2second - p1first);
        norm1.normalise();
        C3dPoint norm2 = crossproduct(vPointMatches[2].p1() - p2first, vPointMatches[2].p2() - p2first);
        norm2.normalise();
        T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin() + 3;
        for(;pPointPair != vPointMatches.end();pPointPair++){
            C3dPoint planeLine1 = pPointPair->p1() - p1first;
            planeLine1.normalise();
            if(dotproduct(planeLine1, norm1) > 0.05)
                break;

            C3dPoint planeLine2 = pPointPair->p2() - p2first;
            planeLine2.normalise();
            if(dotproduct(planeLine2, norm2) > 0.05)
                break;

        }
        if(pPointPair == vPointMatches.end()){
            cout << "Warning: Planar 3d points\n";
        }
    }

}

bool getTandR(const T3dPointMatchVector &vPointMatches, C3dRotation &R, C3dPoint &T, double &c)
{
    THROW( "Needs newmat")

    /* OK but needs NEWMAT or re-writing
    DEBUGONLY(testForPlanarPoints(vPointMatches));

    C3dPoint mean1, mean2;

    const double sizeInv = 1.0 / vPointMatches.size();

    for(T3dPointMatchVector::const_iterator ppPointPair = vPointMatches.begin();ppPointPair != vPointMatches.end();ppPointPair++){
        C3dPointMatch *pPointPair = *ppPointPair;
        const C3dPoint & p1 = *(pPointPair->p1());
        const C3dPoint & p2 = *(pPointPair->p2());
        mean1 += p1;
        mean2 += p2;
    }
    mean1 *= sizeInv;
    mean2 *= sizeInv;

    double var1 = 0, var2 = 0;
    Matrix Sigma(3, 3);
    Sigma = 0;
    for(T3dPointMatchVector::const_iterator ppPointPair = vPointMatches.begin();ppPointPair != vPointMatches.end();ppPointPair++){
        C3dPointMatch *pPointPair = *ppPointPair;
        const C3dPoint & centredPoint1 = (*(pPointPair->p1()) - mean1);
        const C3dPoint & centredPoint2 = (*(pPointPair->p2()) - mean2);
        var1 += centredPoint1.sum_square();
        var2 += centredPoint2.sum_square();
        Sigma += centredPoint2.asCV() * centredPoint1.asCV().t();
    }
    var1 *= sizeInv;
    var2 *= sizeInv;
    Sigma *= sizeInv;
    DiagonalMatrix S(I3);
    if(Sigma.determinant() <= 0){
        return false;
        cout << "Warning: getTandR: large errors (neg determinant)...";
        S(3, 3) = -1;
    }
    DiagonalMatrix D;
    Matrix U;
    Matrix V;
    SVD(Sigma, D, U, V);
    const Matrix & R_temp = U * S * V.t();

    if(!isRotMat(R_temp)) return false;

    R=R_temp;
    //Was c = (D * S).trace() / var1; //Done: swap to get rid of 1 div
    //c = 1. / c;
    double dTrace = (D * S).trace();
    if(0 == dTrace) return false;
    c = var1 / dTrace;
    T = mean2 * c - R * mean1;
    //PRINTMAT(mean2)
    //PRINTMAT(R * mean1)
    return true;*/
}

double getErr(const C3dPointMatch & PM, const C3dRotation &R_hyp, const C3dPoint &T_hyp, double c_hyp)
{
    const C3dPoint & C1 = PM.p1();
    const C3dPoint & C2 = PM.p2();

    //ColumnVector p1calc=(R_hyp * *pC1 + T_hyp/c_hyp);
    //double err=(p1calc - (*pC2)/c_hyp).sum_square();

    //c*points in cam2 frame = R* points in cam1 frame + T
    C3dPoint p1calc=(R_hyp * C1 + T_hyp);
    double err=(p1calc - C2*c_hyp).sum_square();

    return err;
}

/*double getErr(const C3dPointMatch * pPM, const Matrix &R_hyp, const ColumnVector &T_hyp, double c_hyp)
{
    const ColumnVector * pC1 = pPM->p1();
    const ColumnVector * pC2 = pPM->p2();
    ColumnVector p1calc=(R_hyp.i() * *pC2)/c_hyp - T_hyp;
    double err=(p1calc - *pC1).sum_square();
    return err;
}*/

//Return inlier count, or 0 to fail
int getTandR_RANSAC(const C3dPointMatchVector &vPointMatches, C3dRotation &R, C3dPoint &T, double &c)
{
    int count = (int)vPointMatches.size();
    if(IS_DEBUG) CHECK(count<4, "getTandR_RANSAC: Less than 4 points");

    if(count==4)
    {
        if(!getTandR(vPointMatches, R, T, c))
        {
            cout << "getTandR_RANSAC failed: 4 points don't give a good transformation\n";
            return 0;
        }
        return 4;
    }

    int best_good_count=0;

    const int sample_size = 4;
//    const double INLIER_PRIOR_PROB = 0.5;
    int max_samples = 100;//RANSACMaxIters(count, sample_size, INLIER_PRIOR_PROB, 0.98);
    cout << "Use RANSACMaxIters\n";

    T3dPointMatchVector vSelected(count);
    T3dPointMatchVector vBest;

    T3dPointMatchVector randSample(sample_size);

    for(int sample_count = 0; sample_count < max_samples; sample_count++ )
    {
        // choose random <sample_size> (=4) points
        vPointMatches.randomSubset(randSample);

        // find motion hyp.
        C3dPoint T_hyp;
        double c_hyp;
        C3dRotation R_hyp;

        if(!getTandR(randSample, R_hyp, T_hyp, c_hyp)) continue;
        //c*points in cam2 frame = R* points in cam1 frame + T

        const double MAX_SCALE_RANGE = 1000;//Max scale allowed relative to 1st scale
        if(c_hyp > MAX_SCALE_RANGE || c_hyp < 1./MAX_SCALE_RANGE)
        {
            cout << "Scale: " << c_hyp << " out of range\n";
            continue;
        }

        //Now choose a good threshhold based on how well these 4 points fit: (is this invalid when these 4 are in an impossible config??)
        /*double dMaxErr = 0;
        for(int n3dCorr=0; n3dCorr<sample_size; n3dCorr++)
        {
            double err=getErr(randSample[n3dCorr],R_hyp, T_hyp, c_hyp);
            if(dMaxErr<err) dMaxErr = err;
        }

        const double threshold = 0.6*dMaxErr; //todo: probably variable
        cout << threshold << "=RANSAC threshhold\n";
        */
        const double threshold = sqr(5); //RANSAC_3PT_THRESH;
        cout << "missing param: RANSAC_3PT_THRESH\n";

        // for each pair of 3d points see if this R+S+T translates 1 to 2 and thus find
        // the number of in-liers.
        int good_count = 0;
        for(int i = 0; i < count; i++ )
        {
            double err=getErr(vPointMatches[i], R_hyp, T_hyp, c_hyp);

            if(err < threshold)
            {
//                cout << err << "=err\n";
                vSelected[i] = vPointMatches[i];
                good_count++;
            }
            else
                vSelected[i].setZero();
        }

        if( good_count > max<int>( best_good_count, sample_size ) ) //3 points mean nothing much
        {
            // update the current best inlier set
            vBest.clear();
            for (T3dPointMatchVector::iterator pCorr = vSelected.begin(); pCorr < vSelected.end(); pCorr++)
            {
                if(!pCorr->zero()) vBest.push_back(*pCorr);
            }

            R=R_hyp; T=T_hyp; c=c_hyp;
            best_good_count = good_count;

            //const double p = 0.99;
            //max_samples = BoW_DisjointRANSACUpdateNumIters( p, //might want to adjust this to cope with num of points being greater than can actually be
            //    (double)(count - good_count)/count, 7, max_samples );

            cout << "Use RANSACMaxIters\n";
            //int new_max_samples = RANSACMaxItersOld(count, sample_size, best_good_count, p);
            //if (new_max_samples<max_samples) max_samples = new_max_samples;
        }
    }

    if( best_good_count <= sample_size)
    {
        cout << "getTandR_RANSAC failed: No consensus set found\n";
        return 0;
    }

    cout << best_good_count << " inliers found, ";
    cout << max_samples << " max_samples\n";

    if( best_good_count > sample_size )
    {
        if(!getTandR(vBest, R, T, c))
        {
            cout << "getTandR_RANSAC failed: Error calculating transformation from all points\n";
            return 0;
        }
    }
    return best_good_count;
}

inline bool resolveScale1d(const double ALIGNMENT_THETA_MAX_ERR, const C3dPointMatch & pMatch, double &c)
{
    static double TAN_MAX_THETA_ERR_SQ = -1, ALIGNMENT_THETA_MAX_ERR_OLD = -1;
    if(ALIGNMENT_THETA_MAX_ERR_OLD != ALIGNMENT_THETA_MAX_ERR)
    {
        ALIGNMENT_THETA_MAX_ERR_OLD = ALIGNMENT_THETA_MAX_ERR;
        TAN_MAX_THETA_ERR_SQ = sqr(tan(ALIGNMENT_THETA_MAX_ERR));
    }

    // Least-squares c  = pThis_a_lastTime_scaleAmbig . pOther_a / | pThis_a_lastTime_scaleAmbig |^2
    double dotProd = dotproduct(pMatch.p2(), pMatch.p1());
    
    if(dotProd <= 0) 
        return false;

    c = dotProd / pMatch.p2().sum_square(); //No sqrt needed here

    double tan_thetaErr_sq = crossproduct(pMatch.p2(), pMatch.p1()).sum_square();// / sqr(dotprod);

    // return fail if there's clearly a large angular error (????)
    return (tan_thetaErr_sq >= 0 && tan_thetaErr_sq < sqr(dotProd)*TAN_MAX_THETA_ERR_SQ);
}

/*inline bool getTandR_1d(const C3dPointMatch & pMatch, double &c)
{
    static const double MAX_THETA_ERR = 0.2;
    static const double TAN_MAX_THETA_ERR_SQ = tan(MAX_THETA_ERR)*tan(MAX_THETA_ERR);

    const C3dPoint & p1_0 = pMatch.p1(); // point at time a
    const C3dPoint & p2_a_scaleAmbig = pMatch.p2(); // point at time b

    const C3dPoint & p1_a = R0a*(p1_0 - T_0a_0); // Is this calc correct? Changing to a - really improves results!! should be a T_a0_0 I guess
    // Now scale p2_a_scaleAmbig so it is somewhere near p1_a in the direction T_ab_a_dir
    // p2_a_scaleAmbig*c=p2_a    p2_a_scaleAmbig*c = p1_a + T_ab_a_dir*c because T_ab_a_dir is the baseline
    // (p2_a_scaleAmbig - T_ab_0_dir)*c = p1_a
    const C3dPoint & p2_a_lastTime_scaleAmbig = p2_a_scaleAmbig - T_ab_0_dir; //Less outliers with a - here--posibly a scale ambiguity??
    // Least-squares c  = p2_a_lastTime_scaleAmbig . p1_a / | p2_a_lastTime_scaleAmbig |^2

    c = dotproduct(p2_a_lastTime_scaleAmbig, p1_a) / p2_a_lastTime_scaleAmbig.sum_square();

    / * / Alternative: |p2_a_scaleAmbig - T_ab_0_dir|*c = |p1_a|
    c = sqrt(p1_a.sum_square() / p2_a_lastTime_scaleAmbig.sum_square());* /

    double tan_thetaErr_sq = crossproduct(p2_a_lastTime_scaleAmbig, p1_a).sum_square() / sqr(dotproduct(p2_a_lastTime_scaleAmbig, p1_a));

    / *cout << c << "=c\n";
    if(tan_thetaErr_sq > TAN_MAX_THETA_ERR_SQ)
    {
        cout << tan_thetaErr_sq << "=err_sq, max = "<< TAN_MAX_THETA_ERR_SQ << "\n";
        PRINTMAT(p2_a_lastTime_scaleAmbig)
        PRINTMAT(p1_a)
    }* /

    // return fail if there's clearly a large angular error (????)
    return (tan_thetaErr_sq >= 0 && tan_thetaErr_sq < TAN_MAX_THETA_ERR_SQ);
}*/

int getTandR_1d_Mean(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c)
{
    int nSize = (int)vPointMatches.size();
    if(IS_DEBUG) CHECK(nSize < 1, "getTandR_1d: no points");
    //if(IS_DEBUG) CHECK(!zero((1.0 - T_ab_0_dir.sum_square())*0.01), "getTandR_1d: Bad direction--need to normalise");

    double dScale = 0;
    int nCountSuccess = 0;
    for(int i=0; i<nSize; i++)
    {
        bool bSuccess = resolveScale1d(ALIGNMENT_THETA_MAX_ERR, vPointMatches[i], c);
        if(bSuccess)
        {
            dScale += c;
            nCountSuccess++;
        }
    }

    if(nCountSuccess < 3*nSize/4)
        cout << nCountSuccess << "/" << nSize << " good points found\n";

    if(nCountSuccess>0)
        c = dScale/nCountSuccess;

    return nCountSuccess;
}

int getTandR_1d_MeanID(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c)
{
    int nSize = (int)vPointMatches.size();
    if(IS_DEBUG) CHECK(nSize < 1, "getTandR_1d_MeanID: no points");
    //if(IS_DEBUG) CHECK(!zero((1.0 - T_ab_0_dir.sum_square())*0.01), "getTandR_1d: Bad direction--need to normalise");

    double dScale = 0;
    int nCountSuccess = 0;
    for(int i=0; i<nSize; i++)
    {
        bool bSuccess = resolveScale1d(ALIGNMENT_THETA_MAX_ERR, vPointMatches[i], c);
        if(bSuccess)
        {
            dScale += 1.0/c;
            nCountSuccess++;
        }
    }

    if(nCountSuccess < 3*nSize/4)
        cout << nCountSuccess << "/" << nSize << " good points found\n";

    if(nCountSuccess>0)
        c = nCountSuccess/dScale;

    return nCountSuccess;
}

void getScaleVector(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, CDynArray<double> & scales)
{
    int nSize = (int)vPointMatches.size();
    if(IS_DEBUG) CHECK(nSize < 1, "getTandR_1d: no points");

    scales.reserve(nSize);

    double c;
    for(int i=0; i<nSize; i++)
    {
        bool bSuccess = resolveScale1d(ALIGNMENT_THETA_MAX_ERR, vPointMatches[i], c);
        if(bSuccess)
        {
            //cout << c << ' ';
            scales.push_back(c);
        }
    }

    if(scales.size() < 3*nSize/4)
        cout << scales.size() << "/" << nSize << " good points found\n";
}

//50%     60%     70%     80%     90%     95%     98%     99%     99.5%     99.8%     99.9%
enum eAlphaVals { e50, e40, e30, e20, e10, e05, e02, e01, e005, e002, e001 };

//Get from table at http://en.wikipedia.org/wiki/T-distribution (2 sided)
//alpha < 5%/12 => alpha < 0.005
double getTStat(eAlphaVals alpha, const int dof)
{
    double aadLookupT[37][11] = {{ 1.000 , 1.376 , 1.963 , 3.078 , 6.314 , 12.71 , 31.82 , 63.66 , 127.3 , 318.3 , 636.6  } ,
                            { 0.816 , 1.061 , 1.386 , 1.886 , 2.920 , 4.303 , 6.965 , 9.925 , 14.09 , 22.33 , 31.60  } ,
                            { 0.765 , 0.978 , 1.250 , 1.638 , 2.353 , 3.182 , 4.541 , 5.841 , 7.453 , 10.21 , 12.92  } ,
                            { 0.741 , 0.941 , 1.190 , 1.533 , 2.132 , 2.776 , 3.747 , 4.604 , 5.598 , 7.173 , 8.610  } ,
                            { 0.727 , 0.920 , 1.156 , 1.476 , 2.015 , 2.571 , 3.365 , 4.032 , 4.773 , 5.893 , 6.869  } ,
                            { 0.718 , 0.906 , 1.134 , 1.440 , 1.943 , 2.447 , 3.143 , 3.707 , 4.317 , 5.208 , 5.959  } ,
                            { 0.711 , 0.896 , 1.119 , 1.415 , 1.895 , 2.365 , 2.998 , 3.499 , 4.029 , 4.785 , 5.408  } ,
                            { 0.706 , 0.889 , 1.108 , 1.397 , 1.860 , 2.306 , 2.896 , 3.355 , 3.833 , 4.501 , 5.041  } ,
                            { 0.703 , 0.883 , 1.100 , 1.383 , 1.833 , 2.262 , 2.821 , 3.250 , 3.690 , 4.297 , 4.781  } ,
                            { 0.700 , 0.879 , 1.093 , 1.372 , 1.812 , 2.228 , 2.764 , 3.169 , 3.581 , 4.144 , 4.587  } ,
                            { 0.697 , 0.876 , 1.088 , 1.363 , 1.796 , 2.201 , 2.718 , 3.106 , 3.497 , 4.025 , 4.437  } ,
                            { 0.695 , 0.873 , 1.083 , 1.356 , 1.782 , 2.179 , 2.681 , 3.055 , 3.428 , 3.930 , 4.318  } ,
                            { 0.694 , 0.870 , 1.079 , 1.350 , 1.771 , 2.160 , 2.650 , 3.012 , 3.372 , 3.852 , 4.221  } ,
                            { 0.692 , 0.868 , 1.076 , 1.345 , 1.761 , 2.145 , 2.624 , 2.977 , 3.326 , 3.787 , 4.140  } ,
                            { 0.691 , 0.866 , 1.074 , 1.341 , 1.753 , 2.131 , 2.602 , 2.947 , 3.286 , 3.733 , 4.073  } ,
                            { 0.690 , 0.865 , 1.071 , 1.337 , 1.746 , 2.120 , 2.583 , 2.921 , 3.252 , 3.686 , 4.015  } ,
                            { 0.689 , 0.863 , 1.069 , 1.333 , 1.740 , 2.110 , 2.567 , 2.898 , 3.222 , 3.646 , 3.965  } ,
                            { 0.688 , 0.862 , 1.067 , 1.330 , 1.734 , 2.101 , 2.552 , 2.878 , 3.197 , 3.610 , 3.922  } ,
                            { 0.688 , 0.861 , 1.066 , 1.328 , 1.729 , 2.093 , 2.539 , 2.861 , 3.174 , 3.579 , 3.883  } ,
                            { 0.687 , 0.860 , 1.064 , 1.325 , 1.725 , 2.086 , 2.528 , 2.845 , 3.153 , 3.552 , 3.850  } ,
                            { 0.686 , 0.859 , 1.063 , 1.323 , 1.721 , 2.080 , 2.518 , 2.831 , 3.135 , 3.527 , 3.819  } ,
                            { 0.686 , 0.858 , 1.061 , 1.321 , 1.717 , 2.074 , 2.508 , 2.819 , 3.119 , 3.505 , 3.792  } ,
                            { 0.685 , 0.858 , 1.060 , 1.319 , 1.714 , 2.069 , 2.500 , 2.807 , 3.104 , 3.485 , 3.767  } ,
                            { 0.685 , 0.857 , 1.059 , 1.318 , 1.711 , 2.064 , 2.492 , 2.797 , 3.091 , 3.467 , 3.745  } ,
                            { 0.684 , 0.856 , 1.058 , 1.316 , 1.708 , 2.060 , 2.485 , 2.787 , 3.078 , 3.450 , 3.725  } ,
                            { 0.684 , 0.856 , 1.058 , 1.315 , 1.706 , 2.056 , 2.479 , 2.779 , 3.067 , 3.435 , 3.707  } ,
                            { 0.684 , 0.855 , 1.057 , 1.314 , 1.703 , 2.052 , 2.473 , 2.771 , 3.057 , 3.421 , 3.690  } ,
                            { 0.683 , 0.855 , 1.056 , 1.313 , 1.701 , 2.048 , 2.467 , 2.763 , 3.047 , 3.408 , 3.674  } ,
                            { 0.683 , 0.854 , 1.055 , 1.311 , 1.699, 2.045 , 2.462 , 2.756 , 3.038 , 3.396 , 3.659  } ,
                            { 0.683 , 0.854 , 1.055 , 1.310, 1.697 , 2.042 , 2.457 , 2.750 , 3.030 , 3.385 , 3.646  } ,
                            { 0.681 , 0.851 , 1.050 , 1.303 , 1.684 , 2.021 , 2.423 , 2.704 , 2.971 , 3.307 , 3.551  } ,
                            { 0.679 , 0.849 , 1.047 , 1.299 , 1.676 , 2.009 , 2.403 , 2.678 , 2.937 , 3.261 , 3.496  } ,
                            { 0.679 , 0.848 , 1.045 , 1.296 , 1.671 , 2.000 , 2.390 , 2.660 , 2.915 , 3.232 , 3.460  } ,
                            { 0.678 , 0.846 , 1.043 , 1.292 , 1.664 , 1.990 , 2.374 , 2.639 , 2.887 , 3.195 , 3.416  } ,
                            { 0.677 , 0.845 , 1.042 , 1.290 , 1.660 , 1.984 , 2.364 , 2.626 , 2.871 , 3.174 , 3.390  } ,
                            { 0.677 , 0.845 , 1.041 , 1.289 , 1.658 , 1.980 , 2.358 , 2.617 , 2.860 , 3.160 , 3.373  } ,
                            { 0.674 , 0.842 , 1.036 , 1.282 , 1.645 , 1.960 , 2.326 , 2.576 , 2.807 , 3.090 , 3.291 }};


    if(dof <= 30)
        return aadLookupT[dof-1][alpha];
    else
        switch (dof)
        {
        case 40:
            return aadLookupT[30][alpha];
        case 50:
            return aadLookupT[31][alpha];
        case 60:
            return aadLookupT[32][alpha];
        case 80:
            return aadLookupT[33][alpha];
        case 100:
            return aadLookupT[34][alpha];
        case 120:
            return aadLookupT[35][alpha];
        default:
            return aadLookupT[36][alpha];
        }

    /*switch(alpha)
    {
    case 0.05:

    }*/

    THROW("Missing alpha val for computing t stat")
}
//Get from table at http://en.wikipedia.org/wiki/T-distribution (2 sided)
//alpha < 5%/12 => alpha < 0.005
double getTStat(const double alpha, const int dof)
{
    if(IS_DEBUG) CHECK(dof < 1 || alpha <= 0, "getTStat: Values OOB");
    int dof_lo = 10000, dof_hi = 10000;
    if (dof<=30)
        dof_lo=dof_hi=dof;
    else if (dof<=60)
    {
        if(dof % 10 == 0)
            dof_lo=dof_hi=dof;
        else
        {
            dof_lo = (dof/10)*10;
            dof_hi=dof_lo+10;
        }
    }
    else if (dof<=120)
    {
        if(dof % 20 == 0)
            dof_lo=dof_hi=dof;
        else
        {
            dof_lo = (dof/20)*20;
            dof_hi=dof_lo+20;
        }
    }
    //Else equal and will get inf

    eAlphaVals e_lo, e_hi;
    double dWeightHi = 1;

    if(alpha < 0.001)
    {
        e_lo=e_hi = e001;
    }
    else if(alpha < 0.002)
    {
        e_lo=e001, e_hi = e002;
        dWeightHi = (alpha-0.001)*1000;
    }
    else if(alpha < 0.005)
    {
        e_lo=e002, e_hi = e005;
        dWeightHi = (alpha-0.002)*333.333;
    }
    else if(alpha < 0.01)
    {
        e_lo=e005, e_hi = e01;
        dWeightHi = (alpha-0.005)*200;
    }
    else if(alpha < 0.02)
    {
        e_lo=e01, e_hi = e02;
        dWeightHi = (alpha-0.01)*100;
    }
    else if(alpha < 0.05)
    {
        e_lo=e02, e_hi = e05;
        dWeightHi = (alpha-0.03)*33.3333;
    }
    else if(alpha < 0.1)
    {
        e_lo=e05, e_hi = e10;
        dWeightHi = (alpha-0.05)*20;
    }
    else if(alpha < 0.2)
    {
        e_lo=e10, e_hi = e20;
        dWeightHi = (alpha-0.1)*10;
    }
    else if(alpha < 0.2)
    {
        e_lo=e20, e_hi = e30;
        dWeightHi = (alpha-0.2)*10;
    }
    else if(alpha < 0.2)
    {
        e_lo=e30, e_hi = e40;
        dWeightHi = (alpha-0.3)*10;
    }
    else if(alpha < 0.2)
    {
        e_lo=e40, e_hi = e50;
        dWeightHi = (alpha-0.4)*10;
    }
    else
    {
        e_lo=e50, e_hi = e50;
        dWeightHi = 1;
    }

    double dLo = getTStat(e_lo, dof_lo);
    double dHi = getTStat(e_hi, dof_lo);

    return dLo * (1-dWeightHi) + dHi * dWeightHi;
}

//http://en.wikipedia.org/wiki/Grubbs%27_test_for_outliers
void grubbsInliers(CDynArray<double> & scales, double alpha, double & dMean, double & dVar)
{
    std::sort(scales.begin(), scales.end());

    const int MIN_GRUBBS = 3; //need >= 3, 3 seems to help though

    for (int N = scales.size(); ; N = scales.size())
    {
        scales.getMeanVar(dMean, dVar); //need to do this for return, and as var depends on mean

        if(N < MIN_GRUBBS)
        {
            return;
        }

        double dSD = sqrt(dVar);
        double t = getTStat(alpha/(2*N), N-2);

        double dGrubbsThresh = dSD  *  (N-1)*(sqrt(sqr(t) / (N * (N - 2 + sqr(t) ))));

        //Now choose the end with the most deviation from the mean and truncate as many as poss from there
        double deviationLo = dMean - scales.top();
        double deviationHi = scales.back() - dMean;

        if(deviationHi > deviationLo && deviationHi > dGrubbsThresh)
        {
            //Remove from back
            //for(CDynArray<double>::const_iterator pdHi = scales.end()-1; *pdHi-dMean > dGrubbsThresh; pdHi--, Nnew--);
            //Update mean and var (have to recompute var cos mean has changed)
            int nRemove = 0;
            do
            {
                nRemove ++;
                scales.resize(N-nRemove);
                //actual val of dGrubbsThresh would be falling
            } while (scales.back() - dMean > max<double>(dGrubbsThresh, deviationLo));

            if(nRemove*8>N)
                cout << "Removed " << nRemove << " of " << N << " from back\n";

        }
        else if(deviationLo > dGrubbsThresh)
        {
            //Remove from front
            //for(CDynArray<double>::const_iterator pdLo = scales.begin(); dMean-*pdLo > dGrubbsThresh; pdLo++, Nnew--);
            int nRemove = 0;
            do
            {
                nRemove++;
            }
            while (dMean - scales[nRemove] > max<double>(dGrubbsThresh, deviationHi));

            if(nRemove*8>N)
                cout << "Removed " << nRemove << " of " << N << " from front\n";

            scales.popN(nRemove);
        } else //We're done and have a mean+var est.
            break;
    }
}

int saveScales(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, const char * fileName, const char * fileNameInliers, bool bLogs)
{
    CDynArray<double> scales;

    getScaleVector(vPointMatches, ALIGNMENT_THETA_MAX_ERR, scales);

    if(scales.size() == 0)
        return 0;

    if(bLogs)
    {
        double (*flog)(double) = &log;
        scales.apply(*flog);
    }

    ofstream all(fileName), inliers(fileNameInliers);

    std::sort(scales.begin(), scales.end());
    //cout << "All scales: ";
    scales.pp("\n", all);
    //cout << endl;
    double dMean, dVar;
    grubbsInliers(scales, 0.05, dMean, dVar);
    //cout << "Inlier scales: ";
    scales.pp("\n", inliers);

    return scales.size();
}

int getDG(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double & dMean, double & dVar, bool bVerbose)
{
    CDynArray<double> scales;

    getScaleVector(vPointMatches, ALIGNMENT_THETA_MAX_ERR, scales);

    if(scales.size() == 0)
        return 0;

    double (& pflog)(double) = log;
    scales.apply(pflog);

    if(bVerbose)
    {
        std::sort(scales.begin(), scales.end());
        cout << "All scales: "; scales.pp();
    }

    grubbsInliers(scales, 0.05, dMean, dVar);
    
    if(bVerbose)
    {
        cout << "Inlier scales: "; scales.pp();
    }

    //The set has mean dMean and variance dVar
    //The *mean* (d) has variance dVar/(N-3) because [mean - sample mean] / sqrt[sample variance/N] is a T-statistic, so has var v/v-2, so d has variance [sample variance/N]*(N-1)/(N-3) = var/N-3


    if(scales.size() >= 4)
        dVar /= scales.size()-3; //var of T-statistic
    else
        dVar = 1e+5; //infinite

    return scales.size();
}

//Use median of good values as robust stat of translation
int getTandR_1d_Median(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c)
{
    CDynArray<double> scales;

    getScaleVector(vPointMatches, ALIGNMENT_THETA_MAX_ERR, scales);

    if(scales.size() == 0)
        return 0;

    c = median(scales);

    return scales.size();
}

int getTandR_1d_Ransac(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, const double dRANSAC1dPropThresh, const int nRANSAC1dMinInliers, double &c, double & dVar, const bool bVerbose)
{
    int nSize = (int)vPointMatches.size();
    CDynArray<double> scales;     
    
    getScaleVector(vPointMatches, ALIGNMENT_THETA_MAX_ERR, scales);

    if(scales.size() < nRANSAC1dMinInliers)
        return 0;

    int nInliers = RANSAC1d(scales, dRANSAC1dPropThresh, nRANSAC1dMinInliers, !bVerbose, c, dVar, false);

    if(bVerbose)
        cout << nInliers << "/" << nSize << " good points after 1d RANSAC\n";

    return nInliers;
}

int getTandR_1d_MeanLogs(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, double &c, double & dVar)
{
    CDynArray<double> scales;

    getScaleVector(vPointMatches, ALIGNMENT_THETA_MAX_ERR, scales);

    if(scales.size() < 2)
        return 0;

    double dMean = 0;
    dVar = 0;

    for(CDynArray<double>::iterator pd = scales.begin(); pd != scales.end(); pd++)
        dMean += *pd = log(*pd);

    dMean /= scales.size();

    for(CDynArray<double>::iterator pd = scales.begin(); pd != scales.end(); pd++)
        dVar += sqr(*pd - dMean);

    dVar /= (scales.size()-1);

    c = dMean;

    return scales.size();
}

int getTandR_1d_RANSACLogs(const T3dPointMatchVector &vPointMatches, const double ALIGNMENT_THETA_MAX_ERR, const double dRANSAC1dPropThresh, const int nRANSAC1dMinInliers, double &c, double & dVar)
{
    CDynArray<double> scales;

    getScaleVector(vPointMatches, ALIGNMENT_THETA_MAX_ERR, scales);

    if(scales.size() < nRANSAC1dMinInliers)
    {
        return 0;
    }

    for(CDynArray<double>::iterator pd = scales.begin(); pd != scales.end(); pd++)
        *pd = log(*pd);

    int nInliers = RANSAC1d(scales, dRANSAC1dPropThresh, nRANSAC1dMinInliers, false, c, dVar, true);

    return nInliers;
}

/*/ Resolve scale NOT TESTED -- use MEDIAN instead
bool getTandR_1dRANSAC(const C3dPointMatchVector &vPointMatches, const Matrix &R0a, const ColumnVector &T_0a_0, const ColumnVector &T_ab_0_dir, double &c)
{
    int count = (int)vPointMatches.size();
    if(IS_DEBUG) CHECK(count<1, "getTandR_1dRANSAC: No points"); //1 for scale, 1 for trans.

    if(count<=2) //can't really do RANSAC with 2 points (?)
    {
        return getTandR_1d(vPointMatches, R0a, T_0a_0, T_ab_0_dir, c); //, "getTandR_RANSAC: these particular points don't give a good transformation");
    }

    int result = 0, best_good_count=0;

    const int sample_size = 1;
    const double INLIER_PRIOR_PROB = 0.4;
    int max_samples = RANSACMaxIters(count, sample_size, INLIER_PRIOR_PROB, 0.98);

    //Matrix R_hyp;

    T3dPointMatchVector vSelected(count);
    T3dPointMatchVector vBest;
    T3dPointMatchVector randSample(sample_size);

    for(int sample_count = 0; sample_count < max_samples; sample_count++ )
    {
        // choose random <sample_size> (=4) points
        vPointMatches.randomSubset(randSample); //could be faster...

        // find motion hyp.
        //ColumnVector T_hyp;
        double c_hyp;

        if(!getTandR_1d(randSample, R0a, T_0a_0, T_ab_0_dir, c_hyp)) continue;
        //c*points in cam2 frame = R* points in cam1 frame + T

        const double MAX_SCALE_RANGE = 30;//Max scale allowed relative to 1st scale
        if(c_hyp > MAX_SCALE_RANGE || c_hyp < -MAX_SCALE_RANGE)
        {
            cout << "Scale of " << c_hyp << " too large/small compared to original scale of 1\n";
            continue;
        }

        const double threshold = SQR(25);

        // for each pair of 3d points see if this R+S+T translates 1 to 2 and thus find
        // the number of in-liers.
        int good_count = 0;
        for( int i = 0; i < count; i++ )
        {
            double err=getErr(vPointMatches[i], R0a, T_0a_0, c_hyp);
//            double err=getErr(vPointMatches[i], R_hyp, T_hyp, c_hyp);

            if(err < threshold)
            {
//                cout << err << "=err\n";
                vSelected[i] = vPointMatches[i];
                good_count++;
            }
            else
                vSelected[i] = 0;
        }

        if( good_count > MAX( best_good_count, sample_size ) ) //3 points mean nothing much
        {
            // update the current best inlier set
            vBest.resize(0);
            for (T3dPointMatchVector::iterator ppCorr = vSelected.begin(); ppCorr < vSelected.end(); ppCorr++)
            {
                C3dPointMatch * pCorr = *ppCorr;
                if(pCorr) vBest.push_back(pCorr);
            }

            c=c_hyp;
            best_good_count = good_count;

            const double p = 0.99;
            max_samples = RANSACTryReduceIters(count, sample_size, best_good_count, max_samples, p);
        }
    }

    if(IS_DEBUG) CHECK( best_good_count <= sample_size, "getTandR_1dRANSAC: No consensus set found" );

    result = 1;

    cout << best_good_count << " inliers (1d alignment)\n";
    cout << max_samples << " max_samples\n";

    if( best_good_count > sample_size )
        if(IS_DEBUG) CHECK(!getTandR_1d(vBest, R0a, T_0a_0, T_ab_0_dir, c), "getTandR_1dRANSAC: Error calculating trasnformation from all points" );

    return true;
}*/
