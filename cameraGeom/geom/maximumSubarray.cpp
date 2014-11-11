//#include "StdAfx.h"
#include "maximumSubarray.h"
#include <util/convert.h>

void findMaximalSubarray2d(const cv::Mat & image_in, const int nMean, cv::Rect & BB_best)
{
    const int N = image_in.rows, M=image_in.cols;
    int t = 0;
    //int a[100][100];
    std::vector<int> pr(N);
    int S = -MAX_INT, s = 0, k, l, x1 = 0,x2 = 0,y1 = 0,y2 = 0,j;
 
    for( int z = 0; z < N; z++)
    {
        for(int i = 0; i < N; i++) pr[i] = 0;
 
        for(int x = z; x < M; x++)
        {
            t = 0;
            s = 1<<31;
            j = 0;
            k = 0; l = 0;
            for(int i = 0; i < N; i++)
            {
                pr[i] = pr[i] + ((int)image_in.at<uchar>(i,x) - nMean);
                t = t + pr[i];
                if( t > s)
                {
                    s = t;
                    k = i;
                    l = j;
                }
                if( t < 0 )
                {
                    t = 0;
                    j = i + 1;
                }
            }
            if( s > S)
            {
                S = s;
                x1 = x;
                y1 = k;
                x2 = z;
                y2 = l;
            }
        }
    }
 
    //cout << x1 << " " << y1 << " " << x2 << " "  << y2 << endl;
    //cout << S;

    BB_best = cv::Rect(x2,y2,x1-x2,y1-y2);
}

int percentile(const cv::Mat & M, const int nPercentile)
{
    std::vector<int> aVals; aVals.reserve(M.rows*M.cols / 8);
	for(int r=0; r<M.rows; r+=3) 
    	for(int c=0; c<M.cols; c+=3) 
	    	aVals.push_back(M.at<uchar>(r,c));

	std::sort(aVals.begin(), aVals.end());
    int nIdx = (int)((double)nPercentile*0.01*(double)aVals.size());
    if(nIdx >= (int)aVals.size() || nIdx < 0)
        throw "Bad percentile";
	
    return aVals[nIdx];
}

//Choose a better threshold by examining image. Assume about 1 on the breast
void findMaximalSubarray2d_findPercentile(const cv::Mat & image_in, const int nPercentile, cv::Rect & BB_best)
{
    int nThresh = percentile(image_in, nPercentile);

    cout << "Thresh: " << nThresh << endl;

    findMaximalSubarray2d(image_in, nThresh, BB_best);
}

//Choose a better threshold by examining image. Assume about 1 on the breast
void findMaximalSubarray2d_adaptive(const cv::Mat & image_in, cv::Rect & BB_best)
{
    double dMeanBreast = 200;
    double dMin, dMax;
    cv::minMaxLoc(image_in, &dMin, &dMax);

    cv::Mat blockAbove(image_in, cv::Range(0, image_in.rows/10));
    cv::Mat blockLeft(image_in, cv::Range::all(), cv::Range(0, 1+image_in.cols/20));
        
    double dMeanElsewhere = 0.5*(cv::mean(blockAbove)[0] + cv::mean(blockLeft)[0]);

    cout << "Max: " << dMax << endl;
    cout << "MeanElsewhere: " << dMeanElsewhere << endl;

    int nThresh = (int)(0.65*(dMeanBreast + dMeanElsewhere));

    findMaximalSubarray2d(image_in, nThresh, BB_best);
}

//Want to choose a better threshold based on 1st thresh
void findMaximalSubarray2d_2pass(const cv::Mat & image_in, cv::Rect & BB_best)
{
    const int nMean = (int)cv::mean(image_in)[0];
    findMaximalSubarray2d(image_in, nMean, BB_best);

    double dSumInside = cv::sum(image_in(BB_best))[0];
    double dSumOutside = cv::sum(image_in)[0]-dSumInside;

    double dMeanInside=dSumInside / (BB_best.area());
    double dMeanOutside=dSumOutside / (image_in.rows*image_in.cols);

    const int nNewThresh = (int)(0.5*(dMeanInside - dMeanOutside));
    findMaximalSubarray2d(image_in, nNewThresh, BB_best);
}
