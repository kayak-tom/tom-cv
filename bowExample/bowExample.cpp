/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/* ********** Bag-of-Words Example **********8
 *
 * Pass a path to an image directory as a command line argument
 *
 * bowTest function adds up to 1000 images to a BoW database, rebuilds the word bag, and lists the 4 best matches to each image (the first of which should be itself!)
 *
 * imageMatching function finds sets of correspondences between consecutive pairs of images, and also lists matches.
 */
#define NOMINMAX

#include "bow/bagOfWords.h"
#include <boost/filesystem.hpp>
#include <vector>
#include <iomanip>
#include <opencv2/opencv.hpp>
#include <fstream>
#include "geom/geom.h"
#include <boost/smart_ptr.hpp>

#include "geom/geom.h"
#include "bow_imu_params.h"
#include "imageSource/imageSourceFromDir.h"

#include "params/config.h"
#include "image/convert_OpenCV.h"

#include "featureExtract/featureExtractor.h"

#include "time/SpeedTest.h"

using namespace std;

void imageMatching(CBoWIMUParams & BOWIMUPARAMS)
{
	//Load images from a directory (or whatever) specified in config file
    boost::scoped_ptr<CImageSource> pImageLoader( CImageSourceFromDir::makeImageSource(BOWIMUPARAMS.Im) );
	boost::scoped_ptr<CFeatureExtractor> pFeatureExtractor ( CFeatureExtractor::makeFeatureExtractor(BOWIMUPARAMS.Im, BOWIMUPARAMS.Corner, BOWIMUPARAMS.PatchDescriptor, BOWIMUPARAMS.DescriptorSetClustering));

	//Create empty BoW DB with params from config file
    CBoW bow(BOWIMUPARAMS.BOW);

	ofstream logStream;
	logStream.open( "log.txt", fstream::out );

	//cvNamedWindow("Bow1" );
	//cvNamedWindow("Bow2" );

	int iter = 0;
    IplImage * pIm = cvCreateImage(BOWIMUPARAMS.Im.SIZE(), IPL_DEPTH_8U, BOWIMUPARAMS.Im.IM_CHANNELS);
	IplImage * pIm2 = cvCreateImage(BOWIMUPARAMS.Im.SIZE(), IPL_DEPTH_8U, BOWIMUPARAMS.Im.IM_CHANNELS);

    while(iter<15 && pImageLoader->loadImage(iter, pIm) && pImageLoader->loadImage(iter+1, pIm2) ) //reloading--inefficient
	{
		//get the descriptors
		CDescriptorSet * pDesc = pFeatureExtractor->getDescriptors(pIm);
		CDescriptorSet * pDesc2 = pFeatureExtractor->getDescriptors(pIm2);

		//mark the features on image 1.
		markDescriptors( pIm, pDesc );

		//mark the features on image 2
		markDescriptors( pIm2, pDesc2 );

		//Add image to bag-of-words
		bow.addImage(&pDesc, iter); //...and add it to the Bag-of-Words database with id 0

        if(iter % 10 == 0)
            bow.recreateDictionary();

		//work out the correspondences
		const CBoWCorrespondences * pCorr = bow.getCorrespondences(iter, pDesc2, BOWIMUPARAMS.BOWMatching);
		CMask abInliers( pCorr->size() ); //initialise inliers (these are worked out in computeAngles(...)

		//If you have a calibration matrix you can use geom/geom.h functions and BaySAC to compute the essential matrix
		//bool anglesOk = computeAngles( pCorr, abInliers, logStream, BOWIMUPARAMS );

		cout << "Num correspondences " << pCorr->size() << endl;

		//now mark inlier correspondences
		int nInliers = 0;
		CvScalar col = CV_RGB(0,255,0);
		for (int nPoint = 0; nPoint < pCorr->size(); nPoint++)
		{
			//CvScalar col = (abInliers[nPoint]) ? CV_RGB(0,255,0) : CV_RGB(255,0,0);
			CLocation pointLeft  = (*pCorr)[nPoint].Location1();
			CLocation pointRight = (*pCorr)[nPoint].Location2();

			CvPoint cvPointLeft  = cvPoint(pointLeft.x() , pointLeft.y()); //invert points in the y-axis
			CvPoint cvPointRight = cvPoint(pointRight.x(), pointRight.y());//invert points in the y-axis

			if ( abInliers[nPoint] )
			{
				cvCircle( pIm , cvPointRight, 2, CV_RGB( 255, 255, 0 ), 1 );
				cvCircle( pIm2, cvPointLeft , 2, CV_RGB( 255, 255, 0 ), 1 );

				//only draw lines for good points
				cvLine( pIm2, cvPointLeft, cvPointRight, col, 1 );
				nInliers++;
			}
		}

		TBoWMatchVector * pMatches = bow.getMatches(pDesc2);
        cout << "Top 10 matches with latest image:\n";
        int nMatch = 1;
        for(TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end() && nMatch<=10; pMatch++, nMatch++)
            cout << nMatch << ": " << pMatch->id() << " strength = " << pMatch->MatchStrength() << endl;
        delete pMatches;
        
        const int type = (BOWIMUPARAMS.Im.IM_CHANNELS == 3) ? CV_8UC3 : CV_8UC1;
        
        cv::Mat im1 (pIm->height, pIm->width, type, pIm->origin);
        cv::imshow("Bow1", im1);
        cv::Mat im2 (pIm2->height, pIm2->width, type, pIm2->origin);
        cv::imshow("Bow2", im2);

		//cvShowImage( "Bow1", pIm  );
		//cvShowImage( "Bow2", pIm2 );

		cv::waitKey( 1 );

		CDescriptorSet::deleteDS((const CDescriptorSet **)&pDesc2);
		delete pCorr;

        iter++; 
	}
	cvReleaseImage( &pIm  ); 
	cvReleaseImage( &pIm2 ); 
}

//Load every image into BoW database. For every match list 4 closest matches.
void bowTest(CBoWIMUParams & BOWIMUPARAMS)
{
    boost::scoped_ptr<CImageSource> pImageLoader( CImageSourceFromDir::makeImageSource(BOWIMUPARAMS.Im) );
	boost::scoped_ptr<CFeatureExtractor> pFeatureExtractor ( CFeatureExtractor::makeFeatureExtractor(BOWIMUPARAMS.Im, BOWIMUPARAMS.Corner, BOWIMUPARAMS.PatchDescriptor, BOWIMUPARAMS.DescriptorSetClustering));

    CBoW bow(BOWIMUPARAMS.BOW);

    CvPtr<IplImage> pIm ( cvCreateImage(BOWIMUPARAMS.Im.SIZE(), IPL_DEPTH_8U, BOWIMUPARAMS.Im.IM_CHANNELS));

    CStopWatch s; s.startTimer();

	//Load each image in turn and add to BoW
    const int CA = 1;
    BOWIMUPARAMS.BOW.BOWClustering.DESCRIPTORS_PER_WORD = 20/CA;
    BOWIMUPARAMS.BOW.DescriptorBinning.LOWER_BOUND = std::min<int>(CA,4);

    double dClusterTime = 0;
    int nId = 0;
	for(; pImageLoader->loadImage(nId, pIm); nId++ )
	{
		//get descriptors from the image
		CDescriptorSet * pDesc = pFeatureExtractor->getDescriptors(pIm);
		cout << "Adding " << nId << ", " << pDesc->Count() << " descriptors\n";

		/*Can display images with few descriptors:
		if(pDesc->Count() < 50)
		{
			markDescriptors(pIm, pDesc);
			cvShowImage("Features", pIm);
			cvWaitKey(0);
		}*/

		if(nId == CA*1020)
		{
		    bow.recreateDictionary(); //Create codebook once have added CA images
		}

		bow.addImage(&pDesc, nId);
	}
    const int NUM_IMAGES = nId;
    s.stopTimer();
    double dIndexTime = s.getElapsedTime();

    /*s.startTimer();
    bow.recreateDictionary(); //Create codebook only once have added all images
    s.stopTimer();
    double dClusterTime = s.getElapsedTime();*/

    s.startTimer();
    double dScore = 0;
    for(nId = 0; nId < NUM_IMAGES; nId++)
    {
		//if((nId > 2295 && nId < 3000) || (nId > 5431 && nId < 5436)) continue;

    	boost::scoped_ptr<TBoWMatchVector> pMatches ( bow.getMatches(nId, 4) );
    	cout << nId << ": ";
    	int nCorrectMatches=0;
    	for(TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end(); pMatch++)
    	{
    		int nMatchId = pMatch->id();
    		if(nMatchId / 4 == nId/4) nCorrectMatches++;
    		cout << nMatchId << ", ";
    	}
    	cout << "Score=" << nCorrectMatches << endl;
    	dScore += nCorrectMatches;
    }
    dScore /= NUM_IMAGES;
	cout << "Total score=" << dScore << endl;
	s.stopTimer();
	double dQueryTime = s.getElapsedTime();

	cout << "Index time=" << dIndexTime << " per-image: " << dIndexTime/NUM_IMAGES << endl;
	cout << "Cluster time=" << dClusterTime  << endl;
	cout << "Query time=" << dQueryTime << " per-query: " << dQueryTime/NUM_IMAGES << endl;


}

int main(int argc, char* argv[])
{
	CHECK(argc != 2, "Pass 1 argument: path to config file");

	config cfg(argv[1]);
	CBoWIMUParams BOWIMUPARAMS(0, 0);
	BOWIMUPARAMS.init(&cfg); // Initialise parameters

	/* Can also set parameters in code:
	BOWIMUPARAMS.BOW.DescriptorBinning.RADIUS = BOWIMUPARAMS.PatchDescriptor.radius();
	BOWIMUPARAMS.BOW.BOWClustering.LEVELS = 3;
	BOWIMUPARAMS.PatchDescriptor.Patch.ORIENT = true;
	BOWIMUPARAMS.BOW.COMP_METHOD = CBOWParams::eVectorDistFast; //Slightly faster than eNisterDist
	 */
	BOWIMUPARAMS.BOW.BOWClustering.RECLUSTER_FREQUENCY = -1; //Uses specifies when to cluster with recreateDictionary();
	BOWIMUPARAMS.BOW.BOWClustering.CLUSTER_IN_SEPERATE_THREAD = false;

	if(false)
		imageMatching(BOWIMUPARAMS);
	else
		bowTest(BOWIMUPARAMS);

    //BOWIMUPARAMS.printUseSummary();

	return 0;
}
