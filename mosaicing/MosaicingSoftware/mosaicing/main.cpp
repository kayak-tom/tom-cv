/*! \mainpage DTA Image Mosaicer code documentation
 *
 *  \section vs Building the DTA Image Mosaicer
 *
 *  Netbeans 7.0 recommended. See README for build/run instructions.
 *
 */

#include "util/exception.h"


pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "ImageSource.h"
#include "FeatureExtractor.h"
#include "RansacPerspectiveTransformEstimator.h"
#include "TransformEngine.h"

#include "Renderer.h"
#include "FeatheredRenderer.h"
#include "MultiScaleFeatheredRenderer.h"
#include "DijkstraCutRenderer.h"

#include "ShiTomasiFeatureDetector.h"
#include "DoGBlobDetector.h"
#include "SimpleAsymFeatureMatcher.h"
#include "FastAsymFeatureMatcher.h"
#include "Descriptor.h"
#include "GRCException.h"
#include "params/config.h"
#include "params/param.h"
#include "SumSquaredErrorEvaluationFunction.h"

#include "ransac/ransacParams.h"
#include "bow/bagOfWordsParam.h"
#include "image/imageAccess.h"
#include "description/vectorDescriptor.h"
#include "description/patchParams.h"
#include "featureExtract/cornerDetector.h"
#include "featureExtract/featureExtractor.h"
#include "imageSource/imageSource.h"

#include "BFFeatureMatcher.h"
#include "BoWFeatureMatcher.h"
#include "BaySACTransformEstimator.h"
#include "bow/bagOfWords.h"
#include "time/SpeedTest.h"

#include <iostream>

WRAPPERPARAMCLASS(Mosaicing)
	PARAM(VideoSource, 0, 2, 2,"")
	//PARAMSTR(ImageDir)
	//PARAMSTR(VideoFile)
	PARAMSTR(MosaicSaveDir,"")
	PARAMB(MosaicDestination, false,"") //1=save
	PARAM(TransformType, 0, 2, 0,"")
	PARAM(LMIterations, 0, 100, 8,"")
	PARAM(LM, 0, 1, 0, "")
	PARAM(PatchDescriptorType, 0, 3, 3, "")
	PARAM(PatchDescriptorSize, 0, 3, 2, "")
	PARAM(MaxFrames, 0, 1000, 30, "")
	PARAMB(IncrementalRendering, true, "")
	PARAM(FullFrameUpdateFreq, 0, 100, 7, "")
	PARAM(SkipFrames, 0, 1, 0, "") //1=skip
	PARAM(RendererType, 0, 3, 1, "")
	PARAM(FeatherRadius, 1, 100, 16, "")
	PARAM(DijkstraScale, 1, 16, 4, "")
	PARAM(WarpMethod, 0, 1, 0, "")
	PARAM(MosaicX, 32, 10000, 1600, "")
	PARAM(MosaicY, 32, 10000, 1600, "")
	PARAM(EvaluationFunction, 0, 1, 0, "") //1=ssd
	//PARAM(MaxRansacIters, 1, 10000, 100, "")
	PARAM(MaxSearchForTransform, 0, 200, 10, "")
	PARAM(FeatureDetector, 0, 2, 0, "")
	PARAM(NumBlobScales, 1, 10, 3, "")
	//PARAM(CorrespondenceQuality, 0.1, 1, 0.8, "")
	//PARAM(NumFeatures, 10, 10000, 150, "")
	PARAMB(MarkCorrespondences, false, "")
	CHILDCLASS(BOW, "")
	CHILDCLASS(BOWMatching, "")
	CHILDCLASS(PatchDescriptor, "")
	CHILDCLASS(RANSACHomography, "")
	CHILDCLASS(Im, "")
	CHILDCLASS(Corner, "")
	CHILDCLASS(DescriptorSetClustering, "")
	{};

	CNumParam<int> VideoSource; //To enum
	CStringParam MosaicSaveDir;
	CNumParam<bool> MosaicDestination;
	CNumParam<int> TransformType, LMIterations, LM, PatchDescriptorType, PatchDescriptorSize, MaxFrames;
	CNumParam<bool> IncrementalRendering;
	CNumParam<int> FullFrameUpdateFreq;
	CNumParam<bool> SkipFrames;
	CNumParam<int> RendererType;
	CNumParam<int> FeatherRadius;
	CNumParam<int> DijkstraScale;
	CNumParam<int> WarpMethod;
	CNumParam<int> MosaicX;
	CNumParam<int> MosaicY;
	CNumParam<int> EvaluationFunction;
	//CNumParam<int> MaxRansacIters;
	CNumParam<int> MaxSearchForTransform, FeatureDetector, NumBlobScales;
	CNumParam<bool> MarkCorrespondences;
	//CNumParam<double> CorrespondenceQuality;

	MAKECHILDCLASS(BOW)
	MAKECHILDCLASS(BOWMatching)
	MAKECHILDCLASS(PatchDescriptor)
	MAKECHILDCLASS(RANSACHomography)
	MAKECHILDCLASS(Im)
	MAKECHILDCLASS(Corner)
	MAKECHILDCLASS(DescriptorSetClustering)
};

//! Load settings from config and use them to set up mosaicer components, and set them going.
int main(int argc, char* argv[])
{
    //grc::FeatureDetector * featureDetector = 0; //Delete these after (potentially) handling exceptions
    grc::Renderer * mosaicRenderer = 0;
    grc::EvaluationFunction * evaluationFunction = 0;

    try
    {
    	CMosaicingParams PARAMS(0, 0);
        config config(argc > 1 ? argv[1] : "default.cfg");

        PARAMS.init(&config);
        PARAMS.BOW.DescriptorBinning.RADIUS = PARAMS.PatchDescriptor.radius();

       	boost::scoped_ptr<CImageSource> pImageLoader( CImageSource::makeImageSource( PARAMS.Im));

        boost::scoped_ptr<CFeatureExtractor> pFeatureExtractor ( CFeatureExtractor::makeFeatureExtractor(PARAMS.Im, PARAMS.Corner, PARAMS.PatchDescriptor, PARAMS.DescriptorSetClustering));

        grc::Transform::eTransType_s = (grc::Enums::eTransformType)(int)PARAMS.TransformType;
        grc::Transform::warpMethod_s = (grc::Enums::eWarpMethod)(int)PARAMS.WarpMethod == grc::Enums::eWarpBilinear ? CV_INTER_LINEAR : CV_INTER_NN;
        grc::cDescriptor::eDescriptorType_s = (grc::Enums::eDescriptorInvarianceType)(int)PARAMS.PatchDescriptorType;
        grc::cDescriptor::eDescriptorSize_s = (grc::Enums::eDescriptorSize)(int)PARAMS.PatchDescriptorSize;

        grc::ImageSource imSource((grc::Enums::eVideoSource)(int)PARAMS.VideoSource, PARAMS.Im.ImageDir.IMAGE_DIR.asSz(), PARAMS.Im.VideoFile.FILENAME.asSz(), (grc::Enums::eDisplayMode)(int)PARAMS.MosaicDestination, PARAMS.MosaicSaveDir.asSz()); //Pure-Translation-Extended "H:/Projects/External/DTA003 Image mosaicing/Test Data/Rangiora"
        
        IplImage * pIm = cvLoadImage("/home/data/data/data/mosaicing/FixedCamera800/IMG_0004.JPG");

        grc::FeatureMatcher2 * featureMatcher = 0;

        //featureMatcher = new grc::BFFeatureMatcher(imSource, *pFeatureExtractor, PARAMS.BOWMatching);

        CBoW bow(PARAMS.BOW);
        featureMatcher = new grc::BoWFeatureMatcher(imSource, *pFeatureExtractor, bow, PARAMS.BOWMatching);

        grc::BaySACTransformEstimator transformEstimator(*featureMatcher, PARAMS.RANSACHomography, PARAMS.Im.getCamCalibrationMat(),
        		((bool)PARAMS.MarkCorrespondences) ? pImageLoader.get() : 0);

        switch(PARAMS.EvaluationFunction)
        {
        case grc::Enums::eSSDEvaluation:
            evaluationFunction = new grc::SumSquaredErrorEvaluationFunction;
            break;
        }
        if(evaluationFunction == 0 && PARAMS.EvaluationFunction != grc::Enums::eNoEvaluation)
            throw new grc::GRCException("main: No evaluation function initialised");

        CvSize mosaicSize = cvSize(PARAMS.MosaicX,PARAMS.MosaicY);
        grc::TransformEngine engine(transformEstimator, PARAMS.MaxFrames, PARAMS.IncrementalRendering ? 1 : 0, PARAMS.LM * PARAMS.LMIterations,
			mosaicSize, PARAMS.SkipFrames ? grc::Enums::eChooseSequentialSkip : grc::Enums::eChooseSequential, PARAMS.FullFrameUpdateFreq, PARAMS.MaxSearchForTransform);

        switch(PARAMS.RendererType)
        {
        case grc::Enums::eBasicRenderer:
		    mosaicRenderer = new grc::Renderer(imSource, engine, mosaicSize, evaluationFunction);
            break;
        case grc::Enums::eFeatheredRenderer:
            mosaicRenderer = new grc::FeatheredRenderer(imSource, engine, mosaicSize, PARAMS.FeatherRadius, evaluationFunction);
            break;
        case grc::Enums::eMultiScaleRenderer:
		    mosaicRenderer = new grc::MultiScaleFeatheredRenderer(imSource, engine, mosaicSize, PARAMS.FeatherRadius, evaluationFunction);
            break;
        case grc::Enums::eDijkstraRenderer:
		    mosaicRenderer = new grc::DijkstraCutRenderer(imSource, engine, mosaicSize, PARAMS.DijkstraScale, evaluationFunction);
            break;
        }
        if(mosaicRenderer == 0)
            throw new grc::GRCException("main: No renderer initialised");

		//This ensures images are captured by this thread:
        CStopWatch s;
        s.startTimer();

        imSource.doCaptureImages(&engine, mosaicRenderer);

        s.stopTimer();
        std::cout << s.getElapsedTime() << " seconds total" << endl;
        //PARAMS.printUseSummary();
        //PARAMS.printCfgFile();
    }
    catch(grc::GRCException * pEx)
    {
        std::cout << "ERROR: Unhandled exception: " << *(pEx->GetErrorMessage()) << std::endl << "Exiting...";
        cvDestroyAllWindows();
#ifndef __GNUC__
        Sleep(5000);
#endif
        delete pEx;
    }

    //delete featureDetector;
    delete mosaicRenderer;
    delete evaluationFunction;
}

pragma_warning(pop)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
