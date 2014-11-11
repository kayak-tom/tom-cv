#pragma warning (push)
#pragma warning (disable:4996)
#include "bow/bagOfWords.h"
#include <boost/filesystem.hpp>
#include <vector>
#include "cv.h" 
#include "highgui.h" 
#pragma warning (pop)

#include "params/param.h" 
#include "params/config.h" 

#include "featureExtract/featureExtractor.h"
#include "description/descriptor.h"
#include "description/vectorDescriptor.h"
#include "description/patchParams.h"
#include "featureExtract/cornerDetector.h"
#include "imageSource/imageSourceFromDir.h"

using namespace std;
using namespace boost::filesystem;


WRAPPERPARAMCLASS(BoWIO)
	PARAM(START_FRAME, 0, MAX_INT/2, 0) //TODO use some of these for IMU stuff
	PARAM(MAX_TRACK_LEN, 10, 10000, 10000)
	PARAM(MIN_TRACK_LEN_RECONSTRUCT, 0, 16, 2)
	PARAM(MAX_TRIES_FINDING_E, 1, 5, 1)
	PARAM(MIN_INLIERS_NEARBY, 8, 50, 20)
	PARAM(MIN_INLIERS_DISTANT, 8, 100, 28)
	PARAM(NEARBY_TIME, 0, 250, 15)
	PARAM(TARGET_CORRESPONDENCES, 0, 10000, 500)

    PARAMB(CORRECT_RD, true)

    CHILDCLASS(BOW)
	CHILDCLASS(BOWMatching)
	CHILDCLASS(PatchDescriptor)
	CHILDCLASS(Im)
	CHILDCLASS(Corner)
	CHILDCLASS(DescriptorSetClustering)
	{}

	CNumParam<int> START_FRAME, MAX_TRACK_LEN, MIN_TRACK_LEN_RECONSTRUCT, MAX_TRIES_FINDING_E, MIN_INLIERS_NEARBY, MIN_INLIERS_DISTANT, NEARBY_TIME, TARGET_CORRESPONDENCES;
	CNumParam<bool> CORRECT_RD;

	//BoWSLAM requires these parameterisation so instantiate inside
	MAKECHILDCLASS(BOW)
	MAKECHILDCLASS(BOWMatching)
	MAKECHILDCLASS(PatchDescriptor)
//MAKECHILDCLASS(RANSAC)
	//MAKECHILDCLASS(RANSACHomography)
	//MAKECHILDCLASS(LinkSelection)
	MAKECHILDCLASS(Im)
	MAKECHILDCLASS(Corner)
	//MAKECHILDCLASS(ResolveScale)
	MAKECHILDCLASS(DescriptorSetClustering)
};



bool isImage(const string & str)
{
    const char * szFileName = str.c_str();
    return strstr(szFileName, ".jpg") || strstr(szFileName, ".bmp") || strstr(szFileName, ".tif")
               || strstr(szFileName, ".JPG") || strstr(szFileName, ".BMP") || strstr(szFileName, ".TIF");
}

int main(int argc, char* argv[])
{
    //LoadSettingsFromCfg("BoWSLAM/defaultPatch.cfg");
	CBoWIOParams BOWPARAMS(0);
    
    {
        config cfg(argc == 1 ? "BoWSLAM/defaultPatch.cfg" : argv[1]);
	    BOWPARAMS.init(&cfg);
        
        BOWPARAMS.BOW.DescriptorBinning.RADIUS = BOWPARAMS.PatchDescriptor.radius();
    }

    CBoW bow(BOWPARAMS.BOW);

    string strPath, strPath2;

#define IO_DIR "IO"
    create_directory( IO_DIR );

    boost::scoped_ptr<CFeatureExtractor> pFeatureExtractor ;

    remove("matches");
    cout << "Polling " << IO_DIR << " for images or commands...\n";

    for(;;)
    {
        directory_iterator end_itr;
        for ( directory_iterator itr( IO_DIR ); itr != end_itr; ++itr )
        {
            if ( is_regular(itr->status()) )
            {
        	    const path p = itr->path();
                string strImFilename = p.relative_path().file_string();
                if(isImage(strImFilename)) 
                {
                    cout << "Adding image " << strImFilename << endl;
                    if(!pFeatureExtractor)
                    {
                        CImageSourceFromDir(IO_DIR, BOWPARAMS.Im) ;
                        pFeatureExtractor.reset( CFeatureExtractor::makeFeatureExtractor(BOWPARAMS.Im, BOWPARAMS.Corner, BOWPARAMS.PatchDescriptor, BOWPARAMS.DescriptorSetClustering));
                    }

                    char name[1000];
                    sprintf_s(name, 1000, "%s", strImFilename.c_str());
                    char * pcDot = strrchr(name, '.');
                    if(!pcDot)
                    {
                        cout << "Error parsing filename\n";
                        return -1;
                    }
                    *pcDot = 0;
                    char * pcNameId = strrchr(name, '\\');
                    if (pcNameId == 0) 
                        pcNameId = name; //no slash
                    else
                        pcNameId++;

                    while(*pcNameId < '0' || *pcNameId > '9')
                        pcNameId++;

                    if(!pcNameId)
                        cout << "Error: Could not find id in filename\n";
                    else
                    {

                        int nId = atoi(pcNameId);

                        IplImage * pIm = cvLoadImage(strImFilename.c_str());
                        CDescriptorSet * pDesc = 0;
                        
                        if(pIm)
                        {
                            //cvCopy(pIm, rgbImage);
                            pDesc = pFeatureExtractor->getDescriptors(pIm);
                        }

                        if(pDesc && pDesc->Count())
                        {
                            bow.addImage(&pDesc, nId); //...and add it to the Bag-of-Words database with id 0
                            cout << "Added image with id " << nId << endl;
                        }
                        else
                            cout << "Failed to open/find descriptors from this image\n";
                    }
                }
                else if (strstr(strImFilename.c_str(), "recluster"))
                {
                    bow.recreateDictionary(); 
                }
                else if (strstr(strImFilename.c_str(), "quit"))
                {
                    remove(strImFilename.c_str());
                    return 0;
                } 
                else if (strstr(strImFilename.c_str(), "bb"))
                {
                    char ids[1000];
                    const char * pcId1=strstr(strImFilename.c_str(), "bb")+2;
                    sprintf_s(ids, 1000, "%s", pcId1);

                    char * pcId2=strchr(ids, ',');
                    if(!pcId2)
                    {
                        cout << "Error getting correspondences (no comma)" << endl;
                    }
                    else
                    {
                        *pcId2 = 0; pcId2++;
                        if(!*pcId2)
                        {
                            cout << "Error getting correspondences (no second id)" << endl;
                        }
                        else
                        {
                            int id1 = atoi(ids);
                            int id2 = atoi(pcId2);
                            if(id1==id2)
                            {
                                cout << "Ids are the same" << endl;
                            }
                            else
                            {
                                boost::scoped_ptr<CBoWCorrespondences> pCorr(bow.getCorrespondences(id1, id2, BOWPARAMS.BOWMatching));

                                while(exists("correspondences"))
                                {
                                    cout << "Waiting for previous 'correspondences' file to be removed..." << endl;
                                    Sleep(500);
                                }
                                ofstream outFile("correspondences");

                                for(CBoWCorrespondences::const_iterator pC = pCorr->begin(); pC != pCorr->end(); pC++)
                                    outFile << '(' << pC->Location1().dx() << ',' << pC->Location1().dy() << "),(" << pC->Location2().dx() << ',' << pC->Location2().dy() << ")," << pC->priorProb() << endl;
                                
                                outFile.close();

                            }
                        }
                    }
                    remove(strImFilename.c_str());
                } 
                else
                {
                    char name[1000];
                    sprintf_s(name, 1000, "%s", strImFilename.c_str());
                    char * pcDot = strrchr(name, '.');
                    if(pcDot)
                        *pcDot = 0;

                    char * pcNameId = strrchr(name, '\\');
                    if (pcNameId == 0) 
                        pcNameId = name; //no slash
                    else
                        pcNameId++;

                    int nId = atoi(pcNameId);

                    while(exists("matches"))
                    {
                        cout << "Waiting for previous results to be removed..." << endl;
                        Sleep(500);
                    }

                    TBoWMatchVector * pMatches = bow.getMatches(nId);

                    ofstream outFile("matches");

                    for(TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end(); pMatch++)
                        outFile << (pMatch)->id() << ',' << (pMatch)->MatchStrength() << endl;
                    
                    outFile.close();

                    delete pMatches;
                }

                remove(strImFilename.c_str());
            }
        }
        Sleep(10);
    }
}