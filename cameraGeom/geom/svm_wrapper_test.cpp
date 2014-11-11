#include "svm_wrapper.h"
#include <Eigen/Core>
#include <util/random.h>
#include <boost/scoped_ptr.hpp>
#include <boost/filesystem.hpp>

static const char * szDir = "SVMTraining/SVMAutotest";
static const char * szLabel = "SVMAutotest";
/**
 * @class CTestFeature
 * @brief Only even elements are significant
 * Each component is scaled by its idx.
 */
class CTestFeature : public CSVMFeature_base
{
    static const int DIMS = 10;
    static bool s_bSetup;
    Eigen::VectorXd adFeature;
    
public:
    CTestFeature(const Eigen::VectorXd & feature) : CSVMFeature_base((int)feature.rows()), adFeature(feature)
    {
        
        //cout << bPos << ": " << adFeature << endl;
    }
    
    virtual double value_int(const int nIdx)
    {
        return adFeature(nIdx);
    }
};

class CTestFeatureFactory
{
    Eigen::VectorXd aadFeatureMean[2]; //Negative and positive
    
    void setupFeatureMean(const int DIMS)
    {
        aadFeatureMean[false] = Eigen::VectorXd::Random(DIMS);
        //aadFeatureMean[false] *= 0;
        for(int nFeature=0; nFeature < DIMS; nFeature++)
            aadFeatureMean[false](nFeature) *= (1+nFeature);
            
        aadFeatureMean[true] = aadFeatureMean[false];
        for(int nFeature=0; nFeature < DIMS; nFeature+=2)
            aadFeatureMean[true](nFeature)+=2*(1+nFeature);
    }    
public:
    CTestFeatureFactory(const int nDims)
    {
        setupFeatureMean(nDims);
    }
    
    CTestFeature * makeFeature(const bool bPos) const
    {
        Eigen::VectorXd adFeature = aadFeatureMean[bPos];
        for(int nFeature=0; nFeature < adFeature.rows(); nFeature++)
        {
            adFeature(nFeature) += CRandom::Normal(0, 1+nFeature);
        }
        return new CTestFeature(adFeature);
    }
};

typedef std::vector<CTestFeature *> TFeatureVec;

bool label(const int i)
{
    return i%4==0;
}

TFeatureVec trainSVMWrapper(const CTestFeatureFactory & testFeatureFactory, const int nNumTrainingFeatures)
{
    boost::scoped_ptr<CSVMTraining_base> pTrainer(CSVMTraining_base::makeSVMTraining(szDir, szLabel, 1.0f, CSVMTraining_base::eFFS, true));
    
    TFeatureVec aFeatures;
    
    for(int i=0; i<nNumTrainingFeatures; i++)
    {
        CTestFeature * pFeature = testFeatureFactory.makeFeature(label(i));
        pTrainer->addTrainingFeature(pFeature, label(i));
        aFeatures.push_back(pFeature);
    }
    
    //Destructor triggers training
    
    return aFeatures;
}

void testSVMWrapper()
{
    boost::filesystem::remove_all("SVMAutotest"); //clean up old autotests...
    
    const CTestFeatureFactory testFeatureFactory(10);
    
    TFeatureVec aFeatures = trainSVMWrapper(testFeatureFactory, 400);
    
    boost::scoped_ptr<CSVMClassifier_base> pClassifier(CSVMClassifier_base::makeSVMClassifier(szDir, szLabel, CSVMClassifier_base::NO_PRECISION));
    
    double dProbTrainingLabelledCorrectly = 0, dProbTestLabelledCorrectly = 0;
    
    for(int i=0; i<(int)aFeatures.size(); i++)
    {
        CTestFeature * pOldFeature = aFeatures[i];
        
        cout << "Testing " << i << " label=" << label(i) << endl;
        
        const bool bPredictedLabel = pClassifier->binaryClassify(pOldFeature);

        delete pOldFeature;
        
        if(bPredictedLabel == label(i))
            dProbTrainingLabelledCorrectly++;
            
        boost::scoped_ptr<CTestFeature> pNewFeature (testFeatureFactory.makeFeature(label(i)));
        const bool bNewPredictedLabel = pClassifier->binaryClassify(pNewFeature.get());
        if(bNewPredictedLabel == label(i))
            dProbTestLabelledCorrectly++;
            
    }
    
    dProbTrainingLabelledCorrectly /= (double)aFeatures.size();
    dProbTestLabelledCorrectly /= (double)aFeatures.size();
    
    cout << "Proportion of training vectors labelled correctly: " << dProbTrainingLabelledCorrectly << endl;
    cout << "Proportion of test vectors labelled correctly: " << dProbTestLabelledCorrectly << endl;

    CHECK_P(dProbTrainingLabelledCorrectly < 0.9, dProbTrainingLabelledCorrectly, "Looks like SVM training has failed");
    CHECK_P(dProbTestLabelledCorrectly < 0.9, dProbTestLabelledCorrectly, "Looks like SVM training has failed");
}



