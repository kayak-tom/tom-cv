#ifndef CSVMWRAPPER_H
#define CSVMWRAPPER_H

#include <opencv2/core/core.hpp>

/*
 * For each classification problem
 * * There is always a CSVMClassifier_base, which loads settings (SVM state, feature subset, normalising coefficients, sign ambiguity), requests feature components, normalises feature components, classifies.
 * * The prediction is usually done with the CSVMClassifier_base, but this might be turned off the first pass so we can use decisions from the existing classifier.
 * * For each classification we create a CSVMFeature_base, and call result = CSVMClassifier_base.classify(feature)
 * When in training mode (CSVMTraining_base* pTrainer exists), features are generally pulled from training images automatically, so we don't need to save them. ~CSVMTraining_base does the training.
 * * [no If there is no training data then the prediction is based on the training label (hence we can use labelled training data for both) ]
 */
 
class CSVMFeature_base
{
protected:
    cv::Mat feature;
public:
    CSVMFeature_base(const int nDim);
    virtual ~CSVMFeature_base() {}
    
    // Total number of features (from which we will choose a subset)
    int dimension() const;
    
    //Compute, cache, and return a feature value (may compute + cache other feature values as a side effect)
    double value(const int nIdx);

    
    //Compute every coefficient and return feature (used for training)
    const cv::Mat & getEntireFeature(); 
    static float UNINIT();

protected:
    //Compute and return a feature value
    virtual double value_int(const int nIdx) = 0;
    
};

class CSVMClassifier_base
{
public:
    CSVMClassifier_base() {}
    virtual ~CSVMClassifier_base() {}

    static CSVMClassifier_base * makeSVMClassifier(const std::string path, const std::string label, const double dPrecision);
    
    virtual double classify(CSVMFeature_base * pFeature) = 0;
    
    bool binaryClassify(CSVMFeature_base * pFeature) { return classify(pFeature) > 0; }

    static const double NO_PRECISION;
};

class CSVMTraining_base
{
    
public:
    enum eSVMFeatureSelectionMethod {eFFS, eBFS, eLoadFromFile, eNoFS }; 

    CSVMTraining_base() {}
    virtual ~CSVMTraining_base() {}

    static CSVMTraining_base * makeSVMTraining(const std::string path, const std::string label, const float fNegRelativeWeight, const eSVMFeatureSelectionMethod featureSelectionMode, const bool bFilterHyperparams);
    
    virtual void addTrainingFeature(CSVMFeature_base * pFeature, const bool bLabel) = 0;
};

void testSVMWrapper();

#endif // CSVMWRAPPER_H
