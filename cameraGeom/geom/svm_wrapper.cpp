/*
 * Implements:
 * * Train boosting tests to classify the easy cases
 * * Normalise features
 * * Select [feature subset, hyperparameters] by K-fold cross-validation
 * *
 *
 *
 * */
#include <util/convert.h>
#include "svm_wrapper.h"
#include <map>
#include <set>
#include <opencv2/ml/ml.hpp>
#include <boost/foreach.hpp>
#include <util/vectorUtil.h>
#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/bind.hpp>
#include "threadpool.h"
#include <fstream>
#include <Eigen/Core>
#include <opencv2/core/core_c.h>
#include <boost/thread/mutex.hpp>

static const int SVM_TYPE = cv::SVM::NU_SVC;

typedef std::vector<cv::Mat> TLabelledFeatures;

typedef std::pair<std::vector<double>, std::vector<double> > TPRLookup;

cv::Mat vectorToMat(const TLabelledFeatures & aLabelledFeatures)
{
    CHECK(aLabelledFeatures.size()==0, "No vector data");
    cv::Mat M((int)aLabelledFeatures.size(), aLabelledFeatures[0].cols, aLabelledFeatures[0].type(), cv::Scalar());
    for(int i=0; i<M.rows; i++) {
        aLabelledFeatures[i].copyTo(M.row(i));

        /*cout << "M=" << M << endl;
        cout << "M.row(i)=" << M.row(i) << endl;
        cout << "aLabelledFeatures[i]=" << aLabelledFeatures[i] << endl;*/

    }
    CHECK(M.size().area()==0, "No mat data");

    return M;
}


class CBoosterState
{
    int nFeatureIdx;
    double dThreshold;
    bool bRejectAbove;
public:
    CBoosterState() : nFeatureIdx(-1), dThreshold(HUGE), bRejectAbove(false) {}

    CBoosterState(const int nFeatureIdx, const double dThreshold, const bool bRejectAbove) : nFeatureIdx(nFeatureIdx), dThreshold(dThreshold), bRejectAbove(bRejectAbove) {}

    bool rejectAbove() const {
        return bRejectAbove;
    }
    double getThreshold() const {
        return dThreshold;
    }
    int getFeatureIdx() const {
        return nFeatureIdx;
    }
};

std::ostream & operator<<(std::ostream & s, const CBoosterState & info)
{
    s << "CBoosterState: rejectAbove=" << info.rejectAbove() << "getThreshold=" << info.getThreshold() << "getFeatureIdx=" << info.getFeatureIdx() << endl;
    return s;
}

typedef std::vector<CBoosterState> TBoosterStates; //An ordered cascade of boosters

typedef std::vector<int> TFeatureIdxSubset;
typedef std::vector<double> TNormalisingCoefficients; //multiply features by these to get normalised features


class CFeatureSubsetSelecter
{
    TFeatureIdxSubset anFeatureIdxSubset;
    TNormalisingCoefficients adNormalisingMean, adNormalisingScale;
public:
    CFeatureSubsetSelecter(const TFeatureIdxSubset & anFeatureIdxSubset, const TNormalisingCoefficients & adNormalisingMean, const TNormalisingCoefficients & adNormalisingScale) : anFeatureIdxSubset(anFeatureIdxSubset), adNormalisingMean(adNormalisingMean), adNormalisingScale(adNormalisingScale) {}
    CFeatureSubsetSelecter() {}

    void setFeatureIdxSubset(const TFeatureIdxSubset& anFeatureIdxSubset_in) {
        anFeatureIdxSubset = anFeatureIdxSubset_in;
    }

    const TFeatureIdxSubset & getFeatureIdxSubset() const {
        return anFeatureIdxSubset;
    }

    void load(const cv::FileStorage & fs) {
        const bool bVerbose = false;

        cv::Mat substateIndexMat;
        fs["featureSubset"] >> substateIndexMat;
        for(int i=0; i<substateIndexMat.rows; i++)
            anFeatureIdxSubset.push_back(substateIndexMat.at<int>(i));
        fs["normalisingMean"] >> adNormalisingMean;
        fs["normalisingScale"] >> adNormalisingScale;

        if(bVerbose)
        {
            cout << "Loaded feature index subset " << substateIndexMat << endl;
            cout << "Loaded normalisingMean " << cv::Mat(adNormalisingMean) << endl;
            cout << "Loaded normalisingScale " << cv::Mat(adNormalisingScale) << endl;
        }
    }

    void selectAndNormalise(CSVMFeature_base * pFeature, cv::Mat & featureSubset) const {
        const bool bVerbose = false;

        featureSubset.create(cv::Size((int)anFeatureIdxSubset.size(), 1), CV_32FC1);
        int i=0;
        BOOST_FOREACH(int nIdx, anFeatureIdxSubset) {
            featureSubset.at<float>(i) = (float)(adNormalisingScale[nIdx]*(pFeature->value(nIdx) - adNormalisingMean[nIdx]));
            i++;
        }
        if(bVerbose) cout << "Selected "  << featureSubset << " of " << pFeature->getEntireFeature() << endl;
    }

    void save(cv::FileStorage & fs) const {
        cv::Mat substateIndexMat(anFeatureIdxSubset);
        fs << "featureSubset" << substateIndexMat;
        fs << "normalisingMean" << adNormalisingMean;
        fs << "normalisingScale" << adNormalisingScale;
    }

    cv::Mat selectAndNormalise(const cv::Mat & feature) const {
        const bool bVerbose = false;

        CHECK(anFeatureIdxSubset.size() == 0, "Empty feature subset");

        cv::Mat featureSubset(cv::Size((int)anFeatureIdxSubset.size(), 1), feature.type(), cv::Scalar());
        int i=0;
        BOOST_FOREACH(int nIdx, anFeatureIdxSubset) {
            featureSubset.at<float>(i) = (float)(adNormalisingScale[nIdx]*(feature.at<float>(nIdx) - adNormalisingMean[nIdx]));
            i++;
        }
        if(bVerbose) {
            cout << "Selected feature subset " << featureSubset << endl;
            cout << "of feature " << feature << endl;
        }

        CHECK(featureSubset.size().area() == 0, "No data selected");

        return featureSubset;
    }

    void findNormalisingCoeffs(const TLabelledFeatures * aaadFeatures, const int nFeature, double & dMean, double & dScale) const { // http://en.wikipedia.org/wiki/Feature_scaling
        std::vector<double> adVals;

        for(int bLabel=0; bLabel < 2; bLabel++) {
            BOOST_FOREACH(const cv::Mat & feature, aaadFeatures[bLabel]) {
                adVals.push_back(feature.at<float>(nFeature));
            }
        }

        dMean=0;
        double dSD=0;
        meanSD(adVals, dMean, dSD);
        dScale = (dSD > 0) ? 1.0/dSD : 1;
    }

    void findNormalisingCoeffs(const TLabelledFeatures * aaadFeatures) {
        const bool bVerbose = true;

        const int nDims = aaadFeatures[false][0].cols;
        adNormalisingMean = TNormalisingCoefficients(nDims, 0);
        adNormalisingScale = TNormalisingCoefficients(nDims, 0);

        for(int nFeature=0; nFeature < nDims; nFeature++) {
            findNormalisingCoeffs(aaadFeatures, nFeature, adNormalisingMean[nFeature], adNormalisingScale[nFeature]);

            CHECK(adNormalisingScale[nFeature] == 0, "Bad scale");
            CHECKBADNUM(adNormalisingScale[nFeature]);
        }
        if(bVerbose) {
            cout << "normalisingMean=" << cv::Mat(adNormalisingMean) << endl;
            cout << "normalisingScale=" << cv::Mat(adNormalisingScale) << endl;
        }
    }
};

class CSavedSVMState
{
    const std::string path, label;
    std::string fullFilename, svmFilename;
    TBoosterStates boosterStates;
    CFeatureSubsetSelecter featureSubset;
    //TFeatureIdxSubset anFeatureSubset;
    //TNormalisingCoefficients adNormalisingMean, adNormalisingScale;
    double dSignCorrection;
    double dClassificationBoundary;
    
    static double interpPrecisionBoundary(const double dTargetPrecision, const TPRLookup & PRLookup) 
    {
        //Interpolate off 2 closest (should work off the ends as well)
        double dBoundary1=0,dBoundary2=0,dP1=HUGE,dP2=HUGE;
        {
            for(int i=0;i<(int)PRLookup.first.size();i++)
            {
                const double dBoundary = PRLookup.first[i];
                const double dPrecision = PRLookup.second[i];
                if(fabs(dPrecision - dTargetPrecision) < fabs(dP1 - dTargetPrecision))
                {
                    dP2 = dP1; dBoundary2 = dBoundary1;
                    dP1 = dPrecision; dBoundary1 = dBoundary;
                }
                else if(fabs(dPrecision - dTargetPrecision) < fabs(dP2 - dTargetPrecision))
                {
                    dP2 = dPrecision; dBoundary2 = dBoundary;
                }
            }
            
            // Fit straight line:
            // dBoundary1 = m*dP1+c
            // dBoundary2 = m*dP2+c
            
            const double m = (dBoundary2-dBoundary1)/(dP2-dP1);
            const double c = dBoundary1-m*dP1;
            
            CHECK(!zero(c-(dBoundary2-m*dP2)), "Line fit failed (probably a horizontal line)");
            
            return dTargetPrecision*m+c;
        }
    }
public:
    //Load settings:
    CSavedSVMState(const std::string path, const std::string label, const double dPrecision) : path(path), label(label), dClassificationBoundary(0) {
        CHECK_P(!boost::filesystem::exists(getExtraFilename()), getExtraFilename(), "Saved SVM state file doesn't exist");
        cv::FileStorage fs(getExtraFilename().c_str(), cv::FileStorage::READ);
        cv::Mat boosterStatesMat;
        fs["boosterStates"] >> boosterStatesMat;
        for(int i=0; i<boosterStatesMat.rows; i++) {
            const cv::Mat row = boosterStatesMat.row(i);
            boosterStates.push_back(CBoosterState((int)row.at<double>(0), row.at<double>(1), (bool)row.at<double>(2)));
        }

        featureSubset.load(fs);

        fs["signCorrection"] >> dSignCorrection;
        
        TPRLookup PRLookup;
        if(!fs["boundaries"].empty() && dPrecision != CSVMClassifier_base::NO_PRECISION)
        {
            fs["boundaries"] >> PRLookup.first;
            fs["precision"] >> PRLookup.second;
            dClassificationBoundary = interpPrecisionBoundary(dPrecision, PRLookup);
        }
        
        if(dSignCorrection == 0)
            cout << "Boosted classifier has no svm training--all training outliers can be removed by boosting." << endl;
    }

    //Save back settings from training
    CSavedSVMState(const std::string path, const std::string label, const TBoosterStates & boosterStates, const CFeatureSubsetSelecter & featureSubset, cv::SVM & svm, double dSignCorrection, TPRLookup & PRLookup)
        : path(path), label(label), boosterStates(boosterStates), featureSubset(featureSubset), dSignCorrection(dSignCorrection) {

        boost::filesystem::create_directories(path);
        cv::FileStorage fs(getExtraFilename().c_str(), cv::FileStorage::WRITE);

        CHECK_P(!boost::filesystem::exists(getExtraFilename()), getExtraFilename(), "yaml file doesn't exist");

        cv::Mat boosterStatesMat(cv::Size(3, (int)boosterStates.size()), CV_64FC1, cv::Scalar(-1));
        for(int i=0; i<boosterStatesMat.rows; i++) {
            cv::Mat row(cv::Size(3,1), CV_64FC1, cv::Scalar(-1));
            row.at<double>(0) = boosterStates[i].getFeatureIdx();
            row.at<double>(1) = boosterStates[i].getThreshold();
            row.at<double>(2) = boosterStates[i].rejectAbove();
            row.copyTo(boosterStatesMat.row(i));
        }
        fs << "boosterStates" << boosterStatesMat;

        fs << "signCorrection" << dSignCorrection; //0 == boosting only.
        featureSubset.save(fs);
        
        fs << "boundaries" << PRLookup.first;
        fs << "precision" << PRLookup.second;
        
        if(dSignCorrection != 0) {
            svm.save(getSVMFilename().c_str(), label.c_str());
        }
    }

    const CFeatureSubsetSelecter & getFeatureSubset() const {
        return featureSubset;
    }

    const TBoosterStates getBoosterStates() const {
        return boosterStates;
    }

    const std::string getExtraFilename() {
        fullFilename = path + "/savedSVMstate" + label + "_subset.yaml";
        return fullFilename;
    }

    const std::string getSVMFilename() {
        svmFilename = path + "/savedSVMstate" + label + ".yaml";
        return svmFilename;
    }

    double getSignCorrection() const {
        return dSignCorrection;
    }
    
    double getClassificationBoundary() const
    {
        return dClassificationBoundary;
    }

};

///////////////// Feature ////////////////////////////

CSVMFeature_base::CSVMFeature_base(const int nDim) : feature(cv::Size(nDim, 1), CV_32FC1, cv::Scalar(UNINIT()) )
{

}

double CSVMFeature_base::value(const int nIdx)
{
    float & existingVal = feature.at<float>(nIdx);

    if(existingVal == UNINIT())
        existingVal = (float)value_int(nIdx);
    CHECKBADNUM(existingVal);
    return existingVal;
}

// Total number of features (from which we will choose a subset)
int CSVMFeature_base::dimension() const
{
    return feature.cols;
}

float CSVMFeature_base::UNINIT()
{
    return -1e+20f;
}

//Compute every coefficient and return feature (used for training)
const cv::Mat & CSVMFeature_base::getEntireFeature()
{
    for(int i=0; i<dimension(); i++) {
        value(i);
    }

    return feature;
}


/////////////// Boosted filter ///////////////////////

// Only used to remove negative examples. keepPotentialCandidate returns false for candidates which are almost certainly negative.

class CBoostedFilter
{
    CBoosterState state;

    bool keepPotentialCandidate(const double dFeatureVal) const {
        if(state.rejectAbove())
            return (dFeatureVal < state.getThreshold());
        else
            return (dFeatureVal > state.getThreshold());
    }
public:
    CBoostedFilter() {}
    CBoostedFilter(const CBoosterState & state) : state(state) {

    }

    bool keepPotentialCandidate(CSVMFeature_base * pFeature) const {
        const double dFeatureVal = pFeature->value(state.getFeatureIdx());
        return keepPotentialCandidate(dFeatureVal);
    }

    bool keepPotentialCandidate(const cv::Mat & feature) const {
        const double dFeatureVal = feature.at<float>(state.getFeatureIdx());
        CHECKBADNUM(dFeatureVal);
        return keepPotentialCandidate(dFeatureVal);
    }
};

class CBoostedFilters
{
    typedef std::vector<CBoostedFilter> TBoostedFilters;
    TBoostedFilters aFilters;
public:
    CBoostedFilters(const TBoosterStates & aBoosterStates) {
        BOOST_FOREACH(const CBoosterState & boosterState, aBoosterStates) {
            aFilters.push_back(CBoostedFilter(boosterState));
        }
    }

    bool keepPotentialCandidate(CSVMFeature_base * pFeature) const {
        BOOST_FOREACH(const CBoostedFilter & filter, aFilters) {
            if(!filter.keepPotentialCandidate(pFeature))
                return false;
        }

        return true;
    }

};


//////////////// Classifier //////////////////////////

class CSVMClassifier : public CSVMClassifier_base
{
    CSavedSVMState savedState;
    CBoostedFilters boostedFilters;
    CFeatureSubsetSelecter featureSubsetSelecter;
    cv::Mat feature;
    cv::SVM svm;
public:
    CSVMClassifier(const std::string path, const std::string label, const double dPrecision)
        : savedState(path, label, dPrecision),
          boostedFilters(savedState.getBoosterStates()),
          featureSubsetSelecter(savedState.getFeatureSubset()) {
        if(savedState.getSignCorrection() != 0)
            svm.load(savedState.getSVMFilename().c_str(), label.c_str());
    }

    virtual ~CSVMClassifier() {}

    //Return SVM return value
    virtual double classify(CSVMFeature_base * pFeature) {
        const bool bVerbose = false;

        if(!boostedFilters.keepPotentialCandidate(pFeature))
            return -1; //-ve

        if(savedState.getSignCorrection() == 0)
            return 1; //All remaining points are inliers

        featureSubsetSelecter.selectAndNormalise(pFeature, feature);
        const double dSVMVal = svm.predict(feature, true) - savedState.getClassificationBoundary();

        if(bVerbose) cout << "Raw response: " << dSVMVal << endl;

        return savedState.getSignCorrection()*dSVMVal;
    }

};

CSVMClassifier_base * CSVMClassifier_base::makeSVMClassifier(const std::string path, const std::string label, const double dPrecision)
{
    return new CSVMClassifier(path, label, dPrecision);
}

////////////////// Training ///////////////////////////

class CSVMTraining : public CSVMTraining_base
{
    const std::string path, label;

    static const bool bMT = true;
    static const int K = 6; // for k-fold cross-validation (= # threads)
    boost::scoped_ptr<CThreadpool_base> pSVMThreadpool;

    TLabelledFeatures aaadFeatures[2]; //[0] negative and [1] positive examples

    float //fNegativeScore;+ve examples have weight 1, negative examples have weight fNegativeScore so that the appropriate cost is being minimised
    fNegRelativeWeight; //Cost of incorrectly labelling negative examples, relative to +ve weight.

    CvMat * pClassWeights;

    eSVMFeatureSelectionMethod featureSelectionMode;
    const bool bFilterHyperparams;

    boost::mutex mxLockToAdd;

    std::map<double, CBoosterState, std::greater<double> > aBoosterStates;
    typedef std::pair<double, CBoosterState> TBoosterCandidate;

    typedef std::pair<double, bool> TOneFeatureVal;

    static bool sortByOneFeature_asc(const TOneFeatureVal & d1, const TOneFeatureVal & d2) {
        return d1.first < d2.first;
    };
    static bool sortByOneFeature_desc(const TOneFeatureVal & d1, const TOneFeatureVal & d2) {
        return d1.first > d2.first;
    };
    /*
     * A worthwhile boosted classifier will remove (e.g.) at least 10% of -ve's and at most 1 +ve for every 100 -ve's removed (there's usually a lot more - than +).
     * */
    TBoosterCandidate findBoosterState(const int nIdx, const bool bRejectAbove) {
        const bool bVerbose = true;

        std::vector<TOneFeatureVal> aSortedFeatures; //sorted in the order they are to be applied

        for(int bLabel=0; bLabel < 2; bLabel++) {
            BOOST_FOREACH(const cv::Mat & feature, aaadFeatures[bLabel]) {
                aSortedFeatures.push_back(std::pair<double,bool>(feature.at<float>(nIdx), (bool) bLabel));
            }
        }

        std::sort(aSortedFeatures.begin(), aSortedFeatures.end(), bRejectAbove ? sortByOneFeature_desc : sortByOneFeature_asc);
        double dPos=0, dNeg=0;
        double dBestThreshPosBelow = -HUGE;
        double dNumRemovedBelow = -1;
        cout << "Sorted features: " << (aSortedFeatures.size()) << endl;
        for(int i=0; i<((int)aSortedFeatures.size())-1; i++) {
            if(bVerbose && i<10) cout << i << ": " << aSortedFeatures[i].first << "-" << aSortedFeatures[i].second << endl;

            if(aSortedFeatures[i].second)
                dPos++;
            else
                dNeg++;

            //Is a split below here good enough to be worthwhile?
            if(dPos < 0.01*dNeg && aSortedFeatures[i].first != aSortedFeatures[i+1].first) { //Force the split to occur away from integral values
                dBestThreshPosBelow = 0.5*(aSortedFeatures[i].first + aSortedFeatures[i+1].first);
                dNumRemovedBelow = dNeg;
            }
        }
        const double dMinimumPower = 0.1;
        const double dPropOfNegativeExamplesRemoved = dNumRemovedBelow / (double)aaadFeatures[false].size();
        if(dPropOfNegativeExamplesRemoved < dMinimumPower || dNumRemovedBelow < 50) {
            return TBoosterCandidate(0, CBoosterState());
        }
        
        if(bVerbose)cout << "One candidate " << dPropOfNegativeExamplesRemoved << " nIdx=" << nIdx << " dBestThreshPosBelow=" << dBestThreshPosBelow << " bRejectAbove=" << bRejectAbove << endl;
        
        return TBoosterCandidate(dPropOfNegativeExamplesRemoved, CBoosterState(nIdx, dBestThreshPosBelow, bRejectAbove));
    }

    static bool sortBoosterStates(const TBoosterCandidate & d1, const TBoosterCandidate & d2) {
        return d1.first < d2.first;
    }

    TBoosterStates findBoosterStates() {
        TBoosterStates aBoosterStates;

        for(;;) {
            CBoosterState boosterState = findBestBoosterState();
            if(boosterState.getFeatureIdx()<0)
                return aBoosterStates;


            //Otherwise we found a booster state.
            aBoosterStates.push_back(boosterState);
            CBoostedFilter filter(boosterState);
            
            const int nNegFeaturesBefore = (int)aaadFeatures[0].size();
            //Filter and repeat.
            for(int bLabel=0; bLabel < 2; bLabel++) {
                TLabelledFeatures newFeatures;
                BOOST_FOREACH(const cv::Mat & feature, aaadFeatures[bLabel]) {
                    if(filter.keepPotentialCandidate(feature)) {
                        newFeatures.push_back(feature);
                    }
                }
                aaadFeatures[bLabel] = newFeatures;
            }
            
            const int nNegFeaturesAfter = (int)aaadFeatures[0].size();
            
            CHECK_P(nNegFeaturesAfter >= nNegFeaturesBefore, nNegFeaturesAfter, "Boosting failed");
        }
    }

    CBoosterState findBestBoosterState() {
        const bool bVerbose = true;

        //Try each feature in turn. Find high/low percentiles of +ve examples and their percentile amongst negative examples
        TBoosterCandidate bestBoosterCandidate(0, CBoosterState());
        for(int nRejectAbove=0; nRejectAbove<2; nRejectAbove++) {
            for(int nIdx=0; nIdx<aaadFeatures[false][0].cols; nIdx++) {
                const TBoosterCandidate boosterCandidate = findBoosterState(nIdx, (bool)nRejectAbove);
                if(boosterCandidate.first > bestBoosterCandidate.first) {
                    bestBoosterCandidate = boosterCandidate;
                }
            }
        }

        if(bVerbose) {
            cout << "findBestBoosterState, proportion removed: " << bestBoosterCandidate.first << endl;
            cout << bestBoosterCandidate.second << endl;
        }
        return bestBoosterCandidate.second;
    }

    /*void normaliseFeatures() {
        const bool bVerbose = true;
        for(int bLabel=0; bLabel < 2; bLabel++) {
            BOOST_FOREACH(cv::Mat & feature, aaadFeatures[bLabel]) {
                if(bVerbose)
                    cout << "Feature before normalisation " << feature << endl;

                feature -= normalisingMean;
                cv::multiply(feature, normalisingScale, feature);

                if(bVerbose)
                    cout << "Feature after normalisation " << feature << endl;
            }
        }
    }*/
    double propCorrect(const cv::SVM & svm, const cv::Mat & features, const cv::Mat & labels, double & dSignCorrection, double & dPrecision, const double dBoundary) const {
        const bool bVerbose = (dBoundary!=0);

        cv::Mat test_labels = 0*labels;
        //CvMat testLabelsMat = test_labels, featuresMat = features;
        if(bVerbose) cout << "Validating with " << labels.rows << " features" << endl;

        try {

            /*if(dBoundary==0)
            {
                svm.predict(features, test_labels);
            }
            else
            {*/ //This *Sometimes* flips the sign relative to svm.predict(features, test_labels);
            for(int i=0; i<test_labels.rows; i++) {
                cv::Mat row = features.row(i);
                const double dSVMResponse = svm.predict(row, true);
				test_labels.at<float>(i) = (float)(dSVMResponse - dBoundary);

                //cout << "Feature " << row << " response " << dSVMResponse << endl;
            }
            //}
        } catch(const cv::Exception & ex) {
            cout << "Error on svm predict: " << ex.what() << endl;
            cout << "Nu=" << svm.get_params().nu << endl;
            cout << "Gamma=" << svm.get_params().gamma << endl;
            cout << "SVs=" << svm.get_support_vector_count() << endl;

            return 0;
        }

        const bool bClassificationError=false;
        if(bClassificationError) {
            CHECK(dBoundary!=0, "Now using full DF value");
            const int pos_ex = cv::countNonZero(labels == 1);
            const int neg_ex = cv::countNonZero(labels == -1);
            const double dTrivialPropCorrect = std::max<double>(pos_ex, neg_ex)/labels.rows;

            const cv::Mat errors = test_labels - labels;

            const int nNumWrong = cv::countNonZero(errors);
            double dPropCorrect = (1.0-((double)nNumWrong/labels.rows));

            dSignCorrection = 1;
            cout << "dPropCorrect=" << dPropCorrect << " ";
            if(dPropCorrect < 0.5) {
                dPropCorrect = 1-dPropCorrect;
                dSignCorrection = -1;
            }
            if(bVerbose) cout << "dSignCorrection=" << dSignCorrection << endl;

            if(bVerbose) {
                cout << nNumWrong << " incorrectly classified" << endl;
                cout << 100*dPropCorrect << "% correct (trivial=" << (int)100*dTrivialPropCorrect << "%)" << endl;

                if(dPropCorrect < dTrivialPropCorrect) {
                    cout << "Warning: poor performance" << endl;
                    return dTrivialPropCorrect;
                }
            }

            //cout << test_labels << endl;
            return dPropCorrect;
        } else {
            if(bVerbose && dBoundary != 0)
                cout << "Decision boundary = " << dBoundary << endl;

            const double dTotalSuccessRate = computeTotalSuccessRate(labels, test_labels, dSignCorrection, bVerbose);

            const double dBSR = computeBSR(labels, test_labels, dSignCorrection, bVerbose);

            if(bVerbose)
                cout << "BSR=" << dBSR << " dTotalSuccessRate=" << dTotalSuccessRate << endl;
            
            dPrecision = computePrecision(labels, test_labels, dSignCorrection, bVerbose);

            return dTotalSuccessRate;
        }
    }

    double computeTotalSuccessRate(const cv::Mat & labels, const cv::Mat & test_labels, double & dSignCorrection, const bool bVerbose) const {
        double dTotalSuccessRate = -1;
        double dTotalScore=0, dTotalErrorScore=0;
        for(int i=0; i < labels.rows; i++) {
            const bool bGTClass = svmClass(labels.at<float>(i));
            const double dAbsLabel = fabs(cvmGet(pClassWeights, bGTClass, 0));

            const bool bPredictedClass = svmClass(test_labels.at<float>(i));

            if(bGTClass != bPredictedClass) {
                dTotalErrorScore += dAbsLabel;
            }

            dTotalScore += dAbsLabel;
        }

        if(dTotalErrorScore < 0.5*dTotalScore) {
            dSignCorrection = 1;
            dTotalSuccessRate = (dTotalScore-dTotalErrorScore)/dTotalScore; //I think this is precisely what is being minimised in the svm
        } else {
            dSignCorrection = -1;
            dTotalSuccessRate = dTotalErrorScore/dTotalScore; //I think this is precisely what is being minimised in the svm
        }
        if(bVerbose) cout << "dSignCorrection=" << dSignCorrection << endl;

        return dTotalSuccessRate;
    }
    
    double computePrecision(const cv::Mat & labels, const cv::Mat & test_labels, const double dSignCorrection, const bool bVerbose) const {
        
        int nNumActuallyCorrect=0, nNumLabelledCorrect=0, nTotalCorrectExamples=0;
        for(int i=0; i < labels.rows; i++) {
            const double labelGT = labels.at<float>(i);
			const bool bGTClass = svmClass((float)labelGT);
            const bool bPredictedClass = svmClass((float)dSignCorrection*test_labels.at<float>(i));
            if(bPredictedClass)
            {
                nNumLabelledCorrect++;
                if(bGTClass)
                    nNumActuallyCorrect++;
            }
            if(bGTClass)
                nTotalCorrectExamples++;
        }
        const double dPrecision = (double)nNumActuallyCorrect/(double)nNumLabelledCorrect;
        const double dRecall = (double)nNumActuallyCorrect/(double)nTotalCorrectExamples;
        
        if(bVerbose)
            cout << "Precision=" << dPrecision << " Recall=" << dRecall << endl;
        
        return dPrecision;
    }
    

    double computeBSR(const cv::Mat & labels, const cv::Mat & test_labels, const double dSignCorrection, const bool bVerbose) const {
        double dBSR = 0; //2-class balanced success rate, as defined in "A User's Guide to Support Vector Machines"
        //Actually we will play with the weights in labels to adjust the cost of mistakes, so we should also compute error rate in the same way

        for(int bLabel=0; bLabel < 2; bLabel++) {

            double dErrors = 0, dExamples=0;
            for(int i=0; i < labels.rows; i++) {
                const double labelGT = labels.at<float>(i);
                const bool bGTClass = svmClass((float)labelGT);

                if(bGTClass == (bool)bLabel) {
                    const bool bPredictedClass = svmClass((float)dSignCorrection*test_labels.at<float>(i));

                    dExamples++;
                    if(bGTClass != bPredictedClass)
                        dErrors++;
                }
            }
            const double dClassSuccessRate = (dExamples-dErrors)/dExamples;
            dBSR += 0.5*dClassSuccessRate;

            if(bVerbose)
                cout << "Class " << bLabel << " success rate=" << dClassSuccessRate << endl;
        }
        return dBSR;
    }

    float svmScore(const bool bLabel) const {
        return bLabel ? 1.0f : -1.0f; //Unfortunately can't send fNegativeScore to opencv;
    }
    bool svmClass(const float fScore) const {
        return fScore>0;
    }

    class CTrainValidateFeatureSet
    {
        enum eTrainValidate { eTrain, eValidate, NUM_FEATURE_DIVS };
        const CSVMTraining * pTrainer;
        cv::Mat aFeatures[NUM_FEATURE_DIVS], aLabels[NUM_FEATURE_DIVS];
    public:
        CTrainValidateFeatureSet() {}
        CTrainValidateFeatureSet(const CSVMTraining * pTrainer, const TLabelledFeatures * aaadFeatures, const int nSubset, const int K) : pTrainer(pTrainer) {
            const bool bVerbose = false;
            std::vector<cv::Mat> aFeaturesVec[NUM_FEATURE_DIVS];
            std::vector<float> aLabelsVec[NUM_FEATURE_DIVS];
            for(int bLabel=0; bLabel < 2; bLabel++) {
                const TLabelledFeatures & aadFeatures = aaadFeatures[bLabel];

                const int nNumFeatures = (int)aadFeatures.size();
                const int nValidateBlockStart = (nNumFeatures*nSubset)/K;
                const int nValidateBlockEnd = (nNumFeatures*(nSubset+1))/K;

                for(int nFeature=0; nFeature < nNumFeatures; nFeature++) {
                    const eTrainValidate selectTV = (nFeature >= nValidateBlockStart && nFeature < nValidateBlockEnd) ? eValidate : eTrain;

                    const cv::Mat & feature = aadFeatures[nFeature];

                    CHECK_P(feature.size().area() == 0, feature, "Feature has no area");

                    aLabelsVec[selectTV].push_back(pTrainer->svmScore(bLabel));
                    aFeaturesVec[selectTV].push_back(feature);
                }
            }
            for(int eTV = 0; eTV < (int)NUM_FEATURE_DIVS; eTV++) {
                CHECK(aLabelsVec[eTV].size() == 0, "0 labels");
                CHECK(aFeaturesVec[eTV].size() == 0, "0 features");

                aLabels[eTV] = cv::Mat(aLabelsVec[eTV], true);
                aFeatures[eTV] = vectorToMat(aFeaturesVec[eTV]);

                CHECK(aLabels[eTV].size().area() == 0, "0 labels");
                CHECK(aFeatures[eTV].size().area() == 0, "0 features");

                if(bVerbose) {
                    cout << "Features: " << aFeatures[eTV] << endl;
                    cout << "Labels: " << aLabels[eTV] << endl;
                }
            }
        }

        double trainAndValidate(const cv::SVMParams & svmParams, int & nNumSVs) const {
            cv::SVM svm;
            const bool bVerbose = false;
            if(bVerbose) {
                cout << aFeatures[eTrain] << endl;
                cout << aLabels[eTrain] << endl;
            }
            svm.train(aFeatures[eTrain], aLabels[eTrain], cv::Mat(), cv::Mat(), svmParams);
            
            nNumSVs = svm.get_support_vector_count();
            
            return validate(svm);
        }
        
        int dims() const { return aFeatures[eTrain].cols; }
    private:
        double validate(const cv::SVM & svm) const {
            double dSignCorrection = 0, dPrecision=0;
            return pTrainer->propCorrect(svm, aFeatures[eValidate], aLabels[eValidate], dSignCorrection, dPrecision, 0);
        }
    };

    class CSVMParameterisation
    {
        cv::SVMParams svmParams;
        double dCVScore; //Worst 0.5 to best 1
        double dNumSVs;
    public:
        CSVMParameterisation(const double nu=-1, const double gamma=-1, CvMat * pClassWeights = 0) : dCVScore(-1), dNumSVs(-1) {
            svmParams.svm_type = SVM_TYPE;
            svmParams.kernel_type = (gamma > 0) ? cv::SVM::RBF : cv::SVM::LINEAR;
            svmParams.class_weights = pClassWeights;
            svmParams.nu = nu;
            svmParams.C = nu;
            svmParams.gamma = gamma;
        }

        void setCVScore(const double dNewCVScore, const double dNewNumSVs) {
            //CHECK(dNewCVScore < 0, "CV score not set");
            //CHECK(dCVScore >= 0, "CV score already set");
            dCVScore = dNewCVScore;
            dNumSVs = dNewNumSVs;
        }
        
        double getCVScore() const {
            return dCVScore;
        }
        
        double getNumSVs() const {
            return dNumSVs;
        }
        
        const cv::SVMParams& getSvmParams() const {
            return svmParams;
        }

        std::string toString() const {
            std::ostringstream ss;
            ss << "Nu=" << svmParams.nu << "Gamma=" <<svmParams.gamma;
            return ss.str();
        }
    };

    class CComputeMeanCov
    {
        std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > aSamples;
        std::string label;
    public:
        CComputeMeanCov(std::string label) : label(label) {}

        void addSample(cv::Mat & sample_cv) {
            Eigen::VectorXd sample(sample_cv.cols);
            for(int i=0; i<sample_cv.cols; i++)
                sample(i) = sample_cv.at<float>(i);

            aSamples.push_back(sample);
        }

        ~CComputeMeanCov() {
            const int nSamples = (int)aSamples.size();
            if(nSamples==0)
                return;

            const int nDim = (int)aSamples[0].rows();

            Eigen::VectorXd mean = Eigen::VectorXd::Zero(nDim);

            for(int i=0; i<nSamples; i++) {
                mean += aSamples[i];
            }
            mean /= nSamples;

            cout << label << ": Mean=" << mean.transpose() << endl;

            Eigen::MatrixXd var = Eigen::MatrixXd::Zero(nDim, nDim);
            for(int i=0; i<nSamples; i++) {
                Eigen::VectorXd diff = aSamples[i]-mean;
                var += diff*diff.transpose();
            }

            var /= nSamples;

            const Eigen::VectorXd SD = var.diagonal().array().sqrt();

            cout << "SD: " << SD.transpose() << endl;

            //cout << "Var: " << endl << var << endl << endl;
        }

    };

    class CKfoldTrainValidateFeatureSet
    {
        const CSVMTraining * pTrainer;
        std::vector<CTrainValidateFeatureSet> aFeatureDivisions;
        TLabelledFeatures aaadFeatureSubset[2];

        void selectFeatureSubset(const TLabelledFeatures * aaadFeatures, const CFeatureSubsetSelecter & normalisingCoeffs) {
            const bool bVerbose = false;
            for(int bLabel=0; bLabel < 2; bLabel++) {

                //std::vector<std::vector<double> > aFeaturesForMeanSD;
                CComputeMeanCov meanCov(bLabel ? "Positive" : "Negative");

                TLabelledFeatures & aadFeatureSubset = aaadFeatureSubset[bLabel];
                BOOST_FOREACH(const cv::Mat & feature, aaadFeatures[bLabel]) {
                    cv::Mat normalisedFeatureSubset=normalisingCoeffs.selectAndNormalise(feature);
                    aadFeatureSubset.push_back(normalisedFeatureSubset);

                    meanCov.addSample(normalisedFeatureSubset);

                    if(bVerbose) {
                        cout << "aadFeatureSubset[0] " << aadFeatureSubset[0]  << endl;
                        cout << "aadFeatureSubset.back() " << aadFeatureSubset.back()  << endl;
                    }

                    if(aadFeatureSubset.size() > 1) {
                        const double dDist = cv::mean(aadFeatureSubset[0] - aadFeatureSubset.back())[0];
                        if(fabs(dDist) < 1e-8) cout << "Warning: Duplicate training vectors (this is usually ok). Separation=" << dDist << endl;
                        //CHECK_P(fabs(dDist) < 1e-8, dDist, "Duplicate training vectors (this is usually ok)");
                    }
                }
                CHECK(aadFeatureSubset.size() == 0 || aadFeatureSubset.size() != aaadFeatures[bLabel].size(), "Lost size");
            }
        }

    public:
        double trainOnAll(cv::SVM & svm_final, const CSVMParameterisation & bestParameterisation, double & dSignFix, TPRLookup & adPrecision) const {
            //Select all features
            std::vector<cv::Mat> aFeaturesVec;
            std::vector<float> aLabelsVec;
            for(int bLabel=0; bLabel < 2; bLabel++) {
                BOOST_FOREACH(const cv::Mat & feature, aaadFeatureSubset[bLabel]) {
                    aFeaturesVec.push_back(feature);
                    aLabelsVec.push_back(pTrainer->svmScore(bLabel));
                }
            }
            cv::Mat allFeatures = vectorToMat(aFeaturesVec);
            cv::Mat allLabels(aLabelsVec);
            svm_final.train(allFeatures, allLabels, cv::Mat(), cv::Mat(), bestParameterisation.getSvmParams());
            
            
            double dPrecision = -1;
            for(double dBoundary = -1; dBoundary <= 1; dBoundary += 0.1) {
                pTrainer->propCorrect(svm_final, allFeatures, allLabels, dSignFix, dPrecision, dBoundary);
                
                adPrecision.first.push_back(dBoundary);
                adPrecision.second.push_back(dPrecision);
            }
            
            
            for(int i=0; i<aFeaturesVec[0].cols-1; i++) {
                outputDecisionBoundary(bestParameterisation, svm_final, aFeaturesVec[0]*0, i, i+1);
            }

            const double dScore = pTrainer->propCorrect(svm_final, allFeatures, allLabels, dSignFix, dPrecision, 0);
            return dScore;
        }

        void outputDecisionBoundary(const CSVMParameterisation & bestParameterisation, const cv::SVM & svm, cv::Mat testFeature, const int i, const int j) const {
            const std::string folder = pTrainer->getPath() + "/boundaries/";
            boost::filesystem::create_directories(folder);

            cout << "Created directory " << folder << endl;

            std::ostringstream filename;
            filename << bestParameterisation.toString() << "i=" << i << "j=" << j << ".tsv";
            const std::string fullPath = folder + filename.str();

            std::ofstream outputFile(fullPath.c_str());
            cout << "Created output file " << fullPath << endl;

            for(float dMin_i=-2.0f; dMin_i <= 2.0f; dMin_i+=0.04f)
                for(float dMin_j=-2.0f; dMin_j <= 2.0f; dMin_j+=0.04f) {
                    testFeature.at<float>(i) = dMin_i;
                    testFeature.at<float>(j) = dMin_j;
                    const float dResponse = svm.predict(testFeature, true);
                    outputFile << dMin_i << '\t' << dMin_j << '\t' << dResponse << endl;
                }
        }


        CKfoldTrainValidateFeatureSet(const CSVMTraining * pTrainer, const TLabelledFeatures * aaadFeatures, const CFeatureSubsetSelecter & normalisingCoeffs, const int K) : pTrainer(pTrainer) {

            selectFeatureSubset(aaadFeatures, normalisingCoeffs);

            for(int i=0; i<K; i++)
                aFeatureDivisions.push_back(CTrainValidateFeatureSet(pTrainer, aaadFeatureSubset, i, K));
        }

        void trainAndValidate(CSVMParameterisation & svmParams) const {
            const bool bVerbose = false;

            double dKFoldCVScore = 0, dAvNumSVs = 0;
            
            BOOST_FOREACH(const CTrainValidateFeatureSet & featureSet, aFeatureDivisions) {
                int nNumSVs = 0; 
                dKFoldCVScore += featureSet.trainAndValidate(svmParams.getSvmParams(), nNumSVs);
                dAvNumSVs += nNumSVs;
            }
            dKFoldCVScore /= (double)aFeatureDivisions.size();
            dAvNumSVs /= (double)aFeatureDivisions.size();
            
            const double dPenalty = 0.003*aFeatureDivisions[0].dims();
            
            if(bVerbose)
            { 
                cout << aFeatureDivisions.size() << "-fold cross validation score=" << dKFoldCVScore << endl;
                cout << "Penalty = " << dPenalty << endl;
            }
            svmParams.setCVScore(dKFoldCVScore-dPenalty, dAvNumSVs);
        }
    };

    std::vector<CSVMParameterisation> getHyperparamSets() const {
        std::vector<CSVMParameterisation> aParameterisations;

        std::vector<double> adNuVals, adGammaVals;
        if(SVM_TYPE == cv::SVM::NU_SVC)
        {
            adNuVals.push_back(0.001);
            adNuVals.push_back(0.002);
            adNuVals.push_back(0.0033);
            adNuVals.push_back(0.006);
            adNuVals.push_back(0.01);
            adNuVals.push_back(0.02);
            adNuVals.push_back(0.033);
            adNuVals.push_back(0.1); //Make sure we don't end up with far too many SVs
            adNuVals.push_back(0.2);
            adNuVals.push_back(0.3);
            adNuVals.push_back(0.4);
            //adNuVals.push_back(0.5);
        }
        else if(SVM_TYPE == cv::SVM::C_SVC)
        {
            for(double dPow = -5; dPow <= 15; dPow += 2)
                adNuVals.push_back(std::pow(2.0, dPow));
        }
        
        adGammaVals.push_back(-1); //linear
        
        double dLogGammaStart = -12, dLogGammaEnd = 5, dGammaStep = (dLogGammaEnd-dLogGammaStart)/15.01; //16 steps
        for(double dLogGamma = dLogGammaStart; dLogGamma < dLogGammaEnd; dLogGamma+=dGammaStep) { //a bit random results below -1
            adGammaVals.push_back(exp(dLogGamma));
        }

        BOOST_FOREACH(const double gamma, adGammaVals) {
            BOOST_FOREACH(const double nu, adNuVals) {
                aParameterisations.push_back(CSVMParameterisation(nu, gamma, pClassWeights));
            }
        }
        
        return aParameterisations;
    }

    void makeNewSubsets(const TFeatureIdxSubset & bestFeatureSubsetThisSize, std::set<TFeatureIdxSubset> & candidateSubsetsToTry, const bool bFFS) {
        if(bFFS)
            makeNewSubsets_forwards(bestFeatureSubsetThisSize, candidateSubsetsToTry);
        else
            makeNewSubsets_backwards(bestFeatureSubsetThisSize, candidateSubsetsToTry);
    }
    void makeNewSubsets_backwards(const TFeatureIdxSubset & bestFeatureSubsetThisSize, std::set<TFeatureIdxSubset> & candidateSubsetsToTry) const {
        const bool bVerbose = true;

        candidateSubsetsToTry.clear();
        if(bVerbose) cout << "Best=" << toString(bestFeatureSubsetThisSize) << endl;
        //Remove one feature in turn
        BOOST_FOREACH(const int nRemove, bestFeatureSubsetThisSize) {
            TFeatureIdxSubset newSubset;
            BOOST_FOREACH(const int nKeep, bestFeatureSubsetThisSize) {
                if(nKeep != nRemove)
                    newSubset.push_back(nKeep);
            }
            if(bVerbose) cout << "newSubset=" << toString(newSubset) << endl;
            candidateSubsetsToTry.insert(newSubset);
        }
    }

    void makeNewSubsets_forwards(const TFeatureIdxSubset & bestFeatureSubsetThisSize, std::set<TFeatureIdxSubset> & candidateSubsetsToTry) const {
        const bool bVerbose = true;
        const int nDims = aaadFeatures[false][0].cols;

        candidateSubsetsToTry.clear();
        if(bVerbose) cout << "Best=" << toString(bestFeatureSubsetThisSize) << endl;
        //Add one feature in turn
        for(int i=0; i<nDims; i++) {
            if(std::find(bestFeatureSubsetThisSize.begin(), bestFeatureSubsetThisSize.end(), i) == bestFeatureSubsetThisSize.end()) {
                TFeatureIdxSubset newSubset = bestFeatureSubsetThisSize;
                newSubset.push_back(i);
                candidateSubsetsToTry.insert(newSubset);
            }
        }
    }

    static std::string toString(const TFeatureIdxSubset & featureSubset) {
        std::ostringstream ss;
        for(int i=0; i<(int)featureSubset.size(); i++) {
            if(i>0) ss << "-";
            ss << featureSubset[i];
        }
        return ss.str();
    }

    void trainHyperparameters(const CFeatureSubsetSelecter & featureSubset, std::vector<CSVMParameterisation> & aParameterisations, CSVMParameterisation & bestParameterisationForThisSubset) {

        const std::string surfaceDir = path + "/hyperparams";
        boost::filesystem::create_directories(surfaceDir);
        const std::string surfaceName = surfaceDir + "/surface" + toString(featureSubset.getFeatureIdxSubset()) + ".tsv";
        std::ofstream surfaceTSVFile(surfaceName.c_str()); // See ~/pruning/trainSVM/contourPlot.py

        CKfoldTrainValidateFeatureSet trainAndValidateData(this, aaadFeatures, featureSubset, K);

        BOOST_FOREACH(CSVMParameterisation & parameterisation, aParameterisations) {
            TNullaryFnObj fn = boost::bind(&CKfoldTrainValidateFeatureSet::trainAndValidate, &trainAndValidateData, boost::ref(parameterisation));
            pSVMThreadpool->addJob(fn);
        }
        cout << "Training with " << aParameterisations.size() << " hyperparameterisations..." << endl;

        pSVMThreadpool->waitForAll();

        cout << "Done training" << endl;

        BOOST_FOREACH(CSVMParameterisation & parameterisation, aParameterisations) {
            if(parameterisation.getCVScore() > bestParameterisationForThisSubset.getCVScore())
                bestParameterisationForThisSubset = parameterisation;

            surfaceTSVFile << parameterisation.getSvmParams().nu << '\t' << log(parameterisation.getSvmParams().gamma) << '\t' << parameterisation.getCVScore() << '\t' << parameterisation.getNumSVs() << endl;
        }
        cout << "Best parameterisation for this subset has score " << bestParameterisationForThisSubset.getCVScore() << endl;

        //Now retrain the best on all data and save
        /*const double dFinalScoreOneSubset = trainAndValidateData.trainOnAll(svm_final, bestParameterisationForThisSubset, dSignFix);
        cout << "Score on training set: " << dFinalScoreOneSubset << endl; //but we use the CV scores to select the best subset...*/
    }

    void logResults(const TFeatureIdxSubset & bestFeatureSubsetThisSize, const CSVMParameterisation & bestParameterisationForSubsetsThisSize, std::ostream & bestResults) const {
        bestResults << bestFeatureSubsetThisSize.size() << '\t' << toString(bestFeatureSubsetThisSize) << '\t' << bestParameterisationForSubsetsThisSize.getSvmParams().nu << '\t' << bestParameterisationForSubsetsThisSize.getSvmParams().gamma << '\t' << bestParameterisationForSubsetsThisSize.getCVScore() << '\t' << endl;
    }

    //The best hyperparameters are about the same for every subset. After the first run only consider the best K (also should slightly reduce risk of overfitting)
    void filterHyperparameters(std::vector<CSVMParameterisation> & aParameterisations) {
        const int nNumToKeep = K;
        if((int)aParameterisations.size() <= nNumToKeep)
            return;

        typedef std::map<double, CSVMParameterisation, std::greater<double> > TSortedParams;
        TSortedParams aParamScores;

        BOOST_FOREACH(const CSVMParameterisation & parameterisation, aParameterisations) {
            aParamScores[parameterisation.getCVScore()] = parameterisation;
        }

        aParameterisations.clear();

        BOOST_FOREACH(const TSortedParams::value_type & parameterisationPair, aParamScores) {
            aParameterisations.push_back(parameterisationPair.second);
            if((int)aParameterisations.size() >= nNumToKeep)
                break;
        }
    }

    void setupFeatureSelection(eSVMFeatureSelectionMethod & featureSelectionMode, std::set<TFeatureIdxSubset> & candidateSubsetsToTry) const {
        const int nDims = aaadFeatures[false][0].cols;

        std::string featureSet = path + "/featureSet";
        if(boost::filesystem::exists(featureSet)) {
            std::ifstream featureFile(featureSet.c_str());
            TFeatureIdxSubset candidateSubset;
            while(!featureFile.eof()) {
                int nFeature = -1;
                featureFile >> nFeature;
                cout << "Loaded feature " << nFeature << endl;
                if(nFeature<0)
                    break; //doesn't hit eof after end for some reason...
                CHECK_P(nFeature < 0 || nFeature >= nDims, nFeature, "nFeature index OOB on load from featurefile");
                candidateSubset.push_back(nFeature);
            }
            candidateSubsetsToTry.insert(candidateSubset);
            CHECK_P(candidateSubsetsToTry.size()==0 || (int)candidateSubsetsToTry.size() > nDims, candidateSubsetsToTry.size(), "Bad load from featurefile");

            featureSelectionMode=eLoadFromFile;
        } else if(featureSelectionMode == eBFS) {
            TFeatureIdxSubset candidateSubset(nDims);
            for(int i=0; i < nDims; i++)
                candidateSubset[i] = i;

            candidateSubsetsToTry.insert(candidateSubset);
        } else if(featureSelectionMode == eFFS) {
            TFeatureIdxSubset empty;
            makeNewSubsets_forwards(empty, candidateSubsetsToTry);
        }
    }

    /* train nu and gamma and subset (CSVMParameterisation) to maximise k-fold X-validation score */
    void trainSVM(CFeatureSubsetSelecter & featureSubset, cv::SVM & svm_final, double & dSignFix, TPRLookup & PRLookup) {

        const int nDims = aaadFeatures[false][0].cols;

        std::set<TFeatureIdxSubset> candidateSubsetsToTry;
        setupFeatureSelection(featureSelectionMode,  candidateSubsetsToTry);
        std::vector<CSVMParameterisation> aParameterisations = getHyperparamSets();

        CSVMParameterisation bestParameterisationOverall;
        TFeatureIdxSubset bestFeatureSubsetOverall;

        std::string allResultsFile = path + "/" + label + "-allResults.tsv";
        std::ofstream allResults(allResultsFile.c_str());
        CHECK_P(!allResults.is_open(), allResultsFile, "File not open");

        std::string bestResultsFile = path + "/" + label + "-bestResults.tsv";
        std::ofstream bestResults(bestResultsFile.c_str());
        CHECK_P(!bestResults.is_open(), bestResultsFile, "File not open");

        for(int nSubsetSize = nDims; nSubsetSize > 0; nSubsetSize--) { //For each subset

            CSVMParameterisation bestParameterisationForSubsetsThisSize;
            TFeatureIdxSubset bestFeatureSubsetThisSize;

            BOOST_FOREACH(const TFeatureIdxSubset & anFeatureSubset, candidateSubsetsToTry) {
                CSVMParameterisation bestParameterisationForThisSubset(-1,-1);

                featureSubset.setFeatureIdxSubset(anFeatureSubset);

                trainHyperparameters(featureSubset, aParameterisations, bestParameterisationForThisSubset);

                if(bFilterHyperparams && (int)anFeatureSubset.size() > nDims/3) filterHyperparameters(aParameterisations);

                logResults(anFeatureSubset, bestParameterisationForThisSubset, allResults);

                if(bestParameterisationForThisSubset.getCVScore() > bestParameterisationForSubsetsThisSize.getCVScore()) {
                    bestParameterisationForSubsetsThisSize = bestParameterisationForThisSubset;
                    bestFeatureSubsetThisSize = anFeatureSubset;
                }
            }
            logResults(bestFeatureSubsetThisSize, bestParameterisationForSubsetsThisSize, bestResults);

            if(bestParameterisationForSubsetsThisSize.getCVScore() >= bestParameterisationOverall.getCVScore()) {
                bestParameterisationOverall = bestParameterisationForSubsetsThisSize;
                bestFeatureSubsetOverall = bestFeatureSubsetThisSize;
            }

            if(featureSelectionMode==eLoadFromFile)
                break;

            makeNewSubsets(bestFeatureSubsetThisSize, candidateSubsetsToTry, featureSelectionMode == eFFS);
        }

        cout << "Best subset has " << bestFeatureSubsetOverall.size() << " features" << endl;
        featureSubset.setFeatureIdxSubset(bestFeatureSubsetOverall);
        CKfoldTrainValidateFeatureSet bestTrainAndValidateData(this, aaadFeatures, featureSubset, K);

        const double dFinalScoreOverall = bestTrainAndValidateData.trainOnAll(svm_final, bestParameterisationOverall, dSignFix, PRLookup);
        cout << "Score on training set after retrain on all: " << dFinalScoreOverall << endl;
    }

    void computeNegWeight(const int nNumPos, const int nNumNeg) {
        //We want 1*nNumPos = fNegWeightToMakeClassesBalance*nNumNeg (as per "A User's Guide to Support Vector Machines")
        const float fNegWeightToMakeClassesBalance = (float)nNumPos/(float)nNumNeg;
        const float fNegativeScore = -fNegRelativeWeight*fNegWeightToMakeClassesBalance;

        pClassWeights = cvCreateMat(2,1,CV_32FC1);
        cvmSet(pClassWeights, 0, 0, fNegativeScore);
        cvmSet(pClassWeights, 1, 0, svmScore(true));
        //pClassWeights = cv::Mat(2,1,CV_32FC1, svmScore(true));
        //pClassWeights.at<float>(0) = fNegativeScore;

        cout << "fNegativeScore=" << fNegativeScore << endl;
    }

    std::string getPath() const {
        return path;
    }

public:
    CSVMTraining(const std::string path, const std::string label, const float fNegRelativeWeight, const eSVMFeatureSelectionMethod featureSelectionMode, const bool bFilterHyperparams) : path(path), label(label), pSVMThreadpool(CThreadpool_base::makeThreadpool(bMT ? K : 1)), fNegRelativeWeight(fNegRelativeWeight), pClassWeights(0), featureSelectionMode(featureSelectionMode), bFilterHyperparams(bFilterHyperparams) {
        boost::filesystem::create_directories(path);

        CHECK(fNegRelativeWeight <= 0, "Bad fNegRelativeWeight");
    }

    /*
     * Repeat:
     *   1: compute best boosting test
     *   2: filter features
     *
     * 3: compute normalising coeffs
     *
     * 4: normalise
     *
     * 5: Train backward or forward feature selection + hyperparameters to max. CV score
     * For each subset:
     *   select subset
     *   For each parameter set:
     *     getCVScore(subset, parameters)
     *
     * 6: resolve return sign
     *
     * save SVM with best CV score
     */
    virtual ~CSVMTraining() {

        int nNumPos=(int)aaadFeatures[true].size(), nNumNeg = (int)aaadFeatures[false].size();
        cout << "Training from " << nNumPos << " positive and " << nNumNeg << " negative examples" << endl;
        
        if(nNumPos < 20 || nNumNeg < 20)
        {
            cout << "INSUFFICIENT TRAINING DATA" << endl;
            return;
        }
        
        TBoosterStates aBoosterStates = findBoosterStates();
        //features in aaadFeatures are now filtered

        CFeatureSubsetSelecter featureSubset;
        cv::SVM svm;
        double dSignFix = 0;
        nNumPos=(int)aaadFeatures[true].size(), nNumNeg = (int)aaadFeatures[false].size();
        cout << nNumPos << " positive and " << nNumNeg << " negative examples after boosting" << endl;
        
        TPRLookup PRLookup;

        if(nNumPos > 0 && nNumNeg > 0) {
            //We've still got something to train from...
            computeNegWeight(nNumPos, nNumNeg);

            featureSubset.findNormalisingCoeffs(aaadFeatures);

            trainSVM(featureSubset, svm, dSignFix, PRLookup);
        } else {
            cout << "Boosting left no training data. Save boosting-only classifier." << endl;
        }
        CSavedSVMState savedState(path, label, aBoosterStates, featureSubset, svm, dSignFix, PRLookup);

        if(pClassWeights) cvReleaseMat(&pClassWeights);
    }
    
    static bool equal(const cv::Mat & M1, const cv::Mat & M2)
    {
        for(int i=0;i<M1.rows;i++)
            if(M1.at<float>(i) != M2.at<float>(i) )
                return false;
                
        return true;
    }

    virtual void addTrainingFeature(CSVMFeature_base * pFeature, const bool bLabel) {
        boost::mutex::scoped_lock lock(mxLockToAdd);

        const bool bVerbose = false;
        
        const cv::Mat & entireFeature=pFeature->getEntireFeature();
        
        BOOST_FOREACH(const cv::Mat & sample, aaadFeatures[bLabel])
        {
            if(equal(sample, entireFeature))
            {
                if(bVerbose)
                    cout << "Duplicate feature " << bLabel << " " << aaadFeatures[bLabel].back() << endl;
                
                return;
            }             
        }
        
        aaadFeatures[bLabel].push_back(entireFeature);
        
        if(bVerbose)
            cout << "Added feature " << bLabel << " " << aaadFeatures[bLabel].back() << endl;
    }
};

CSVMTraining_base * CSVMTraining_base::makeSVMTraining(const std::string path, const std::string label, const float fNegRelativeWeight, const eSVMFeatureSelectionMethod featureSelectionMode, const bool bFilterHyperparams)
{
    return new CSVMTraining(path, label, fNegRelativeWeight, featureSelectionMode, bFilterHyperparams);
}

const double CSVMClassifier_base::NO_PRECISION = -1;