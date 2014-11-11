/* 
 * File:   houghForE.cpp
 * Author: tom
 * 
 * Created on 13 September 2011, 3:28 PM
 */

#include "houghForE.h"
#include "geom/geom.h"
#include <Eigen/Core>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "geom/geom_eigen.h"
#include "boost/math/distributions/binomial.hpp"
#include "util/random.h"
#include "models.h"

template<typename TFloat, int ROWS = 3 > //Can set ROWS=4 for possible speedup
        class CPFForE_int : public CHoughForE {
    typedef Eigen::Matrix<TFloat, 1, ROWS> TFullVec;
    typedef Eigen::Matrix<TFloat, 3, 1 > T3Vec;
    typedef Eigen::Matrix<TFloat, ROWS, 3 > TRMat;

    static const int NUM_HYPS = 100;

    static const TFloat binomialPMF(double n, double p, double N) {
        if (n < 0 || n > N)
            return 0;

        boost::math::binomial_distribution<TFloat> dist(N, p);

        return boost::math::pdf(dist, n);
    }
    
    static const TFloat binomialCMF(double n, double p, double N) {
        if (n < 0)
            return 0;

        //boost::math::binomial_distribution<TFloat> dist(N, p);
        //return boost::math::cdf(dist, n);
        
        TFloat dProb = 0;
        for(int i=1; i<=n; i++)
            dProb += binomialPMF(n-i,p,N);
        
        dProb /= N;
        return dProb;
    }
    //Represents a hypervolume of essential matrices

    class CEssentialMat_int {
        C3dRotation q;
        C3dPoint t;

        TRMat E;

        TFloat dVotes, dProb, dCumulativeProb;

        T3Vec vt;
    public:

        void init_test(const TRMat & E_in)
        {
            E= E_in;
        }
        void init(const C3dRotation & q_old, const C3dPoint & t_old, const TFloat dSD)
        {
            q = C3dRotation(q_old[0]+CRandom::Normal(0, dSD), 
                            q_old[1]+CRandom::Normal(0, dSD),
                            q_old[2]+CRandom::Normal(0, dSD),
                            q_old[3]+CRandom::Normal(0, dSD));
                    
            t = t_old;
            t.addNoise(dSD);
            
            setupE();
        }
        
        CEssentialMat_int() {
            q.setRandom();
            t.setRandomNormal();
            
            setupE();
        }
        
        void setupE()
        {
            reset();
                    
            t.normalise();
            t.asVector(vt);

            makeE(q, t, E);
        }

        void reset() {
            dVotes = 0;
        }
        
        TFloat samsonsErr(const T3Vec &x, const TFullVec &xp) const
        {
            T3Vec xpE = xp * E;
            T3Vec Ex = E * x;

            TFloat numerator = sqr(xpE.dot(x));
            TFloat denom = xpE.segment(0, 2).squaredNorm() + Ex.segment(0, 2).squaredNorm();

            return numerator / denom;
        }

        void observePoint(const T3Vec &x, const TFullVec &xp, const TFloat dMaxAngularError) {
            //x,xp are normalised. The residual should be angle from solution plane squared
            TFloat dErr = samsonsErr(x, xp);

            dVotes += dErr < dMaxAngularError ? 1 : 0;

            //cout << dErr << endl;
        }
        
        void setModel(C3x3MatModel & model) const
        {
            for(int r=0;r<3;r++)
                for(int c=0;c<3;c++)
                    model(r,c) = E(r,c);
        }

        TFloat votes() const {
            return dVotes;
        }

        void pp() const {
            cout << "Votes: " << dVotes << " Prob: " << dProb << ", t=" << t << ", q=" << q << endl;
        }
        
        TFloat SSDFromBest(const Eigen::Matrix3d & E_exact) const {
            return EssentialMat_SSD(E, E_exact, false);
        }
        
        void pp(const Eigen::Matrix3d & E_exact) const {
            pp();
            cout << "SSD=" << SSDFromBest(E_exact) << endl;
            //cout << E << endl;
        }

        TFloat computeProb(const TFloat dTotalVotes, const TFloat NUM_POINTS) {
            TFloat P_votes_if_correct = 1;
            TFloat P_votes_if_incorrect = 1;
            if(false)
            {
                const TFloat P = dTotalVotes/(NUM_HYPS * NUM_POINTS);
                P_votes_if_correct = binomialCMF(dVotes, P, NUM_POINTS);
                P_votes_if_incorrect = binomialPMF(dVotes, P, NUM_POINTS);
            }
            
            P_votes_if_correct = dVotes;

            dProb = P_votes_if_correct / P_votes_if_incorrect;

            return dProb;
        }

        TFloat normaliseProb(const TFloat dTotalProb_inv, const TFloat dCumulativeProb_in) {

            dProb *= dTotalProb_inv;
            //cout << dVotes << endl;
            //pp();
            dCumulativeProb = dCumulativeProb_in + dProb;
            
            return dProb;
        }
        
        TFloat cumulativeProb() const { return dCumulativeProb; }
        
        void makeNewHyp(CEssentialMat_int & newHyp, const double dRadius) const
        {
            newHyp.init(q, t, dRadius);
        }

    };


    static const int MAX_ITERS = 10;
    CEssentialMat_int aaHypotheses[MAX_ITERS][NUM_HYPS];

    class CRotPred {
    public:

        bool operator()(const TRMat& lhs, const TRMat& rhs) const {
            return lhs(0, 0) < rhs(0, 0);
        }
    };

public:

    CPFForE_int() {
    }

    int findRTHough(const T2dPoints & points0, const T2dPoints & points1, const Eigen::Matrix3d & E_exact, const CMask & mask_exact, C3x3MatModel & model, CMask & mask) {

        for (int i = 0; i < NUM_HYPS; i++)
            aaHypotheses[0][i].reset();

        const TFloat NORMAL_LOCALISATION_ERR_SQ=sqr(0.0005);
        const int NUM_POINTS = points0.size();

        std::vector<T3Vec> ax;
        std::vector<TFullVec> axp;
        ax.reserve(NUM_POINTS);
        axp.reserve(NUM_POINTS);
        
        for (int i = 0; i < NUM_POINTS; i++) {
            T3Vec x;
            
            x(0) = points0[i].getX();
            x(1) = points0[i].getY();
            x(2) = 1;
            ax.push_back(x);
            
            TFullVec xp;
            xp.setZero();
            xp(0) = points1[i].getX();
            xp(1) = points1[i].getY();
            xp(2) = 1;
            
            axp.push_back(xp);
        }
        
        CEssentialMat_int E_test;
        bool bTest = true;
        if(bTest)
            E_test.init_test(E_exact);
        
        TFloat dMaxAngularError_sq = sqr(0.1); //Should set this so that on 1st iteration, about 2000 votes (100^0.8 * 50) for 50 points, 100 hyps.
        //TFloat dMaxAngularError_sq = sqr(0.15); //Should set this so that on 1st iteration, about 12500 votes (1000^0.8 * 50) for 50 points, 1000 hyps.
        
        for (int nIter = 0; nIter < MAX_ITERS; nIter++) {
            
            CEssentialMat_int * aHypotheses = aaHypotheses[nIter];
            
            const int NUM_POINTS = points0.size();

            for (int i = 0; i < NUM_POINTS; i++) {

                for (int nHyp = 0; nHyp < NUM_HYPS; nHyp++) {
                    aHypotheses[nHyp].observePoint(ax[i], axp[i], dMaxAngularError_sq);
                }
                
                if(bTest)
                    E_test.observePoint(ax[i], axp[i], dMaxAngularError_sq);
            }
            int nBestHyp = -1;
            double dBestVotes = 0, dTotalVotes = 0;
            for (int i = 0; i < NUM_HYPS; i++) {
                if (aHypotheses[i].votes() > dBestVotes) {
                    dBestVotes = aHypotheses[i].votes();
                    nBestHyp = i;
                }
                //cout << i << " " << aHypotheses[i].votes() << endl;
                //aHypotheses[i].pp();

                dTotalVotes += aHypotheses[i].votes();
            }
            
            aHypotheses[nBestHyp].setModel( model );
            for (int i = 0; i < (int)points0.size(); i++) {
                mask[i] = aHypotheses[nBestHyp].samsonsErr(ax[i], axp[i]) < dMaxAngularError_sq ? true : false;
            }

            aHypotheses[nBestHyp].pp(E_exact);
            mask.pp();
            
            /* //POOR indicator of best bin
             * TFloat dClosest = HUGE;
            int nCorrectHyp = -1;
            for (int i = 0; i < NUM_HYPS; i++) {
                TFloat dDist = aHypotheses[i].SSDFromBest(E_exact);
                if(dDist < dClosest)
                {
                    dClosest = dDist;
                    nCorrectHyp = i;
                }
            }*/
            
            cout << "Best bin has " << aHypotheses[nBestHyp].votes() << " votes" << endl;
            //cout << "Correct bin has " << aHypotheses[nCorrectHyp].votes() << " votes" << endl;
            
            if(bTest)
            {
                cout << "E_test bin has " << E_test.votes() << " votes" << endl;
                int nInliersInBest = 0;
                int nIntersection = 0;
                for (int i = 0; i < (int)points0.size(); i++) {
                    if(mask[i] && (E_test.samsonsErr(ax[i], axp[i]) < dMaxAngularError_sq))
                        nIntersection++;
                    
                    if(mask_exact[i] && mask[i])
                        nInliersInBest++;
                    
                }

                cout << "Intersection size: " << nIntersection << endl;
                cout << "InliersInBest: " << nInliersInBest << " " << ((TFloat)nInliersInBest/aHypotheses[nBestHyp].votes()) << endl;
                
                E_test.reset();
            }
            
            cout << "Average votes: " << dTotalVotes / NUM_HYPS << endl;
            
            if(nIter == MAX_ITERS-1 || NORMAL_LOCALISATION_ERR_SQ == dMaxAngularError_sq)
                break;            

            TFloat dTotalProb = 0;//, dNumInliers = 0.5 * (TFloat) NUM_POINTS;

            for (int i = 0; i < NUM_HYPS; i++) {
                dTotalProb += aHypotheses[i].computeProb(dTotalVotes, NUM_POINTS);
            }

            TFloat dTotalProb_inv = 1.0 / dTotalProb, dTotalProbSq = 0;
            dTotalProb = 0;
            for (int i = 0; i < NUM_HYPS; i++) {
                TFloat dProb = aHypotheses[i].normaliseProb(dTotalProb_inv, dTotalProb);
                dTotalProb += dProb;
                
                dTotalProbSq += sqr(dProb);
            }
            
            //Resample:
            int nOldHyp = 0;
            static const TFloat NUM_HYPS_inv = 1.0/NUM_HYPS;
            for (int i = 0; i < NUM_HYPS; i++) {
                double dThresh = (i+1)*NUM_HYPS_inv - 1e-8;
                
                while(aHypotheses[nOldHyp].cumulativeProb() < dThresh)
                    nOldHyp++;
                
                aHypotheses[nOldHyp].makeNewHyp(aaHypotheses[nIter+1][i], sqrt(dMaxAngularError_sq));
            }
            
            const TFloat dEffectiveNumParticles = 1.0/dTotalProbSq; //2* because otehrwise underestimates by about 2x
            
            TFloat reductionInVolume = dEffectiveNumParticles*NUM_HYPS_inv;
            TFloat reductionInRadius = pow(reductionInVolume, 2.0/5.0);
            
            dMaxAngularError_sq *= reductionInRadius;
            if(dMaxAngularError_sq < NORMAL_LOCALISATION_ERR_SQ)
                dMaxAngularError_sq = NORMAL_LOCALISATION_ERR_SQ; //About 3px, shouldn't go lower than this
        }

        return mask.countInliers();
    }
};

template<typename TFloat, int ROWS = 3 > //Can set ROWS=4 for possible speedup
        class CHoughForE_int : public CHoughForE {
    typedef Eigen::Matrix<TFloat, ROWS, 1 > TFullVec;
    typedef Eigen::Matrix<TFloat, 3, 1 > T3Vec;
    typedef Eigen::Matrix<TFloat, ROWS, 3 > TRMat;

    static const int NUM_ROTATIONS = 1000;
    static const int NUM_TRANSLATIONS = 1000;

    TRMat aRotations[NUM_ROTATIONS];
    T3Vec aTranslations[NUM_TRANSLATIONS];

    class CTransPred {
    public:

        bool operator()(const T3Vec& lhs, const T3Vec& rhs) const {
            return lhs(0) < rhs(0);
        }
    };

    class CRotPred {
    public:

        bool operator()(const TRMat& lhs, const TRMat& rhs) const {
            return lhs(0, 0) < rhs(0, 0);
        }
    };

public:

    CHoughForE_int() {
        for (int i = 0; i < NUM_ROTATIONS; i++) {
            aRotations[i].setZero();

            C3dRotation q;
            q.setRandom();
            Eigen::Matrix3d R;
            q.asMat(R);

            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                    aRotations[i](r, c) = R(r, c);
        }

        for (int i = 0; i < NUM_TRANSLATIONS; i++) {
            aTranslations[i].setRandom();
            aTranslations[i] /= aTranslations[i].norm();
        }

        std::sort(aRotations, aRotations + NUM_ROTATIONS, CRotPred());
        std::sort(aTranslations, aTranslations + NUM_TRANSLATIONS, CTransPred());
    }

    int findRTHough(const T2dPoints & points0, const T2dPoints & points1, const Eigen::Matrix3d & E_exact, const CMask & mask_exact, C3x3MatModel & model, CMask & mask) {
        cv::Mat im(NUM_ROTATIONS, NUM_TRANSLATIONS, CV_8UC1);
        im.setTo(cv::Scalar(0));

        TFullVec xp;
        xp.setZero();
        xp(2) = 1;
        T3Vec x;
        x(2) = 1;
        for (int i = 0; i < points0.size(); i++) {
            x(0) = points0[i].getX();
            x(1) = points0[i].getY();
            xp(0) = points1[i].getX();
            xp(1) = points1[i].getY();

            x /= x.norm();
            xp /= xp.norm();

            for (int nRot = 0; nRot < NUM_ROTATIONS; nRot++) {
                T3Vec Rxp = aRotations[nRot] * xp; //NB transpose 
                T3Vec t_perp = Rxp.cross(x);

                //Easier alg: iterate around a circle, quantise components (?)
                for (int nTrans = 0; nTrans < NUM_TRANSLATIONS; nTrans++) {
                    TFloat resid = fabs(aTranslations[nTrans].dot(t_perp));
                    int nResid = (int)(10.0 * sqr(resid));
                    im.at<uchar > (nRot, nTrans) += nResid;
                }
            }
        }

        cv::imshow("hough", im);

        cv::waitKey(0);
        return 0;
    }
};

CHoughForE * CHoughForE::makeHough() {
    return new CPFForE_int<double>;
    return new CHoughForE_int<double>;
}
