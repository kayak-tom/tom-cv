#include "LevenbergMarquardt.h"
#include <cmath>
#include <iostream>

//#include "SpeedTest.h"

using namespace std;

namespace grc {

LevenbergMarquardt::~LevenbergMarquardt()
{
    cvReleaseMat(&J_vec_);
    cvReleaseMat(&J_);
    cvReleaseMat(&J_inv_);
    cvReleaseMat(&delta_);
    cvReleaseMat(&delta_scaled_);
    cvReleaseMat(&gradient_vec_);
};

double LevenbergMarquardt::maxAbs(const CvMat * M) 
{
	double g_min, g_max;
    cvMinMaxLoc(M, &g_min, &g_max);
    if(g_max+g_min > 0)
        return g_max;

    return -g_min;
}

double LevenbergMarquardt::maxDiagVal(const CvMat * M)  
{
	double maxDiag = -HUGE;
    for (int i = 0; i < M->cols; i++) {
        double diagVal = cvmGet(M, i,i);
		if (diagVal > maxDiag) {
			maxDiag = diagVal;
		}
	}
    return maxDiag;
}

double LevenbergMarquardt::sumSquare(const CvMat * M) 
{
    double SS = 0;
	double * data = M->data.db;
    int matSize = M->rows*M->cols;
	for (int i = 0; i < matSize; i++) {
        double elem = data[i];
        SS += elem*elem;
	}
    return SS;
}

double LevenbergMarquardt::frobeniusNorm(const CvMat * M)  
{
    return sqrt(sumSquare(M));
}

void LevenbergMarquardt::addLMDampingFactor(CvMat * M, double u) 
{
	double * data = M->data.db;
    int matSize = M->cols;
    int step = matSize+1;
	for (int i = 0; i < (matSize*matSize); i+= step) {
        data[i] += u;
	}
}

LevenbergMarquardt::LevenbergMarquardt(TransformErrFunction &errFn, const int maxLMIterations) : errFn_(errFn), maxLMIterations_(maxLMIterations)
{
    const bool VERBOSE = false;

    const double e1 = 1e-15;
	const double e2 = 1e-15;
	const double e3 = 1e-15;
	const double t = 1e-3;

	double v = 2;

    //TransformErrFunction::cParamVals * originalParams = errFn.getVals();

    TransformErrFunction::cParamVals * params = errFn.getVals();
    numParams_ = params->size();

	J_ = cvCreateMat(numParams_, numParams_, CV_64FC1);
	J_inv_ = cvCreateMat(numParams_, numParams_, CV_64FC1);
	J_vec_ = cvCreateMat(numParams_, 1, CV_64FC1);
    delta_ = cvCreateMat(numParams_, 1, CV_64FC1);
    delta_scaled_ = cvCreateMat(numParams_, 1, CV_64FC1);
    gradient_vec_ = cvCreateMat(numParams_, 1, CV_64FC1);

	double initErr = errFn_.evaluateReprojErr();

    computeJacobian();
	//cout << J << endl;

	cout << fixed;
	cout.precision(8);

    cvScale(J_vec_, gradient_vec_, -initErr);
    double gradientSize = maxAbs(gradient_vec_);
	bool stop = (gradientSize <= e1);
	if (stop) {
		cout << "Stopping because gradient is very small" << endl;
	}
	double maxDiag = maxDiagVal(J_);
	double damp = t*maxDiag;

    int k = 0;
	for (; !stop && k < maxLMIterations_; k++) {
        double rho = -1.0;
		while (!stop && rho <= 0.0 ) {
			rho = -1.0;
			//Matrix J_inv_ = (J_ + u*I).i()*g;
           
            addLMDampingFactor(J_, damp);
            //CStopWatch s;
            //s.startTimer();
            cvInv(J_, J_inv_, CV_LU); // CV_SVD_SYM is ok as is positively defined (I think) + symetric
            //Takes 0.0018s for CV_LU or 0.004s for CV_SVD_SYM (96 params)
            //s.stopTimer();
            //cout << "Inversion took " <<s.getElapsedTime() << " secs, " << numParams_ << " params\n";

            cvMatMul(J_inv_, gradient_vec_, delta_); //delta is the innovation vector

			if ( frobeniusNorm(delta_) <= e2*params->frobeniusNorm() && k > 1 ) {
				stop = true;
				cout << scientific;
				cout << "Stopping because delta_ is very small " << frobeniusNorm(delta_) << " <= " << e2*params->frobeniusNorm() << endl;
				cout << fixed;

			} else {
				TransformErrFunction::cParamVals * newParams = new TransformErrFunction::cParamVals(params, delta_); //initErr = errFn_.evaluateReprojErr();
                // rho is the improvement
                errFn_.assignAllVals(newParams);
                double newErr = errFn_.evaluateReprojErr();
                
                //delta_scaled_ = delta_*damp + gradient_vec_
                cvScaleAdd(delta_, cvScalar(damp), gradient_vec_, delta_scaled_);

                rho = (initErr - newErr)/cvDotProduct(delta_, delta_scaled_);
				if (rho > 0.0 ) { //Yes this step is good so apply it
                    initErr = newErr;
                    delete params; params = newParams;
					
                    computeJacobian();

					if ( VERBOSE ) {
						cout << "Iteration " << k << " error = " << initErr << endl;
					}
					cvScale(J_vec_, gradient_vec_, -initErr);

                    if (maxAbs(gradient_vec_) <= e1) {
						stop = true;
						cout << "Stopping because gradient is very small" << endl;
						if ( VERBOSE )
						{
							//cout << "Final residuals:" << endl;
							//for (int i=0;i<x.Nrows();i++)
							//	cout << i << ' ' << x.element(i) << '\t' << evaluate(params).element(i) << '\t' << initErr.element(i) << endl;
						}

					} else if (initErr <= e3) {
						stop = true;
						cout << "Stopping because error is very small" << endl;
					}
					double um = 1-pow(2.0*rho-1.0, 3);
					if (um < 0.33) {
						um = 0.33;
					}
					damp *= um;
					v = 2.0;
				} else {// this step is bad so don't apply it
                    delete newParams;
                    errFn_.assignAllVals(params); // reset parameters
					damp *= v;
					v *= 2.0;
				}
			}
		}
	}
    ////Output the actual changes
    //params->printDiff(originalParams);
    //delete originalParams;

    std::cout << "LM error=" << initErr << " of initial error after " << k << " iters; " << params->size() << " parameters adjusted\n";

    delete params;
}


	void LevenbergMarquardt::computeJacobian() {

        for(int i = 0; i < numParams_; i++)
        {
            double df_dxi = errFn_.pd(i);
            cvmSet(J_vec_, i, 0, df_dxi);
        }
        cvGEMM(J_vec_, J_vec_, 1.0, NULL, 0.0, J_, CV_GEMM_B_T);
	}
}