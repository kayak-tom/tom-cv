#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include "util/convert.h"
#include "levMarNumerical.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <math.h>

template<class TSparseMat>
void denseToSparse(const Eigen::MatrixXd & D, TSparseMat & S) {
    S.resize((int) D.rows(), (int) D.cols());
    S.reserve(20 * (int) S.cols()); //each row contains 20 elements

    for (int j = 0; j < D.cols(); ++j) {
        for (int i = 0; i < D.rows(); i++) // with increasing i
            if (D(i, j) != 0)
                S.insert(i, j) = D(i, j);
    }
    S.makeCompressed();
}

template<class TSparseMat>
void pp_sparse_warning(const TSparseMat & S) {
    if (S.nonZeros() == 0)
        return;

    double dPropFill = S.nonZeros() / ((double) S.cols() * (double) S.rows());
    if (dPropFill > 0.15) {
        cout << S.nonZeros() << " nonzeros, " << 100.0 * dPropFill << "% fill " << (int) (S.cols() * dPropFill) << " nonzeros per row" << endl;
    }

    //spy(S);
}

void pp(const std::string label, const Eigen::MatrixXd & J, const int nSubmatSize)
{
    if(J.size() == 0)
    {
        cout << " [empty]" << endl;
        return;
    }
    
    int nNonZeros = 0;
    for(int r=0;r<J.rows(); r++)
        for(int c=0;c<J.cols(); c++)
        {
            const double dFabs = fabs(J(r,c));
            if(dFabs > 1e-12)
            {
                nNonZeros++;
            }
        }
    int nR=-1, nC=-1;
    cout << label << ": size " << J.rows() << "x" << J.cols() << "=" << J.size() << " Nonzeros: " << nNonZeros << " Min: " << J.minCoeff(&nR, &nC) << "[" << nR << "," << nC << "]" << " Max: " << J.maxCoeff(&nR, &nC) << "[" << nR << "," << nC << "]" << endl;

    const Eigen::MatrixXd topLeft = J.block(0,0,std::min<int>((int)J.rows(), nSubmatSize), std::min<int>((int)J.cols(), nSubmatSize));
    if(J.cols() == 1)
    {
        if(topLeft.rows() < 10)
        {
            cout << topLeft.transpose() << endl;
        }
        else
        {
            for(int i=0; i < topLeft.rows(); i++)
                cout << "[" << i << "] " << topLeft(i) << " ";
            cout << endl;
        }
    }
    else
        cout << topLeft << endl;
}
void CLevMar::computeDerivatives(const Eigen::VectorXd &x) {
    boost::timer time; //for mat solve timing dense/sparse

    if (bMakeSparseJ)
        THROW("Not implemented yet");
    
    const int nInputs = function.inputs(), nValues = function.values();
    
    Eigen::VectorXd Dresiduals(nValues);
    
    for (int nParam = 0; nParam < nInputs; nParam++) {
        function.analyticDeriv(x, Dresiduals, nParam);
        J.col(nParam) = Dresiduals;
    }

    grad = J.transpose() * residual;

    if(nVerbose>0) REPEAT(1, pp_sparse_warning(spJ));

    if (nVerbose>=0 && nLMIter == 1 && time.elapsed() > 0.1)
        REPEAT(100, cout << "Compute J: " << time.elapsed() << endl);
}

CLMFunction::eLMSuccessStatus CLevMar::computeDerivativesNumerically(const Eigen::VectorXd &x, const bool bForwards) {
    boost::timer time; //for mat solve timing dense/sparse

    if (IS_DEBUG) {
        Eigen::VectorXd resids_temp = residual;
        if(robustFunction(x, resids_temp, false, -1) == CLMFunction::eLMFail)
            return CLMFunction::eLMFail;

        if (resids_temp != residual) {
            if(residual.size() < 1000)
            {
                cout << residual.transpose() << endl << endl;
                cout << resids_temp.transpose() << endl << endl;
                cout << (resids_temp - residual).transpose() << endl << endl;
            }
            cout << "diff=" << (resids_temp - residual).norm() << endl;
            THROW("computeDerivativesNumerically not supplied with derivatives needed (is some renormalisation of x going on?)");
        }
    }

    if (bMakeSparseJ) {
        CHECK(spJ.Options & Eigen::RowMajorBit, "Sparse fill will fail");
        spJ.resize((int) function.values(), (int) function.inputs());
        spJ.reserve(20 * (int) spJ.cols()); //each row contains 20 elements
    }

    //double dMaxDelta = std::max<double>(dMinDelta, 0.005);
    double delta = bForwards ? dMinDelta : -dMinDelta;

    Eigen::VectorXd resids_plus = residual;
    Eigen::VectorXd x_plus = x;
    Eigen::VectorXd diff;

    const int nInputs = function.inputs(), nValues = function.values();
    for (int nParam = 0; nParam < nInputs; nParam++) {

        const double delta_inv = 1.0 / delta;
        if (nParam > 0)
            x_plus(nParam - 1) -= delta;

        x_plus(nParam) += delta;
        if(robustFunction(x_plus, resids_plus, false, nParam) == CLMFunction::eLMFail)
            return CLMFunction::eLMFail;

        if (IS_DEBUG && nParam == 12) //Check param updates working properly
        {
            Eigen::VectorXd resids_temp = residual;
            robustFunction(x_plus, resids_temp, false);

            CHECK((resids_temp - resids_plus).squaredNorm() > 1e-10, "Param incremental updating failed.");

            robustFunction(x, residual, false);
        }

        if (bMakeSparseJ) {
            spJ.startVec(nParam);

            diff = resids_plus - residual;
            for (int nResid = 0; nResid < nValues; nResid++) {
                if (diff(nResid) != 0) {
                    double dDeriv = diff(nResid) * delta_inv;
                    spJ.insertBack(nResid, nParam) = dDeriv; //Much faster than insert
                    resids_plus(nResid) = residual(nResid);
                    //cout << nResid << " " << nParam << ": " << dDeriv << endl;
                    //spJT.insert(nParam, nResid) = dDeriv;
                }
            }
        } else {
            J.col(nParam) = (resids_plus - residual).array() * delta_inv;

            resids_plus = residual;
        }

        //cout << "Max deriv:" << fjac.col(nParam).maxCoeff() << " using delta " << delta << endl;
    }
    //cout << fjac << endl << endl;

    if (bMakeSparseJ) {
        spJ.makeCompressed();
    } else if (method == eLMSparse)
        denseToSparse(J, spJ);

    if (spJ.rows() > 0)
        grad = spJ.transpose() * residual;
    else
        grad = J.transpose() * residual;

    if(nVerbose>0) REPEAT(1, pp_sparse_warning(spJ));

    if (nVerbose>=0 && nLMIter == 1 && time.elapsed() > 0.1)
    {
        REPEAT(100, cout << "Compute J: " << time.elapsed() << endl);
    }
    
    //PPMAT(J);
    return CLMFunction::eLMSuccess;
}

CLevMar::eIterState CLevMar::dampMarquardt(const double lambda, TParamVector & JTJ_diag) {
    JTJ_diag *= (1 + lambda);
    double dMinDiag = JTJ_diag.minCoeff();
    if (dMinDiag <= 0) {
        double dMean = fabs(JTJ_diag.mean());
        JTJ_diag.array() += lambda*dMean; //Whether marquardt or not, to stop it going singular
        if(lambda*dMean == 0)
        {
            cout << "JTJ has a zero diagonal (probably an error, but probably exactly at the minimum => converge)" << endl;
            return eConverged;
        }
    }
    return eDescending;
}

CLevMar::eIterState CLevMar::dampMarquardt(const double lambda, Eigen::MatrixXd & JTJ) {
    JTJ.diagonal().array() *= (1 + lambda);
    double dMinDiag = JTJ.diagonal().array().minCoeff();
    if (dMinDiag <= 0) {
        double dMean = fabs(JTJ.diagonal().array().mean());
        JTJ.diagonal().array() += lambda*dMean; //Whether marquardt or not, to stop it going singular
        if(lambda*dMean == 0)
        {
            cout << "JTJ has a zero diagonal (probably an error, but probably exactly at the minimum => converge)" << endl;
            return eConverged;
        }
    }
    return eDescending;
}

CLevMar::eIterState CLevMar::dampMarquardt(const double lambda, Eigen::SparseMatrix<double> & JTJ) {
    double dSumDiagonalSq = 0;
    for (int k = 0; k < JTJ.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(JTJ, k); it; ++it) {
            if (it.row() == it.col()) {
                dSumDiagonalSq += sqr(it.value());
                it.valueRef() *= (1 + lambda);
                if (it.value() <= 0)
                    it.valueRef() = 1; //doesn't matter as long as it is nonzero
            }
        }
    }
    
    return (dSumDiagonalSq == 0) ? eConverged : eDescending;
}

CLevMar::CLevMar(CLMFunction & f, const bool bVerbose, const double dMinDelta, const eMethod method) : function(f), nVerbose(bVerbose ? 1 : 0), dMinDelta(dMinDelta), dMinStepLength(0.00001), pseudoHuber_t(-1), method(method), residual(function.values()) {
    paramUpdateVec.resize(function.inputs());
    paramUpdateVec.setConstant(1); //initial delta
    bMakeSparseJ = (method == eLMSparse || method == eCGD);
}

void CLevMar::setMinStepLength(const double dMinStepLength_in)
{
    dMinStepLength = dMinStepLength_in;
}

void CLevMar::setPseudoHuberSD(const double t)
{
    pseudoHuber_t = t;
}

const Eigen::VectorXd & CLevMar::residuals() const {
    return residual;
}

const double CLevMar::RMSError() const
{
    CHECK(residual.size() == 0, "No residuals (yet?)");
    return sqrt(residual.squaredNorm() / residual.size());
}

CLMFunction::eLMSuccessStatus CLevMar::robustFunction(const Eigen::VectorXd & params, Eigen::VectorXd &resids, bool bVerbose, const int nParamChanged)
{
    //return function.function(params, resids, bVerbose, -1);
    if(pseudoHuber_t<=0)
        return function.function(params, resids, bVerbose, nParamChanged);
    
    const double b_sq = sqr(pseudoHuber_t);
    
    if(nParamChanged == -1)
    {
        CLMFunction::eLMSuccessStatus result = function.function(params, resids, bVerbose, -1);
        
        //Reweight residuals
        //const Eigen::ArrayXd d_abs = resids.array().abs() + 1e-12;

        //C(delta) = 2*b^2*(sqrt(1+(delta/b)^2) - 1);
        //cout << "Before " << resids.transpose() << endl;
        resids = (2 * b_sq * (sqrt(1 + (resids.array() / pseudoHuber_t).eval().square()) - 1)).sqrt().eval();
        //cout << "Weights: " << weights.transpose() << endl;
        //resids.array() /= weights;
        //cout << "After " << resids.transpose() << endl;

        return result;
    }
    else //only a few residuals will change
    {
        const double NO_RESID = -99999;
        Eigen::VectorXd tempResids = Eigen::VectorXd::Constant(resids.rows(), NO_RESID);
        CLMFunction::eLMSuccessStatus result =  function.function(params, tempResids, bVerbose, nParamChanged);
        
        for(int nResid = 0; nResid < resids.rows(); nResid++)
        {
            if(tempResids(nResid) != NO_RESID)
                resids(nResid) = sqrt(2 * b_sq * (sqrt(1 + sqr(tempResids(nResid) / pseudoHuber_t)) - 1));
        }

        return result;
    }
}

#ifdef _WIN32
#  ifndef _DEBUG
#    pragma optimize("x", on)
#  endif
#endif
CLevMar::eIterState CLevMar::levMarIter(Eigen::VectorXd & params, double & lambda, const bool bMarquardt, double & dErr) {

    boost::timer totaltime; //for total LM time
    const double eps = sqr(dMinStepLength);

    if (method == eLMSparse) {

        if (!bMakeSparseJ)
            denseToSparse(J, spJ);

        if(nVerbose>0) REPEAT(1, pp_sparse_warning(spJ));

        spJTJ = spJ.transpose() * spJ;

        if(nVerbose>0) REPEAT(1, pp_sparse_warning(spJ));
    }

    for (;;) {
        boost::timer time; //for mat solve timing dense/sparse

        if (method == eLM) {
            JTJ = J.transpose() * J; //TODO:Move outside loop

            if (nVerbose>=0 && nLMIter == 1 && time.elapsed() > 0.1) {
                cout << "Compute JTJ: " << time.elapsed() << "s" << endl;
                time.restart();
            }

            if (bMarquardt)
            {
                if(dampMarquardt(lambda, JTJ) == eConverged)
                {
                    //if (nVerbose>=0)
                    cout << "J with zero diagonal\n" << J << endl; //Usually an error
                    return eConverged;
                }
            }       
            else
                JTJ.diagonal().array() += lambda;

            if(IS_DEBUG) CHECK(JTJ.diagonal().array().minCoeff() <= 0, "LM damping failed to 'nonzero the diagonal");

            if(nVerbose>1) {
                cout << "Iter " << nLMIter << " lambda=" << lambda << ":\n";
                pp("J", J, 1000);
                pp("JTJ", JTJ, 1000);
                pp("residual", residual, 1000);
                pp("param update", paramUpdateVec);
            }
            
            paramUpdateVec = JTJ.llt().solve(grad); //generally faster than .inverse() (todo: write fixed size version and restore?)
        } else if (method == eLMSparse) {
            if(dampMarquardt(lambda, spJTJ) == eConverged)
                return eConverged;
                

#if EIGEN_VERSION_AT_LEAST(3,0,9)
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
            cg.compute(spJTJ);
            paramUpdateVec = cg.solve(grad);
            if (cg.iterations() == cg.maxIterations()) {
				if(nVerbose>=0) 
				{
					std::cout << "Max # iterations hit: " << cg.iterations();
					std::cout << " estimated error: " << cg.error() << std::endl;
				}
            }
#else
            THROW("Need Eigen 3.1");
#endif
        } else {
            //LM Diag
            const int nInputs = function.inputs();
            TParamVector JTJ_diag(nInputs);
            for (int i = 0; i < nInputs; i++) {
                JTJ_diag(i) = J.col(i).squaredNorm();
            }
            if(dampMarquardt(lambda, JTJ_diag) == eConverged)
                return eConverged;

            paramUpdateVec = grad.array() / JTJ_diag.array();
        }

        if (nVerbose>=0 && nLMIter == 1 && time.elapsed() > 0.1) {
            cout << "Solve JTJ: " << time.elapsed() << endl;
        }

        if (IS_DEBUG && std::isnan(paramUpdateVec.sum())) {
            
            if(method == eLMSparse)
            {
                cout << "spJ: " << spJ << endl << endl;
                cout << "spJTJ: " << spJTJ << endl << endl;
            }
            else
            {
                cout << "J: " << J << endl << endl;
                cout << "JTJ: " << JTJ << endl << endl;
            }

            cout << "Residuals: " << residuals().transpose() << endl << endl;
            cout << "paramUpdateVec: " << paramUpdateVec.transpose() << endl << endl;
            cout << "Lambda: " << lambda << endl;
            
            THROW("Error--paramUpdateVec is nan");
        }

        const double dUpdateStepLengthSq = paramUpdateVec.squaredNorm();

        if (nVerbose>0)
            cout << "Step size: " << sqrt(dUpdateStepLengthSq) << " dErr= " << dErr << " step: " << paramUpdateVec.transpose().segment(0, std::min<int>((int) paramUpdateVec.size(), 10)) << "..." << endl;

        const TParamVector paramsNew = params - paramUpdateVec;

        const CLMFunction::eLMSuccessStatus successStatus = robustFunction(paramsNew, residual, nVerbose>0); 
        
        if(nVerbose>0)
        {
            pp("Residuals", residual, 100);
        }

        const double dErrNew = (successStatus == CLMFunction::eLMSuccess) ? residual.squaredNorm() : dErr; //Force step back on failure
        if (dErrNew < dErr) {
            lambda *= 0.25;
            dErr = dErrNew;
            params = paramsNew;

            if (dUpdateStepLengthSq < eps) //Step is tiny. 
            {
                if(nVerbose>0)
                    cout << "Lambda " << lambda << " - converged on small step of length " << sqrt(dUpdateStepLengthSq) << "\n";

                return eConverged;
            }
            return eDescending;

        } else {
            lambda *= 4;
        }

        if (dUpdateStepLengthSq < eps) {
            if(nVerbose>0)
                cout << "Lambda " << lambda << " - converged on small step of length " << sqrt(dUpdateStepLengthSq) << "\n";
            return eConverged; // derivatives wrong near minimum, and damped until lambda is massive
        } 
        
        
        if(lambda > 100000) { //If we end up here then fn is not really approx quadratic. Does happen though, eg. cane reconstruction
            if(nVerbose>0)
            {
                cout << "Converge on large lambda: Param Update Length = " << sqrt(dUpdateStepLengthSq);
                cout << ", min delta = " << dMinDelta << " Iter = " << nLMIter << endl;

                PPMAT(paramUpdateVec);
                PPMAT(params);
                if(!bMakeSparseJ)
                {
                    PPMAT(grad);
                    PPMAT(J);
                    PPMAT(JTJ);
                }
            }
            return eConverged;
        }
    }
    
    if (nVerbose>=0 && totaltime.elapsed() > 0.1) {
        cout << "LM iter total: " << totaltime.elapsed() << "s" << endl;
        totaltime.restart();
    }

    return eDescending;
}

CLevMar::eIterState CLevMar::minLineSearch(TParamVector & params, const TParamVector & delXn, double & alpha, double & dInitErr) {
    //double dInitErr = residual.squaredNorm();
    double dStepLength = delXn.norm();

    //Interval bracketing
    //double adAlpha[3] = {0,0,0};
    //double adVals[3] = {dInitErr, HUGE, HUGE};

    SVal LB(0, dInitErr);
    SVal UB;
    std::vector<SVal> aMidpoints;

    //First increase alpha until we have bracketed minimum

    for (;;) {
        TParamVector params_temp = params + alpha*delXn;
        double dNewErr = function.sumSquare(params_temp, residual, false);
        if (dNewErr > dInitErr) //We have found the far side of the interval
        {
            //We have bracketed solution
            UB = SVal(alpha, dNewErr);
            break;
        } else {
            aMidpoints.push_back(SVal(alpha, dNewErr));
        }
        alpha *= 2;
    }

    SVal smallest(HUGE, HUGE);

    SVal secondSmallest(HUGE, HUGE);

    //Now set smallest

    // choose 2 intermediate points
    if (aMidpoints.size() > 0) {
        for (int i = 0; i < (int) aMidpoints.size(); i++) {
            if (aMidpoints[i].val < smallest.val) {
                secondSmallest = smallest;
                smallest = aMidpoints[i];
            } else if (aMidpoints[i].val < secondSmallest.val) {
                secondSmallest = aMidpoints[i];
            }
        }

        //need another one. Split the first interval
        if (aMidpoints.size() == 1) {
            alpha = smallest.alpha * 0.5;
            TParamVector params_temp = params + alpha*delXn;
            double dNewErr = function.sumSquare(params_temp, residual, false);

            secondSmallest = SVal(alpha, dNewErr);
            if (secondSmallest.val < smallest.val)
                std::swap(secondSmallest, smallest);
        }

        //Now choose 2 better endpoints
        if (smallest.alpha < secondSmallest.alpha) //largest is near the UB of the interval
            UB = secondSmallest; //Move UB down
        else
            LB = secondSmallest; //Move LB up
    } else {
        //Need to compute them:
        //Drop alpha until find a point within interval
        for (;;) {
            alpha = UB.alpha * 0.5; //Todo use the fact that the OF is > 0

            TParamVector params_temp = params + alpha*delXn;
            double dNewErr = function.sumSquare(params_temp, residual, false);

            if (dNewErr > dInitErr)//Still a UB
            {
                UB = SVal(alpha, dNewErr);
            } else //we have a midpoint
            {
                smallest = SVal(alpha, dNewErr);
                break;
            }
        }
    }

    const int MAX_ITERS = 10;
    for (int nIter = 0; nIter < MAX_ITERS; nIter++) {
        //Split largest interval
        if (smallest.alpha - LB.alpha > UB.alpha - smallest.alpha) {
            //Largest interval is below smallest
            alpha = 0.5 * (smallest.alpha + LB.alpha);
        } else {
            alpha = 0.5 * (UB.alpha + smallest.alpha);
        }


        TParamVector params_temp = params + alpha*delXn;
        double dNewErr = function.sumSquare(params_temp, residual, false);

        secondSmallest = SVal(alpha, dNewErr);
        if (secondSmallest.val < smallest.val)
            std::swap(secondSmallest, smallest);

        //Now choose 2 better endpoints
        if (smallest.alpha < secondSmallest.alpha) //largest is near the UB of the interval
            UB = secondSmallest; //Move UB down
        else
            LB = secondSmallest; //Move LB up

        if ((UB.alpha - LB.alpha) * dStepLength < dMinDelta * 2)
            break;
    }

    params = params + alpha * delXn;

    dInitErr = function.sumSquare(params, residual, nVerbose>0); //To make sure residuals are reset

    double dStepSize = dStepLength * alpha;
    if (dStepSize < dMinDelta)
        return eConverged;
    else
        return eDescending;
}

void CLevMar::checkGrad(const TParamVector & grad) const {

    int nPosGradients = 0;
    for (int i = 0; i < grad.rows(); i++) {
        if (grad(i) >= 0)
            nPosGradients++;
    }
    if ((nPosGradients - grad.rows() / 2) / (double) grad.rows() > 0.25)
        cout << nPosGradients << " / " << function.inputs() << " positive gradients\n";
}



CLevMar::eIterState CLevMar::CGDIter(TParamVector & params, const bool bReset, double & alpha, double & dErr) {//Notation from http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
    //checkGrad(grad);

    TParamVector gradXn = -grad;

    if (delXn.size() != function.inputs() || bReset) {
        gradXn_take_1 = Eigen::VectorXd::Zero(function.inputs());
        delXn = gradXn_take_1;
    }

    const double beta_num = gradXn.transpose() * (gradXn - gradXn_take_1);
    const double beta_denom = gradXn_take_1.squaredNorm();
    double beta = 0;
    if (beta_denom > 0 && beta_num > 0 /*&& nResetCount > 0*/)
        beta = beta_num / beta_denom;

    delXn = (gradXn + beta * delXn).eval();


    if (nVerbose>0) {
        cout << "beta " << beta << endl;
        cout << "Max: " << delXn.maxCoeff() << endl;
        cout << "Min: " << delXn.minCoeff() << endl;
        cout << "RMS: " << sqrt(delXn.squaredNorm() / delXn.size()) << endl;
    }
    delXn /= delXn.norm();

    gradXn_take_1 = -grad;

    return minLineSearchNR(params, delXn, dErr);
    //return minLineSearch(params, delXn, alpha, dErr);
}

//Also updates params, and evaluated residual at minimum
//See numerical recipes section 5.7
void CLevMar::testLineSearchDelta(TParamVector & params, const TParamVector & delXn, double & dErr) {

    Eigen::VectorXd params_plus(function.inputs());

    for (double delta = dMinDelta; delta < 10000*dMinDelta; delta*=1.25) {

        double f_zero = dErr;
        params_plus = params - delta*delXn;
        double f_minus1 = function.sumSquare(params_plus, residual);

        params_plus += 2 * delta*delXn;
        double f_plus1 = function.sumSquare(params_plus, residual);

        //Take a factor delta out of each to avoid the delta^2
        //const double Df = 0.5 * (f_plus1 - f_minus1);
        double D2f = ((f_plus1 - f_zero) - (f_zero - f_minus1)) / sqr(delta);
        
        cout << delta << ": " << D2f << endl;
    }
}
    

CLevMar::eIterState CLevMar::minLineSearchNR(TParamVector & params, const TParamVector & delXn, double & dErr) {
    
    //testLineSearchDelta(params, delXn, dErr);
    
    const int MAX_NR_ITERS = 5;
    Eigen::VectorXd params_plus(function.inputs());

    const double dNorm = delXn.norm();
    double dTotalUpdate = 0;

    const double delta = pow(dMinDelta, 0.75); //Heuristic appears to work well. Numerical recipies says macheps^0.25

    int nNRIter = 0;
    double dInitUpdate = -1;
    
    for (; nNRIter < MAX_NR_ITERS; nNRIter++) {
        /*double f_minus1 = dErr;
        params_plus = params + delta*delXn;
        double f_zero = function.sumSquare(params_plus, residual);

        params_plus += delta*delXn;
        double f_plus1 = function.sumSquare(params_plus, residual);*/

        double f_zero = dErr;
        params_plus = params - delta*delXn;
        double f_minus1 = function.sumSquare(params_plus, residual, (nVerbose>0) && (nNRIter == 0));

        params_plus += 2 * delta*delXn;
        double f_plus1 = function.sumSquare(params_plus, residual);

        //Take a factor delta out of each to avoid the delta^2
        const double Df = 0.5 * (f_plus1 - f_minus1);
        double D2f = ((f_plus1 - f_zero) - (f_zero - f_minus1)) / delta;

        if(IS_DEBUG) CHECK(std::isnan(Df) || std::isinf(Df), "Inf or nan computing Df");
        if(IS_DEBUG) CHECK(std::isnan(D2f) || std::isinf(D2f), "Inf or nan computing D2f");

        if (Df == 0) {
            //This is actually good convergence I think
            /*
            cout << "Converge on Df=0" << endl;
            
            params_plus -= 20 * delta*delXn;
            double f_minus = function.sumSquare(params_plus, residual);
            if(f_minus < f_zero)
                cout << "Would not have converged on larger delta" << endl;
             */
            nNRIter = MAX_NR_ITERS;
            break;
        }

        if (D2f == 0) {
            D2f = 1;
            cout << "Set D2f=1" << endl;
        }

        if (nVerbose>0) {
            cout << "f_minus1=" << f_minus1 << endl;
            cout << "f_zero=" << f_zero << endl;
            cout << "f_plus1=" << f_plus1 << endl;

            cout << "Df=" << Df / delta << endl;
            cout << "D2f=" << D2f / delta << endl;
        }
        double update = Df / D2f;
        if (std::isnan(update) || std::isinf(update)) {
            cout << "Converge on NaN/inf update" << endl;
            return eConverged;
        }

        if (nNRIter == 0 && Df >= 0) {
            if(nVerbose>0) cout << "Converging when fn does not decrease in gradient direction\n";
            nNRIter = MAX_NR_ITERS;
            break;
        } else if (dTotalUpdate - update < 0) {
            if(nVerbose>0) cout << "dTotalUpdate - update < 0, choosing damped update...\n";
            update = (nNRIter == 0) ? -1 : (dTotalUpdate * 0.5);
        }

        for (;;) {
            params_plus = params - update*delXn;

            dErr = function.sumSquare(params_plus, residual);
            if(IS_DEBUG) CHECK(std::isnan(dErr) || std::isinf(dErr), "Inf or nan computing residual");

            if(nVerbose>0) cout << "NR iter " << nNRIter << " Update " << update << " err=" << dErr << endl;

            if (dErr < f_zero) {
                if(nNRIter == 0)
                    dInitUpdate = fabs(update);
                
                dTotalUpdate -= update;
                params = params_plus;
                break;
            }

            if (fabs(update) * dNorm < delta)//function.inputs() *
            {
                /*if (nNRIter == 0)
                {
                    checkGrad(grad);
                    cout << "Warning: Failed to find small step which reduces cost" << endl;
                }*/
                nNRIter = MAX_NR_ITERS;
                break;
            }

            update *= 0.5; //Step backwards until hit a minimum
        }
        if(fabs(update) < 0.1*dInitUpdate) //Numerical recipes: don't want to spend too long doing line search
            break;
    }

    if (nNRIter >= MAX_NR_ITERS)
        dErr = function.sumSquare(params, residual);

    if (dTotalUpdate * dNorm < delta)
        return eConverged;
    else
        return eDescending;
}


double CLevMar::minimise(TParamVector & params, const int MAX_ITERS) {
    if(IS_DEBUG) CHECK(params.size() != function.inputs(), "Param vector has wrong size");

    static int s_nTotalRuns = 0;
    s_nTotalRuns++;
    static int s_nTotalIters = 0;
    static double s_dTotalReduction = 0;

    if(robustFunction(params, residual, nVerbose>0) == CLMFunction::eLMFail)
    {
        if(nVerbose>=0)
        {
            cout << "LM failed on first function call--initial parameters are invalid" << endl;
        }
        return HUGE;
    }        

    double dErr = residual.squaredNorm();
    double dErr_init = dErr;
    if (std::isnan(dErr)) {
        cout << "Params: " << params.transpose() << endl;
        cout << "Residuals: " << residual.transpose() << endl;
        THROW("Error is nan")
    }

    const bool bMarquardt = true;
    double lambda = bMarquardt ? 0.5 : 100;
    double alpha = 1;

    if (!bMakeSparseJ)
    {
        //Need dense J as well
        J.resize(function.values(), function.inputs());
        int nJSize = ((int) J.size() * 8) / 1024;
        if (nJSize > 100000)
            cout << "J size=" << nJSize << "kb" << endl;
    }

    CLevMar::eIterState iterState = eDescending;
    double dPreviousErr = dErr;
    for (nLMIter = 0; nLMIter < MAX_ITERS && iterState == eDescending; nLMIter++) { //Error will not actually go to 0
        s_nTotalIters++;
        
        dPreviousErr = dErr;

        if (function.useAnalyticDerivatives())
        {
            computeDerivatives(params);
        }
        else
        {
            if(computeDerivativesNumerically(params) == CLMFunction::eLMFail)
            {
                if(nVerbose>=0)
                    cout << "LM failed to compute derivatives--we have stepped to the edge of the allowed parameter space" << endl;
                return HUGE;
            }
        }
        
        if (method == eCGD) {
            bool bReset = nLMIter % 5 == 0;
            iterState = CGDIter(params, bReset, alpha, dErr);
            if (iterState == eConverged && !bReset) {
                iterState = CGDIter(params, true, alpha, dErr);
                if (nVerbose>0 && iterState != eConverged)
                    cout << "Continued descent after reset" << endl;
            }
        } else {
            iterState = levMarIter(params, lambda, bMarquardt, dErr);
        }

        if(nVerbose>0)
            cout << "Iter=" << nLMIter << ", lambda=" << lambda << ", err=" << dErr << endl;

        //cout << "params = " << params.transpose() << endl;
    }

    double dPropInitialError = 0;
    if (dErr_init > 0.00001)
        dPropInitialError = dErr / dErr_init;

    s_dTotalReduction += dPropInitialError;

	const double dPropFinal = (dPreviousErr - dErr)/(((dErr_init - dErr) != 0) ? (dErr_init - dErr) : 1);
    const bool bOutput = (nVerbose>0) || (iterState == eDescending && dPropFinal > 0.01 && nVerbose>=0);

    if (bOutput) {
        if (iterState == eConverged)
            cout << "Convergence after " << nLMIter << " iterations, ";
        else
            cout << "No convergence after " << nLMIter << " iterations, ";

        cout << "error = " << dErr;

		if (dErr_init > 0.00001 || nVerbose>0) {
            cout << " = " << dPropInitialError << " of initial error";
            cout << " " << dPropFinal << " of the reduction was on the final iteration";
        }

        cout << endl;
        cout << "RMS error = " << sqrt(dErr / residual.size()) << endl;
        if (params.size() < 25) cout << "params = " << params.transpose() << endl;
    }

    if (nVerbose>=0 && s_nTotalRuns % 2000 == 0) {
        cout << "Runs: " << s_nTotalRuns << endl;
        cout << "Iterations per run: " << (double) s_nTotalIters / s_nTotalRuns << endl;
        cout << "Reduction in error per run: " << (double) s_dTotalReduction / s_nTotalRuns << endl;
    }

    return dErr;
}

//Simple gradient descent to find a maximum

double CLevMar::maximiseGD(TParamVector & params) {
    return optimiseGD(params, false);
}

//Simple gradient descent to find a minimum

double CLevMar::minimiseGD(TParamVector & params) {
    return optimiseGD(params, true);
}


double CLevMar::optimiseGD(TParamVector & params, const bool bMinimise) {
    if(IS_DEBUG) CHECK(function.values() != 1, "For GD, function should have 1 value for now")
    TResidVector residual(function.values());

    robustFunction(params, residual, nVerbose>0);

    double dScale = bMinimise ? 1 : -1;

    double dVal = dScale * residual(0);
    const double dInitVal = dVal;
    //const double eps = sqr(0.000001);

    double gamma = 2;

    int nIter = 0;
    const int MAX_ITERS = 150;
    for (; nIter < MAX_ITERS; nIter++) {
        TJMatrix J(function.values(), function.inputs());
        CHECK(std::isnan(params.sum()), "params is nan");
        computeDerivativesNumerically(params);
        CHECK(std::isnan(J.sum()), "J is nan");
        
        Eigen::VectorXd Dir = J.transpose();
        if(J.norm() == 0)
            cout << "Warning: zero derivatives" << endl;
        else
            Dir /= J.norm();


        //Now line-search (todo: improve...):
        TParamVector bestParams = params * 0, newParams = params - gamma * Dir;

        bool bImprovement = false;

        for (;;) {
            newParams = params - gamma * Dir;
            CHECK(std::isnan(newParams.sum()), "params is nan");
            
            robustFunction(newParams, residual, false);
            
            const double dNewVal = dScale * residual(0);

            if (dVal > dNewVal) {
                dVal = dNewVal;
                bestParams = newParams;
                gamma *= 2;
                if(nVerbose>0) cout << "New val " << dNewVal << endl;
                bImprovement = true;
            } else {
                gamma *= 0.5;
                break;
            }
        }//while(residual(0) < dVal)//Keep increasing gamma to find lower minima to

        if (!bImprovement) {
            for (;;) {
                newParams = params - gamma * Dir;
                robustFunction(newParams, residual, false);
                const double dNewVal = dScale * residual(0);

                if (dVal > dNewVal) {
                    dVal = dNewVal;
                    bestParams = newParams;
                    bImprovement = true;
                    if(nVerbose>0) cout << "New val " << dNewVal << endl;
                    break;
                } else {
                    gamma *= 0.5;
                    if (gamma < dMinDelta)
                        break;
                }
            }//while(residual(0) < dVal)//Keep increasing gamma to find lower minima to
        }

        if (!bImprovement) {
            cout << "Converged in " << nIter << " iterations\n";
            cout << "Improvement of " << (dInitVal - dVal) / fabs(dInitVal) << "\n";
            return dVal;
        } else {
            params = bestParams;
        }
    }

    cout << "Improvement of " << (dInitVal - dVal) / fabs(dInitVal) << "\n";
    cout << "No convergence after " << nIter << " iterations\n";
    return dVal;
}


void CLevMar::testSmoothness(const TParamVector & params, const double dDelta, const int nParam) {
    double dParamValStart = params(nParam) / 5;
    double dParamValEnd = params(nParam)*5 + dDelta * 10;

    TParamVector params_test = params;
    int i = 0, nNumTurningPoints = 0;
    double dOldVal = 0;

    enum eMode {
        eIncreasing, eDecreasing, eFixed
    };
    eMode oldMode = eFixed;

    for (double dParam = dParamValStart; dParam < dParamValEnd && i < 1000; dParam += dDelta, i++) {
        try {
            params_test(nParam) = dParam;
            double dVal = function.sumSquare(params_test, residual);

            eMode mode;
            if (dVal > dOldVal)
                mode = eIncreasing;
            else if (dVal == dOldVal)
                mode = eFixed;
            else
                mode = eDecreasing;

            if (i >= 1 && mode != oldMode)
                nNumTurningPoints++;

            //cout << mode << " " << dVal << endl;

            oldMode = mode;
            dOldVal = dVal;
        } catch (...) {
            break;
        }
    }

    cout << "Param " << nParam << ", Delta=" << dDelta << ", turning points=" << nNumTurningPoints << endl;
}

void CLevMar::testSmoothness(const TParamVector & params) {
    for (double dDelta = 1e-10; dDelta < 10; dDelta *= 4)
        for (int nParam = 0; nParam < std::min<int>(10, function.inputs()); nParam++) {
            testSmoothness(params, dDelta, nParam);
        }
}

//Compute delta minimising the difference between forward and backward derivatives
double CLevMar::normalisedDiff(const TParamVector & params)
{
    computeDerivativesNumerically(params, true);
    TJTJMatrix J_forwards = J;
    computeDerivativesNumerically(params, false);
    return sqrt((J-J_forwards).squaredNorm()/(J.rows()*J.cols()));
}

//Compute delta minimising the difference between forward and backward derivatives
double CLevMar::computeOptimalDelta(const TParamVector & params)
{
    J.resize(function.values(), function.inputs());
    
    double dBestDelta = dMinDelta, dMinDiff = HUGE;

    robustFunction(params, residual, false);
    
    for(dMinDelta = 1e-10; dMinDelta < 1e-2; dMinDelta *= 1.25)
    {
        double dDiff = normalisedDiff(params);
        
        if(dMinDiff > dDiff)
        {
            dMinDiff = dDiff;
            dBestDelta = dMinDelta;
        }
        cout << dMinDelta << " " << dDiff << endl;
    }
    dMinDelta = dBestDelta;
    
    return dBestDelta;
}

Eigen::VectorXd CLevMar::minimiseFunction(CLMFunction & lmFunction, const bool bVerbose)
{
    CLevMar LM(lmFunction, bVerbose);
    Eigen::VectorXd params = lmFunction.init();
    LM.minimise(params);
    return params;
}

/*void CLMFunction::setLargeResids(Eigen::VectorXd &resids, const int nParamChanged, const bool bVerbose) { 
    if(nVerbose>0) 
        cout << "Setting large residuals to force step back" << endl;
    
    if(nParamChanged >= 0)
        THROW("WARNING: setLargeResids while computing numerical derivatives (we have either stepped up against some threshold, or have started at an invalid position, or have failed to step back from a bad value)");
        
    resids.setConstant(10000);
}*/
