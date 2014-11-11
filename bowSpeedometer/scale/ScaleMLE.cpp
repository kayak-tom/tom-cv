/*
 * ScaleMLE.cpp
 *
 *  Created on: 23 Apr 2010
 *      Author: tom
 */

#include "ScaleMLE.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/QR>

#include "util/dynArray.h"
#include <algorithm>

using namespace Eigen;

const int NUM_PARAMS = 2;
const int NUM_MODEL_VARS = 2;

typedef Eigen::Matrix<double, NUM_PARAMS, 1> TParamVector;
typedef Eigen::Matrix<double, NUM_MODEL_VARS,1> TResidVector;
typedef Eigen::Matrix<double, NUM_MODEL_VARS, NUM_PARAMS> TJMatrix;
typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;

CScaleMLE::CScaleMLE() {
	// TODO Auto-generated constructor stub

}

//We are minimising the 2 derivatives of the LL simultaneously

//Got to start very close to minimum because -> infinity the derivative also goes to 0.

//TODO: Weighted??
void getMeanVar(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, double & dRange, bool bVerbose)
{
	dMean = dVar = 0;
	double dMin = HUGE, dMax=0;

	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pd = aLengths.begin(); pd != aLengths.end(); pd++)
	{
		if(pd->hasSLAMLength())
		{
			double dLength = pd->getSLAMLength();
			dMean += dLength;
			if(dMin > dLength) dMin = dLength;
			if(dMax < dLength) dMax = dLength;

			cout << "M:"<<dLength << " V:" << pd->getVar() << endl;
		}
	}

	dRange = dMax - dMin;

	dMean /= aLengths.size();

	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pd = aLengths.begin(); pd != aLengths.end(); pd++)
	{
		if(pd->hasSLAMLength())
		{
			double diff = (pd->getSLAMLength() - dMean);
			dVar += diff*diff;
		}
	}

	dVar /= (aLengths.size());
}

double getLogLikelihood(const CBoWSpeedo::CBoWObject::TLengths & aLengths, const double m, const double s_2)
{
	//const double s = sqrt(s_2);
	double LL = 0;

	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pMeasurement = aLengths.begin(); pMeasurement != aLengths.end(); pMeasurement++)
	{
		if(pMeasurement->hasSLAMLength())
		{
			const double x = pMeasurement->getSLAMLength();
			const double si_2 = pMeasurement->getVar();

			const double K2_inv = 1.0/(si_2 + s_2);

			LL -= 0.5*(sqr(x-m)*K2_inv + /*log(2Pi) + */ log(si_2 + s_2)); //Easy to optimise out log...
		}
	}

	return LL;
}

void getDLogLikelihood(TResidVector & resid, const CBoWSpeedo::CBoWObject::TLengths & aLengths, const double m, const double s_2)
{

	resid.setZero();

	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pMeasurement = aLengths.begin(); pMeasurement != aLengths.end(); pMeasurement++)
	{
		if(pMeasurement->hasSLAMLength())
		{
			const double x = pMeasurement->getSLAMLength();
			const double si_2 = pMeasurement->getVar();

			const double K2_inv = 1.0/(si_2 + s_2);

			resid(0) += (x-m)*K2_inv; //This is dLL/dm
			resid(1) += 0.5*(sqr(x-m) - (si_2 + s_2))*sqr(K2_inv);
			//s * (sqr(x-m)*sqr(K2_inv) - K2_inv); //This is dLL/ds
		}
	}
}

//Use gradient-descent. TODO: 2nd order...
bool CScaleMLE::getMLEScaleParams_GD(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose)
{
	double dRange_unused=-1;
	getMeanVar(aLengths, dMean, dVar, dRange_unused, bVerbose);

	const double EPS = 1e-4;
	const int MAX_ITERS = 200;
	const double DROP=0.5;
	double LAMBDA = 2;

	if(dVar < EPS)
		dVar = EPS;

	double LL = getLogLikelihood(aLengths, dMean, dVar);

	TParamVector grad; grad.setConstant(1);

	for(int iter = 0;
			iter < MAX_ITERS && LAMBDA > EPS && grad.squaredNorm() > EPS;
			iter++)
	{
		getDLogLikelihood(grad, aLengths, dMean, dVar);

		/*cout << "Gradient: " << grad.transpose() << endl;

		double LL_plus_dm = getLogLikelihood(aLengths, dMean+EPS, dVar);
		cout << "Grad m: " << (LL_plus_dm-LL)/EPS;
		double LL_plus_ds = getLogLikelihood(aLengths, dMean, dVar+EPS);
		cout << ", grad s: " << (LL_plus_ds-LL)/EPS << endl;*/

		for(int iter2 = 0; LAMBDA > EPS; iter2++)
		{
			double dMean_temp = dMean + grad(0) * LAMBDA;
			double dVar_temp = dVar + grad(1) * LAMBDA;

			if(dVar_temp < EPS)
			{
				dVar_temp = EPS; //Keep positive
				if(LAMBDA < EPS && fabs(grad(0) < EPS))
				{
					cout << "Variance going negative\n";
					return false;
				}
			}

			{
				double LL_new = getLogLikelihood(aLengths, dMean_temp, dVar_temp);
				if(LL <= LL_new)
				{
					dMean = dMean_temp;
					dVar = dVar_temp;
					LL=LL_new;
					break;
				}
			}
			LAMBDA *= DROP;
		}

	}
	return LAMBDA < EPS || grad.squaredNorm() < EPS;
}

//this is the 2nd derivatives
void getDerivs(TJMatrix & J, const CBoWSpeedo::CBoWObject::TLengths & aLengths, const double m, const double s_2)
{
	J.setZero();

	//const double s = sqrt(s_2);
	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pMeasurement = aLengths.begin(); pMeasurement != aLengths.end(); pMeasurement++)
	{
		if(pMeasurement->hasSLAMLength())
		{
			const double x = pMeasurement->getSLAMLength();
			const double si_2 = pMeasurement->getVar();

			const double K2_inv = 1.0/(si_2 + s_2);

			J(0,0) += -K2_inv;
			const double numerator = (0.5*(si_2 + s_2) - sqr(m-x));
			J(1,1) += numerator * cube(K2_inv);

			const double dDmDs = (m-x) * sqr(K2_inv);
			J(0,1) += dDmDs;
		}
    }
    J(1,0) = J(0,1);
}

//Use gradient-descent. 2nd order...
bool CScaleMLE::getMLEScaleParams_2ndOrderGD(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose)
{
	bool res = getMLEScaleParams_2ndOrderGD_int(aLengths, dMean, dVar, bVerbose);
	if(res)
	{
		cout << "Mean: " << dMean << endl;
		cout << "Var: " << dVar << endl;

		double dLL = getLogLikelihood(aLengths, dMean, dVar);
		double dLL1 = getLogLikelihood(aLengths, dMean+0.1, dVar);
		double dLL2 = getLogLikelihood(aLengths, dMean-0.1, dVar);
		double dLL3 = getLogLikelihood(aLengths, dMean, dVar+0.1);
		double dLL4 = (dVar-0.1>0) ? getLogLikelihood(aLengths, dMean, dVar-0.1) : -HUGE;

		CHECK(dLL<dLL1,"Gradient-descent failed")
		CHECK(dLL<dLL2,"Gradient-descent failed")
		CHECK(dLL<dLL3,"Gradient-descent failed")
		CHECK(dLL<dLL4,"Gradient-descent failed")
		if(dVar < 0)
		{
			cout << "Setting res false as variance is negative\n";
			res = false;
		}
	}

	cout << "\nSURFACE=[[";
	for(double dM=0; dM < 6; dM += 0.2)
	{
		if(dM > 0)
			cout << "]; [";
		for(double dS=-2; dS < 6; dS += 0.2)
		{
			double dLL = getLogLikelihood(aLengths, dM, dS);
			if(dS > -2)
				cout << ',';

			cout << dLL;
		}
	}
	cout << "]];\n" << endl;

	double dMtemp, dVtemp;

	getBayesianScaleParams_NG(aLengths, dMtemp, dVtemp, bVerbose);
	return res;
}

bool CScaleMLE::getMLEScaleParams_2ndOrderGD_int(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose)
{
	double dRange_unused=-1;
	getMeanVar(aLengths, dMean, dVar, dRange_unused, bVerbose);
	dVar = sqr(dRange_unused);
	cout << "Init mean: " << dMean << endl;
	cout << "Init var: " << dVar << endl;

	const double EPS = 1e-4;
	const int MAX_ITERS = 20;
	const double DROP=0.5;
	double LAMBDA = 0.9;

	double improvement = HUGE;

	if(dVar < EPS)
		dVar = EPS;

	double LL = getLogLikelihood(aLengths, dMean, dVar);

	TParamVector grad; grad.setConstant(1);
	bool bConverged = false;

	for(int iter = 0;
			iter < MAX_ITERS && !bConverged;
			iter++)
	{
		getDLogLikelihood(grad, aLengths, dMean, dVar);

		/*{
			TParamVector gradM;
			getDLogLikelihood(gradM, aLengths, dMean+EPS, dVar);
			gradM -= grad;
			gradM /= EPS;

			cout << gradM << "=dJdm\n";

			getDLogLikelihood(gradM, aLengths, dMean, dVar+EPS);
			gradM -= grad;
			gradM /= EPS;

			cout << gradM << "=dJds\n";
		}*/

		TJMatrix H;
		getDerivs(H, aLengths, dMean, dVar);

		//cout << H << "=H\n";

		Eigen::SelfAdjointEigenSolver<TJMatrix> eig(H);
		//cout << eig.eigenvalues() << "=EVs\n";

		if(eig.eigenvalues().maxCoeff() > -EPS)
		{
			TJMatrix Lambda; Lambda.setZero();
			Lambda.diagonal() = eig.eigenvalues();

			//TJMatrix test = eig.eigenvectors() * Lambda * eig.eigenvectors().inverse();
			//cout << test << "check = H\n";

			for(int i=0;i<2;i++)
				if(Lambda(i,i) > -EPS)
					Lambda(i,i) = -EPS;

			TJMatrix Hmod = eig.eigenvectors() * Lambda * eig.eigenvectors().inverse();

			{
				Eigen::SelfAdjointEigenSolver<TJMatrix> eig2(Hmod);
				if(IS_DEBUG) CHECK(eig2.eigenvalues().minCoeff() > -EPS*0.99, "Failed to make neagitev-definite")
			}

			grad = Hmod.inverse() * grad;
		}
		else
			grad = H.inverse() * grad;

		double LAMBDA_ThisIter = LAMBDA;

		for(int iter2 = 0; LAMBDA > EPS; iter2++)
		{
			double dMean_temp = dMean - grad(0) * LAMBDA_ThisIter;
			double dVar_temp = dVar - grad(1) * LAMBDA_ThisIter;

			/*if(dVar_temp < EPS)
			{
				dVar_temp = EPS; //Keep positive
				if(LAMBDA_ThisIter < EPS && fabs(grad(0) < EPS))
				{
					cout << "Variance going negative\n";
					return false;
				}
			}*/

			{
				double LL_new = getLogLikelihood(aLengths, dMean_temp, dVar_temp); //Probably can do without evaluating...
				if(LL <= LL_new)
				{
					dMean = dMean_temp;
					dVar = dVar_temp;
					improvement = LL_new - LL;
					LL=LL_new;
					break;
				}
			}
			LAMBDA_ThisIter *= DROP;
		}
		bConverged = grad.squaredNorm() < EPS || improvement < EPS;
	}
	if(!bConverged)
	{
		cout << "Failed mean: " << dMean << endl;
		cout << "Failed var: " << dVar << endl;
		cout << "Convergence failed; improvement: " << improvement << "\nGradient " << grad.transpose() << endl;
	}

	return bConverged;
}

bool CScaleMLE::getMLEScaleParams(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose)
{
	dMean = 0;
	dVar = 1;

	double dRange = 0;

	getMeanVar(aLengths, dMean, dVar, dRange, bVerbose);

	double dVarUB = sqr(dRange); //Detect when dVar is going to infinity;

	if(bVerbose)
		cout << "Mean: " << dMean << " var: " << dVar << " range: " << dRange << endl;

	double LAMBDA = 0.2;
	const double DROP=0.2, BOOST=5;

	/*{
		TResidVector resid_old;
		for(double dV=0.1; dV<2; dV+=0.2)
			for(double dM=-4; dM<4; dM+=0.25)
			{
				getDLogLikelihood(resid_old, aLengths, dM, dV);
				cout << "Mean: " << dM << " var: " << dV;
				cout << "Residual: " << resid_old.transpose() << " = " << resid_old.squaredNorm() << endl;
			}
		cout << "#################################################################\n";
		for(double dV=0.01; dV<6; dV+=0.04)
		{
			getDLogLikelihood(resid_old, aLengths, 2, dV);
			cout << "Mean: " << 2 << " var: " << dV;
			cout << "Residual: " << resid_old.transpose() << " = " << resid_old.squaredNorm() << endl;
		}
	}*/

	TResidVector resid_old;
	getDLogLikelihood(resid_old, aLengths, dMean, dVar);
	double dOldSquaredErr = resid_old.squaredNorm();
	const double EPS = 1e-4;
	const int MAX_ITERS = 200;

	for(int iter = 0; iter < MAX_ITERS 
		&& dOldSquaredErr > sqr(EPS) 
		&& dVar < dVarUB; 
		iter++)
	{
		TJMatrix J;
		getDerivs(J, aLengths, dMean, dVar);

		if(bVerbose)
		{
			cout << "Residual: " << resid_old.transpose() << endl;
			cout << "Mean: " << dMean << " var: " << dVar << endl;
		}

		for(int iter2 = 0; ; iter2++)
		{
			TJTJMatrix JTJ = J.transpose() * J;
			JTJ.diagonal().array() += LAMBDA;

			TParamVector delta = JTJ.inverse() * J.transpose() * resid_old;

			const double dMean_temp = dMean - delta(0);
			const double dVar_temp = dVar - delta(1);

			if(dVar_temp > 0) 
			{
				TResidVector resid_new;
				getDLogLikelihood(resid_new, aLengths, dMean_temp, dVar_temp);
				double dNewSquaredErr = resid_new.squaredNorm();
				if(dNewSquaredErr < dOldSquaredErr)
				{
					dMean = dMean_temp;
					dVar = dVar_temp;
					resid_old = resid_new;
					dOldSquaredErr = dNewSquaredErr;
					LAMBDA *= DROP;
					break;
				}
			}
			LAMBDA *= BOOST;
		}

		if(dVar <= 0)//Otherwise might want to RESET, going to -inf
		{
			cout << "Error, Variance has gone negative";
			return false;
		}
	}
	bool bSuccess = dOldSquaredErr < sqr(EPS);

	if(bVerbose || !bSuccess)
	{
		cout << (bSuccess ? "Success" : "Fail") << endl;
		cout << "Mean: " << dMean << endl;
		cout << "Var: " << dVar << endl;
	}

	return bSuccess;
}

void CScaleMLE::testScaleMLE()
{
	return;

	CBoWSpeedo::CBoWObject::TLengths aLengths;
	for(double dM = 0; dM < 7; dM++)
	{
		CBoWSpeedo::CBoWObject::CLengthAndWeight LW(dM, CScale(dM, 0.1*(sqr(dM)+1)), -1);
		aLengths.push_back(LW);

		double dMean, dVar;
		//getMLEScaleParams(aLengths, dMean, dVar, true);
		getMLEScaleParams_2ndOrderGD(aLengths, dMean, dVar, true);
		cout << "Mean: " << dMean << endl;
		cout << "Var: " << dVar << endl;

		double dLL = getLogLikelihood(aLengths, dMean, dVar);
		double dLL1 = getLogLikelihood(aLengths, dMean+0.1, dVar);
		double dLL2 = getLogLikelihood(aLengths, dMean-0.1, dVar);
		double dLL3 = getLogLikelihood(aLengths, dMean, dVar+0.1);
		double dLL4 = (dVar-0.1>0) ? getLogLikelihood(aLengths, dMean, dVar-0.1) : -HUGE;

		CHECK(dLL<dLL1,"Gradient-descent failed")
		CHECK(dLL<dLL2,"Gradient-descent failed")
		CHECK(dLL<dLL3,"Gradient-descent failed")
		CHECK(dLL<dLL4,"Gradient-descent failed")
	}
}

/*class NGParams
{
	double lambda, gamma, alpha, beta;
public:
	NGParams(double lambda, double gamma, double alpha, double beta)
	 : lambda(lambda), gamma(gamma), alpha(alpha), beta(beta)
	{}

	void observe(double mu, double tau)
	{
		lambda = (lambda*gamma + mu) / (gamma+1);
		gamma++;
		alpha += 0.5;
		beta += gamma*sqr(mu-lambda) / (2 * (gamma+1));
	}
};

//See http://en.wikipedia.org/wiki/Normal-gamma_distribution and http://en.wikipedia.org/wiki/Conjugate_prior
bool CScaleMLE::getBayesianScaleParams_NG(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose)
{
	//Make prior
	NGParams ngParams();

	//Observe 1 measurement
	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pMeasurement = aLengths.begin(); pMeasurement != aLengths.end(); pMeasurement++)
	{
		const double mu = pMeasurement->getLength();
		const double tau = 1.0/pMeasurement->getVar();

		ngParams.observe(mu, tau);
	}

	dMean = ngParams.mean();
	dVar = ngParams.var();
}*/

//See http://en.wikipedia.org/wiki/Normal-gamma_distribution and http://en.wikipedia.org/wiki/Conjugate_prior
//Seperated, rev1002
bool CScaleMLE::getBayesianScaleParams_NG(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose)
{
	//Make prior: params for Normal distn for mean
	double mu0 = 0, tau_forMean0 = 0; //Improper--could be eps

	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pMeasurement = aLengths.begin(); pMeasurement != aLengths.end(); pMeasurement++)
	{
		if(pMeasurement->hasSLAMLength())
		{
			//Observe 1 measurement
			const double mu = pMeasurement->getSLAMLength();
			const double tau = pMeasurement->getPrecision();

			//Weerahandi-1995 and Wikipedia are sources for calculations
			mu0 = (tau_forMean0*mu0 + tau*mu) / (tau_forMean0 + tau);
			tau_forMean0 = tau_forMean0 + tau; //Strictly converging????? This is the confidence in the MEAN

			if(bVerbose)
			{
				cout << "Observation " << mu << ", var " << 1.0/tau << endl;
				cout << "Expected measurement " << mu0 << ", var " << 1.0/tau_forMean0 << endl;
				if(IS_DEBUG) CHECK(isnan(mu0) || isnan(tau_forMean0) , "NAN calculating distn params")
			}
		}
	}

	//Gamma prior for tau
	//alpha=k, beta = 1/theta for normal parameterisation (k, theta)
	//Mean alpha/beta and var alpha/beta^2.
	double alpha = 1, beta=0.001;
	//Could be 0,0 also...

	for(CBoWSpeedo::CBoWObject::TLengths::const_iterator pMeasurement = aLengths.begin(); pMeasurement != aLengths.end(); pMeasurement++)
	{
		if(pMeasurement->hasSLAMLength())
		{
			//Observe 1 measurement
			const double mu = pMeasurement->getSLAMLength();
			const double tau = pMeasurement->getPrecision();

			//Seperate this because mu parameterisation independent of these
			alpha = alpha + 0.5;
			beta = beta + sqr(mu0 - mu)*0.5;

			if(bVerbose)
			{
				cout << "Observation " << mu << ", var " << 1.0/tau << endl;
				cout << "Expected measurement var " << alpha/beta << ", var " << alpha/sqr(beta) << endl;
				if(IS_DEBUG) CHECK(isnan(mu0) || isnan(tau_forMean0) || isnan(alpha) || isnan(beta), "NAN calculating distn params")
			}
		}
	}
	dMean = mu0;
	dVar = alpha/beta;

	return true;
}
