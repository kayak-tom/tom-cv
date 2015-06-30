/*
 * fullRelPose.cpp
 *
 *  Created on: 14 Jan 2010
 *      Author: tom
 */

#include "fullRelPose.h"
#include <ostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;

double CScale::s_MaxBadness = 100.01;

void CFullRelPose::writeEdge(std::ostream & toroGraphFile, const bool b2d) const
{
	if(IS_DEBUG) CHECK(!hasPosition(), "Node is disconnected?? Or too bad");

	const double dLength = length();
	//double dVar = sqr(dLength/4);// scale.scaleVarHacked();
	double dVar = variance();
	dVar = std::max<double>(1e-5, dVar);

	normalisedPose.scale(dLength).write(toroGraphFile, b2d);

	Vector3d t, x_axis; x_axis.setZero();
	normalisedPose.t.asVector(t);

	if(b2d)
	{
		//double inf_ff, inf_fs, inf_ss, inf_rr, inf_fr, inf_sr;
		Vector2d t2d(t(0), t(2));
		Vector2d t2d_perp(t(2), -t(0));

		Matrix2d Q, LAMBDA;
		LAMBDA.setZero();

		Q << t2d, t2d_perp;

		LAMBDA(0,0) = dVar; //dVar may be zero
		LAMBDA(1,1) = dVar * sqr(normalisedPose.SD.cameraMotionAngleSD());

		CHECK(std::isnan(LAMBDA.sum()), "writeEdge: nan")

		const Matrix2d C = Q * LAMBDA * Q.transpose();
		const Matrix2d & I = C.inverse();
		//cout << I << endl;

		//THROW("Not complete...")
		double dRotInf = 1.0/sqr(normalisedPose.SD.relOrientationSD());
		toroGraphFile << I(0,0) << ' ' << I(1,0) << ' ' << I(1,1) << ' ' << dRotInf << " 0 0";

		return;
	}

	if(t(0) < 0.9) //check t isn't x axis aligned (x is arbitrary)
		x_axis(0) = 1;
	else
		x_axis(1) = 1;

	Matrix3d Q;

	if(t.sum() == 0) //Pure rotation. Should be arbitrary
	{
		Q.setIdentity();
	}
	else
	{
		Vector3d tperp1 = t.cross(x_axis);
		tperp1.normalize();

		Vector3d tperp2 = tperp1.cross(t);

		Q << t, tperp1, tperp2;
	}

	CHECK(std::isnan(Q.sum()), "writeEdge: nan")

	Matrix3d LAMBDA; LAMBDA.setZero();
	LAMBDA(0,0) = dVar; //dVar may be zero
	LAMBDA(1,1) = LAMBDA(2,2) = dVar * sqr(normalisedPose.SD.cameraMotionAngleSD());

	CHECK(std::isnan(LAMBDA.sum()), "writeEdge: nan")

	const Matrix3d C = Q * LAMBDA * Q.transpose();

	//std::cout << "Covariance: " << std::endl << C << std::endl;
	CHECK(std::isnan(C.sum()), "writeEdge: nan")

	Matrix3d Crpy; Crpy.setZero();
	Crpy.diagonal().setConstant(sqr(normalisedPose.SD.relOrientationSD()));
	Matrix<double, 6, 6> Cfull; Cfull << C, Matrix3d::Zero(), Matrix3d::Zero(), Crpy;

	CHECK(std::isnan(Cfull.sum()), "writeEdge: nan")

	if(Cfull.trace() < 0.0001)
	{
		std::cout << "Warning: Covariance matrix near-singular, adjusting diagonal\n";
		Cfull.diagonal().array() += 0.0001;
	}

	const Matrix<double, 6, 6> & Ifull = Cfull.inverse();
	//std::cout << "Information: " << std::endl << Ifull << std::endl;

	//cout << "Warning: Not inverting information matrix\n"; Yes the information mat does work a little better.
	//Possibly ok for 1 edge...if(IS_DEBUG) CHECK(Cfull.determinant()<0.0001, "CFullRelPose::writeEdge: Singular covariance matrix");

	CHECK(std::isnan(Ifull.sum()), "writeEdge: nan")

	for (int nRow = 0; nRow <6; nRow++)
	{
		for(int nCol = nRow; nCol < 6; nCol++)
			toroGraphFile << Ifull(nRow, nCol) << ' ';
	}
}

std::ostream& operator<<(std::ostream& s, const CScale & scale)
{
	if(scale.hasScale())
		s << "D=" << scale.D << ", G^2=" << scale.G_sq << /*", mean=" << scale.scaleMean() <<*/  ", median=" << scale.scaleMedian() << ", var=" << scale.scaleVar();
	else
		s << "[Scale not initialised]";
	return s;
}

std::ostream& operator<<(std::ostream& s, const CRelScale & relscale)
{
	if(relscale.notTooBad())
	{
		CScale scale; scale.setOrigin(true); scale = scale + relscale;
		s << "Relative scale: " << scale;
	}
	else
		s << "[Relative scale not initialised]";

	return s;
}

std::ostream& operator<<(std::ostream& s, const CFullRelPose & X)
{
	s << "Position: " << X.position(false) << std::endl;
	s << "Scale: " << X.getScale() << std::endl;
	return s;
}
