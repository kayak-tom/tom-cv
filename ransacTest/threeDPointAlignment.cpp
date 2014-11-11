/*
 * threeDPointAlignment.cpp
 *
 *  Created on: 18 Jan 2010
 *      Author: tom
 */

#include <iostream>
#include "ransac/ransac.h"
#include "ransac/findBestMatchingDisjoint.h"
#include "ransac/essentialMat_Eigen.h"
#include "ransac/openCV8PtEssentialMatModified.h"
#include <boost/smart_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include "ransac/ransacParams.h"
#include "geom/geom.h"
#include "geom/geom_eigen.h"
#include "time/SpeedTest.h"
#include "util/stats.h"
#include <fstream>
#include "geom/pointAlignment.h"

using namespace Eigen;
using namespace std;
using namespace boost;

void testAlignment(const int nIters, const int nPoints, int dDepth, double dSD, double dScale, ostream & results)
{
	//Simulate some 3d structure, project into 3 cameras, given transformations between points compute scale factor
	ARRAY(C3dPoint, aWorldPoints, nPoints);
	ARRAY(C3dPoint, aReconWorldPoints1, nPoints);
	ARRAY(C3dPoint, aReconWorldPoints2, nPoints);

	CStats estimatedScale;
	CStats estimatedVar;

	for(int i=0; i<nIters; i++)
	{

		C3dPoint t1, t2;
		t1.addNoise(1), t2.addNoise(1);

		//double dScale = t1.length() / t2.length();
		C3dRotation R1, R2;

		CCamera aP[] = {CCamera(), R1 | t1, R2 | t2};

		ARRAY(C2dPoint, aImPoints1, nPoints);
		ARRAY(C2dPoint, aImPoints2, nPoints);
		ARRAY(C2dPoint, aImPoints3, nPoints);

		C2dPoint * aImPoints[3] = {PTR(aImPoints1), PTR(aImPoints2), PTR(aImPoints3)};

		for(int i=0;i<nPoints; i++)
		{
			aWorldPoints[i] = C3dPoint(CRandom::Uniform(-1., 1.), CRandom::Uniform(-1., 1.), CRandom::Uniform(-1. + dDepth, 1. + dDepth));
			//Make depth actually dDepth
			aWorldPoints[i] *= dDepth/aWorldPoints[i].length();
		}

		for(int nCam=0;nCam<3; nCam++)
			for(int i=0;i<nPoints; i++)
			{
				aImPoints[nCam][i] = aP[nCam] * aWorldPoints[i];
				aImPoints[nCam][i].addNoise(dSD);
			}

		C3dPointMatchVector matches;

		for(int i=0;i<nPoints; i++)
		{
			aReconWorldPoints1[i] = reconstruct(aP[0], aP[1], aImPoints[0][i], aImPoints[1][i]);
			//cout << aReconWorldPoints1[i].depth(aP[0]) << endl;
			aReconWorldPoints2[i] = dScale*reconstruct(aP[0], aP[2], aImPoints[0][i], aImPoints[2][i]);
			matches.push_back(C3dPointMatch(aReconWorldPoints1[i], aReconWorldPoints2[i]));
		}

		char fileName[1000], fileNameInliers[1000];
		sprintf_s(fileName, 1000, "Points=%d_Depth=%d_all.dat",  nPoints, dDepth);
		sprintf_s(fileNameInliers, 1000, "Points=%d_Depth=%d_inliers.dat",  nPoints, dDepth);

		//saveScales(matches, 0.05, fileName, fileNameInliers, true);

		double dScaleEst = -1, dVarEst = -1;
		getDG(matches, 0.05, dScaleEst, dVarEst, false);
		cout << dScaleEst << "== scale estimate grubbs, dVar=" << dVarEst << ", exp(mean)=median=" << exp(dScaleEst) << endl;

		estimatedScale.add(dScaleEst);
		estimatedVar.add(dVarEst);
	}
	estimatedScale.writeTSVdata(results);
	estimatedVar.writeTSVdata(results);
	results << endl;

	return;

	/*getTandR_1d_Median(matches, 0.05, dScaleEst);
	cout << dScaleEst << "== scale estimate median" << endl;

	getTandR_1d_Mean(matches, 0.05, dScaleEst);
	cout << dScaleEst << "== scale estimate mean" << endl;
	getTandR_1d_Ransac(matches, 0.05, 0.1, 1, dScaleEst, dVarEst);
	cout << dScaleEst << "== scale estimate 1d ransac, var=" << dVarEst << endl;
	getTandR_1d_MeanID(matches, 0.05, dScaleEst);
	cout << dScaleEst << "== scale estimate 1d mean ID" << endl;
	
	getTandR_1d_MeanLogs(matches, 0.05, dScaleEst, dVarEst);
	cout << dScaleEst << "== scale estimate 1d mean Logs, dVar=" << dVarEst << ", exp(mean)=median=" << exp(dScaleEst) << endl;
	getTandR_1d_RANSACLogs(matches, 0.05, 0.1, 1, dScaleEst, dVarEst);
	cout << dScaleEst << "== scale estimate 1d RANSAC Logs, dVar=" << dVarEst << ", exp(mean)=median=" << exp(dScaleEst) << endl;

	cout << dScale << "== true scale" << endl;*/
}

void testAlignment()
{
	ofstream results("1dAlignResults.tsv");

	results << "NumPoints\tSD\tScale\tDepth\t";
	CStats::writeTSVheader("Scale", results);
	CStats::writeTSVheader("Var", results);
	results << endl;

	double dSD = 0.01; //radians
	double dScale = 4;
	for(int nPoints = 5; nPoints < 20; nPoints++)
		for(int dDepth = 100; dDepth < 200; dDepth+=25)
		{
			results << nPoints << '\t' << dSD << '\t' << dScale << '\t' << dDepth << '\t';
			testAlignment(2000, nPoints, dDepth, dSD, dScale, results);
			return;
		}

	results.close();
}
