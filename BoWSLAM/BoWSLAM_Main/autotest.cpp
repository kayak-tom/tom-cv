/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "autotest.h"
#include "ransac/ransac.h"
#include "ransac/essentialMat_Eigen.h"
#include "ransac/inlierCounter.h"
#include "geom/geom.h"
#include "util/opencv.h"
#include "image/HSI.h"
#include <algorithm>
#include <queue>
#include "util/random.h"
#include "time/SpeedTest.h"
#include "util/smallHashTable.h"
#include "util/fastSet.h"
#include "util/set2.h"
#include "scale/EdgeScaleObserver.h"
#include "bow/scales.h"
#include "scale/ScaleMLE.h"

#include "geom/pointAlignment.h" //and expose a couple more fns for testing:

double getErr(const C3dPointMatch * pPM, const C3dRotation &R_hyp, const C3dPoint &T_hyp, double c_hyp);
bool getTandR_1d(const C3dPointMatch * pMatch, const C3dRotation &R0a, const C3dPoint &T_0a_0, const C3dPoint &T_ab_0_dir, double &c);

#include "util/simpleBinaryTree.h"

using namespace std;

/*bool fpeq(const Matrix &M1, const Matrix &M2)
{
	return zero((M1-M2).sum_square());
}*/
template<typename FLOAT, bool DIV>
void timeDiv()
{
	CStopWatch s;
	s.startTimer();
	FLOAT d = 0;
	for(int i=1; i<200000000; i++)
		d += DIV ? FLOAT(1.0)/FLOAT(i) : FLOAT(0.123)*FLOAT(i);
	s.stopTimer();
	cout << sizeof(FLOAT) << " result " << d << " time " << s.getElapsedTime() << endl;
}

template<typename T>
void timeAlignmentCopy()
{
	T * pX = new T[1000000/sizeof(T)], *pY = pX;//+(CRandom::Uniform(0.1,1.73) == 0.5 ? 3 : 29);

	CStopWatch s;
	s.startTimer();

	int sum = 0;
	for(int i=0; i<100000; i++, pY += 4/sizeof(T))
		sum += *(reinterpret_cast<int *>(pY));
	s.stopTimer();

	delete [] pX;

	cout << "Size: " << sizeof(T) << " sum: " << sum << " time: " << s.getElapsedTime() << endl;
}

double E_sumSquaredResiduals(const Eigen::Matrix3d & E, const T2dPoints & m0in, const T2dPoints & m1in, const CMask & acMask)
{
	//THROW("This measure is a bit meaningless, reprojection error would be better") OK for debug tho
	//Matrix E = cvMatToMatrix(cvF);
	const TPointVec2d & m0 = m0in;
	const TPointVec2d & m1 = m1in;
	double dSum = 0;

	for(int i=0;i<acMask.size(); i++)
	{
		if(acMask[i])
		{
			Eigen::Vector3d v0, v1;
			m0[i].asVector(v0);
			m1[i].asVector(v1);
			dSum += (v1.transpose() * E * v0).squaredNorm();
		}
	}

	return dSum;
}

double E_sumSquaredResiduals(const double * pdE, const T2dPoints & m0in, const T2dPoints & m1in, const CMask & acMask)
{
	Eigen::Matrix3d E(pdE); E.transposeInPlace();
	return E_sumSquaredResiduals(E, m0in, m1in, acMask);
}


bool fpeq(double a, double b)
{
	return zero(a-b);
}

bool fpne(double a, double b)
{
	return !zero(a-b);
}

/*bool isId(const Matrix &M)
{
	IdentityMatrix I(M.ncols());
	return zero((M-I).sum_square());
}*/

bool coloursEq(int r1, int g1, int b1, int r2, int g2, int b2)
{
	return (abs(r1-r2) + abs(g1-g2) + abs(b1-b2)) < 15;
}



#define ATEST3(t, l) CHECK(!(t), "AUTOTEST FAILED: " __FILE__ # l " \"" #t "\"");
#define ATEST2(t, l) ATEST3(t, l);
#define ATEST(t) ATEST2(t, __LINE__);

//Test HSI/HSV (?) conversions are invertable
void testColConvertHSI(int h1, int s1, int i1)
{
	int r, g, b, h, s, i, r2, g2, b2;
	CHsiInt::RGB(h1, s1, i1, r, g, b);
	CHsiInt::HSI(r, g, b, h, s, i);
	CHsiInt::RGB(h, s, i, r2, g2, b2);
	ATEST(coloursEq(r, g, b, r2, g2, b2)) //HSI space is redundant, but we only ever get HSI colours from the HSI function so this is ok
	//ATEST(coloursEq(h, s, i, h1, s1, i1))
}

void testColConvertRGB(int r1, int g1, int b1)
{
	int r, g, b, h, s, i;
	CHsiInt::HSI(r1, g1, b1, h, s, i);
	CHsiInt::RGB(h, s, i, r, g, b);
	ATEST(coloursEq(r, g, b, r1, g1, b1))
}

void testVectors();
void testFindBestMatching();

#define MATRIX(M, r, c) double M ## _data[r*c]; LaGenMatDouble M(M ## _data, r, c, false);
#define VECTOR(M, r) double M ## _data[r]; LaVectorDouble M(M ## _data, r);
#define ROWVECTOR(M, r) double M ## _data[r]; LaRowVectorDouble M(M ## _data, r);

//Check SVs are 2 equal and 0
bool checkSVs(const C3x3MatModel & E)
{
	CCamera aPp[4];
	getCamsFromE(E.asDouble9(), aPp); //Does the SVD check...
	return true;
	//This is one from MATLAB--check we also find it has dodgy SVs (no it doesn't):
	/*double pdE2[9] =	   {-0.5301  ,	   -0.1615   ,	    0.2346   ,	   -0.0983  ,	    0.5385  ,	   -0.3300 ,	    0.4355 ,	   -0.1978 ,	    0.0244 };*/
/*
	const LaGenMatDouble E((double*)(void*)pdE.asDouble9(), 3, 3, false);
	LaGenMatDouble Ecopy = E;
//	cout << "E:" << E << "\n";

	LaVectorDouble S(3);
	LaGenMatDouble U_dontCalc(0, 0); //Todo: don't calc (try taking transpose of everything and not calc'ing VT)
	LaGenMatDouble V_trans(3, 3);

	LaSVD_IP(Ecopy, S, U_dontCalc, V_trans);

 //cout << "\nU=\n" << U_dontCalc;
 //cout << "\nS=\n" << S;
 //cout << "\nV=\n" << V_trans;
 //cout << "\n";

	return zero(S(2)) && zero(S(0)-S(1));*/
}

double checkE(const C3x3MatModel & pdE, const T2dPoints & m0, const T2dPoints & m1, const CMask & acMask)
{
	cout << "Todo: E check\n";
	return 0;
/*
	Eigen::Matrix E((double*)(void*)pdE.asDouble9(), 3, 3, false);
	int nPoints = m0.size();

	ROWVECTOR(p1, 3);
	VECTOR(p2, 3);
	MATRIX(temp, 1, 3)
	MATRIX(resid, 1, 1);

	double dResid=0, dBadResid = 0;

	const bool BAD_RESID = false;

	//Copy points into vectors
	if(acMask.countInliers())
	{
		for(int i=0; i<nPoints; i++)
		{
			if(acMask[i])
			{
				p1(0) = m0[i].getX();
				p1(1) = m0[i].getY();
				p1(2) = 1;
				p2(0) = m1[i].getX();
				p2(1) = m1[i].getY();
				p2(2) = 1;

				//LaGenMatDouble resid = p1*E*p2;
				Blas_Mat_Mat_Mult(p1,E, temp);
				Blas_Mat_Mat_Mult(temp, p2, resid);
				//cout << "Residual = " << resid << "\n";
				dResid += *resid_data;

				if(BAD_RESID)
				{
					p1(0) = m1[i].getX();
					p1(1) = m1[i].getY();
					p2(0) = m0[i].getX();
					p2(1) = m0[i].getY();
					//LaGenMatDouble resid = p1*E*p2;
					Blas_Mat_Mat_Mult(p1,E, temp);
					Blas_Mat_Mat_Mult(temp, p2, resid);
					//cout << "Bad residual = " << resid << "\n";
					dBadResid += *resid_data;
				}
			}
		}
	}
	else
	{
		//Copy points into vectors
		for(int i=0; i<nPoints; i++)
		{
			p1(0) = m0[i].getX();
			p1(1) = m0[i].getY();
			p1(2) = 1;
			p2(0) = m1[i].getX();
			p2(1) = m1[i].getY();
			p2(2) = 1;

			//LaGenMatDouble resid = p1*E*p2;
			Blas_Mat_Mat_Mult(p1,E, temp);
			Blas_Mat_Mat_Mult(temp, p2, resid);
			//cout << "Residual = " << resid << "\n";
			dResid += *resid_data;

			if(BAD_RESID)
			{

				p1(0) = m1[i].getX();
				p1(1) = m1[i].getY();
				p2(0) = m0[i].getX();
				p2(1) = m0[i].getY();
				//LaGenMatDouble resid = p1*E*p2;
				Blas_Mat_Mat_Mult(p1,E, temp);
				Blas_Mat_Mat_Mult(temp, p2, resid);
				//cout << "Bad residual = " << resid << "\n";
				dBadResid += *resid_data;
			}
		}
	}
	//cout << "Total residual = " << dResid << ", bad = " << dBadResid << "\n";

	if(BAD_RESID)
	{
		return fabs(dResid/dBadResid);
	}
	return fabs(dResid);*/
}

//Autotests:

/*void testEigenSpeed()
{
    double test2_Q1x[5] = {14  ,  79   ,  4   , 68 ,   39};
    double test2_Q1y[5] = {42  ,  96  ,  85   , 76    ,66};

    double test2_Q2x[5] = {71  ,   5   , 69  ,   3 ,   77};
    double test2_Q2y[5] = {3  ,  10  ,  32  ,  44  ,  80};

    CvPoint2D64f m0[5];
    CvPoint2D64f m1[5];

    int nPoints=5;
    for(int i=0;i<nPoints;i++)
    {
    	m0[i].x = test2_Q1x[i];
    	m0[i].y = test2_Q1y[i];
    	m1[i].x = test2_Q2x[i];
    	m1[i].y = test2_Q2y[i];
    }
    double pdEssentialMat2[10*9];

    int x=0;
    CStopWatch s;
    s.startTimer();
    for(m1[0].y=-50000; m1[0].y<50000; m1[0].y++)
    	x += calcEssentialMat_5point_Eigen( m0, m1, pdEssentialMat2, 0,0, false, true);
    s.stopTimer();
    cout << "Time (100000*)=" << s.getElapsedTime() << endl;
}*/

template<typename T>
class testStdSet : public std::set<T, std::less<T>, std::allocator<T> >
{
public:
	bool exists(const T & t) const { return std::set<T>::find(t) != std::set<T>::end(); }
};

#define TFS(TSet) testFastSet<TSet>(false /*Set back to true for better test*/, #TSet)

template<typename TSetType> void testFastSet(const bool bMillion = false, const char * szName = 0)
{
	CStopWatch s;
	s.startTimer();
	{
		TSetType fastSet;
		fastSet.insert(1);
		fastSet.insert(1);
		ATEST(fastSet.exists(1))
		ATEST(!fastSet.exists(2))
		ATEST(fastSet.size() == 1)
		fastSet.erase(1);
		ATEST(fastSet.size() == 0)
		fastSet.insert(1);
		fastSet.insert(2);
		fastSet.insert(3);
		fastSet.erase(2);
		ATEST(fastSet.size() == 2)
		fastSet.insert(1);
		fastSet.insert(2);
		fastSet.insert(3);
		fastSet.erase(2);
		ATEST(fastSet.size() == 2)
		fastSet.insert(4);
		fastSet.insert(5);
		fastSet.insert(6);
		ATEST(fastSet.size() == 5)
		ATEST(fastSet.exists(1))
		ATEST(!fastSet.exists(2))

		fastSet.erase(2);
		ATEST(fastSet.size() == 5)
		for(int i=2;i<=100;i++)
		{
			fastSet.insert(i-1);
			fastSet.insert(i);
			ATEST(fastSet.exists(i))
		}
		ATEST(fastSet.size() == 100)

		int nCount = 0;
		for(typename TSetType::const_iterator p = fastSet.begin(); p != fastSet.end(); ++p)
			nCount++;

		ATEST(nCount == 100);

		for(int i=1;i<=100;i++)
		{
			fastSet.erase(i);
			fastSet.insert(i+1);
		}
		ATEST(fastSet.size() == 1)

		ATEST(fastSet.exists(101))

		nCount = 0;
		for(typename TSetType::const_iterator p = fastSet.begin(); p != fastSet.end(); ++p)
		{
			nCount++;
			ATEST(*p == 101);
		}
		ATEST(nCount == 1);

		ATEST(fastSet.exists(101));

		TSetType fastSet2;

		int nMax = bMillion ? 100000 : 1000;

		for(int i=0; i< nMax; i++)
		{
			if(i%2 == 0)
			{
				fastSet2.insert(CRandom::Uniform(0, i));
			}
			else
			{
				fastSet2.erase(CRandom::Uniform(0, 1000000));
				if(bMillion && CRandom::Uniform(0, 10000) == 0)
				{
					int nLastElement = -1;
					for(typename TSetType::const_iterator p = fastSet2.begin(); p != fastSet2.end(); ++p)
					{
						int element = *p;
						ATEST(nLastElement < element);
					}
				}
			}
		}
		int nLastElement = -1;
		CDynArray<int> arr;
		for(typename TSetType::const_iterator p = fastSet2.begin(); p != fastSet2.end(); ++p)
		{
			int element = *p;
			if(bMillion)
				ATEST(nLastElement < element);
			//fastSet2.erase(element); good way to break it!!
			arr.push_back(element);
		}
		ATEST( (int)arr.size() == (int)fastSet2.size() );

		for(CDynArray<int>::iterator p = arr.begin(); p < arr.end(); p++)
			fastSet2.erase(*p);

		ATEST(fastSet2.size() == 0);
	}
	s.stopTimer(); 

	if (szName && bMillion)
	{
		cout << szName;
		cout << " took " << s.getElapsedTime() << endl;
	}
}

template<typename TSetType> void testFastSetLocations()
{
	TSetType locFastSet;
	for(int i=1;i<=100;i++)
	{
		locFastSet.insert(CLocation(i,i));
		ATEST(locFastSet.size() == i)
		ATEST(locFastSet.exists(CLocation(i,i)));

		locFastSet.insert(CLocation(i,i));
		ATEST(locFastSet.size() == i)
		ATEST(locFastSet.exists(CLocation(i,i)));
	}
	for(int i=1;i<=100;i++)
	{
		locFastSet.erase(CLocation(i,i));
		ATEST(locFastSet.size() == 100-i)
		ATEST(!locFastSet.exists(CLocation(i,i)));
	}
}

void testLognormal()
{
	///Lognormal distn tests:
	CScale initScale, uninitScale;
	ATEST(!(initScale.isMoreAccurateThan(uninitScale)));
	ATEST(!initScale.hasScale());
	initScale.setOrigin(true);
	//cout << initScale << endl;
	ATEST(initScale.hasScale());
	CRelScale uninformativeScale(CRelScale::eUninformativeScale, 1);
	ATEST(uninformativeScale.notTooBad())
	CScale scale1 = initScale + uninformativeScale;
	//cout << scale1 << endl;
	ATEST(scale1.hasScale());
	ATEST(scale1.scaleEstimate() == 1);
	CScale scale2 = scale1 + uninformativeScale;
	//cout << scale2 << endl;
	ATEST(scale2.hasScale());
	ATEST(fpeq(scale2.scaleEstimate(), 1));

	ATEST(scale2.isMoreAccurateThan(uninitScale));
	ATEST(scale1.isMoreAccurateThan(scale2));
	//ATEST(scale2 > scale1);
	//ATEST(!(scale2 > scale2));
	ATEST(!(scale2.isMoreAccurateThan(scale2)));
	ATEST(!(scale1.isMoreAccurateThan(scale1)));

	const double dScale = 2, dVar=1;
	double d, g;
//	CLogNormal::wikipediaParamEstimation(dScale, dVar, d, g);
	CLogNormal::medianParamEstimation(dScale, dVar, d, g);

	CRelScale scaleDoubling(d, g);
	//cout << scaleDoubling << "scaleDoubling" << endl;
	CScale twoTimes = initScale + scaleDoubling;
	//cout << twoTimes << "= twoTimes" << endl;
	ATEST(fpeq(twoTimes.scaleEstimate(), 2));
	CScale fourTimes = twoTimes + scaleDoubling;
	//cout << fourTimes << "= fourTimes" << endl;
	ATEST(fpeq(fourTimes.scaleEstimate(), 4));
	CScale fourTimesBigVar = fourTimes + uninformativeScale;
	//cout << fourTimesBigVar << "= fourTimesBigVar" << endl;
	ATEST(fpeq(fourTimesBigVar.scaleEstimate(), 4));
	CRelScale scaleHalving = scaleDoubling.inverse();
	//cout << scaleHalving << "scaleHalving" << endl;
	CScale twoTimesBiggerVar = fourTimesBigVar + scaleHalving;
	//cout << twoTimesBiggerVar << "= twoTimesBiggerVar" << endl;
	ATEST(fpeq(twoTimesBiggerVar.scaleEstimate(), 2));
	CScale halfTimesBiggerVar = twoTimesBiggerVar + scaleHalving + scaleHalving;
	//cout << halfTimesBiggerVar << "= halfTimesBiggerVar" << endl;
	ATEST(fpeq(halfTimesBiggerVar.scaleEstimate(), 0.5));
	ATEST(twoTimesBiggerVar.isMoreAccurateThan(halfTimesBiggerVar))
	ATEST(twoTimes.isMoreAccurateThan(twoTimesBiggerVar) )
	ATEST(halfTimesBiggerVar.isMoreAccurateThan(uninitScale))
}

bool test()
{
	cout << "Running Autotests...";

	//testTiming();

    try
    {
    	/*timeDiv<float, true>();
    	timeDiv<double, true>();
    	timeDiv<float, true>();
    	timeDiv<double, true>();
    	timeDiv<float, false>();
    	timeDiv<double, false>();
    	timeDiv<float, false>();
    	timeDiv<double, false>();

		timeAlignmentCopy<char>();
		timeAlignmentCopy<int>();
		timeAlignmentCopy<char>();
		timeAlignmentCopy<int>();
    	THROW("quit");*/

    	float f_pos=123.4, f_neg=-123.4;
    	ATEST(doubleToInt(f_pos) == (int) f_pos)
    	ATEST(doubleToInt(f_neg) == (int) f_neg)

    	double d_pos=123.4, d_neg=-123.4;
		ATEST(doubleToInt(d_pos) == (int) d_pos)
		ATEST(doubleToInt(d_neg) == (int) d_neg)

    	srand(-100);
#if EXACT_HACK > 1 //If we've enabled this hack, need to be using synth data
    	ATEST(bSyntheticData)
#endif

    	testVectors();

    	C3dPoint axis(1,2,3);
		C3dRotation R(axis, 0.3), I;

		ATEST(fpeq(R.angle(), 0.3))
		ATEST(fpeq(R.t().angle(), 0.3))
		ATEST(fpeq((R*R).angle(), 0.6))
		ATEST(fpeq((R*R.t()).angle(), 0.0))
		ATEST(fpeq((R.t()*R.t()).angle(), 0.6))
		ATEST(fpeq(I.angle(), 0.0))
		ATEST(fpeq(I.t().angle(), 0.0))
		ATEST(fpeq((R*I).angle(), 0.3))
		ATEST(fpeq((I*R).angle(), 0.3))

/*    	C3dRotation Ry=getRotMat(0, 0.3, 0);
    	C3dRotation Rz=getRotMat(0, 0, 0.3);

		ATEST(fpeq(quatAngle(rotMatToQuat(I3)), 0));

		ATEST(fpeq(Rx.angle(), 0.3));
		ATEST(fpeq(Ry.angle(), 0.3));
		ATEST(fpeq(Rz.angle(), 0.3));

		cout << Rx << endl;
		ColumnVector c_testConversion(3); c_testConversion << 0 << 0 << 1;
		C3dPoint p_testConversion(0,0,1);
		PRINTMAT(Rx.asMat())
		PRINTMAT(Rx.asMat()*c_testConversion)
		PRINTMAT(Rx * p_testConversion)

		PRINTMAT(getRotMat(0.3, 0, 0))
		PRINTMAT(getRotMat(0.3, 0, 0)*c_testConversion)

		ATEST(fpeq(Rx.asMat(), getRotMat(0.3, 0, 0)))

		C3dRotation Rpi=getRotMat(3.14159, 0, 0);
		ATEST(fpeq(Rpi.angle() , 3.14159));

		Matrix Rmat = getRotMat(3.14159, 0.4, 0.7);
		C3dRotation R=Rmat;
		ATEST(isRotMat(R.asMat()));
		ATEST(isRotMat(Rz.asMat()));
		ATEST(isRotMat(Ry.asMat()));
		ATEST(isRotMat(Rx.asMat()));
		ATEST(isRotMat(Rpi.asMat()));
		ATEST(!isRotMat(-Rpi.asMat()));

		ATEST(fpeq(R.asMat(), Rmat))*/

		//Now test 3d alignment:
		const double dActualScale = 17.3;
		//C3dPoint actualT(123, -3, 0.04);
		//C3dRotation & actualR = R;
		C3dPointMatchVector vPointMatches;
		C3dPointMatchVector vPointMatchesRev;

		for(double x=13; x<15;x++)
			for(double y=-10; y<20; y += 10)
				for(double z=704; z<705;z+=0.3)
				{
					C3dPoint pPoint1(x,y,z);
					C3dPoint pPoint2 = pPoint1 * dActualScale;
					//PRINTMAT(*pPoint1)
					//PRINTMAT(*pPoint2)
					vPointMatches.push_back(C3dPointMatch(pPoint1, pPoint2));
					vPointMatchesRev.push_back(C3dPointMatch(pPoint2, pPoint1));
				}

		double dCalcScale = -1, dCalcVar = -1;
		ATEST(vPointMatches.size() == getTandR_1d_Mean(vPointMatches, 0.1, dCalcScale));
		ATEST(fpeq(dCalcScale, 1.0/dActualScale) )

		ATEST(vPointMatches.size() == getTandR_1d_Median(vPointMatches, 0.1, dCalcScale));
		ATEST(fpeq(dCalcScale, 1.0/dActualScale) )

		ATEST(vPointMatches.size() == getTandR_1d_Ransac(vPointMatches, 0.1, 0.1, 5, dCalcScale, dCalcVar, false));
		ATEST(fpeq(dCalcScale, 1.0/dActualScale) )

		ATEST(vPointMatches.size() == getTandR_1d_Mean(vPointMatchesRev, 0.1, dCalcScale));
		ATEST(fpeq(dCalcScale, dActualScale) )

		ATEST(vPointMatches.size() == getTandR_1d_Median(vPointMatchesRev, 0.1, dCalcScale));
		ATEST(fpeq(dCalcScale, dActualScale) )

		ATEST(vPointMatches.size() == getTandR_1d_Ransac(vPointMatchesRev, 0.1, 0.1, 5, dCalcScale, dCalcVar, false));
		ATEST(fpeq(dCalcScale, dActualScale) )
		//double dScale; //The scale is relative to a 3d reconstruction assuming a baseline of Tdir => Tactual = Tdir/dScale

		C3dPoint T_from3d;
		C3dRotation R_from3d;
		/*ATEST(getTandR(vPointMatches, R_from3d, T_from3d, dScale));

		ATEST((R_from3d == actualR));
		ATEST((T_from3d == actualT));
		ATEST(fpeq(dScale*dActualScale, 1));
		//ATEST(isRotMat(R_from3d));

		getTandR_RANSAC(vPointMatches, R_from3d, T_from3d, dScale);

		ATEST((R_from3d == actualR));
		ATEST((T_from3d == actualT));
		ATEST(fpeq(dScale*dActualScale, 1));
		//ATEST(isRotMat(R_from3d));
///////Now test the reverse...
		ATEST(getTandR(vPointMatchesRev, R_from3d, T_from3d, dScale));

		ATEST((R_from3d == actualR.t()));
//		PRINTMAT(T_from3d)
//		PRINTMAT(actualT)
//		PRINTMAT(R_from3d*T_from3d)
//		PRINTMAT((-1./dScale)*actualR*T_from3d);
//		cout << (T_from3d- (-1./dScale)*actualR*T_from3d).sum_square();
		ATEST((actualT == actualR*T_from3d*(-1./dScale)));
		ATEST(fpeq(dScale,dActualScale));
		//ATEST(isRotMat(R_from3d));

		getTandR_RANSAC(vPointMatchesRev, R_from3d, T_from3d, dScale);

		ATEST((R_from3d == actualR.t()));
		ATEST((actualT == actualR*T_from3d*(-1./dScale)));
		ATEST(fpeq(dScale, dActualScale));
		//ATEST(isRotMat(R_from3d));

		double err=getErr(vPointMatchesRev[0], R_from3d, T_from3d, dScale);
		ATEST(err>=0);
		C3dPoint * pv1 = const_cast<C3dPoint *>(vPointMatchesRev[0]->p1());

		C3dPoint makeErr(0.05, 0.05, 0.05);
		*pv1 = *pv1 + makeErr;
		double errBigger=getErr(vPointMatchesRev[0], R_from3d, T_from3d, dScale);
		ATEST(err<errBigger);
		ATEST(err<=errBigger);
		*pv1 = *pv1 - makeErr;
		*pv1 = *pv1 * 1.1;
		errBigger=getErr(vPointMatchesRev[0], R_from3d, T_from3d, dScale);
		ATEST(err<errBigger);
		*pv1 /= 1.1;

//////Now introduce outliers and test RANSAC
		for(double z=704; z<706;z+=0.3) //add a few dodgy points
		{
			C3dPoint * pPoint1 = new C3dPoint(z , z , z);
			C3dPoint * pPoint2 = new C3dPoint(22+13*z , -z , z*z*.1+1);
			//PRINTMAT(*pPoint1)
			//PRINTMAT(*pPoint2)
			C3dPointMatch * pMatch = new C3dPointMatch(pPoint1, pPoint2);
			vPointMatches.push_back(pMatch);

			vCleanup.push_back(pPoint1);
			vCleanup.push_back(pPoint2);
		}
		getTandR_RANSAC(vPointMatches, R_from3d, T_from3d, dScale);

		// if these fail we're admitting far-out outliers
		ATEST((R_from3d == actualR));
		ATEST((T_from3d == actualT));
		ATEST(fpeq(dScale*dActualScale, 1));
		//ATEST(isRotMat(R_from3d));

		//Now test 1d alignment:
		ATEST(getTandR_1d_Mean(vPointMatches, actualR, Origin, dActualScale*T_from3d, dScale));

		//ATEST((R_from3d == actualR));
		ATEST(fpeq(dScale*dActualScale, 1));
		//ATEST(isRotMat(R_from3d));

		//Now test 1d alignment:
		ATEST(getTandR_1d_Median(vPointMatches, actualR, Origin, dActualScale*T_from3d, dScale));

		//ATEST((R_from3d == actualR));
		ATEST(fpeq(dScale*dActualScale, 1));*/
		//ATEST(isRotMat(R_from3d));

		//Just 1 point...
		// ATEST(getTandR_1d(vPointMatches[0], actualR, Origin, dActualScale*T_from3d, dScale));

		//ATEST(fpeq(R_from3d, actualR));
		//ATEST(fpeq(dScale*dActualScale, 1));
		//ATEST(isRotMat(R_from3d));

		//1d ransac TODO
		/*
		ATEST(getTandR_1dRANSAC(vPointMatches, actualR, T_from3d*0, dActualScale*T_from3d, dScale));

		ATEST(fpeq(R_from3d, actualR));
		ATEST(fpeq(dScale*dActualScale, 1));
		ATEST(isRotMat(R_from3d));

*/

		//Test HSI=HSV (not HSL) conversions are invertable

		testColConvertHSI(50, 255, 128);
		testColConvertHSI(0, 25, 50);
		testColConvertHSI(255, 25, 50);
		testColConvertHSI(50, 255, 128);
		testColConvertHSI(0, 125, 150);
		testColConvertHSI(255, 225, 200);

		testColConvertRGB(0, 0, 0);
		testColConvertRGB(255, 255, 255);
		testColConvertRGB(100, 100, 100);
		testColConvertRGB(255, 0, 0);
		testColConvertRGB(0, 255, 245);
		testColConvertRGB(0, 0, 255);
		testColConvertRGB(0, 255, 0);
		testColConvertRGB(255, 255, 0);
		testColConvertRGB(255, 0, 245);
		testColConvertRGB(100, 100, 200);
		testColConvertRGB(100, 200, 100);
		testColConvertRGB(200, 100, 100);

		//Test radial distortion:
#define eigAccess(M, x, y) eig ## M(x-1, y-1)
#define M(x, y) eigAccess(M, x, y)

/*	    Eigen::Matrix eigM(3,6);
	    M.row(2) = 1;
	    Matrix calibratedPoints1(2, 6);

	    //Choose 3 image coords in straight line near each axis
		M(1,1)=0;
		M(2,1)=23;
		M(1,2)=0;
		M(2,2)=100;
		M(1,3)=0;
		M(2,3)=123;

		M(1,4)=0;
		M(2,4)=10;
		M(1,5)=100;
		M(2,5)=10;
		M(1,6)=200;
		M(2,6)=10;

		Matrix K(3,3);
		double * adRDcoeffs = 0;
		int nRDcoeffs = 0;

		camParams::HEIGHT = 426;
		getCalibrationMat(K, &adRDcoeffs, nRDcoeffs);

		Matrix K_inv = K.i().rows(1,2);*/

		/*CInvCamCalibMatrix K_inv;
		 *

	    calibratedPoints1 = K_inv.asMat()*M;

	    //Correct RD for calibrated points only--ONLY FOR C2dPoints now
	    //correctRD(calibratedPoints1, adRDcoeffs, nRDcoeffs);

	    //(c-a)*|b-a|/|c-a| should be left/above b
	    ColumnVector c_take_a = calibratedPoints1.column(3) - calibratedPoints1.column(1);
	    ColumnVector b_take_a = calibratedPoints1.column(2) - calibratedPoints1.column(1);
	    ColumnVector b_straightline = calibratedPoints1.column(1) + c_take_a*sqrt(b_take_a.sum_square()/c_take_a.sum_square());

	    / *double distInside = b_straightline(1)-calibratedPoints1(1,2);
	    ATEST(distInside > 0 || nRDcoeffs==0)* /

	    c_take_a = calibratedPoints1.column(6) - calibratedPoints1.column(4);
	    b_take_a = calibratedPoints1.column(5) - calibratedPoints1.column(4);
	    b_straightline = calibratedPoints1.column(4) + c_take_a*sqrt(b_take_a.sum_square()/c_take_a.sum_square());

	    / *distInside = b_straightline(2)-calibratedPoints1(2, 5);
	    ATEST(distInside > 0 || nRDcoeffs==0)* /
*/
	    //5-point AT
	    T2dPoints m0(5);
	    T2dPoints m1(5);
	    //double pdE1[9*10], pdE2[9*10];
	    //double * pdEssentialMat1 = pdE1;
	    //double * pdEssentialMat2 = pdE2;

	    double test_Q1x[5] = {0.6160,    0.8308    ,0.9172 ,   0.7537    ,0.0759};
	    double test_Q1y[5] = {0.4733 ,   0.5853   , 0.2858  ,  0.3804   , 0.0540};

	    double test_Q2x[5] = {0.7792  ,  0.5688  ,  0.3371   , 0.3112  ,  0.6020};
	    double test_Q2y[5] = {0.9340   , 0.4694 ,   0.1622    ,0.5285 ,   0.2630};

	    int nPoints=5;
	    for(int i=0;i<nPoints;i++)
	    {
	    	m0[i] = CSimple2dPoint(test_Q1x[i], test_Q1y[i]);
	    	m1[i] = CSimple2dPoint(test_Q2x[i], test_Q2y[i]);
	    }


	    TSubSet anHypSet(5, 0);
	    for(int i=0;i<5;i++)
	    	anHypSet[i] = i;

//#ifdef TODO_FIX_ATESTS

	    C5ptEssentialMat hypothesise(m0, m1);
	    const T3x3MatModels * models = dynamic_cast<const T3x3MatModels *>(hypothesise.getModels(anHypSet));
	    ATEST(models);

	    //COpenCV8PtEssentialMatModified refiner(m0, m1, 1.5, false);
	    CEssentialMatInlierCounter inlierCounter(0.001, m0, m1);

	    int numMats1 = models->numModels();
	    CMask emptyMask(5);

	    ATEST(numMats1 == 4);
	    for(int i=0;i<numMats1;i++)
	    {
	    	const C3x3MatModel & model = dynamic_cast<const C3x3MatModel &>( models->getData(i) );
	    	ATEST(E_sumSquaredResiduals(model.asDouble9(), m0, m1, emptyMask ) < 0.00002);
	    	//cout << "todo: restore this\n";

	    	ATEST(checkSVs(model));
	    	//pdEssentialMat1 += 9;

	    	//refiner.fitModel(emptyMask, models.getData(i), false)
	    	CSimpleTerminator terminator(m0.size());
	    	int nInlierCount = 0;
	    	emptyMask.setZero();
	    	inlierCounter.countInliers(model, &terminator, 0, 0, emptyMask, 0, nInlierCount, 1.0);

	    	ATEST(emptyMask.countInliers() == m0.size());
	    	ATEST(nInlierCount == m0.size());
	    }

	    //######## New exact data thru matlab
	    double test2_Q1x[5] = {14  ,  79   ,  4   , 68 ,   39};
	    double test2_Q1y[5] = {42  ,  96  ,  85   , 76    ,66};

	    double test2_Q2x[5] = {71  ,   5   , 69  ,   3 ,   77};
	    double test2_Q2y[5] = {3  ,  10  ,  32  ,  44  ,  80};

	    nPoints=5;
	    for(int i=0;i<nPoints;i++)
	    {
	    	m0[i] = CSimple2dPoint(test2_Q1x[i], test2_Q1y[i]);
	    	m1[i] = CSimple2dPoint(test2_Q2x[i], test2_Q2y[i]);
	    }

	    const T3x3MatModels * models2 = dynamic_cast<const T3x3MatModels *>(hypothesise.getModels(anHypSet));
	    ATEST(models2)

	    numMats1 = models2->numModels();
	    ATEST(numMats1 == 4);
	    for(int i=0;i<numMats1;i++)
	    {
	    	const C3x3MatModel & model = dynamic_cast<const C3x3MatModel &>( models2->getData(i) );
	    	ATEST(checkSVs(model));

	    	CSimpleTerminator terminator(m0.size());
	    	int nInlierCount = 0;
	    	emptyMask.setZero();
	    	/*const C3x3MatModel & modelOld = dynamic_cast<const C3x3MatModel &>( models.getData(i) );
	    	inlierCounter.countInliers(modelOld, &terminator, 0, 0, emptyMask, 0, nInlierCount);
	     	ATEST(emptyMask.countInliers() == 0); //No inliers to original models*/

	    	inlierCounter.countInliers(model, &terminator, 0, 0, emptyMask, 0, nInlierCount, 1.0);
	    	ATEST(emptyMask.countInliers() == m0.size());
	    	ATEST(nInlierCount == m0.size());
	    	emptyMask.setZero();

	    	//Check BGC doesn't break it
	    	inlierCounter.countInliers(model, &terminator, 0, 4 /*BGC*/, emptyMask, 0, nInlierCount, 1.0);
	    	ATEST(emptyMask.countInliers() == m0.size());
	    	ATEST(nInlierCount == m0.size());
	    }

//#endif
	    //############

	    //If we have 6 points or more, soln is a lin. comb. of just a few cols of E

//	    give this a rest and implement relative positions. Only recompute when a new position is calculated, propogate changes (if there has actually been a change). Should never have to propagate very far.

	    //Compare essential matrix recovery methods:
	    //Simulate 2 calibrated images:
	    CPointVec2d points1(0), points2(0);

	    double adArrLikelihood[1000];
		C3dPoint t(0.01 , 0.02 , -0.03);

		//R=getRotMat(0.3, 0.4, 0.5);
		CCamera P = R0 | Origin;
		CCamera Pp = R | t;

		//Matrix E = (Xmat(t) * Pp.rotation().asMat());

		const double CUBE_RAD = 0.2;
		const double CUBE_STEP = 0.04;
		const double CUBE_DEPTH = 0.6 ;

		int numPoints = 0;
		for(double z=CUBE_DEPTH-CUBE_RAD; z<=CUBE_DEPTH+CUBE_RAD; z+=CUBE_STEP)
			for(double x=-CUBE_RAD; x<=CUBE_RAD; x+=CUBE_STEP)
				for(double y=-CUBE_RAD; y<=CUBE_RAD; y+=CUBE_STEP)
				{
					C3dPoint Q(x,y,z); //3d scene point
 					C2dPoint p = Q.photo(P);// p/=p(3);
					C2dPoint pp = Q.photo(Pp);// pp/=pp(3);

					/*double reproj = epipolarResid(pp, E, p);
					ATEST( zero(reproj) ) //this failing shows the autotest is wrong*/

					///////////////// Check reconstruction is working:
					C3dPoint Q_recon=reconstruct(P, Pp, p, pp);

					ATEST((Q_recon-Q).sum_square()/Q.sum_square() < 0.01);

					////////////////////
					//if(bClip) if(clip(L1) || clip(L2)) continue; //Apply clipping frame



					ATEST(Q.testInFront(P, Pp))
					ATEST(Q.depth(P) > 0)
					ATEST(Q.depth(Pp) > 0)

					ATEST(testPair(P, Pp, p, pp))

					//int nLikelihoodInv = 1;
					//pCorr->push_back(new CCorrespondence(L1,L2, 1, nLikelihoodInv));

					points1.push_back(p);
					points2.push_back(pp);
					adArrLikelihood[numPoints] = 0.5;
					numPoints++;
					//also fill array of locs TODO

					/*for(int i=1;i<numPoints;i++)
					{
						if(points1.column(i) == points1.column(numPoints))
							cout << "Col " << i << " matches " << numPoints << ": " << p << endl;
					}*/
				}

		/*CvMat * cvF=cvCreateMat(3,3,CV_64FC1);
		CvMat * cvArrStatus = 0;//cvCreateMat(1,numPoints,CV_8UC1);
		int numOutliers = 0;


		eSamplingMethod eSM = eSM_RANSAC;
		//RANSAC_VERIFICATION_METHOD = 0; cout << "TO DO";
		//cout << points1;

		for(int nOutlierMode=0; nOutlierMode<2; nOutlierMode++)
		{
			CvMat * cvPoints1 = pointVectorToCvMat(points1);
			CvMat * cvPoints2 = pointVectorToCvMat(points2);
			eFAlgorithm aAlgs[3] =  { e5point , e7point, e8point  };

			for(int nHypAlgIdx = 0; nHypAlgIdx<1; nHypAlgIdx++)
			{
				eFAlgorithm eFHypAlg = aAlgs[nHypAlgIdx];
				for(int nRefineAlgIdx = 2; nRefineAlgIdx<3; nRefineAlgIdx += 2)
				{
					eFAlgorithm eFRefineAlg = aAlgs[nRefineAlgIdx];
					for(int numPointsHyp = (int)eFHypAlg; numPointsHyp<(int)eFHypAlg+1; numPointsHyp++)
					{
						for(int bGB = 0; bGB<1; bGB++)
						{
							if(bGB && eFRefineAlg != e5point && eFHypAlg != e5point) continue;

							int nInliers=-1;

							//if(CRandom::Uniform(10)!=0) continue; //Save time...

							bool bSuccessFindingE=BoW_DisjointFindFundamentalMat(cvPoints1, cvPoints2, cvF, CV_FM_RANSAC | CV_FM_8POINT, RANSAC_E_THRESH , RANSAC_E_PROB, cvArrStatus, &nInliers, true, eBestMatchOnly, numPointsHyp, eFHypAlg, eFRefineAlg, bGB, eSM, RANSAC_8PT_SV_CUTOFF, (eVerificationMethod) (RANSAC_VERIFICATION_METHOD), 0, RANSAC_TOPDOWN_REFINEMENT, RANSAC_MAX_ITERS, adArrLikelihood, 3);

							ATEST(bSuccessFindingE)
							if(numOutliers==0)
							{
								if(!(nInliers <= numPoints && nInliers > (int)(numPoints*0.98)))
									cout << nInliers << " inliers found\n";

								ATEST(nInliers <= numPoints && nInliers > (int)(numPoints*0.98))
							}
							else
							{
								bool bInlierTest = nInliers <= numPoints && nInliers > (int)((numPoints-numOutliers)*0.97);
								if(!bInlierTest)
									cout << nInliers << " found: ";
								ATEST(bInlierTest)
							}
						}
					}
				}
			}
			cvReleaseMat(&cvPoints1);
			cvReleaseMat(&cvPoints2);
			numOutliers = 250;
			//points1.submatrix(1, 1, 100, 100+numOutliers) *= 0.7; //NB these are systematically wrong
			for(int i=0;i<100;i++)
				points1[i] *= 0.7;
			//points1.submatrix(1, 1, 100, 100+numOutliers/2) *= 2.7;
		}
		*/

		//testFindBestMatching();

		testFastSet<CFastSet<int, 3> >(false);
		testFastSet<CSmallHashTable<int, 23, 1, intHash<23> > >(false);
		testFastSetLocations<CFastSet<CLocation, 96> >();
		testFastSetLocations<CSmallHashTable<CLocation, 13, 2, locationHash<13> > >();

		testFastSet<CSimpleSet<int> >(false);
		testFastSetLocations<CSimpleSet<CLocation> >();
		TFS(testStdSet<int>);
		TFS(set2<int>);

		typedef set2<int, std::less<int> > TSetBF128;
		TFS(TSetBF128);
		typedef set2<int, std::less<int> > TSetBF1024;
		TFS(TSetBF128);
		typedef set2<int, std::less<int> > TSet16;
		TFS(TSet16);
		typedef set2<int, std::less<int> > TSet128;
		TFS(TSet128);
		typedef set2<int, std::less<int> > TSet1024;
		TFS(TSet1024);
		typedef set2<int, std::less<int> > TSetNF128;
		TFS(TSetNF128);
		typedef set2<int, std::less<int> > TSetNF1024;
		TFS(TSetNF1024);

		/*TFS(set2<int, std::less<int>, CBitfieldIndividualPoolAllocator< T, 1024 > >);
		TFS(set2<int, std::less<int>, CIndividualPoolAllocator< T, 128 > >);
		TFS(set2<int, std::less<int>, CIndividualPoolAllocator< T, 1024 > >);
		TFS(set2<int, std::less<int>, CIndividualPool_NoFree_Allocator< T, 128 > >);
		TFS(set2<int, std::less<int>, CIndividualPool_NoFree_Allocator< T, 1024 > >);
		TFS(set2<int, std::less<int>, CIndividualPoolAllocator< T, 16 > >);*/

		testLognormal();

		//Normal distns calc
		double dMean=0, dVar=0;
		CVariablesFromMultiDistnsObserved observer;
		observer.addDist(23,34);
		observer.computeMeanVar(dMean, dVar);
		ATEST(dMean == 23);
		ATEST(dVar == 34);
		observer.addDist(23,34);
		observer.addDist(23,34);
		observer.computeMeanVar(dMean, dVar);
		ATEST(dMean == 23);
		ATEST(dVar < 34);
		observer.reset();
		observer.addDist(1,4);
		observer.addDist(1,9);
		observer.computeMeanVar(dMean, dVar);
		ATEST(dMean == 1);
		ATEST(dVar < 4);
		observer.reset();
		observer.addDist(1,4);
		observer.addDist(1,9);
		observer.computeMeanVar(dMean, dVar);
		ATEST(dMean == 1);
		ATEST(dVar < 4);
		observer.reset();
		observer.addDist(1,9);
		observer.addDist(2,9);
		observer.addDist(3,9);
		observer.computeMeanVar(dMean, dVar);
		ATEST(dMean == 2);
		ATEST(dVar == 3);

		observer.reset();
		observer.addDist(10,11);
		observer.addDist(12,13);
		observer.computeMeanVar(dMean, dVar);
		double dMean2=0, dVar2=0;
		CVariablesFromMultiDistnsObserved::combineTwoObservations(10,11,12,13, dMean2, dVar2);
		ATEST(zero(dMean - dMean2));
		ATEST(zero(dVar - dVar2));
//return true;
		/*CScaleFromObservingAndMeasuringKnownObjs scaleObserver;
		scaleObserver.addObjectMeasurement(5, 10, 324);
		double dScale = 0;
		scaleObserver.computeMeanVar(dScale, dVar);
		ATEST(dScale == 2);

		scaleObserver.addObjectMeasurement(5, 10, 324);
		scaleObserver.computeMeanVar(dScale, dVar);
		ATEST(dScale == 2);

		scaleObserver.addObjectMeasurement(23, 46, 100);
		scaleObserver.computeMeanVar(dScale, dVar);
		ATEST(dScale == 2);

		scaleObserver.addObjectMeasurement(23, 45, 100);
		scaleObserver.computeMeanVar(dScale, dVar);
		ATEST(dScale < 2);

		scaleObserver.addObjectMeasurement(23, 46, 80);
		scaleObserver.computeMeanVar(dScale, dVar);
		ATEST(dScale < 2);

		scaleObserver.addObjectMeasurement(23, 47, 80);
		scaleObserver.computeMeanVar(dScale, dVar);
		ATEST(dScale > 2);
		//ATEST(dVar < 80);

		scaleObserver.addObjectMeasurement(5, 10,324);*/


  		CScaleMLE::testScaleMLE();


		//TODO......
	}
    catch(CException pEx)
    {
		cout << "AUTOTEST FAILED\n";
        cout << pEx.GetErrorMessage() << '\n';
        cout.flush();
        return false;
    }
    cout << "Autotest succeeded\n";
    return true;
}

void testVectors()
{
	C3dPoint p1(1,2,3);
	C3dPoint p2=p1+p1*3;
	//cout << p2;
	ATEST(p1+p2 == p2+p1)

	//Matrix R=getRotMat(0.3, 0.4, 0.5);
	C3dPoint axis(23,45,-1);
	C3dRotation R2(axis, 1.1);

	//cout << R;
	//cout << R2 << endl;

	C3dRotation Id;
	ATEST(R2*R2.t() == Id)
	ATEST(R2.t()*R2 == Id)
	ATEST(Id * R2 == R2)
	ATEST(R2 * Id == R2)
	ATEST(Id * p2 == p2)

	//cout << R2*p2;
	ATEST(R2.t().t() == R2)

	ATEST(R2 * (R2.t() * p2) == p2)

	cout << R2.t()*R2*p2;

	C3dPoint p3 = p2;
	ATEST(p2 == p3)
	p2.rotate(R2);

	ATEST(p2 == R2 * p3)
	//cout << endl << p2;
	p2.rotateInv(R2);

	//cout << endl << p2 << endl << p3 << endl;
	ATEST(p2 == p3)

	C3dRotation R3(-axis, 0.1), R4, R5(-axis, 0.01);
	for(int i=0;i<10;i++)
	{
		R4 = R4 * R3;
		R4 = R4 * R5;
	}

	ATEST(zero(1.1-R4.angle()))
}

/*#define EUCLIDNORM(an1, an2, i, len) \
	d += sqr(an1[i]-an2[i]); \
	EUCLIDNORM(an1, an2, i+1, len)*/
//template<int len, int i> inline void euclidNormExpanded(int * an1, int *an2, double & d);

//template<int len> double euclidNormExpand(int * an1, int *an2, int);

template<int len, int i> class euclidNormExpanded;

template<int len>
class euclidNormExpanded<len, 0>
{
public:
	inline static void addNext(int *, int *, int &) {}
};

template<int len, int i> class euclidNormExpanded
{
	friend class euclidNormExpanded<len, i+1>;
public:
	inline static void addNext(int * an1, int *an2, int & d)
	{
		d += sqr(an1[len-i]-an2[len-i]);

		if (i>0)
		{
			euclidNormExpanded<len, i-1>::addNext(an1, an2, d);
		}
	}
};

template<int len> int euclidNormExpand(int * an1, int *an2, int)
{
	int d=0;
	euclidNormExpanded<len, len>::addNext(an1, an2, d);
	return d;
}
template<int len> int euclidNormExpand2(int * an1, int *an2, int)
{
	int d=0;
	euclidNormExpanded<len/4, len/4>::addNext(an1, an2, d);
	an1 += len/4; an2 += len/4;
	euclidNormExpanded<len/4, len/4>::addNext(an1, an2, d);
	an1 += len/4; an2 += len/4;
	euclidNormExpanded<len/4, len/4>::addNext(an1, an2, d);
	an1 += len/4; an2 += len/4;
	euclidNormExpanded<len/4, len/4>::addNext(an1, an2, d);
	return d;
}

void testNorms()
{
#define len 172
	int * an1 = new int[len+200];
	int * an2 = new int[len+200];
	int dist = 0;
	for(int i=0;i<len+200; i++)
		{ an1[i] = i; an2[i] = an1[i]-10; }

	CStopWatch s;

#define TIMENORM(NORM) \
	dist=0;\
	s.startTimer(); \
	for(int i=0;i<10; i++) \
	    dist += NORM(an1+((i*20)%200), an2+(i%200), len); \
	s.stopTimer(); \
	if(dist) \
		cout << # NORM << s.getElapsedTime() << ' ' << dist << endl;

	TIMENORM(L1Dist)
	TIMENORM(L1DistSlower)
	TIMENORM(euclidDist)
	TIMENORM(MaxDist)
	TIMENORM(MaxDistSlower)
	TIMENORM((euclidDistUnroll<int, 2, len>))
	TIMENORM((euclidDistUnroll<int, 3, len>))
	TIMENORM((euclidDistUnroll<int, 4, len>))
	//TIMENORM((euclidDistUnroll4<int, len>))
	TIMENORM((euclidDistUnroll<int, 5, len>))
	TIMENORM((euclidDistUnroll<int, 6, len>))
	TIMENORM((euclidDistUnroll<int, 7, len>))
	TIMENORM((euclidDistUnroll<int, 8, len>))
    TIMENORM((euclidDistUnroll<int, 9, len>))
	TIMENORM(euclidNormExpand<len>)
    TIMENORM(euclidNormExpand2<len>)
}
