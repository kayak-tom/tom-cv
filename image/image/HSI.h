/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

//Declaration of class CHsi
// HSI.h file
//
#ifndef __HSI_H__
#define __HSI_H__

#include "util/exception.h"
#include "util/convert.h"

#define HSI_Max(x,y) ((x)>(y) ? (x) : (y))
#define HSI_Min(x,y) ((x)<(y) ? (x) : (y))

//Struct really
class CHSI
{
public:
    int H, S, I;

    CHSI() : H(-1), S(-1), I(-1) {}
    CHSI(int H, int S, int I) : H(H), S(S), I(I) {}

    static double distance(const CHSI & col1, const CHSI & col2);

};

class CRedGreenBlue
{
public:
    int R, G, B;

    CRedGreenBlue() : R(-1), G(-1), B(-1) {}
    CRedGreenBlue(int R, int G, int B) : R(R), G(G), B(B) {}
    CRedGreenBlue(const unsigned char * pBGR) : R(pBGR[2]), G(pBGR[1]), B(pBGR[0]) {}
};

class CXYZ
{
public:
    int X, Y, Z;

    CXYZ() : X(-1), Y(-1), Z(-1) {}
    CXYZ(int X, int Y, int Z) : X(X), Y(Y), Z(Z) {}
    CXYZ(const CHSI & colHSI) : X(-1), Y(-1), Z(-1)
    {
        Z=colHSI.I;

        //Convert polar H, S to cartesian:
        double dHueAngle = colHSI.H * (2 * M_PI / 255.0);
        double dHueCircleRadius = colHSI.S*fabs((double)(colHSI.I - 128))/255.0;
        X = doubleToInt(dHueCircleRadius*sin(dHueAngle));
        Y = doubleToInt(dHueCircleRadius*cos(dHueAngle));
    }
};

class CHsi
{

public:
  CHsi();
  virtual ~CHsi();
  void HSI(int r,int g,int b,\
          double &Hue, double &Saturation, double &Intensity);
private:
  int nImax,nImin,nSum,nDifference;

};

#ifdef RGB
#undef RGB
#endif

class CHsiInt
{
public:
  static void HSI( int r, int g, int b,
           int &Hue, int &Saturation, int &Intensity);
  static void RGB( int h, int s, int i,
           int &R, int &G, int &B);
  static void RGBfromHSV( int h, int s, int i,
           int &nR, int &nG, int &nB);


  static void HSI( const CRedGreenBlue & inRGB,
          CHSI & outHSI) { return HSI(inRGB.R, inRGB.G, inRGB.B, outHSI.H, outHSI.S, outHSI.I ); }

  static void RGB( const CHSI & inHSI,
          CRedGreenBlue & outRGB) { RGB( inHSI.H, inHSI.S, inHSI.I, outRGB.R, outRGB.G, outRGB.B ); }
  static void RGBfromHSV(const CHSI & inHSI,
          CRedGreenBlue & outRGB )  { RGBfromHSV( inHSI.H, inHSI.S, inHSI.I, outRGB.R, outRGB.G, outRGB.B ); }
};

std::ostream& operator<<(std::ostream& s, const CRedGreenBlue & X);


#endif// __HSI_H__
