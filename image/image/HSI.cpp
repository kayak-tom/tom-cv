/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

//Implementation of class CHsi
#include "HSI.h"
#include "util/exception.h"
#include "util/convert.h"

double CHSI::distance(const CHSI & col1, const CHSI & col2)
{
    CXYZ xyz1(col1);
    CXYZ xyz2(col2);
    return (sqr(xyz1.X-xyz2.X) + sqr(xyz1.Y-xyz2.Y) + sqr(xyz1.Z-xyz2.Z));
}


CHsi::CHsi()
{

}

CHsi::~CHsi()
{

}

void CHsi::HSI(int r, int g, int b,
           double &Hue, double &Saturation, double &Intensity)
{

  if( (r<0 && g<0 && b<0) || (r>255 || g>255 || b>255) )
  {
    Hue=Saturation=Intensity=0;
    return;
  }

  if(g==b)
  {
    if(b<255)
    {
      b=b+1;
    }
    else
    {
      b=b-1;
    }
  }

  nImax = HSI_Max(r,b);
  nImax = HSI_Max(nImax,g);
  nImin = HSI_Min(r,b);
  nImin = HSI_Min(nImin,g);
  nSum = nImin+nImax;
  nDifference =nImax-nImin;

  Intensity = (float)nSum/2;

  if(Intensity<128)
  {
    Saturation=(255*((float)nDifference/nSum));
  }
  else
  {
    Saturation=(float)(255*((float)nDifference/(510-nSum)));
  }

  if(Saturation!=0)
  {
    if(nImax == r)
    {
      Hue=(60*((float)g-(float)b)/nDifference);
    }
    else if(nImax == g)
    {
      Hue=(60*((float)b-(float)r)/nDifference+120);
    }
    else if(nImax == b)
    {
      Hue=(60*((float)r-(float)g)/nDifference+240);
    }

    if(Hue<0)
    {
      Hue=(60*((float)b-(float)r)/nDifference+120);
    }
  }
  else
  {
    Hue=0;
  }
}

#define SETCOL(c) \
if(t ## c < 256/6) \
    c = p + ((q-p)*6*t ## c)/255; \
else if (t ## c < 256/2) \
    c = q; \
else if (t ## c < (2*256)/3) \
    c = p + ((q-p)*6*((2*256)/3 - t ## c))/255; \
else \
    c = p;

void CHsiInt::RGB(const int h,const int s,const int i,
         int &r, int &g, int &b)
{
    // Use http://en.wikipedia.org/wiki/HSI_color_space#Conversion_from_RGB_to_HSL_or_HSV
    if(s==0) r=g=b=i;

    int p, q, tr, tg, tb;
    if(i < 128)
    {
        q = (i*(255+s))/255;
    }
    else
    {
        q = (i+s) - (i*s)/255;
        if(q>255) q=255;
        if(q<0) q=0;
    }
    p = 2*i - q;
    tr = h + 255/3; if(tr>255) tr -= 255;
    tg = h;
    tb = h - 255/3; if(tb<0) tb += 255;

/*    if(tr < 256/6)
        r = p + ((q-p)*6*tr)/255;
    else if (tr < 256/2)
        r = q;
    else if (tr < (2*256)/3)
        r = p + ((q-p)*6*((2*256)/3 - tr))/255;
    else
        r = p;*/
    SETCOL(r)
    SETCOL(g)
    SETCOL(b)
}

void CHsiInt::RGBfromHSV( int h, int s, int i,
         int &nR, int &nG, int &nB)
{

    // ######################################################################
    // T. Nathan Mundhenk
    // mundhenk@usc.edu
    // C/C++ Macro HSV to RGB TODO: To int
    static const double intTo01 = 1.0/255.0;
    const double H=h*intTo01*6.0, S=s*intTo01, V=i*intTo01;
    double R,G,B;

    if( V == 0 )
    { R = 0; G = 0; B = 0; }
    else if( S == 0 )
    {
      R = V;
      G = V;
      B = V;
    }
    else
    {
      const double hf = H ; // / (256.0/6.0); //Was /60
      const double floor_hf = floor( hf );
      const int    i  = (int) floor_hf;
      const double f  = hf - floor_hf;
      const double pv  = V * ( 1 - S );
      const double qv  = V * ( 1 - S * f );
      const double tv  = V * ( 1 - S * ( 1 - f ) );
      switch( i )
        {

    //Red is the dominant color

        case 0:
          R = V;
          G = tv;
          B = pv;
          break;

    //Green is the dominant color

        case 1:
          R = qv;
          G = V;
          B = pv;
          break;
        case 2:
          R = pv;
          G = V;
          B = tv;
          break;

    //Blue is the dominant color

        case 3:
          R = pv;
          G = qv;
          B = V;
          break;
        case 4:
          R = tv;
          G = pv;
          B = V;
          break;

    //Red is the dominant color

        case 5:
          R = V;
          G = pv;
          B = qv;
          break;

    //Just in case we overshoot on our math by a little, we put these here. Since its a switch it won't slow us down at all to put these here.

        case 6:
          R = V;
          G = tv;
          B = pv;
          break;
        case -1:
          R = V;
          G = pv;
          B = qv;
          break;

    //The color is not defined, we should throw an error.

        default:
          THROW("i Value error in Pixel conversion");
          break;
        }
    }
    R *= 255.99;
    G *= 255.99;
    B *= 255.99;

    nR=doubleToInt(R);
    nG=doubleToInt(G);
    nB=doubleToInt(B);
}

void CHsiInt::HSI( int r, int g,  int b,
           int &Hue,  int &Saturation,  int &Intensity)
{
  int nImax,nImin,nSum,nDifference;

  //if( (r<0 && g<0 && b<0) || (r>255 || g>255 || b>255) )
  //{
  //  Hue=Saturation=Intensity=0;
  //  return;
  //}

  if(g==b) //hack to make g!=b
  {
    if(b<255)
    {
      b++;
    }
    else
    {
      b--;
    }
  }

  nImax = HSI_Max(r,b);
  nImax = HSI_Max(nImax,g);
  nImin = HSI_Min(r,b);
  nImin = HSI_Min(nImin,g);
  nSum = nImin+nImax;
  nDifference =nImax-nImin;

  Intensity = nSum/2; //dubious

  if(Intensity<128) //Definitely HSL colour space, agrees with Wikipedia
  {
    Saturation=(255*nDifference)/nSum;
  }
  else
  {
    Saturation=(255*nDifference)/(510-nSum);
  }

  //const int SEGMENTSIZE = 60; for hue in 0..360
  const int SEGMENTSIZE = 42; //256/6;

  if(Saturation!=0)
  {
    if(nImax == r)
    {
      Hue=(SEGMENTSIZE*(g-b))/nDifference;
      if (Hue < 0 ) Hue += 256;
    }
    else if(nImax == g)
    {
      Hue=(SEGMENTSIZE*(b-r))/nDifference + 2*SEGMENTSIZE;
    }
    else if(nImax == b)
    {
      Hue=(SEGMENTSIZE*(r-g))/nDifference + 4*SEGMENTSIZE;
    }

    if(Hue<0)
    {
      Hue=(SEGMENTSIZE*(b-r))/nDifference + 2*SEGMENTSIZE;
    }
  }
  else
  {
    Hue=0;
  }

  return;
}

std::ostream& operator<<(std::ostream& s, const CRedGreenBlue & X)
{
   s << "" << X.R << ", " << X.G << ", "  << X.B << " " << std::flush;
   return s;
}
