/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "convert.h"
#include "random.h"

int intLookup::anSQRT[intLookup::nLookupSize] = {};
bool intLookup::bSetup = false;
COUNTHITS(int intLookup::nHit = 0; int intLookup::nMiss = 0)
double intLookup::adUcharToDouble[256] = {};
float intLookup::afUcharToFloat[256] = {};

unsigned CRandom::g_seed = 0;
const double CRandom::NEXT_NORMAL_NOT_SET = 100000;
double CRandom::s_dNextNormal = CRandom::NEXT_NORMAL_NOT_SET;
const float CRandom::fNEXT_NORMAL_NOT_SET = 100000;
float CRandom::s_fNextNormal = CRandom::fNEXT_NORMAL_NOT_SET;

double factorial(double n)
{
    if (n<=1)
        return 1.0;
    return n*factorial(n-1.0);
}

double nCr(double n,double r)//could be implemented efficiently...
{
    if(r>n) return 0;
    return factorial(n)/(factorial(r)*factorial(n-r));
}

void breakInCpp()
{
    cout << "Break in convert.cpp" << endl;
    //raise(SIGINT);
    static int dummy=0; 
    dummy++;

#ifdef WIN32
    //__debugbreak();
#endif
}

