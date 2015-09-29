/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */
/* INCLUDE THIS FILE BEFORE ANY OTHER INCLUDES!! */
#pragma once  

#ifndef EXCEPTION_H  
#define EXCEPTION_H

#define NOMINMAX // otherwise windows defines min, max, which conflicts with Eigen

#ifndef __GNUC__
#pragma warning( disable : 4127 ) //conditional expression is constant
#pragma warning( disable : 4800 ) //performance warning converting value to bool

#endif
 
#ifdef __GNUC__ 
#  ifndef __OPTIMIZE__
#    ifndef _DEBUG
#      define _DEBUG //if -O0 set.
#    endif
#  endif
#endif

#if defined(_DEBUG) || defined(BOWDLL_EXPORTS)
#  define IS_DEBUG true
#  define DEBUGONLY(...) __VA_ARGS__
#else
#  ifdef __EXCEPTIONS
//#    undef __EXCEPTIONS //STL assertions
#  endif
#  define IS_DEBUG false
#  define DEBUGONLY(...)
#  define EIGEN_NO_DEBUG //Disable some Eigen assertions

#  ifndef NDEBUG
#    define NDEBUG //Disable STL assertions
#  endif

#  define EIGEN_VECTORIZE

#  ifndef BOOST_DISABLE_ASSERTS
#    define BOOST_DISABLE_ASSERTS //Array bound checking on scoped_array. 
#  endif
#endif

inline bool UNINIT(const void * p) //heuristic indicating pointer that's been overwritten
{
    return (p) < (const void *)1000000;
}

#include <iostream>
#include <signal.h>
#include <string> //use strings in exception handling
#include <stdlib.h>
#include <stdexcept>
#include <string.h>

//Actually breakInCpp is useful in Windows as well #ifndef _WIN32
#define BREAK_IN_CPP
//#endif

void breakInCpp(); //Dummy function so a breakpoint can be placed in a cpp file for netbeans.

class CException : public std::runtime_error
{
protected:
    
public:
    CException(const std::string & errorString, const char * szFnName = 0) : std::runtime_error(std::string("ERROR: ") + (szFnName ? szFnName : "") + " " + errorString)
    {

        std::cout << getErrorString() << std::endl;
        std::cerr << getErrorString() << std::endl; //because we sometimes direct cout only to file.

#ifdef BREAK_IN_CPP
        breakInCpp();
#endif
    }
    CException() : std::runtime_error("") {}
    //CException(const CException & ex) : szError(ex.GetErrorMessage()) {}
    /*void operator=(const CException & ex)
    {
        *const_cast<char **>( &szError ) = const_cast<char *>(ex.GetErrorMessage());
    };*/

    const char * GetErrorMessage() const { return what(); };
    std::string getErrorString(/*const char * szFnName=0*/) const
    {
        return what();
    }

    operator bool() const { return strlen(what()) > 0; }
};

#ifdef __GNUC__
#  define THROW_int(s) throw CException(s, __PRETTY_FUNCTION__);
#else
#  define THROW_int(s) throw CException(s);
#endif

#define ex_OP(z) ex_OP2(__FILE__, __LINE__, z)
#define ex_OP2(x,y,z) ex_OP3(x, y, z)
#define ex_OP3(x,y,z) x ":" #y ": " + std::string(z)

#define CHECK_PROFILING // rand(); // Define this to a non-pure function so that the callgrind profile shows up how often CHECK conditions are been hit, in case there are any which are slowing the program down.

#define THROW(s) THROW_int(ex_OP(s))
#define CHECK(x, s) { if(x) THROW(s)}
#define CHECK_P(x, expr, s) {CHECK_PROFILING; if(x) { std::cout << "CHECK_P: " << #expr << "=" << expr << std::endl; std::cerr << "CHECK_P: " << #expr << "=" << expr << std::endl;  THROW(s); }}
#define CHECKNAN(x) {CHECK_PROFILING; CHECK(std::isnan((double)x) || std::isinf((double)x), #x " is inf or nan");}
#define CHECKBADNUM(x) {CHECKNAN(x); const double xabs=std::fabs((double)x); CHECK_P((xabs>0&&xabs<1e-50) || xabs>1e+100, x, #x " is *probably* uninitialised (denormal, or >>HUGE)")};
#define CHECKNOTNULL(x) {CHECK_PROFILING; CHECK(!x, #x " is null");}
#define CHECKPROBABILITY(x) { CHECKNAN(x); /*do not CHECKBADNUM because probs can get tiny */ CHECK(x<0 || x>1, #x " is not a probability"); }
#define CHECKEQUAL(x,y) { if(!zero(x-y)) { std::cout << #x "=" << x << " != " #y "=" << y << std::endl; THROW("Arguments not equal"); } }
#define CHECKOPTIONAL(x) CHECKNOTNULL(x) //for boost::optional when we're expecting a return
#define CHECK_MAT_INIT(M) CHECK(M.size().area() == 0, "cv::Mat " #M " is uninitialised")

/*#ifdef __GNUC__
#define CHECKOOB(n, N) { if(n<0 || n>=(typeof(n)) N) { cout <<  #n " is OOB: " #n " = " << n << ", upper bound " << #N " = " << N << endl; THROW("Index OOB");  } } 
#else This one should work with C++11... */
#define CHECKOOB(n, N) { if((int)n<0 || (int)n>=(int) N) { cout <<  #n " is OOB: " #n " = " << n << ", upper bound " << #N " = " << N << endl; THROW("Index OOB");  } } 
//#endif
#define CHECKEXISTS(vec, el) {  if(std::find(vec.begin(), vec.end(), el) != vec.end()) { cout << #el << " = " << el << " already exists in " << #vec << endl; THROW("Element already exists"); } } 

#define pragma_warning(x) //Unfortunately have to exclude each warning in __GNUC__ individually

#ifndef __GNUC__

#  define ARRAY(T, atName, length) boost::scoped_array<T> atName(new T[length]) //__GNUC__ allows int anArr[nLength] when nLength is const
#  define PTR(array) array.get()

#else

#  define sprintf_s(txt,len,format,...) sprintf(txt,format,##__VA_ARGS__)

#  define ARRAY(T, atName, length) T atName[length] //__GNUC__ allows int anArr[nLength] when nLength is const
#  define PTR(array) array

#endif

using std::cout;
using std::endl;

#define REPEAT(times, ...) { static int nRepeatCount = times; if(nRepeatCount>0) { nRepeatCount--; {__VA_ARGS__;} }}
#define PERIODIC(freq, ...) { static int nRepeatCount = 0; nRepeatCount++; if(nRepeatCount % freq == 0)  {__VA_ARGS__;}}
class CVALGRIND_CHECK_INIT {
public: static void checkInit(const double d) { if(d == 123456789) cout << "Detecting use of uninitialised variables in valgrind (this actually prints sometimes?!): d=" << d << endl; }
};
#define VALGRIND_CHECK_INIT(d) CVALGRIND_CHECK_INIT::checkInit(d);

//#define CHECK_SIZES_EQUAL(T1,T2) _x_CCASSERT_LINE_CAT(sizeof(T1)==sizeof(T2), T1, T2);
//#define eSTR(x) # x
//#define eXSTR(x) eSTR(x)
//#define CHECK_SIZES_EQUAL_int(T1, T2) CHECK_SIZES_EQUAL(T1, T2, "x" eXSTR(__LINE__));
#define CHECK_SIZES_EQUAL(T1, T2, id) \
typedef char Sizes_must_be_equal_##id[2*((sizeof(T1)==sizeof(T2)))-1];

#define CHECK_SIZES_EQUAL_RT(T1, T2, id) CHECK(sizeof(T1)!=sizeof(T2), "Sizes of T1 and T2 must be equal");


#endif // EXCEPTION_H
