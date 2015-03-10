#ifndef CTHREADPOOL_H
#define CTHREADPOOL_H

//This file is generally not required for geometry library (remove if you're using an old/incompatible compiler)
#ifndef __GNUC__
#  include <functional>
#else
#  include <tr1/functional>
#endif

typedef void TNullaryFn(void);
typedef std::tr1::function<TNullaryFn> TNullaryFnObj;

class CThreadpool_base
{
public:
    CThreadpool_base() {}
    virtual ~CThreadpool_base() { }
    
    static CThreadpool_base * makeThreadpool(const int nNumThreads);
    
    virtual void addJob(TNullaryFnObj & function) = 0;
    
    //Block until everything has finished
    virtual void waitForAll(TNullaryFnObj mainThreadFn = (TNullaryFnObj)(nothingFn)) = 0;
    
    virtual int getNumThreads() const = 0;

    static void nothingFn() {}
};

#endif // CTHREADPOOL_H
