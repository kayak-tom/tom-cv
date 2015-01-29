#include <util/exception.h>
#include "threadpool.h"
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

static const char * szEXITING = "Exiting from worker thread", * szUnknown = "Unknown exception caught from job run by threadpool";

/**
 * @class CThreadpool
 * @brief Locking with semaphores. 1 thread is a special case: the functions are called by waitForAll(). This makes 
 * valgrind/callgrind output much nicer, and means that OpenCV showim functions will work when MT is turned off.
 */
class CThreadpool : public CThreadpool_base
{
    boost::interprocess::interprocess_semaphore semaphore_workers, semaphore_waitForAll;
    boost::mutex mxLockVector;
    boost::thread_group aThreads;
    
    std::vector<TNullaryFnObj> aJobs;
    
    bool bExiting;
    
    const char * szErrorMessage;
    
    //Return false if there's no jobs left
    bool onejobHandler()
    {
        //If there's a queued job then do it, else return 'done'
        TNullaryFnObj fn;
        {
            boost::mutex::scoped_lock scopedLock(mxLockVector);
            if(aJobs.size() == 0)
                return false;
            
            fn = aJobs.back();
            aJobs.pop_back();
        }
        fn();
        
        return true;
    }
    
    //Return true when we should exit
    void doJobs()
    {
        try
        {
            for(;;)
            {
                if(!onejobHandler())
                    return;
            }
        }
        catch(const CException & ex)
        {
            szErrorMessage = ex.GetErrorMessage();
            if(!szErrorMessage)
                szErrorMessage = szEXITING;
        }
        catch(const std::exception & ex)
        {
            cout << "ERROR: std::exception caught in threadpool: " << ex.what() << endl;
            szErrorMessage = "std::exception"; 
        }
        catch(...)
        {
            szErrorMessage = szUnknown;
        }
        
        if(szErrorMessage)
        {
            boost::mutex::scoped_try_lock scopedLock(mxLockVector);
            if(scopedLock)
                aJobs.clear(); //kill everything else quickly
        }
        //semaphore will be posted after return
    }
    
    void jobHandler()
    {
        for(;;)
        {
            semaphore_workers.wait();
            
            if(bExiting)
                return;
            
            doJobs();
            
            semaphore_waitForAll.post();
        }
    }

    int getActualNumThreads() const  { return (int)aThreads.size(); } 
    
public:
    CThreadpool(const int nNumThreads) : semaphore_workers(0), semaphore_waitForAll(0), bExiting(false), szErrorMessage(0)
    {
        CHECK(nNumThreads < 1, "Set at least 1 thread (1 is a special case and won't actually create any threads)");
        
        if(nNumThreads > 1)
            for(int i=0;i<nNumThreads;i++)
                aThreads.create_thread(boost::bind(&CThreadpool::jobHandler, this));
    }
    
    virtual ~CThreadpool() 
    {
        bExiting=true;
        for(int i=0; i < getActualNumThreads(); i++)
            semaphore_workers.post();
            
        aThreads.join_all();
    }
    
    virtual void addJob(TNullaryFnObj & function)
    {
        boost::mutex::scoped_lock scopedLock(mxLockVector);
        aJobs.push_back(function);
    }
    
    virtual void waitForAll(TNullaryFnObj mainThreadFn)
    {
        if(getActualNumThreads() == 0 || aJobs.size() <= 1) //special case
        {
            doJobs();
            mainThreadFn();
        }
        else
        {
            const int nNumThreadsToDispatch = std::min<int>(getActualNumThreads(), (int)aJobs.size());
            
            for(int i=0;i<nNumThreadsToDispatch;i++)
                semaphore_workers.post();
                
            try
            {
                mainThreadFn();
            }
            catch(CException & ex) {
                if(ex.GetErrorMessage())
                    szErrorMessage = ex.GetErrorMessage();
                else
                    szErrorMessage = szEXITING;
            }
            catch(...) {
                szErrorMessage = szUnknown;
            }

            for(int i=0;i<nNumThreadsToDispatch;i++)
                semaphore_waitForAll.wait();
        }
        
        if(szErrorMessage)
        {
            if(szErrorMessage == szEXITING)
            {
                std::cout << "Exit triggered from a worker thread" << std::endl;
                throw CException();
            }
            std::cerr << "ERROR thrown in worker thread: " << szErrorMessage << endl;
            THROW("ERROR thrown in threadpool worker thread");
        }
    }
    
    virtual int getNumThreads() const  { return ((int)aThreads.size() < 1) ? 1 : (int)aThreads.size(); } 

};

/*Use condition variables rather than semaphores
 * 
 * Currently unreliable and no speed difference.
 *
class CThreadpool_condition : public CThreadpool_base
{
    boost::interprocess::interprocess_semaphore semaphore_waitForAll;
    boost::condition_variable cond;
    boost::mutex mxConditionVar;    
    bool bDispatch;
    
    boost::mutex mxLockVector;
    boost::thread_group aThreads;
    
    std::vector<TNullaryFnObj> aJobs;
    
    bool bExiting;
    
    //Return false if there's no jobs left
    bool onejobHandler()
    {
        //If there's a queued job then do it, else return 'done'
        TNullaryFnObj fn;
        {
            boost::mutex::scoped_lock scopedLock(mxLockVector);
            if(aJobs.size() == 0)
                return false;
            
            fn = aJobs.back();
            aJobs.pop_back();
        }
        fn();
        
        return true;
    }
    
    //Return true when we should exit
    void doJobs()
    {
        for(;;)
        {
            if(!onejobHandler())
                return;
        }
    }
    
    void jobHandler()
    {
        for(;;)
        {
            
            for(;;)
            {
                boost::unique_lock<boost::mutex> lock(mxConditionVar);
                cond.wait(lock); //occasionally this can release without being dispatched
                if(bDispatch)
                    break;
            }
            
            if(bExiting)
                return;
            
            doJobs();
            
            semaphore_waitForAll.post();
        }
    }

    int getActualNumThreads() const  { return (int)aThreads.size(); } 
    
public:
    CThreadpool_condition(const int nNumThreads) : semaphore_waitForAll(0), bExiting(false), bDispatch(false)
    {
        CHECK(nNumThreads < 1, "Set at least 1 thread (1 is a special case and won't actually create any threads)");
        
        if(nNumThreads > 1)
            for(int i=0;i<nNumThreads;i++)
                aThreads.create_thread(boost::bind(&CThreadpool_condition::jobHandler, this));
    }
    
    virtual ~CThreadpool_condition() 
    {
        bExiting=true;
        bDispatch=true;
        cond.notify_all();
            
        aThreads.join_all();
    }
    
    virtual void addJob(TNullaryFnObj & function)
    {
        aJobs.push_back(function);
    }
    
    virtual void waitForAll()
    {
        if(getActualNumThreads() == 0) //special case
            doJobs();
        else
        {
            bDispatch=true;
            cond.notify_all();
            
            //for(int i=0;i<getActualNumThreads();i++)
              //  semaphore_workers.post();

            for(int i=0;i<getActualNumThreads();i++)
                semaphore_waitForAll.wait();
                
            bDispatch=false;
        }
    }
    
    virtual int getNumThreads() const  { return ((int)aThreads.size() < 1) ? 1 : (int)aThreads.size(); } 

};*/

#define SAFECOUT(x) if(bVerbose) { \
                    boost::unique_lock<boost::mutex> lock(mxcout); \
                    cout << x << endl; \
                }

class CThreadpool_condition : public CThreadpool_base
{
    static const bool bVerbose = false;
    
    //boost::interprocess::interprocess_semaphore semaphore_waitForAll;
    boost::condition_variable cond, cond_wakeupParent;
    boost::mutex mxConditionVar, mxcout, mx_wakeupParent;
    bool bDispatch;
    int nNumDispatched; //protect by mxConditionVar
    
    //boost::mutex mxLockVector;
    boost::thread_group aThreads;
    
    std::vector<TNullaryFnObj> aJobs;
    
    bool bExiting;
    
    //Return true when we should exit
    void doJobs()
    {
        boost::unique_lock<boost::mutex> lock(mxConditionVar);
        while(!bDispatch)
        {
            SAFECOUT("Waiting for dispatch");
            cond.wait(lock); //occasionally this can release without being dispatched
            
            if(bExiting)
                return;
        }
        
        nNumDispatched++; //count from here (while mxConditionVar is locked) to make sure we only count the threads that are actually dispatched
        SAFECOUT("Dispatch nNumDispatched=" << nNumDispatched);
        
        if(aJobs.size() == 0)
        {
            SAFECOUT("No jobs on first check");
            //semaphore_waitForAll.post();
            doWakeup();
            return;
        }

        TNullaryFnObj fn = aJobs.back();
        aJobs.pop_back();
        
        lock.unlock();
        fn();
        ///...to here but then we lose the lock
        
        doRemainingJobs();
    }
    void doWakeup() 
    {
        SAFECOUT("Return after jobs finished, nNumDispatched=" << nNumDispatched);
        
        nNumDispatched--; //we have a lock on mxConditionVar protecting this BUT cond_wakeupParent is polling it
        if(!bExiting) bDispatch = false; //stop any subsequent stray dispatches from getting this far
        
        if(IS_DEBUG || bVerbose)
        {
            CHECK_P(nNumDispatched < 0, nNumDispatched, "Error counting threads");
            CHECK(mxConditionVar.try_lock(), "mxConditionVar should already be locked")
        }
        
        if(nNumDispatched==0)
        {
            {
                boost::lock_guard<boost::mutex> lock_wakeupParent(mx_wakeupParent); 
                if(!bExiting) bDispatch = false; //trying to speedup...
            }
            SAFECOUT("Return after jobs finished, nNumDispatched=" << nNumDispatched);
            cond_wakeupParent.notify_one();
        }
    }
    void doRemainingJobs()
    {
        for(;;)
        {
            //If there's a queued job then do it, else return
            boost::mutex::scoped_lock scopedLock(mxConditionVar);
            SAFECOUT(aJobs.size() << " jobs");
            if(aJobs.size() == 0)
            {
                //semaphore_waitForAll.post();
                doWakeup();
                return;
            }
            TNullaryFnObj fn = aJobs.back();
            aJobs.pop_back();

            scopedLock.unlock();
            
            fn();
        }
    }
    
    void jobHandler()
    {
        SAFECOUT("Started thread");
        
        for(;;)
        {
            SAFECOUT("Waiting for more jobs");
            doJobs();
            
            if(bExiting)
            {
                SAFECOUT("Exiting thread");
                return;
            }
        }
    }

    int getActualNumThreads() const  { return (int)aThreads.size(); } 
    
public:
    CThreadpool_condition(const int nNumThreads) : bDispatch(false), nNumDispatched(0), bExiting(false)
    {
        CHECK(nNumThreads < 1, "Set at least 1 thread (1 is a special case and won't actually create any threads)");
        
        if(nNumThreads > 1)
        {
            for(int i=0;i<nNumThreads;i++)
            {
                aThreads.create_thread(boost::bind(&CThreadpool_condition::jobHandler, this));
                SAFECOUT("Created thread " << i);
            }
        }
    }
    
    virtual ~CThreadpool_condition() 
    {
        {
            boost::unique_lock<boost::mutex> lockDispatch(mxConditionVar);
            SAFECOUT("Deleting");

            bExiting=true;
            bDispatch=true;
        }
        
        cond.notify_all();
            
        SAFECOUT("Done notify all");
        
        aThreads.join_all();
        
        SAFECOUT("Done join all");
    }
    
    virtual void addJob(TNullaryFnObj & function)
    {
        boost::lock_guard<boost::mutex> lockDispatch(mxConditionVar);
        aJobs.push_back(function);
        SAFECOUT("Added job: # jobs=" << aJobs.size());
    }
    
    virtual void waitForAll(TNullaryFnObj mainThreadFn)
    {
        mainThreadFn();
        if(getActualNumThreads() == 0 || aJobs.size() <= 1) //special case
        {
            SAFECOUT("Executing jobs in one thread. # jobs=" << aJobs.size());
            nNumDispatched=1;
            doRemainingJobs();
        }
        else
        {
            SAFECOUT("Executing jobs in all threads. # jobs=" << aJobs.size());
            {
                boost::lock_guard<boost::mutex> lockDispatch(mxConditionVar);
                bDispatch=true;
                //nNumDispatched = getActualNumThreads(); 
            }
            cond.notify_all();

            /*for(int i=0;i<getActualNumThreads();i++)
            {
                SAFECOUT("Waiting for semaphore... i=" << i);
                semaphore_waitForAll.wait();
            }   */
            
            boost::unique_lock<boost::mutex> lock_wakeupParent(mx_wakeupParent); //Different mutex because we want wakeup to be fast
            while(bDispatch /* ensure at least one has been woken up (to do the jobs) */ || nNumDispatched > 0)
            {
                SAFECOUT("Waiting for all threads to finish, nNumDispatched=" << nNumDispatched << " bDispatch=" << bDispatch);
                cond_wakeupParent.wait(lock_wakeupParent); //occasionally this can release without being dispatched
            }
            
            bDispatch=false;
            CHECK(aJobs.size() > 0, "Failed to complete jobs");
        }
    }
    
    virtual int getNumThreads() const  { return ((int)aThreads.size() < 1) ? 1 : (int)aThreads.size(); } 
};


/**
 * @class CThreadpool_RW
 * @brief Locking with notify condition variable. Jobs start as soon as they arrive. 1 thread is a special case: the functions are called by waitForAll()
 * 
 * This makes valgrind/callgrind output much nicer, and means that OpenCV showim functions will work. when MT is turned off.
 *
class CThreadpool_RW : public CThreadpool_base
{
    boost::condition_variable cond;
    boost::interprocess::interprocess_semaphore semaphore_waitForAll;
    boost::mutex mxCond;

    boost::thread_group aThreads;
    
    std::vector<TNullaryFnObj> aJobs;
    
    bool bExiting;
    
    //Return false if there's no jobs left
    bool onejobHandler()
    {
        //If there's a queued job then do it, else return 'done'
        TNullaryFnObj fn;
        {
            boost::mutex::scoped_lock scopedLock(mxLockVector);
            if(aJobs.size() == 0)
                return false;
            
            fn = aJobs.back();
            aJobs.pop_back();
        }
        fn();
        
        return true;
    }
    
    //Return true when we should exit
    void doJobs()
    {
        for(;;)
        {
            if(!onejobHandler())
                return;
        }
    }
    
    void jobHandler()
    {
        for(;;)
        {
            boost::lock<boost::mutex> lock(workerLock);
            
            if(bExiting)
                return;
            
            doJobs();
            
            //finishedLock.unlock();
            semaphore_waitForAll.post();
        }
    }

    int getActualNumThreads() const  { return (int)aThreads.size(); } 
    
public:
    CThreadpool_RW(const int nNumThreads) : semaphore_waitForAll(0), bExiting(false)
    {
        CHECK(nNumThreads < 1, "Set at least 1 thread (1 is a special case and won't actually create any threads)");

        workerLock.lock();
        
        if(nNumThreads > 1)
            for(int i=0;i<nNumThreads;i++)
                aThreads.create_thread(boost::bind(&CThreadpool_RW::jobHandler, this));
    }
    
    virtual ~CThreadpool_RW() 
    {
        bExiting=true;

        workerLock.unlock();
        aThreads.join_all();
    }
    
    virtual void addJob(TNullaryFnObj & function)
    {
        aJobs.push_back(function);
    }
    
    virtual void waitForAll()
    {
        if(getActualNumThreads() == 0) //special case
            doJobs();
        else
        {
            //finishedLock.lock();//unlocked by the *first* worker to finish
            //todo not what we want
            
            workerLock.unlock();
            
            //boost::mutex::scoped_lock lockUntilFinished(finishedLock); //force wait for jobs finished here
            for(int i=0;i<getActualNumThreads();i++)
                semaphore_waitForAll.wait();
            
            workerLock.lock();
        }
    }
    
    virtual int getNumThreads() const  { return ((int)aThreads.size() < 1) ? 1 : (int)aThreads.size(); } 

};*/


CThreadpool_base * CThreadpool_base::makeThreadpool(const int nNumThreads)
{
    if(nNumThreads>0)
        return new CThreadpool(nNumThreads);
    else 
        //return new CThreadpool_RW(-nNumThreads);
        return new CThreadpool_condition(-nNumThreads);
}   
