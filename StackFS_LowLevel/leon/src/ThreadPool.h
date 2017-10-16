//Reference: https://13abyknight.wordpress.com/2013/03/20/a-simple-thread-pool-c-implementation-on-linux/
#ifndef __THREADPOOL_H_
#define __THREADPOOL_H_

#include <iostream>
#include <queue>
#include <pthread.h>

using namespace std;

class job
{
public:
    job(int id) : jobID(id) {}
    virtual ~job(){}
    //will be overrided by specific JOB
    void virtual working() = 0;
    static int finished_jobs;
    static pthread_mutex_t jobLock;
    int jobID;
};
 
class thread_pool
{
public:
    static pthread_mutex_t jobQueue_lock;
    static pthread_cond_t jobQueue_cond;
    thread_pool(){ thread_pool(2); }
    thread_pool(int num) : numOfThreads(num) {}
    virtual ~thread_pool() { while(!jobQueue.empty()) jobQueue.pop(); };
    void initThreads(pthread_t *);
    void assignJob(job *_job_);
    bool loadJob(job *&_job_);
    static void *threadExecute(void *);
private:
    queue<job*> jobQueue;
    int numOfThreads;
};

#endif