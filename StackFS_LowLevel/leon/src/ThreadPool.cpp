
//Reference : https://13abyknight.wordpress.com/2013/03/20/a-simple-thread-pool-c-implementation-on-linux/
#include <iostream>
#include <queue>
#include <pthread.h>
#include "ThreadPool.h"
 
using namespace std;
  
int job::finished_jobs = 0;
 
pthread_mutex_t job::jobLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t thread_pool::jobQueue_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t thread_pool::jobQueue_cond = PTHREAD_COND_INITIALIZER;
 
void thread_pool::initThreads(pthread_t *threads)
{
     
    for(int i = 0; i < numOfThreads; i++)
    {
        pthread_create(&threads[i], NULL, &thread_pool::threadExecute, (void *)this);
        cout << "Thread:" << i << " is alive now!\n";
    }
}
 
void thread_pool::assignJob(job* _job_)
{
    pthread_mutex_lock(&jobQueue_lock);
    jobQueue.push(_job_);
    pthread_mutex_unlock(&jobQueue_lock);
    pthread_cond_signal(&jobQueue_cond);
}
 
bool thread_pool::loadJob(job*& _job_)
{
    pthread_mutex_lock(&jobQueue_lock);
    while(jobQueue.empty())
        pthread_cond_wait(&jobQueue_cond, &jobQueue_lock);
    _job_ = jobQueue.front();
    jobQueue.pop();
    pthread_mutex_unlock(&jobQueue_lock);
    return true;
}
 
void *thread_pool::threadExecute(void *param)
{
    thread_pool *p = (thread_pool *)param;
    job *oneJob = NULL;
    while(p->loadJob(oneJob))
    {
        if(oneJob)
            oneJob->working();
        delete oneJob;
        oneJob = NULL;
    }
}

