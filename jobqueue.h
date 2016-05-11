/**
 * This file implements a multithreaded job queue. Jobs are pushed
 * onto a queue by the main program. Each thread (or worker) removes a
 * job from the queue, executes it, and then goes back for
 * another. When all jobs are finished, control returns to the main
 * function. 
 *
 * Alan R. Rogers 2013-6-28
 */

/*
 * Copyright (c) 2013, Alan R. Rogers
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef LDPSIZ_JOBQUEUE
#define LDPSIZ_JOBQUEUE

typedef struct JobQueue JobQueue;

JobQueue   *JobQueue_new(int nthreads);
void        JobQueue_addJob(JobQueue * jq, int (*jobfun) (void *),
                            void *param);
void        JobQueue_addTask(JobQueue * jq, int (*taskfun) (void *),
                             void *param);
void        JobQueue_waitOnJobs(JobQueue * jq);
void        JobQueue_free(JobQueue * jq);
#endif