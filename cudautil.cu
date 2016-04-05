/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
/*#include <nvToolsExt.h>*/

#include "cudautil.h"

#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>
#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h> /* for MAXHOSTNAMELEN, if available */
#endif
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

static void *CUDA_malloc(size_t nBytes) {
#ifdef __linux__
    uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
    uint64_t nPageSize = 512;
#endif
    void *blk;
    if (posix_memalign(&blk,nPageSize,nBytes)) blk = NULL;
    return blk;
    }

static void CUDA_free(void *data) {
    free(data);
    }

static void *CUDA_gpu_malloc(size_t nBytes) {
    void *blk = NULL;
    CUDA_CHECK(cudaMalloc,(&blk, nBytes));
    return blk;
    }

static void CUDA_gpu_free(void *blk) {
    CUDA_CHECK(cudaFree,(blk));
    }

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line) {
    fprintf(stderr,"%s error %d in %s(%d)\n%s\n", fname, rc, file, line, cudaGetErrorString(rc));
    exit(1);
    }

#ifdef _MSC_VER
double CUDA_getTime() {
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    return clock / 10000000.0;
    }
#else
double CUDA_getTime() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
    }
#endif

// This is cheeserific!
static int recovery_epoch = 0;
static pthread_mutex_t recovery_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_barrier_t recovery_barrier;

/*
** This function allocates buffers, and can restart outstanding requests after a device reset.
*/
static CUDAwqNode *setup_node(CUDACTX cuda,CUDAwqNode *work) {
    int hostBufSize = cuda->inCudaBufSize > cuda->outCudaBufSize ? cuda->inCudaBufSize : cuda->outCudaBufSize;
    CUDA_CHECK(cudaHostRegister,(work->pHostBufFromGPU, hostBufSize, cudaHostRegisterPortable));
    CUDA_CHECK(cudaHostRegister,(work->pHostBufToGPU,   hostBufSize, cudaHostRegisterPortable));
    work->pCudaBufIn = CUDA_gpu_malloc(cuda->inCudaBufSize);
    assert(work->pCudaBufIn!=NULL);
    work->pCudaBufOut = CUDA_gpu_malloc(cuda->outCudaBufSize);
    assert(work->pCudaBufOut!=NULL);
    CUDA_CHECK(cudaStreamCreate,( &work->stream ));
#ifdef USE_CUDA_EVENTS
    CUDA_CHECK(cudaEventCreateWithFlags,( &work->event, cudaEventDisableTiming ));
#endif
    return work;
    }

static int setup_cuda(CUDACTX cuda) {
    CUDAwqNode *work;
    int nQueued = 0;

    // Free element -- just allocate buffers and create the stream and event
    for(work=cuda->wqFree; work!=NULL; work=work->next) {
        setup_node(cuda,work);
        }

    // During a restart we have to restart any active requests as well.
    // At normal startup this work queue will be empty.
    for(work=cuda->wqCuda; work!=NULL; work=work->next) {
        setup_node(cuda,work);
        }

    if (cuda->nodePP) {
        work = setup_node(cuda,cuda->nodePP);
        }
    if (cuda->nodePC) {
        work = setup_node(cuda,cuda->nodePC);
        }

    for(work=cuda->wqCuda; work!=NULL; work=work->next) {
        (*work->initFcn)(work->ctx,work); // This restarts the work
	    work->startTime = CUDA_getTime();
        ++nQueued;
        }

    return nQueued;
    }

static void setup_ewald(CUDACTX cuda) {
    CUDA_CHECK(cudaStreamCreate,( &cuda->streamEwald ));
#ifdef USE_CUDA_EVENTS
    CUDA_CHECK(cudaEventCreateWithFlags,( &cuda->eventEwald, cudaEventDisableTiming ));
#endif
    }

extern "C"
void *CUDA_initialize(int nCores,int iCore) {
    int nDevices;

    CUDA_CHECK(cudaGetDeviceCount,(&nDevices))
    if (nDevices == 0 ) return NULL;
    CUDA_CHECK(cudaSetDevice,(iCore % nDevices));
    CUDACTX ctx = reinterpret_cast<CUDACTX>(malloc(sizeof(struct cuda_ctx)));
    assert(ctx!=NULL);
    ctx->nCores = nCores;
    ctx->iCore = iCore;
    ctx->epoch = 0;
    ctx->nWorkQueueSize = 0;
    ctx->nWorkQueueBusy = 0;
    ctx->wqCuda = ctx->wqFree = NULL;
    ctx->nodePP = NULL;
    ctx->nodePC = NULL;
    ctx->ewIn = NULL;
    ctx->ewt = NULL;

    if (ctx->iCore == 0) pthread_barrier_init(&recovery_barrier,NULL,ctx->nCores);
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    if (gethostname(ctx->hostname,MAXHOSTNAMELEN))
#endif
        strcpy(ctx->hostname,"unknown");
    CUDA_CHECK(cudaGetDeviceProperties,(&ctx->prop,iCore % nDevices));

    return ctx;
    }

cudaError_t cuda_setup_ewald(CUDACTX cuda);

void CUDA_attempt_recovery(CUDACTX cuda,cudaError_t errorCode) {
    int rc;
    // Decide which thread will initiate the recovery
    pthread_mutex_lock(&recovery_mutex);
    if (cuda->epoch == recovery_epoch) { // Looks like it's all up to me!
        int nEquals = (80 - 21 - strlen(cuda->hostname)) / 2;
        char sEquals[41];
        if (nEquals < 10) nEquals = 10;
        sEquals[nEquals] = 0;
        while(nEquals) sEquals[--nEquals] = '=';
        printf("%s CUDA ERROR ON NODE %s %s\n", sEquals, cuda->hostname, sEquals);
        printf("%2d: cuda error %d: %s\n", cuda->iCore, errorCode, cudaGetErrorString(errorCode));
        time_t now;
        time(&now);
        struct tm* tm_info = localtime(&now);
        strftime(sEquals, 26, "%Y:%m:%d %H:%M:%S", tm_info);
        puts(sEquals);
        ++recovery_epoch;
        }
    pthread_mutex_unlock(&recovery_mutex);

    // One thread will perform the device reset
    rc = pthread_barrier_wait(&recovery_barrier);
    if (rc==PTHREAD_BARRIER_SERIAL_THREAD) {
        printf("%2d: calling cudaDeviceReset()\n", cuda->iCore);
        cudaDeviceReset();
        sleep(1);
	    }

    cuda->epoch = recovery_epoch;

    pthread_barrier_wait(&recovery_barrier);
    setup_ewald(cuda); // Create the necessary event/stream for Ewald setup

    rc = pthread_barrier_wait(&recovery_barrier);
    if (rc==PTHREAD_BARRIER_SERIAL_THREAD) {
        CUDA_CHECK(cuda_setup_ewald,(cuda)); // This MUST work now
        }
    // Now all threads must restart any work in progress
    rc = pthread_barrier_wait(&recovery_barrier);
    int nRestart = setup_cuda(cuda);
    pthread_mutex_lock(&recovery_mutex);
    printf("%2d: restarted %d request%c\n", cuda->iCore, nRestart, nRestart==1 ? ' ':'s');
    pthread_mutex_unlock(&recovery_mutex);

    rc = pthread_barrier_wait(&recovery_barrier);
    if (rc==PTHREAD_BARRIER_SERIAL_THREAD) {
        printf(
            "========================================"
            "========================================\n");
        }
    }

extern "C"
void CUDA_checkForRecovery(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    if (cuda->epoch != recovery_epoch) {
        CUDA_attempt_recovery(cuda,cudaSuccess);
        }
    }

extern "C"
int CUDA_flushDone(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work, **last = &cuda->wqCuda;
    // Someone else has started the recovery process -- we just participate
    if (cuda->epoch != recovery_epoch) {
        CUDA_attempt_recovery(cuda,cudaSuccess);
        }
    while( (work=*last) !=NULL ) {
#ifdef USE_CUDA_EVENTS
        cudaError_t rc = cudaEventQuery(work->event);
#else
        cudaError_t rc = cudaStreamQuery(work->stream);
#endif
        if (rc==cudaSuccess) {
            assert(work->checkFcn != NULL); /* Only one cudaSuccess per customer! */
            int rc = (*work->checkFcn)(work->ctx,work);
            assert(rc == 0);
            work->checkFcn = NULL; /* Make sure we don't do this multiple times! */
            *last = work->next;
            work->next = cuda->wqFree;
            cuda->wqFree = work;
            --cuda->nWorkQueueBusy;
            continue;
            }
        else if (rc!=cudaErrorNotReady) {
            CUDA_attempt_recovery(cuda,rc);
            }
        else if (work->startTime != 0) {
            double seconds = CUDA_getTime() - work->startTime;
            if (seconds>=1.0) {
                fprintf(stderr,"%s: cudaXxxxxQuery has returned cudaErrorNotReady for %f seconds\n",cuda->hostname,seconds);
                work->startTime = 0;
                CUDA_attempt_recovery(cuda,cudaErrorLaunchTimeout);
                break;
                }
            }
        last = &work->next;
        }
    return cuda->wqCuda != NULL;
    }

extern "C"
int CUDA_queue(void *vcuda, void *ctx,
    int (*initWork)(void *ctx,void *work),
    int (*checkWork)(void *ctx,void *work),
    int (*doneWork)(void *ctx)) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;
    CUDA_flushDone(vcuda);
    if (cuda->wqFree == NULL || initWork==NULL) return 0;
    work = cuda->wqFree;
    cuda->wqFree = work->next;
    ++cuda->nWorkQueueBusy;
    work->ctx = ctx;
    work->checkFcn = checkWork;
    work->initFcn = initWork;
    work->next = cuda->wqCuda;
    cuda->wqCuda = work;
    work->startTime = CUDA_getTime();
    cudaError_t rc = static_cast<cudaError_t>((*work->initFcn)(work->ctx,work));
    if ( rc != cudaSuccess) CUDA_attempt_recovery(cuda,rc);
    return 1;
    }

extern "C"
void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize, int outCudaBufSize) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;
    cuda->inCudaBufSize = inCudaBufSize;
    cuda->outCudaBufSize = outCudaBufSize;
    cuda->nWorkQueueSize = cudaSize < 6 ? 6 : cudaSize;
    int hostBufSize = cuda->inCudaBufSize > cuda->outCudaBufSize ? cuda->inCudaBufSize : cuda->outCudaBufSize;
    while(cudaSize--) {
        work = reinterpret_cast<CUDAwqNode *>(malloc(sizeof(CUDAwqNode)));
        work->pHostBufFromGPU = CUDA_malloc(hostBufSize);
        assert(work->pHostBufFromGPU!=NULL);
        work->pHostBufToGPU = CUDA_malloc(hostBufSize);
        assert(work->pHostBufToGPU!=NULL);
        work->ctx = NULL;
        work->checkFcn = NULL;
        work->startTime = 0;
        work->next = cuda->wqFree;
        cuda->wqFree = work;
        }
    cuda->nWorkQueueBusy = 0;
    setup_cuda(cuda);
    setup_ewald(cuda);
    }

extern "C"
void CUDA_finish(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;

    while( (work=cuda->wqFree) != NULL ) {
        cuda->wqFree = work->next;
        CUDA_free(work->pHostBufFromGPU);
        CUDA_free(work->pHostBufToGPU);
        CUDA_gpu_free(work->pCudaBufIn);
        CUDA_gpu_free(work->pCudaBufOut);
        free(work);
        }
    cudaDeviceReset();
    free(cuda);
    }
