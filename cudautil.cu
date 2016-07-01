/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cuda.h>

/*#include <nvToolsExt.h>*/
#ifdef USE_NVML
#include <nvidia/gdk/nvml.h>
#endif
#include <signal.h>

#include "cudautil.h"

#include <assert.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
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
#ifdef _MSC_VER
    blk = _aligned_malloc(nBytes, nPageSize);
#else
    if (posix_memalign(&blk, nPageSize, nBytes)) blk = NULL;
#endif
    char *p = reinterpret_cast<char *>(blk);
    char *e = p + nBytes;
    for(; p<e; p+= nPageSize) *p = 0;

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

CUDAwqNode *getNode(CUDACTX cuda) {
    CUDAwqNode *work;
    CUDA_flushDone(cuda);

    if (OPA_Queue_is_empty(&cuda->wqFree)) return NULL;
    OPA_Queue_dequeue(&cuda->wqFree, work, CUDAwqNode, q.hdr);
    work->pwqDone = &cuda->wqDone;

    ++cuda->nWorkQueueBusy;
    work->ctx = cuda;
    work->initFcn = NULL;
    work->doneFcn = NULL;
    work->kernelName = "unknown";

    work->ppSizeIn = 0;
    work->ppSizeOut = 0;
    work->ppnBuffered = 0;
    work->ppnBlocks = 0;
    return work;
    }

// This is cheeserific!
static int recovery_epoch = 0;
static pthread_mutex_t recovery_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_barrier_t recovery_barrier;

/*
** This function allocates buffers, and can restart outstanding requests after a device reset.
*/
static CUDAwqNode *setup_node(CUDACTX cuda,CUDAwqNode *work) {
    int hostBufSize = work->inBufferSize > work->outBufferSize ? work->inBufferSize : work->outBufferSize;
    CUDA_CHECK(cudaHostRegister,(work->pHostBufFromGPU, hostBufSize, cudaHostRegisterPortable));
    CUDA_CHECK(cudaHostRegister,(work->pHostBufToGPU,   hostBufSize, cudaHostRegisterPortable));
    work->pCudaBufIn = CUDA_gpu_malloc(work->inBufferSize);
    assert(work->pCudaBufIn!=NULL);
    work->pCudaBufOut = CUDA_gpu_malloc(work->outBufferSize);
    assert(work->pCudaBufOut!=NULL);
#ifdef USE_SINGLE_STREAM
    work->stream = cuda->stream;
    work->eventCopyDone = cuda->eventCopyDone;
    work->eventKernelDone = cuda->eventKernelDone;
#else
    CUDA_CHECK(cudaStreamCreate,( &work->stream ));
#ifdef USE_CUDA_EVENTS
    CUDA_CHECK(cudaEventCreateWithFlags,( &work->event, cudaEventDisableTiming ));
#endif
    CUDA_CHECK(cudaEventCreateWithFlags,( &work->eventCopyDone, cudaEventDisableTiming ));
    CUDA_CHECK(cudaEventCreateWithFlags,( &work->eventKernelDone, cudaEventDisableTiming ));
#endif
    return work;
    }

#if 0
static void setup_list(CUDACTX cuda,OPA_Queue_info_t *q) {
    CUDAwqNode *work;
    OPA_Queue_info_t t;
    OPA_Queue_init(&t);
    while (!OPA_Queue_is_empty(q)) {
        OPA_Queue_dequeue(q, work, CUDAwqNode, q.hdr);
        OPA_Queue_enqueue(&t, work, CUDAwqNode, q.hdr);
        }
    while (!OPA_Queue_is_empty(&t)) {
        OPA_Queue_dequeue(&t, work, CUDAwqNode, q.hdr);
        setup_node(cuda,work);
        OPA_Queue_enqueue(q, work, CUDAwqNode, q.hdr);
        }
    }

static int setup_cuda(CUDACTX cuda) {
    // Free element -- just allocate buffers and create the stream and event
    setup_list(cuda,&cuda->wqFree);

    if (cuda->nodePP) setup_node(cuda,cuda->nodePP);
    if (cuda->nodePC) setup_node(cuda,cuda->nodePC);
    if (cuda->nodeEwald) setup_node(cuda,cuda->nodeEwald);

    // During a restart we have to restart any active requests as well.
    // At normal startup this work queue will be empty.
    CUDAwqNode *work;
    int nQueued = 0;
    for(work=cuda->wqCudaBusy; work!=NULL; work=work->q.next) {
        setup_node(cuda,work);
        (*work->initFcn)(work->ctx,work); // This restarts the work
        work->startTime = CUDA_getTime();
        ++nQueued;
        }
    return nQueued;
    }
#endif

static void setup_ewald(CUDACTX cuda) {
//    CUDA_CHECK(cudaStreamCreate,( &cuda->streamEwald ));
#ifdef USE_CUDA_EVENTS
    CUDA_CHECK(cudaEventCreateWithFlags,( &cuda->eventEwald, cudaEventDisableTiming ));
#endif
    }

extern "C"
void *CUDA_initialize(int nCores, int iCore, OPA_Queue_info_t *queueWORK, OPA_Queue_info_t *queueREGISTER) {
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
    ctx->nKernelLaunches = 0;

#ifdef USE_SINGLE_STREAM
    CUDA_CHECK(cudaStreamCreate,( &ctx->stream ));
    CUDA_CHECK(cudaEventCreateWithFlags,( &ctx->eventCopyDone, cudaEventDisableTiming ));
    CUDA_CHECK(cudaEventCreateWithFlags,( &ctx->eventKernelDone, cudaEventDisableTiming ));
#endif
    ctx->wqCudaBusy = NULL;
    ctx->nwqCudaBusy = 0;
    OPA_Queue_init(&ctx->wqFree);
    OPA_Queue_init(&ctx->wqDone);
    ctx->queueWORK = queueWORK;
    ctx->queueREGISTER = queueREGISTER;
    ctx->nodePP = NULL;
    ctx->nodePC = NULL;
    ctx->nodeEwald = NULL;
    ctx->ewIn = NULL;
    ctx->ewt = NULL;

    if (ctx->iCore == 0) pthread_barrier_init(&recovery_barrier,NULL,ctx->nCores);
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    if (gethostname(ctx->hostname,MAXHOSTNAMELEN))
#endif
        strcpy(ctx->hostname,"unknown");
    CUDA_CHECK(cudaGetDeviceProperties,(&ctx->prop,iCore % nDevices));

#ifdef USE_NVML
    nvmlReturn_t result;
    result = nvmlInit();
    if (result != NVML_SUCCESS) printf("Failed to initialize NVML: %s\n", nvmlErrorString(result));
    assert(result == NVML_SUCCESS);
#endif

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
    cuda->epoch = recovery_epoch;

#ifdef ATTEMPT_RESET
    // One thread will perform the device reset
    rc = pthread_barrier_wait(&recovery_barrier);
    if (rc==PTHREAD_BARRIER_SERIAL_THREAD) {
        printf("%2d: calling cudaDeviceReset()\n", cuda->iCore);
        cudaDeviceReset();
        sleep(1);
	    }


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
#else
    raise(SIGINT); /* If the debugger is running we might get a traceback. */
#endif
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

static void CUDA_finishWork(CUDACTX cuda,CUDAwqNode *work) {
    int rc = (*work->doneFcn)(work->ctx,work);
    assert(rc == 0);
    work->doneFcn = NULL; /* Make sure we don't do this multiple times! */
    OPA_Queue_enqueue(&cuda->wqFree, work, CUDAwqNode, q.hdr);
    --cuda->nWorkQueueBusy;
    }

extern "C" void CUDA_startWork(void *vcuda,OPA_Queue_info_t *queueWORK) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    assert( cuda->nwqCudaBusy==0 || cuda->wqCudaBusy!=NULL);
    while (cuda->nwqCudaBusy < 3 && !OPA_Queue_is_empty(queueWORK)) {
        CUDAwqNode *node;
        OPA_Queue_dequeue(queueWORK, node, CUDAwqNode, q.hdr);
        node->startTime = CUDA_getTime();
        node->q.next = cuda->wqCudaBusy;
        cuda->wqCudaBusy = node;
        ++cuda->nwqCudaBusy;
        cudaError_t rc = static_cast<cudaError_t>((*node->initFcn)(node->ctx,node));
        ++cuda->nKernelLaunches;
        if ( rc != cudaSuccess) CUDA_attempt_recovery(cuda,rc);
        }
    }

extern "C" void CUDA_registerBuffers(void *vcuda, OPA_Queue_info_t *queueWORK) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    while (!OPA_Queue_is_empty(queueWORK)) {
	CUDAwqNode *node;
	OPA_Queue_dequeue(queueWORK, node, CUDAwqNode, q.hdr);
	setup_node(cuda, node);
	OPA_Queue_enqueue(node->pwqDone, node, CUDAwqNode, q.hdr);
    }
}

extern "C"
int CUDA_flushDone(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work, **last;
    // Someone else has started the recovery process -- we just participate
    if (cuda->epoch != recovery_epoch) {
        CUDA_attempt_recovery(cuda,cudaSuccess);
        }
    last = &cuda->wqCudaBusy;
    while( (work = *last) != NULL) {
#ifdef USE_CUDA_EVENTS
        cudaError_t rc = cudaEventQuery(work->event);
#else
        cudaError_t rc = cudaStreamQuery(work->stream);
#endif
        if (rc==cudaSuccess) {
            assert(work->doneFcn != NULL); /* Only one cudaSuccess per customer! */
            *last = work->q.next;
            --cuda->nwqCudaBusy;
            OPA_Queue_enqueue(work->pwqDone, work, CUDAwqNode, q.hdr);
            continue;
            }
        else if (rc!=cudaErrorNotReady) {
            fprintf(stderr,"%s: dimBlock=%d,%d,%d  dimGrid=%d,%d,%d\n",
                cuda->hostname, work->dimBlock.x,work->dimBlock.y,work->dimBlock.z,work->dimGrid.x,work->dimGrid.y,work->dimGrid.z);
            CUDA_attempt_recovery(cuda,rc);
            }
        else if (work->startTime != 0) {
            double seconds = CUDA_getTime() - work->startTime;
            if (seconds>=2.0) {
                rc = cudaEventQuery(work->eventCopyDone);
                const char *done1 = (rc==cudaSuccess?"yes":"no");
                rc = cudaEventQuery(work->eventKernelDone);
                const char *done2 = (rc==cudaSuccess?"yes":"no");
                fprintf(stderr,"%s: cudaStreamQuery for kernel %s has returned cudaErrorNotReady for %f seconds, Copy=%s Kernel=%s Copy=no\n",
                    cuda->hostname, work->kernelName, seconds, done1, done2);
                fprintf(stderr,"%s: dimBlock=%d,%d,%d  dimGrid=%d,%d,%d Launches=%llu\n",
                    cuda->hostname, work->dimBlock.x,work->dimBlock.y,work->dimBlock.z,
                    work->dimGrid.x,work->dimGrid.y,work->dimGrid.z,cuda->nKernelLaunches);
                if (work->dumpFcn) (*work->dumpFcn)(work);
                work->startTime = 0;
                CUDA_attempt_recovery(cuda,cudaErrorLaunchTimeout);
                break;
                }
            }
        last = &work->q.next;
        }
	while (!OPA_Queue_is_empty(&cuda->wqDone)) {
	    OPA_Queue_dequeue(&cuda->wqDone, work, CUDAwqNode, q.hdr);
        CUDA_finishWork(cuda,work);
        }
    return cuda->nWorkQueueBusy > 0;
    }

extern "C"
int CUDA_queue(void *vcuda, void *ctx,
    int (*initWork)(void *ctx,void *work),
    int (*checkWork)(void *ctx,void *work)) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;
    CUDA_flushDone(vcuda);

    if (OPA_Queue_is_empty(&cuda->wqFree) || initWork==NULL) return 0;
	OPA_Queue_dequeue(&cuda->wqFree, work, CUDAwqNode, q.hdr);

    ++cuda->nWorkQueueBusy;
    work->ctx = ctx;
    work->doneFcn = checkWork;
    work->initFcn = initWork;
    work->dumpFcn = NULL;
    work->startTime = CUDA_getTime();
    work->q.next = cuda->wqCudaBusy;
    cuda->wqCudaBusy = work;
    ++cuda->nwqCudaBusy;

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
        work->inBufferSize = cuda->inCudaBufSize;
        work->outBufferSize = cuda->outCudaBufSize;
        work->pwqDone = &cuda->wqDone;
        work->pHostBufFromGPU = CUDA_malloc(hostBufSize);
        assert(work->pHostBufFromGPU!=NULL);
        work->pHostBufToGPU = CUDA_malloc(hostBufSize);
        assert(work->pHostBufToGPU!=NULL);
        work->ctx = NULL;
        work->initFcn = NULL;
        work->doneFcn = NULL;
        work->startTime = 0;
//        OPA_Queue_enqueue(&cuda->wqFree, work, CUDAwqNode, q.hdr);
        work->pwqDone = &cuda->wqFree;
        OPA_Queue_enqueue(cuda->queueREGISTER, work, CUDAwqNode, q.hdr);
        }
    cuda->nWorkQueueBusy = 0;
//    setup_cuda(cuda);
    setup_ewald(cuda);
    }

extern "C"
void CUDA_finish(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;
    while (!OPA_Queue_is_empty(&cuda->wqFree)) {
        OPA_Queue_dequeue(&cuda->wqFree, work, CUDAwqNode, q.hdr);
        CUDA_free(work->pHostBufFromGPU);
        CUDA_free(work->pHostBufToGPU);
        CUDA_gpu_free(work->pCudaBufIn);
        CUDA_gpu_free(work->pCudaBufOut);
        free(work);
        }
    cudaDeviceReset();
    free(cuda);
    }
