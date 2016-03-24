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

static void *CUDA_malloc(size_t nBytes) {
    void *blk = valloc(nBytes);
//    void *blk = NULL;
//    cudaMallocHost( &blk, nBytes);
    return blk;
    }

static void CUDA_free(void *data) {
//    cudaFreeHost(data);
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

// This is cheeserific!
static int recovery_epoch = 0;
static pthread_mutex_t recovery_mutex = PTHREAD_MUTEX_INITIALIZER;

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
    CUDA_CHECK(cudaEventCreateWithFlags,( &work->event, cudaEventDisableTiming ));
    cudaEventQuery( work->event );
    CUDA_CHECK(cudaEventQuery,( work->event ));
    CUDA_CHECK(cudaGetLastError,());
    CUDA_CHECK(cudaEventQuery,( work->event ));
    return work;
    }


static void setup_cuda(CUDACTX cuda) {
    CUDAwqNode *work;
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
        }

    }

extern "C"
void *CUDA_initialize(int iCore) {
    int nDevices;

    CUDA_CHECK(cudaGetDeviceCount,(&nDevices))
    if (nDevices == 0 ) return NULL;
    CUDA_CHECK(cudaSetDevice,(iCore % nDevices));
    CUDACTX ctx = reinterpret_cast<CUDACTX>(malloc(sizeof(struct cuda_ctx)));
    assert(ctx!=NULL);
    ctx->epoch = 0;
    ctx->nWorkQueueSize = 0;
    ctx->nWorkQueueBusy = 0;
    ctx->wqCuda = ctx->wqFree = NULL;
    ctx->nodePP = NULL;
    ctx->nodePC = NULL;
    ctx->ewIn = NULL;
    ctx->ewt = NULL;

#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    if (gethostname(ctx->hostname,MAXHOSTNAMELEN))
#endif
        strcpy(ctx->hostname,"unknown");


    CUDA_CHECK(cudaGetDeviceProperties,(&ctx->prop,iCore % nDevices));

    return ctx;
    }

void cuda_setup_ewald(CUDACTX cuda);

void CUDA_attempt_recovery(CUDACTX cuda) {
    cudaDeviceSynchronize();
    pthread_mutex_lock(&recovery_mutex);
    // Looks like it's all up to me!
    if (cuda->epoch == recovery_epoch) {
        printf("%s: ************* CUDA attempting recovery\n", cuda->hostname);
        ++recovery_epoch;
        sleep(1);
        cudaDeviceSynchronize();
        printf("Reset!\n");
        cudaDeviceReset();
        cuda_setup_ewald(cuda);
        sleep(1);
        }
    cudaDeviceSynchronize();
    cudaGetLastError();
    CUDA_CHECK(cudaGetDeviceProperties,(&cuda->prop,0));
    cuda->epoch = recovery_epoch;
    pthread_mutex_unlock(&recovery_mutex);
    setup_cuda(cuda);
    }

extern "C"
int CUDA_flushDone(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work, **last = &cuda->wqCuda;
    if (cuda->epoch != recovery_epoch) {
        CUDA_attempt_recovery(cuda);
        }
    while( (work=*last) !=NULL ) {
        cudaError_t rc = cudaEventQuery(work->event);
// To test the "reset"
//        if (time(NULL)%5 == 1 ) CUDA_attempt_recovery(cuda);
//        else
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
            fprintf(stderr,"%s: cudaEventQuery error %d: %s\n", cuda->hostname, rc, cudaGetErrorString(rc));
            exit(1);
//            CUDA_attempt_recovery(cuda);
            }
        else if (work->startTime != 0) {
            time_t now;
            double seconds;
            time(&now);
            seconds = difftime(now,work->startTime);
            if (seconds>=30) {
                fprintf(stderr,"%s: cudaEventQuery has returned cudaErrorNotReady for %f seconds\n",cuda->hostname,seconds);
                work->startTime = 0;
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
    if ( (*initWork)(work->ctx,work) ) {
	    time(&work->startTime);
	    work->next = cuda->wqCuda;
	    cuda->wqCuda = work;
	    }
    else {
        assert(0);
        }
    return 1;
    }

extern "C"
void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize, int outCudaBufSize) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;
    cuda->inCudaBufSize = inCudaBufSize;
    cuda->outCudaBufSize = outCudaBufSize;
    cuda->nWorkQueueSize = cudaSize;
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
