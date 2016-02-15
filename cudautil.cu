/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
/*#include <nvToolsExt.h>*/

#include "cudautil.h"

#include <assert.h>
#include <stdio.h>

static void *CUDA_malloc(size_t nBytes) {
    void *blk = NULL;
    cudaMallocHost( &blk, nBytes);
    return blk;
    }

static void CUDA_free(void *data) {
    cudaFreeHost(data);
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

extern "C"
void *CUDA_initialize(int iCore) {
    int nDevices;

    CUDA_CHECK(cudaGetDeviceCount,(&nDevices))
    if (nDevices == 0 ) return NULL;
    CUDA_CHECK(cudaSetDevice,(iCore % nDevices));
    CUDACTX ctx = reinterpret_cast<CUDACTX>(malloc(sizeof(struct cuda_ctx)));
    assert(ctx!=NULL);
    ctx->wqCuda = ctx->wqFree = NULL;
    ctx->nodePP = NULL;
    ctx->nodePC = NULL;

    CUDA_CHECK(cudaGetDeviceProperties,(&ctx->prop,iCore % nDevices));

    return ctx;
    }

extern "C"
int CUDA_flushDone(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work, **last = &cuda->wqCuda;
    while( (work=*last) !=NULL ) {
        cudaError_t rc = cudaEventQuery(work->event);
        if (rc==cudaSuccess) {
            if ( (*work->checkFcn)(work->ctx,work) == 0) {
                *last = work->next;
                work->next = cuda->wqFree;
                cuda->wqFree = work;
                continue;
                }
            }
        else if (rc!=cudaErrorNotReady) {
            fprintf(stderr,"cudaEventQuery error %d: %s\n", rc, cudaGetErrorString(rc));
            exit(1);
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
    if ( (*initWork)(ctx,work) ) {
	    work->ctx = ctx;
	    work->checkFcn = checkWork;
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
    int hostBufSize = inCudaBufSize > outCudaBufSize ? inCudaBufSize : outCudaBufSize;
    cuda->inCudaBufSize = inCudaBufSize;
    cuda->outCudaBufSize = outCudaBufSize;
    while(cudaSize--) {
        work = reinterpret_cast<CUDAwqNode *>(malloc(sizeof(CUDAwqNode)));
        work->pHostBuf = CUDA_malloc(hostBufSize);
        assert(work->pHostBuf!=NULL);
        work->pCudaBufIn = CUDA_gpu_malloc(inCudaBufSize);
        assert(work->pCudaBufIn!=NULL);
        work->pCudaBufOut = CUDA_gpu_malloc(outCudaBufSize);
        assert(work->pCudaBufOut!=NULL);
        CUDA_CHECK(cudaEventCreateWithFlags,( &work->event, cudaEventDisableTiming ));
        CUDA_CHECK(cudaStreamCreate,( &work->stream ));
        work->ctx = NULL;
        work->checkFcn = NULL;
        work->next = cuda->wqFree;
        cuda->wqFree = work;
        }
    }

extern "C"
void CUDA_finish(void *vcuda) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;

    while( (work=cuda->wqFree) != NULL ) {
        cuda->wqFree = work->next;
        CUDA_free(work->pHostBuf);
        CUDA_gpu_free(work->pCudaBufIn);
        CUDA_gpu_free(work->pCudaBufOut);
        free(work);
        }
    cudaDeviceReset();
    free(cuda);
    }
