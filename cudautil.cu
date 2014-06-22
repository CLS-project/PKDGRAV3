/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cudautil.h"

#include <assert.h>
#include <stdio.h>

extern "C"
void *CUDA_malloc(size_t nBytes) {
    void *blk = NULL;
    cudaHostAlloc( &blk, nBytes, 0 );
    return blk;
    }

extern "C"
void CUDA_free(void *data) {
    cudaFreeHost(data);
    }

extern "C"
void *CUDA_gpu_malloc(size_t nBytes) {
    void *blk = NULL;
    CUDA_CHECK(cudaMalloc,(&blk, nBytes));
    return blk;
    }

extern "C"
void CUDA_gpu_free(void *blk) {
    CUDA_CHECK(cudaFree,(blk));
    }

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line) {
    fprintf(stderr,"%s error %d in %s(%d)\n%s\n", fname, rc, file, line, cudaGetErrorString(rc));
    exit(1);
    }

extern "C" void CUDAsetupPP(void);
extern "C" void CUDAsetupPC(void);
#if 0

extern "C"
void *CUDA_initialize(int iProc,int iWorkQueueSize, size_t tileSize, size_t ParticlesSize, size_t OutSize) {
    int nDevices, i;

    CUDA_CHECK(cudaGetDeviceCount,(&nDevices))
    if (nDevices == 0 ) return NULL;
    CUDA_CHECK(cudaSetDevice,(iProc % nDevices));
    /*CUDA_CHECK(cudaSetDevice,(iProc%1));*/



    CUDACTX ctx = reinterpret_cast<CUDACTX>(malloc(sizeof(struct cuda_ctx)));
    assert(ctx!=NULL);
    ctx->iWorkQueueSize = iWorkQueueSize;
    ctx->in = NULL;
    ctx->block = NULL;

    CUDA_CHECK(cudaGetDeviceProperties,(&ctx->prop,iProc % nDevices));

    for(i=0; i<iWorkQueueSize; ++i) {

        gpuInput *in = reinterpret_cast<gpuInput *>(malloc(sizeof(gpuInput)));
        assert(in!=NULL);
        in->next = ctx->in;
        ctx->in = in;
        CUDA_CHECK(cudaMalloc,(&in->in, ParticlesSize));
        CUDA_CHECK(cudaHostAlloc,(reinterpret_cast<void **>(&in->cpuIn), ParticlesSize, 0));
        CUDA_CHECK(cudaEventCreate,( &in->event ));

        gpuBlock *blk = reinterpret_cast<gpuBlock *>(malloc(sizeof(gpuBlock)));
        assert(blk!=NULL);
        blk->next = ctx->block;
        ctx->block = blk;

        CUDA_CHECK(cudaMalloc,(&blk->gpuBlk, tileSize));
        CUDA_CHECK(cudaMalloc,(reinterpret_cast<void **>(&blk->gpuResults), OutSize));
        CUDA_CHECK(cudaHostAlloc,(reinterpret_cast<void **>(&blk->cpuResults), OutSize, 0));
        CUDA_CHECK(cudaEventCreate,( &blk->event ));
        CUDA_CHECK(cudaStreamCreate,( &blk->stream ));
        }


    CUDAsetupPP();
    CUDAsetupPC();

    return ctx;
    }

extern "C"
void CUDA_finish(void *vctx) {
    CUDACTX ctx = reinterpret_cast<CUDACTX>(vctx);

    while(ctx->in != NULL) {
        gpuInput *in = ctx->in;
        ctx->in = in->next;
        cudaFree(in->in);
        cudaEventDestroy( in->event );
        free(in);
        }

    while(ctx->block != NULL) {
        gpuBlock *blk = ctx->block;
        ctx->block = blk->next;
        cudaFree(blk->gpuBlk);
        cudaFree(reinterpret_cast<void *>(blk->gpuResults));
        cudaFreeHost(reinterpret_cast<void *>(blk->cpuResults));
        cudaEventDestroy( blk->event );
        cudaStreamDestroy( blk->stream );
        free(blk);
        }

    free(ctx);
    }
#endif

extern "C"
void *CUDA_initialize(int iCore) {
    int nDevices;

    CUDA_CHECK(cudaGetDeviceCount,(&nDevices))
    if (nDevices == 0 ) return NULL;
    CUDA_CHECK(cudaSetDevice,(iCore % nDevices));
    CUDACTX ctx = reinterpret_cast<CUDACTX>(malloc(sizeof(struct cuda_ctx)));
    assert(ctx!=NULL);
    ctx->wqSize = 0;
    ctx->wqCuda = ctx->wqFree = NULL;

    CUDA_CHECK(cudaGetDeviceProperties,(&ctx->prop,iCore % nDevices));



#if 0
    for(i=0; i<iWorkQueueSize; ++i) {

        gpuInput *in = reinterpret_cast<gpuInput *>(malloc(sizeof(gpuInput)));
        assert(in!=NULL);
        in->next = ctx->in;
        ctx->in = in;
        CUDA_CHECK(cudaMalloc,(&in->in, ParticlesSize));
        CUDA_CHECK(cudaHostAlloc,(reinterpret_cast<void **>(&in->cpuIn), ParticlesSize, 0));
        CUDA_CHECK(cudaEventCreate,( &in->event ));

        gpuBlock *blk = reinterpret_cast<gpuBlock *>(malloc(sizeof(gpuBlock)));
        assert(blk!=NULL);
        blk->next = ctx->block;
        ctx->block = blk;

        CUDA_CHECK(cudaMalloc,(&blk->gpuBlk, tileSize));
        CUDA_CHECK(cudaMalloc,(reinterpret_cast<void **>(&blk->gpuResults), OutSize));
        CUDA_CHECK(cudaHostAlloc,(reinterpret_cast<void **>(&blk->cpuResults), OutSize, 0));
        CUDA_CHECK(cudaEventCreate,( &blk->event ));
        CUDA_CHECK(cudaStreamCreate,( &blk->stream ));
        }
#endif
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
	    work->doneFcn = doneWork;
	    work->next = cuda->wqCuda;
	    cuda->wqCuda = work;
	    }
    else {
        assert(0);
        }
    return 1;
    }
extern "C"
void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(vcuda);
    CUDAwqNode *work;
    while(cudaSize--) {
        work = reinterpret_cast<CUDAwqNode *>(malloc(sizeof(CUDAwqNode)));
        work->pHostBuf = CUDA_malloc(inCudaBufSize);
        work->pCudaBuf = CUDA_gpu_malloc(inCudaBufSize);
        CUDA_CHECK(cudaEventCreate,( &work->event ));
        CUDA_CHECK(cudaStreamCreate,( &work->stream ));
        work->ctx = NULL;
        work->checkFcn = NULL;
        work->doneFcn = NULL;
        work->next = cuda->wqFree;
        cuda->wqFree = work;
//        ++mdl->wqCudaMaxSize;
        }
    }
