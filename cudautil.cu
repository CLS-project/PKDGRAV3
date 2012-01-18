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

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line) {
    fprintf(stderr,"%s error %d in %s(%d)\n%s\n", fname, rc, file, line, cudaGetErrorString(rc));
    exit(1);
    }

extern "C"
void *CUDA_initialize(int iWorkQueueSize, size_t tileSize, size_t ParticlesSize, size_t OutSize) {
    int nDevices, i;

    CUDA_CHECK(cudaGetDeviceCount,(&nDevices))
    if (nDevices == 0 ) return NULL;
    CUDA_CHECK(cudaSetDevice,(nDevices-1));

    CUDACTX ctx = reinterpret_cast<CUDACTX>(malloc(sizeof(struct cuda_ctx)));
    assert(ctx!=NULL);
    ctx->iWorkQueueSize = iWorkQueueSize;
    ctx->in = NULL;
    ctx->block = NULL;

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
