#ifndef CUDAUTIL_H
#define CUDAUTIL_H

#ifdef __cplusplus
extern "C" {
#endif
#ifdef USE_CUDA
    void *CUDA_malloc(size_t nBytes);
    void CUDA_free(void *data);
    void *CUDA_initialize(int iProc,int iWorkQueueSize, size_t tileSize, size_t ParticlesSize, size_t OutSize);
    void CUDA_finish(void *vctx);
#else
#define CUDA_malloc(n) malloc(n)
#define CUDA_free(p) free(p)
#endif
#ifdef __cplusplus
    }
#endif

#ifdef __CUDACC__

#include "pkd.h"

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line);

#define CUDA_CHECK(f,a) {cudaError_t rc = (f)a; if (rc!=cudaSuccess) CUDA_Abort(rc,#f,__FILE__,__LINE__);}

typedef struct gpu_input {
    struct gpu_input *next;
    PINFOIN *cpuIn;
    PINFOIN *in;       // Input particles
    cudaEvent_t event; // Input particles have been copied
    int nRefs;
    } gpuInput;

typedef struct gpu_block {
    struct gpu_block *next;
    union {
	ILP_BLK *gpuBlk;       // Input block
	ILC_BLK *gpuBlkILC;    // Input block
	};
    PINFOOUT *gpuResults;  // Results in GPU memory
    PINFOOUT *cpuResults;  // Results in CPU memory
    cudaEvent_t event;     // Results have been copied back
    cudaStream_t stream;   // execution stream
    } gpuBlock;

typedef struct cuda_ctx {
    int iWorkQueueSize;
    gpuInput   *in;
    gpuBlock   *block;
    } *CUDACTX;

#endif




#endif
