#ifndef CUDAUTIL_H
#define CUDAUTIL_H

#ifdef __cplusplus
extern "C" {
#endif
#ifdef USE_CUDA
    void *CUDA_malloc(size_t nBytes);
    void CUDA_free(void *data);
    void *CUDA_gpu_malloc(size_t nBytes);
    void CUDA_gpu_free(void *data);
    void *CUDA_initialize(int iCore);
    void CUDA_finish(void *vctx);
    void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize);

    int CUDA_queue(void *vcuda, void *ctx,
	int (*initWork)(void *ctx,void *work),
	int (*checkWork)(void *ctx,void *work),
	int (*doneWork)(void *ctx));
    int CUDA_flushDone(void *vcuda);

#else
#include "simd.h"
#define CUDA_malloc SIMD_malloc
#define CUDA_free SIMD_free
#endif
#ifdef __cplusplus
    }
#endif

#ifdef __CUDACC__

#define CUDA_NUM_BANKS 16
#define CUDA_LOG_NUM_BANKS 4

#include "ilp.h"
#include "ilc.h"
#include "basetype.h"

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

typedef struct cuda_wq_node {
    /* We can put this on different types of queues */
    struct cuda_wq_node *next;
    void *ctx;
    int (*checkFcn)(void *,void *);
    int (*doneFcn)(void *);
    void *pHostBuf;
    void *pCudaBuf;
    cudaEvent_t event;     // Results have been copied back
    cudaStream_t stream;   // execution stream
    } CUDAwqNode;

typedef struct cuda_ctx {
    int iWorkQueueSize;
    gpuInput   *in;
    gpuBlock   *block;
    struct cudaDeviceProp prop;
    int wqSize;
    CUDAwqNode *wqCuda;
    CUDAwqNode *wqFree;
    } *CUDACTX;

#endif




#endif
