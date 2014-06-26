#ifndef CUDAUTIL_H
#define CUDAUTIL_H
#include "basetype.h"
#include "ilp.h"

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
    void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize, int outCudaBufSiz);

    int CUDA_queue(void *vcuda, void *ctx,
	int (*initWork)(void *ctx,void *work),
	int (*checkWork)(void *ctx,void *work),
	int (*doneWork)(void *ctx));
    int CUDA_flushDone(void *vcuda);
    int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile);
    int CUDA_queuePC(void *cudaCtx,workParticle *wp, ILCTILE tile);
    void CUDA_sendWork(void *cudaCtx);

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
#define cfOffset(i) ((i) + ((i)>>CUDA_LOG_NUM_BANKS))


template <typename T,unsigned int blockSize>
inline __device__ T warpReduce(volatile T * data, int tid) {
    T t = data[tid];
    if (tid<blockSize/2) {
	if (blockSize >= 64) data[tid] = t = t + data[tid + 32]; 
	if (blockSize >= 32) data[tid] = t = t + data[tid + 16]; 
	if (blockSize >= 16) data[tid] = t = t + data[tid +  8]; 
	if (blockSize >= 8)  data[tid] = t = t + data[tid +  4]; 
	if (blockSize >= 4)  data[tid] = t = t + data[tid +  2]; 
	if (blockSize >= 2)  data[tid] = t = t + data[tid +  1]; 
	}
    return t;
    }

template <typename T,unsigned int blockSize>
inline __device__ T warpReduce(volatile T * data, int tid, T t) {
    data[tid] = t;
    if (tid<blockSize/2) {
	if (blockSize >= 64) data[tid] = t = t + data[tid + 32]; 
	if (blockSize >= 32) data[tid] = t = t + data[tid + 16]; 
	if (blockSize >= 16) data[tid] = t = t + data[tid +  8]; 
	if (blockSize >= 8)  data[tid] = t = t + data[tid +  4]; 
	if (blockSize >= 4)  data[tid] = t = t + data[tid +  2]; 
	if (blockSize >= 2)  data[tid] = t = t + data[tid +  1]; 
	}
    return t;
    }

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStore(volatile T * data, int tid,T *result) {
    T t = warpReduce<T,blockSize>(data,tid);
    if (tid==0) *result = t;
    }

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStore(volatile T * data, int tid,T t,T *result) {
    t = warpReduce<T,blockSize>(data,tid,t);
    if (tid==0) *result = t;
    }

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line);

#define CUDA_CHECK(f,a) {cudaError_t rc = (f)a; if (rc!=cudaSuccess) CUDA_Abort(rc,#f,__FILE__,__LINE__);}
#define CUDA_PP_MAX_BUFFERED 128

typedef struct cuda_wq_node {
    /* We can put this on different types of queues */
    struct cuda_wq_node *next;
    void *ctx;
    int (*checkFcn)(void *,void *);
    void *pHostBuf;
    void *pCudaBufIn;
    void *pCudaBufOut;
    cudaEvent_t event;       // Results have been copied back
    cudaStream_t stream;     // execution stream
    workParticle *ppWP[CUDA_PP_MAX_BUFFERED];
    int ppNI[CUDA_PP_MAX_BUFFERED];
    int ppSizeIn; // Number of bytes consumed in the buffer
    int ppSizeOut; // Number of bytes consumed in the buffer
    int ppnBlocks;
    int ppnBuffered;
    } CUDAwqNode;

typedef struct cuda_ctx {
    struct cudaDeviceProp prop;
    CUDAwqNode *wqCuda;
    CUDAwqNode *wqFree;
    CUDAwqNode *nodePP; // We are building a PP request
    CUDAwqNode *nodePC; // We are building a PC request
    int iWorkQueueSize;
    int inCudaBufSize, outCudaBufSize;
    } *CUDACTX;



#endif




#endif
