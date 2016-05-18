#ifndef CUDAUTIL_H
#define CUDAUTIL_H
#include "basetype.h"
#include "ilp.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifdef USE_CUDA
    void CUDA_nvtxRangePush(char *name);
    void CUDA_nvtxRangePop();
    void *CUDA_initialize(int nCores, int iCore);
    void CUDA_finish(void *vctx);
    void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize, int outCudaBufSiz);

    int CUDA_queue(void *vcuda, void *ctx,
	int (*initWork)(void *ctx,void *work),
	int (*checkWork)(void *ctx,void *work));
    int CUDA_flushDone(void *vcuda);
    int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile, int bGravStep);
    int CUDA_queuePC(void *cudaCtx,workParticle *wp, ILCTILE tile, int bGravStep);
    void CUDA_sendWork(void *cudaCtx);
    void CUDA_checkForRecovery(void *vcuda);
    void pkdAccumulateCUDA(void *vpkd,workEwald *we,gpuEwaldOutput *fromGPU);

#else
#include "simd.h"
#define CUDA_malloc SIMD_malloc
#define CUDA_free SIMD_free
#endif
#ifdef __cplusplus
    }
#endif

#ifdef __CUDACC__
#include "sm_30_intrinsics.h"

#define CUDA_NUM_BANKS 16
#define CUDA_LOG_NUM_BANKS 4
#define cfOffset(i) ((i) + ((i)>>CUDA_LOG_NUM_BANKS))

template <typename T,unsigned int blockSize>
inline __device__ T warpReduce(/*volatile T * data, int tid,*/ T t) {
    if (blockSize >= 32) t += __shfl_xor(t,16);
    if (blockSize >= 16) t += __shfl_xor(t,8);
    if (blockSize >= 8)  t += __shfl_xor(t,4);
    if (blockSize >= 4)  t += __shfl_xor(t,2);
    if (blockSize >= 2)  t += __shfl_xor(t,1);
    return t;
    }

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStore(volatile T * data, int tid,T *result) {
    T t = warpReduce<T,blockSize>(data[tid]);
    if (tid==0) *result = t;
    }

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStore(int tid,T t,T *result) {
    t = warpReduce<T,blockSize>(t);
    if (tid==0) *result = t;
    }

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line);

#define CUDA_CHECK(f,a) {cudaError_t rc = (f)a; if (rc!=cudaSuccess) CUDA_Abort(rc,#f,__FILE__,__LINE__);}
#define CUDA_RETURN(f,a) {cudaError_t rc = (f)a; if (rc!=cudaSuccess) return rc;}

#define CUDA_PP_MAX_BUFFERED 128

typedef struct cuda_wq_node {
    /* We can put this on different types of queues */
    struct cuda_wq_node *next;
    void *ctx;
    int (*initFcn)(void *,void *);
    int (*checkFcn)(void *,void *);
    void *pHostBufToGPU, *pHostBufFromGPU;
    void *pCudaBufIn;
    void *pCudaBufOut;
    double startTime;
#ifdef USE_CUDA_EVENTS
    cudaEvent_t event;       // Results have been copied back
#endif
    cudaStream_t stream;     // execution stream
    workParticle *ppWP[CUDA_PP_MAX_BUFFERED];
    int ppNI[CUDA_PP_MAX_BUFFERED];
    int ppSizeIn; // Number of bytes consumed in the buffer
    int ppSizeOut; // Number of bytes consumed in the buffer
    int ppnBlocks;
    int ppnBuffered;
    int bGravStep;
    union {
	struct {
	    size_t nBufferIn;
	    size_t nBufferOut;
	    int nGrid;
	    } pppc;
	};
    } CUDAwqNode;

typedef struct cuda_ctx {
    struct cudaDeviceProp prop;
    CUDAwqNode *wqCuda;
    CUDAwqNode *wqFree;
    CUDAwqNode *nodePP; // We are building a PP request
    CUDAwqNode *nodePC; // We are building a PC request
    int nWorkQueueSize, nWorkQueueBusy;
    int inCudaBufSize, outCudaBufSize;
    int epoch;
    int nCores, iCore;

    struct EwaldVariables *ewIn;
    EwaldTable *ewt;
    cudaEvent_t eventEwald;       // Results have been copied back
    cudaStream_t streamEwald;     // execution stream

    char hostname[256];
    } *CUDACTX;

void CUDA_attempt_recovery(CUDACTX cuda,cudaError_t errorCode);
double CUDA_getTime();
#endif

#endif
