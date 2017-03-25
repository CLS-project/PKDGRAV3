#ifndef CUDAUTIL_H
#define CUDAUTIL_H
#include "opa_queue.h"
#include "basetype.h"
#include "ilp.h"

//#define CUDA_STREAMS 16

#ifdef __cplusplus
extern "C" {
#endif
#ifdef USE_CUDA
    void CUDA_nvtxRangePush(char *name);
    void CUDA_nvtxRangePop();
    void *CUDA_initialize(int nCores, int iCore, OPA_Queue_info_t *queueWORK, OPA_Queue_info_t *queueREGISTER);
    void CUDA_finish(void *vctx);
    void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize, int outCudaBufSiz);

    int CUDA_queue(void *vcuda, void *ctx,
	int (*initWork)(void *ctx,void *work),
	int (*checkWork)(void *ctx,void *work));
    int CUDA_flushDone(void *vcuda);
    int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile, int bGravStep);
    int CUDA_queuePC(void *cudaCtx,workParticle *wp, ILCTILE tile, int bGravStep);
    int CUDA_queueEwald(void *cudaCtx,workParticle *work);
    void CUDA_sendWork(void *cudaCtx);
    void CUDA_flushEwald(void *cudaCtx);
    void CUDA_checkForRecovery(void *vcuda);
    void CUDA_startWork(void *vcuda,OPA_Queue_info_t *queueWORK);
    void CUDA_registerBuffers(void *vcuda, OPA_Queue_info_t *queueWORK);
#else
#if !defined(__CUDACC__)
#include "simd.h"
#define CUDA_malloc SIMD_malloc
#define CUDA_free SIMD_free
#endif
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

#define CUDA_WP_MAX_BUFFERED 128

typedef struct cuda_wq_node {
    union {
	OPA_Queue_element_hdr_t hdr;
	struct cuda_wq_node *next;
	} q;
    OPA_Queue_info_t *pwqDone; // Where to return this
    /* We can put this on different types of queues */
//    struct cuda_wq_node *next;
    void *ctx;
    int (*initFcn)(void *ctx,void *work);
    int (*doneFcn)(void *ctx,void *work);
    void (*dumpFcn)(struct cuda_wq_node *work);
    const char *kernelName;
    int inBufferSize, outBufferSize;
    void *pHostBufToGPU, *pHostBufFromGPU;
    void *pCudaBufIn, *pCudaBufOut;
    double startTime;
#ifdef USE_CUDA_EVENTS
    cudaEvent_t event;       // Results have been copied back
#endif
    cudaEvent_t eventCopyDone;
    cudaEvent_t eventKernelDone;
    cudaStream_t stream;     // execution stream
    dim3 dimBlock, dimGrid;
    workParticle *ppWP[CUDA_WP_MAX_BUFFERED];
    int ppNI[CUDA_WP_MAX_BUFFERED];
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
	struct {
	    int nParticles;
	    } ewald;
	struct {
	    struct EwaldVariables *ewIn;
	    EwaldTable *ewt;
	    } initEwald;
        };
    } CUDAwqNode;

#ifdef CUDA_STREAMS
typedef struct cuda_stream {
    union {
	struct cuda_stream *next;
	} q;
    cudaEvent_t eventCopyDone;
    cudaEvent_t eventKernelDone;
    cudaStream_t stream;     // execution stream
    CUDAwqNode *node;
    } cudaStream;
#endif

typedef struct cuda_ctx {
    struct cudaDeviceProp prop;
#ifdef CUDA_STREAMS
    cudaStream *wqCudaFree;  // Private to this CUDACTX
    cudaStream *wqCudaBusy;  // Private to this CUDACTX
#else
    CUDAwqNode *wqCudaFree;  // Private to this CUDACTX
    CUDAwqNode *wqCudaBusy;  // Private to this CUDACTX
#endif
    OPA_Queue_info_t wqFree; // We can receive from another thread
    OPA_Queue_info_t wqDone; // We can receive from another thread
    OPA_Queue_info_t *queueWORK;
    OPA_Queue_info_t *queueREGISTER;
    CUDAwqNode *nodePP; // We are building a PP request
    CUDAwqNode *nodePC; // We are building a PC request
    CUDAwqNode *nodeEwald; // We are building an Ewald request
    int nwqCudaBusy;
    int nWorkQueueSize, nWorkQueueBusy;
    int inCudaBufSize, outCudaBufSize;
    int epoch;
    int nCores, iCore;
    uint64_t nKernelLaunches;

    struct EwaldVariables *ewIn;
    EwaldTable *ewt;
#ifndef CUDA_STREAMS
    cudaEvent_t eventEwald;       // Results have been copied back
    cudaStream_t streamEwald;     // execution stream
#endif
    char hostname[256];
    } *CUDACTX;

void CUDA_attempt_recovery(CUDACTX cuda,cudaError_t errorCode);
double CUDA_getTime();
CUDAwqNode *getNode(CUDACTX cuda);
#endif

#endif
