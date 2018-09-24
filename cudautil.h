/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CUDAUTIL_H
#define CUDAUTIL_H
#include "opa_queue.h"
#include "basetype.h"

//#define CUDA_STREAMS 16

#ifdef __cplusplus
extern "C" {
#endif
#ifdef USE_CUDA
    void *CUDA_initialize(int nCores, int iCore, OPA_Queue_info_t *queueWORK, OPA_Queue_info_t *queueREGISTER);
    void CUDA_finish(void *vctx);
    void CUDA_SetQueueSize(void *vcuda,int cudaSize, int inCudaBufSize, int outCudaBufSiz);

    int CUDA_queue(void *vcuda, void *ctx,
	int (*initWork)(void *ctx,void *work),
	int (*checkWork)(void *ctx,void *work));
    int CUDA_flushDone(void *vcuda);
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

#define CUDA_DEVICE __device__
inline __device__ bool testz(bool p) { return !p; }
inline __device__ float maskz_mov(bool p,float a) { return p ? a : 0.0f; }

template <typename T,unsigned int blockSize>
inline __device__ T warpReduce(/*volatile T * data, int tid,*/ T t) {
#if CUDART_VERSION >= 9000
    if (blockSize >= 32) t += __shfl_xor_sync(0xffffffff,t,16);
    if (blockSize >= 16) t += __shfl_xor_sync(0xffffffff,t,8);
    if (blockSize >= 8)  t += __shfl_xor_sync(0xffffffff,t,4);
    if (blockSize >= 4)  t += __shfl_xor_sync(0xffffffff,t,2);
    if (blockSize >= 2)  t += __shfl_xor_sync(0xffffffff,t,1);
#else
    if (blockSize >= 32) t += __shfl_xor(t,16);
    if (blockSize >= 16) t += __shfl_xor(t,8);
    if (blockSize >= 8)  t += __shfl_xor(t,4);
    if (blockSize >= 4)  t += __shfl_xor(t,2);
    if (blockSize >= 2)  t += __shfl_xor(t,1);
#endif
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

#ifdef __cplusplus
class WorkNode {
    OPA_Queue_element_hdr_t hdr;
    OPA_Queue_info_t *pwqDone; // Where to return this
public:
    void queueTo(OPA_Queue_info_t *q,        // Send to this queue for processing
    	         OPA_Queue_info_t *qBack=0); // Return to this queue when complete
    void sendBack(); // Return this work node back to the requested queue
    };

class CudaWorkNode : public WorkNode {
    int inBufferSize, outBufferSize;
    void *pHostBufToGPU, *pHostBufFromGPU;
    void *pCudaBufIn, *pCudaBufOut;
    cudaStream_t stream;     // execution stream
    dim3 dimBlock, dimGrid;

public:
    explicit CudaWorkNode(int inSize,int outSize=0);
    virtual ~CudaWorkNode();
    virtual int init(void *context) = 0; // Start the work (mdl thread)
    virtual int done(void *context) = 0; // finish the work (mdl thread)
    virtual int save(void *context) = 0; // save the results (caller)
    void create_buffers();
    void destroy_buffers();
    };
#endif

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
