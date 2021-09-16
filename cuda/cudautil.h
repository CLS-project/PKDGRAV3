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
#ifdef USE_CUDA
#include "opa_queue.h"
#include "basetype.h"

#define CUDA_STREAMS 8
#define CUDA_WP_MAX_BUFFERED 128

#ifdef __cplusplus
#include "mdl.h"
#include "mdlcuda.h"
#include <vector>

class cudaDataMessage : public mdl::cudaMessage {
protected:
    // This should be handled in a better way, but for now let the cheese happen.
    const size_t requestBufferSize = 2*1024*1024;
    const size_t resultsBufferSize = 2*1024*1024;
    void *pHostBufIn, *pHostBufOut;
protected:
    virtual void launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) override = 0;
    virtual void finish() = 0;
public:
    explicit cudaDataMessage();
    virtual ~cudaDataMessage();
    };

class MessageEwald : public cudaDataMessage {
protected:
    class CudaClient &cuda;
    std::vector<workParticle *> ppWP; // [CUDA_WP_MAX_BUFFERED]
    virtual void launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
    int nParticles, nMaxParticles;
public:
    explicit MessageEwald(class CudaClient &cuda);
    bool queue(workParticle *work); // Add work to this message; return false if this is not possible
    };

class MessageEwaldSetup : public mdl::cudaMessage {
protected:
    virtual void launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) override;
public:
    explicit MessageEwaldSetup(struct EwaldVariables * const ew, EwaldTable * const ewt);
protected:
    struct EwaldVariables * const ewIn;
    EwaldTable * const ewt;
    std::vector<momFloat> dLx, dLy, dLz;
    std::vector<int> ibHole;
    };

template<class TILE,int nIntPerTB, int nIntPerWU>
class MessagePPPC : public cudaDataMessage {
protected:
    mdl::messageQueue<MessagePPPC> &freeQueue;
    bool bGravStep;
    size_t requestBufferCount, resultsBufferCount;
    int nInteractionBlocks, nGrid;
    struct workInformation {
    	workParticle *wp;
	size_t nInteractions;
	workInformation(workParticle *wp, size_t nInteractions) : wp(wp), nInteractions(nInteractions) {}
	};
    std::vector<workInformation> work; // [CUDA_WP_MAX_BUFFERED]
    virtual void launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
public:
    explicit MessagePPPC(mdl::messageQueue<MessagePPPC> &freeQueue);
	// : cuda(cuda), requestBufferCount(0), resultsBufferCount(0), nInteractionBlocks(0) {
	// work.reserve(CUDA_WP_MAX_BUFFERED);
	// }

    void clear();
    bool queue(workParticle *wp, TILE *tile, bool bGravStep);
    MessagePPPC<TILE,nIntPerTB,nIntPerWU> &prepare();
    };
typedef MessagePPPC<ilpTile,128,128> MessagePP;
typedef MessagePPPC<ilcTile,128,32>  MessagePC;

class CudaClient {
    friend class MessageEwald;
protected:
    // Message on this queue need to have finish() called on them.
    // The finish() routine will add them back to the correct "free" list.
    mdl::messageQueue<cudaDataMessage> doneCuda;
    mdl::messageQueue<MessageEwald> freeEwald;
    MessageEwald *ewald;
    mdl::messageQueue<MessagePP> freePP;
    MessagePP *pp;
    mdl::messageQueue<MessagePC> freePC;
    MessagePC *pc;

    int nEwhLoop;
    std::list<MessageEwald> free_Ewald, busy_Ewald;
    mdl::mdlClass &mdl;
protected:
    template<class MESSAGE> void flush(MESSAGE * &M);
    template<class MESSAGE,class QUEUE,class TILE>
    int queue(MESSAGE * &m, QUEUE &Q, workParticle *wp, TILE *tile, bool bGravStep);
public:
    explicit CudaClient(mdl::mdlClass &mdl);
    void flushCUDA();
    int  queuePP(workParticle *wp, ilpTile *tile, bool bGravStep);
    int  queuePC(workParticle *wp, ilcTile *tile, bool bGravStep);
    int  queueEwald(workParticle *wp);
    void setupEwald(struct EwaldVariables * const ew, EwaldTable * const ewt);
    };
#endif

#ifdef __cplusplus
extern "C" {
#endif
#ifdef USE_CUDA
    void *CudaClientInitialize(void *vmdl);
    void CudaClientFlush(void *vcudaClient);
    int CudaClientQueueEwald(void *vcudaClient, workParticle *work);
    int CudaClientQueuePP(void *vcudaClient, workParticle *work, struct ilpTile *tile, int bGravStep);
    int CudaClientQueuePC(void *vcudaClient, workParticle *work, struct ilcTile *tile, int bGravStep);
#else
#if !defined(__CUDACC__)
#include "core/simd.h"
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

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line);

#define CUDA_CHECK(f,a) {cudaError_t rc = (f)a; if (rc!=cudaSuccess) CUDA_Abort(rc,#f,__FILE__,__LINE__);}
#define CUDA_RETURN(f,a) {cudaError_t rc = (f)a; if (rc!=cudaSuccess) return rc;}

double CUDA_getTime();

#endif
#endif
#endif
