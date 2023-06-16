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

/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#include "pkd_config.h"
#include <stdio.h>
#include "cooperative_groups/memcpy_async.h"

//****************************************************************************************************
// The strategy is as follows:
// - We want each thread to handle one source interaction (but multiple sink particles)
// - We copy the interaction into shared memory for speed and to reduce register pressure
//   - LIMIT: P-C 64k shared / (29 x 4 bytes) = 512 threads = 59'392 bytes of shared memory
// - We want 2048 threads resident, so we need to spread particles amongst 4 threads
//   - Thread block: (32,4,1) = 128 threads, the "work unit"
//   - Resident blocks: 2048 / 128 = 16 (but subject to potential register limits)
//****************************************************************************************************

inline __device__ bool testz(bool p) { return !p; }
inline __device__ float maskz_mov(bool p,float a) { return p ? a : 0.0f; }

#include "basetype.h"
#include "cudapppc.h"
#include "gravity/pp.h"
#include "gravity/pc.h"
#include "gravity/ilp.h"
#include "gravity/ilc.h"
#include "cuda/reduce.h"

#define WIDTH 32

/// Evaluate a single P-P interaction
template<bool bGravStep>
__device__ __forceinline__ auto evalInteraction(const gpu::ppInput pp,const gpu::Blk<WIDTH,ilpTile> *__restrict__ blk,int i) {
    return EvalPP<float,bool>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
                              blk->dx[i], blk->dy[i], blk->dz[i], blk->fourh2[i], blk->m[i],
                              pp.ax, pp.ay, pp.az, pp.dImaga);
}

/// Evaluate a single P-C interaction
template<bool bGravStep>
__device__ __forceinline__ auto evalInteraction(const gpu::ppInput pp,const gpu::Blk<WIDTH,ilcTile> *__restrict__ blk,int i) {
    return EvalPC<float,bool,bGravStep>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
                                        blk->dx[i],blk->dy[i],blk->dz[i],blk->m[i],blk->u[i],
                                        blk->xxxx[i],blk->xxxy[i],blk->xxxz[i],blk->xxyz[i],blk->xxyy[i],
                                        blk->yyyz[i],blk->xyyz[i],blk->xyyy[i],blk->yyyy[i],
                                        blk->xxx[i],blk->xyy[i],blk->xxy[i],blk->yyy[i],blk->xxz[i],blk->yyz[i],blk->xyz[i],
                                        blk->xx[i],blk->xy[i],blk->xz[i],blk->yy[i],blk->yz[i],
#ifdef USE_DIAPOLE
                                        blk->x[i],blk->y[i],blk->z[i],
#endif
                                        pp.ax, pp.ay, pp.az, pp.dImaga);
}

// Reduction of a P-P or P-C interaction
template <bool bGravStep,class RESULT>
__device__ __forceinline__ void reduceInteraction(gpu::ppResult &out, RESULT result) {
    warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.ax,&out.ax);
    warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.ay,&out.ay);
    warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.az,&out.az);
    warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.pot,&out.fPot);
    if (bGravStep) {
        warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.ir,&out.dirsum);
        warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.norm,&out.normsum);
    }
}

/// CUDA Kernel to evaluate a block of P-P or P-C interactions
template <bool bGravStep,class INPUT,class BLK,class OUTPUT>
__global__ void cudaInteract(
    const gpu::ppWorkUnit *__restrict__ work,
    const INPUT *__restrict__ pPart, // e.g., gpu::ppInput
    const BLK *__restrict__ gblk,
    OUTPUT *out) { // e.g., gpu::ppResult
    extern __shared__ char cache[];
    auto sblk = reinterpret_cast<BLK *>(cache);
    auto block = cooperative_groups::this_thread_block();

    work += blockIdx.x;
    gblk += blockIdx.x;
    pPart += work->iP;
    out += work->iP;

    // Load the interactions into shared memory
    cooperative_groups::memcpy_async(block, sblk, gblk, sizeof(*gblk));
    cooperative_groups::wait(block); // Joins all threads, waits for all copies to complete

    decltype(evalInteraction<bGravStep>(pPart[0],sblk,0)) result {0,0,0,0,0,0};
    for (auto iP=threadIdx.y; iP<work->nP; iP += blockDim.y) {
        if (threadIdx.x < work->nI) {
            result = evalInteraction<bGravStep>(pPart[iP],sblk,threadIdx.x);
        }
        reduceInteraction<bGravStep>(out[iP],result);
    }
}

void pkdParticleWorkDone(workParticle *wp);

template<class TILE,int N>
void MessagePPPC<TILE,N>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) {
    typedef gpu::Blk<N,TILE> BLK;
    auto *pCudaOutput = reinterpret_cast<gpu::ppResult *>(pCudaBufOut);

    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, this->pHostBufIn, this->requestBufferCount, cudaMemcpyHostToDevice, stream));

    // The interation blocks
    auto *__restrict__ blkCuda = reinterpret_cast<BLK *>(pCudaBufIn);
    // The particle information
    auto *__restrict__ partCuda = reinterpret_cast<gpu::ppInput *>(blkCuda + this->nTotalInteractionBlocks);
    // The interaction block descriptors
    auto *__restrict__ wuCuda = reinterpret_cast<gpu::ppWorkUnit *>(partCuda + this->nTotalParticles);

    cudaMemsetAsync(pCudaBufOut,0,this->resultsBufferCount,stream);

    dim3 dimBlock( N, 8, 1 );
    dim3 dimGrid( this->nGrid, 1,1);
    if (bGravStep) {
        cudaInteract<true>
        <<<dimGrid, dimBlock, sizeof(BLK), stream>>>
        (wuCuda,partCuda,blkCuda,pCudaOutput );
    }
    else {
        cudaInteract<false>
        <<<dimGrid, dimBlock, sizeof(BLK), stream>>>
        (wuCuda,partCuda,blkCuda,pCudaOutput );
    }

    CUDA_CHECK(cudaMemcpyAsync,(this->pHostBufOut, pCudaBufOut, this->resultsBufferCount, cudaMemcpyDeviceToHost, stream) );
}

template<class TILE,int N>
void MessagePPPC<TILE,N>::finish() {
    auto *pR = reinterpret_cast<gpu::ppResult *>(this->pHostBufOut);

    for ( auto &w : this->work ) {
        auto nP = w.wp->nP;
        auto *pInfoOut = w.wp->pInfoOut;
        for (auto ip=0; ip<nP; ++ip) {
            pInfoOut[ip].a[0]    += pR[ip].ax;
            pInfoOut[ip].a[1]    += pR[ip].ay;
            pInfoOut[ip].a[2]    += pR[ip].az;
            pInfoOut[ip].fPot    += pR[ip].fPot;
            pInfoOut[ip].dirsum  += pR[ip].dirsum;
            pInfoOut[ip].normsum += pR[ip].normsum;
        }
        pR += this->align_nP(nP);
        pkdParticleWorkDone(w.wp);
    }
    this->clear();
    freeQueue.enqueue(*this);
}

/*****************************************************************************\
*   Explicit instantiation for methods we need elsewhere
\*****************************************************************************/

template void MessagePP::finish();
template void MessagePC::finish();
template void MessagePP::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessagePC::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);
