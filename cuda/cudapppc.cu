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
#include "cudapppc.h"
#include "gravity/pp.h"
#include "gravity/pc.h"
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

#define SYNC_RATE 16  // Must be: 1, 2, 4, 8, 16
#define WIDTH 32

#include "gravity/ilp.h"
#include "gravity/ilc.h"
#include "basetype.h"

/*
** The following are the basically the same as ILP_BLK and ILC_BLK,
** but we need to be able to alter their sizes.
*/

template <int n,class TILE> struct Blk;
template <int n> struct Blk<n,ilpTile> : BlockPP<n> {
    // float dx[n], dy[n], dz[n];    /* Offset from ilp->cx, cy, cz */
    // float m[n];             /* Mass */
    // float fourh2[n];        /* Softening: calculated */
};
template <int n> struct Blk<n,ilcTile> : BlockPC<n> {
    // float dx[n],dy[n],dz[n];
    // float xxxx[n],xxxy[n],xxxz[n],xxyz[n],xxyy[n],yyyz[n],xyyz[n],xyyy[n],yyyy[n];
    // float xxx[n],xyy[n],xxy[n],yyy[n],xxz[n],yyz[n],xyz[n];
    // float xx[n],xy[n],xz[n],yy[n],yz[n];
    // float x[n],y[n],z[n];
    // float m[n],u[n];
};

__device__ __forceinline__ auto evalInteraction(const ppInput pp,const Blk<WIDTH,ilpTile> *__restrict__ blk,int i) {
    return EvalPP<float,bool>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
                              blk->dx[i], blk->dy[i], blk->dz[i], blk->fourh2[i], blk->m[i],
                              pp.ax, pp.ay, pp.az, pp.dImaga);
}

__device__ __forceinline__ auto evalInteraction(const ppInput pp,const Blk<WIDTH,ilcTile> *__restrict__ blk,int i) {
    return EvalPC<float,bool,true>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
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

template <bool bGravStep,class BLK>
__global__ void cudaInteract(
    const ppWorkUnit *__restrict__ work,
    const ppInput *__restrict__ pPart,
    const BLK *__restrict__ gblk,
    ppResult *out) {
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

    decltype(evalInteraction(pPart[0],sblk,0)) result {0,0,0,0,0,0};
    for (auto iP=threadIdx.y; iP<work->nP; iP += blockDim.y) {
        if (threadIdx.x < work->nI) {
            result = evalInteraction(pPart[iP],sblk,threadIdx.x);
        }
        warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.ax,&out[iP].ax);
        warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.ay,&out[iP].ay);
        warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.az,&out[iP].az);
        warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.pot,&out[iP].fPot);
        if (bGravStep) {
            warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.ir,&out[iP].dirsum);
            warpReduceAndStoreAtomic<float,32>(threadIdx.x,result.norm,&out[iP].normsum);
        }
    }
}

void pkdParticleWorkDone(workParticle *wp);

/*****************************************************************************\
*   CudaClient interface (new!)
\*****************************************************************************/

int CudaClient::queuePP(workParticle *work, ilpTile &tile, bool bGravStep) {
    return queue(pp,freePP,work,tile,bGravStep);
}

int CudaClient::queuePC(workParticle *work, ilcTile &tile, bool bGravStep) {
    return queue(pc,freePC,work,tile,bGravStep);
}

template<class MESSAGE,class QUEUE,class TILE>
int CudaClient::queue(MESSAGE *&M,QUEUE &Q, workParticle *work, TILE &tile, bool bGravStep) {
    if (M) { // If we are in the middle of building data for a kernel
        if (M->queue(work,tile,bGravStep)) return work->nP; // Successfully queued
        flush(M); // The buffer is full, so send it
    }
    mdl.flushCompletedCUDA();
    if (Q.empty()) return 0; // No buffers so the CPU has to do this part
    M = & Q.dequeue();
    if (M->queue(work,tile,bGravStep)) return work->nP; // Successfully queued
    return 0; // Not sure how this would happen, but okay.
}

template<class MESSAGE>
void CudaClient::flush(MESSAGE *&M) {
    if (M) {
        mdl.enqueue(M->prepare());
        M = nullptr;
    }
}

/*****************************************************************************\
*   DATA LAYOUT
*
*   The memory block sent to the GPU has three distinct sections:
*   1. An array of interaction blocks
*   2. An array of particles
*   3. An array of interaction block descriptors
*
*   When an interaction list is queued, we:
*   1. Make sure that the additional interaction blocks and the pending
*      descriptors and particles will fit, returning false if they won't.
*   2. The individual blocks are copied into the output buffer.
*   3. The workParticle structure is saved for later use
\*****************************************************************************/

template<class TILE>
MessagePPPC<TILE>::MessagePPPC(mdl::messageQueue<MessagePPPC> &freeQueue)
    : freeQueue(freeQueue), requestBufferCount(0), resultsBufferCount(0), nTotalInteractionBlocks(0), nTotalParticles(0) {
    work.reserve(CUDA_WP_MAX_BUFFERED);
}

// This function empties the message for subsequent queue requests
template<class TILE>
void MessagePPPC<TILE>::clear() {
    work.clear();
    requestBufferCount = resultsBufferCount = 0;
    nTotalInteractionBlocks = nTotalParticles = 0;
}

template<int n,class TILE>
int copyBLKs2(Blk<n,TILE> *out, TILE &in) {
    assert(n==ILP_PART_PER_BLK);
    auto nIlp = in.count();
    int i, nBlk = (nIlp+n-1) / n;
    for (i=0; i<nBlk; ++i) memcpy(&out[i],&in[i],sizeof(out[i]));
    return nBlk;
}

static double getFlops(workParticle *wp,ilpTile &tile) {
    return COST_FLOP_PP*wp->nP*tile.count();
}
static double getFlops(workParticle *wp,ilcTile &tile) {
    return COST_FLOP_PC*wp->nP*tile.count();
}

// Queue all of the work of a tile by adding it to our buffer
template<class TILE>
bool MessagePPPC<TILE>::queue(workParticle *wp, TILE &tile, bool bGravStep) {
    typedef Blk<WIDTH,TILE> BLK;
    if (work.size() == CUDA_WP_MAX_BUFFERED) return false;  // Too many work packages so send the work
    this->bGravStep = bGravStep;
    const auto nP = wp->nP;                                 // Queue this many particles
    const auto nI = tile.count();                           // ... operating on this many interactions

    if (inputSize<BLK>(nP,nI) > requestBufferSize - 1024) return false;     // Refuse if this tile won't fit in this buffer
    if (outputSize<BLK>(nP,nI) > requestBufferSize) return false;           // Refuse if this response won't fit
    auto blk = reinterpret_cast<Blk<WIDTH,TILE> *>(pHostBufIn);             // Copy the blocks to the input buffer
    nTotalInteractionBlocks += copyBLKs2(blk+nTotalInteractionBlocks,tile); // (the ILP tile can now be freed/reused)
    nTotalParticles += align_nP(nP);

    ++wp->nRefs;
    work.emplace_back(wp,tile.count());
    wp->dFlopSingleGPU += getFlops(wp,tile);
    return true;
}

// Final preparation before sending this message to the GPU thread.
// We need to setup the descriptors and add the particles.
template<class TILE>
MessagePPPC<TILE> &MessagePPPC<TILE>::prepare() {
    typedef Blk<WIDTH,TILE> BLK;
    // The interation blocks -- already copied to the host memory
    auto *__restrict__ blkHost = reinterpret_cast<BLK *>(pHostBufIn);
    // The particle information
    auto *__restrict__ partHost = reinterpret_cast<ppInput *>(blkHost + nTotalInteractionBlocks);
    // The interaction block descriptors
    auto *__restrict__ wuHost = reinterpret_cast<ppWorkUnit *>(partHost + nTotalParticles);

    uint32_t iP = 0;
    for ( auto &w : work ) {
        int nI = w.nInteractions;
        auto nP = w.wp->nP;
        auto *pInfoIn = w.wp->pInfoIn;
        auto nBlocks = (nI + BLK::width - 1) / BLK::width;

        // Generate a interaction block descriptor for each block
        for (auto j=0; j<nBlocks; ++j) {
            wuHost->nP = nP;
            wuHost->iP = iP;
            wuHost->nI = nI > BLK::width ? BLK::width : nI;
            nI -= wuHost->nI;
            ++wuHost;
        }

        // Copy in nP particles
        for (auto j=0; j<nP; ++j) {
            partHost[j].dx =  pInfoIn[j].r[0];
            partHost[j].dy =  pInfoIn[j].r[1];
            partHost[j].dz =  pInfoIn[j].r[2];
            partHost[j].ax =  pInfoIn[j].a[0];
            partHost[j].ay =  pInfoIn[j].a[1];
            partHost[j].az =  pInfoIn[j].a[2];
            partHost[j].fSoft2 = pInfoIn[j].fSmooth2;
            /*partHost[j].dImaga = 0;*/
        }
        nP = align_nP(nP);
        partHost += nP;
        iP += nP;
    }
    requestBufferCount = reinterpret_cast<char *>(wuHost) - reinterpret_cast<char *>(pHostBufIn);
    resultsBufferCount = iP * sizeof(ppResult);
    nGrid = nTotalInteractionBlocks;

    assert(requestBufferCount <= requestBufferSize);
    assert(resultsBufferCount <= resultsBufferSize);
    return *this;
}

template<class TILE>
void MessagePPPC<TILE>::launch(mdl::Device &device,cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) {
    typedef Blk<WIDTH,TILE> BLK;
    auto *pCudaOutput = reinterpret_cast<ppResult *>(pCudaBufOut);

    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, pHostBufIn, requestBufferCount, cudaMemcpyHostToDevice, stream));

    // The interation blocks
    auto *__restrict__ blkCuda = reinterpret_cast<BLK *>(pCudaBufIn);
    // The particle information
    auto *__restrict__ partCuda = reinterpret_cast<ppInput *>(blkCuda + nTotalInteractionBlocks);
    // The interaction block descriptors
    auto *__restrict__ wuCuda = reinterpret_cast<ppWorkUnit *>(partCuda + nTotalParticles);

    cudaMemsetAsync(pCudaBufOut,0,resultsBufferCount,stream);

    dim3 dimBlock( BLK::width, 8, 1 );
    dim3 dimGrid( nGrid, 1,1);
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

    CUDA_CHECK(cudaMemcpyAsync,(pHostBufOut, pCudaBufOut, resultsBufferCount, cudaMemcpyDeviceToHost, stream) );
}

template<class TILE>
void MessagePPPC<TILE>::finish() {
    auto *pR = reinterpret_cast<ppResult *>(pHostBufOut);

    for ( auto &w : work ) {
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
        pR += align_nP(nP);
        pkdParticleWorkDone(w.wp);
    }
    clear();
    freeQueue.enqueue(*this);
}

/*****************************************************************************\
*   Explicit instantiation for methods we need elsewhere
\*****************************************************************************/

template MessagePP::MessagePPPC(mdl::messageQueue<MessagePPPC> &free);
template MessagePC::MessagePPPC(mdl::messageQueue<MessagePPPC> &free);
template void MessagePP::finish();
template void MessagePC::finish();
template void MessagePP::launch(mdl::Device &device,cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessagePC::launch(mdl::Device &device,cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut);
template void CudaClient::flush(MessagePP *&M);
template void CudaClient::flush(MessagePC *&M);
