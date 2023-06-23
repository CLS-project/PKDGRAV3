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
inline __device__ float mask_mov(float src, bool p, float a) { return p ? a : src; }

#include "basetype.h"
#include "cudapppc.h"
#include "gravity/pp.h"
#include "gravity/pc.h"
#include "gravity/ilp.h"
#include "gravity/ilc.h"
#include "cudautil.h"
#include "cuda/reduce.h"
#include "SPH/SPHOptions.h"
#include <queue>

#define WIDTH 32

__constant__ SPHOptionsGPU SPHoptions;

/// Evaluate a single P-P interaction
template<bool bGravStep>
__device__ __forceinline__ auto evalInteraction(const gpu::ppInput pp,const gpu::Blk<WIDTH,ilpTile> *__restrict__ blk,int i) {
    return EvalPP<float,bool>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
                              blk->dx[i], blk->dy[i], blk->dz[i], blk->fourh2[i], blk->m[i],
                              pp.ax, pp.ay, pp.az, pp.dImaga);
}

template<bool bGravStep>
__device__ __forceinline__ auto evalInteraction(const gpu::denInput pp,const gpu::denBlk<WIDTH> *__restrict__ blk,int i) {
    return EvalDensity<float,bool>(pp.dx, pp.dy, pp.dz, pp.fBall, pp.iMat,
                                   blk->dx[i], blk->dy[i], blk->dz[i], blk->m[i], blk->iMat[i],
                                   SPHoptions.kernelType, SPHoptions.doInterfaceCorrection);
}

template<bool bGravStep>
__device__ __forceinline__ auto evalInteraction(const gpu::denCorrInput pp,const gpu::denCorrBlk<WIDTH> *__restrict__ blk,int i) {
    return EvalDensityCorrection<float,bool>(pp.dx, pp.dy, pp.dz, pp.fBall,
            blk->dx[i], blk->dy[i], blk->dz[i], blk->T[i], blk->P[i], blk->expImb2[i],
            SPHoptions.kernelType);
}

template<bool bGravStep>
__device__ __forceinline__ auto evalInteraction(const gpu::sphForceInput pp,const gpu::sphForceBlk<WIDTH> *__restrict__ blk,int i) {
    return EvalSPHForces<float,bool>(pp.dx, pp.dy, pp.dz, pp.fBall, pp.Omega, pp.vx, pp.vy, pp.vz, pp.rho, pp.P, pp.c,
                                     blk->dx[i], blk->dy[i], blk->dz[i], blk->m[i], blk->fBall[i], blk->Omega[i], blk->vx[i],
                                     blk->vy[i], blk->vz[i], blk->rho[i], blk->P[i], blk->c[i], blk->rung[i],
                                     SPHoptions.kernelType, SPHoptions.epsilon, SPHoptions.alpha, SPHoptions.beta,
                                     SPHoptions.EtaCourant, SPHoptions.a, SPHoptions.H, SPHoptions.useIsentropic);
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
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.ax,&out.ax);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.ay,&out.ay);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.az,&out.az);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.pot,&out.fPot);
    if (bGravStep) {
        warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.ir,&out.dirsum);
        warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.norm,&out.normsum);
    }
}

// Reduction of a density interaction
template <bool bGravStep,class RESULT>
__device__ __forceinline__ void reduceInteraction(gpu::denResult &out, RESULT result) {
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.rho,&out.rho);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.drhodfball,&out.drhodfball);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.nden,&out.nden);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.dndendfball,&out.dndendfball);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.nSmooth,&out.nSmooth);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.imbalanceX,&out.imbalanceX);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.imbalanceY,&out.imbalanceY);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.imbalanceZ,&out.imbalanceZ);
}

// Reduction of a density correction interaction
template <bool bGravStep,class RESULT>
__device__ __forceinline__ void reduceInteraction(gpu::denCorrResult &out, RESULT result) {
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.corrT,&out.corrT);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.corrP,&out.corrP);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.corr,&out.corr);
}

// Reduction of an SPH force interaction
template <bool bGravStep,class RESULT>
__device__ __forceinline__ void reduceInteraction(gpu::sphForceResult &out, RESULT result) {
    if (threadIdx.x==0) {
        if (out.dtEst == 0.0f) {
            atomicAdd(&out.dtEst,1e14f);
        }
    }
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.uDot,&out.uDot);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.ax,&out.ax);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.ay,&out.ay);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.az,&out.az);
    warpReduceAndStoreAtomicAdd<float,32>(threadIdx.x,result.divv,&out.divv);
    warpReduceAndStoreAtomicMin<float,32>(threadIdx.x,result.dtEst,&out.dtEst);
    warpReduceAndStoreAtomicMax<float,32>(threadIdx.x,result.maxRung,&out.maxRung);
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

    decltype(evalInteraction<bGravStep>(pPart[0],sblk,0)) result;
    result.zero();
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

template<int N>
void MessageDen<N>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) {
    typedef gpu::denBlk<N> BLK;
    auto *pCudaOutput = reinterpret_cast<gpu::denResult *>(pCudaBufOut);

    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, this->pHostBufIn, this->requestBufferCount, cudaMemcpyHostToDevice, stream));

    // The interation blocks
    auto *__restrict__ blkCuda = reinterpret_cast<BLK *>(pCudaBufIn);
    // The particle information
    auto *__restrict__ partCuda = reinterpret_cast<gpu::denInput *>(blkCuda + this->nTotalInteractionBlocks);
    // The interaction block descriptors
    auto *__restrict__ wuCuda = reinterpret_cast<gpu::ppWorkUnit *>(partCuda + this->nTotalParticles);

    cudaMemsetAsync(pCudaBufOut,0,this->resultsBufferCount,stream);

    dim3 dimBlock( N, 8, 1 );
    dim3 dimGrid( this->nGrid, 1,1);
    cudaInteract<false>
    <<<dimGrid, dimBlock, sizeof(BLK), stream>>>
    (wuCuda,partCuda,blkCuda,pCudaOutput );

    CUDA_CHECK(cudaMemcpyAsync,(this->pHostBufOut, pCudaBufOut, this->resultsBufferCount, cudaMemcpyDeviceToHost, stream) );
}

template<int N>
void MessageDen<N>::finish() {
    auto *pR = reinterpret_cast<gpu::denResult *>(this->pHostBufOut);

    for ( auto &w : this->work ) {
        auto nP = w.wp->nP;
        auto *pInfoOut = w.wp->pInfoOut;
        for (auto ip=0; ip<nP; ++ip) {
            pInfoOut[ip].rho         += pR[ip].rho;
            pInfoOut[ip].drhodfball  += pR[ip].drhodfball;
            pInfoOut[ip].nden        += pR[ip].nden;
            pInfoOut[ip].dndendfball += pR[ip].dndendfball;
            pInfoOut[ip].nSmooth     += pR[ip].nSmooth;
            pInfoOut[ip].imbalanceX  += pR[ip].imbalanceX;
            pInfoOut[ip].imbalanceY  += pR[ip].imbalanceY;
            pInfoOut[ip].imbalanceZ  += pR[ip].imbalanceZ;
        }
        pR += this->align_nP(nP);
        // pkdParticleWorkDone(w.wp);
        this->wps->push(w.wp);
    }
    this->clear();
    freeQueue.enqueue(*this);
}

template<int N>
void MessageDenCorr<N>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) {
    typedef gpu::denCorrBlk<N> BLK;
    auto *pCudaOutput = reinterpret_cast<gpu::denCorrResult *>(pCudaBufOut);

    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, this->pHostBufIn, this->requestBufferCount, cudaMemcpyHostToDevice, stream));

    // The interation blocks
    auto *__restrict__ blkCuda = reinterpret_cast<BLK *>(pCudaBufIn);
    // The particle information
    auto *__restrict__ partCuda = reinterpret_cast<gpu::denCorrInput *>(blkCuda + this->nTotalInteractionBlocks);
    // The interaction block descriptors
    auto *__restrict__ wuCuda = reinterpret_cast<gpu::ppWorkUnit *>(partCuda + this->nTotalParticles);

    cudaMemsetAsync(pCudaBufOut,0,this->resultsBufferCount,stream);

    dim3 dimBlock( N, 8, 1 );
    dim3 dimGrid( this->nGrid, 1,1);
    cudaInteract<false>
    <<<dimGrid, dimBlock, sizeof(BLK), stream>>>
    (wuCuda,partCuda,blkCuda,pCudaOutput );

    CUDA_CHECK(cudaMemcpyAsync,(this->pHostBufOut, pCudaBufOut, this->resultsBufferCount, cudaMemcpyDeviceToHost, stream) );
}

template<int N>
void MessageDenCorr<N>::finish() {
    auto *pR = reinterpret_cast<gpu::denCorrResult *>(this->pHostBufOut);

    for ( auto &w : this->work ) {
        auto nP = w.wp->nP;
        auto *pInfoOut = w.wp->pInfoOut;
        for (auto ip=0; ip<nP; ++ip) {
            pInfoOut[ip].corrT += pR[ip].corrT;
            pInfoOut[ip].corrP += pR[ip].corrP;
            pInfoOut[ip].corr  += pR[ip].corr;
        }
        pR += this->align_nP(nP);
        pkdParticleWorkDone(w.wp);
    }
    this->clear();
    freeQueue.enqueue(*this);
}

template<int N>
void MessageSPHForce<N>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) {
    typedef gpu::sphForceBlk<N> BLK;
    auto *pCudaOutput = reinterpret_cast<gpu::sphForceResult *>(pCudaBufOut);

    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, this->pHostBufIn, this->requestBufferCount, cudaMemcpyHostToDevice, stream));

    // The interation blocks
    auto *__restrict__ blkCuda = reinterpret_cast<BLK *>(pCudaBufIn);
    // The particle information
    auto *__restrict__ partCuda = reinterpret_cast<gpu::sphForceInput *>(blkCuda + this->nTotalInteractionBlocks);
    // The interaction block descriptors
    auto *__restrict__ wuCuda = reinterpret_cast<gpu::ppWorkUnit *>(partCuda + this->nTotalParticles);

    cudaMemsetAsync(pCudaBufOut,0,this->resultsBufferCount,stream);

    dim3 dimBlock( N, 8, 1 );
    dim3 dimGrid( this->nGrid, 1,1);
    cudaInteract<false>
    <<<dimGrid, dimBlock, sizeof(BLK), stream>>>
    (wuCuda,partCuda,blkCuda,pCudaOutput );

    CUDA_CHECK(cudaMemcpyAsync,(this->pHostBufOut, pCudaBufOut, this->resultsBufferCount, cudaMemcpyDeviceToHost, stream) );
}

template<int N>
void MessageSPHForce<N>::finish() {
    auto *pR = reinterpret_cast<gpu::sphForceResult *>(this->pHostBufOut);

    for ( auto &w : this->work ) {
        auto nP = w.wp->nP;
        auto *pInfoOut = w.wp->pInfoOut;
        for (auto ip=0; ip<nP; ++ip) {
            pInfoOut[ip].uDot   += pR[ip].uDot;
            pInfoOut[ip].a[0]   += pR[ip].ax;
            pInfoOut[ip].a[1]   += pR[ip].ay;
            pInfoOut[ip].a[2]   += pR[ip].az;
            pInfoOut[ip].divv   += pR[ip].divv;
            pInfoOut[ip].dtEst   = std::min(pInfoOut[ip].dtEst,pR[ip].dtEst);
            pInfoOut[ip].maxRung = std::max(pInfoOut[ip].maxRung,pR[ip].maxRung);
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
template void MessageDen<32>::finish();
template void MessageDenCorr<32>::finish();
template void MessageSPHForce<32>::finish();
template void MessagePP::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessagePC::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessageDen<32>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessageDenCorr<32>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessageSPHForce<32>::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut);

/*
** This has to live here currently, but will be moved out later.
*/
void CudaClient::setupSPHOptions(SPHOptionsGPU *const SPHoptions) {
    mdl::cudaMessageQueue wait;
    for (auto i=0; i<cuda.numDevices(); ++i) {
        auto m = new MessageSPHOptionsSetup(SPHoptions,i);
        cuda.enqueue(*m,wait);
    }
    for (auto i=0; i<cuda.numDevices(); ++i) {
        auto &m = wait.wait();
        delete &m;
    }
}

MessageSPHOptionsSetup::MessageSPHOptionsSetup(SPHOptionsGPU *const SPHoptions, int iDevice)
    : mdl::cudaMessage(iDevice), SPHoptionsIn(SPHoptions) {}

void MessageSPHOptionsSetup::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) {
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (SPHoptions, SPHoptionsIn, sizeof(SPHoptions), 0, cudaMemcpyHostToDevice, stream));
}