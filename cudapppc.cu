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
#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#include <stdio.h>
#include "cudapppc.h"
#include "pp.h"
#include "pc.h"

#define SYNC_RATE 16  // Must be: 1, 2, 4, 8, 16
#define WIDTH 32

#define TB_THREADS 128
#define WARPS (TB_THREADS/32)

#define PP_WU 128
#define PC_WU 32

#include "ilp.h"
#include "ilc.h"
#include "basetype.h"

/*
** The following are the basically the same as ILP_BLK and ILC_BLK,
** but we need to be able to alter their sizes.
*/

template <int n,class TILE> struct Blk;
template <int n> struct Blk<n,ilpTile> {
    float dx[n], dy[n], dz[n];    /* Offset from ilp->cx, cy, cz */
    float m[n];             /* Mass */
    float fourh2[n];        /* Softening: calculated */
    };
template <int n> struct Blk<n,ilcTile> {
    float dx[n],dy[n],dz[n];
    float xxxx[n],xxxy[n],xxxz[n],xxyz[n],xxyy[n],yyyz[n],xyyz[n],xyyy[n],yyyy[n];
    float xxx[n],xyy[n],xxy[n],yyy[n],xxz[n],yyz[n],xyz[n];
    float xx[n],xy[n],xz[n],yy[n],yz[n];
    float x[n],y[n],z[n];
    float m[n],u[n];
    };

// One of these entries for each interaction block
struct ppWorkUnit {
    uint32_t nP;   // Number of particles
    uint32_t nI;   // Number of interactions in the block
    uint32_t iP;   // Index of first particle
    uint32_t iO;   // Index of the output block
    };


struct ppInput {
    float dx, dy, dz;
    float ax, ay, az;
    float fSoft2;
    float dImaga;
    };

/* Each thread block outputs this for each particle */
struct __align__(32) ppResult {
    float ax;
    float ay;
    float az;
    float fPot;
    float dirsum;
    float normsum;
    };

#define NP_ALIGN (128/sizeof(ppResult))
#define NP_ALIGN_MASK (NP_ALIGN-1)

/*
** Occupancy (theoretical):
**   Compute 2.x:  8 (blocks) * PP_THREADS (128) = 1024 (2/3 MAX)
**   Compute 3.x: 16 (blocks) * PP_THREADS (128) = 2048 (MAX!)
**
** To reduce memory, we syncthreads() and flush the particles
** results every SYNC_RATE particles => 8 seems a good choice.
**1
** Shared memory requirements
**  - Particles         32 * nSyncRate (16)                =  512
**  - warp reduction     4 * nWarps (4) * 32 (threads)     =  512
**  - warp results      24 * nSyncRate (16)  * nWarps (4)  = 1536
** TOTAL 2560 * 16 blocks = 40 KB
**
** nvcc -DHAVE_CONFIG_H --ptxas-options=-v -c  -I. -arch=sm_20 cudapp.cu
** ptxas info    : 11 bytes gmem, 8 bytes cmem[14]
** ptxas info    : Compiling entry function '_Z6cudaPPILi4ELi16EEvPK10ppWorkUnitPK7ppInputPK7ILP_BLKP8ppResult' for 'sm_20'
** ptxas info    : Function properties for _Z6cudaPPILi4ELi16EEvPK10ppWorkUnitPK7ppInputPK7ILP_BLKP8ppResult
**     0 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
** ptxas info    : Used 35 registers, 2560 bytes smem, 64 bytes cmem[0], 12 bytes cmem[16]
**
** Interaction lists limits.
**   Compute 2.x: 65535 * PP_THREADS (128) = 8 million
**   Compute 3.x: effectively unlimited
**
*/

// A good number for nWarps is 4 giving 128 threads per thread block, nSyncRate=8
// Each thread block outputs ay,ay,az,fPot,dirsum,normsum for each particle
template <int nWarps,int nWarpsPerWU,int nSyncRate,int bGravStep>
__global__ void cudaInteract(
    const ppWorkUnit * __restrict__ work,
    const ppInput * __restrict__ pPart,
    const Blk<WIDTH,ilpTile> * __restrict__ blk,
    ppResult *out) {
    int i, iSync;
    int iWork, iI, iWarp;

    if (nWarpsPerWU==1) {           // blockDim.z == nWarps, blockDim.y == 1, blockDim.x == 32
        iWork = blockIdx.x * nWarps + threadIdx.z; // Work and corresponding blk
        iI = threadIdx.x; // Index into blk
        iWarp = threadIdx.z;
        }
    else if (nWarps==nWarpsPerWU) { // blockDim.z == 1, blockDim.y == nWarps, blockDim.x == 32
        iWork = blockIdx.x; // Index of work and blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Index of interaction
        iWarp = threadIdx.y;
        }
    else {                          // blockDim.z == 2, blockDim.y == 2, blockDim.x == 32
        // Calculate our interaction and particle group
        iWork = blockIdx.x*blockDim.z + threadIdx.z; // Work and corresponding blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Thread working on blk
        iWarp = threadIdx.y + blockDim.y*threadIdx.z;
        }
    int iTinW = iI % 32;

    uint32_t nP = work[iWork].nP; // Number of particles
    pPart += work[iWork].iP; // First particle
    uint32_t nI = work[iWork].nI; // Number of interactions
//    blk += work[iWork].iB*blockDim.y + threadIdx.y; // blk[threadIdx.x] is our interaction
    blk += iWork*blockDim.y + threadIdx.y; // blk[threadIdx.x] is our interaction
    out += work[iWork].iO;   // Result for each particle

    __shared__ union {
        ppInput P[nWarps/nWarpsPerWU][nSyncRate];
        float   W[nWarps/nWarpsPerWU][nSyncRate*sizeof(ppInput)/sizeof(float)];
        } Particles;

    __shared__ float wX[nSyncRate][nWarps];
    __shared__ float wY[nSyncRate][nWarps];
    __shared__ float wZ[nSyncRate][nWarps];

    __shared__ float wPot[nSyncRate][nWarps];
    __shared__ float wDirsum[nSyncRate][nWarps];
    __shared__ float wNormsum[nSyncRate][nWarps];


    // Load the interaction. It is blocked for performance.
    float iX,iY,iZ,iM,ifourh2;
    if (iI < nI) {
        iX = blk->dx[threadIdx.x];
        iY = blk->dy[threadIdx.x];
        iZ = blk->dz[threadIdx.x];
        iM = blk->m[threadIdx.x];
        ifourh2 = blk->fourh2[threadIdx.x];
        }

    // Apply the particles, nSyncRate at a time
    for(iSync=0; iSync<nP; iSync += nSyncRate) {
        int iEnd = nP - iSync;
        if (iEnd > nSyncRate) iEnd=nSyncRate;
        // Preload the bucket of particles - this is a memcpy
        if (iI < iEnd*sizeof(ppInput) / sizeof(float)) {
            Particles.W[threadIdx.z][iI] = (reinterpret_cast<const float *>(pPart+iSync))[iI];
            }
        if (nWarpsPerWU>1) __syncthreads();
        if (iI < iEnd && bGravStep) { 
            float ax = Particles.P[threadIdx.z][iI].ax;
            float ay = Particles.P[threadIdx.z][iI].ay;
            float az = Particles.P[threadIdx.z][iI].az;
            ax = ax*ax + ay*ay + az*az;
            if (ax > 0.0f) ax = rsqrtf(ax);
            Particles.P[threadIdx.z][iI].dImaga = ax;
            }
        if (nWarpsPerWU>1) __syncthreads();

        for( i=0; i<iEnd; ++i) {
            float ax=0.0f, ay=0.0f, az=0.0f, fPot=0.0f, dir=0.0f, norm=0.0f;
            if (iI < nI) {
                float Px = Particles.P[threadIdx.z][i].dx;
                float Py = Particles.P[threadIdx.z][i].dy;
                float Pz = Particles.P[threadIdx.z][i].dz;
		float fSoft2 = Particles.P[threadIdx.z][i].fSoft2;
		float iax = Particles.P[threadIdx.z][i].ax;
		float iay = Particles.P[threadIdx.z][i].ay;
		float iaz = Particles.P[threadIdx.z][i].az;
		float imaga = Particles.P[threadIdx.z][i].dImaga;
		EvalPP<float,bool,true>(
		    Px, Py, Pz, fSoft2, 
		    iX, iY, iZ, ifourh2, iM,
		    ax, ay, az, fPot,
		    iax, iay, iaz, imaga,
		    dir, norm);
                }
            // Horizontal add within each warp -- no sychronization required
            warpReduceAndStore<float,32>(iTinW,ax,       &wX[i][iWarp]);
            warpReduceAndStore<float,32>(iTinW,ay,       &wY[i][iWarp]);
            warpReduceAndStore<float,32>(iTinW,az,       &wZ[i][iWarp]);
            warpReduceAndStore<float,32>(iTinW,fPot,     &wPot[i][iWarp]);
            if (bGravStep) {
                warpReduceAndStore<float,32>(iTinW,dir,  &wDirsum[i][iWarp]);
                warpReduceAndStore<float,32>(iTinW,norm, &wNormsum[i][iWarp]);
                }
            }

        if (nWarpsPerWU>1) __syncthreads();
        // Assuming four warps & SYNC of 8, the cache looks like this:
        //                0     1     2     3       4
        // Cache: 1x4x8   P0    P0    P0    P0      P1 P1 P1 P1 P2 ... P6 P7 P7 P7 P7
        // Cache: 4x1x8   P0,I0 P0,I1 P0,I2 P0,I3   P1,I0 ... 
        // Cache: 2x2x8   P0,I0 P0,I0 P0,I1 P0,I1   P1,I0
        // every set of 4 threads does another reduce. As long as the
        // number of warps is a power of 2 <= the warp size (32), we
        // can do this step without any further synchronization.
        int nOut = iEnd * nWarpsPerWU; // Normally 64
        if (iI<nOut) {
            int iP    = (iI & ~(nWarpsPerWU-1)) * nWarps/nWarpsPerWU + threadIdx.z*nWarpsPerWU; // 0,4,8,...
            int iWarp = iI &  (nWarpsPerWU-1); // 0 .. 3
            int iOut  = iI / nWarpsPerWU + iSync;
         
            warpReduceAndStore<float,nWarpsPerWU>(      &wX[0][0]+iP,iWarp,&out[iOut].ax);
            warpReduceAndStore<float,nWarpsPerWU>(      &wY[0][0]+iP,iWarp,&out[iOut].ay);
            warpReduceAndStore<float,nWarpsPerWU>(      &wZ[0][0]+iP,iWarp,&out[iOut].az);
            warpReduceAndStore<float,nWarpsPerWU>(    &wPot[0][0]+iP,iWarp,&out[iOut].fPot);
            if (bGravStep) {
                warpReduceAndStore<float,nWarpsPerWU>( &wDirsum[0][0]+iP,iWarp,&out[iOut].dirsum);
                warpReduceAndStore<float,nWarpsPerWU>(&wNormsum[0][0]+iP,iWarp,&out[iOut].normsum);
                }
            }
        }
    }

template <int nWarps,int nWarpsPerWU,int nSyncRate,int bGravStep>
__global__ void cudaInteract(
    const ppWorkUnit * __restrict__ work,
    const ppInput * __restrict__ pPart,
    const Blk<WIDTH,ilcTile> * __restrict__ blk,
    ppResult *out) {
    int i, iSync;
    int iWork, iI, iWarp;

    if (nWarpsPerWU==1) {           // blockDim.z == nWarps, blockDim.y == 1, blockDim.x == 32
        iWork = blockIdx.x * nWarps + threadIdx.z; // Work and corresponding blk
        iI = threadIdx.x; // Index into blk
        iWarp = threadIdx.z;
        }
    else if (nWarps==nWarpsPerWU) { // blockDim.z == 1, blockDim.y == nWarps, blockDim.x == 32
        iWork = blockIdx.x; // Index of work and blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Index of interaction
        iWarp = threadIdx.y;
        }
    else {
        // Calculate our interaction and particle group
        iWork = blockIdx.x*blockDim.z + threadIdx.z; // Work and corresponding blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Thread working on blk
        int iAll = iI + threadIdx.z*blockDim.y*blockDim.x;
        iWarp = iAll / 32;
        }
    int iTinW = iI % 32;

    int nP = work[iWork].nP; // Number of particles
    pPart += work[iWork].iP; // First particle
    int nI = work[iWork].nI; // Number of interactions
//    blk += work[iWork].iB*blockDim.y + threadIdx.y; // blk[threadIdx.x] is our interaction
    blk += iWork*blockDim.y + threadIdx.y; // blk[threadIdx.x] is our interaction
    out += work[iWork].iO;   // Result for each particle

    __shared__ union {
        ppInput P[nWarps/nWarpsPerWU][nSyncRate];
        float   W[nWarps/nWarpsPerWU][nSyncRate*sizeof(ppInput)/sizeof(float)];
        } Particles;

    __shared__ float wX[nSyncRate][nWarps];
    __shared__ float wY[nSyncRate][nWarps];
    __shared__ float wZ[nSyncRate][nWarps];

    __shared__ float wPot[nSyncRate][nWarps];
    __shared__ float wDirsum[nSyncRate][nWarps];
    __shared__ float wNormsum[nSyncRate][nWarps];


    // Load the interaction. It is blocked for performance.
    float Idx,Idy,Idz;
    float Ixxxx,Ixxxy,Ixxxz,Ixxyz,Ixxyy,Iyyyz,Ixyyz,Ixyyy,Iyyyy;
    float Ixxx,Ixyy,Ixxy,Iyyy,Ixxz,Iyyz,Ixyz;
    float Ixx,Ixy,Ixz,Iyy,Iyz;
#ifdef USE_DIAPOLE
    float Ix,Iy,Iz;
#endif
    float Im,Iu;
    if (iI < nI) {
        Idx = blk->dx[threadIdx.x];
        Idy = blk->dy[threadIdx.x];
        Idz = blk->dz[threadIdx.x];
        Ixxxx = blk->xxxx[threadIdx.x];
        Ixxxy = blk->xxxy[threadIdx.x];
        Ixxxz = blk->xxxz[threadIdx.x];
        Ixxyz = blk->xxyz[threadIdx.x];
        Ixxyy = blk->xxyy[threadIdx.x];
        Iyyyz = blk->yyyz[threadIdx.x];
        Ixyyz = blk->xyyz[threadIdx.x];
        Ixyyy = blk->xyyy[threadIdx.x];
        Iyyyy = blk->yyyy[threadIdx.x];
        Ixxx = blk->xxx[threadIdx.x];
        Ixyy = blk->xyy[threadIdx.x];
        Ixxy = blk->xxy[threadIdx.x];
        Iyyy = blk->yyy[threadIdx.x];
        Ixxz = blk->xxz[threadIdx.x];
        Iyyz = blk->yyz[threadIdx.x];
        Ixyz = blk->xyz[threadIdx.x];
        Ixx = blk->xx[threadIdx.x];
        Ixy = blk->xy[threadIdx.x];
        Ixz = blk->xz[threadIdx.x];
        Iyy = blk->yy[threadIdx.x];
        Iyz = blk->yz[threadIdx.x];
#ifdef USE_DIAPOLE
        Ix = blk->x[threadIdx.x];
        Iy = blk->y[threadIdx.x];
        Iz = blk->z[threadIdx.x];
#endif
        Im = blk->m[threadIdx.x];
        Iu = blk->u[threadIdx.x];
        }
    for(iSync=0; iSync<nP; iSync += nSyncRate) {
        int iEnd = nP - iSync;
        if (iEnd > nSyncRate) iEnd=nSyncRate;
        // Preload the bucket of particles - this is a memcpy
        if (iI < iEnd*sizeof(ppInput) / sizeof(float)) {
            Particles.W[threadIdx.z][iI] = (reinterpret_cast<const float *>(pPart+iSync))[iI];
            }
        if (nWarpsPerWU>1) __syncthreads();
        if (iI < iEnd && bGravStep) { 
            float ax = Particles.P[threadIdx.z][iI].ax;
            float ay = Particles.P[threadIdx.z][iI].ay;
            float az = Particles.P[threadIdx.z][iI].az;
            Particles.P[threadIdx.z][iI].dImaga = ax*ax + ay*ay + az*az;
            if (Particles.P[threadIdx.z][iI].dImaga > 0.0f)
                Particles.P[threadIdx.z][iI].dImaga = rsqrtf(Particles.P[threadIdx.z][iI].dImaga);
            }
        if (nWarpsPerWU>1) __syncthreads();

        for( i=0; i<iEnd; ++i) {
            float ax=0.0f, ay=0.0f, az=0.0f, fPot=0.0f, dirsum=0.0f, norm=0.0f;
            if (iI < nI) {
                float Pdx = Particles.P[threadIdx.z][i].dx;
                float Pdy = Particles.P[threadIdx.z][i].dy;
                float Pdz = Particles.P[threadIdx.z][i].dz;
		float Pax = Particles.P[threadIdx.z][i].ax;
		float Pay = Particles.P[threadIdx.z][i].ay;
		float Paz = Particles.P[threadIdx.z][i].az;
		float fSoft2 = Particles.P[threadIdx.z][i].fSoft2;
		float Pimaga = Particles.P[threadIdx.z][i].dImaga;
		EvalPC<float,bool,true>(
		    Pdx, Pdy, Pdz,fSoft2,
		    Idx, Idy, Idz, Im, Iu,
		    Ixxxx, Ixxxy, Ixxxz, Ixxyz, Ixxyy, Iyyyz, Ixyyz, Ixyyy, Iyyyy,
		    Ixxx, Ixyy, Ixxy, Iyyy, Ixxz, Iyyz, Ixyz, Ixx, Ixy, Ixz, Iyy, Iyz,
#ifdef USE_DIAPOLE
		    Ix, Iy, Iz,
#endif
		    ax, ay, az, fPot,
		    Pax, Pay, Paz, Pimaga,
		    dirsum, norm);
                }
            // Horizontal add within each warp -- no sychronization required
            warpReduceAndStore<float,32>(iTinW,ax,     &wX[i][iWarp]);
            warpReduceAndStore<float,32>(iTinW,ay,     &wY[i][iWarp]);
            warpReduceAndStore<float,32>(iTinW,az,     &wZ[i][iWarp]);
            warpReduceAndStore<float,32>(iTinW,fPot,   &wPot[i][iWarp]);
            if (bGravStep) {
                warpReduceAndStore<float,32>(iTinW,dirsum, &wDirsum[i][iWarp]);
                warpReduceAndStore<float,32>(iTinW,norm,&wNormsum[i][iWarp]);
                }
	}

        if (nWarpsPerWU>1) __syncthreads();
        // Assuming four warps & SYNC of 8, the cache looks like this:
        //                0     1     2     3       4
        // Cache: 1x4x8   P0    P0    P0    P0      P1 P1 P1 P1 P2 ... P6 P7 P7 P7 P7
        // Cache: 4x1x8   P0,I0 P0,I1 P0,I2 P0,I3   P1,I0 ... 
        // Cache: 2x2x8   P0,I0 P0,I0 P0,I1 P0,I1   P1,I0
        // every set of 4 threads does another reduce. As long as the
        // number of warps is a power of 2 <= the warp size (32), we
        // can do this step without any further synchronization.
        int nOut = iEnd * nWarpsPerWU; // Normally 64
        if (iI<nOut) {
            int iP    = (iI & ~(nWarpsPerWU-1)) * nWarps/nWarpsPerWU + threadIdx.z; // 0,4,8,...
            int iWarp = iI &  (nWarpsPerWU-1); // 0 .. 3
            int iOut  = iI / nWarpsPerWU + iSync;
         
            warpReduceAndStore<float,nWarpsPerWU>(      &wX[0][iP],iWarp,&out[iOut].ax);
            warpReduceAndStore<float,nWarpsPerWU>(      &wY[0][iP],iWarp,&out[iOut].ay);
            warpReduceAndStore<float,nWarpsPerWU>(      &wZ[0][iP],iWarp,&out[iOut].az);
            warpReduceAndStore<float,nWarpsPerWU>(    &wPot[0][iP],iWarp,&out[iOut].fPot);
            if (bGravStep) {
                warpReduceAndStore<float,nWarpsPerWU>( &wDirsum[0][iP],iWarp,&out[iOut].dirsum);
                warpReduceAndStore<float,nWarpsPerWU>(&wNormsum[0][iP],iWarp,&out[iOut].normsum);
                }
            }
        }
    }

extern "C"
void pkdParticleWorkDone(workParticle *wp);

/*****************************************************************************\
*   CudaClient interface (new!)
\*****************************************************************************/

extern "C"
int CudaClientQueuePP(void *vcudaClient, workParticle *work, struct ilpTile *tile, int bGravStep) {
    auto cuda = reinterpret_cast<CudaClient *>(vcudaClient);
    return cuda->queuePP(work,tile,bGravStep);
    }
int CudaClient::queuePP(workParticle *work, ilpTile *tile, bool bGravStep) {
    return queue(pp,freePP,work,tile,bGravStep);
    }

extern "C"
int CudaClientQueuePC(void *vcudaClient, workParticle *work, struct ilcTile *tile, int bGravStep) {
    auto cuda = reinterpret_cast<CudaClient *>(vcudaClient);
    return cuda->queuePC(work,tile,bGravStep);
    }
int CudaClient::queuePC(workParticle *work, ilcTile *tile, bool bGravStep) {
    return queue(pc,freePC,work,tile,bGravStep);
    }

template<class MESSAGE,class QUEUE,class TILE>
int CudaClient::queue(MESSAGE * &M,QUEUE &Q, workParticle *work, TILE *tile, bool bGravStep) {
    if (M) { // If we are in the middle of building data for a kernel
	if (M->queue(work,tile,bGravStep)) return work->nP; // Successfully queued
	flush(M); // The buffer is full, so send it
	}
    if (Q.empty()) return 0; // No buffers so the CPU has to do this part
    M = & Q.dequeue();
    if (M->queue(work,tile,bGravStep)) return work->nP; // Successfully queued
    return 0; // Not sure how this would happen, but okay.
    }

template<class MESSAGE>
void CudaClient::flush(MESSAGE * &M) {
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
*   2. An array of interaction block descriptors
*   3. An array of particles
*
*   When an interaction list is queued, we:
*   1. Make sure that the additional interaction blocks and the pending
*      descriptors and particles will fit, returning false if they won't.
*   2. The individual blocks are copied into the output buffer.
*   3. The workParticle structure is saved for later use
\*****************************************************************************/

template<class TILE,int nIntPerTB, int nIntPerWU>
MessagePPPC<TILE,nIntPerTB,nIntPerWU>::MessagePPPC(mdl::messageQueue<MessagePPPC> &freeQueue)
    : freeQueue(freeQueue), requestBufferCount(0), resultsBufferCount(0), nInteractionBlocks(0) {
    work.reserve(CUDA_WP_MAX_BUFFERED);
    }

// This function empties the message for subsequent queue requests
template<class TILE,int nIntPerTB, int nIntPerWU>
void MessagePPPC<TILE,nIntPerTB,nIntPerWU>::clear() {
    work.clear();
    requestBufferCount = resultsBufferCount = 0;
    nInteractionBlocks = 0;
    }

template<int n,class TILE>
int copyBLKs2(Blk<n,TILE> *out, TILE *in,const int nIlp) {
    assert(n==ILP_PART_PER_BLK);
    int i, nBlk = (nIlp+n-1) / n;
    for(i=0; i<nBlk; ++i) memcpy(&out[i],&in->blk[i],sizeof(out[i]));
    return nBlk;
    }

// Add the interactions to this list
template<class TILE,int nIntPerTB, int nIntPerWU>
bool MessagePPPC<TILE,nIntPerTB,nIntPerWU>::queue(workParticle *wp, TILE *tile, bool bGravStep) {
    if (work.size() == CUDA_WP_MAX_BUFFERED) return false; // Too many work packages
    this->bGravStep = bGravStep;

    typedef Blk<WIDTH,TILE> BLK;

    const int nBlkPerWU = nIntPerWU / WIDTH;
    const int nP = wp->nP;
    const int nPaligned = (nP+NP_ALIGN_MASK) & ~NP_ALIGN_MASK;
    const int nBlocks = tile->lstTile.nBlocks + (tile->lstTile.nInLast?1:0);
    const int nBlocksAligned = (nBlocks + nBlkPerWU - 1) & ~(nBlkPerWU - 1);
    const int nInteract = tile->lstTile.nBlocks*ILP_PART_PER_BLK + tile->lstTile.nInLast;
    const int nWork = nBlocksAligned / nBlkPerWU;
    const int nBytesIn = nPaligned * sizeof(ppInput) + nBlocksAligned*sizeof(BLK) + nWork*sizeof(ppWorkUnit);
    const int nBytesOut = nP * sizeof(ppResult) * nWork;
    assert(nWork*nBlkPerWU == nBlocksAligned);

    if ( requestBufferCount + nBytesIn + 8*sizeof(ppWorkUnit) > requestBufferSize || resultsBufferCount + nBytesOut > resultsBufferSize) return false;
    requestBufferCount += nBytesIn;
    resultsBufferCount += nBytesOut;

    // Copy in the interactions. The ILP tiles can then be freed/reused.
    auto blk = reinterpret_cast<Blk<WIDTH,TILE> *>(pHostBufIn);

    copyBLKs2(blk+nInteractionBlocks,tile,nInteract);
    nInteractionBlocks += nBlocksAligned;
    work.emplace_back(wp,nInteract);
    ++wp->nRefs;

    return true;
    }

// Final preparation before sending this message to the GPU thread.
// We need to setup the descriptors and add the particles.
template<class TILE,int nIntPerTB, int nIntPerWU>
MessagePPPC<TILE,nIntPerTB,nIntPerWU> & MessagePPPC<TILE,nIntPerTB,nIntPerWU>::prepare() {
    typedef Blk<WIDTH,TILE> BLK;
    int iI=0, iP=0, iO=0;
    const int nBlkPer = nIntPerWU / WIDTH;
    const int nWork = nInteractionBlocks / nBlkPer;

    // The interation blocks -- already copied to the host memory
    auto * __restrict__ blkHost = reinterpret_cast<BLK*>(pHostBufIn);

    // The interaction block descriptors
    auto * __restrict__ wuHost = reinterpret_cast<ppWorkUnit *>(blkHost + nInteractionBlocks);

    // The particle information
    auto * __restrict__ partHost = reinterpret_cast<ppInput *>(wuHost + ((nWork+7)&~7));

    size_t nOutputBytes = 0;
    for ( auto &w : work ) {
        auto nP = w.wp->nP;
        auto *pInfoIn = w.wp->pInfoIn;
        int nPaligned = (nP+NP_ALIGN_MASK) & ~NP_ALIGN_MASK;
        int nInteract = w.nInteractions;
        int nBlocks = (nInteract+nIntPerWU-1) / nIntPerWU;

        // Generate a interaction block descriptor for each block
        for(auto j=0; j<nBlocks; ++j) {
            wuHost->nP = nP;
            wuHost->iP = iP;
            wuHost->nI = nInteract > nIntPerWU ? nIntPerWU : nInteract;
            wuHost->iO = iO;
	    //wuHost->iB = iI;
            iO += nP;
            nInteract -= wuHost->nI;
            ++wuHost;
            ++iI;
            nOutputBytes += nP * sizeof(ppResult);
            }
        assert(nInteract==0);
        // Copy in nP particles
        for(auto j=0; j<nP; ++j) {
            partHost[j].dx =  pInfoIn[j].r[0];
            partHost[j].dy =  pInfoIn[j].r[1];
            partHost[j].dz =  pInfoIn[j].r[2];
            partHost[j].ax =  pInfoIn[j].a[0];
            partHost[j].ay =  pInfoIn[j].a[1];
            partHost[j].az =  pInfoIn[j].a[2];
            partHost[j].fSoft2 = pInfoIn[j].fSmooth2;
            /*partHost[j].dImaga = 0;*/
            }
        partHost += nPaligned;
        iP += nPaligned;
	}
    assert(resultsBufferCount == nOutputBytes);
    assert(iI == nInteractionBlocks/nBlkPer);
    /* Pad the work out so all work units have valid data */
    const int nWUPerTB = nIntPerTB/nIntPerWU;
    while((iI&(nWUPerTB-1)) != 0) {
        wuHost->nP = 0;
        wuHost->iP = 0;
        wuHost->nI = 0;
        wuHost->iO = iO;
        ++wuHost;
        ++iI;
        }

    assert((iI & (nWUPerTB-1)) == 0);
    nGrid = iI/nWUPerTB;
    requestBufferCount = reinterpret_cast<char *>(partHost) - reinterpret_cast<char *>(pHostBufIn);
    assert(requestBufferCount <= requestBufferSize);
    assert(resultsBufferCount <= resultsBufferSize);

    return *this;
    }

template<class TILE,int nIntPerTB, int nIntPerWU>
void MessagePPPC<TILE,nIntPerTB,nIntPerWU>::launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) {
    typedef Blk<WIDTH,TILE> BLK;
    const int nBlkPer = nIntPerWU / WIDTH;
    const int nWork = nInteractionBlocks / nBlkPer;

    // The interation blocks -- already copied to the host memory
    auto * __restrict__ blkCuda = reinterpret_cast<BLK*>(pCudaBufIn);

    // The interaction block descriptors
    auto * __restrict__ wuCuda = reinterpret_cast<ppWorkUnit *>(blkCuda + nInteractionBlocks);

    // The particle information
    auto * __restrict__ partCuda = reinterpret_cast<ppInput *>(wuCuda + ((nWork+7)&~7));

    auto *pCudaOutput = reinterpret_cast<ppResult *>(pCudaBufOut);

    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, pHostBufIn, requestBufferCount, cudaMemcpyHostToDevice, stream));
    dim3 dimBlock( WIDTH, nIntPerWU/WIDTH, nIntPerTB/nIntPerWU );
    dim3 dimGrid( nGrid, 1,1);

    if (bGravStep) {
        cudaInteract<WARPS,nIntPerWU/32,SYNC_RATE*nIntPerWU/nIntPerTB,1>
            <<<dimGrid, dimBlock, 0, stream>>>
            (wuCuda,partCuda,blkCuda,pCudaOutput );
        }
    else {
        cudaInteract<WARPS,nIntPerWU/32,SYNC_RATE*nIntPerWU/nIntPerTB,0>
            <<<dimGrid, dimBlock, 0, stream>>>
            (wuCuda,partCuda,blkCuda,pCudaOutput );
        }

    CUDA_CHECK(cudaMemcpyAsync,(pHostBufOut, pCudaBufOut, resultsBufferCount, cudaMemcpyDeviceToHost, stream) );
    }

template<class TILE,int nIntPerTB, int nIntPerWU>
void MessagePPPC<TILE,nIntPerTB,nIntPerWU>::finish() {
    auto *pR = reinterpret_cast<ppResult *>(pHostBufOut);

    for ( auto &w : work ) {
        auto nP = w.wp->nP;
        auto *pInfoOut = w.wp->pInfoOut;

        int nWork = (w.nInteractions + nIntPerWU - 1) / nIntPerWU;
        for(auto iw=0; iw<nWork; ++iw) {
            for(auto ip=0; ip<nP; ++ip) {
                pInfoOut[ip].a[0]    += pR[ip].ax;
                pInfoOut[ip].a[1]    += pR[ip].ay;
                pInfoOut[ip].a[2]    += pR[ip].az;
                pInfoOut[ip].fPot    += pR[ip].fPot;
                pInfoOut[ip].dirsum  += pR[ip].dirsum;
                pInfoOut[ip].normsum += pR[ip].normsum;
                }
            pR += nP;
            }
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
template void MessagePP::launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut);
template void MessagePC::launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut);
template void CudaClient::flush(MessagePP * &M);
template void CudaClient::flush(MessagePC * &M);
