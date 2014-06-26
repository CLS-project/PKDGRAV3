/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include "cudautil.h"

#define SYNC_RATE 16  // Must be power of 2 <= warp size or trouble
#define SYNC_MASK (SYNC_RATE-1)
#define WARPS 4

/*
** Occupancy (theoretical):
**   Compute 2.x:  8 (blocks) * PP_THREADS (128) = 1024 (2/3 MAX)
**   Compute 3.x: 16 (blocks) * PP_THREADS (128) = 2048 (MAX!)
**
** To reduce memory, we syncthreads() and flush the particles
** results every SYNC_RATE particles => 8 seems a good choice.
**1
** Shared memory requirements
**  - Particles         32 * nSyncRate (8)                 =  256
**  - warp reduction     4 * nWarps (4) * 32 (threads)     =  512
**  - warp results      24 * nSyncRate (8)  * nWarps (4)   =  768
** TOTAL 1536 * 16 blocks = 24 KB
** Same nSyncRate=16: 2560 * 16 blocks = 40 KB

**
** Interaction lists limits.
**   Compute 2.x: 65535 * PP_THREADS (128) = 8 million
**   Compute 3.x: effectively unlimited
**
*/

/* Each thread block outputs this for each particle */
struct __align__(32) ppResult {
    float ax;
    float ay;
    float az;
    float fPot;
    float dirsum;
    float normsum;
    };

template<int n>
struct ppBLK {
    float dx[n],dy[n],dz[n];/* Offset from ilp->cx, cy, cz */
    float m[n],fourh2[n];   /* Mass & Softening */
    };

// One of these entries for each interaction block
struct ppWorkUnit {
    int nP;   // Number of particles
    int iP;   // Index of first particle
    int nI;   // Number of interactions in the block
    int iO;   // Index of the output block
    };


struct ppInput {
    float dx, dy, dz;
    float ax, ay, az;
    float fSoft2;
    float dImaga;
    };

#define NP_ALIGN (128/sizeof(ppResult))
#define NP_ALIGN_MASK (NP_ALIGN-1)

// A good number for nWarps is 4 giving 128 threads per thread block, nSyncRate=8
// Each thread block outputs ay,ay,az,fPot,dirsum,normsum for each particle
template <int nWarps,int nSyncRate>
__global__ void cudaPP(
    const ppWorkUnit * __restrict__ work,
    const ppInput * __restrict__ pPart,
    const ILP_BLK * __restrict__ blk,
    ppResult *out) {
    int i, iSync;

    int iWarp = threadIdx.x / 32;
    int iTinW = threadIdx.x % 32;


    // We do pPart[0..nPart-1]
    int nPart = work[blockIdx.x].nP;
    pPart += work[blockIdx.x].iP;

    // Our interaction block is quite simply found
    int nI = work[blockIdx.x].nI;
    blk += blockIdx.x;

    // Results for this block of interactions goes here
    out += work[blockIdx.x].iO;


    __shared__ union {
        ppInput P[nSyncRate];
        float   W[nSyncRate*sizeof(ppInput)/sizeof(float)];
        } Particles;

    __shared__ float reduce[nWarps][32];

    __shared__ float wX[nSyncRate][nWarps];
    __shared__ float wY[nSyncRate][nWarps];
    __shared__ float wZ[nSyncRate][nWarps];

    __shared__ float wPot[nSyncRate][nWarps];
    __shared__ float wDirsum[nSyncRate][nWarps];
    __shared__ float wNormsum[nSyncRate][nWarps];


    // Load the interaction. It is blocked for performance.
    float iX,iY,iZ,iM,ifourh2;
    if (threadIdx.x < nI) {
        iX = blk->dx.f[threadIdx.x];
        iY = blk->dy.f[threadIdx.x];
        iZ = blk->dz.f[threadIdx.x];
        iM = blk->m.f[threadIdx.x];
        ifourh2 = blk->fourh2.f[threadIdx.x];
        }
    for(iSync=0; iSync<nPart; iSync += nSyncRate) {
        int iEnd = nPart - iSync;
        if (iEnd > nSyncRate) iEnd=nSyncRate;
        // Preload the bucket of particles - this is a memcpy
        i = threadIdx.x;
        if (threadIdx.x < iEnd*sizeof(ppInput) / sizeof(float)) {
            Particles.W[i] = (reinterpret_cast<const float *>(pPart+iSync))[i];
            }
        if (i < iEnd ) { 
            float ax = Particles.P[i].ax;
            float ay = Particles.P[i].ay;
            float az = Particles.P[i].az;
            Particles.P[i].dImaga = ax*ax + ay*ay + az*az;
            if (Particles.P[i].dImaga > 0.0f) Particles.P[i].dImaga = rsqrtf(Particles.P[i].dImaga);
            }
        __syncthreads();

        for( i=0; i<iEnd; ++i) {
            float ax=0.0f, ay=0.0f, az=0.0f, fPot=0.0f, dirsum=0.0f, normsum=0.0f;
            if (threadIdx.x < nI) {
                float fourh2,dir,dir2,dir3;
                float dx = iX + Particles.P[i].dx;
                float dy = iY + Particles.P[i].dy;
                float dz = iZ + Particles.P[i].dz;
                float d2 = dx*dx + dy*dy + dz*dz;
                if (d2 != 0.0f ) { /* Ignore self interactions */
                    fourh2 = ifourh2;
                    if (d2 > fourh2) fourh2 = d2;
                    dir = rsqrtf(fourh2);
                    dir2 = dir*dir;
                    dir3 = dir2*dir;
                    if (d2 < fourh2) {
                        /*
                        ** This uses the Dehnen K1 kernel function now, it's fast!
                        */
                        dir2 *= d2;
                        dir2 = 1.0f - dir2;
                        dir *= 1.0f + dir2*(0.5f + dir2*(3.0f/8.0f + dir2*(45.0f/32.0f)));
                        dir3 *= 1.0f + dir2*(1.5f + dir2*(135.0f/16.0f));
                        }
                    dir3 *= -iM;
                    ax = dx * dir3;
                    ay = dy * dir3;
                    az = dz * dir3;
                    fPot = -iM*dir;
                    /*
                    ** Calculations for determining the timestep.
                    */
                    float adotai;
                    adotai = Particles.P[i].ax*ax + Particles.P[i].ay*ay + Particles.P[i].az*az;
                    if (adotai > 0.0f && d2 >= Particles.P[i].fSoft2 ) {
                        adotai *= Particles.P[i].dImaga;
                        dirsum = dir*adotai*adotai;
                        normsum = adotai*adotai;
                        }
                    }
                }
            // Horizontal add within each warp -- no sychronization required
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,ax,     &wX[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,ay,     &wY[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,az,     &wZ[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,fPot,   &wPot[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,dirsum, &wDirsum[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,normsum,&wNormsum[i][iWarp]);
            }

        __syncthreads();
        // Assuming four warps & SYNC of 8, the cache looks like this:
        // Cache: P0 P0 P0 P0 P1 P1 P1 P1 P2 ... P6 P7 P7 P7 P7
        // every set of 4 threads does another reduce. As long as the
        // number of warps is a power of 2 <= the warp size (32), we
        // can do this step without any further synchronization.

        int nOut = iEnd * nWarps; // Normally 32
        if (threadIdx.x<nOut) {
            int iP    = threadIdx.x & ~(nWarps-1); // 0,4,8,...
            int iWarp = threadIdx.x &  (nWarps-1); // 0 .. 3
            int iOut  = threadIdx.x / nWarps + iSync;
            warpReduceAndStore<float,nWarps>(&wX[0][0]+iP,iWarp,&out[iOut].ax);
            warpReduceAndStore<float,nWarps>(&wY[0][0]+iP,iWarp,&out[iOut].ay);
            warpReduceAndStore<float,nWarps>(&wZ[0][0]+iP,iWarp,&out[iOut].az);
            warpReduceAndStore<float,nWarps>(&wPot[0][0]+iP,iWarp,&out[iOut].fPot);
            warpReduceAndStore<float,nWarps>(&wDirsum[0][0]+iP,iWarp,&out[iOut].dirsum);
            warpReduceAndStore<float,nWarps>(&wNormsum[0][0]+iP,iWarp,&out[iOut].normsum);
            }
        }
    }

extern "C"
void pkdParticleWorkDone(workParticle *wp);

extern "C"
int CUDAcheckWorkPP( void *vpp, void *vwork ) {
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    ppResult *pR       = reinterpret_cast<ppResult *>(work->pHostBuf);
    int ib, ig, ip;


    for( ib=0; ib<work->ppnBuffered; ++ib) {
        workParticle *wp = work->ppWP[ib];
        PINFOOUT *pInfoOut = wp->pInfoOut;
        int nGridBlocks = (work->ppNI[ib] + 32*WARPS - 1) / (32*WARPS);
        for(ig=0; ig<nGridBlocks; ++ig) {
            for(ip=0; ip<wp->nP; ++ip) {
                pInfoOut[ip].a[0]    += pR[ip].ax;
                pInfoOut[ip].a[1]    += pR[ip].ay;
                pInfoOut[ip].a[2]    += pR[ip].az;
                pInfoOut[ip].fPot    += pR[ip].fPot;
                pInfoOut[ip].dirsum  += pR[ip].dirsum;
                pInfoOut[ip].normsum += pR[ip].normsum;
                }
            pR += wp->nP;
            }
        pkdParticleWorkDone(wp);

        }
    return 0;
    }


extern "C"
void CUDAsetupPP(void) {
    //cudaFuncSetCacheConfig(cudaPP,cudaFuncCachePreferL1);
    }

static CUDAwqNode *getNode(CUDACTX cuda) {
    CUDAwqNode *work;
    if (cuda->wqFree == NULL) return NULL;
    work = cuda->wqFree;
    cuda->wqFree = work->next;
    work->ctx = cuda;
    work->checkFcn = CUDAcheckWorkPP;
    work->next = cuda->wqCuda;
    work->ppSizeIn = 0;
    work->ppSizeOut = 0;
    work->ppnBuffered = 0;
    work->ppnBlocks = 0;
    return work;
    }


extern "C"
void CUDA_sendWorkPP(void *cudaCtx) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDAwqNode *work = cuda->nodePP;
    if (work != NULL) {
        int i, j;
        int iI=0, iP=0, iO=0;
        int nBufferOut = 0;

        // The interation blocks -- already copied to the host memory
        ILP_BLK * __restrict__ blkHost = reinterpret_cast<ILP_BLK*>(work->pHostBuf);
        ILP_BLK * __restrict__ blkCuda = reinterpret_cast<ILP_BLK*>(work->pCudaBufIn);
        
        // The interaction block descriptors
        ppWorkUnit * __restrict__ wuHost = reinterpret_cast<ppWorkUnit *>(blkHost + work->ppnBlocks);
        ppWorkUnit * __restrict__ wuCuda = reinterpret_cast<ppWorkUnit *>(blkCuda + work->ppnBlocks);

        // The particle information
        ppInput * __restrict__ partHost = reinterpret_cast<ppInput *>(wuHost + ((work->ppnBlocks+7)&~7));
        ppInput * __restrict__ partCuda = reinterpret_cast<ppInput *>(wuCuda + ((work->ppnBlocks+7)&~7));

        for( i=0; i<work->ppnBuffered; ++i) {
            int nP = work->ppWP[i]->nP;
            PINFOIN *pInfoIn = work->ppWP[i]->pInfoIn;

            int nPaligned = (nP+NP_ALIGN_MASK) & ~NP_ALIGN_MASK;
            int nInteract = work->ppNI[i];
            int nBlocks = (nInteract+ILP_PART_PER_BLK-1) / ILP_PART_PER_BLK;

            // Generate a interaction block descriptor for each block
            for(j=0; j<nBlocks; ++j) {
                wuHost->nP = nP;
                wuHost->iP = iP;
                wuHost->nI = nInteract > ILP_PART_PER_BLK ? ILP_PART_PER_BLK : nInteract;
                wuHost->iO = iO;
                iO += nP;
                nInteract -= wuHost->nI;
                ++wuHost;
                ++iI;
                nBufferOut += nP * sizeof(ppResult);
                }
            assert(nInteract==0);

            // Copy in nP particles
            for(j=0; j<nP; ++j) {
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
        assert(iI == work->ppnBlocks);

        ppResult *pCudaBufOut = reinterpret_cast<ppResult *>(work->pCudaBufOut);
        int nBufferIn = reinterpret_cast<char *>(partHost) - reinterpret_cast<char *>(work->pHostBuf);
        CUDA_CHECK(cudaMemcpyAsync,(blkCuda, blkHost, nBufferIn, cudaMemcpyHostToDevice, work->stream));
        dim3 dimBlock( 32*WARPS, 1 );
        dim3 dimGrid( iI, 1,1);
        cudaPP<WARPS,SYNC_RATE><<<dimGrid, dimBlock, 0, work->stream>>>(wuCuda,partCuda,blkCuda,pCudaBufOut );

        CUDA_CHECK(cudaMemcpyAsync,(blkHost, work->pCudaBufOut, nBufferOut, cudaMemcpyDeviceToHost, work->stream));
        CUDA_CHECK(cudaEventRecord,(work->event,work->stream));

        work->next = cuda->wqCuda;
        cuda->wqCuda = work;

        cuda->nodePP = NULL;
        }
    }

extern "C"
int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDAwqNode *work = cuda->nodePP;
    int nP = wp->nP;
    int nPaligned = (nP+NP_ALIGN_MASK) & ~NP_ALIGN_MASK;
    int nBlocks = tile->lstTile.nBlocks + (tile->lstTile.nInLast?1:0);
    int nInteract = tile->lstTile.nBlocks*ILP_PART_PER_BLK + tile->lstTile.nInLast;
    int nBytesIn = nPaligned * sizeof(ppInput) + nBlocks*sizeof(ILP_BLK) + nBlocks*sizeof(ppWorkUnit);
    int i;

//    int nGridBlocks = (nInteract+(32*WARPS)-1) / (32*WARPS); // number of ppBLK[]

    // Figure out the total amount of space we need, and see if there is enough
    // in the current block. If not, send it to the GPU and get another.
    // Space: nPaligned * 4 (float) * 7 (coordinates)
    //        nBlocks * sizeof(ILP_BLK)
    //
    if (work!=NULL && work->ppSizeIn + nBytesIn + 8*sizeof(ppWorkUnit) > cuda->inCudaBufSize) {
        CUDA_sendWorkPP(cuda);
        work = NULL;
        assert(cuda->nodePP==NULL);
        }

    // If we don't have a PP work element, try to grab one
    if (work==NULL) {
        cuda->nodePP = work = getNode(cuda);
        if (work==NULL) return 0;
        work->ppSizeIn = 0;
        work->ppSizeOut = 0;
        work->ppnBuffered = 0;
        work->ppnBlocks = 0;
        }
    if (work==NULL) return 0;

    // Copy in the interactions. The ILP tiles can then be freed/reused.
    ILP_BLK * __restrict__ blk = reinterpret_cast<ILP_BLK*>(work->pHostBuf);
    for(i=0; i<nBlocks; ++i) blk[work->ppnBlocks++] = tile->blk[i];

    work->ppSizeIn += nBytesIn;
    work->ppNI[work->ppnBuffered] = nInteract;
    work->ppWP[work->ppnBuffered] = wp;
    ++wp->nRefs;

    if ( ++work->ppnBuffered == CUDA_PP_MAX_BUFFERED) CUDA_sendWorkPP(cuda);

    return 1;
    }
