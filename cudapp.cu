/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include "cudautil.h"

#define MAX_BUCKET 32 // Larger buckets just become a different kernel
#define SYNC_RATE 8  // Must be power of 2 <= warp size or trouble
#define SYNC_MASK (SYNC_RATE-1)
#define WARPS 4

/*
** Occupancy (theoretical):
**   Compute 2.x:  8 (blocks) * PP_THREADS (128) = 1024 (2/3 MAX)
**   Compute 3.x: 16 (blocks) * PP_THREADS (128) = 2048 (MAX!)
**
** To reduce memory, we syncthreads() and flush the particles
** results every SYNC_RATE particles => 8 seems a good choice.
**
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

/* Each thread block outputs this for each nSyncRate particles */
template <int nSyncRate>
struct __align__(128) ppResults {
    float ax[nSyncRate];
    float ay[nSyncRate];
    float az[nSyncRate];
    float fPot[nSyncRate];
    float dirsum[nSyncRate];
    float normsum[nSyncRate];
    };

// A good number for nWarps is 4 giving 128 threads per thread block
// nMaxBucket=32 and nSyncRate=8
// Each thread block outputs ay,ay,az,fPot,dirsum,normsum for each particle
template <int nWarps,int nMaxBucket,int nSyncRate>
__global__ void cudaPP(
    int nPart,
    const float * __restrict__ pX, const float * __restrict__ pY, const float * __restrict__ pZ, // Particles
    const float * __restrict__ pAX, const float * __restrict__ pAY, const float * __restrict__ pAZ, // Particles
    const float * __restrict__ pSoft2,
    int nInteract, const ILP_BLK * __restrict__ blk,ppResults<nSyncRate> *out) { // interactions
    int iWarp = threadIdx.x / 32;
    int iTinW = threadIdx.x % 32;
    int nParticleBlocks = (nPart + nSyncRate - 1) / nSyncRate;
    int iOutBlock = nParticleBlocks * blockIdx.x;
    int i, iSync;

    __shared__ float X[nSyncRate];
    __shared__ float Y[nSyncRate];
    __shared__ float Z[nSyncRate];
    __shared__ float AX[nSyncRate];
    __shared__ float AY[nSyncRate];
    __shared__ float AZ[nSyncRate];
    __shared__ float Soft2[nSyncRate];
    __shared__ float ImagA[nSyncRate];

    __shared__ float reduce[nWarps][32];

    __shared__ float wX[nSyncRate][nWarps];
    __shared__ float wY[nSyncRate][nWarps];
    __shared__ float wZ[nSyncRate][nWarps];

    __shared__ float wPot[nSyncRate][nWarps];
    __shared__ float wDirsum[nSyncRate][nWarps];
    __shared__ float wNormsum[nSyncRate][nWarps];


    // Load the interaction. It is blocked for performance.
    // nWarps and ILP_PART_PER_BLK are constant, so only one code path here.
    int ii = threadIdx.x + blockIdx.x * gridDim.x; // Our interaction
    int iix, iiy;
    if (nWarps*32 == ILP_PART_PER_BLK) {
        iiy = blockIdx.x;
        iix = threadIdx.x;
        }
    else {
        iiy = ii / ILP_PART_PER_BLK;
        iix = ii % ILP_PART_PER_BLK;
        }
    float iX = blk[iiy].dx.f[iix];
    float iY = blk[iiy].dy.f[iix];
    float iZ = blk[iiy].dz.f[iix];
    float iM = blk[iiy].m.f[iix];
    float ifourh2 = blk[iiy].fourh2.f[iix];

    for(iSync=0; iSync<nPart; iSync += nSyncRate) {
        int iEnd = nPart - iSync;
        if (iEnd > nSyncRate) iEnd=nSyncRate;

        // Preload the bucket of particles
        if (threadIdx.x<iEnd) {
            i = threadIdx.x;
            X[i] = pX[iSync + i];
            Y[i] = pY[iSync + i];
            Z[i] = pZ[iSync + i];
            AX[i] = pAX[iSync + i];
            AY[i] = pAY[iSync + i];
            AZ[i] = pAZ[iSync + i];
            Soft2[i] = pSoft2[iSync + i];
            ImagA[i] = AX[i]*AX[i] + AY[i]*AY[i] + AZ[i]*AZ[i];
            if (ImagA[i] > 0.0f) ImagA[i] = rsqrtf(ImagA[i]);
            }
        __syncthreads();

        for( i=0; i<iEnd; ++i) {
            float ax=0.0f, ay=0.0f, az=0.0f, fPot=0.0f, dirsum=0.0f, normsum=0.0f;
            float dx, dy, dz, d2, fourh2,dir,dir2,dir3;
            if (ii<nInteract) {
                dx = iX + X[i];
                dy = iY + Y[i];
                dz = iZ + Z[i];
                d2 = dx*dx + dy*dy + dz*dz;
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
                    adotai = AX[i]*ax + AY[i]*ay + AZ[i]*az;
                    if (adotai > 0.0f && d2 >= Soft2[i] ) {
                        adotai *= ImagA[i];
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
        // __shared__ float wX[nSyncRate][nWarps];

        int nOut = iEnd * nWarps; // Normally 32
        if (threadIdx.x<nOut) {
            int iP    = threadIdx.x & ~(nWarps-1); // 0,4,8,...
            int iWarp = threadIdx.x &  (nWarps-1); // 0 .. 3
            int iOut  = threadIdx.x / nWarps;
            warpReduceAndStore<float,nWarps>(&wX[0][0]+iP,iWarp,&out[iOutBlock].ax[iOut]);
            warpReduceAndStore<float,nWarps>(&wY[0][0]+iP,iWarp,&out[iOutBlock].ay[iOut]);
            warpReduceAndStore<float,nWarps>(&wZ[0][0]+iP,iWarp,&out[iOutBlock].az[iOut]);
            warpReduceAndStore<float,nWarps>(&wPot[0][0]+iP,iWarp,&out[iOutBlock].fPot[iOut]);
            warpReduceAndStore<float,nWarps>(&wDirsum[0][0]+iP,iWarp,&out[iOutBlock].dirsum[iOut]);
            warpReduceAndStore<float,nWarps>(&wNormsum[0][0]+iP,iWarp,&out[iOutBlock].normsum[iOut]);
            }
	++out;
        }
    }

extern "C"
void pkdParticleWorkDone(workParticle *wp);

extern "C"
int CUDAcheckWorkPP( void *vpp, void *vwork ) {
//    CUDACTX cuda = reinterpret_cast<CUDACTX>(vpp);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    char *pHostBuf = reinterpret_cast<char *>(work->pHostBuf);
    int ib, ig, ip;

    ppResults<SYNC_RATE> *pR       = reinterpret_cast<ppResults<SYNC_RATE> *>(pHostBuf);


    for( ib=0; ib<work->ppnBuffered; ++ib) {
        workParticle *wp = work->ppWP[ib];
	PINFOOUT *pInfoOut = wp->pInfoOut;
	int nGridBlocks = (work->ppNI[ib] + 32*WARPS - 1) / (32*WARPS);
	int nPartBlocks = (wp->nP + SYNC_RATE - 1) / SYNC_RATE;

    static int bDo  =1;
    if (bDo) {
	printf("got %g %g %g\n", pR[0].ax[0], pR[0].ay[0], pR[0].az[0] );
	bDo =0;
	}



	for(ig=0; ig<nGridBlocks; ++ig) {
	    for(ip=0; ip<wp->nP; ++ip) {
		

		pInfoOut[ip].a[0]    += pR[ip/SYNC_RATE].ax[ip%SYNC_RATE];
		pInfoOut[ip].a[1]    += pR[ip/SYNC_RATE].ay[ip%SYNC_RATE];
		pInfoOut[ip].a[2]    += pR[ip/SYNC_RATE].az[ip%SYNC_RATE];
		pInfoOut[ip].fPot    += pR[ip/SYNC_RATE].fPot[ip%SYNC_RATE];
		pInfoOut[ip].dirsum  += pR[ip/SYNC_RATE].dirsum[ip%SYNC_RATE];
		pInfoOut[ip].normsum += pR[ip/SYNC_RATE].normsum[ip%SYNC_RATE];


		}
	    pR += nPartBlocks;
	    }
        pkdParticleWorkDone(wp);

        }



#if 0
    PINFOOUT *pInfoOut = pp->pInfoOut;
    float *pHostBuf = reinterpret_cast<float *>(work->pHostBuf);

    int i,j;


    const int nP = pp->work->nP;
    const int nInteract = pp->nBlocks*ILP_PART_PER_BLK + pp->nInLast;
    const int nGridBlocks = (nInteract + 32*4 - 1) / (32*4);

    printf("done work CUDA PP, %d particles, %d interactions a[0]=%g\n", nP, nInteract, pHostBuf[0]);



    for( j=0; j<nP; j++ ) {
        for( i=0; i<nGridBlocks; i++) {
            pInfoOut[j].a[0]    += pHostBuf[0];
            pInfoOut[j].a[1]    += pHostBuf[1*nP*nGridBlocks];
            pInfoOut[j].a[2]    += pHostBuf[2*nP*nGridBlocks];
            pInfoOut[j].fPot    += pHostBuf[3*nP*nGridBlocks];
            pInfoOut[j].dirsum  += pHostBuf[4*nP*nGridBlocks];
            pInfoOut[j].normsum += pHostBuf[5*nP*nGridBlocks];
            ++pHostBuf;
            }
        }
#endif
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
    work->ppSize = 0;
    work->ppnBuffered = 0;
    return work;
    }


extern "C"
void CUDA_sendWorkPP(void *cudaCtx) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDAwqNode *work = cuda->nodePP;
    if (work != NULL) {
        char *pHostBuf = reinterpret_cast<char *>(work->pHostBuf);
        char *pCudaBufIn = reinterpret_cast<char *>(work->pCudaBufIn);
        ppResults<SYNC_RATE> *pCudaBufOut = reinterpret_cast<ppResults<SYNC_RATE> *>(work->pCudaBufOut);
        int nBufferIn = cuda->inCudaBufSize - work->ppSize;
        int nBufferOut = 0;
        int i;

        CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, pHostBuf, nBufferIn, cudaMemcpyHostToDevice, work->stream));
        // Launch a kernel for each buffered work package. Could be separate streams, but problematic.
        // Why problematic? When do I do the return transfer without chaining events?
        for( i=0; i<work->ppnBuffered; ++i) {
            int nP = work->ppWP[i]->nP;
            int nPaligned = (nP+31) & ~31;
            int nInteract = work->ppNI[i];
            int nBlocks = (nInteract+ILP_PART_PER_BLK-1) / ILP_PART_PER_BLK;
            float *pX       = reinterpret_cast<float *>(pCudaBufIn);
            float *pY       = pX  + nPaligned;
            float *pZ       = pY  + nPaligned;
            float *pAX      = pZ  + nPaligned;
            float *pAY      = pAX + nPaligned;
            float *pAZ      = pAY + nPaligned;
            float *pSmooth2 = pAZ + nPaligned;
            ILP_BLK * __restrict__ blk = reinterpret_cast<ILP_BLK*>(pSmooth2 + nPaligned);
            pCudaBufIn = reinterpret_cast<char *>(blk + nBlocks );
            int nGridBlocks = (nInteract + 32*WARPS - 1) / (32*WARPS);
            dim3 dimBlock( 32*WARPS, 1 );
            dim3 dimGrid( nGridBlocks, 1,1);
            cudaPP<WARPS,MAX_BUCKET,SYNC_RATE><<<dimGrid, dimBlock, 0, work->stream>>>(
                nP, pX, pY, pZ, pAX, pAY, pAZ, pSmooth2,
                nInteract, blk,pCudaBufOut );
	    int nOutBlocks = nGridBlocks * (nP+SYNC_RATE-1)/SYNC_RATE;
	    pCudaBufOut += nOutBlocks;
            nBufferOut += sizeof(ppResults<SYNC_RATE>) * nOutBlocks;
            }
        assert(nBufferIn == reinterpret_cast<char *>(pCudaBufIn) - reinterpret_cast<char *>(work->pCudaBufIn));
        assert(nBufferOut <= cuda->inCudaBufSize);
        CUDA_CHECK(cudaMemcpyAsync,(pHostBuf, work->pCudaBufOut, nBufferOut, cudaMemcpyDeviceToHost, work->stream));
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
    PINFOIN *pInfoIn = wp->pInfoIn;
    int nP = wp->nP;
    int nPaligned = (nP+31) & ~31;
    int nBlocks = tile->lstTile.nBlocks + (tile->lstTile.nInLast?1:0);
    int nInteract = tile->lstTile.nBlocks*ILP_PART_PER_BLK + tile->lstTile.nInLast;
    int nBytesIn = nPaligned * sizeof(float) * 7 + nBlocks*sizeof(ILP_BLK);
    int i;


    // Figure out the total amount of space we need, and see if there is enough
    // in the current block. If not, send it to the GPU and get another.
    // Space: nPaligned * 4 (float) * 7 (coordinates)
    //        nBlocks * sizeof(ILP_BLK)
    //
    if (work!=NULL && nBytesIn > work->ppSize) {
        CUDA_sendWorkPP(cuda);
        work = NULL;
        assert(cuda->nodePP==NULL);
        }

    // If we don't have a PP work element, try to grab one
    if (work==NULL) {
        cuda->nodePP = work = getNode(cuda);
        if (work==NULL) return 0;
        work->ppSize = cuda->inCudaBufSize;
        work->ppnBuffered = 0;
        }
    if (work==NULL) return 0;

    char *pHostBuf = reinterpret_cast<char *>(work->pHostBuf) + cuda->inCudaBufSize - work->ppSize;

    float *pX       = reinterpret_cast<float *>(pHostBuf);
    float *pY       = pX  + nPaligned;
    float *pZ       = pY  + nPaligned;
    float *pAX      = pZ  + nPaligned;
    float *pAY      = pAX + nPaligned;
    float *pAZ      = pAY + nPaligned;
    float *pSmooth2 = pAZ + nPaligned;
    ILP_BLK * __restrict__ blk = reinterpret_cast<ILP_BLK*>(pSmooth2 + nPaligned);
    for(i=0; i<nP; ++i) {
        pX[i] =  pInfoIn[i].r[0];
        pY[i] =  pInfoIn[i].r[1];
        pZ[i] =  pInfoIn[i].r[2];
        pAX[i] =  pInfoIn[i].a[0];
        pAY[i] =  pInfoIn[i].a[1];
        pAZ[i] =  pInfoIn[i].a[2];
        pSmooth2[i] = pInfoIn[i].fSmooth2;
        }
    for(i=0; i<nBlocks; ++i) blk[i] = tile->blk[i];

//    printf("CUDA: buffering %d bytes\n", nBytesIn);
    work->ppSize -= nBytesIn;
    work->ppNI[work->ppnBuffered] = nInteract;
    work->ppWP[work->ppnBuffered] = wp;
    ++wp->nRefs;

    if ( ++work->ppnBuffered == CUDA_PP_MAX_BUFFERED) CUDA_sendWorkPP(cuda);

    return 1;
    }
