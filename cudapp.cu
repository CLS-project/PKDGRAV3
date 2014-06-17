/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
//#include "pkd.h"
#include "ilp.h"
#include "cudautil.h"

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

/* With 512 threads we can have better occupancy as 1536 can be active. */
#define PP_THREADS 256
#define PP_BLKS_PER_THREAD (PP_THREADS/ILP_PART_PER_BLK)

#define NBLOCKS 1

__global__ void cudaPP( int nP, PINFOIN *in, int nPart, ILP_BLK *blk, PINFOOUT *out ) {
    int bid, tid, pid, wid;
    float d2, dir, dir2, dir3, fourh2, dx, dy, dz, p, ds, ns;
    __shared__ float ax[32];
    __shared__ float ay[32];
    __shared__ float az[32];
    __shared__ float fPot[32];
    __shared__ float dirsum[32];
    __shared__ float normsum[32];
    float m, dimaga, adotai;
    int i = blockIdx.y;
    float *a = in[i].a;

    bid = blockIdx.x * PP_BLKS_PER_THREAD + threadIdx.y;
    tid = threadIdx.x + threadIdx.y*ILP_PART_PER_BLK;
    pid = tid + blockIdx.x * PP_THREADS;
    if (tid<32) ax[tid] = ay[tid] = az[tid] = fPot[tid] = dirsum[tid] = normsum[tid] = 0.0f;
    __syncthreads();

    blk += bid;
    wid = tid & 31;

    dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0.0f) {
        dimaga = rsqrtf(dimaga);
        }
    while(pid < nPart) {
        dx = blk->dx.f[threadIdx.x] + in[i].r[0];
        dy = blk->dy.f[threadIdx.x] + in[i].r[1];
        dz = blk->dz.f[threadIdx.x] + in[i].r[2];
        d2 = dx*dx + dy*dy + dz*dz;

        m = blk->m.f[threadIdx.x];
        fourh2 = blk->fourh2.f[threadIdx.x];
        if (d2 != 0.0f ) { /* Ignore self interactions */
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

            dir3 *= -m;
            dx *= dir3;
            dy *= dir3;
            dz *= dir3;

            atomicAdd(&ax[wid],dx);
            atomicAdd(&ay[wid],dy);
            atomicAdd(&az[wid],dz);
            atomicAdd(&fPot[wid],-m*dir);

            /*
            ** Calculations for determining the timestep.
            */
            adotai = a[0]*dx + a[1]*dy + a[2]*dz;
            if (adotai > 0.0f && d2 >= in[i].fSmooth2) {
                adotai *= dimaga;
                atomicAdd(&dirsum[wid],dir*adotai*adotai);
                atomicAdd(&normsum[wid],adotai*adotai);
                }
            }
        blk += PP_BLKS_PER_THREAD * NBLOCKS;
        pid += PP_THREADS * NBLOCKS;
        }

    /*
    ** Now reduce the result
    */

#define R(S,T,O) S[tid] = T = T + S[tid + O]

    __syncthreads();

#ifndef __DEVICE_EMULATION__
    if (tid < 16)
#endif
        {
        // now that we are using warp-synchronous programming (below)
        // we need to declare our shared memory volatile so that the compiler
        // doesn't reorder stores to it and induce incorrect behavior.
        volatile float * vax = ax;
        volatile float * vay = ay;
        volatile float * vaz = az;
        volatile float * vPot = fPot;
        volatile float * vdirsum = dirsum;
        volatile float * vnormsum = normsum;
        
        dx = vax[tid];  dy = vay[tid]; dz = vaz[tid]; p = vPot[tid]; ds = vdirsum[tid]; ns = vnormsum[tid];
        { R(vax,dx,16); R(vay,dy,16); R(vaz,dz,16); R(vPot,p,16); R(vdirsum,ds,16); R(vnormsum,ns,16); EMUSYNC; }
        { R(vax,dx, 8); R(vay,dy, 8); R(vaz,dz, 8); R(vPot,p, 8); R(vdirsum,ds, 8); R(vnormsum,ns, 8); EMUSYNC; }
        { R(vax,dx, 4); R(vay,dy, 4); R(vaz,dz, 4); R(vPot,p, 4); R(vdirsum,ds, 4); R(vnormsum,ns, 4); EMUSYNC; }
        { R(vax,dx, 2); R(vay,dy, 2); R(vaz,dz, 2); R(vPot,p, 2); R(vdirsum,ds, 2); R(vnormsum,ns, 2); EMUSYNC; }
        { R(vax,dx, 1); R(vay,dy, 1); R(vaz,dz, 1); R(vPot,p, 1); R(vdirsum,ds, 1); R(vnormsum,ns, 1); EMUSYNC; }

        if (tid==0) {
            out[i+nP*blockIdx.x].a[0] = ax[0];
            out[i+nP*blockIdx.x].a[1] = ay[0];
            out[i+nP*blockIdx.x].a[2] = az[0];
            out[i+nP*blockIdx.x].fPot = fPot[0];
            out[i+nP*blockIdx.x].dirsum = dirsum[0];
            out[i+nP*blockIdx.x].normsum= normsum[0];
            }
        }
    }

/*
** Queue the given work element for execution on the GPU and return 1.
** If we cannot, then we return 0.
*/
extern "C"
int CUDAinitWorkPP( void *vpp ) {
    workPP *pp = reinterpret_cast<workPP *>(vpp);
    int nP = pp->work->nP;
    ILPTILE tile = pp->tile;
    int nPart = pp->nBlocks*ILP_PART_PER_BLK + pp->nInLast;
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pp->work->cudaCtx);
    gpuInput *in;
    gpuBlock *blk;
    int j;

    if (ctx->in==NULL || ctx->block==NULL) return 0; /* good luck */

    PINFOOUT *pInfoOut = pp->pInfoOut;
    for( j=0; j<nP; j++ ) {
        pInfoOut[j].a[0] = 0;
        pInfoOut[j].a[1] = 0;
        pInfoOut[j].a[2] = 0;
        pInfoOut[j].fPot = 0;
        pInfoOut[j].dirsum = 0;
        pInfoOut[j].normsum = 0;
        }

    // Grab a block of memory
    blk = ctx->block;
    ctx->block = blk->next;
    pp->gpu_memory = blk;

    // The particles need only be copied once. We can use the stream from the first block.
    if ((in=reinterpret_cast<gpuInput *>(pp->work->gpu_memory)) == NULL) {
        in = ctx->in;
        ctx->in = in->next;
        pp->work->gpu_memory = in;
        in->nRefs = 0;
        memcpy(in->cpuIn,pp->work->pInfoIn,sizeof(PINFOIN)*nP);
        CUDA_CHECK(cudaMemcpyAsync,(in->in, in->cpuIn, sizeof(PINFOIN)*nP, cudaMemcpyHostToDevice, blk->stream));;
        CUDA_CHECK(cudaEventRecord,(in->event,blk->stream));
        }
    in->nRefs++;

    int nBlocks = pp->nBlocks + (pp->nInLast ? 1 : 0);
    pp->nCudaBlks = (nBlocks+PP_BLKS_PER_THREAD-1)/(PP_BLKS_PER_THREAD);
    if (pp->nCudaBlks > NBLOCKS ) pp->nCudaBlks = NBLOCKS;

    // cuda global arrays
    dim3 threads(ILP_PART_PER_BLK,PP_BLKS_PER_THREAD);
    dim3 numBlocks(pp->nCudaBlks, nP );
    int bytes = sizeof(ILP_BLK) * (nBlocks);

    // copy data directly to device memory
    CUDA_CHECK(cudaMemcpyAsync,(blk->gpuBlk, tile->blk, bytes, cudaMemcpyHostToDevice, blk->stream));
    // We need the particles, so we need to wait for that transfer to complete
    CUDA_CHECK(cudaStreamWaitEvent,(blk->stream, in->event, 0));
    cudaPP<<<numBlocks, threads, 0, blk->stream>>>( nP, in->in, nPart, blk->gpuBlk, blk->gpuResults);
    CUDA_CHECK(cudaMemcpyAsync,(blk->cpuResults, blk->gpuResults, nP*numBlocks.x*sizeof(PINFOOUT),
            cudaMemcpyDeviceToHost, blk->stream));
    CUDA_CHECK(cudaEventRecord,(blk->event,blk->stream));

    return 1;
    }

extern "C"
int CUDAcheckWorkPP( void *vpp ) {
    workPP *pp = reinterpret_cast<workPP *>(vpp);
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pp->work->cudaCtx);
    PINFOOUT *pInfoOut = pp->pInfoOut;
    int nP = pp->work->nP;
    cudaError_t rc;
    gpuBlock *blk;
    gpuInput *in;
    int i,j;

    blk = reinterpret_cast<gpuBlock *>(pp->gpu_memory);
    in = reinterpret_cast<gpuInput *>(pp->work->gpu_memory);
    rc = cudaEventQuery(blk->event);
    if (rc==cudaErrorNotReady) return 1;
    else if (rc!=cudaSuccess) {
        fprintf(stderr,"cudaEventQuery error %d: %s\n", rc, cudaGetErrorString(rc));
        exit(1);
        }
    for( j=0; j<nP; j++ ) {
        for( i=0; i<pp->nCudaBlks; i++) {
            pInfoOut[j].a[0] += blk->cpuResults[j+nP*i].a[0];
            pInfoOut[j].a[1] += blk->cpuResults[j+nP*i].a[1];
            pInfoOut[j].a[2] += blk->cpuResults[j+nP*i].a[2];
            pInfoOut[j].fPot += blk->cpuResults[j+nP*i].fPot;
            pInfoOut[j].dirsum += blk->cpuResults[j+nP*i].dirsum;
            pInfoOut[j].normsum += blk->cpuResults[j+nP*i].normsum;
            }
        }

    // Return the memory
    if ( --in->nRefs == 0 ) {
        in->next = ctx->in;
        ctx->in = in;
        pp->work->gpu_memory = NULL;
        }

    blk->next = ctx->block;
    ctx->block = blk;
    pp->gpu_memory = NULL;

    return 0;
    }


extern "C"
void CUDAsetupPP(void) {
    cudaFuncSetCacheConfig(cudaPP,cudaFuncCachePreferL1);
    }
