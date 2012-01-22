/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "pkd.h"
#include "ilp.h"
#include "cudautil.h"

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

#define PP_THREADS 1024
#define PP_BLKS_PER_THREAD (PP_THREADS/ILP_PART_PER_BLK)

__device__ float MWS(float m1, float h12, float m2, float h22) {
    float tmp = h12*h22;
    if (m1 == 0.0) return(h22);
    if (m2 == 0.0) return(h12);
    if (tmp > 0.0) return((m1+m2)*tmp/(h22*m1+h12*m2));
    else return(0.0);
    }

__global__ void cudaPP( int nP, PINFOIN *in, int nPart, ILP_BLK *blk, PINFOOUT *out ) {
    int bid = blockIdx.x * PP_BLKS_PER_THREAD + threadIdx.y;
    int tid = threadIdx.x + threadIdx.y*ILP_PART_PER_BLK;
    int pid = tid + blockIdx.x * PP_THREADS;
    int wid = tid & 31;
    float d2, dir, dir2, fourh2, dx, dy, dz, p;
    float fSoft, fMass;
    __shared__ float ax[32];
    __shared__ float ay[32];
    __shared__ float az[32];
    __shared__ float fPot[32];

    if (tid<32) ax[tid] = ay[tid] = az[tid] = fPot[tid] = 0.0;
    __syncthreads();

    if ( pid < nPart ) {
        blk += bid;

        dx = blk->dx.f[threadIdx.x] + in[blockIdx.y].r[0];
        dy = blk->dy.f[threadIdx.x] + in[blockIdx.y].r[1];
        dz = blk->dz.f[threadIdx.x] + in[blockIdx.y].r[2];
        d2 = dx*dx + dy*dy + dz*dz;
        fSoft = in[blockIdx.y].fSoft;
        fMass = in[blockIdx.y].fMass;
        fourh2 = MWS(fMass,4.0f*fSoft*fSoft,blk->m.f[threadIdx.x],blk->fourh2.f[threadIdx.x]);
        if (d2 == 0.0 ) dir2 = 0.0; /* Ignore self interactions */
        else if (d2 <= fourh2) {
            /*
            ** This uses the Dehnen K1 kernel function now, it's fast!
            */
            dir = rsqrtf(fourh2);
            dir2 = dir*dir;
            d2 *= dir2;
            dir2 *= dir;
            d2 = 1.0f - d2;
            dir *= 1.0f + d2*(0.5f + d2*(3.0f/8.0f + d2*(45.0f/32.0f)));
            dir2 *= 1.0f + d2*(1.5f + d2*(135.0f/16.0f));
            }
        else {
            dir = rsqrtf(d2);
            dir2 = dir*dir*dir;
            }

        dir2 *= -blk->m.f[threadIdx.x];
        dx *= dir2;
        dy *= dir2;
        dz *= dir2;
        atomicAdd(&ax[wid],dx);
        atomicAdd(&ay[wid],dy);
        atomicAdd(&az[wid],dz);
        atomicAdd(&fPot[wid],-blk->m.f[threadIdx.x]*dir);

        /*
        ** Calculations for determining the timestep.
        */
        float adotai = in[blockIdx.y].a[0]*dx + in[blockIdx.y].a[1]*dy + in[blockIdx.y].a[2]*dz;
        if (adotai > 0.0 && d2 >= in[blockIdx.y].fSmooth2) {
//	    adotai *= dimaga;
//	    *dirsum += dir*adotai*adotai;
//	    *normsum += adotai*adotai;
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
            dx = vax[tid];  dy = vay[tid]; dz = vaz[tid]; p = fPot[tid];
            { R(vax,dx,16); R(vay,dy,16); R(vaz,dz,16); R(fPot,p,16); EMUSYNC; }
            { R(vax,dx, 8); R(vay,dy, 8); R(vaz,dz, 8); R(fPot,p, 8); EMUSYNC; }
            { R(vax,dx, 4); R(vay,dy, 4); R(vaz,dz, 4); R(fPot,p, 4); EMUSYNC; }
            { R(vax,dx, 2); R(vay,dy, 2); R(vaz,dz, 2); R(fPot,p, 2); EMUSYNC; }
            { R(vax,dx, 1); R(vay,dy, 1); R(vaz,dz, 1); R(fPot,p, 1); EMUSYNC; }
            }
        if (tid==0) {
            out[blockIdx.y+nP*blockIdx.x].a[0] = ax[0];
            out[blockIdx.y+nP*blockIdx.x].a[1] = ay[0];
            out[blockIdx.y+nP*blockIdx.x].a[2] = az[0];
            out[blockIdx.y+nP*blockIdx.x].fPot = fPot[0];
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
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pp->work->pkd->cudaCtx);
    gpuInput *in;
    gpuBlock *blk;

    if (ctx->in==NULL || ctx->block==NULL) return 0; /* good luck */

    //cudaFuncSetCacheConfig(cudaPP,cudaFuncCachePreferL1);


    int j;
    PINFOOUT *pInfoOut = pp->pInfoOut;
    for( j=0; j<nP; j++ ) {
        pInfoOut[j].a[0] = 0;
        pInfoOut[j].a[1] = 0;
        pInfoOut[j].a[2] = 0;
        pInfoOut[j].fPot = 0;
#if 0
        pp->i = j;
        extern int CPUdoWorkPP(void *vpp);
        CPUdoWorkPP(pp);
        pInfoOut[j].a[0] = -pInfoOut[j].a[0];
        pInfoOut[j].a[1] = -pInfoOut[j].a[1];
        pInfoOut[j].a[2] = -pInfoOut[j].a[2];
        pInfoOut[j].fPot = -pInfoOut[j].fPot;
#endif
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

    // cuda global arrays
    int nBlocks = pp->nBlocks + (pp->nInLast ? 1 : 0);
    dim3 threads(ILP_PART_PER_BLK,PP_BLKS_PER_THREAD);
    dim3 numBlocks((nBlocks+PP_BLKS_PER_THREAD-1)/PP_BLKS_PER_THREAD, nP);
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
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pp->work->pkd->cudaCtx);
    PINFOOUT *pInfoOut = pp->pInfoOut;
    int nP = pp->work->nP;
    int numBlocks = (pp->nBlocks + (pp->nInLast ? 1 : 0) + PP_BLKS_PER_THREAD - 1)/PP_BLKS_PER_THREAD;
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
        for( i=0; i<numBlocks; i++) {
            pInfoOut[j].a[0] += blk->cpuResults[j+nP*i].a[0];
            pInfoOut[j].a[1] += blk->cpuResults[j+nP*i].a[1];
            pInfoOut[j].a[2] += blk->cpuResults[j+nP*i].a[2];
            pInfoOut[j].fPot += blk->cpuResults[j+nP*i].fPot;
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
