/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "pkd.h"
#include "ilc.h"
#include "cudautil.h"

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

/* With 512 threads we can have better occupancy as 1536 can be active. */
#define PC_THREADS 128
#define PC_BLKS_PER_THREAD (PC_THREADS/ILC_PART_PER_BLK)

__global__ void cudaPC( int nP, PINFOIN *in, int nPart, ILC_BLK *blk, PINFOOUT *out ) {
    int bid = blockIdx.x * PC_BLKS_PER_THREAD + threadIdx.y;
    int tid = threadIdx.x + threadIdx.y*ILC_PART_PER_BLK;
    int pid = tid + blockIdx.x * PC_THREADS;
    int wid = tid & 31;
    const float onethird = 1.0f/3.0f;
    float d2, dir, dx, dy, dz, p;
    float u,g0,g2,g3,g4;
    float tax, tay, taz;
    float x,y,z;
//    float adotai;
    float tx,ty,tz;
    float xx,xy,xz,yy,yz,zz;
    float xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    int i = blockIdx.y;

    __shared__ float ax[32];
    __shared__ float ay[32];
    __shared__ float az[32];
    __shared__ float fPot[32];

    blk += bid;

//    for( i=0; i<nP; ++i ) {

    if (tid<32) ax[tid] = ay[tid] = az[tid] = fPot[tid] = 0.0f;
    __syncthreads();

    if ( pid < nPart ) {

        dx = blk->dx.f[threadIdx.x] + in[i].r[0];
        dy = blk->dy.f[threadIdx.x] + in[i].r[1];
        dz = blk->dz.f[threadIdx.x] + in[i].r[2];
        d2 = dx*dx + dy*dy + dz*dz;
        dir = rsqrtf(d2);
	    u = blk->u.f[threadIdx.x]*dir;
	    g2 = 3.0f*dir*u*u;
	    g3 = 5.0f*g2*u;
	    g4 = 7.0f*g3*u;
	    /*
	    ** Calculate the funky distance terms.
	    */
	    x = dx*dir;
	    y = dy*dir;
	    z = dz*dir;
	    xx = 0.5f*x*x;
	    xy = x*y;
	    xz = x*z;
	    yy = 0.5f*y*y;
	    yz = y*z;
	    zz = 0.5f*z*z;
	    xxx = x*(onethird*xx - zz);
	    xxz = z*(xx - onethird*zz);
	    yyy = y*(onethird*yy - zz);
	    yyz = z*(yy - onethird*zz);
	    xx -= zz;
	    yy -= zz;
	    xxy = y*xx;
	    xyy = x*yy;
	    xyz = xy*z;
	    /*
	    ** Now calculate the interaction up to Hexadecapole order.
	    */
	    tx = g4*(blk->xxxx.f[threadIdx.x]*xxx + blk->xyyy.f[threadIdx.x]*yyy + blk->xxxy.f[threadIdx.x]*xxy + blk->xxxz.f[threadIdx.x]*xxz + blk->xxyy.f[threadIdx.x]*xyy + blk->xxyz.f[threadIdx.x]*xyz + blk->xyyz.f[threadIdx.x]*yyz);
	    ty = g4*(blk->xyyy.f[threadIdx.x]*xyy + blk->xxxy.f[threadIdx.x]*xxx + blk->yyyy.f[threadIdx.x]*yyy + blk->yyyz.f[threadIdx.x]*yyz + blk->xxyy.f[threadIdx.x]*xxy + blk->xxyz.f[threadIdx.x]*xxz + blk->xyyz.f[threadIdx.x]*xyz);
	    tz = g4*(-blk->xxxx.f[threadIdx.x]*xxz - (blk->xyyy.f[threadIdx.x] + blk->xxxy.f[threadIdx.x])*xyz - blk->yyyy.f[threadIdx.x]*yyz + blk->xxxz.f[threadIdx.x]*xxx + blk->yyyz.f[threadIdx.x]*yyy - blk->xxyy.f[threadIdx.x]*(xxz + yyz) + blk->xxyz.f[threadIdx.x]*xxy + blk->xyyz.f[threadIdx.x]*xyy);
	    g4 = 0.25f*(tx*x + ty*y + tz*z);
	    xxx = g3*(blk->xxx.f[threadIdx.x]*xx + blk->xyy.f[threadIdx.x]*yy + blk->xxy.f[threadIdx.x]*xy + blk->xxz.f[threadIdx.x]*xz + blk->xyz.f[threadIdx.x]*yz);
	    xxy = g3*(blk->xyy.f[threadIdx.x]*xy + blk->xxy.f[threadIdx.x]*xx + blk->yyy.f[threadIdx.x]*yy + blk->yyz.f[threadIdx.x]*yz + blk->xyz.f[threadIdx.x]*xz);
	    xxz = g3*(-(blk->xxx.f[threadIdx.x] + blk->xyy.f[threadIdx.x])*xz - (blk->xxy.f[threadIdx.x] + blk->yyy.f[threadIdx.x])*yz + blk->xxz.f[threadIdx.x]*xx + blk->yyz.f[threadIdx.x]*yy + blk->xyz.f[threadIdx.x]*xy);
	    g3 = onethird*(xxx*x + xxy*y + xxz*z);
	    xx = g2*(blk->xx.f[threadIdx.x]*x + blk->xy.f[threadIdx.x]*y + blk->xz.f[threadIdx.x]*z);
	    xy = g2*(blk->yy.f[threadIdx.x]*y + blk->xy.f[threadIdx.x]*x + blk->yz.f[threadIdx.x]*z);
	    xz = g2*(-(blk->xx.f[threadIdx.x] + blk->yy.f[threadIdx.x])*z + blk->xz.f[threadIdx.x]*x + blk->yz.f[threadIdx.x]*y);
	    g2 = 0.5f*(xx*x + xy*y + xz*z);
	    g0 = dir * blk->m.f[threadIdx.x];
        atomicAdd(&fPot[wid],-(g0 + g2 + g3 + g4));
	    g0 += 5.0f*g2 + 7.0f*g3 + 9.0f*g4;
	    tax = dir*(xx + xxx + tx - x*g0);
	    tay = dir*(xy + xxy + ty - y*g0);
	    taz = dir*(xz + xxz + tz - z*g0);
	    /*
	    ** Calculations for determining the timestep.
	    */
//	    adotai = pInfoIn[i].a[0]*tax + pInfoIn[i].a[1]*tay + pInfoIn[i].a[2]*taz;
//	    if (adotai > 0) {
//		adotai *= dimaga;
//		dirsum += dir*adotai*adotai;
//		normsum += adotai*adotai;

        atomicAdd(&ax[wid],tax);
        atomicAdd(&ay[wid],tay);
        atomicAdd(&az[wid],taz);

        /*
        ** Calculations for determining the timestep.
        */
//PP        float adotai = in[blockIdx.y].a[0]*dx + in[blockIdx.y].a[1]*dy + in[blockIdx.y].a[2]*dz;
//        if (adotai > 0.0f && d2 >= in[blockIdx.y].fSmooth2) {
//	    adotai *= dimaga;
//	    *dirsum += dir*adotai*adotai;
//	    *normsum += adotai*adotai;
//	    }

        /*
        ** Now reduce the result
        */
        }

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
        out[i+nP*blockIdx.x].a[0] = ax[0];
        out[i+nP*blockIdx.x].a[1] = ay[0];
        out[i+nP*blockIdx.x].a[2] = az[0];
        out[i+nP*blockIdx.x].fPot = fPot[0];
        out[i+nP*blockIdx.x].dirsum = 0;
        out[i+nP*blockIdx.x].normsum = 0;
        }
//    }

    }

/*
** Queue the given work element for execution on the GPU and return 1.
** If we cannot, then we return 0.
*/
extern "C"
int CUDAinitWorkPC( void *vpp) {
    workPC *pc = reinterpret_cast<workPC *>(vpp);
    int nP = pc->work->nP;
    ILCTILE tile = pc->tile;
    int nPart = pc->nBlocks*ILC_PART_PER_BLK + pc->nInLast;
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pc->work->pkd->cudaCtx);
    gpuInput *in;
    gpuBlock *blk;

    if (ctx->in==NULL || ctx->block==NULL) return 0; /* good luck */

    //cudaFuncSetCacheConfig(cudaPC,cudaFuncCachePreferL1);

    int j;
    PINFOOUT *pInfoOut = pc->pInfoOut;
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
    pc->gpu_memory = blk;

    // The particles need only be copied once. We can use the stream from the first block.
    if ((in=reinterpret_cast<gpuInput *>(pc->work->gpu_memory)) == NULL) {
        in = ctx->in;
        ctx->in = in->next;
        pc->work->gpu_memory = in;
        in->nRefs = 0;
        memcpy(in->cpuIn,pc->work->pInfoIn,sizeof(PINFOIN)*nP);
        CUDA_CHECK(cudaMemcpyAsync,(in->in, in->cpuIn, sizeof(PINFOIN)*nP, cudaMemcpyHostToDevice, blk->stream));;
        CUDA_CHECK(cudaEventRecord,(in->event,blk->stream));
        }
    in->nRefs++;

    // cuda global arrays

    int nBlocks = pc->nBlocks + (pc->nInLast ? 1 : 0);
    dim3 threads(ILC_PART_PER_BLK,PC_BLKS_PER_THREAD);
    dim3 numBlocks((nBlocks+PC_BLKS_PER_THREAD-1)/PC_BLKS_PER_THREAD, nP);
    int bytes = sizeof(ILC_BLK) * (nBlocks);

    // copy data directly to device memory
    CUDA_CHECK(cudaMemcpyAsync,(blk->gpuBlkILC, tile->blk, bytes, cudaMemcpyHostToDevice, blk->stream));
    // We need the particles, so we need to wait for that transfer to complete
    CUDA_CHECK(cudaStreamWaitEvent,(blk->stream, in->event, 0));
    cudaPC<<<numBlocks, threads, 0, blk->stream>>>( nP, in->in, nPart, blk->gpuBlkILC, blk->gpuResults);
    CUDA_CHECK(cudaMemcpyAsync,(blk->cpuResults, blk->gpuResults, nP*numBlocks.x*sizeof(PINFOOUT),
            cudaMemcpyDeviceToHost, blk->stream));
    CUDA_CHECK(cudaEventRecord,(blk->event,blk->stream));

    return 1;
    }

extern "C"
int CUDAcheckWorkPC( void *vpc ) {
    workPC *pc = reinterpret_cast<workPC *>(vpc);
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pc->work->pkd->cudaCtx);
    PINFOOUT *pInfoOut = pc->pInfoOut;
    int nP = pc->work->nP;
    int numBlocks = (pc->nBlocks + (pc->nInLast ? 1 : 0) + PC_BLKS_PER_THREAD - 1)/PC_BLKS_PER_THREAD;
    cudaError_t rc;
    gpuBlock *blk;
    gpuInput *in;
    int i,j;

    blk = reinterpret_cast<gpuBlock *>(pc->gpu_memory);
    in = reinterpret_cast<gpuInput *>(pc->work->gpu_memory);
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
            pInfoOut[j].dirsum += blk->cpuResults[j+nP*i].dirsum;
            pInfoOut[j].normsum+= blk->cpuResults[j+nP*i].normsum;
            }
        }

    // Return the memory
    if ( --in->nRefs == 0 ) {
        in->next = ctx->in;
        ctx->in = in;
        pc->work->gpu_memory = NULL;
        }

    blk->next = ctx->block;
    ctx->block = blk;
    pc->gpu_memory = NULL;

    return 0;
    }

extern "C"
void CUDAsetupPC(void) {
    cudaFuncSetCacheConfig(cudaPC,cudaFuncCachePreferL1);
    }
