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

__global__ void cudaPC( int nP, PINFOIN *in, int nPart, ILC_BLK *blk, PINFOOUT *out ) {
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
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


    __shared__ float ax[32];
    __shared__ float ay[32];
    __shared__ float az[32];
    __shared__ float fPot[32];

    if (tid<32) ax[tid] = ay[tid] = az[tid] = fPot[tid] = 0.0;
    __syncthreads();

    if ( pid < nPart ) {
        blk += blockIdx.x;

        dx = blk->dx.f[tid] + in[blockIdx.y].r[0];
        dy = blk->dy.f[tid] + in[blockIdx.y].r[1];
        dz = blk->dz.f[tid] + in[blockIdx.y].r[2];
        d2 = dx*dx + dy*dy + dz*dz;
        dir = rsqrtf(d2);
	    u = blk->u.f[tid]*dir;
	    g0 = dir;
	    g2 = 3*dir*u*u;
	    g3 = 5*g2*u;
	    g4 = 7*g3*u;
	    /*
	    ** Calculate the funky distance terms.
	    */
	    x = dx*dir;
	    y = dy*dir;
	    z = dz*dir;
	    xx = 0.5*x*x;
	    xy = x*y;
	    xz = x*z;
	    yy = 0.5*y*y;
	    yz = y*z;
	    zz = 0.5*z*z;
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
	    tx = g4*(blk->xxxx.f[tid]*xxx + blk->xyyy.f[tid]*yyy + blk->xxxy.f[tid]*xxy + blk->xxxz.f[tid]*xxz + blk->xxyy.f[tid]*xyy + blk->xxyz.f[tid]*xyz + blk->xyyz.f[tid]*yyz);
	    ty = g4*(blk->xyyy.f[tid]*xyy + blk->xxxy.f[tid]*xxx + blk->yyyy.f[tid]*yyy + blk->yyyz.f[tid]*yyz + blk->xxyy.f[tid]*xxy + blk->xxyz.f[tid]*xxz + blk->xyyz.f[tid]*xyz);
	    tz = g4*(-blk->xxxx.f[tid]*xxz - (blk->xyyy.f[tid] + blk->xxxy.f[tid])*xyz - blk->yyyy.f[tid]*yyz + blk->xxxz.f[tid]*xxx + blk->yyyz.f[tid]*yyy - blk->xxyy.f[tid]*(xxz + yyz) + blk->xxyz.f[tid]*xxy + blk->xyyz.f[tid]*xyy);
	    g4 = 0.25*(tx*x + ty*y + tz*z);
	    xxx = g3*(blk->xxx.f[tid]*xx + blk->xyy.f[tid]*yy + blk->xxy.f[tid]*xy + blk->xxz.f[tid]*xz + blk->xyz.f[tid]*yz);
	    xxy = g3*(blk->xyy.f[tid]*xy + blk->xxy.f[tid]*xx + blk->yyy.f[tid]*yy + blk->yyz.f[tid]*yz + blk->xyz.f[tid]*xz);
	    xxz = g3*(-(blk->xxx.f[tid] + blk->xyy.f[tid])*xz - (blk->xxy.f[tid] + blk->yyy.f[tid])*yz + blk->xxz.f[tid]*xx + blk->yyz.f[tid]*yy + blk->xyz.f[tid]*xy);
	    g3 = onethird*(xxx*x + xxy*y + xxz*z);
	    xx = g2*(blk->xx.f[tid]*x + blk->xy.f[tid]*y + blk->xz.f[tid]*z);
	    xy = g2*(blk->yy.f[tid]*y + blk->xy.f[tid]*x + blk->yz.f[tid]*z);
	    xz = g2*(-(blk->xx.f[tid] + blk->yy.f[tid])*z + blk->xz.f[tid]*x + blk->yz.f[tid]*y);
	    g2 = 0.5*(xx*x + xy*y + xz*z);
	    g0 *= blk->m.f[tid];
        atomicAdd(&fPot[wid],-(g0 + g2 + g3 + g4));
	    g0 += 5*g2 + 7*g3 + 9*g4;
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
//        if (adotai > 0.0 && d2 >= in[blockIdx.y].fSmooth2) {
//	    adotai *= dimaga;
//	    *dirsum += dir*adotai*adotai;
//	    *normsum += adotai*adotai;
//	    }

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
int CUDAinitWorkPC( void *vpp) {
    workPC *pc = reinterpret_cast<workPC *>(vpp);
    int nP = pc->work->nP;
    ILCTILE tile = pc->tile;
    int nPart = pc->nBlocks*ILC_PART_PER_BLK + pc->nInLast;
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pc->work->pkd->cudaCtx);
    gpuInput *in;
    gpuBlock *blk;

    if (ctx->in==NULL || ctx->block==NULL) return 0; /* good luck */

    int j;
    PINFOOUT *pInfoOut = pc->pInfoOut;
    for( j=0; j<nP; j++ ) {
        pInfoOut[j].a[0] = 0;
        pInfoOut[j].a[1] = 0;
        pInfoOut[j].a[2] = 0;
        pInfoOut[j].fPot = 0;
#if 0
        pc->i = j;
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
    int threads = ILC_PART_PER_BLK;
    dim3 numBlocks(pc->nBlocks + (pc->nInLast ? 1 : 0), nP);
    int bytes = sizeof(ILC_BLK) * (numBlocks.x);

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
int CUDAcheckWorkPC( void *vpp ) {
    workPP *pc = reinterpret_cast<workPP *>(vpp);
    CUDACTX ctx = reinterpret_cast<CUDACTX>(pc->work->pkd->cudaCtx);
    PINFOOUT *pInfoOut = pc->pInfoOut;
    int nP = pc->work->nP;
    int numBlocks = pc->nBlocks + (pc->nInLast ? 1 : 0);
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
