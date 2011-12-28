/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "pkd.h"
#include "ilp.h"

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

typedef struct {
    cudaStream_t stream;
    PINFOIN  *g_in;
    ILP_BLK  *g_blk;
    PINFOOUT *g_out;
    PINFOOUT *aout;
    PINFOIN *pInfo;
    } PKDCUDA;

extern "C"
bool isPow2(unsigned int x)
{
    return ((x&(x-1))==0);
}

__device__ float MWS(float m1, float h12, float m2, float h22) {
    float tmp = h12*h22;
    if (m1 == 0.0) return(h22);
    if (m2 == 0.0) return(h12);
    if (tmp > 0.0) return((m1+m2)*tmp/(h22*m1+h12*m2));
    else return(0.0);
    }

//static __constant__ PINFOIN in[64];

__global__ void cudaPP( int nP, PINFOIN *in, int nPart, ILP_BLK *blk, PINFOOUT *out ) {
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;
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
        blk += blockIdx.x;

        dx = blk->dx.f[tid] + in[blockIdx.y].r[0];
        dy = blk->dy.f[tid] + in[blockIdx.y].r[1];
        dz = blk->dz.f[tid] + in[blockIdx.y].r[2];
        d2 = dx*dx + dy*dy + dz*dz;
        fSoft = in[blockIdx.y].fSoft;
        fMass = in[blockIdx.y].fMass;
        fourh2 = MWS(fMass,4.0f*fSoft*fSoft,blk->m.f[tid],blk->fourh2.f[tid]);
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

        dir2 *= -blk->m.f[tid];
        dx *= dir2;
        dy *= dir2;
        dz *= dir2;
        atomicAdd(&ax[wid],dx);
        atomicAdd(&ay[wid],dy);
        atomicAdd(&az[wid],dz);
        atomicAdd(&fPot[wid],-blk->m.f[tid]*dir);

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

#if 1
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
#endif
        if (tid==0) {
            out[blockIdx.y+nP*blockIdx.x].a[0] = ax[0];
            out[blockIdx.y+nP*blockIdx.x].a[1] = ay[0];
            out[blockIdx.y+nP*blockIdx.x].a[2] = az[0];
            out[blockIdx.y+nP*blockIdx.x].p    = fPot[0];
            }
        }
    }

extern "C" ILP_BLK * pkdGravCudaPPAllocateBlk() {
    ILP_BLK *blk;
    cudaHostAlloc( (void**)&blk, sizeof(ILP_BLK) * ILP_BLK_PER_TILE, 0 ? cudaHostAllocWriteCombined : 0 );
//    cudaHostAlloc( (void**)&blk, sizeof(ILP_BLK) * ILP_BLK_PER_TILE, cudaHostAllocMapped | cudaHostAllocWriteCombined);
    return blk;
    }

extern "C" void pkdGravCudaPPFreeBlk(ILP_BLK *blk) {
    cudaFreeHost( blk );
    }

extern "C" void *pkdGravCudaPPAllocate() {
    PKDCUDA *ctx = reinterpret_cast<PKDCUDA *>(malloc(sizeof(PKDCUDA)));
    size_t bytes = ILP_BLK_PER_TILE * sizeof(ILP_BLK);

    cudaStreamCreate(&ctx->stream);
    cudaMalloc((void**) &ctx->g_blk, bytes);
    cudaMalloc((void**) &ctx->g_out, ILP_BLK_PER_TILE * PKD_GROUP_SIZE * sizeof(PINFOOUT) );
    cudaMalloc((void**) &ctx->g_in, PKD_GROUP_SIZE * sizeof(PINFOIN) );

    cudaHostAlloc( (void**)&ctx->pInfo, sizeof(PINFOIN) * PKD_GROUP_SIZE, 0 ? cudaHostAllocWriteCombined : 0 );
    cudaHostAlloc( (void**)&ctx->aout, sizeof(PINFOOUT) * ILP_BLK_PER_TILE * PKD_GROUP_SIZE,0 );

    return ctx;
    }

extern "C" void pkdGravCudaPPFree(void *cudaCtx) {
    PKDCUDA *ctx = reinterpret_cast<PKDCUDA *>(cudaCtx);
    cudaFree(ctx->g_in);
    cudaFree(ctx->g_out);
    cudaFree(ctx->g_blk);
    cudaStreamSynchronize(ctx->stream);
    cudaStreamDestroy(ctx->stream);
    free(ctx);
    }

extern "C"
    void pkdGravCudaPP( PKD pkd, void *cudaCtx, int nP, PARTICLE **pPart, PINFOIN *pInfoIn, ILP ilp ) {
    PKDCUDA *ctx = reinterpret_cast<PKDCUDA *>(cudaCtx);
    int i, j;
    ILPTILE tile;

//    printf( "ilp->first=%p ilp->first->blk=%p\n", ilp->first, ilp->first->blk);


    // Copy the particles
    memcpy(ctx->pInfo,pInfoIn,sizeof(PINFOIN)*nP);
    cudaMemcpyAsync(ctx->g_in, ctx->pInfo, sizeof(PINFOIN)*nP, cudaMemcpyHostToDevice, ctx->stream);
//    cudaMemcpyAsync(ctx->g_in, pInfoIn, sizeof(PINFOIN)*nP, cudaMemcpyHostToDevice, ctx->stream);

    //cudaMemcpyToSymbol(in, pInfoIn, sizeof(PINFOIN)*nP);

	ILP_LOOP(ilp,tile) {
        // cuda global arrays
        int threads = ILP_PART_PER_BLK;
        dim3 numBlocks((tile->nPart + threads - 1) / threads, nP);
        int bytes = sizeof(ILP_BLK) * (numBlocks.x);

        // copy data directly to device memory
        cudaMemcpyAsync(ctx->g_blk, tile->blk, bytes, cudaMemcpyHostToDevice, ctx->stream);
        cudaPP<<<numBlocks, threads, 0, ctx->stream>>>( nP, ctx->g_in, tile->nPart, ctx->g_blk, ctx->g_out);

//        ILP_BLK *in;
//        cudaHostGetDevicePointer(&in, tile->blk, 0);
//        cudaPP<<<numBlocks, threads, 0, ctx->stream>>>( nP, ctx->g_in, tile->nPart, in, ctx->g_out);

        cudaMemcpyAsync(ctx->aout, ctx->g_out, nP*numBlocks.x*sizeof(PINFOOUT), cudaMemcpyDeviceToHost, ctx->stream);

        cudaStreamSynchronize(ctx->stream);
#if 1
        for( j=0; j<nP; j++ ) {
            float *a = pkdAccel(pkd,pPart[j]);
            for( i=0; i<numBlocks.x; i++) {
                a[0] += ctx->aout[j+nP*i].a[0];
                a[1] += ctx->aout[j+nP*i].a[1];
                a[2] += ctx->aout[j+nP*i].a[2];
                }
            }
#endif
	    }
    }
