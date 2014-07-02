/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include "basetype.h"
#include "moments.h"
#include "cudautil.h"


#define MAX_TOTAL_REPLICAS (7*7*7)
#define SCAN_SIZE 512 // Must be larger than MAX_TOTAL_REPLICAS and a power of two.

__constant__ struct EwaldVariables ew;
__constant__ float hx[MAX_TOTAL_REPLICAS];
__constant__ float hy[MAX_TOTAL_REPLICAS];
__constant__ float hz[MAX_TOTAL_REPLICAS];
__constant__ float hCfac[MAX_TOTAL_REPLICAS];
__constant__ float hSfac[MAX_TOTAL_REPLICAS];

/*
** nvcc -DHAVE_CONFIG_H --ptxas-options=-v -c  -I. -arch=sm_35 cudaewald.cu
** ptxas info    : 0 bytes gmem, 7464 bytes cmem[3]
** ptxas info    : Compiling entry function '_Z9cudaEwaldPdS_S_S_' for 'sm_35'
** ptxas info    : Function properties for _Z9cudaEwaldPdS_S_S_
**     8 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
** ptxas info    : Used 82 registers, 2176 bytes smem, 352 bytes cmem[0], 716 bytes cmem[2]
** ptxas info    : Function properties for __internal_trig_reduction_slowpathd
**     40 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
*/


/*
** threadIdx.x: all work on the same particle -- this is the warp size, i.e., 32
** blockIdx.x:  different particles. If y=z=1, then x can be anything, otherwise
**              the total number of particles is a block of x*y*z
**
** We are allowed 16 resident blocks so this corresponds to 512 threads per SM.
** This is fine because we are actually close to shared memory limited:
**   bValid[512] = 2052 bytes
**   We have 48k per SM so around 23 active thread blocks (but we use only 16)
**   compute 5.0: 32 thread blocks, but 64K shared gives us ~ 31 resident.
*/

__global__ void cudaEwald(double *X,double *Y,double *Z,double *pPot) {
    double x, y, z;
    double r2,dir,dir2,a;
    double rx, ry, rz;
    double g0,g1,g2,g3,g4,g5,alphan;
    int i, ix, iy, iz, bInHole;
    int pidx = threadIdx.x + 32*blockIdx.x;
    double tax = 0.0, tay = 0.0, taz = 0.0, tpot=0.0;

    rx = X[pidx] - ew.r[0];
    ry = Y[pidx] - ew.r[1];
    rz = Z[pidx] - ew.r[2];

    for(ix=-3; ix<=3; ++ix) {
        for(iy=-3; iy<=3; ++iy) {
            for(iz=-3; iz<=3; ++iz) {
                if (ix==0 && iy==0 & iz==0) continue;
                bInHole = (abs(ix) <= ew.nReps && abs(iy) <= ew.nReps && abs(iz) <= ew.nReps);
 
                x = rx + ew.Lbox * ix;
                y = ry + ew.Lbox * iy;
                z = rz + ew.Lbox * iz;
                r2 = x*x + y*y + z*z;
                if (r2 >= ew.fEwCut2 && !bInHole) continue;

                dir = rsqrt(r2);
                dir2 = dir*dir;
                a = exp(-r2*ew.alpha2);
                a *= ew.ka*dir2;
                if (bInHole) g0 = -erf(ew.alpha*r2*dir);
                else         g0 = erfc(ew.alpha*r2*dir);
                g0 *= dir;
                g1 = g0*dir2 + a;
                alphan = 2*ew.alpha2;
                g2 = 3*g1*dir2 + alphan*a;
                alphan *= 2*ew.alpha2;
                g3 = 5*g2*dir2 + alphan*a;
                alphan *= 2*ew.alpha2;
                g4 = 7*g3*dir2 + alphan*a;
                alphan *= 2*ew.alpha2;
                g5 = 9*g4*dir2 + alphan*a;

                double onethird = 1.0/3.0;
                double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
                double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
                double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;

                xx = 0.5*x*x;
                xxx = onethird*xx*x;
                xxy = xx*y;
                xxz = xx*z;
                yy = 0.5*y*y;
                yyy = onethird*yy*y;
                xyy = yy*x;
                yyz = yy*z;
                zz = 0.5*z*z;
                zzz = onethird*zz*z;
                xzz = zz*x;
                yzz = zz*y;
                xy = x*y;
                xyz = xy*z;
                xz = x*z;
                yz = y*z;
                Q2mirx = ew.mom.xx*x + ew.mom.xy*y + ew.mom.xz*z;
                Q2miry = ew.mom.xy*x + ew.mom.yy*y + ew.mom.yz*z;
                Q2mirz = ew.mom.xz*x + ew.mom.yz*y + ew.mom.zz*z;
                Q3mirx = ew.mom.xxx*xx + ew.mom.xxy*xy + ew.mom.xxz*xz + ew.mom.xyy*yy + ew.mom.xyz*yz + ew.mom.xzz*zz;
                Q3miry = ew.mom.xxy*xx + ew.mom.xyy*xy + ew.mom.xyz*xz + ew.mom.yyy*yy + ew.mom.yyz*yz + ew.mom.yzz*zz;
                Q3mirz = ew.mom.xxz*xx + ew.mom.xyz*xy + ew.mom.xzz*xz + ew.mom.yyz*yy + ew.mom.yzz*yz + ew.mom.zzz*zz;
                Q4mirx = ew.mom.xxxx*xxx + ew.mom.xxxy*xxy + ew.mom.xxxz*xxz + ew.mom.xxyy*xyy + ew.mom.xxyz*xyz +
                    ew.mom.xxzz*xzz + ew.mom.xyyy*yyy + ew.mom.xyyz*yyz + ew.mom.xyzz*yzz + ew.mom.xzzz*zzz;
                Q4miry = ew.mom.xxxy*xxx + ew.mom.xxyy*xxy + ew.mom.xxyz*xxz + ew.mom.xyyy*xyy + ew.mom.xyyz*xyz +
                    ew.mom.xyzz*xzz + ew.mom.yyyy*yyy + ew.mom.yyyz*yyz + ew.mom.yyzz*yzz + ew.mom.yzzz*zzz;
                Q4mirz = ew.mom.xxxz*xxx + ew.mom.xxyz*xxy + ew.mom.xxzz*xxz + ew.mom.xyyz*xyy + ew.mom.xyzz*xyz +
                    ew.mom.xzzz*xzz + ew.mom.yyyz*yyy + ew.mom.yyzz*yyz + ew.mom.yzzz*yzz + ew.mom.zzzz*zzz;
                Q4x = ew.Q4xx*x + ew.Q4xy*y + ew.Q4xz*z;
                Q4y = ew.Q4xy*x + ew.Q4yy*y + ew.Q4yz*z;
                Q4z = ew.Q4xz*x + ew.Q4yz*y + ew.Q4zz*z;
                Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew.Q3x*x + ew.Q3y*y + ew.Q3z*z) + ew.Q4;
                Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
                Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
                Qta = g1*ew.mom.m - g2*ew.Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
                tpot -= g0*ew.mom.m - g1*ew.Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
                tax += g2*(Q2mirx - ew.Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
                tay += g2*(Q2miry - ew.Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
                taz += g2*(Q2mirz - ew.Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
                }
            }
	}

    // the H-Loop
    float fx=rx, fy=ry, fz=rz;
    float fax=0, fay=0, faz=0;
    for( i=0; i<ew.nEwhLoop; ++i) {
	float hdotx,s,c,t;
	hdotx = hx[i]*fx + hy[i]*fy + hz[i]*fz;
	sincosf(hdotx,&s,&c);
	tpot += hCfac[i]*c + hSfac[i]*s;
	t = hCfac[i]*s - hSfac[i]*c;
	fax += hx[i]*t;
	fay += hy[i]*t;
	faz += hz[i]*t;
	}
    X[pidx] = tax + fax;
    Y[pidx] = tay + fay;
    Z[pidx] = taz + faz;
    pPot[pidx] = tpot;
    }

extern "C"
void cudaEwaldInit(struct EwaldVariables *ewIn, EwaldTable *ewt ) {
    CUDA_CHECK(cudaMemcpyToSymbol,(ew,ewIn,sizeof(ew),0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hx,ewt->hx.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hy,ewt->hy.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hz,ewt->hz.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hCfac,ewt->hCfac.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hSfac,ewt->hSfac.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    }

extern "C"
int CUDAinitWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    //CUDACTX ctx = reinterpret_cast<CUDACTX>(e->cudaCtx);
    double *pHostBuf = reinterpret_cast<double *>(work->pHostBuf);
    double *pCudaBufIn = reinterpret_cast<double *>(work->pCudaBufIn);
    double *X, *Y, *Z;
    double *cudaX, *cudaY, *cudaZ, *cudaPot;
    int align, i;

    align = (e->nP+31)&~31; /* Warp align the memory buffers */
    X       = pHostBuf + 0*align;
    Y       = pHostBuf + 1*align;
    Z       = pHostBuf + 2*align;
    cudaX   = pCudaBufIn + 0*align;
    cudaY   = pCudaBufIn + 1*align;
    cudaZ   = pCudaBufIn + 2*align;
    cudaPot = pCudaBufIn + 3*align;

    dim3 dimBlock( 32, 1 );
    dim3 dimGrid( align/32, 1,1 );

    for(i=0; i<e->nP; ++i) {
        workParticle *wp = e->ppWorkPart[i];
        PARTICLE *p = wp->pPart[e->piWorkPart[i]];
	X[i] = p->r[0];
	Y[i] = p->r[1];
	Z[i] = p->r[2];
	}

    // copy data directly to device memory
    CUDA_CHECK(cudaMemcpyAsync,(pCudaBufIn, pHostBuf, align*3*sizeof(double),
	    cudaMemcpyHostToDevice, work->stream));
    cudaEwald<<<dimGrid, dimBlock, 0, work->stream>>>(cudaX,cudaY,cudaZ,cudaPot);
    CUDA_CHECK(cudaMemcpyAsync,(pHostBuf, pCudaBufIn, align*4*sizeof(double),
            cudaMemcpyDeviceToHost, work->stream));
    CUDA_CHECK(cudaEventRecord,(work->event,work->stream));

    return 1;
    }

extern "C"
void pkdAccumulateCUDA(void * pkd,workEwald *we,double *pax,double *pay,double *paz,double *pot);


extern "C"
int CUDAcheckWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    double *pHostBuf = reinterpret_cast<double *>(work->pHostBuf);
    double *X, *Y, *Z, *pPot;
    int align;

    align = (e->nP+31)&~31; /* As above! Warp align the memory buffers */
    X       = pHostBuf + 0*align;
    Y       = pHostBuf + 1*align;
    Z       = pHostBuf + 2*align;
    pPot    = pHostBuf + 3*align;
    pkdAccumulateCUDA(e->pkd,e,X,Y,Z,pPot);
    free(e->ppWorkPart);
    free(e->piWorkPart);
    free(e);
    return 0;

    }
