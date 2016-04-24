/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <stdio.h>
#include "basetype.h"
#include "moments.h"
#include "cudautil.h"

#define ALIGN 64
#define MASK (ALIGN-1)

#define MAX_TOTAL_REPLICAS (7*7*7)

__constant__ struct EwaldVariables ew;
__constant__ float hx[MAX_TOTAL_REPLICAS];
__constant__ float hy[MAX_TOTAL_REPLICAS];
__constant__ float hz[MAX_TOTAL_REPLICAS];
__constant__ float hCfac[MAX_TOTAL_REPLICAS];
__constant__ float hSfac[MAX_TOTAL_REPLICAS];

__constant__ double Lx[MAX_TOTAL_REPLICAS];
__constant__ double Ly[MAX_TOTAL_REPLICAS];
__constant__ double Lz[MAX_TOTAL_REPLICAS];
__constant__ int bHole[MAX_TOTAL_REPLICAS];


/*
** nvcc -DHAVE_CONFIG_H --ptxas-options=-v -c  -I. -arch=sm_35 cudaewald.cu
** ptxas info    : 0 bytes gmem, 16948 bytes cmem[3]
** ptxas info    : Compiling entry function '_Z9cudaEwaldPdS_S_S_S_S_S_S_' for 'sm_35'
** ptxas info    : Function properties for _Z9cudaEwaldPdS_S_S_S_S_S_S_
**    32 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
** ptxas info    : Used 79 registers, 384 bytes cmem[0], 668 bytes cmem[2]
*/


/*
** threadIdx.x: all work on the same particle -- this is the warp size, i.e., 32
** blockIdx.x:  different particles. If y=z=1, then x can be anything, otherwise
**              the total number of particles is a block of x*y*z
*/

__global__ void cudaEwald(double *X,double *Y,double *Z,
    double *Xout, double *Yout, double *Zout, double *pPot,double *pdFlop) {
    int pidx = threadIdx.x + ALIGN*blockIdx.x;

    const double onethird = 1.0/3.0;
    double g0,g1,g2,g3,g4,g5,alphan;
    int i, bInHole;
    double tax, tay, taz, dPot, dFlop=0.0;
    const double rx = X[pidx] - ew.r[0];
    const double ry = Y[pidx] - ew.r[1];
    const double rz = Z[pidx] - ew.r[2];

    // the H-Loop
    float fx=rx, fy=ry, fz=rz;
    float fax=0, fay=0, faz=0, fPot=0;
    for( i=0; i<ew.nEwhLoop; ++i) {
	float hdotx,s,c,t;
	hdotx = hx[i]*fx + hy[i]*fy + hz[i]*fz;
	sincosf(hdotx,&s,&c);
	fPot += hCfac[i]*c + hSfac[i]*s;
	t = hCfac[i]*s - hSfac[i]*c;
	fax += hx[i]*t;
	fay += hy[i]*t;
	faz += hz[i]*t;
	}
    tax = fax;
    tay = fay;
    taz = faz;
    dPot = fPot;

    for(i=0; i<MAX_TOTAL_REPLICAS; ++i) {
        bInHole = bHole[i];
        const double x = rx + Lx[i];
        const double y = ry + Ly[i];
        const double z = rz + Lz[i];
        double r2 = x*x + y*y + z*z;
        if (r2 >= ew.fEwCut2 && !bInHole) continue;
        if (r2 < ew.fInner2) { /* Once, at most per particle */
            /*
             * For small r, series expand about
             * the origin to avoid errors caused
             * by cancellation of large terms.
             */
            alphan = ew.ka;
            r2 *= ew.alpha2;
            g0 = alphan*((1.0/3.0)*r2 - 1.0);
            alphan *= 2*ew.alpha2;
            g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
            alphan *= 2*ew.alpha2;
            g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
            alphan *= 2*ew.alpha2;
            g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
            alphan *= 2*ew.alpha2;
            g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
            alphan *= 2*ew.alpha2;
            g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
            }
        else {
            const double dir = rsqrt(r2);
            const double dir2 = dir*dir;
            const double a = exp(-r2*ew.alpha2) * ew.ka*dir2;
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
            }

        dPot -= g0*ew.mom.m - g1*ew.Q2;

        double Q4mirx, Q4miry, Q4mirz;
        double Q3mirx, Q3miry, Q3mirz;

        const  double xx = 0.5*x*x;
        Q3mirx = ew.mom.xxx*xx;
        Q3miry = ew.mom.xxy*xx;
        Q3mirz = ew.mom.xxz*xx;
        const  double xxx = onethird*xx*x;
        Q4mirx = ew.mom.xxxx*xxx;
        Q4miry = ew.mom.xxxy*xxx;
        Q4mirz = ew.mom.xxxz*xxx;
        const  double xxy = xx*y;
        Q4mirx += ew.mom.xxxy*xxy;
        Q4miry += ew.mom.xxyy*xxy;
        Q4mirz += ew.mom.xxyz*xxy;
        const  double xxz = xx*z;
        Q4mirx += ew.mom.xxxz*xxz;
        Q4miry += ew.mom.xxyz*xxz;
        Q4mirz += ew.mom.xxzz*xxz;

        const  double yy = 0.5*y*y;
        Q3mirx += ew.mom.xyy*yy;
        Q3miry += ew.mom.yyy*yy;
        Q3mirz += ew.mom.yyz*yy;
        const  double xyy = yy*x;
        Q4mirx += ew.mom.xxyy*xyy;
        Q4miry += ew.mom.xyyy*xyy;
        Q4mirz += ew.mom.xyyz*xyy;
        const  double yyy = onethird*yy*y;
        Q4mirx += ew.mom.xyyy*yyy;
        Q4miry += ew.mom.yyyy*yyy;
        Q4mirz += ew.mom.yyyz*yyy;
        const  double yyz = yy*z;
        Q4mirx += ew.mom.xyyz*yyz;
        Q4miry += ew.mom.yyyz*yyz;
        Q4mirz += ew.mom.yyzz*yyz;

        const  double xy = x*y;
        Q3mirx += ew.mom.xxy*xy;
        Q3miry += ew.mom.xyy*xy;
        Q3mirz += ew.mom.xyz*xy;
        const  double xyz = xy*z;
        Q4mirx += ew.mom.xxyz*xyz;
        Q4miry += ew.mom.xyyz*xyz;
        Q4mirz += ew.mom.xyzz*xyz;

        const  double zz = 0.5*z*z;
        Q3mirx += ew.mom.xzz*zz;
        Q3miry += ew.mom.yzz*zz;
        Q3mirz += ew.mom.zzz*zz;
        const  double xzz = zz*x;
        Q4mirx += ew.mom.xxzz*xzz;
        Q4miry += ew.mom.xyzz*xzz;
        Q4mirz += ew.mom.xzzz*xzz;
        const  double yzz = zz*y;
        Q4mirx += ew.mom.xyzz*yzz;
        Q4miry += ew.mom.yyzz*yzz;
        Q4mirz += ew.mom.yzzz*yzz;
        const  double zzz = onethird*zz*z;
        Q4mirx += ew.mom.xzzz*zzz;
        Q4miry += ew.mom.yzzz*zzz;
        Q4mirz += ew.mom.zzzz*zzz;

        tax += g4*Q4mirx;
        tay += g4*Q4miry;
        taz += g4*Q4mirz;
        const double Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
        dPot -= g4*Q4mir;

        const  double xz = x*z;
        Q3mirx += ew.mom.xxz*xz;
        Q3miry += ew.mom.xyz*xz;
        Q3mirz += ew.mom.xzz*xz;

        const  double yz = y*z;
        Q3mirx += ew.mom.xyz*yz;
        Q3miry += ew.mom.yyz*yz;
        Q3mirz += ew.mom.yzz*yz;

        const double Q4x = ew.Q4xx*x + ew.Q4xy*y + ew.Q4xz*z;
        const double Q4y = ew.Q4xy*x + ew.Q4yy*y + ew.Q4yz*z;
        const double Q4z = ew.Q4xz*x + ew.Q4yz*y + ew.Q4zz*z;
        const double Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
        dPot -= g3*Q3mir;
        tax += g3*(Q3mirx - Q4x);
        tay += g3*(Q3miry - Q4y);
        taz += g3*(Q3mirz - Q4z);

        const double Q2mirx = ew.mom.xx*x + ew.mom.xy*y + ew.mom.xz*z;
        const double Q2miry = ew.mom.xy*x + ew.mom.yy*y + ew.mom.yz*z;
        const double Q2mirz = ew.mom.xz*x + ew.mom.yz*y + ew.mom.zz*z;
        const double Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew.Q3x*x + ew.Q3y*y + ew.Q3z*z) + ew.Q4;
        dPot -= g2*Q2mir;
        tax += g2*(Q2mirx - ew.Q3x);
        tay += g2*(Q2miry - ew.Q3y);
        taz += g2*(Q2mirz - ew.Q3z);

        const double Qta = g1*ew.mom.m - g2*ew.Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
        tax -= x*Qta;
        tay -= y*Qta;
        taz -= z*Qta;
        dFlop += COST_FLOP_EWALD;
	}

/*    dFlop += COST_FLOP_HLOOP * ew.nEwhLoop;*/ /* Accounted for outside */
    Xout[pidx] = tax;
    Yout[pidx] = tay;
    Zout[pidx] = taz;
    pPot[pidx] = dPot;
    pdFlop[pidx] = dFlop;
    }

/* If this returns an error, then the caller must attempt recovery or abort */
cudaError_t cuda_setup_ewald(CUDACTX cuda) {
    if (cuda->ewIn && cuda->ewt) {
        double start = CUDA_getTime();
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(ew,cuda->ewIn,sizeof(ew),0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(hx,cuda->ewt->hx.f,sizeof(float)*cuda->ewIn->nEwhLoop,0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(hy,cuda->ewt->hy.f,sizeof(float)*cuda->ewIn->nEwhLoop,0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(hz,cuda->ewt->hz.f,sizeof(float)*cuda->ewIn->nEwhLoop,0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(hCfac,cuda->ewt->hCfac.f,sizeof(float)*cuda->ewIn->nEwhLoop,0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(hSfac,cuda->ewt->hSfac.f,sizeof(float)*cuda->ewIn->nEwhLoop,0,cudaMemcpyHostToDevice,cuda->streamEwald));
// Time(%)      Time     Calls       Avg       Min       Max  Name
// 14.93%  1.47255s       413  3.5655ms  2.6458ms  3.9733ms  cudaEwald(double*, double*, double*, double*, double*, double*, double*, d
        double dLx[MAX_TOTAL_REPLICAS];
        double dLy[MAX_TOTAL_REPLICAS];
        double dLz[MAX_TOTAL_REPLICAS];
        int ibHole[MAX_TOTAL_REPLICAS];
        int i=0, ix, iy, iz;
        for(ix=-3; ix<=3; ++ix) {
            for(iy=-3; iy<=3; ++iy) {
                for(iz=-3; iz<=3; ++iz) {
                    ibHole[i] = (abs(ix) <= cuda->ewIn->nReps && abs(iy) <= cuda->ewIn->nReps && abs(iz) <= cuda->ewIn->nReps);
                    dLx[i] = cuda->ewIn->Lbox * ix;
                    dLy[i] = cuda->ewIn->Lbox * iy;
                    dLz[i] = cuda->ewIn->Lbox * iz;
                    ++i;
                    }
                }
            }
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(Lx,dLx,sizeof(Lx),0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(Ly,dLy,sizeof(Ly),0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(Lz,dLz,sizeof(Lz),0,cudaMemcpyHostToDevice,cuda->streamEwald));
        CUDA_RETURN(cudaMemcpyToSymbolAsync,(bHole,ibHole,sizeof(bHole),0,cudaMemcpyHostToDevice,cuda->streamEwald));


#ifdef USE_CUDA_EVENTS
        CUDA_RETURN(cudaEventRecord,(cuda->eventEwald,cuda->streamEwald));
#endif
        cudaError_t rc;
        do {
#ifdef USE_CUDA_EVENTS
            rc = cudaEventQuery(cuda->eventEwald);
#else
            rc = cudaStreamQuery(cuda->streamEwald);
#endif
            switch(rc) {
            case cudaSuccess:
            case cudaErrorNotReady:
                break;
            default:
                return rc;
                }
            if (CUDA_getTime() - start > 1.0) {
                return cudaErrorLaunchTimeout;
                }
            } while (rc!=cudaSuccess);
        }
    return cudaSuccess;
    }

extern "C"
void cudaEwaldInit(void *cudaCtx, struct EwaldVariables *ewIn, EwaldTable *ewt ) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    cuda->ewIn = ewIn;
    cuda->ewt = ewt;
    if (cuda->iCore==0) {
        cudaError_t ec = cuda_setup_ewald(cuda);
        if (ec != cudaSuccess) CUDA_attempt_recovery(cuda,ec);
        }
    }

/* If this returns an error, then the caller must attempt recovery or abort */
extern "C"
int CUDAinitWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    double *pHostBufFromGPU  = reinterpret_cast<double *>(work->pHostBufFromGPU);
    double *pHostBufToGPU    = reinterpret_cast<double *>(work->pHostBufToGPU);
    double *pCudaBufIn = reinterpret_cast<double *>(work->pCudaBufIn);
    double *pCudaBufOut = reinterpret_cast<double *>(work->pCudaBufOut);
    double *X, *Y, *Z;
    double *cudaX, *cudaY, *cudaZ, *cudaXout, *cudaYout, *cudaZout, *cudaPot, *cudaFlop;
    int align, i;

    align = (e->nP+MASK)&~MASK; /* Warp align the memory buffers */
    X       = pHostBufToGPU + 0*align;
    Y       = pHostBufToGPU + 1*align;
    Z       = pHostBufToGPU + 2*align;
    cudaX   = pCudaBufIn + 0*align;
    cudaY   = pCudaBufIn + 1*align;
    cudaZ   = pCudaBufIn + 2*align;
    cudaXout= pCudaBufOut + 0*align;
    cudaYout= pCudaBufOut + 1*align;
    cudaZout= pCudaBufOut + 2*align;
    cudaPot = pCudaBufOut + 3*align;
    cudaFlop= pCudaBufOut + 4*align;

    dim3 dimBlock( ALIGN, 1 );
    dim3 dimGrid( align/ALIGN, 1,1 );
    for(i=0; i<e->nP; ++i) {
        const workParticle *wp = e->ppWorkPart[i];
	const int wi = e->piWorkPart[i];
	const PINFOIN *in = &wp->pInfoIn[wi];
	X[i] = wp->c[0] + in->r[0];
	Y[i] = wp->c[1] + in->r[1];
	Z[i] = wp->c[2] + in->r[2];
	}
    for(;i<align;++i) X[i]=Y[i]=Z[i] = 100;

    // copy data directly to device memory
    CUDA_RETURN(cudaMemcpyAsync,(pCudaBufIn, pHostBufToGPU, align*3*sizeof(double),
	    cudaMemcpyHostToDevice, work->stream));
    cudaEwald<<<dimGrid, dimBlock, 0, work->stream>>>(cudaX,cudaY,cudaZ,cudaXout,cudaYout,cudaZout,cudaPot,cudaFlop);
    CUDA_RETURN(cudaMemcpyAsync,(pHostBufFromGPU, pCudaBufOut, align*5*sizeof(double),
            cudaMemcpyDeviceToHost, work->stream));
#ifdef USE_CUDA_EVENTS
    CUDA_RETURN(cudaEventRecord,(work->event,work->stream));
#endif

    return cudaSuccess;
    }

extern "C"
void pkdAccumulateCUDA(void * pkd,workEwald *we,double *pax,double *pay,double *paz,double *pot,double *pdFlop);


extern "C"
int CUDAcheckWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    double *pHostBuf = reinterpret_cast<double *>(work->pHostBufFromGPU);
    double *X, *Y, *Z, *pPot, *pdFlop;
    int align;

    align = (e->nP+MASK)&~MASK; /* As above! Warp align the memory buffers */
    X       = pHostBuf + 0*align;
    Y       = pHostBuf + 1*align;
    Z       = pHostBuf + 2*align;
    pPot    = pHostBuf + 3*align;
    pdFlop  = pHostBuf + 4*align;
    pkdAccumulateCUDA(e->pkd,e,X,Y,Z,pPot,pdFlop);
    free(e->ppWorkPart);
    free(e->piWorkPart);
    free(e);
    return 0;

    }
