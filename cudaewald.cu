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

#define MAX_TOTAL_REPLICAS (7*7*7)

__constant__ struct EwaldVariables ew;
__constant__ float hx[MAX_TOTAL_REPLICAS];
__constant__ float hy[MAX_TOTAL_REPLICAS];
__constant__ float hz[MAX_TOTAL_REPLICAS];
__constant__ float hCfac[MAX_TOTAL_REPLICAS];
__constant__ float hSfac[MAX_TOTAL_REPLICAS];

__constant__ momFloat Lx[MAX_TOTAL_REPLICAS];
__constant__ momFloat Ly[MAX_TOTAL_REPLICAS];
__constant__ momFloat Lz[MAX_TOTAL_REPLICAS];
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

__global__ void cudaEwald(gpuEwaldInput *onGPU,gpuEwaldOutput *outGPU) {
    const momFloat f_1_2 = 1.0 / 2.0;
    const momFloat f_1_3 = 1.0 / 3.0;
    const momFloat f_1_4 = 1.0 / 4.0;
    const momFloat f_1_5 = 1.0 / 5.0;
    const momFloat f_1_7 = 1.0 / 7.0;
    const momFloat f_1_9 = 1.0 / 9.0;
    const momFloat f_1_11 = 1.0 / 11.0;
    const momFloat f_1_13 = 1.0 / 13.0;
    momFloat g0, g1, g2, g3, g4, g5, alphan;
    int i, bInHole;
    momFloat tax, tay, taz, dPot, dFlop=0.0;
    const momFloat rx = onGPU[blockIdx.x].X[threadIdx.x] - ew.r[0];
    const momFloat ry = onGPU[blockIdx.x].Y[threadIdx.x] - ew.r[1];
    const momFloat rz = onGPU[blockIdx.x].Z[threadIdx.x] - ew.r[2];

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
        const momFloat x = rx + Lx[i];
        const momFloat y = ry + Ly[i];
        const momFloat z = rz + Lz[i];
        momFloat r2 = x*x + y*y + z*z;
        if (r2 >= ew.fEwCut2 && !bInHole) continue;
        if (r2 < ew.fInner2) { /* Once, at most per particle */
            /*
             * For small r, series expand about
             * the origin to avoid errors caused
             * by cancellation of large terms.
             */
            alphan = ew.ka;
            r2 *= ew.alpha2;
            g0 = alphan*(f_1_3*r2 - 1);
            alphan *= 2*ew.alpha2;
	    g1 = alphan*(f_1_5*r2 - f_1_3);
            alphan *= 2*ew.alpha2;
	    g2 = alphan*(f_1_7*r2 - f_1_5);
            alphan *= 2*ew.alpha2;
	    g3 = alphan*(f_1_9*r2 - f_1_7);
            alphan *= 2*ew.alpha2;
	    g4 = alphan*(f_1_11*r2 - f_1_9);
            alphan *= 2*ew.alpha2;
	    g5 = alphan*(f_1_13*r2 - f_1_11);
            }
        else {
	    const momFloat dir = rsqrt(r2);
	    const momFloat dir2 = dir*dir;
	    const momFloat a = exp(-r2*ew.alpha2) * ew.ka*dir2;
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

        momFloat Q4mirx, Q4miry, Q4mirz;
        momFloat Q3mirx, Q3miry, Q3mirz;

	const  momFloat xx = f_1_2*x*x;
        Q3mirx = ew.mom.xxx*xx;
        Q3miry = ew.mom.xxy*xx;
        Q3mirz = ew.mom.xxz*xx;
	const  momFloat xxx = f_1_3*xx*x;
        Q4mirx = ew.mom.xxxx*xxx;
        Q4miry = ew.mom.xxxy*xxx;
        Q4mirz = ew.mom.xxxz*xxx;
        const  momFloat xxy = xx*y;
        Q4mirx += ew.mom.xxxy*xxy;
        Q4miry += ew.mom.xxyy*xxy;
        Q4mirz += ew.mom.xxyz*xxy;
        const  momFloat xxz = xx*z;
        Q4mirx += ew.mom.xxxz*xxz;
        Q4miry += ew.mom.xxyz*xxz;
        Q4mirz += ew.mom.xxzz*xxz;

	const  momFloat yy = f_1_2*y*y;
        Q3mirx += ew.mom.xyy*yy;
        Q3miry += ew.mom.yyy*yy;
        Q3mirz += ew.mom.yyz*yy;
        const  momFloat xyy = yy*x;
        Q4mirx += ew.mom.xxyy*xyy;
        Q4miry += ew.mom.xyyy*xyy;
        Q4mirz += ew.mom.xyyz*xyy;
	const  momFloat yyy = f_1_3*yy*y;
        Q4mirx += ew.mom.xyyy*yyy;
        Q4miry += ew.mom.yyyy*yyy;
        Q4mirz += ew.mom.yyyz*yyy;
        const  momFloat yyz = yy*z;
        Q4mirx += ew.mom.xyyz*yyz;
        Q4miry += ew.mom.yyyz*yyz;
        Q4mirz += ew.mom.yyzz*yyz;

        const  momFloat xy = x*y;
        Q3mirx += ew.mom.xxy*xy;
        Q3miry += ew.mom.xyy*xy;
        Q3mirz += ew.mom.xyz*xy;
        const  momFloat xyz = xy*z;
        Q4mirx += ew.mom.xxyz*xyz;
        Q4miry += ew.mom.xyyz*xyz;
        Q4mirz += ew.mom.xyzz*xyz;

	const  momFloat zz = f_1_2*z*z;
        Q3mirx += ew.mom.xzz*zz;
        Q3miry += ew.mom.yzz*zz;
        Q3mirz += ew.mom.zzz*zz;
        const  momFloat xzz = zz*x;
        Q4mirx += ew.mom.xxzz*xzz;
        Q4miry += ew.mom.xyzz*xzz;
        Q4mirz += ew.mom.xzzz*xzz;
        const  momFloat yzz = zz*y;
        Q4mirx += ew.mom.xyzz*yzz;
        Q4miry += ew.mom.yyzz*yzz;
        Q4mirz += ew.mom.yzzz*yzz;
	const  momFloat zzz = f_1_3*zz*z;
        Q4mirx += ew.mom.xzzz*zzz;
        Q4miry += ew.mom.yzzz*zzz;
        Q4mirz += ew.mom.zzzz*zzz;

        tax += g4*Q4mirx;
        tay += g4*Q4miry;
        taz += g4*Q4mirz;
	const momFloat Q4mir = f_1_4*(Q4mirx*x + Q4miry*y + Q4mirz*z);
        dPot -= g4*Q4mir;

        const  momFloat xz = x*z;
        Q3mirx += ew.mom.xxz*xz;
        Q3miry += ew.mom.xyz*xz;
        Q3mirz += ew.mom.xzz*xz;

        const  momFloat yz = y*z;
        Q3mirx += ew.mom.xyz*yz;
        Q3miry += ew.mom.yyz*yz;
        Q3mirz += ew.mom.yzz*yz;

        const momFloat Q4x = ew.Q4xx*x + ew.Q4xy*y + ew.Q4xz*z;
        const momFloat Q4y = ew.Q4xy*x + ew.Q4yy*y + ew.Q4yz*z;
        const momFloat Q4z = ew.Q4xz*x + ew.Q4yz*y + ew.Q4zz*z;
        const momFloat Q3mir = f_1_3*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
        dPot -= g3*Q3mir;
        tax += g3*(Q3mirx - Q4x);
        tay += g3*(Q3miry - Q4y);
        taz += g3*(Q3mirz - Q4z);

        const momFloat Q2mirx = ew.mom.xx*x + ew.mom.xy*y + ew.mom.xz*z;
        const momFloat Q2miry = ew.mom.xy*x + ew.mom.yy*y + ew.mom.yz*z;
        const momFloat Q2mirz = ew.mom.xz*x + ew.mom.yz*y + ew.mom.zz*z;
        const momFloat Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew.Q3x*x + ew.Q3y*y + ew.Q3z*z) + ew.Q4;
        dPot -= g2*Q2mir;
        tax += g2*(Q2mirx - ew.Q3x);
        tay += g2*(Q2miry - ew.Q3y);
        taz += g2*(Q2mirz - ew.Q3z);

        const momFloat Qta = g1*ew.mom.m - g2*ew.Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
        tax -= x*Qta;
        tay -= y*Qta;
        taz -= z*Qta;
        dFlop += COST_FLOP_EWALD;
	}

/*    dFlop += COST_FLOP_HLOOP * ew.nEwhLoop;*/ /* Accounted for outside */
    outGPU[blockIdx.x].X[threadIdx.x] = tax;
    outGPU[blockIdx.x].Y[threadIdx.x] = tay;
    outGPU[blockIdx.x].Z[threadIdx.x] = taz;
    outGPU[blockIdx.x].Pot[threadIdx.x]  = dPot;
    outGPU[blockIdx.x].Flop[threadIdx.x] = dFlop;
    }

static int doneEwaldFcn( void *vpp, void *vwork ) {
    return 0;
    }

static int initEwaldFcn(void *ve, void *vwork) {
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (ew, work->initEwald.ewIn, sizeof(ew), 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (hx, work->initEwald.ewt->hx.f, sizeof(float)*work->initEwald.ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (hy, work->initEwald.ewt->hy.f, sizeof(float)*work->initEwald.ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (hz, work->initEwald.ewt->hz.f, sizeof(float)*work->initEwald.ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (hCfac, work->initEwald.ewt->hCfac.f, sizeof(float)*work->initEwald.ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (hSfac, work->initEwald.ewt->hSfac.f, sizeof(float)*work->initEwald.ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, work->stream));
    // Time(%)      Time     Calls       Avg       Min       Max  Name
    // 14.93%  1.47255s       413  3.5655ms  2.6458ms  3.9733ms  cudaEwald(double*, double*, double*, double*, double*, double*, double*, d
    momFloat dLx[MAX_TOTAL_REPLICAS];
    momFloat dLy[MAX_TOTAL_REPLICAS];
    momFloat dLz[MAX_TOTAL_REPLICAS];
    int ibHole[MAX_TOTAL_REPLICAS];
    int i = 0, ix, iy, iz;
    for (ix = -3; ix <= 3; ++ix) {
	for (iy = -3; iy <= 3; ++iy) {
	    for (iz = -3; iz <= 3; ++iz) {
		ibHole[i] = (abs(ix) <= work->initEwald.ewIn->nReps && abs(iy) <= work->initEwald.ewIn->nReps && abs(iz) <= work->initEwald.ewIn->nReps);
		dLx[i] = work->initEwald.ewIn->Lbox * ix;
		dLy[i] = work->initEwald.ewIn->Lbox * iy;
		dLz[i] = work->initEwald.ewIn->Lbox * iz;
		++i;
	    }
	}
    }
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (Lx, dLx, sizeof(Lx), 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (Ly, dLy, sizeof(Ly), 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (Lz, dLz, sizeof(Lz), 0, cudaMemcpyHostToDevice, work->stream));
    CUDA_RETURN(cudaMemcpyToSymbolAsync, (bHole, ibHole, sizeof(bHole), 0, cudaMemcpyHostToDevice, work->stream));
    return cudaSuccess;
    }

#if 0
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
        momFloat dLx[MAX_TOTAL_REPLICAS];
        momFloat dLy[MAX_TOTAL_REPLICAS];
        momFloat dLz[MAX_TOTAL_REPLICAS];
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
#endif

extern "C"
void cudaEwaldInit(void *cudaCtx, struct EwaldVariables *ewIn, EwaldTable *ewt ) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    cuda->ewIn = ewIn;
    cuda->ewt = ewt;
    if (cuda->iCore==0) {
	assert(OPA_Queue_is_empty(&cuda->wqDone));
	CUDAwqNode *node = getNode(cuda);
        if (node==NULL) return;
	assert(node);
	node->initEwald.ewIn = cuda->ewIn;
	node->initEwald.ewt = cuda->ewt;
	node->initFcn = initEwaldFcn;
	node->dumpFcn = NULL;
	node->doneFcn = doneEwaldFcn;
	OPA_Queue_enqueue(cuda->queueWORK, node, CUDAwqNode, q.hdr);
	while (OPA_Queue_is_empty(&cuda->wqDone)) {}
        OPA_Queue_dequeue(&cuda->wqDone, node, CUDAwqNode, q.hdr);
	OPA_Queue_enqueue(&cuda->wqFree, node, CUDAwqNode, q.hdr);
	--cuda->nWorkQueueBusy;
//	cudaError_t ec = cuda_setup_ewald(cuda);
//        if (ec != cudaSuccess) CUDA_attempt_recovery(cuda,ec);
        }
    }

extern "C" void pkdParticleWorkDone(workParticle *work);

extern "C"
int CUDAdoneEwald( void *cudaCtx, void *vwork ) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDAwqNode *node = reinterpret_cast<CUDAwqNode *>(vwork);
    gpuEwaldOutput *fromGPU = reinterpret_cast<gpuEwaldOutput *>(node->pHostBufFromGPU);
    int iNode, i;
    int iResult = 0;
    for( iNode=0; iNode<node->ppnBuffered; ++iNode) {
	workParticle *wp = node->ppWP[iNode];
        for(i=0; i<wp->nP; ++i) {
            int ij = iResult / EWALD_ALIGN;
            int ii = iResult % EWALD_ALIGN;
            ++iResult;
            PINFOOUT *out = &wp->pInfoOut[i];
            out->a[0] += fromGPU[ij].X[ii];
            out->a[1] += fromGPU[ij].Y[ii];
            out->a[2] += fromGPU[ij].Z[ii];
            out->fPot += fromGPU[ij].Pot[ii];
            wp->dFlopSingleGPU += COST_FLOP_HLOOP*cuda->ewIn->nEwhLoop;
            wp->dFlopDoubleGPU += fromGPU[ij].Flop[ii];
            }
        pkdParticleWorkDone(wp);
        }
    assert(iResult == node->ewald.nParticles);
    return 0;
    }

extern "C"
int CUDAlaunchEwald( void *ve, void *vwork ) {
    CUDAwqNode *node = reinterpret_cast<CUDAwqNode *>(vwork);

    int align = (node->ewald.nParticles+EWALD_MASK)&~EWALD_MASK; /* Warp align the memory buffers */
    int ngrid = align/EWALD_ALIGN;
    dim3 dimBlock( EWALD_ALIGN, 1 );
    dim3 dimGrid( ngrid, 1,1 );
    gpuEwaldInput *toGPU = reinterpret_cast<gpuEwaldInput *>(node->pHostBufToGPU);
    gpuEwaldInput *onGPU = reinterpret_cast<gpuEwaldInput *>(node->pCudaBufIn);
    gpuEwaldOutput *outGPU = reinterpret_cast<gpuEwaldOutput *>(node->pCudaBufOut);
    gpuEwaldOutput *fromGPU  = reinterpret_cast<gpuEwaldOutput *>(node->pHostBufFromGPU);
    CUDA_RETURN(cudaMemcpyAsync,(onGPU, toGPU, ngrid * sizeof(gpuEwaldInput),
            cudaMemcpyHostToDevice, node->stream));
    cudaEwald<<<dimGrid, dimBlock, 0, node->stream>>>(onGPU,outGPU);
    CUDA_RETURN(cudaMemcpyAsync,(fromGPU, outGPU, ngrid * sizeof(gpuEwaldOutput),
            cudaMemcpyDeviceToHost, node->stream));
    return cudaSuccess;
    }

extern "C"
void CUDA_flushEwald(void *cudaCtx) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDAwqNode *node = cuda->nodeEwald;
    if (node && node->ppnBuffered) {
        node->startTime = CUDA_getTime();
        int ij = node->ewald.nParticles / EWALD_ALIGN;
        int ii = node->ewald.nParticles % EWALD_ALIGN;
        gpuEwaldInput *toGPU = reinterpret_cast<gpuEwaldInput *>(node->pHostBufToGPU);
        while(ii & EWALD_MASK) {
            toGPU[ij].X[ii] = toGPU[ij].X[0];
            toGPU[ij].Y[ii] = toGPU[ij].Y[0];
            toGPU[ij].Z[ii] = toGPU[ij].Z[0];
            ++ii;
            }
        OPA_Queue_enqueue(cuda->queueWORK, node, CUDAwqNode, q.hdr);
        cuda->nodeEwald = NULL;
        }
    }
extern "C"
int CUDA_queueEwald(void *cudaCtx,workParticle *work) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDAwqNode *node = cuda->nodeEwald;
    int i;

    // Make sure that there is enough buffer memory for this set of particles
    if (node!=NULL) {
        int nNewTotal = node->ewald.nParticles + work->nP;
        int align = (nNewTotal+EWALD_MASK)&~EWALD_MASK; /* Warp align the memory buffers */
        int nTiles = align / EWALD_ALIGN;
        if (nTiles*sizeof(gpuEwaldInput) > cuda->inCudaBufSize 
            || nTiles*sizeof(gpuEwaldOutput) > cuda->outCudaBufSize ) {
            CUDA_flushEwald(cuda);
            node = cuda->nodeEwald; // Should be NULL at this point
            }
        }

    if (node==NULL) {
        cuda->nodeEwald = node = getNode(cuda);
        if (node==NULL) return 0; // The CPU will handle this set of particles
        node->ewald.nParticles = 0;
        node->ppnBuffered = 0;
        node->ctx = cudaCtx;
        node->doneFcn = CUDAdoneEwald;
        node->initFcn = CUDAlaunchEwald;
        node->kernelName = "Ewald";
        }
    gpuEwaldInput *toGPU = reinterpret_cast<gpuEwaldInput *>(node->pHostBufToGPU);
    for( i=0; i<work->nP; i++ ) {
	const PINFOIN *in = &work->pInfoIn[i];
        int ij = node->ewald.nParticles / EWALD_ALIGN;
        int ii = node->ewald.nParticles % EWALD_ALIGN;
        toGPU[ij].X[ii] = work->c[0] + in->r[0];
	toGPU[ij].Y[ii] = work->c[1] + in->r[1];
	toGPU[ij].Z[ii] = work->c[2] + in->r[2];
        ++node->ewald.nParticles;
	}

    node->ppNI[node->ppnBuffered] = 0; // Unused
    node->ppWP[node->ppnBuffered] = work;
    ++work->nRefs;
    if ( ++node->ppnBuffered == CUDA_WP_MAX_BUFFERED) {
        CUDA_flushEwald(cuda);
        }

    return work->nP; // Return the number of particles sucessfully queued
    }
