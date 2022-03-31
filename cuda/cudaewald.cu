/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#define OLD_CUDA

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <stdio.h>
#include "basetype.h"
#include "gravity/moments.h"
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
    dPot = fPot + ew.mom.m*ew.k1;

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
    outGPU[blockIdx.x].FlopDouble[threadIdx.x] = dFlop;
    }

void pkdParticleWorkDone(workParticle *work);

/*****************************************************************************\
*   CudaClient interface (new!)
\*****************************************************************************/

// The Ewald routines need a table of contant values related to the replicates and moment
// This function will queue a message to the CUDA/MDL thread to setup the values.
void CudaClient::setupEwald(struct EwaldVariables * const ew, EwaldTable * const ewt) {
    nEwhLoop = ew->nEwhLoop;
    // The Ewald table needs to be sent once to every GPU on the system
    if (mdl.Core()==0) {
        mdl::cudaMessageQueue wait;
        for(auto i=0; i<mdl.numGPUs(); ++i) {
            auto m = new MessageEwaldSetup(ew,ewt,i);
            mdl.enqueue(*m,wait);
        }
        for(auto i=0; i<mdl.numGPUs(); ++i) {
            auto &m = wait.wait();
            delete &m;
        }
    }
    mdl.ThreadBarrier();
}

// Contruct the message with Ewald tables. We can just tuck away the pointers as we have to wait.
MessageEwaldSetup::MessageEwaldSetup(struct EwaldVariables *const ew, EwaldTable *const ewt, int iDevice)
    : mdl::cudaMessage(iDevice), ewIn(ew), ewt(ewt) {}

// This copies all of the variables to the device.
void MessageEwaldSetup::launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) {
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (ew,   ewIn, sizeof(ew), 0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (hx,   ewt->hx.f,   sizeof(float)*ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (hy,   ewt->hy.f,   sizeof(float)*ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (hz,   ewt->hz.f,   sizeof(float)*ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (hCfac,ewt->hCfac.f,sizeof(float)*ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (hSfac,ewt->hSfac.f,sizeof(float)*ewIn->nEwhLoop, 0, cudaMemcpyHostToDevice, stream));
    dLx.reserve(MAX_TOTAL_REPLICAS);     dLx.clear();
    dLy.reserve(MAX_TOTAL_REPLICAS);     dLy.clear();
    dLz.reserve(MAX_TOTAL_REPLICAS);     dLz.clear();
    ibHole.reserve(MAX_TOTAL_REPLICAS);  ibHole.clear();

    for (auto ix = -3; ix <= 3; ++ix) {
	for (auto iy = -3; iy <= 3; ++iy) {
	    for (auto iz = -3; iz <= 3; ++iz) {
		ibHole.push_back(abs(ix) <= ewIn->nReps && abs(iy) <= ewIn->nReps && abs(iz) <= ewIn->nReps);
		dLx.push_back(ewIn->Lbox * ix);
		dLy.push_back(ewIn->Lbox * iy);
		dLz.push_back(ewIn->Lbox * iz);
	    }
	}
    }
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (Lx,   dLx.data(),   dLx.size()*sizeof(decltype(dLx)::value_type),  0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (Ly,   dLy.data(),   dLy.size()*sizeof(decltype(dLy)::value_type),  0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (Lz,   dLz.data(),   dLz.size()*sizeof(decltype(dLz)::value_type),  0, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (bHole,ibHole.data(),ibHole.size()*sizeof(decltype(ibHole)::value_type), 0, cudaMemcpyHostToDevice, stream));
    }

extern "C"
int CudaClientQueueEwald(void *vcudaClient, workParticle *work) {
    auto cuda = reinterpret_cast<CudaClient *>(vcudaClient);
    return cuda->queueEwald(work);
    }

int CudaClient::queueEwald(workParticle *work) {
    if (ewald) {
    	if (ewald->queue(work)) return work->nP; // Sucessfully queued
    	mdl.enqueue(*ewald); // Full, so send it to the GPU
    	ewald = nullptr;
	}
    mdl.flushCompletedCUDA();
    if (freeEwald.empty()) return 0; // No buffers so the CPU has to do this part
    ewald = & freeEwald.dequeue();
    if (ewald->queue(work)) return work->nP; // Sucessfully queued
    return 0; // Not sure how this would happen, but okay.
    }

MessageEwald::MessageEwald(class CudaClient &cuda) : cuda(cuda) {
    nParticles = 0;
    nMaxParticles = requestBufferSize / sizeof(gpuEwaldInput) * EWALD_ALIGN;
    ppWP.reserve(CUDA_WP_MAX_BUFFERED);
    }

bool MessageEwald::queue(workParticle *work) {
    if (ppWP.size() == CUDA_WP_MAX_BUFFERED) return false; // Too many work packages
    if ( nParticles + work->nP > nMaxParticles) return false; // Not enough space
    auto toGPU = reinterpret_cast<gpuEwaldInput *>(pHostBufIn);
    for( auto i=0; i<work->nP; i++ ) {
	const PINFOIN *in = &work->pInfoIn[i];
        int ij = nParticles / EWALD_ALIGN;
        int ii = nParticles % EWALD_ALIGN;
        toGPU[ij].X[ii] = work->c[0] + in->r[0];
	toGPU[ij].Y[ii] = work->c[1] + in->r[1];
	toGPU[ij].Z[ii] = work->c[2] + in->r[2];
        ++nParticles;
	}
    ppWP.push_back(work);
    ++work->nRefs;

    return true;
    }

void MessageEwald::launch(cudaStream_t stream,void *pCudaBufIn, void *pCudaBufOut) {
    int align = (nParticles+EWALD_MASK)&~EWALD_MASK; /* Warp align the memory buffers */
    int ngrid = align/EWALD_ALIGN;
    dim3 dimBlock( EWALD_ALIGN, 1 );
    dim3 dimGrid( ngrid, 1,1 );
    gpuEwaldInput *toGPU = reinterpret_cast<gpuEwaldInput *>(pHostBufIn);
    gpuEwaldInput *onGPU = reinterpret_cast<gpuEwaldInput *>(pCudaBufIn);
    gpuEwaldOutput *outGPU = reinterpret_cast<gpuEwaldOutput *>(pCudaBufOut);
    gpuEwaldOutput *fromGPU  = reinterpret_cast<gpuEwaldOutput *>(pHostBufOut);
    CUDA_CHECK(cudaMemcpyAsync,(onGPU, toGPU, ngrid * sizeof(gpuEwaldInput), cudaMemcpyHostToDevice, stream));
    cudaEwald<<<dimGrid, dimBlock, 0, stream>>>(onGPU,outGPU);
    CUDA_CHECK(cudaMemcpyAsync,(fromGPU, outGPU, ngrid * sizeof(gpuEwaldOutput), cudaMemcpyDeviceToHost, stream));
    }

void MessageEwald::finish() {
    auto fromGPU = reinterpret_cast<gpuEwaldOutput *>(pHostBufOut);
    int iResult = 0;
    for( auto & wp : ppWP) {
        for(auto i=0; i<wp->nP; ++i) {
            int ij = iResult / EWALD_ALIGN;
            int ii = iResult % EWALD_ALIGN;
            ++iResult;
            PINFOOUT *out = &wp->pInfoOut[i];
            out->a[0] += fromGPU[ij].X[ii];
            out->a[1] += fromGPU[ij].Y[ii];
            out->a[2] += fromGPU[ij].Z[ii];
            out->fPot += fromGPU[ij].Pot[ii];
            wp->dFlopSingleGPU += COST_FLOP_HLOOP * cuda.nEwhLoop;
            wp->dFlopDoubleGPU += fromGPU[ij].FlopDouble[ii];
            }
        pkdParticleWorkDone(wp);
        }
    assert(iResult == nParticles);
    nParticles = 0;
    ppWP.clear();
    cuda.freeEwald.enqueue(*this);
    }
