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

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <float.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <assert.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include "pkd.h"
#include "moments.h"
#include "meval.h"
#include "qeval.h"
#include "ewald.h"
#include "grav.h"
#include "cudautil.h"
#include "cudapppc.h"
#include "cudaewald.h"
#include "clutil.h"
#include "SPHOptions.h"

#if 1
#if defined(USE_SIMD) && defined(__SSE2__)
/* Caution: This uses v/sqrt(v) so v cannot be zero! */
static inline float asqrtf(float v) {
    __m128 r2 = _mm_set_ss(v);
    __m128 r = _mm_rsqrt_ps(r2);
    r = _mm_mul_ss(r,_mm_sub_ss(_mm_set_ss(3.0/2.0),_mm_mul_ss(_mm_mul_ss(r,r),_mm_mul_ss(r2,_mm_set_ss(0.5)))));
    r = _mm_mul_ss(r,r2);
    v = _mm_cvtss_f32(r);
    return v;
    }
static inline float rsqrtf(float v) {
    __m128 r2 = _mm_set_ss(v);
    __m128 r = _mm_rsqrt_ps(r2);
    r = _mm_mul_ss(r,_mm_sub_ss(_mm_set_ss(3.0/2.0),_mm_mul_ss(_mm_mul_ss(r,r),_mm_mul_ss(r2,_mm_set_ss(0.5)))));
    v =_mm_cvtss_f32(r);
    return v;
    }
#else
static inline float asqrtf(float v) {
    return sqrtf(v);
    }
static inline float rsqrtf(float v) {
    return 1.0f / sqrtf(v);
    }
#endif
#endif

#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    }

/*
** This is called after work has been done for this particle group.
** If everyone has finished, then the particle is updated.
*/
void pkdParticleWorkDone(workParticle *wp) {
    PKD pkd = wp->ctx;
    int i,gid;
    PARTICLE *p;
    double r[3];
    vel_t *v,v2;
    float *a;
    float fx, fy, fz, m;
    float maga, dT, dtGrav;
    unsigned char uNewRung;

    if ( --wp->nRefs == 0 ) {
	float fiDelta = 1.0/wp->ts->dDelta;
	float fEta = wp->ts->dEta;
	float fiAccFac = 1.0 / wp->ts->dAccFac;
	pkd->dFlop += wp->dFlop;
	pkd->dFlopSingleCPU += wp->dFlopSingleCPU;
	pkd->dFlopDoubleCPU += wp->dFlopDoubleCPU;
	pkd->dFlopSingleGPU += wp->dFlopSingleGPU;
	pkd->dFlopDoubleGPU += wp->dFlopDoubleGPU;
	for( i=0; i<wp->nP; i++ ) {
	    p = wp->pPart[i];

        if (wp->SPHoptions->doDensity) {
        //SPHFIELDS *pSph = pkdSph(pkd,p);
        pkdSetDensity(pkd,p,wp->pInfoOut[i].rho);
        pkdSetBall(pkd,p,wp->pInfoOut[i].fBall);
        //pSph->Omega = 1.0f + wp->pInfoOut[i].fBall/(3.0f * wp->pInfoOut[i].rho)*wp->pInfoOut[i].drhodfball;
        //pSph->divv = wp->pInfoOut[i].divv;
        //pSph->thetaDot = wp->pInfoOut[i].thetaDot;
        }

        if (wp->SPHoptions->doGravity) {
	    pkdGetPos1(pkd,p,r);
	    m = pkdMass(pkd,p);
	    if (pkd->oFieldOffset[oAcceleration]) {
		a = pkdAccel(pkd,p);
		a[0] = wp->pInfoOut[i].a[0];
		a[1] = wp->pInfoOut[i].a[1];
		a[2] = wp->pInfoOut[i].a[2];
		}
	    if (pkd->oFieldOffset[oPotential]) {
		float *pPot = pkdPot(pkd,p);
		*pPot = wp->pInfoOut[i].fPot;
		}
	    if (pkd->ga != NULL) {
		gid = pkdGetGroup(pkd,p);
		if (gid && wp->pInfoOut[i].fPot < pkd->ga[gid].minPot) {
		    pkd->ga[gid].minPot = wp->pInfoOut[i].fPot;
		    pkd->ga[gid].iMinPart = wp->iPart[i];
		    }
		}
	    a = wp->pInfoOut[i].a;
	    pkd->dEnergyU += 0.5 * m * wp->pInfoOut[i].fPot;
	    pkd->dEnergyW += m*(r[0]*a[0] + r[1]*a[1] + r[2]*a[2]);
	    pkd->dEnergyF[0] += m*a[0];
	    pkd->dEnergyF[1] += m*a[1];
	    pkd->dEnergyF[2] += m*a[2];

	    if (wp->ts->bGravStep) {
		float dirsum = wp->pInfoOut[i].dirsum;
		float normsum = wp->pInfoOut[i].normsum;

		/*
		** If this is the first time through, the accelerations will have 
		** all been zero resulting in zero for normsum (and nan for dtGrav).
		** We repeat this process again, so dtGrav will be correct.
		*/
		if (normsum > 0.0) {
		    /*
		    ** Use new acceleration here!
		    */
		    fx = wp->pInfoOut[i].a[0];
		    fy = wp->pInfoOut[i].a[1];
		    fz = wp->pInfoOut[i].a[2];
		    maga = fx*fx + fy*fy + fz*fz;
		    if (maga>0.0f) maga = asqrtf(maga);
		    dtGrav = maga*dirsum/normsum;
		    }
		else dtGrav = 0.0;
		dtGrav += wp->ts->dPreFacRhoLoc*wp->pInfoIn[i].fDensity;
		dtGrav = (wp->pInfoOut[i].rhopmax > dtGrav?wp->pInfoOut[i].rhopmax:dtGrav);
		if (dtGrav > 0.0) {
		    dT = fEta * rsqrtf(dtGrav*wp->ts->dRhoFac);
            dT = fmin(dT,wp->pInfoOut[i].dtEst);
		    uNewRung = pkdDtToRungInverse(dT,fiDelta,wp->ts->uMaxRung-1);
		    }
		else uNewRung = 0; /* Assumes current uNewRung is outdated -- not ideal */
		} /* end of wp->bGravStep */
	    else {
		/*
		** We are doing eps/a timestepping.
		*/
		fx = wp->pInfoOut[i].a[0];
		fy = wp->pInfoOut[i].a[1];
		fz = wp->pInfoOut[i].a[2];
		maga = fx*fx + fy*fy + fz*fz;
		if (maga > 0) {
		    float imaga = rsqrtf(maga) * fiAccFac;
		    dT = fEta*asqrtf(pkdSoft(pkd,p)*imaga);
            dT = fmin(dT,wp->pInfoOut[i].dtEst);
		    uNewRung = pkdDtToRungInverse(dT,fiDelta,wp->ts->uMaxRung-1);
		    }
		else uNewRung = 0;
		}
	    /*
	    ** Here we must make sure we do not try to take a larger opening
	    ** timestep than the current active particles involved in the
	    ** gravity calculation!
	    */
	    if (uNewRung < wp->ts->uRungLo) uNewRung = wp->ts->uRungLo;
	    else if (uNewRung > wp->ts->uRungHi) uNewRung = wp->ts->uRungHi; 
	    if (!pkd->bNoParticleOrder) p->uNewRung = uNewRung;
	    /*
	    ** Now we want to kick the particle velocities with a closing kick 
	    ** based on the old rung and an opening kick based on the new rung.
	    ** However, we are not always allowed to decrease uNewRung by an
	    ** arbitrary amount, as this depends on where we are in the 
	    ** substep hierarchy.
	    ** This also requires the various timestep factors to have been 
	    ** pre-calculated for this dTime and all possible rungs at play.
	    ** Note that the only persistent info needed is the old rung which 
	    ** gets updated at the end of this so that p->uRung indicates 
	    ** which particles are active in this time step. 
	    ** p->uRung = uNewRung could be done here, but we let the outside 
	    ** code handle this.
	    */
	    if (pkd->oFieldOffset[oVelocity]) {
		v = pkdVel(pkd,p);
		if (wp->kick && wp->kick->bKickClose) {
		    v[0] += wp->kick->dtClose[p->uRung]*wp->pInfoOut[i].a[0];
		    v[1] += wp->kick->dtClose[p->uRung]*wp->pInfoOut[i].a[1];
		    v[2] += wp->kick->dtClose[p->uRung]*wp->pInfoOut[i].a[2];
		    }
		v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
		/*
		** Now calculate the kinetic energy term.
		*/
		pkd->dEnergyT += 0.5*m*v2;
		/* L is calculated with respect to the origin (0,0,0) */
		pkd->dEnergyL[0] += m*(r[1]*v[2] - r[2]*v[1]);
		pkd->dEnergyL[1] += m*(r[2]*v[0] - r[0]*v[2]);
		pkd->dEnergyL[2] += m*(r[0]*v[1] - r[1]*v[0]);

		if (wp->kick && wp->kick->bKickOpen) {
		    p->uRung = uNewRung;
		    ++pkd->nRung[p->uRung];
		    v[0] += wp->kick->dtOpen[p->uRung]*wp->pInfoOut[i].a[0];
		    v[1] += wp->kick->dtOpen[p->uRung]*wp->pInfoOut[i].a[1];
		    v[2] += wp->kick->dtOpen[p->uRung]*wp->pInfoOut[i].a[2];		    
		    /*
		    ** On KickOpen we also always check for intersection with the lightcone
		    ** surface over the entire next timestep of the particle (not a half
		    ** timestep as is usual for kicking (we are drifting afterall).
		    */
		    if (wp->lc->dLookbackFac > 0) {
			pkdProcessLightCone(pkd,p,wp->pInfoOut[i].fPot,wp->lc->dLookbackFac,wp->lc->dLookbackFacLCP,
				wp->lc->dtLCDrift[p->uRung],wp->lc->dtLCKick[p->uRung],
				wp->lc->dBoxSize,wp->lc->bLightConeParticles);
			}
		    }
		}
        } else {
            ++pkd->nRung[p->uRung];
        }
	    }
	free(wp->pPart);
	free(wp->iPart);
	free(wp->pInfoIn);
	free(wp->pInfoOut);
	free(wp);
	}
    }

int CPUdoWorkPP(void *vpp) {
    workPP *pp = vpp;
    workParticle *wp = pp->work;
    ILPTILE tile = pp->tile;
    ILP_BLK *blk = tile->blk;
    PINFOIN *pPart = &pp->work->pInfoIn[pp->i];
    PINFOOUT *pOut = &pp->pInfoOut[pp->i];
    int nBlocks = tile->lstTile.nBlocks;
    int nInLast = tile->lstTile.nInLast;

    pOut->a[0] = 0.0;
    pOut->a[1] = 0.0;
    pOut->a[2] = 0.0;
    pOut->fPot = 0.0;
    pOut->dirsum = 0.0;
    pOut->normsum = 0.0;
    pkdGravEvalPP(pPart,nBlocks,nInLast,blk,pOut);
    wp->dFlopSingleCPU += COST_FLOP_PP*(tile->lstTile.nBlocks*ILP_PART_PER_BLK  + tile->lstTile.nInLast);
    if ( ++pp->i == pp->work->nP ) return 0;
    else return 1;
    }

int CPUdoWorkDensity(void *vpp) {
    workPP *pp = vpp;
    workParticle *wp = pp->work;
    ILPTILE tile = pp->tile;
    ILP_BLK *blk = tile->blk;
    PINFOIN *pPart = &pp->work->pInfoIn[pp->i];
    PINFOOUT *pOut = &pp->pInfoOut[pp->i];
    int nBlocks = tile->lstTile.nBlocks;
    int nInLast = tile->lstTile.nInLast;
    SPHOptions *SPHoptions = wp->SPHoptions;

    pOut->rho = 0.0;
    pOut->drhodfball = 0.0;
    pOut->nden = 0.0;
    pOut->dndendfball = 0.0;
    pOut->nSmooth = 0.0;
    pkdDensityEval(pPart,nBlocks,nInLast,blk,pOut,SPHoptions);
    //wp->dFlopSingleCPU += COST_FLOP_PP*(tile->lstTile.nBlocks*ILP_PART_PER_BLK  + tile->lstTile.nInLast);
    if ( ++pp->i == pp->work->nP ) return 0;
    else return 1;
    }

int CPUdoWorkSPHForces(void *vpp) {
    workPP *pp = vpp;
    workParticle *wp = pp->work;
    ILPTILE tile = pp->tile;
    ILP_BLK *blk = tile->blk;
    PINFOIN *pPart = &pp->work->pInfoIn[pp->i];
    PINFOOUT *pOut = &pp->pInfoOut[pp->i];
    int nBlocks = tile->lstTile.nBlocks;
    int nInLast = tile->lstTile.nInLast;
    SPHOptions *SPHoptions = wp->SPHoptions;

    pOut->thetaDot = 0.0;
    pOut->a[0] = 0.0;
    pOut->a[1] = 0.0;
    pOut->a[2] = 0.0;
    pOut->divv = 0.0;
    pOut->dtEst = HUGE_VAL;
    pkdSPHForcesEval(pPart,nBlocks,nInLast,blk,pOut,SPHoptions);
    //wp->dFlopSingleCPU += COST_FLOP_PP*(tile->lstTile.nBlocks*ILP_PART_PER_BLK  + tile->lstTile.nInLast);
    if ( ++pp->i == pp->work->nP ) return 0;
    else return 1;
    }

int doneWorkPP(void *vpp) {
    workPP *pp = vpp;
    int i;

    for(i=0; i<pp->work->nP; ++i) {
	pp->work->pInfoOut[i].a[0] += pp->pInfoOut[i].a[0];
	pp->work->pInfoOut[i].a[1] += pp->pInfoOut[i].a[1];
	pp->work->pInfoOut[i].a[2] += pp->pInfoOut[i].a[2];
	pp->work->pInfoOut[i].fPot += pp->pInfoOut[i].fPot;
	pp->work->pInfoOut[i].dirsum += pp->pInfoOut[i].dirsum;
	pp->work->pInfoOut[i].normsum += pp->pInfoOut[i].normsum;
	}
    lstFreeTile(&pp->ilp->lst,&pp->tile->lstTile);
    pkdParticleWorkDone(pp->work);
    free(pp->pInfoOut);
    free(pp);
    return 0;
    }

int doneWorkDensity(void *vpp) {
    workPP *pp = vpp;
    int i;

    for(i=0; i<pp->work->nP; ++i) {
	pp->work->pInfoOut[i].rho += pp->pInfoOut[i].rho;
	pp->work->pInfoOut[i].drhodfball += pp->pInfoOut[i].drhodfball;
    pp->work->pInfoOut[i].nden += pp->pInfoOut[i].nden;
	pp->work->pInfoOut[i].dndendfball += pp->pInfoOut[i].dndendfball;
    pp->work->pInfoOut[i].nSmooth += pp->pInfoOut[i].nSmooth;
	}
    lstFreeTile(&pp->ilp->lst,&pp->tile->lstTile);
    pkdParticleWorkDone(pp->work);
    free(pp->pInfoOut);
    free(pp);
    return 0;
    }

int doneWorkSPHForces(void *vpp) {
    workPP *pp = vpp;
    int i;

    for(i=0; i<pp->work->nP; ++i) {
    pp->work->pInfoOut[i].thetaDot += pp->pInfoOut[i].thetaDot;
    pp->work->pInfoOut[i].a[0] += pp->pInfoOut[i].a[0];
	pp->work->pInfoOut[i].a[1] += pp->pInfoOut[i].a[1];
	pp->work->pInfoOut[i].a[2] += pp->pInfoOut[i].a[2];
    pp->work->pInfoOut[i].divv += pp->pInfoOut[i].divv;
    pp->work->pInfoOut[i].dtEst = fmin(pp->work->pInfoOut[i].dtEst,pp->pInfoOut[i].dtEst);
	}
    lstFreeTile(&pp->ilp->lst,&pp->tile->lstTile);
    pkdParticleWorkDone(pp->work);
    free(pp->pInfoOut);
    free(pp);
    return 0;
    }

static void queuePP( PKD pkd, workParticle *wp, ILP ilp, int bGravStep ) {
    ILPTILE tile;
    workPP *pp;
    ILP_LOOP(ilp,tile) {
#ifdef USE_CUDA
	if (CudaClientQueuePP(pkd->cudaClient,wp,tile,bGravStep)) continue;
#endif
	pp = malloc(sizeof(workPP));
	assert(pp!=NULL);
	pp->pInfoOut = malloc(sizeof(PINFOOUT) * wp->nP);
	assert(pp->pInfoOut!=NULL);
	pp->work = wp;
	pp->ilp = ilp;
	pp->tile = tile;
	pp->i = 0;
	tile->lstTile.nRefs++;
	wp->nRefs++;
	mdlAddWork(pkd->mdl,pp,NULL,NULL,CPUdoWorkPP,doneWorkPP);
	}
    }

static void queueDensity( PKD pkd, workParticle *wp, ILP ilp, int bGravStep ) {
    ILPTILE tile;
    workPP *pp;
    // initialize kernel mass deviation so the loop runs at least once
    int loopcount = 0;
    while(1) {// start while loop
    loopcount++;
    // Zero density and density derivative
    for( int i=0; i<wp->nP; i++ ) {
        wp->pInfoOut[i].rho = 0.0f;
        wp->pInfoOut[i].drhodfball = 0.0f;
        wp->pInfoOut[i].nden = 0.0f;
        wp->pInfoOut[i].dndendfball = 0.0f;
        wp->pInfoOut[i].nSmooth = 0.0f;
    }

    ILP_LOOP(ilp,tile) {
#ifdef USE_CUDA
	assert(0);
#endif
	pp = malloc(sizeof(workPP));
	assert(pp!=NULL);
	pp->pInfoOut = malloc(sizeof(PINFOOUT) * wp->nP);
	assert(pp->pInfoOut!=NULL);
	pp->work = wp;
	pp->ilp = ilp;
	pp->tile = tile;
	pp->i = 0;
	tile->lstTile.nRefs++;
	wp->nRefs++;
	mdlAddWork(pkd->mdl,pp,NULL,NULL,CPUdoWorkDensity,doneWorkDensity);
	}

    // we have to wait for all work to have finished here
    // all this only makes sense if we are not asynchronous
    assert(wp->nRefs == 1);

    float maxkerneldeviation = 0.0f;
    // calculate maximum kernel mass deviation
    for (int i=0; i<wp->nP; i++) {
        float kerneldeviation = 0.0f;
        if (wp->SPHoptions->useNumDen) {
            kerneldeviation = 4.0f/3.0f*M_PI*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoOut[i].nden - wp->SPHoptions->fKernelTarget;
        } else {
            kerneldeviation = 4.0f/3.0f*M_PI*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoOut[i].rho - wp->SPHoptions->fKernelTarget;
        }
        kerneldeviation = (kerneldeviation > 0) ? kerneldeviation : -kerneldeviation;
        maxkerneldeviation = (kerneldeviation > maxkerneldeviation) ? kerneldeviation : maxkerneldeviation;
    }

    /*
    ** decide if loop has to continue
    ** if true, calculate new fBall for all particles
    ** else, exit loop
    */
    if (maxkerneldeviation/wp->SPHoptions->fKernelTarget > 1e-6f) {
        // do another loop
        for (int i=0; i<wp->nP; i++) {
            float prefac = 4.0f/3.0f*M_PI;
            float fBall = wp->pInfoIn[i].fBall;
            float fx, dfdx;
            if (wp->SPHoptions->useNumDen) {
                fx = prefac * fBall * fBall * fBall * wp->pInfoOut[i].nden - wp->SPHoptions->fKernelTarget;
                dfdx = prefac * 3.0f * fBall * fBall * wp->pInfoOut[i].nden + prefac * fBall * fBall * fBall * wp->pInfoOut[i].dndendfball;
            } else {
                fx = prefac * fBall * fBall * fBall * wp->pInfoOut[i].rho - wp->SPHoptions->fKernelTarget;
                dfdx = prefac * 3.0f * fBall * fBall * wp->pInfoOut[i].rho + prefac * fBall * fBall * fBall * wp->pInfoOut[i].drhodfball;
            }
            float newfBall = wp->pInfoIn[i].fBall - fx / dfdx;
            if (newfBall < 0.5f * wp->pInfoIn[i].fBall) {
                wp->pInfoIn[i].fBall = 0.5f * wp->pInfoIn[i].fBall;
            } else {
                wp->pInfoIn[i].fBall = newfBall;
            }
        }
    } else {
        // finish
        // printf("loopcount = %d, deviation = %.6e\n",loopcount,maxkerneldeviation/wp->SPHoptions->fKernelTarget);
        break;
    }

    }// end while loop

    for( int i=0; i<wp->nP; i++ ) {
        // save the new fBall for each particle
        wp->pInfoOut[i].fBall = wp->pInfoIn[i].fBall;
        // wp->pInfoOut[i].fBall = wp->pInfoOut[i].nSmooth;
    }
    }

static void queueSPHForces( PKD pkd, workParticle *wp, ILP ilp, int bGravStep ) {
    ILPTILE tile;
    workPP *pp;
    ILP_LOOP(ilp,tile) {
#ifdef USE_CUDA
	assert(0);
#endif
	pp = malloc(sizeof(workPP));
	assert(pp!=NULL);
	pp->pInfoOut = malloc(sizeof(PINFOOUT) * wp->nP);
	assert(pp->pInfoOut!=NULL);
	pp->work = wp;
	pp->ilp = ilp;
	pp->tile = tile;
	pp->i = 0;
	tile->lstTile.nRefs++;
	wp->nRefs++;
	mdlAddWork(pkd->mdl,pp,NULL,NULL,CPUdoWorkSPHForces,doneWorkSPHForces);
	}
    }

int CPUdoWorkPC(void *vpc) {
    workPC *pc = vpc;
    workParticle *wp = pc->work;
    ILCTILE tile = pc->tile;
    ILC_BLK *blk = tile->blk;
    PINFOIN *pPart = &pc->work->pInfoIn[pc->i];
    PINFOOUT *pOut = &pc->pInfoOut[pc->i];
    int nBlocks = tile->lstTile.nBlocks;
    int nInLast = tile->lstTile.nInLast;

    pOut->a[0] = 0.0;
    pOut->a[1] = 0.0;
    pOut->a[2] = 0.0;
    pOut->fPot = 0.0;
    pOut->dirsum = 0.0;
    pOut->normsum = 0.0;
    pkdGravEvalPC(pPart,nBlocks,nInLast,blk,pOut);
    wp->dFlopSingleCPU += COST_FLOP_PC*(tile->lstTile.nBlocks*ILC_PART_PER_BLK  + tile->lstTile.nInLast);
    if ( ++pc->i == pc->work->nP ) return 0;
    else return 1;
    }

int doneWorkPC(void *vpc) {
    workPC *pc = vpc;
    int i;

    for(i=0; i<pc->work->nP; ++i) {
	pc->work->pInfoOut[i].a[0] += pc->pInfoOut[i].a[0];
	pc->work->pInfoOut[i].a[1] += pc->pInfoOut[i].a[1];
	pc->work->pInfoOut[i].a[2] += pc->pInfoOut[i].a[2];
	pc->work->pInfoOut[i].fPot += pc->pInfoOut[i].fPot;
	pc->work->pInfoOut[i].dirsum += pc->pInfoOut[i].dirsum;
	pc->work->pInfoOut[i].normsum += pc->pInfoOut[i].normsum;
	}

    lstFreeTile(&pc->ilc->lst,&pc->tile->lstTile);
    pkdParticleWorkDone(pc->work);
    free(pc->pInfoOut);
    free(pc);
    return 0;
    }

static void queuePC( PKD pkd,  workParticle *wp, ILC ilc, int bGravStep ) {
    ILCTILE tile;
    workPC *pc;

    ILC_LOOP(ilc,tile) {
#ifdef USE_CUDA
        if (CudaClientQueuePC(pkd->cudaClient,wp,tile,bGravStep)) continue;
#endif
	pc = malloc(sizeof(workPC));
	assert(pc!=NULL);
	pc->pInfoOut = malloc(sizeof(PINFOOUT) * wp->nP);
	assert(pc->pInfoOut!=NULL);
	pc->work = wp;
	pc->ilc = ilc;
	pc->tile = tile;
	pc->i = 0;
	tile->lstTile.nRefs++;
	wp->nRefs++;
	mdlAddWork(pkd->mdl,pc,NULL,NULL,CPUdoWorkPC,doneWorkPC);
	}
    }

static void queueEwald( PKD pkd, workParticle *wp ) {
    int i;
#ifdef USE_CUDA
    int nQueued = CudaClientQueueEwald(pkd->cudaClient,wp);
#else
    int nQueued = 0;
#endif
    ++wp->nRefs;
    for( i=nQueued; i<wp->nP; ++i) {
	PINFOIN *in = &wp->pInfoIn[i];
	PINFOOUT *out = &wp->pInfoOut[i];
	double r[3];
	r[0] = wp->c[0] + in->r[0];
	r[1] = wp->c[1] + in->r[1];
	r[2] = wp->c[2] + in->r[2];
	wp->dFlop += pkdParticleEwald(pkd,r,out->a,&out->fPot,&wp->dFlopSingleCPU,&wp->dFlopDoubleCPU);
	}
    pkdParticleWorkDone(wp);
    }

/*
** This version of grav.c does all the operations inline, including
** v_sqrt's and such.
** Returns nActive.
*/

int pkdGravInteract(PKD pkd,
    struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
    KDN *pBucket,LOCR *pLoc,ILP ilp,ILC ilc,
    float dirLsum,float normLsum,int bEwald,double *pdFlop,
    SMX smx,SMF *smf,int iRoot1,int iRoot2,SPHOptions *SPHoptions) {
    PARTICLE *p;
    KDN *pkdn = pBucket;
    double r[3], kdn_r[3];
    vel_t *v;
    double vx,vy,vz;
    float *ap;
    float d2,dir,dir2;
    float fMass,fSoft;
    float fx,fy,fz;
    double dx,dy,dz;
    float rhopmax,rhopmaxlocal;
    float summ;
    float fBall;
    ILPTILE tile;
    int i,n,nSoft,nActive;
    float fourh2;
    int nP;

    pkdNodeGetPos(pkd,pkdn,kdn_r);

    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;

    /* Collect the bucket particle information */
    workParticle *wp = malloc(sizeof(workParticle));
    assert(wp!=NULL);
    wp->nRefs = 1; /* I am using it currently */
    wp->ctx = pkd;
    wp->dFlop = 0.0;
    wp->dFlopSingleCPU = wp->dFlopSingleGPU = 0.0;
    wp->dFlopDoubleCPU = wp->dFlopDoubleGPU = 0.0;
    wp->nP = 0;

    nP = pkdn->pUpper - pkdn->pLower + 1;
    wp->pPart = malloc(sizeof(PARTICLE *) * nP); assert(wp->pPart != NULL);
    wp->iPart = malloc(sizeof(uint32_t) * nP); assert(wp->iPart != NULL);
    wp->pInfoIn = malloc(sizeof(PINFOIN) * nP); assert(wp->pInfoIn != NULL);
    wp->pInfoOut = malloc(sizeof(PINFOOUT) * nP); assert(wp->pInfoOut != NULL);
    wp->c[0] = ilp->cx; assert(wp->c[0] == ilc->cx);
    wp->c[1] = ilp->cy; assert(wp->c[1] == ilc->cy);
    wp->c[2] = ilp->cz; assert(wp->c[2] == ilc->cz);

    wp->dirLsum = dirLsum; // not used here...
    wp->normLsum = normLsum; // not used here...

    wp->ts = ts;
    wp->lc = lc;
    wp->kick = kick;

    wp->SPHoptions = SPHoptions;

    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	if (!pkdIsRungRange(p,ts->uRungLo,ts->uRungHi)) continue;
	pkdGetPos1(pkd,p,r);
	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	v = pkdVel(pkd,p);

	nP = wp->nP++;
	wp->pPart[nP] = p;
	wp->iPart[nP] = i;
	wp->pInfoIn[nP].r[0]  = (float)(r[0] - ilp->cx);
	wp->pInfoIn[nP].r[1]  = (float)(r[1] - ilp->cy);
	wp->pInfoIn[nP].r[2]  = (float)(r[2] - ilp->cz);
	if (pkd->oFieldOffset[oAcceleration]) {
	    ap = pkdAccel(pkd,p);
	    wp->pInfoIn[nP].a[0]  = ap[0];
	    wp->pInfoIn[nP].a[1]  = ap[1];
	    wp->pInfoIn[nP].a[2]  = ap[2];
	    ap[0] = ap[1] = ap[2] = 0.0;
	    }
	else {
	    wp->pInfoIn[nP].a[0]  = 0;
	    wp->pInfoIn[nP].a[1]  = 0;
	    wp->pInfoIn[nP].a[2]  = 0;
	    }

    //SPHFIELDS *pSph = pkdSph(pkd,p);
    float dtPredDrift = getDtPredDrift(kick,p->bMarked,ts->uRungLo,p->uRung);
    wp->pInfoIn[nP].fBall = pkdBall(pkd,p);
    wp->pInfoIn[nP].Omega = 1.0f; //pSph->Omega;
    wp->pInfoIn[nP].v[0] = v[0];
    wp->pInfoIn[nP].v[1] = v[1];
    wp->pInfoIn[nP].v[2] = v[2];
    wp->pInfoIn[nP].rho = pkdDensity(pkd,p);
    wp->pInfoIn[nP].P = 0.0f;                       /* should be calculated by the EOS, nyi */
    wp->pInfoIn[nP].c = 0.0f;                       /* should be calculated by the EOS, nyi */
    wp->pInfoIn[nP].species = pkdSpecies(pkd,p);

    wp->pInfoOut[nP].rho = 0.0f;
    wp->pInfoOut[nP].drhodfball = 0.0f;
    wp->pInfoOut[nP].nden = 0.0f;
    wp->pInfoOut[nP].dndendfball = 0.0f;
    wp->pInfoOut[nP].fBall = 0.0f;
    wp->pInfoOut[nP].nSmooth = 0.0f;
    wp->pInfoOut[nP].thetaDot = 0.0f;
    wp->pInfoOut[nP].divv = 0.0f;
    wp->pInfoOut[nP].dtEst = HUGE_VAL;

	wp->pInfoOut[nP].a[0] = 0.0f;
	wp->pInfoOut[nP].a[1] = 0.0f;
	wp->pInfoOut[nP].a[2] = 0.0f;
	wp->pInfoOut[nP].fPot = 0.0f;
	wp->pInfoOut[nP].dirsum = dirLsum;
	wp->pInfoOut[nP].normsum = normLsum;
	wp->pInfoOut[nP].rhopmax = 0.0f;

	/*
	** Calculate local density and kernel smoothing length for dynamical time-stepping
	*/
    if (SPHoptions->doGravity) {
	if (ts->bGravStep) {
	    /*
	    ** Calculate local density using smooth; this is fast because the particles are
	    ** likely to be cached already because they will be on the P-P list.
	    */
	    smf->pfDensity = &wp->pInfoIn[nP].fDensity;
	    fBall = smSmoothSingle(smx,smf,p,iRoot1,iRoot2);
	    wp->pInfoIn[nP].fSmooth2 = fBall * fBall;
	    }
	else {
	    /*
	    ** We are not using GravStep!
	    */
	    wp->pInfoIn[nP].fSmooth2 = 0.0;
	    }
	}
    }

    nActive += wp->nP;

    /*
    ** Evaluate the local expansion.
    */
    if (pLoc && SPHoptions->doGravity) {
	for( i=0; i<wp->nP; i++ ) {
	    momFloat ax,ay,az, dPot;
	    double *c = wp->c;
	    float *in = wp->pInfoIn[i].r;
	    //pkdGetPos1(wp->pPart[i]->r,r);
	    r[0] = c[0] + in[0];
	    r[1] = c[1] + in[1];
	    r[2] = c[2] + in[2];

	    /*
	    ** Evaluate local expansion.
	    */
	    dx = r[0] - kdn_r[0];
	    dy = r[1] - kdn_r[1];
	    dz = r[2] - kdn_r[2];
	    dPot = 0;
	    ax = 0;
	    ay = 0;
	    az = 0;

	    momEvalLocr(pLoc,dx,dy,dz,&dPot,&ax,&ay,&az);

	    wp->pInfoOut[i].fPot += dPot;
	    wp->pInfoOut[i].a[0] += ax;
	    wp->pInfoOut[i].a[1] += ay;
	    wp->pInfoOut[i].a[2] += az;
	    }
	}

    if (SPHoptions->doGravity) {
    /*
    ** Evaluate the P-C interactions
    */
    assert(ilc->cx==ilp->cx && ilc->cy==ilp->cy && ilc->cz==ilp->cz );
    queuePC( pkd,  wp, ilc, ts->bGravStep );

    /*
    ** Evaluate the P-P interactions
    */
    queuePP( pkd, wp, ilp, ts->bGravStep );
    }

    if (SPHoptions->doDensity) {
    /*
    ** Evaluate the Density on the P-P interactions
    */
    queueDensity( pkd, wp, ilp, ts->bGravStep );
    }

    if (SPHoptions->doSPHForces) {
    /*
    ** Evaluate the SPH forces on the P-P interactions
    */
    queueSPHForces( pkd, wp, ilp, ts->bGravStep );
    }

    if (SPHoptions->doGravity) {
    /*
    ** Calculate the Ewald correction for this particle, if it is required.
    */
    if (bEwald) {
	queueEwald( pkd,  wp );
	}
    }

#ifdef TIMESTEP_CRITICAL
    if (SPHoptions->doGravity) {
    for( i=0; i<wp->nP; i++ ) {
	double *c = wp->c;
	float *in = wp->pInfoIn[i].r;
	r[0] = c[0] + in[0];
	r[1] = c[1] + in[1];
	r[2] = c[2] + in[2];
        //p = wp->pPart[i];
	//pkdGetPos1(p->r,r);
        v = pkdVel(pkd,p);
        /*
        ** Set value for time-step, note that we have the current ewald acceleration
        ** in this now as well!
        */
        if (ts->bGravStep && wp->ts->iTimeStepCrit == 1) {
	    /*
	    ** GravStep if iTimeStepCrit =
	    ** 0: Mean field regime for dynamical time (normal/standard setting)
	    ** 1: Gravitational scattering regime for dynamical time with eccentricity correction
	    */
	    rhopmax = 0.0;
	    ILP_LOOP(ilp,tile) {
		int blk,prt;
		for( blk=0; blk<=tile->lstTile.nBlocks; ++blk ) {
		    n = (blk==tile->lstTile.nBlocks ? tile->lstTile.nInLast : ilp->lst.nPerBlock);
		    for (prt=0; prt<n; ++prt) {
			if (p->iOrder < wp->ts->nPartColl || tile->xtr[blk].iOrder.i[prt] < wp->ts->nPartColl) {
			    fx = r[0] - ilp->cx + tile->blk[blk].dx.f[prt];
			    fy = r[1] - ilp->cy + tile->blk[blk].dy.f[prt];
			    fz = r[2] - ilp->cz + tile->blk[blk].dz.f[prt];
			    d2 = fx*fx + fy*fy + fz*fz;
			    if (p->iOrder == tile->xtr[blk].iOrder.i[prt]) continue;
			    fourh2 = softmassweight(fMass,4*fSoft*fSoft,
				tile->blk[blk].m.f[prt],tile->blk[blk].fourh2.f[prt]);
			    if (d2 > fourh2) {
				SQRT1(d2,dir);
				dir2 = dir*dir*dir;
				}
			    else {
				/*
				** This uses the Dehnen K1 kernel function now, it's fast!
				*/
				SQRT1(fourh2,dir);
				dir2 = dir*dir;
				d2 *= dir2;
				dir2 *= dir;
				d2 = 1 - d2;
				dir *= 1.0f + d2*(0.5f + d2*(3.0f/8.0f + d2*(45.0f/32.0f)));
				dir2 *= 1.0f + d2*(1.5f + d2*(135.0f/16.0f));
				}
			    summ = fMass+tile->blk[blk].m.f[prt];
			    rhopmaxlocal = summ*dir2;
			    vx = v[0] - tile->xtr[blk].vx.d[prt];
			    vy = v[1] - tile->xtr[blk].vy.d[prt];
			    vz = v[2] - tile->xtr[blk].vz.d[prt];
			    rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,
				fx,fy,fz,vx,vy,vz,wp->ts->dEccFacMax);
			    rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;
			    }
			}
		    }
		}
	    wp->pInfoOut[i].rhopmax = rhopmax;
	    }
        } /* end of i-loop cells & particles */
    }
#endif

    pkdParticleWorkDone(wp);

    *pdFlop += nActive*(ilpCount(pkd->ilp)*COST_FLOP_PP + ilcCount(pkd->ilc)*COST_FLOP_PC) + nSoft*15;
    return(nActive);
    }

#ifdef TIMESTEP_CRITICAL
/*
** Gravitational scattering regime (iTimeStepCrit=1)
*/
double pkdRho1(double rhopmaxlocal, double summ, double dir, double x, double y, double z, double vx, double vy, double vz, double EccFacMax) {

    double Etot, L2, ecc, eccfac, v2;
    /*
    ** Etot and L are normalized by the reduced mass 
    */
    v2 = vx*vx + vy*vy + vz*vz;
    Etot = 0.5*v2 - summ*dir;
    L2 = (y*vz - z*vy)*(y*vz - z*vy) + (z*vx - x*vz)*(z*vx - x*vz) + (x*vy - y*vx)*(x*vy - y*vx);
    ecc = 1+2*Etot*L2/(summ*summ);
    ecc = (ecc <= 0)?0:sqrt(ecc);
    eccfac = (1 + 2*ecc)/fabs(1-ecc);
    eccfac = (eccfac > EccFacMax)?EccFacMax:eccfac;
    if (eccfac > 1.0) rhopmaxlocal *= eccfac;
    return rhopmaxlocal;
    }
#endif
