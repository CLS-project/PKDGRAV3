#ifdef HAVE_CONFIG_H
#include "config.h"
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

#ifdef USE_SIMD
static const struct CONSTS {
    vfloat minSoftening;
    vfloat zero;
    vfloat onequarter;
    vfloat onethird;
    vfloat half;
    vfloat one;
    vfloat threehalves;
    vfloat three;
    vfloat four;
    vfloat five;
    vfloat seven;
    vfloat nine;
    vfloat R3_8;
    vfloat R45_32;
    vfloat R135_16;
    } consts = {
        {SIMD_CONST(1e-18f)},
        {SIMD_CONST(0.0f)},
	{SIMD_CONST(0.25f)},
	{SIMD_CONST(1.0f/3.0f)},
	{SIMD_CONST(0.5f)},
	{SIMD_CONST(1.0f)},
	{SIMD_CONST(1.5f)},
	{SIMD_CONST(3.0f)},
	{SIMD_CONST(4.0f)},
	{SIMD_CONST(5.0f)},
	{SIMD_CONST(7.0f)},
	{SIMD_CONST(9.0f)},
	{SIMD_CONST(3.0f/8.0f)},
	{SIMD_CONST(45.0f/32.0f)},
	{SIMD_CONST(135.0f/16.0f)},
    };
#endif

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
void pkdParticleWorkDone(workParticle *work) {
    PKD pkd = work->ctx;
    int i;
    PARTICLE *p;
    double r[3];
    vel_t *v,v2;
    float *a;
    float fx, fy, fz, m;
    float maga, dT, dtGrav;
    unsigned char uNewRung;

    if ( --work->nRefs == 0 ) {
	float fiDelta = 1.0/pkd->param.dDelta;
	float fEta = pkd->param.dEta;
	float fiAccFac = 1.0 / work->dAccFac;
	pkd->dFlop += work->dFlop;
	pkd->dFlopSingleCPU += work->dFlopSingleCPU;
	pkd->dFlopDoubleCPU += work->dFlopDoubleCPU;
	pkd->dFlopSingleGPU += work->dFlopSingleGPU;
	pkd->dFlopDoubleGPU += work->dFlopDoubleGPU;
	for( i=0; i<work->nP; i++ ) {
	    p = work->pPart[i];
	    pkdGetPos1(pkd,p,r);
	    m = pkdMass(pkd,p);
	    if (pkd->oAcceleration) {
		a = pkdAccel(pkd,p);
		a[0] = work->pInfoOut[i].a[0];
		a[1] = work->pInfoOut[i].a[1];
		a[2] = work->pInfoOut[i].a[2];
		}
	    if (pkd->oPotential) {
		float *pPot = pkdPot(pkd,p);
		*pPot = work->pInfoOut[i].fPot;
		}
	    a = work->pInfoOut[i].a;
	    pkd->dEnergyU += 0.5 * m * work->pInfoOut[i].fPot;
	    pkd->dEnergyW += m*(r[0]*a[0] + r[1]*a[1] + r[2]*a[2]);
	    pkd->dEnergyF[0] += m*a[0];
	    pkd->dEnergyF[1] += m*a[1];
	    pkd->dEnergyF[2] += m*a[2];

	    if (work->bGravStep) {
		float dirsum = work->pInfoOut[i].dirsum;
		float normsum = work->pInfoOut[i].normsum;

		/*
		** If this is the first time through, the accelerations will have 
		** all been zero resulting in zero for normsum (and nan for dtGrav).
		** We repeat this process again, so dtGrav will be correct.
		*/
		if (normsum > 0.0) {
		    /*
		    ** Use new acceleration here!
		    */
		    fx = work->pInfoOut[i].a[0];
		    fy = work->pInfoOut[i].a[1];
		    fz = work->pInfoOut[i].a[2];
		    maga = fx*fx + fy*fy + fz*fz;
		    if (maga>0.0f) maga = asqrtf(maga);
		    dtGrav = maga*dirsum/normsum;
		    }
		else dtGrav = 0.0;
		dtGrav += pkd->param.dPreFacRhoLoc*work->pInfoIn[i].fDensity;
		dtGrav = (work->pInfoOut[i].rhopmax > dtGrav?work->pInfoOut[i].rhopmax:dtGrav);
		if (dtGrav > 0.0) {
		    dT = fEta * rsqrtf(dtGrav*work->dRhoFac);
		    uNewRung = pkdDtToRungInverse(dT,fiDelta,pkd->param.iMaxRung-1);
		    }
		else uNewRung = 0; /* Assumes current uNewRung is outdated -- not ideal */
		} /* end of work->bGravStep */
	    else {
		/*
		** We are doing eps/a timestepping.
		*/
		fx = work->pInfoOut[i].a[0];
		fy = work->pInfoOut[i].a[1];
		fz = work->pInfoOut[i].a[2];
		maga = fx*fx + fy*fy + fz*fz;
		if (maga > 0) {
		    float imaga = rsqrtf(maga) * fiAccFac;
		    dT = fEta*asqrtf(pkdSoft(pkd,p)*imaga);
		    uNewRung = pkdDtToRungInverse(dT,fiDelta,pkd->param.iMaxRung-1);
		    }
		else uNewRung = 0;
		}
	    /*
	    ** Here we must make sure we do not try to take a larger opening
	    ** timestep than the current active particles involved in the
	    ** gravity calculation!
	    */
	    if (uNewRung < work->uRungLo) uNewRung = work->uRungLo;
	    else if (uNewRung > work->uRungHi) uNewRung = work->uRungHi; 
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
	    if (pkd->oVelocity) {
		v = pkdVel(pkd,p);
		if (work->bKickClose) {
		    v[0] += work->dtClose[p->uRung]*work->pInfoOut[i].a[0];
		    v[1] += work->dtClose[p->uRung]*work->pInfoOut[i].a[1];
		    v[2] += work->dtClose[p->uRung]*work->pInfoOut[i].a[2];
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

		if (work->bKickOpen) {
		    p->uRung = uNewRung;
		    ++pkd->nRung[p->uRung];
		    v[0] += work->dtOpen[p->uRung]*work->pInfoOut[i].a[0];
		    v[1] += work->dtOpen[p->uRung]*work->pInfoOut[i].a[1];
		    v[2] += work->dtOpen[p->uRung]*work->pInfoOut[i].a[2];		    
		    }
		}
	    }
	free(work->pPart);
	free(work->pInfoIn);
	free(work->pInfoOut);
	free(work);
	}
    }

void pkdGravEvalPP(PINFOIN *pPart, int nBlocks, int nInLast, ILP_BLK *blk,  PINFOOUT *pOut ) {
    int nLeft;
    int j;
#if defined(USE_SIMD_PP)
    v_sf t1, t2, t3, pd2;
    v_sf pax, pay, paz, pfx, pfy, pfz, pdx, pdy, pdz;
    v_sf piax, piay, piaz;
    v_sf ppot /*,pmass,p4soft2*/;
    v_sf padotai,pimaga,psmooth2,pirsum,pnorms;
#else
    float d2,dx,dy,dz,fourh2,dir,dir2,adotai;
    int n;
#endif
    float ax,ay,az,fPot,dirsum,normsum;
    float tax,tay,taz;
    float dimaga;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    /*float fMass = pPart->.fMass;*/
    /*float fSoft = pPart->.fSoft;*/
    float fsmooth2 = pPart->fSmooth2;
    float *a =pPart->a;
    int nSoft, nIntr;

    dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0) {
	dimaga = 1.0f/sqrtf(dimaga);
	}

#ifdef USE_SIMD_PP
    pimaga = SIMD_SPLAT(dimaga);

    /*
    ** This is a little trick to speed up the calculation. By setting
    ** unused entries in the list to have a zero mass, the resulting
    ** forces are zero. Setting the distance to a large value avoids
    ** softening the non-existent forces which is slightly faster.
    */
    for( j = nInLast; j&SIMD_MASK; j++) {
	blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
	blk[nBlocks].m.f[j] = 0.0f;
	blk[nBlocks].fourh2.f[j] = 1e-18f;
	}

    /*
    ** The list sets mass to zero for unused entries which results
    ** in zero forces. Be careful if that is changed.
    */
    pax = consts.zero.p;
    pay = consts.zero.p;
    paz = consts.zero.p;
    ppot= consts.zero.p;
    pirsum = consts.zero.p;
    pnorms = consts.zero.p;

    piax    = SIMD_SPLAT(a[0]);
    piay    = SIMD_SPLAT(a[1]);
    piaz    = SIMD_SPLAT(a[2]);
    pfx     = SIMD_SPLAT(fx);
    pfy     = SIMD_SPLAT(fy);
    pfz     = SIMD_SPLAT(fz);
    /*pmass   = SIMD_SPLAT(fMass);*/
    /*p4soft2 = SIMD_SPLAT(4.0*fSoft*fSoft);*/
    psmooth2= SIMD_SPLAT(fsmooth2);

    nSoft = 0;
    nIntr = nBlocks * ILP_PART_PER_BLK + nInLast;
    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = ((nLeft ? ILP_PART_PER_BLK : nInLast) + SIMD_MASK) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    v_sf pfourh2, td2, pir, pir2;
	    v_bool vcmp;
	    int msk;

	    pdx = SIMD_ADD(blk->dx.p[j],pfx);
	    pdy = SIMD_ADD(blk->dy.p[j],pfy);
	    pdz = SIMD_ADD(blk->dz.p[j],pfz);
	    pd2 = SIMD_MADD(pdz,pdz,SIMD_MADD(pdy,pdy,SIMD_MUL(pdx,pdx)));

	    pfourh2 = blk->fourh2.p[j];
	    pfourh2 = SIMD_MAX(pfourh2,consts.minSoftening.p); /* There is always a self interaction */
	    vcmp = SIMD_CMP_LT(pd2,pfourh2);
	    td2 = SIMD_MAX(pd2,pfourh2);
	    msk = SIMD_ALL_ZERO(vcmp);  /* zero means nothing is softened - optimization */

	    td2 = SIMD_MAX(consts.minSoftening.p,td2);
	    pir = SIMD_RSQRT_EXACT(td2);
	    pir2 = SIMD_MUL(pir,pir);
	    td2 = SIMD_MUL(pir2,pd2); /* for SOFTENED */
	    pir2 = SIMD_MUL(pir2,pir);

	    /* pir and pir2 are valid now for both softened and unsoftened particles */
	    /* Now we apply the fix to softened particles only */
	    if (msk) {
		++nSoft;
		td2 = SIMD_SUB(consts.one.p,td2);
		td2 = SIMD_AND(vcmp,td2);
		t1 = SIMD_MADD(consts.R45_32.p, td2, consts.R3_8.p);
		t1 = SIMD_MADD(t1, td2, consts.half.p);
		t1 = SIMD_MADD(t1, td2, consts.one.p);
		t2 = SIMD_MADD(consts.R135_16.p, td2, consts.threehalves.p);
		t2 = SIMD_MADD(t2, td2, consts.one.p);
		pir = SIMD_MUL(pir,t1);
		pir2 = SIMD_MUL(pir2,t2);
		}
	    pir2 = SIMD_MUL(pir2,blk->m.p[j]);

	    t1 = SIMD_NMADD(pdx,pir2,consts.zero.p);
	    t2 = SIMD_NMADD(pdy,pir2,consts.zero.p);
	    t3 = SIMD_NMADD(pdz,pir2,consts.zero.p);
	    ppot = SIMD_NMADD(blk->m.p[j],pir,ppot);

	    /* Time stepping criteria stuff */
	    padotai = SIMD_MADD(piaz,t3,SIMD_MADD(piay,t2,SIMD_MUL(piax,t1)));
	    vcmp = SIMD_AND(SIMD_CMP_GT(padotai,consts.zero.p),SIMD_CMP_GE(pd2,psmooth2));
	    padotai= SIMD_AND(padotai,vcmp);
	    padotai= SIMD_MUL(padotai,pimaga);
	    td2 = SIMD_MUL(padotai,padotai);
	    pirsum = SIMD_MADD(pir,td2,pirsum);
	    pnorms = SIMD_ADD(pnorms,td2);

	    pax = SIMD_ADD(pax,t1);
	    pay = SIMD_ADD(pay,t2);
	    paz = SIMD_ADD(paz,t3);
	    }
	}
    ax = SIMD_HADD(pax);
    ay = SIMD_HADD(pay);
    az = SIMD_HADD(paz);
    fPot = SIMD_HADD(ppot);
    dirsum = SIMD_HADD(pirsum);
    normsum = SIMD_HADD(pnorms);
#else
    /*
    ** DO NOT MODIFY THE CODE BELOW UP TO THE #endif!
    ** This code MUST match the SIMD code above.
    */
    ax = 0.0;
    ay = 0.0;
    az = 0.0;
    fPot= 0.0;
    dirsum = 0.0;
    normsum = 0.0;

    nSoft = 0;
    nIntr = nBlocks * ILP_PART_PER_BLK + nInLast;
    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = (nLeft ? ILP_PART_PER_BLK : nInLast);
	for (j=0; j<n; ++j) {
	    dx = fx + blk->dx.f[j];
	    dy = fy + blk->dy.f[j];
	    dz = fz + blk->dz.f[j];
	    fourh2 = blk->fourh2.f[j];
	    d2 = dx*dx + dy*dy + dz*dz;
	    if (d2==0.0) continue; //dir2 = 0.0;
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {
		++nSoft;
		/*
		** This uses the Dehnen K1 kernel function now, it's fast!
		*/
		SQRT1(fourh2,dir);
		dir2 = dir*dir;
		tax = d2 * dir2;
		dir2 *= dir;
		tax = 1 - tax;
		dir *= 1.0 + tax*(0.5 + tax*(3.0/8.0 + tax*(45.0/32.0)));
		dir2 *= 1.0 + tax*(1.5 + tax*(135.0/16.0));
		}
	    dir2 *= -blk->m.f[j];
	    tax = dx*dir2;
	    tay = dy*dir2;
	    taz = dz*dir2;
	    fPot -= blk->m.f[j]*dir;
	    /*
	    ** Calculations for determining the timestep.
	    */
	    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
	    if (adotai > 0 && d2 >= fsmooth2) {
		adotai *= dimaga;
		dirsum += dir*adotai*adotai;
		normsum += adotai*adotai;
		}
	    ax += tax;
	    ay += tay;
	    az += taz;
	    }
	}
#endif
    pOut->a[0] += ax;
    pOut->a[1] += ay;
    pOut->a[2] += az;
    pOut->fPot += fPot;
    pOut->dirsum += dirsum;
    pOut->normsum += normsum;
    }


int CPUdoWorkPP(void *vpp) {
    workPP *pp = vpp;
    workParticle *work = pp->work;
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
    work->dFlopSingleCPU += COST_FLOP_PP*(tile->lstTile.nBlocks*ILP_PART_PER_BLK  + tile->lstTile.nInLast);
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

static void queuePP( PKD pkd, workParticle *work, ILP ilp, int bGravStep ) {
    ILPTILE tile;
    workPP *pp;

    ILP_LOOP(ilp,tile) {
#ifdef USE_CUDA
	if (CUDA_queuePP(pkd->mdl->cudaCtx,work,tile,bGravStep)) continue;
#endif
	pp = malloc(sizeof(workPP));
	assert(pp!=NULL);
	pp->pInfoOut = malloc(sizeof(PINFOOUT) * work->nP);
	assert(pp->pInfoOut!=NULL);
	pp->work = work;
	pp->ilp = ilp;
	pp->tile = tile;
	pp->i = 0;
	tile->lstTile.nRefs++;
	work->nRefs++;
	mdlAddWork(pkd->mdl,pp,NULL,NULL,CPUdoWorkPP,doneWorkPP);
	}
    }

void pkdGravEvalPC(PINFOIN *pPart, int nBlocks, int nInLast, ILC_BLK *blk,  PINFOOUT *pOut ) {

#if defined(USE_SIMD_PC)
    v_sf u,g0,g1,g2,g3,g4;
    v_sf x,y,z;
    v_sf tx,ty,tz;
    v_sf xx,xy,xz,yy,yz,zz;
    v_sf xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    v_sf t1, t2, t3;
    v_sf pax, pay, paz;
    v_sf pdx, pdy, pdz;
    v_sf pfx, pfy, pfz;
    v_sf piax, piay, piaz;
    v_sf ppot /*,pmass,p4soft2*/;
    v_sf padotai,pimaga,pirsum,pnorms;
#else
    const float onethird = 1.0f/3.0f;
    float u,g0,g1,g2,g3,g4;
    float x,y,z;
    float adotai;
    float tx,ty,tz;
    float xx,xy,xz,yy,yz,zz;
    float xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    float dir;
#endif
    float dimaga, tax, tay, taz;
    int j, n, nLeft;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    float *a =pPart->a;
    float ax,ay,az,fPot,dirsum,normsum;

    dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0.0) {
        dimaga = 1.0f / sqrtf(dimaga);
        }

#if defined(USE_SIMD_PC)
    pax = consts.zero.p;
    pay = consts.zero.p;
    paz = consts.zero.p;
    ppot= consts.zero.p;
    pirsum = consts.zero.p;
    pnorms = consts.zero.p;
    pfx     = SIMD_SPLAT(fx);
    pfy     = SIMD_SPLAT(fy);
    pfz     = SIMD_SPLAT(fz);
    piax    = SIMD_SPLAT(a[0]);
    piay    = SIMD_SPLAT(a[1]);
    piaz    = SIMD_SPLAT(a[2]);
    /*pmass   = SIMD_SPLAT(fMass);*/
    /*p4soft2 = SIMD_SPLAT(4.0*fSoft*fSoft);*/
    pimaga  = SIMD_SPLAT(dimaga);

    /* Pad the last value if necessary */
    for( j = nInLast; j&SIMD_MASK; j++) {
	blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
	blk[nBlocks].m.f[j] = 0.0f;
	blk[nBlocks].u.f[j] = 0.0f;
	}

    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = ((nLeft ? ILC_PART_PER_BLK : nInLast) + SIMD_MASK) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    v_sf pir, pd2;
	    v_bool vcmp;

	    pdx = SIMD_ADD(blk->dx.p[j],pfx);
	    pdy = SIMD_ADD(blk->dy.p[j],pfy);
	    pdz = SIMD_ADD(blk->dz.p[j],pfz);
	    pir = SIMD_RSQRT_EXACT(SIMD_MADD(pdz,pdz,SIMD_MADD(pdy,pdy,SIMD_MUL(pdx,pdx))));
	    u = SIMD_MUL(blk->u.p[j],pir);
	    g0 = pir;
	    g1 = SIMD_MUL(g0,u);
	    g2 = SIMD_MUL(consts.three.p,SIMD_MUL(g1,u));
	    g3 = SIMD_MUL(consts.five.p,SIMD_MUL(g2,u));
	    g4 = SIMD_MUL(consts.seven.p,SIMD_MUL(g3,u));
	    /*
	    ** Calculate the funky distance terms.
	    */
	    x = SIMD_MUL(pdx,pir);
	    y = SIMD_MUL(pdy,pir);
	    z = SIMD_MUL(pdz,pir);
	    xx = SIMD_MUL(consts.half.p,SIMD_MUL(x,x));
	    xy = SIMD_MUL(x,y);
	    xz = SIMD_MUL(x,z);
	    yy = SIMD_MUL(consts.half.p,SIMD_MUL(y,y));
	    yz = SIMD_MUL(y,z);
	    zz = SIMD_MUL(consts.half.p,SIMD_MUL(z,z));
	    xxx = SIMD_MSUB(consts.onethird.p,xx,zz);
	    xxx = SIMD_MUL(x,xxx);
	    xxz = SIMD_NMADD(consts.onethird.p,zz,xx);
	    xxz = SIMD_MUL(z,xxz);
	    yyy = SIMD_MSUB(consts.onethird.p,yy,zz);
	    yyy = SIMD_MUL(y,yyy);
	    yyz = SIMD_NMADD(consts.onethird.p,zz,yy);
	    yyz = SIMD_MUL(z,yyz);
	    xx = SIMD_SUB(xx,zz);
	    yy = SIMD_SUB(yy,zz);
	    xxy = SIMD_MUL(y,xx);
	    xyy = SIMD_MUL(x,yy);
	    xyz = SIMD_MUL(xy,z);
	    /*
	    ** Now calculate the interaction up to Hexadecapole order.
	    */
	    tx = SIMD_MUL(blk->xxxx.p[j],xxx);
	    tx = SIMD_MADD(blk->xyyy.p[j],yyy,tx);
	    tx = SIMD_MADD(blk->xxxy.p[j],xxy,tx);
	    tx = SIMD_MADD(blk->xxxz.p[j],xxz,tx);
	    tx = SIMD_MADD(blk->xxyy.p[j],xyy,tx);
	    tx = SIMD_MADD(blk->xxyz.p[j],xyz,tx);
	    tx = SIMD_MADD(blk->xyyz.p[j],yyz,tx);
	    tx = SIMD_MUL(tx,g4);
	    ty = SIMD_MUL(blk->xyyy.p[j],xyy);
	    ty = SIMD_MADD(blk->xxxy.p[j],xxx,ty);
	    ty = SIMD_MADD(blk->yyyy.p[j],yyy,ty);
	    ty = SIMD_MADD(blk->yyyz.p[j],yyz,ty);
	    ty = SIMD_MADD(blk->xxyy.p[j],xxy,ty);
	    ty = SIMD_MADD(blk->xxyz.p[j],xxz,ty);
	    ty = SIMD_MADD(blk->xyyz.p[j],xyz,ty);
	    ty = SIMD_MUL(ty,g4);
	    tz = SIMD_MUL(blk->xxxz.p[j],xxx);
	    tz = SIMD_NMADD(blk->xxxx.p[j],xxz,tz);
	    tz = SIMD_NMADD(SIMD_ADD(blk->xyyy.p[j],blk->xxxy.p[j]),xyz,tz);
	    tz = SIMD_NMADD(blk->yyyy.p[j],yyz,tz);
	    tz = SIMD_MADD(blk->yyyz.p[j],yyy,tz);
	    tz = SIMD_NMADD(blk->xxyy.p[j],SIMD_ADD(xxz,yyz),tz);
	    tz = SIMD_MADD(blk->xxyz.p[j],xxy,tz);
	    tz = SIMD_MADD(blk->xyyz.p[j],xyy,tz);
	    tz = SIMD_MUL(tz,g4);
	    g4 = SIMD_MUL(consts.onequarter.p,SIMD_MADD(tx,x,SIMD_MADD(ty,y,SIMD_MUL(tz,z))));
	    xxx = SIMD_MUL(g3,SIMD_MADD(blk->xxx.p[j],xx,SIMD_MADD(blk->xyy.p[j],yy,
			SIMD_MADD(blk->xxy.p[j],xy,SIMD_MADD(blk->xxz.p[j],xz,SIMD_MUL(blk->xyz.p[j],yz))))));
	    xxy = SIMD_MUL(g3,SIMD_MADD(blk->xyy.p[j],xy,SIMD_MADD(blk->xxy.p[j],xx,SIMD_MADD(blk->yyy.p[j],yy,
			    SIMD_MADD(blk->yyz.p[j],yz,SIMD_MUL(blk->xyz.p[j],xz))))));
	    xxz = SIMD_MUL(g3,SIMD_NMADD(SIMD_ADD(blk->xxx.p[j],blk->xyy.p[j]),xz,
		    SIMD_NMADD(SIMD_ADD(blk->xxy.p[j],blk->yyy.p[j]),yz,
			SIMD_MADD(blk->xxz.p[j],xx,SIMD_MADD(blk->yyz.p[j],yy,SIMD_MUL(blk->xyz.p[j],xy))))));
	    g3 = SIMD_MUL(consts.onethird.p,SIMD_MADD(xxx,x,SIMD_MADD(xxy,y,SIMD_MUL(xxz,z))));
	    xx = SIMD_MUL(g2,SIMD_MADD(blk->xx.p[j],x,SIMD_MADD(blk->xy.p[j],y,SIMD_MUL(blk->xz.p[j],z))));
	    xy = SIMD_MUL(g2,SIMD_MADD(blk->yy.p[j],y,SIMD_MADD(blk->xy.p[j],x,SIMD_MUL(blk->yz.p[j],z))));
	    xz = SIMD_MUL(g2,SIMD_NMADD(SIMD_ADD(blk->xx.p[j],blk->yy.p[j]),z,SIMD_MADD(blk->xz.p[j],x,SIMD_MUL(blk->yz.p[j],y))));
	    g2 = SIMD_MUL(consts.half.p,SIMD_MADD(xx,x,SIMD_MADD(xy,y,SIMD_MUL(xz,z))));
	    g0 = SIMD_MUL(g0,blk->m.p[j]);
	    ppot = SIMD_SUB(ppot,SIMD_ADD(SIMD_ADD(g0,g2),SIMD_ADD(g3,g4)));
	    g0 = SIMD_MADD(consts.five.p,g2,SIMD_MADD(consts.seven.p,g3,SIMD_MADD(consts.nine.p,g4,g0)));
#ifdef USE_DIAPOLE
	    yy = SIMD_MUL(g1,blk->x.p[j]);
	    yz = SIMD_MUL(g1,blk->y.p[j]);
	    zz = SIMD_MUL(g1,blk->z.p[j]);
	    g1 = SIMD_MADD(yy,x,SIMD_MADD(yz,y,SIMD_MUL(zz,z)));
	    ppot = SIMD_SUB(ppot,g1);
	    g0 = SIMD_MADD(consts.three.p,g1,g0);
#else
	    yy = consts.zero.p;
	    yz = consts.zero.p;
	    zz = consts.zero.p;
#endif
	    t1 = SIMD_MUL(pir,SIMD_NMADD(x,g0,SIMD_ADD(SIMD_ADD(yy,xx),SIMD_ADD(xxx,tx))));
	    t2 = SIMD_MUL(pir,SIMD_NMADD(y,g0,SIMD_ADD(SIMD_ADD(yz,xy),SIMD_ADD(xxy,ty))));
	    t3 = SIMD_MUL(pir,SIMD_NMADD(z,g0,SIMD_ADD(SIMD_ADD(zz,xz),SIMD_ADD(xxz,tz))));

	    /* Time stepping criteria stuff */
	    padotai = SIMD_MADD(piaz,t3,SIMD_MADD(piay,t2,SIMD_MUL(piax,t1)));
	    vcmp = SIMD_CMP_GT(padotai,consts.zero.p);
	    padotai= SIMD_AND(padotai,vcmp);
	    padotai= SIMD_MUL(padotai,pimaga);
	    pd2 = SIMD_MUL(padotai,padotai);
	    pirsum = SIMD_MADD(pir,pd2,pirsum);
	    pnorms = SIMD_ADD(pnorms,pd2);

	    pax = SIMD_ADD(pax,t1);
	    pay = SIMD_ADD(pay,t2);
	    paz = SIMD_ADD(paz,t3);
	    }
	}
    ax = SIMD_HADD(pax);
    ay = SIMD_HADD(pay);
    az = SIMD_HADD(paz);
    fPot = SIMD_HADD(ppot);
    dirsum = SIMD_HADD(pirsum);
    normsum = SIMD_HADD(pnorms);
#else
/*
** DO NOT MODIFY THE CODE BELOW UP TO THE #endif!
** This code MUST match the SIMD code above.
*/

    ax = 0.0;
    ay = 0.0;
    az = 0.0;
    fPot= 0.0;
    dirsum = 0.0;
    normsum = 0.0;

    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = (nLeft ? ILC_PART_PER_BLK : nInLast);
	for (j=0; j<n; ++j) {
	    float dx = blk->dx.f[j] + fx;
	    float dy = blk->dy.f[j] + fy;
	    float dz = blk->dz.f[j] + fz;
	    float d2 = dx*dx + dy*dy + dz*dz;
	    SQRT1(d2,dir);
	    u = blk->u.f[j]*dir;
	    g0 = dir;
	    g1 = g0*u;
	    g2 = 3*g1*u;
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
	    tx = g4*(blk->xxxx.f[j]*xxx + blk->xyyy.f[j]*yyy + blk->xxxy.f[j]*xxy + blk->xxxz.f[j]*xxz + blk->xxyy.f[j]*xyy + blk->xxyz.f[j]*xyz + blk->xyyz.f[j]*yyz);
	    ty = g4*(blk->xyyy.f[j]*xyy + blk->xxxy.f[j]*xxx + blk->yyyy.f[j]*yyy + blk->yyyz.f[j]*yyz + blk->xxyy.f[j]*xxy + blk->xxyz.f[j]*xxz + blk->xyyz.f[j]*xyz);
	    tz = g4*(-blk->xxxx.f[j]*xxz - (blk->xyyy.f[j] + blk->xxxy.f[j])*xyz - blk->yyyy.f[j]*yyz + blk->xxxz.f[j]*xxx + blk->yyyz.f[j]*yyy - blk->xxyy.f[j]*(xxz + yyz) + blk->xxyz.f[j]*xxy + blk->xyyz.f[j]*xyy);
	    g4 = 0.25*(tx*x + ty*y + tz*z);
	    xxx = g3*(blk->xxx.f[j]*xx + blk->xyy.f[j]*yy + blk->xxy.f[j]*xy + blk->xxz.f[j]*xz + blk->xyz.f[j]*yz);
	    xxy = g3*(blk->xyy.f[j]*xy + blk->xxy.f[j]*xx + blk->yyy.f[j]*yy + blk->yyz.f[j]*yz + blk->xyz.f[j]*xz);
	    xxz = g3*(-(blk->xxx.f[j] + blk->xyy.f[j])*xz - (blk->xxy.f[j] + blk->yyy.f[j])*yz + blk->xxz.f[j]*xx + blk->yyz.f[j]*yy + blk->xyz.f[j]*xy);
	    g3 = onethird*(xxx*x + xxy*y + xxz*z);
	    xx = g2*(blk->xx.f[j]*x + blk->xy.f[j]*y + blk->xz.f[j]*z);
	    xy = g2*(blk->yy.f[j]*y + blk->xy.f[j]*x + blk->yz.f[j]*z);
	    xz = g2*(-(blk->xx.f[j] + blk->yy.f[j])*z + blk->xz.f[j]*x + blk->yz.f[j]*y);
	    g2 = 0.5*(xx*x + xy*y + xz*z);
	    g0 *= blk->m.f[j];
	    fPot -= g0 + g2 + g3 + g4;
	    g0 += 5*g2 + 7*g3 + 9*g4;
#ifdef USE_DIAPOLE
	    yy = g1*blk->x.f[j];
	    yz = g1*blk->y.f[j];
	    zz = g1*blk->z.f[j];
	    g1 = (yy*x + yz*y + zz*z);
	    fPot -= g1;
	    g0 += 3*g1;
#else
	    yy = 0.0f;
	    yz = 0.0f;
	    zz = 0.0f;
#endif
	    tax = dir*(yy + xx + xxx + tx - x*g0);
	    tay = dir*(yz + xy + xxy + ty - y*g0);
	    taz = dir*(zz + xz + xxz + tz - z*g0);
	    /*
	    ** Calculations for determining the timestep.
	    */
	    adotai = pPart->a[0]*tax + pPart->a[1]*tay + pPart->a[2]*taz;
	    if (adotai > 0) {
		adotai *= dimaga;
		dirsum += dir*adotai*adotai;
		normsum += adotai*adotai;
		}
	    ax += tax;
	    ay += tay;
	    az += taz;
	    }
	}
#endif
    pOut->a[0] += ax;
    pOut->a[1] += ay;
    pOut->a[2] += az;
    pOut->fPot += fPot;
    pOut->dirsum += dirsum;
    pOut->normsum += normsum;
    }

int CPUdoWorkPC(void *vpc) {
    workPC *pc = vpc;
    workParticle *work = pc->work;
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
    work->dFlopSingleCPU += COST_FLOP_PC*(tile->lstTile.nBlocks*ILC_PART_PER_BLK  + tile->lstTile.nInLast);
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

static void queuePC( PKD pkd,  workParticle *work, ILC ilc, int bGravStep ) {
    ILCTILE tile;
    workPC *pc;

    ILC_LOOP(ilc,tile) {
#ifdef USE_CUDA
	if (CUDA_queuePC(pkd->mdl->cudaCtx,work,tile,bGravStep)) continue;
#endif
	pc = malloc(sizeof(workPC));
	assert(pc!=NULL);
	pc->pInfoOut = malloc(sizeof(PINFOOUT) * work->nP);
	assert(pc->pInfoOut!=NULL);
	pc->work = work;
	pc->ilc = ilc;
	pc->tile = tile;
	pc->i = 0;
	tile->lstTile.nRefs++;
	work->nRefs++;
	mdlAddWork(pkd->mdl,pc,NULL,NULL,CPUdoWorkPC,doneWorkPC);
	}
    }

int CPUdoWorkEwald(void *ve) {
    workEwald *ew = ve;
    PKD pkd = (PKD)ew->pkd;
    double r[3];
    while(ew->nP--) {
	workParticle *wp = ew->ppWorkPart[ew->nP];
	int wi = ew->piWorkPart[ew->nP];
	//PARTICLE *p = wp->pPart[wi];
	PINFOIN *in = &wp->pInfoIn[wi];
	PINFOOUT *out = &wp->pInfoOut[wi];
	//pkdGetPos1(p->r,r);
	r[0] = wp->c[0] + in->r[0];
	r[1] = wp->c[1] + in->r[1];
	r[2] = wp->c[2] + in->r[2];
	wp->dFlop += pkdParticleEwald(pkd,r,out->a,&out->fPot,&wp->dFlopSingleCPU,&wp->dFlopDoubleCPU);
	pkdParticleWorkDone(wp);
	}
    return 0;
    }

int doneWorkEwald(void *ve) {
    workEwald *ew = ve;
    free(ew->ppWorkPart);
    free(ew->piWorkPart);
    free(ew);
    return 0;
    }

void pkdGravStartEwald(PKD pkd) {
    pkd->ewWork = malloc(sizeof(workEwald));
    assert(pkd->ewWork!=NULL);
    pkd->ewWork->pkd = pkd;
    pkd->ewWork->nP = 0;
    pkd->ewWork->ppWorkPart = malloc(MAX_EWALD_PARTICLES * sizeof(workParticle *));
    pkd->ewWork->piWorkPart = malloc(MAX_EWALD_PARTICLES * sizeof(int));
    assert(pkd->ewWork->ppWorkPart!=NULL);
    assert(pkd->ewWork->piWorkPart!=NULL);
    }

void pkdGravFinishEwald(PKD pkd) {
    /* Finish any Ewald work */
    if (pkd->ewWork->nP) {
#if defined(USE_CL)
	CL_queue(pkd->mdl->clCtx,pkd->ewWork,CLinitWorkEwald,CLcheckWorkEwald,doneWorkEwald);
#elif defined(USE_CUDA)
	mdlAddWork(pkd->mdl,pkd->ewWork,CUDAinitWorkEwald,CUDAcheckWorkEwald,CPUdoWorkEwald,doneWorkEwald);
#else
	mdlAddWork(pkd->mdl,pkd->ewWork,NULL,NULL,CPUdoWorkEwald,doneWorkEwald);
#endif
	}
    pkd->ewWork = NULL;
    }

static void queueEwald( PKD pkd, workParticle *work ) {
    int i;
    for( i=0; i<work->nP; i++ ) {
	if (pkd->ewWork->nP == MAX_EWALD_PARTICLES) {
	    pkdGravFinishEwald(pkd);
	    pkdGravStartEwald(pkd);
	    }
	pkd->ewWork->ppWorkPart[pkd->ewWork->nP] = work;
	pkd->ewWork->piWorkPart[pkd->ewWork->nP] = i;
	++pkd->ewWork->nP;
	++work->nRefs;
	}
    }

/*
** This version of grav.c does all the operations inline, including
** v_sqrt's and such.
** Returns nActive.
*/
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    int bKickClose,int bKickOpen,double dTime,vel_t *dtClose,vel_t *dtOpen,double dAccFac,
    KDN *pBucket,LOCR *pLoc,ILP ilp,ILC ilc,
    float dirLsum,float normLsum,int bEwald,int bGravStep,double *pdFlop,
    double dRhoFac,SMX smx,SMF *smf,int iRoot1,int iRoot2) {
    PARTICLE *p;
    KDN *pkdn = pBucket;
    double r[3], kdn_r[3];
    vel_t *v;
    double vx,vy,vz;
    float *ap;
    float d2,dir,dir2;
    float fMass,fSoft;
    float fx,fy,fz;
    double dx,dy,dz,dPot,ax,ay,az;
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
    workParticle *work = malloc(sizeof(workParticle));
    assert(work!=NULL);
    /* This is the maximum number of particles -- there may be fewer of course */
    nP = pkdn->pUpper - pkdn->pLower + 1;
    work->pPart = malloc(sizeof(PARTICLE *) * nP); assert(work->pPart != NULL);
    work->pInfoIn = malloc(sizeof(PINFOIN) * nP); assert(work->pInfoIn != NULL);
    work->pInfoOut = malloc(sizeof(PINFOOUT) * nP); assert(work->pInfoOut != NULL);
    work->c[0] = ilp->cx; assert(work->c[0] == ilc->cx);
    work->c[1] = ilp->cy; assert(work->c[1] == ilc->cy);
    work->c[2] = ilp->cz; assert(work->c[2] == ilc->cz);

    work->nRefs = 1; /* I am using it currently */
    work->dFlop = 0.0;
    work->dFlopSingleCPU = work->dFlopSingleGPU = 0.0;
    work->dFlopDoubleCPU = work->dFlopDoubleGPU = 0.0;
    work->nP = 0;
    work->dRhoFac = dRhoFac;
    work->ctx = pkd;
    work->bGravStep = bGravStep;
    work->uRungLo = uRungLo;
    work->uRungHi = uRungHi;
    work->bKickClose = bKickClose;
    work->bKickOpen = bKickOpen;
    /*
    ** We copy the pointers here assuming that the storage for them lasts at least as long as
    ** the work structure. Depending on how this is called it could create problems if the 
    ** work is flushed out somewhere else, as is the case when using CUDA.
    */
    work->dTime = dTime;  /* maybe we want the look-back factor. */
    work->dtClose = dtClose;
    work->dtOpen = dtOpen;
    work->dAccFac = dAccFac;
#ifdef USE_CUDA
    work->cudaCtx = pkd->cudaCtx;
#endif

    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;

	pkdGetPos1(pkd,p,r);
	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	v = pkdVel(pkd,p);

	nP = work->nP++;
	work->pPart[nP] = p;

	work->pInfoIn[nP].r[0]  = (float)(r[0] - ilp->cx);
	work->pInfoIn[nP].r[1]  = (float)(r[1] - ilp->cy);
	work->pInfoIn[nP].r[2]  = (float)(r[2] - ilp->cz);
	if (pkd->oAcceleration) {
	    ap = pkdAccel(pkd,p);
	    work->pInfoIn[nP].a[0]  = ap[0];
	    work->pInfoIn[nP].a[1]  = ap[1];
	    work->pInfoIn[nP].a[2]  = ap[2];
	    ap[0] = ap[1] = ap[2] = 0.0;
	    }
	else {
	    work->pInfoIn[nP].a[0]  = 0;
	    work->pInfoIn[nP].a[1]  = 0;
	    work->pInfoIn[nP].a[2]  = 0;
	    }

	work->pInfoOut[nP].a[0] = 0.0f;
	work->pInfoOut[nP].a[1] = 0.0f;
	work->pInfoOut[nP].a[2] = 0.0f;
	work->pInfoOut[nP].fPot = 0.0f;
	work->pInfoOut[nP].dirsum = dirLsum;
	work->pInfoOut[nP].normsum = normLsum;
	work->pInfoOut[nP].rhopmax = 0.0f;

	/*
	** Calculate local density and kernel smoothing length for dynamical time-stepping
	*/
	if (bGravStep) {
	    /*
	    ** Calculate local density using smooth; this is fast because the particles are
	    ** likely to be cached already because they will be on the P-P list.
	    */
	    smf->pfDensity = &work->pInfoIn[nP].fDensity;
	    fBall = smSmoothSingle(smx,smf,p,iRoot1,iRoot2);
	    work->pInfoIn[nP].fSmooth2 = fBall * fBall;
	    }
	else {
	    /*
	    ** We are not using GravStep!
	    */
	    work->pInfoIn[nP].fSmooth2 = 0.0;
	    }
	}

    nActive += work->nP;

    /*
    ** Evaluate the local expansion.
    */
    if (pLoc) {
	for( i=0; i<work->nP; i++ ) {
	    double *c = work->c;
	    float *in = work->pInfoIn[i].r;
	    //pkdGetPos1(work->pPart[i]->r,r);
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

	    work->pInfoOut[i].fPot += dPot;
	    work->pInfoOut[i].a[0] += ax;
	    work->pInfoOut[i].a[1] += ay;
	    work->pInfoOut[i].a[2] += az;
	    }
	}

    /*
    ** Evaluate the P-C interactions
    */
    assert(ilc->cx==ilp->cx && ilc->cy==ilp->cy && ilc->cz==ilp->cz );
    queuePC( pkd,  work, ilc, bGravStep );

    /*
    ** Evaluate the P-P interactions
    */
    queuePP( pkd, work, ilp, bGravStep );

    /*
    ** Calculate the Ewald correction for this particle, if it is required.
    */
    if (bEwald) {
	queueEwald( pkd,  work );
	}

    for( i=0; i<work->nP; i++ ) {
	double *c = work->c;
	float *in = work->pInfoIn[i].r;
	r[0] = c[0] + in[0];
	r[1] = c[1] + in[1];
	r[2] = c[2] + in[2];
        //p = work->pPart[i];
	//pkdGetPos1(p->r,r);
        v = pkdVel(pkd,p);
        /*
        ** Set value for time-step, note that we have the current ewald acceleration
        ** in this now as well!
        */
#ifdef TIMESTEP_CRITICAL
        if (bGravStep && pkd->param.iTimeStepCrit == 1) {
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
			if (p->iOrder < pkd->param.nPartColl || tile->xtr[blk].iOrder.i[prt] < pkd->param.nPartColl) {
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
				fx,fy,fz,vx,vy,vz,pkd->param.dEccFacMax);
			    rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;
			    }
			}
		    }
		}
	    work->pInfoOut[i].rhopmax = rhopmax;
	    }
#endif
        } /* end of i-loop cells & particles */

    pkdParticleWorkDone(work);

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
