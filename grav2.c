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


#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    }

/*
** This is called after work has been done for this particle group.
** If everyone has finished, then the particle is updated.
*/
static void workDone(workParticle *work) {
    int i;
    PARTICLE *p;
    float *pPot, *a;

    if ( --work->nRefs == 0 ) {
	for( i=0; i<work->nP; i++ ) {
	    p = work->pPart[i];
	    a = pkdAccel(work->pkd,p);
	    pPot = pkdPot(work->pkd,p);
	    a[0] = work->pInfoOut[i].a[0];
	    a[1] = work->pInfoOut[i].a[1];
	    a[2] = work->pInfoOut[i].a[2];
	    *pPot = work->pInfoOut[i].fPot;
	    if (work->bGravStep) {
		float fx, fy, fz;
		float maga, dT, dtGrav;
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
		    maga = sqrt(fx*fx + fy*fy + fz*fz);
		    dtGrav = maga*dirsum/normsum;
		    }
		else dtGrav = 0.0;
		dtGrav += work->pkd->param.dPreFacRhoLoc*p->fDensity;
		dtGrav = (work->pInfoOut[i].rhopmax > dtGrav?work->pInfoOut[i].rhopmax:dtGrav);
		if (dtGrav > 0.0) {
		    dT = work->pkd->param.dEta/sqrt(dtGrav*work->dRhoFac);
		    p->uNewRung = pkdDtToRung(dT,work->pkd->param.dDelta,work->pkd->param.iMaxRung-1);
		    }
		else p->uNewRung = 0; /* Assumes current uNewRung is outdated -- not ideal */
		}
	    }
	free(work->pPart);
	free(work->pInfoIn);
	free(work->pInfoOut);
	free(work);
	}
    }


int CPUdoWorkPP(void *vpp) {
    workPP *pp = vpp;
    int nLeft, n;
    int j;
    ILP_BLK *blk;
#if defined(USE_SIMD_PP)
    v_sf t1, t2, t3, pd2;
    v_sf pax, pay, paz, pfx, pfy, pfz, pdx, pdy, pdz;
    v_sf piax, piay, piaz;
    v_sf ppot, pmass, p4soft2;
    v_sf padotai,pimaga,psmooth2,pirsum,pnorms;
#else
    float d2,dx,dy,dz,fourh2,dir,dir2,adotai;
#endif
    float ax,ay,az,fPot,dirsum,normsum;
    float tax,tay,taz;
    float dimaga;

    int i = pp->i;
    ILPTILE tile = pp->tile;
    float fx = pp->work->pInfoIn[i].r[0];
    float fy = pp->work->pInfoIn[i].r[1];
    float fz = pp->work->pInfoIn[i].r[2];
    float fMass = pp->work->pInfoIn[i].fMass;
    float fSoft = pp->work->pInfoIn[i].fSoft;
    float fsmooth2 = pp->work->pInfoIn[i].fSmooth2;
    float *a =pp->work->pInfoIn[i].a;

    dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0) {
	dimaga = 1.0/sqrt(dimaga);
	}
#ifdef USE_SIMD_PP
    pimaga = SIMD_SPLAT(dimaga);

    /*
    ** This is a little trick to speed up the calculation. By setting
    ** unused entries in the list to have a zero mass, the resulting
    ** forces are zero. Setting the distance to a large value avoids
    ** softening the non-existent forces which is slightly faster.
    */
    n = pp->nBlocks;
    j = pp->nInLast;
    for( blk=tile->blk+n; j&SIMD_MASK; j++) {
	blk->dx.f[j] = blk->dy.f[j] = blk->dz.f[j] = 1e18f;
	blk->m.f[j] = 0.0f;
	blk->fourh2.f[j] = 1e-18f;
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
    pmass   = SIMD_SPLAT(fMass);
    p4soft2 = SIMD_SPLAT(4.0*fSoft*fSoft);
    psmooth2= SIMD_SPLAT(fsmooth2);

    blk = tile->blk;
    for( nLeft=pp->nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = ((nLeft ? ILP_PART_PER_BLK : pp->nInLast) + SIMD_MASK) >> SIMD_BITS;
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

    blk = tile->blk;
    for( nLeft=pp->nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = (nLeft ? ILP_PART_PER_BLK : pp->nInLast);
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
    pp->pInfoOut[i].a[0] = ax;
    pp->pInfoOut[i].a[1] = ay;
    pp->pInfoOut[i].a[2] = az;
    pp->pInfoOut[i].fPot = fPot;
    pp->pInfoOut[i].dirsum = dirsum;
    pp->pInfoOut[i].normsum = normsum;

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
    workDone(pp->work);
    free(pp->pInfoOut);
    free(pp);
    return 0;
    }

static void queuePP( PKD pkd, workParticle *work, ILP ilp ) {
    int i;
    ILPTILE tile;
    workPP *pp;

    ILP_LOOP(ilp,tile) {
	pp = malloc(sizeof(workPP));
	assert(pp!=NULL);
	pp->pInfoOut = malloc(sizeof(PINFOOUT) * work->nP);
	assert(pp->pInfoOut!=NULL);
	pp->work = work;
	pp->ilp = ilp;
	pp->tile = tile;
	pp->nBlocks = tile->lstTile.nBlocks;
	pp->nInLast = tile->lstTile.nInLast;
	pp->i = 0;
	tile->lstTile.nRefs++;
	work->nRefs++;
#ifdef USE_CUDA
	pp->gpu_memory = NULL;
	mdlAddWork(pkd->mdl,pp,CUDAinitWorkPP,CUDAcheckWorkPP,CPUdoWorkPP,doneWorkPP);
#else
	mdlAddWork(pkd->mdl,pp,NULL,NULL,CPUdoWorkPP,doneWorkPP);
#endif


	}
    }

int CPUdoWorkPC(void *vpc) {
    workPC *pc = vpc;

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
    v_sf ppot, pmass, p4soft2;
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
    ILC_BLK *blk;
    int j, n, nLeft;

    int i = pc->i;
    ILCTILE ctile = pc->tile;
    float fx = pc->work->pInfoIn[i].r[0];
    float fy = pc->work->pInfoIn[i].r[1];
    float fz = pc->work->pInfoIn[i].r[2];
    float fMass = pc->work->pInfoIn[i].fMass;
    float fSoft = pc->work->pInfoIn[i].fSoft;
    float fsmooth2 = pc->work->pInfoIn[i].fSmooth2;
    float *a =pc->work->pInfoIn[i].a;
    float ax,ay,az,fPot,dirsum,normsum;

    dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0.0) {
        dimaga = 1.0 / sqrtf(dimaga);
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
    pmass   = SIMD_SPLAT(fMass);
    p4soft2 = SIMD_SPLAT(4.0*fSoft*fSoft);
    pimaga  = SIMD_SPLAT(dimaga);

    /* Pad the last value if necessary */
    n = pc->nBlocks;
    j = pc->nInLast;
    for( blk=ctile->blk+n; j&SIMD_MASK; j++) {
	blk->dx.f[j] = blk->dy.f[j] = blk->dz.f[j] = 1e18f;
	blk->m.f[j] = 0.0f;
	blk->u.f[j] = 0.0f;
	}

    blk = ctile->blk;
    for( nLeft=pc->nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = ((nLeft ? ILC_PART_PER_BLK : pc->nInLast) + SIMD_MASK) >> SIMD_BITS;
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
	    xxx = SIMD_MUL(x,SIMD_SUB(SIMD_MUL(consts.onethird.p,xx),zz));
	    xxz = SIMD_MUL(z,SIMD_SUB(xx,SIMD_MUL(consts.onethird.p,zz)));
	    yyy = SIMD_MUL(y,SIMD_SUB(SIMD_MUL(consts.onethird.p,yy),zz));
	    yyz = SIMD_MUL(z,SIMD_SUB(yy,SIMD_MUL(consts.onethird.p,zz)));
	    xx = SIMD_SUB(xx,zz);
	    yy = SIMD_SUB(yy,zz);
	    xxy = SIMD_MUL(y,xx);
	    xyy = SIMD_MUL(x,yy);
	    xyz = SIMD_MUL(xy,z);
	    /*
	    ** Now calculate the interaction up to Hexadecapole order.
	    */

	    tx = SIMD_MUL(blk->xxxx.p[j],xxx);


	    tx = SIMD_ADD(tx,SIMD_MUL(blk->xyyy.p[j],yyy));
	    tx = SIMD_MUL(blk->xxxy.p[j],xxy);
	    tx = SIMD_ADD(tx,SIMD_MUL(blk->xxxz.p[j],xxz));
	    tx = SIMD_ADD(tx,SIMD_MUL(blk->xxyy.p[j],xyy));
	    tx = SIMD_ADD(tx,SIMD_MUL(blk->xxyz.p[j],xyz));
	    tx = SIMD_ADD(tx,SIMD_MUL(blk->xyyz.p[j],yyz));
	    tx = SIMD_MUL(tx,g4);

	    ty = SIMD_MUL(blk->xyyy.p[j],xyy);
	    ty = SIMD_ADD(ty,SIMD_MUL(blk->xxxy.p[j],xxx));
	    ty = SIMD_ADD(ty,SIMD_MUL(blk->yyyy.p[j],yyy));
	    ty = SIMD_ADD(ty,SIMD_MUL(blk->yyyz.p[j],yyz));
	    ty = SIMD_ADD(ty,SIMD_MUL(blk->xxyy.p[j],xxy));
	    ty = SIMD_ADD(ty,SIMD_MUL(blk->xxyz.p[j],xxz));
	    ty = SIMD_ADD(ty,SIMD_MUL(blk->xyyz.p[j],xyz));
	    ty = SIMD_MUL(ty,g4);
	    tz = SIMD_MUL(blk->xxxz.p[j],xxx);
	    tz = SIMD_SUB(tz,SIMD_MUL(blk->xxxx.p[j],xxz));
	    tz = SIMD_SUB(tz,SIMD_MUL(SIMD_ADD(blk->xyyy.p[j],blk->xxxy.p[j]),xyz));
	    tz = SIMD_SUB(tz,SIMD_MUL(blk->yyyy.p[j],yyz));
	    tz = SIMD_ADD(tz,SIMD_MUL(blk->yyyz.p[j],yyy));
	    tz = SIMD_SUB(tz,SIMD_MUL(blk->xxyy.p[j],SIMD_ADD(xxz,yyz)));
	    tz = SIMD_ADD(tz,SIMD_MUL(blk->xxyz.p[j],xxy));
	    tz = SIMD_ADD(tz,SIMD_MUL(blk->xyyz.p[j],xyy));
	    tz = SIMD_MUL(tz,g4);
//	    tx = SIMD_MUL(g4,SIMD_MADD(blk->xxxx.p[j],xxx,SIMD_MADD(blk->xyyy.p[j],yyy,
//			SIMD_MADD(blk->xxxy.p[j],xxy,SIMD_MADD(blk->xxxz.p[j],xxz,
//				SIMD_MADD(blk->xxyy.p[j],xyy,SIMD_MADD(blk->xxyz.p[j],xyz,SIMD_MUL(blk->xyyz.p[j],yyz))))))));
//	    ty = SIMD_MUL(g4,SIMD_MADD(blk->xyyy.p[j],xyy,SIMD_MADD(blk->xxxy.p[j],xxx,
//			SIMD_MADD(blk->yyyy.p[j],yyy,SIMD_MADD(blk->yyyz.p[j],yyz,SIMD_MADD(blk->xxyy.p[j],xxy,
//				    SIMD_MADD(blk->xxyz.p[j],xxz,SIMD_MUL(blk->xyyz.p[j],xyz))))))));
//	    tz = SIMD_MUL(g4,SIMD_NMADD(blk->xxxx.p[j],xxz,SIMD_NMADD(SIMD_ADD(blk->xyyy.p[j],blk->xxxy.p[j]),xyz,
//			SIMD_NMADD(blk->yyyy.p[j],yyz,SIMD_NMADD(blk->xxyy.p[j],SIMD_ADD(xxz,yyz),
//				SIMD_MADD(blk->xxxz.p[j],xxx,SIMD_MADD(blk->yyyz.p[j],yyy,SIMD_MADD(blk->xxyz.p[j],xxy,SIMD_MUL(blk->xyyz.p[j],xyy)))))))));
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

    blk = ctile->blk;
    for( nLeft=pc->nBlocks; nLeft >= 0; --nLeft,++blk ) {
	n = (nLeft ? ILC_PART_PER_BLK : pc->nInLast);
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
	    adotai = pc->work->pInfoIn[i].a[0]*tax + pc->work->pInfoIn[i].a[1]*tay + pc->work->pInfoIn[i].a[2]*taz;
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
    pc->pInfoOut[i].a[0] = ax;
    pc->pInfoOut[i].a[1] = ay;
    pc->pInfoOut[i].a[2] = az;
    pc->pInfoOut[i].fPot = fPot;
    pc->pInfoOut[i].dirsum = dirsum;
    pc->pInfoOut[i].normsum = normsum;
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
    workDone(pc->work);
    free(pc->pInfoOut);
    free(pc);
    return 0;
    }

static void queuePC( PKD pkd,  workParticle *work, ILC ilc ) {
    int i;
    ILCTILE tile;
    workPC *pc;

    ILC_LOOP(ilc,tile) {
	pc = malloc(sizeof(workPC));
	assert(pc!=NULL);
	pc->pInfoOut = malloc(sizeof(PINFOOUT) * work->nP);
	assert(pc->pInfoOut!=NULL);
	pc->work = work;
	pc->ilc = ilc;
	pc->tile = tile;
	pc->nBlocks = tile->lstTile.nBlocks;
	pc->nInLast = tile->lstTile.nInLast;
	pc->i = 0;
	tile->lstTile.nRefs++;
	work->nRefs++;

#ifdef USE_CUDA
	pc->gpu_memory = NULL;
	mdlAddWork(pkd->mdl,pc,CUDAinitWorkPC,CUDAcheckWorkPC,CPUdoWorkPC,doneWorkPC);
#else
	mdlAddWork(pkd->mdl,pc,NULL,NULL,CPUdoWorkPC,doneWorkPC);
#endif
	}
    }


/*
** This version of grav.c does all the operations inline, including
** v_sqrt's and such.
** Returns nActive.
*/
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,KDN *pBucket,LOCR *pLoc,ILP ilp,ILC ilc,
    float dirLsum,float normLsum,int bEwald,int bGravStep,int nGroup,double *pdFlop,double *pdEwFlop,double dRhoFac,
    SMX smx,SMF *smf) {
    PARTICLE *p,*pj;
    KDN *pkdn = pBucket;
    double *v, *vTmp;
    double vx,vy,vz;
    float *a;
    float d2,dir,dir2;
    float fMass,fSoft;
    float fx,fy,fz;
    double dx,dy,dz,dPot,ax,ay,az;
    float dtGrav,dT;
    float rhopmax,rhopmaxlocal,fsmooth2;
    float summ;
    ILPTILE tile;
    int i,j,n,nSoft,nActive,nLeft;
    float fourh2;
    ILPCHECKPT checkPt;
    int nP;

    assert(pkd->oPotential);
    assert(pkd->oAcceleration);

    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;

    /* We need to add these particles to the P-P interaction. Note that self-interactions are ignored. */
    ilpCheckPt(ilp,&checkPt);

    /* Collect the bucket particle information */
    workParticle *work = malloc(sizeof(workParticle));
    assert(work!=NULL);
    /* This is the maximum number of particles -- there may be fewer of course */
    nP = pkdn->pUpper - pkdn->pLower + 1;
    work->pPart = malloc(sizeof(PARTICLE *) * nP); assert(work->pPart != NULL);
    work->pInfoIn = malloc(sizeof(PINFOIN) * nP); assert(work->pInfoIn != NULL);
    work->pInfoOut = malloc(sizeof(PINFOOUT) * nP); assert(work->pInfoOut != NULL);
    work->nRefs = 1; /* I am using it currently */
    work->nP = 0;
    work->dRhoFac = dRhoFac;
    work->pkd = pkd;
    work->bGravStep = bGravStep;
#ifdef USE_CUDA
    work->cudaCtx = pkd->cudaCtx;
    work->gpu_memory = NULL;
#endif
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	v = pkdVel(pkd,p);

	/* Beware of self-interaction - must result in zero acceleration */
	ilpAppend(ilp,p->r[0],p->r[1],p->r[2],fMass,4*fSoft*fSoft,p->iOrder,v[0],v[1],v[2]);

	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;

	nP = work->nP++;
	work->pPart[nP] = p;

	a = pkdAccel(pkd,p);
	work->pInfoIn[nP].r[0]  = p->r[0] - ilp->cx;
	work->pInfoIn[nP].r[1]  = p->r[1] - ilp->cy;
	work->pInfoIn[nP].r[2]  = p->r[2] - ilp->cz;
	work->pInfoIn[nP].a[0]  = a[0];
	work->pInfoIn[nP].a[1]  = a[1];
	work->pInfoIn[nP].a[2]  = a[2];
	work->pInfoIn[nP].fMass = fMass;
	work->pInfoIn[nP].fSoft = fSoft;

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
	    smSmoothSingle(smx,smf,p);
	    work->pInfoIn[nP].fSmooth2 = p->fBall * p->fBall;
	    }
	else {
	    /*
	    ** We are not using GravStep!
	    */
	    work->pInfoIn[nP].fSmooth2 = 0.0;
	    }
	}
    assert(work->nP<=nGroup);

    nActive += work->nP;

    /*
    ** Evaluate the local expansion.
    */
    if (pLoc) {
	for( i=0; i<work->nP; i++ ) {
	    /*
	    ** Evaluate local expansion.
	    */
	    dx = work->pPart[i]->r[0] - pkdn->r[0];
	    dy = work->pPart[i]->r[1] - pkdn->r[1];
	    dz = work->pPart[i]->r[2] - pkdn->r[2];
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
    queuePC( pkd,  work, ilc );

    /*
    ** Evaluate the P-P interactions
    */
    queuePP( pkd, work, ilp );

    /*
    ** Calculate the Ewald correction for this particle, if it is required.
    */
    if (bEwald) {
	for( i=0; i<work->nP; i++ ) {
	    p = work->pPart[i];
	    a = pkdAccel(pkd,p);
	    *pdEwFlop += pkdParticleEwald(pkd,uRungLo,uRungHi,p,
		work->pInfoOut[i].a, &work->pInfoOut[i].fPot );
	    }
	}

    for( i=0; i<work->nP; i++ ) {
        p = work->pPart[i];
        v = pkdVel(pkd,p);
        a = pkdAccel(pkd,p);
        /*
        ** Set value for time-step, note that we have the current ewald acceleration
        ** in this now as well!
        */
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
			    fx = p->r[0] - ilp->cx + tile->blk[blk].dx.f[prt];
			    fy = p->r[1] - ilp->cy + tile->blk[blk].dy.f[prt];
			    fz = p->r[2] - ilp->cz + tile->blk[blk].dz.f[prt];
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
				dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
				dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
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
        } /* end of i-loop cells & particles */

    workDone(work);

    ilpRestore(ilp,&checkPt);
    *pdFlop += nActive*(ilpCount(pkd->ilp)*40 + ilcCount(pkd->ilc)*200) + nSoft*15;
    return(nActive);
    }

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
