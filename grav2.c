#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
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

#ifdef USE_SIMD
static const struct CONSTS {
    v4 zero;
    v4 onequarter;
    v4 onethird;
    v4 half;
    v4 one;
    v4 threehalves;
    v4 three;
    v4 four;
    v4 five;
    v4 seven;
    v4 nine;
    v4 R3_8;
    v4 R45_32;
    v4 R135_16;
    } consts = {
        {SIMD_CONST(0.0)},
	{SIMD_CONST(0.25)},
	{SIMD_CONST(1.0f/3.0f)},
	{SIMD_CONST(0.5)},
	{SIMD_CONST(1.0)},
	{SIMD_CONST(1.5)},
	{SIMD_CONST(3.0)},
	{SIMD_CONST(4.0)},
	{SIMD_CONST(5.0)},
	{SIMD_CONST(7.0)},
	{SIMD_CONST(9.0)},
	{SIMD_CONST(3.0/8.0)},
	{SIMD_CONST(45.0/32.0)},
	{SIMD_CONST(135.0/16.0)},
    };
#endif


#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    }

void pkdGravTilePP(PKD pkd,ILPTILE tile,
    float fx, float fy, float fz,
    float fMass, float fSoft, float fsmooth2, float *a,
    float *ax, float *ay, float *az,
    float *fPot, float *dirsum, float *normsum) {
    int nLeft, nSIMD, n;
    int j;
    ILP_BLK *blk;
#if defined(USE_SIMD_PP)
    v4sf t1, t2, t3, pd2;
    v4sf pax, pay, paz, pfx, pfy, pfz, pdx, pdy, pdz;
    v4sf piax, piay, piaz;
    v4sf ppot, pmass, p4soft2;
    v4sf padotai,pimaga,psmooth2,pirsum,pnorms;
#else
    float d2,dx,dy,dz,fourh2,dir,dir2,tax,tay,taz,adotai;
    int nSoft;
#endif
    float dimaga;

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
    n = tile->nPart/ILP_PART_PER_BLK;
    j = tile->nPart - n*ILP_PART_PER_BLK;
    for( blk=tile->blk+n; j&SIMD_MASK; j++) {
	blk->dx.f[j] = blk->dy.f[j] = blk->dz.f[j] = 1e18f;
	blk->m.f[j] = 0.0f;
	blk->fourh2.f[j] = 1e-18f;
	}
    nSIMD = (tile->nPart+SIMD_MASK) >> SIMD_BITS;

    /*
    ** The list sets mass to zero for unused entries which results
    ** in zero forces. Be careful if that is changed.
    */
    pax = SIMD_LOADS(*ax);
    pay = SIMD_LOADS(*ay);
    paz = SIMD_LOADS(*az);
    ppot= SIMD_LOADS(*fPot);
    pirsum = SIMD_LOADS(*dirsum);
    pnorms = SIMD_LOADS(*normsum);

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
    for( nLeft=nSIMD; nLeft > 0; nLeft -= n ) {
	n = nLeft > (ILP_PART_PER_BLK/SIMD_WIDTH) ? (ILP_PART_PER_BLK/SIMD_WIDTH) : nLeft;
	for (j=0; j<n; ++j) {
	    v4sf pfourh2, td2, pir, pir2;
	    v4bool vcmp;
	    int msk;

	    pdx = SIMD_ADD(blk->dx.p[j],pfx);
	    pdy = SIMD_ADD(blk->dy.p[j],pfy);
	    pdz = SIMD_ADD(blk->dz.p[j],pfz);
	    pd2 = SIMD_MADD(pdz,pdz,SIMD_MADD(pdy,pdy,SIMD_MUL(pdx,pdx)));
	    vcmp = SIMD_CMP_GT(blk->fourh2.p[j],consts.zero.p);
	    msk = SIMD_ALL_ZERO(vcmp); /* softenings are not zero */
	    if(msk) {
		t1 = SIMD_MUL(SIMD_ADD(pmass,blk->m.p[j]),SIMD_MUL(p4soft2,blk->fourh2.p[j]));
		t2 = SIMD_ADD(SIMD_MUL(blk->fourh2.p[j],pmass),SIMD_MUL(p4soft2,blk->m.p[j]));
#if defined(__SSE2__) || defined(__ALTIVEC__)
		pfourh2 = SIMD_RE_EXACT(t2);
		pfourh2 = SIMD_MUL(pfourh2,t1);
#else
		pfourh2 = SIMD_DIV(t1,t2);
#endif
		vcmp = SIMD_CMP_LT(pd2,pfourh2);
		td2 = SIMD_MAX(pd2,pfourh2);
		msk = SIMD_ALL_ZERO(vcmp);  /* zero means nothing is softened - optimization */
		}
	    else {
		td2 = pd2;
		}

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

	    t1 = SIMD_NMSUB(pdx,pir2,consts.zero.p);
	    t2 = SIMD_NMSUB(pdy,pir2,consts.zero.p);
	    t3 = SIMD_NMSUB(pdz,pir2,consts.zero.p);
	    ppot = SIMD_NMSUB(blk->m.p[j],pir,ppot);

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
	blk++;
	}
    *ax = SIMD_HADD(pax);
    *ay = SIMD_HADD(pay);
    *az = SIMD_HADD(paz);
    *fPot = SIMD_HADD(ppot);
    *dirsum = SIMD_HADD(pirsum);
    *normsum = SIMD_HADD(pnorms);
#else
    /*
    ** DO NOT MODIFY THE CODE BELOW UP TO THE #endif!
    ** This code MUST match the SIMD code above.
    */
    blk = tile->blk;
    for( nLeft=tile->nPart; nLeft > 0; nLeft -= n ) {
	n = nLeft > ILP_PART_PER_BLK ? ILP_PART_PER_BLK : nLeft;
	for (j=0; j<n; ++j) {
	    dx = fx + blk->dx.f[j];
	    dy = fy + blk->dy.f[j];
	    dz = fz + blk->dz.f[j];
	    fourh2 = softmassweight(fMass,4*fSoft*fSoft,
		blk->m.f[j],blk->fourh2.f[j]);
	    d2 = dx*dx + dy*dy + dz*dz;
	    if (d2==0.0) dir2 = 0.0;
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
		++nSoft;
		}
	    dir2 *= blk->m.f[j];
	    tax = -dx*dir2;
	    tay = -dy*dir2;
	    taz = -dz*dir2;
	    *fPot -= blk->m.f[j]*dir;
	    /*
	    ** Calculations for determining the timestep.
	    */
	    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
	    if (adotai > 0 && d2 >= fsmooth2) {
		adotai *= dimaga;
		*dirsum += dir*adotai*adotai;
		*normsum += adotai*adotai;
		}
	    *ax += tax;
	    *ay += tay;
	    *az += taz;
	    }
	blk++;
	}
#endif
    }

static void evaluatePP( PKD pkd, int nP, PARTICLE **pPart, PINFOIN *pInfoIn, ILP ilp ) {
    int i;
    float *a, *pPot;
    float dirsum, normsum;
    ILPTILE tile;

    dirsum = normsum = 0.0;

    for( i=0; i<nP; i++ ) {
	a = pkdAccel(pkd,pPart[i]);
	pPot = pkdPot(pkd,pPart[i]);
	ILP_LOOP(ilp,tile) {
	    pkdGravTilePP(pkd,tile,
		pInfoIn[i].r[0], pInfoIn[i].r[1], pInfoIn[i].r[2],
		pInfoIn[i].fMass,pInfoIn[i].fSoft,pInfoIn[i].fSmooth2,
		pInfoIn[i].a,a+0,a+1,a+2,pPot,&dirsum,&normsum);
	    }
	}
    }

/*
** This version of grav.c does all the operations inline, including
** v_sqrt's and such.
** Returns nActive.
*/
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,KDN *pBucket,FLOCR *pLoc,ILP ilp,ILC ilc,
		    float dirLsum,float normLsum,int bEwald,double *pdFlop,double *pdEwFlop,double dRhoFac,
		    SMX smx,SMF *smf) {
    PARTICLE *p,*pj;
    KDN *pkdn = pBucket;
    double *v, *vTmp;
    double vx,vy,vz;
    pBND bnd;
    float *a, *pPot;
    float ax,ay,az,fPot;
    float d2,dir,dir2;
    float fMass,fSoft;
    float fx, fy, fz;
    float dtGrav,dT;
    float maga,dimaga,dirsum,normsum;
    float rhopmax,rhopmaxlocal,fsmooth2;
    float summ;
#if defined(USE_SIMD_PC)
    v4sf u,g0,g2,g3,g4;
    v4sf x,y,z;
    v4sf tx,ty,tz;
    v4sf xx,xy,xz,yy,yz,zz;
    v4sf xxx,xxz,yyy,yyz,xxy,xyy,xyz;
#else
    const float onethird = 1.0f/3.0f;
    float u,g0,g2,g3,g4;
    float x,y,z;
    float tax,tay,taz,adotai;
    float tx,ty,tz;
    float xx,xy,xz,yy,yz,zz;
    float xxx,xxz,yyy,yyz,xxy,xyy,xyz;
#endif
    ILPTILE tile;
    ILCTILE ctile;
    int i,j,nSoft,nActive;
#if defined(USE_SIMD)
    v4sf t1, t2, t3;
    v4sf pax, pay, paz;
    v4sf pdx, pdy, pdz;
    v4sf pfx, pfy, pfz;
    v4sf piax, piay, piaz;
    v4sf ppot, pmass, p4soft2;
    v4sf padotai,pimaga,pirsum,pnorms;
#endif
    float fourh2;
    ILPCHECKPT checkPt;

    assert(pkd->oPotential);
    assert(pkd->oAcceleration);

    pkdNodeBnd(pkd, pkdn, &bnd);
    normsum = dirsum = maga = 0.0;

    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;


    /* We need to add these particles to the P-P interaction. Note that self-interactions are ignored. */
    ilpCheckPt(ilp,&checkPt);

    /* Collect the bucket particle information */
    int nP;
    PINFOIN pInfoIn[PKD_GROUP_SIZE];
    PARTICLE *pPart[PKD_GROUP_SIZE];
    nP = 0;
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;
	pPart[nP] = p;

	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	pPot = pkdPot(pkd,p);
	v = pkdVel(pkd,p);
	a = pkdAccel(pkd,p);
	pInfoIn[nP].r[0]  = p->r[0] - ilp->cx;
	pInfoIn[nP].r[1]  = p->r[1] - ilp->cy;
	pInfoIn[nP].r[2]  = p->r[2] - ilp->cz;
	pInfoIn[nP].a[0]  = a[0];
	pInfoIn[nP].a[1]  = a[1];
	pInfoIn[nP].a[2]  = a[2];
	pInfoIn[nP].fMass = fMass;
	pInfoIn[nP].fSoft = fSoft;

	a[0] = a[1] = a[2] = *pPot = 0.0;

	/*
	** Calculate local density and kernel smoothing length for dynamical time-stepping
	*/
	if (pkd->param.bGravStep) {
	    /*
	    ** Calculate local density using smooth; this is fast because the particles are
	    ** likely to be cached already because they will be on the P-P list.
	    */
	    smSmoothSingle(smx,smf,p);
	    pInfoIn[nP].fSmooth2 = p->fBall * p->fBall;
	    }
	else {
	    /*
	    ** We are not using GravStep!
	    */
	    pInfoIn[nP].fSmooth2 = 0.0;
	    }

	/* Beware of self-interaction - must result in zero acceleration */
	ilpAppend(ilp,p->r[0],p->r[1],p->r[2],fMass,4*fSoft*fSoft,p->iOrder,v[0],v[1],v[2]);
	++nP;
	}
    assert(nP<=PKD_GROUP_SIZE);

    nActive += nP;

    /*
    ** Evaluate the local expansion.
    */
    if (pLoc) {
	for( i=0; i<nP; i++ ) {
	    a = pkdAccel(pkd,pPart[i]);
	    pPot = pkdPot(pkd,pPart[i]);
	    /*
	    ** Evaluate local expansion.
	    */
	    fx = pPart[i]->r[0] - pkdn->r[0];
	    fy = pPart[i]->r[1] - pkdn->r[1];
	    fz = pPart[i]->r[2] - pkdn->r[2];
	    momEvalFlocr(pLoc,pkdn->bMax,fx,fy,fz,pPot,a+0,a+1,a+2);
	    a = pkdAccel(pkd,pPart[i]);
	    }
	}

    /*
    ** Evaluate the P-C interactions
    */
    for( i=0; i<nP; i++ ) {
	pInfoIn[i].r[0]  = pPart[i]->r[0] - ilc->cx;
	pInfoIn[i].r[1]  = pPart[i]->r[1] - ilc->cy;
	pInfoIn[i].r[2]  = pPart[i]->r[2] - ilc->cz;
	}

    /* Padd the ILC if required */
#if defined(USE_SIMD_PC)
    ILC_LOOP(ilc,ctile) {
	for( j=ctile->nCell; j&SIMD_MASK; ++j) { 
	    ctile->dx.f[j] = ctile->dy.f[j] = ctile->dz.f[j] = 1e18f;
	    ctile->m.f[j] = 0.0f;
	    ctile->u.f[j] = 0.0f;
	    }
	}
#endif

    for( i=0; i<nP; i++ ) {
	p = pPart[i];
	a = pkdAccel(pkd,pPart[i]);
	pPot = pkdPot(pkd,pPart[i]);

	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
	fx = pInfoIn[i].a[0];
	fy = pInfoIn[i].a[1];
	fz = pInfoIn[i].a[2];

	maga = fx*fx + fy*fy + fz*fz;
	dimaga = maga;
	if (dimaga > 0) {
	    dimaga = 1.0/sqrt(dimaga);
	    }
#ifdef USE_SIMD
	pimaga = SIMD_SPLAT(dimaga);
#endif
	dirsum = dirLsum;
	normsum = normLsum;

	fx = pInfoIn[i].r[0];
	fy = pInfoIn[i].r[1];
	fz = pInfoIn[i].r[2];

//	ilc->cx = p->r[0]; /* => cx += fx */
//	ilc->cy = p->r[1];
//	ilc->cz = p->r[2];
//	ilcCompute(ilc,fx,fy,fz);
#if defined(USE_SIMD_PC)
	pax = SIMD_LOADS(ax);
	pay = SIMD_LOADS(ay);
	paz = SIMD_LOADS(az);
	ppot= SIMD_LOADS(fPot);
	pirsum = SIMD_LOADS(dirsum);
	pnorms = SIMD_LOADS(normsum);
	pfx     = SIMD_SPLAT(fx);
	pfy     = SIMD_SPLAT(fy);
	pfz     = SIMD_SPLAT(fz);
	piax    = SIMD_SPLAT(pInfoIn[i].a[0]);
	piay    = SIMD_SPLAT(pInfoIn[i].a[1]);
	piaz    = SIMD_SPLAT(pInfoIn[i].a[2]);
	pmass   = SIMD_SPLAT(pInfoIn[i].fMass);
	p4soft2 = SIMD_SPLAT(4.0*pInfoIn[i].fSoft*pInfoIn[i].fSoft);

	ILC_LOOP(ilc,ctile) {
	    uint32_t n = (ctile->nCell+SIMD_MASK) >> SIMD_BITS;
	    for (j=0;j<n;++j) {
		v4sf pir, pd2;
		v4bool vcmp;

		pdx = SIMD_ADD(ctile->dx.p[j],pfx);
		pdy = SIMD_ADD(ctile->dy.p[j],pfy);
		pdz = SIMD_ADD(ctile->dz.p[j],pfz);
		pir = SIMD_RSQRT_EXACT(SIMD_MADD(pdz,pdz,SIMD_MADD(pdy,pdy,SIMD_MUL(pdx,pdx))));
		u = SIMD_MUL(ctile->u.p[j],pir);
		g0 = pir;
		g2 = SIMD_MUL(SIMD_MUL(consts.three.p,pir),SIMD_MUL(u,u));
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
		tx = SIMD_MUL(g4,SIMD_MADD(ctile->xxxx.p[j],xxx,SIMD_MADD(ctile->xyyy.p[j],yyy,
		    SIMD_MADD(ctile->xxxy.p[j],xxy,SIMD_MADD(ctile->xxxz.p[j],xxz,
		    SIMD_MADD(ctile->xxyy.p[j],xyy,SIMD_MADD(ctile->xxyz.p[j],xyz,SIMD_MUL(ctile->xyyz.p[j],yyz))))))));
		ty = SIMD_MUL(g4,SIMD_MADD(ctile->xyyy.p[j],xyy,SIMD_MADD(ctile->xxxy.p[j],xxx,
		    SIMD_MADD(ctile->yyyy.p[j],yyy,SIMD_MADD(ctile->yyyz.p[j],yyz,SIMD_MADD(ctile->xxyy.p[j],xxy,
		    SIMD_MADD(ctile->xxyz.p[j],xxz,SIMD_MUL(ctile->xyyz.p[j],xyz))))))));
		tz = SIMD_MUL(g4,SIMD_NMSUB(ctile->xxxx.p[j],xxz,SIMD_NMSUB(SIMD_ADD(ctile->xyyy.p[j],ctile->xxxy.p[j]),xyz,
			SIMD_NMSUB(ctile->yyyy.p[j],yyz,SIMD_NMSUB(ctile->xxyy.p[j],SIMD_ADD(xxz,yyz),
			SIMD_MADD(ctile->xxxz.p[j],xxx,SIMD_MADD(ctile->yyyz.p[j],yyy,SIMD_MADD(ctile->xxyz.p[j],xxy,SIMD_MUL(ctile->xyyz.p[j],xyy)))))))));
		g4 = SIMD_MUL(consts.onequarter.p,SIMD_MADD(tx,x,SIMD_MADD(ty,y,SIMD_MUL(tz,z))));
		xxx = SIMD_MUL(g3,SIMD_MADD(ctile->xxx.p[j],xx,SIMD_MADD(ctile->xyy.p[j],yy,
		    SIMD_MADD(ctile->xxy.p[j],xy,SIMD_MADD(ctile->xxz.p[j],xz,SIMD_MUL(ctile->xyz.p[j],yz))))));
		xxy = SIMD_MUL(g3,SIMD_MADD(ctile->xyy.p[j],xy,SIMD_MADD(ctile->xxy.p[j],xx,SIMD_MADD(ctile->yyy.p[j],yy,
		    SIMD_MADD(ctile->yyz.p[j],yz,SIMD_MUL(ctile->xyz.p[j],xz))))));
		xxz = SIMD_MUL(g3,SIMD_NMSUB(SIMD_ADD(ctile->xxx.p[j],ctile->xyy.p[j]),xz,
			    SIMD_NMSUB(SIMD_ADD(ctile->xxy.p[j],ctile->yyy.p[j]),yz,
			    SIMD_MADD(ctile->xxz.p[j],xx,SIMD_MADD(ctile->yyz.p[j],yy,SIMD_MUL(ctile->xyz.p[j],xy))))));
		g3 = SIMD_MUL(consts.onethird.p,SIMD_MADD(xxx,x,SIMD_MADD(xxy,y,SIMD_MUL(xxz,z))));
		xx = SIMD_MUL(g2,SIMD_MADD(ctile->xx.p[j],x,SIMD_MADD(ctile->xy.p[j],y,SIMD_MUL(ctile->xz.p[j],z))));
		xy = SIMD_MUL(g2,SIMD_MADD(ctile->yy.p[j],y,SIMD_MADD(ctile->xy.p[j],x,SIMD_MUL(ctile->yz.p[j],z))));
		xz = SIMD_MUL(g2,SIMD_NMSUB(SIMD_ADD(ctile->xx.p[j],ctile->yy.p[j]),z,SIMD_MADD(ctile->xz.p[j],x,SIMD_MUL(ctile->yz.p[j],y))));
		g2 = SIMD_MUL(consts.half.p,SIMD_MADD(xx,x,SIMD_MADD(xy,y,SIMD_MUL(xz,z))));
		g0 = SIMD_MUL(g0,ctile->m.p[j]);
		ppot = SIMD_SUB(ppot,SIMD_ADD(SIMD_ADD(g0,g2),SIMD_ADD(g3,g4)));
		g0 = SIMD_MADD(consts.five.p,g2,SIMD_MADD(consts.seven.p,g3,SIMD_MADD(consts.nine.p,g4,g0)));
		t1 = SIMD_MUL(pir,SIMD_NMSUB(x,g0,SIMD_ADD(xx,SIMD_ADD(xxx,tx))));
		t2 = SIMD_MUL(pir,SIMD_NMSUB(y,g0,SIMD_ADD(xy,SIMD_ADD(xxy,ty))));
		t3 = SIMD_MUL(pir,SIMD_NMSUB(z,g0,SIMD_ADD(xz,SIMD_ADD(xxz,tz))));

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
	ILC_LOOP(ilc,ctile) {
	    for (j=0;j<ctile->nCell;++j) {
		float dx = ctile->dx.f[j] + fx;
		float dy = ctile->dy.f[j] + fy;
		float dz = ctile->dz.f[j] + fz;
		float d2 = dx*dx + dy*dy + dz*dz;
		SQRT1(d2,dir);
		u = ctile->u.f[j]*dir;
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
		tx = g4*(ctile->xxxx.f[j]*xxx + ctile->xyyy.f[j]*yyy + ctile->xxxy.f[j]*xxy + ctile->xxxz.f[j]*xxz + ctile->xxyy.f[j]*xyy + ctile->xxyz.f[j]*xyz + ctile->xyyz.f[j]*yyz);
		ty = g4*(ctile->xyyy.f[j]*xyy + ctile->xxxy.f[j]*xxx + ctile->yyyy.f[j]*yyy + ctile->yyyz.f[j]*yyz + ctile->xxyy.f[j]*xxy + ctile->xxyz.f[j]*xxz + ctile->xyyz.f[j]*xyz);
		tz = g4*(-ctile->xxxx.f[j]*xxz - (ctile->xyyy.f[j] + ctile->xxxy.f[j])*xyz - ctile->yyyy.f[j]*yyz + ctile->xxxz.f[j]*xxx + ctile->yyyz.f[j]*yyy - ctile->xxyy.f[j]*(xxz + yyz) + ctile->xxyz.f[j]*xxy + ctile->xyyz.f[j]*xyy);
		g4 = 0.25*(tx*x + ty*y + tz*z);
		xxx = g3*(ctile->xxx.f[j]*xx + ctile->xyy.f[j]*yy + ctile->xxy.f[j]*xy + ctile->xxz.f[j]*xz + ctile->xyz.f[j]*yz);
		xxy = g3*(ctile->xyy.f[j]*xy + ctile->xxy.f[j]*xx + ctile->yyy.f[j]*yy + ctile->yyz.f[j]*yz + ctile->xyz.f[j]*xz);
		xxz = g3*(-(ctile->xxx.f[j] + ctile->xyy.f[j])*xz - (ctile->xxy.f[j] + ctile->yyy.f[j])*yz + ctile->xxz.f[j]*xx + ctile->yyz.f[j]*yy + ctile->xyz.f[j]*xy);
		g3 = onethird*(xxx*x + xxy*y + xxz*z);
		xx = g2*(ctile->xx.f[j]*x + ctile->xy.f[j]*y + ctile->xz.f[j]*z);
		xy = g2*(ctile->yy.f[j]*y + ctile->xy.f[j]*x + ctile->yz.f[j]*z);
		xz = g2*(-(ctile->xx.f[j] + ctile->yy.f[j])*z + ctile->xz.f[j]*x + ctile->yz.f[j]*y);
		g2 = 0.5*(xx*x + xy*y + xz*z);
		g0 *= ctile->m.f[j];
		fPot -= g0 + g2 + g3 + g4;
		g0 += 5*g2 + 7*g3 + 9*g4;
		tax = dir*(xx + xxx + tx - x*g0);
		tay = dir*(xy + xxy + ty - y*g0);
		taz = dir*(xz + xxz + tz - z*g0);
		/*
		** Calculations for determining the timestep.
		*/
		adotai = pInfoIn[i].a[0]*tax + pInfoIn[i].a[1]*tay + pInfoIn[i].a[2]*taz;
		if (adotai > 0) {
		    adotai *= dimaga;
		    dirsum += dir*adotai*adotai;
		    normsum += adotai*adotai;
		    }
		ax += tax;
		ay += tay;
		az += taz;
		}
	    } /* end of cell list gravity loop */

#endif
	a[0] += ax;
	a[1] += ay;
	a[2] += az;
	mdlCacheCheck(pkd->mdl);
	}

    /*
    ** Evaluate the P-P interactions
    */
    for( i=0; i<nP; i++ ) {
	pInfoIn[i].r[0]  = pPart[i]->r[0] - ilp->cx;
	pInfoIn[i].r[1]  = pPart[i]->r[1] - ilp->cy;
	pInfoIn[i].r[2]  = pPart[i]->r[2] - ilp->cz;
	}

#ifdef USE_CUDA
    pkdGravCudaPP( pkd, pkd->cudaCtx, nP, pPart, pInfoIn, ilp );
#else
    evaluatePP( pkd, nP, pPart, pInfoIn, ilp );
#endif

    for( i=0; i<nP; i++ ) {
	p = pPart[i];
	a = pkdAccel(pkd,pPart[i]);
	pPot = pkdPot(pkd,pPart[i]);

	/*
	** Now finally calculate the Ewald correction for this particle, if it is
	** required.
	*/
	if (bEwald) {
	    *pdEwFlop += pkdParticleEwald(pkd,uRungLo,uRungHi,p);
	    }
	/*
	** Set value for time-step, note that we have the current ewald acceleration
	** in this now as well!
	*/
	if (pkd->param.bGravStep) {
	    /*
	    ** If this is the first time through, the accelerations will have 
	    ** all been zero resulting in zero for normsum (and nan for dtGrav).
	    ** We repeat this process again, so dtGrav will be correct.
	    */
	    if (normsum > 0.0) {
		/*
		** Use new acceleration here!
		*/
		fx = pInfoIn[i].a[0];
		fy = pInfoIn[i].a[1];
		fz = pInfoIn[i].a[2];
		maga = sqrt(fx*fx + fy*fy + fz*fz);
		dtGrav = maga*dirsum/normsum;
		}
	    else dtGrav = 0.0;
	    dtGrav += pkd->param.dPreFacRhoLoc*p->fDensity;
	    if (pkd->param.iTimeStepCrit == 1) {
	      /*
	      ** GravStep if iTimeStepCrit =
	      ** 0: Mean field regime for dynamical time (normal/standard setting)
	      ** 1: Gravitational scattering regime for dynamical time with eccentricity correction
	      */
	      rhopmax = 0.0;
	      ILP_LOOP(ilp,tile) {
		for (j=0;j<tile->nPart;++j) {
		  int blk = j / ILP_PART_PER_BLK;
		  int prt = j - blk * ILP_PART_PER_BLK;
		  if (p->iOrder < pkd->param.nPartColl || tile->iOrder.i[j] < pkd->param.nPartColl) {
		    fx = p->r[0] - ilp->cx + tile->blk[blk].dx.f[prt];
		    fy = p->r[1] - ilp->cy + tile->blk[blk].dy.f[prt];
		    fz = p->r[2] - ilp->cz + tile->blk[blk].dz.f[prt];
		    d2 = fx*fx + fy*fy + fz*fz;
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
		    vx = v[0] - tile->vx.f[j];
		    vy = v[1] - tile->vy.f[j];
		    vz = v[2] - tile->vz.f[j];
		    rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,
			tile->blk[blk].dx.f[prt],tile->blk[blk].dy.f[prt],tile->blk[blk].dz.f[prt],
			vx,vy,vz,pkd->param.dEccFacMax);
		    rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;
		  }
		}
	      }
	      dtGrav = (rhopmax > dtGrav?rhopmax:dtGrav);
	    }
	    if (dtGrav > 0.0) {
		dT = pkd->param.dEta/sqrt(dtGrav*dRhoFac);
		p->uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
		}
	    else p->uNewRung = 0; /* Assumes current uNewRung is outdated -- not ideal */
	    }
	/*
	** Restore the ilp to the same state that GravInteract was called.
	*/
	} /* end of i-loop cells & particles */

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
