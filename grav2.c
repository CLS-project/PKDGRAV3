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
	{{0.0,     0.0,     0.0,     0.0}},
	{{0.25,    0.25,    0.25,    0.25}},
	{{1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f}},
	{{0.5,     0.5,     0.5,     0.5}},
	{{1.0,     1.0,     1.0,     1.0}},
	{{1.5,     1.5,     1.5,     1.5}},
	{{3.0,     3.0,     3.0,     3.0}},
	{{4.0,     4.0,     4.0,     4.0}},
	{{5.0,     5.0,     5.0,     5.0}},
	{{7.0,     7.0,     7.0,     7.0}},
	{{9.0,     9.0,     9.0,     9.0}},
	{{3.0/8.0, 3.0/8.0, 3.0/8.0, 3.0/8.0}},
	{{45.0/32.0, 45.0/32.0, 45.0/32.0, 45.0/32.0}},
	{{135.0/16.0, 135.0/16.0, 135.0/16.0, 135.0/16.0}},
    };
#endif


#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
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
    float fMassTmp,fSoftTmp;
    float fx, fy, fz;
    float dtGrav,dT;
    float adotai,maga,dimaga,dirsum,normsum;
    float tax,tay,taz;
    float rholoc,rhopmax,rhopmaxlocal,fsmooth2,fSoftMedian,fEps,fEps2;
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
    float tx,ty,tz;
    float xx,xy,xz,yy,yz,zz;
    float xxx,xxz,yyy,yyz,xxy,xyy,xyz;
#endif
    ILPTILE tile;
    ILCTILE ctile;
    int i,j,nSoft,nActive;
    float rMax;
#if defined(USE_SIMD_PP)
    v4sf t1, t2, t3;
    v4sf pax, pay, paz;
    v4sf piax, piay, piaz;
    v4sf ppot, pmass, p4soft2;
    v4sf padotai,pimaga,psmooth2,pirsum,pnorms;
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
    rMax = pkd->param.nPartRhoLoc/(pkdn->pUpper-pkdn->pLower+1)
	   *(bnd.fMax[0]*bnd.fMax[0]
	     + bnd.fMax[1]*bnd.fMax[1]
	     + bnd.fMax[2]*bnd.fMax[2]);
    /*
    ** Save the ilp list so that we can restore it for each particle. We will be 
    ** adding the local bucket interactions within this loop, and so modifying the 
    ** ilp.
    */
    ilpCheckPt(ilp,&checkPt);
    p = pkdParticle(pkd,pkdn->pLower);
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;
	pPot = pkdPot(pkd,p);
	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	v = pkdVel(pkd,p);
	a = pkdAccel(pkd,p);
	++nActive;
	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
	fx = a[0];
	fy = a[1];
	fz = a[2];
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
	if (pLoc) {
	    /*
	    ** Evaluate local expansion.
	    */
	    fx = p->r[0] - pkdn->r[0];
	    fy = p->r[1] - pkdn->r[1];
	    fz = p->r[2] - pkdn->r[2];
	    momEvalFlocr(pLoc,pkdn->bMax,fx,fy,fz,&fPot,&ax,&ay,&az);
	    }

	fx = p->r[0] - ilc->cx;
	fy = p->r[1] - ilc->cy;
	fz = p->r[2] - ilc->cz;

	ilc->cx = p->r[0]; /* => cx += fx */
	ilc->cy = p->r[1];
	ilc->cz = p->r[2];
	ilcCompute(ilc,fx,fy,fz);
#if defined(USE_SIMD_PC)
	pax = SIMD_LOADS(ax);
	pay = SIMD_LOADS(ay);
	paz = SIMD_LOADS(az);
	ppot= SIMD_LOADS(fPot);
	pirsum = SIMD_LOADS(dirsum);
	pnorms = SIMD_LOADS(normsum);

	piax    = SIMD_SPLAT(a[0]);
	piay    = SIMD_SPLAT(a[1]);
	piaz    = SIMD_SPLAT(a[2]);
	pmass   = SIMD_SPLAT(fMass);
	p4soft2 = SIMD_SPLAT(4.0*fSoft*fSoft);

	ILC_LOOP(ilc,ctile) {
	    uint32_t n = (ctile->nCell+ILP_ALIGN_MASK) >> ILP_ALIGN_BITS;
	    for (j=0;j<n;++j) {
		v4sf pir, pd2;
		v4bool vcmp;
		pir = SIMD_RSQRT_EXACT(ctile->d2.p[j]);
		u = SIMD_MUL(ctile->u.p[j],pir);
		g0 = pir;
		g2 = SIMD_MUL(SIMD_MUL(consts.three.p,pir),SIMD_MUL(u,u));
		g3 = SIMD_MUL(consts.five.p,SIMD_MUL(g2,u));
		g4 = SIMD_MUL(consts.seven.p,SIMD_MUL(g3,u));
		/*
		** Calculate the funky distance terms.
		*/
		x = SIMD_MUL(ctile->dx.p[j],pir);
		y = SIMD_MUL(ctile->dy.p[j],pir);
		z = SIMD_MUL(ctile->dz.p[j],pir);
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
		SQRT1(ctile->d2.f[j],dir);
		u = ctile->u.f[j]*dir;
		g0 = dir;
		g2 = 3*dir*u*u;
		g3 = 5*g2*u;
		g4 = 7*g3*u;
		/*
		** Calculate the funky distance terms.
		*/
		x = ctile->dx.f[j]*dir;
		y = ctile->dy.f[j]*dir;
		z = ctile->dz.f[j]*dir;
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
		adotai = a[0]*tax + a[1]*tay + a[2]*taz;
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
	mdlCacheCheck(pkd->mdl);
	/*
	** Part 1: Calculate distance between particle and each interaction
	*/
	for (j=pkdn->pLower;j<=pkdn->pUpper;++j) {
	    if (j == i) continue;
	    pj = pkdParticle(pkd,j);
	    if ( !pkdIsSrcActive(pj,uRungLo,uRungHi) ) continue;
	    fMassTmp = pkdMass(pkd,pj);
	    fSoftTmp = pkdSoft(pkd,pj);
	    vTmp = pkdVel(pkd,pj);
	    ilpAppend(ilp,pj->r[0],pj->r[1],pj->r[2],fMassTmp,4*fSoftTmp*fSoftTmp,pj->iOrder,vTmp[0],vTmp[1],vTmp[2]);
	    }

	fx = p->r[0] - ilp->cx;
	fy = p->r[1] - ilp->cy;
	fz = p->r[2] - ilp->cz;

	ilp->cx = p->r[0]; /* => cx += fx */
	ilp->cy = p->r[1];
	ilp->cz = p->r[2];
	ilpCompute(ilp,fx,fy,fz);
	/*
	** Calculate local density and kernel smoothing length for dynamical time-stepping
	*/
	if (pkd->param.bGravStep) {
	    /*
	    ** Calculate local density using smooth; this is fast because the particles are
	    ** likely to be cached already because they will be on the P-P list.
	    */
	    smSmoothSingle(smx,smf,p);
	    fsmooth2 = p->fBall * p->fBall;
#ifdef USE_SIMD
	    psmooth2 = SIMD_SPLAT(fsmooth2);
#endif
	    }
	 else {
	   /*
	   ** We are not using GravStep!
	   */
	    fsmooth2 = 0.0;
#ifdef USE_SIMD_PP
	    psmooth2 = consts.zero.p;
#endif	    
	    }
	
#ifdef USE_SIMD_PP
	 /*
	 ** The list sets mass to zero for unused entries which results
	 ** in zero forces. Be careful if that is changed.
	 */
	pax = SIMD_LOADS(ax);
	pay = SIMD_LOADS(ay);
	paz = SIMD_LOADS(az);
	ppot= SIMD_LOADS(fPot);
	pirsum = SIMD_LOADS(dirsum);
	pnorms = SIMD_LOADS(normsum);

	piax    = SIMD_SPLAT(a[0]);
	piay    = SIMD_SPLAT(a[1]);
	piaz    = SIMD_SPLAT(a[2]);
	pmass   = SIMD_SPLAT(fMass);
	p4soft2 = SIMD_SPLAT(4.0*fSoft*fSoft);

	ILP_LOOP(ilp,tile) {
	    uint32_t n = (tile->nPart+ILP_ALIGN_MASK) >> ILP_ALIGN_BITS;

	    for ( j=0; j<n; j++ ) {
		v4sf pfourh2, pd2, pir, pir2;
		v4bool vcmp;
		int msk;

		vcmp = SIMD_CMP_GT(tile->fourh2.p[j],consts.zero.p);
		msk = SIMD_ALL_ZERO(vcmp); /* softenings are not zero */
		if(msk) {
		    t1 = SIMD_MUL(SIMD_ADD(pmass,tile->m.p[j]),SIMD_MUL(p4soft2,tile->fourh2.p[j]));
		    t2 = SIMD_ADD(SIMD_MUL(tile->fourh2.p[j],pmass),SIMD_MUL(p4soft2,tile->m.p[j]));
#if defined(__SSE2__) || defined(__ALTIVEC__)
		    pfourh2 = SIMD_RE_EXACT(t2);
		    pfourh2 = SIMD_MUL(pfourh2,t1);
#else
		    pfourh2 = SIMD_DIV(t1,t2);
#endif
		    vcmp = SIMD_CMP_LT(tile->d2.p[j],pfourh2);
		    pd2 = SIMD_MAX(tile->d2.p[j],pfourh2);
		    msk = SIMD_ALL_ZERO(vcmp);  /* zero means nothing is softened - optimization */
		} else {
		    pd2 = tile->d2.p[j];
		}

		pir = SIMD_RSQRT_EXACT(pd2);
		pir2 = SIMD_MUL(pir,pir);
		pd2 = SIMD_MUL(pir2,tile->d2.p[j]); /* for SOFTENED */
		pir2 = SIMD_MUL(pir2,pir);

		/* pir and pir2 are valid now for both softened and unsoftened particles */
		/* Now we apply the fix to softened particles only */
		if (msk) {
		    pd2 = SIMD_SUB(consts.one.p,pd2);
		    pd2 = SIMD_AND(vcmp,pd2);
		    t1 = SIMD_MADD(consts.R45_32.p, pd2, consts.R3_8.p);
		    t1 = SIMD_MADD(t1, pd2, consts.half.p);
		    t1 = SIMD_MADD(t1, pd2, consts.one.p);
		    t2 = SIMD_MADD(consts.R135_16.p, pd2, consts.threehalves.p);
		    t2 = SIMD_MADD(t2, pd2, consts.one.p);
		    pir = SIMD_MUL(pir,t1);
		    pir2 = SIMD_MUL(pir2,t2);
		    }
		pir2 = SIMD_MUL(pir2,tile->m.p[j]);

		t1 = SIMD_NMSUB(tile->dx.p[j],pir2,consts.zero.p);
		t2 = SIMD_NMSUB(tile->dy.p[j],pir2,consts.zero.p);
		t3 = SIMD_NMSUB(tile->dz.p[j],pir2,consts.zero.p);
		ppot = SIMD_NMSUB(tile->m.p[j],pir,ppot);

		/* Time stepping criteria stuff */
		padotai = SIMD_MADD(piaz,t3,SIMD_MADD(piay,t2,SIMD_MUL(piax,t1)));
		vcmp = SIMD_AND(SIMD_CMP_GT(padotai,consts.zero.p),SIMD_CMP_GE(tile->d2.p[j],psmooth2));
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
	ILP_LOOP(ilp,tile) {
	    for (j=0;j<tile->nPart;++j) {
		d2 = tile->d2.f[j];
		fourh2 = softmassweight(fMass,4*fSoft*fSoft,
					tile->m.f[j],tile->fourh2.f[j]);
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
		    ++nSoft;
		    }
		dir2 *= tile->m.f[j];
		tax = -tile->dx.f[j]*dir2;
		tay = -tile->dy.f[j]*dir2;
		taz = -tile->dz.f[j]*dir2;
		fPot -= tile->m.f[j]*dir;
		/*
		** Calculations for determining the timestep.
		*/
		adotai = a[0]*tax + a[1]*tay + a[2]*taz;
		if (adotai > 0 && tile->d2.f[j] >= fsmooth2) {
		    adotai *= dimaga;
		    dirsum += dir*adotai*adotai;
		    normsum += adotai*adotai;
		    }
		ax += tax;
		ay += tay;
		az += taz;
		}
	    } /* end of particle list gravity loop */
#endif
	/*
	** Finally set new acceleration and potential.
	** Note that after this point we cannot use the new timestepping criterion since we
	** overwrite the acceleration.
	*/
	*pPot = fPot;
	a[0] = ax;
	a[1] = ay;
	a[2] = az;
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
		fx = a[0];
		fy = a[1];
		fz = a[2];
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
		  if (p->iOrder < pkd->param.nPartColl || tile->iOrder.i[j] < pkd->param.nPartColl) {
		    d2 = tile->d2.f[j];
		    fourh2 = softmassweight(fMass,4*fSoft*fSoft,
					tile->m.f[j],tile->fourh2.f[j]);
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
		    summ = fMass+tile->m.f[j];
		    rhopmaxlocal = summ*dir2;
		    vx = v[0] - tile->vx.f[j];
		    vy = v[1] - tile->vy.f[j];
		    vz = v[2] - tile->vz.f[j];
		    rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,tile->dx.f[j],tile->dy.f[j],tile->dz.f[j],vx,vy,vz,pkd->param.dEccFacMax);
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
	ilpRestore(ilp,&checkPt);
	} /* end of i-loop cells & particles */
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
