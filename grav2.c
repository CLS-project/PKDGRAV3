#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *grav2_module_id = "$Id$";

#include <math.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include "pkd.h"
#include "moments.h"
#include "meval.h"
#include "qeval.h"
#include "ewald.h"
#include "grav.h"


#ifdef USE_SIMD
static const struct CONSTS {
    v4sf zero;
    v4sf half;
    v4sf one;
    v4sf threehalves;
    v4sf three;
    v4sf four;
    v4sf R3_8;
    v4sf R45_32;
    v4sf R135_16;
    } consts = {
	{0.0,     0.0,     0.0,     0.0},
    {0.5,     0.5,     0.5,     0.5},
    {1.0,     1.0,     1.0,     1.0},
    {1.5,     1.5,     1.5,     1.5},
    {3.0,     3.0,     3.0,     3.0},
    {4.0,     4.0,     4.0,     4.0},
    {3.0/8.0, 3.0/8.0, 3.0/8.0, 3.0/8.0},
    {45.0/32.0, 45.0/32.0, 45.0/32.0, 45.0/32.0},
    {135.0/16.0, 135.0/16.0, 135.0/16.0, 135.0/16.0},
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
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,KDN *pBucket,LOCR *pLoc,ILP ilp,ILC ilc,
		    double dirLsum,double normLsum,int bEwald,double *pdFlop,double *pdEwFlop,double dRhoFac) {
    PARTICLE *p,*pj;
    KDN *pkdn = pBucket;
    const double onethird = 1.0/3.0;
    float *a, *pPot;
    momFloat ax,ay,az,fPot;
    double x,y,z,d2,dir,dir2;
    FLOAT fMass,fSoft;
    FLOAT fMassTmp,fSoftTmp;
    float fx, fy, fz;
    double dtGrav;
    momFloat adotai,maga,dimaga,dirsum,normsum;
    momFloat tax,tay,taz,tmon;
    double rholoc,dirDTS,dsmooth2,fSoftMedian,fEps,fEps2;
#ifndef USE_SIMD_MOMR
    double g2,g3,g4;
    double xx,xy,xz,yy,yz,zz;
    double xxx,xxz,yyy,yyz,xxy,xyy,xyz;
#else
    int nCellILC;
#endif
    double tx,ty,tz;
    ILPTILE tile;
    ILCTILE ctile;
    int nPartX;
    int i,j,nSP,nSoft,nActive;
    float rMax;
#if defined(USE_SIMD_PP)
    v4sf t1, t2, t3;
    v4sf pax, pay, paz;
    v4sf piax, piay, piaz;
    v4sf ppot, pmass, p4soft2;
    v4sf padotai,pimaga,psmooth2,pirsum,pnorms;
#else
    double fourh2;
#endif
    ILPCHECKPT checkPt;

    assert(pkd->oPotential);
    assert(pkd->oAcceleration);

#ifdef USE_SIMD_MOMR
    nCellILC = nCell;
    momPadSIMDMomr( &nCellILC, ilc );
#endif
    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;
    //rMax = 0.0;
    rMax = pkd->param.nPartRhoLoc/(pkdn->pUpper-pkdn->pLower+1)
	   *(pkdn->bnd.fMax[0]*pkdn->bnd.fMax[0]
	     + pkdn->bnd.fMax[1]*pkdn->bnd.fMax[1]
	     + pkdn->bnd.fMax[2]*pkdn->bnd.fMax[2]);
    /*
    ** Save the ilp list so that we can restore it for each particle. We will be 
    ** adding the local bucket interactions within this loop, and so modifying the 
    ** ilp.
    */
    ilpCheckPt(ilp,&checkPt);
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;
	pPot = pkdPot(pkd,p);
	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	a = pkdAccel(pkd,p);
	++nActive;
	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
	dsmooth2 = 0;
	tx = a[0];
	ty = a[1];
	tz = a[2];
	dimaga = tx*tx + ty*ty + tz*tz;
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
	    x = p->r[0] - pkdn->r[0];
	    y = p->r[1] - pkdn->r[1];
	    z = p->r[2] - pkdn->r[2];
	    momEvalLocr(pLoc,x,y,z,&fPot,&ax,&ay,&az);
	    }

#ifdef USE_SIMD_MOMR
	momEvalSIMDMomr( nCellILC, ilc, p->r, p->a,
			 &ax, &ay, &az, &fPot, &dirsum, &normsum );
#else
	ILC_LOOP(ilc,ctile) {
	    for (j=0;j<ctile->nCell;++j) {
		x = p->r[0] - ctile->d[j].x.f;
		y = p->r[1] - ctile->d[j].y.f;
		z = p->r[2] - ctile->d[j].z.f;
		d2 = x*x + y*y + z*z;
		SQRT1(d2,dir);
		dirDTS = dir;
		dir2 = dir*dir;
		g2 = 3*dir*dir2*dir2;
		g3 = 5*g2*dir2;
		g4 = 7*g3*dir2;
		/*
		** Calculate the funky distance terms.
		*/
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
		tx = g4*(ctile->d[j].xxxx.f*xxx + ctile->d[j].xyyy.f*yyy + ctile->d[j].xxxy.f*xxy + ctile->d[j].xxxz.f*xxz + ctile->d[j].xxyy.f*xyy + ctile->d[j].xxyz.f*xyz + ctile->d[j].xyyz.f*yyz);
		ty = g4*(ctile->d[j].xyyy.f*xyy + ctile->d[j].xxxy.f*xxx + ctile->d[j].yyyy.f*yyy + ctile->d[j].yyyz.f*yyz + ctile->d[j].xxyy.f*xxy + ctile->d[j].xxyz.f*xxz + ctile->d[j].xyyz.f*xyz);
		tz = g4*(-ctile->d[j].xxxx.f*xxz - (ctile->d[j].xyyy.f + ctile->d[j].xxxy.f)*xyz - ctile->d[j].yyyy.f*yyz + ctile->d[j].xxxz.f*xxx + ctile->d[j].yyyz.f*yyy - ctile->d[j].xxyy.f*(xxz + yyz) + ctile->d[j].xxyz.f*xxy + ctile->d[j].xyyz.f*xyy);
		g4 = 0.25*(tx*x + ty*y + tz*z);
		xxx = g3*(ctile->d[j].xxx.f*xx + ctile->d[j].xyy.f*yy + ctile->d[j].xxy.f*xy + ctile->d[j].xxz.f*xz + ctile->d[j].xyz.f*yz);
		xxy = g3*(ctile->d[j].xyy.f*xy + ctile->d[j].xxy.f*xx + ctile->d[j].yyy.f*yy + ctile->d[j].yyz.f*yz + ctile->d[j].xyz.f*xz);
		xxz = g3*(-(ctile->d[j].xxx.f + ctile->d[j].xyy.f)*xz - (ctile->d[j].xxy.f + ctile->d[j].yyy.f)*yz + ctile->d[j].xxz.f*xx + ctile->d[j].yyz.f*yy + ctile->d[j].xyz.f*xy);
		g3 = onethird*(xxx*x + xxy*y + xxz*z);
		xx = g2*(ctile->d[j].xx.f*x + ctile->d[j].xy.f*y + ctile->d[j].xz.f*z);
		xy = g2*(ctile->d[j].yy.f*y + ctile->d[j].xy.f*x + ctile->d[j].yz.f*z);
		xz = g2*(-(ctile->d[j].xx.f + ctile->d[j].yy.f)*z + ctile->d[j].xz.f*x + ctile->d[j].yz.f*y);
		g2 = 0.5*(xx*x + xy*y + xz*z);
		tmon = ctile->d[j].m.f*dir;
		dir2 *= tmon + 5*g2 + 7*g3 + 9*g4;
		fPot -= tmon + g2 + g3 + g4;
		tax = xx + xxx + tx - x*dir2;
		tay = xy + xxy + ty - y*dir2;
		taz = xz + xxz + tz - z*dir2;
		adotai = a[0]*tax + a[1]*tay + a[2]*taz;
		if (adotai > 0) {
		    adotai *= dimaga;
		    dirsum += dirDTS*adotai*adotai;
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
	    if (pkdIsSrcActive(pj,0,MAX_RUNG)) {
		fMassTmp = pkdMass(pkd,pj);
		fSoftTmp = pkdSoft(pkd,pj);
		ilpAppend(ilp,pj->r[0],pj->r[1],pj->r[2],
		    fMassTmp, 4*fSoftTmp*fSoftTmp,
		    pj->iOrder, pj->v[0], pj->v[1], pj->v[2]);
		}
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
	    ** Calculate local density only in the case of more than 1 neighbouring particle!
	    */
	    rholoc = 0.0;
	    nPartX = ilpCount(ilp);
	    if (nPartX > 1) {
		nSP = (nPartX < pkd->param.nPartRhoLoc)?nPartX:pkd->param.nPartRhoLoc;
		dsmooth2 = ilpSelect(ilp,nSP,&rMax);
#ifdef USE_SIMD
		psmooth2 = SIMD_SPLAT(dsmooth2);
#endif
                fSoftMedian = ilpSelectMass(ilp,nSP/2+1,nSP);
		SQRT1(dsmooth2,dir);
		dir2 = dir * dir;
		for (j=0;j<nSP;++j) {
		    /*
		    ** We keep this test for masses above zero for safety. 
		    ** Tracer particles should be implemented by having bSrcActive = 0!
		    */
		    assert(ilp->first->s.m.f[j] > 0.0);
		    if (ilp->first->s.fourh2.f[j] > 0.0) {
			fEps = fSoftMedian/ilp->first->s.fourh2.f[j];
			if (fEps > 1.0) fEps = 1.0;
			}
		    else fEps = 1.0;
		    fEps2 = fEps*fEps;
		    d2 = ilp->first->s.d2.f[j]*dir2*fEps2;
		    d2 = (1-d2);
		    if (d2 < 0) d2 = 0.0;
		    rholoc += d2*ilp->first->s.m.f[j]*fEps2*fEps;
		    }
		rholoc = 1.875*M_1_PI*rholoc*dir2*dir; /* F1 Kernel (15/8) */
		}
	    assert(rholoc >= 0.0);
	    }
	else {
	    dsmooth2 = 0.0;
#ifdef USE_SIMD_PP
	    psmooth2 = consts.zero;
#endif	    
	    }
	

#ifdef USE_SIMD_PP
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

		vcmp = SIMD_CMP_GT(tile->d.fourh2.p[j],consts.zero);
		msk = SIMD_ALL_ZERO(vcmp); /* softenings are not zero */
		if(msk) {
		    t1 = SIMD_MUL(SIMD_ADD(pmass,tile->d.m.p[j]),SIMD_MUL(p4soft2,tile->d.fourh2.p[j]));
		    t2 = SIMD_ADD(SIMD_MUL(tile->d.fourh2.p[j],pmass),SIMD_MUL(p4soft2,tile->d.m.p[j]));
#if defined(__SSE2__) || defined(__ALTIVEC__)
		    pfourh2 = SIMD_RE_EXACT(t2);
		    pfourh2 = SIMD_MUL(pfourh2,t1);
#else
		    pfourh2 = SIMD_DIV(t1,t2);
#endif
		    vcmp = SIMD_CMP_LT(tile->d.d2.p[j],pfourh2);
		    pd2 = SIMD_MAX(tile->d.d2.p[j],pfourh2);
		    msk = SIMD_ALL_ZERO(vcmp);  /* zero means nothing is softened - optimization */
		} else {
		    pd2 = tile->d.d2.p[j];
		}

		pir = SIMD_RSQRT_EXACT(pd2);
		pir2 = SIMD_MUL(pir,pir);
		pd2 = SIMD_MUL(pir2,tile->d.d2.p[j]); /* for SOFTENED */
		pir2 = SIMD_MUL(pir2,pir);

		/* pir and pir2 are valid now for both softened and unsoftened particles */
		/* Now we apply the fix to softened particles only */
		if (msk) {
		    pd2 = SIMD_SUB(consts.one,pd2);
		    pd2 = SIMD_AND(vcmp,pd2);
		    t1 = SIMD_MADD(consts.R45_32, pd2, consts.R3_8);
		    t1 = SIMD_MADD(t1, pd2, consts.half);
		    t1 = SIMD_MADD(t1, pd2, consts.one);
		    t2 = SIMD_MADD(consts.R135_16, pd2, consts.threehalves);
		    t2 = SIMD_MADD(t2, pd2, consts.one);
		    pir = SIMD_MUL(pir,t1);
		    pir2 = SIMD_MUL(pir2,t2);
		    }
		pir2 = SIMD_MUL(pir2,tile->d.m.p[j]);

		t1 = SIMD_NMSUB(tile->d.dx.p[j],pir2,consts.zero);
		t2 = SIMD_NMSUB(tile->d.dy.p[j],pir2,consts.zero);
		t3 = SIMD_NMSUB(tile->d.dz.p[j],pir2,consts.zero);

		/* Time stepping criteria stuff */
		padotai = SIMD_MADD(piaz,t3,SIMD_MADD(piay,t2,SIMD_MUL(piax,t1)));
		vcmp = SIMD_AND(SIMD_CMP_GT(padotai,consts.zero),SIMD_CMP_GE(tile->d.d2.p[j],psmooth2));
		padotai= SIMD_AND(padotai,vcmp);
		padotai= SIMD_MUL(padotai,pimaga);
		pd2 = SIMD_MUL(padotai,padotai);
		pirsum = SIMD_MADD(pir,pd2,pirsum);
		pnorms = SIMD_ADD(pnorms,pd2);

		ppot = SIMD_NMSUB(tile->d.m.p[j],pir,ppot);
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
	ILP_LOOP(ilp,tile) {
	    for (j=0;j<tile->nPart;++j) {
		double d2DTS;
		d2 = tile->d.d2.f[j];
		d2DTS = d2;
		fourh2 = softmassweight(fMass,4*fSoft*fSoft,
					tile->d.m.f[j],tile->d.fourh2.f[j]);

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

		dir2 *= tile->d.m.f[j];
		tax = -tile->d.dx.f[j]*dir2;
		tay = -tile->d.dy.f[j]*dir2;
		taz = -tile->d.dz.f[j]*dir2;
		adotai = a[0]*tax + a[1]*tay + a[2]*taz;
		if (adotai > 0 && d2DTS >= dsmooth2) {
		    adotai *= dimaga;
		    dirsum += dir*adotai*adotai;
		    normsum += adotai*adotai;
		    }
		fPot -= tile->d.m.f[j]*dir;
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
	    double dT;

	    /*
	    ** If this is the first time through, the accelerations will have 
	    ** all been zero resulting in zero for normsum (and nan for dtGrav).
	    ** We repeat this process again, so dtGrav will be correct.
	    */
	    if ( normsum > 0.0 ) {
		/*
		** Use new acceleration here!
		*/
		tx = a[0];
		ty = a[1];
		tz = a[2];
		maga = sqrt(tx*tx + ty*ty + tz*tz);
		dtGrav = maga*dirsum/normsum;
		}
	    else dtGrav = 0.0;
	    dtGrav += pkd->param.dPreFacRhoLoc*rholoc;
	    if ( dtGrav > 0.0 ) {
		dT = pkd->param.dEta/sqrt(dtGrav*dRhoFac);
		p->uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
		}
	    else p->uNewRung = 0;
	    p->fDensity = rholoc;
	    }
	/*
	** Restore the ilp to the same state that GravInteract was called.
	*/
	ilpRestore(ilp,&checkPt);
	} /* end of i-loop cells & particles */

    *pdFlop += nActive*(ilpCount(pkd->ilp)*40 + ilcCount(pkd->ilc)*200) + nSoft*15;
    return(nActive);
    }
