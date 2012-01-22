#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef PROFILE_GRAVWALK
#include "VtuneApi.h"
#endif
#ifdef __SSE__
#include <xmmintrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <string.h>
#include "pkd.h"
#include "walk.h"
#include "grav.h"
#include "smooth.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "moments.h"
#include "cl.h"

#ifdef USE_DEHNEN_THETA
double brent(double x1, double x2, double (*my_f)(double v, void *params),void *params) {
    static double EPS = 1e-12;
    int iter;
    double a=x1, b=x2, c, d, e, min1, min2;
    double fa, fb, fc, p, q, r, s, tol1, xm;

    fa = my_f(x1,params);
    fb = my_f(x2,params);

    if ( fb*fa > 0.0) abort();
    fc = fb;
    c = d = e = 0.0;
    for(iter=1;iter<=100;iter++) {
	if (fb*fc > 0.0) {
	    c = a;
	    fc = fa;
	    e = d = b-a;
	    }
	if (fabs(fc) < fabs(fb)) {
	    a = b;
	    b = c;
	    c = a;
	    fa = fb;
	    fb = fc;
	    fc = fa;
	    }
	tol1 = 2.0*EPS*fabs(b)+0.5*1e-9;
	xm = 0.5*(c-b);
	if (fabs(xm) <= tol1 || fb == 0.0) return b;
	if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
	    s = fb/fa;
	    if (a==c) {
		p = 2.0*xm*s;
		q = 1.0-s;
		}
	    else {
		q = fa/fc;
		r = fb/fc;
		p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
		q = (q-1.0)*(r-1.0)*(s-1.0);
		}
	    if (p>0.0) q = -q;
	    p = fabs(p);
	    min1 = 3.0*xm*q - fabs(tol1*q);
	    min2 = fabs(e*q);
	    if (2.0*p < (min1 < min2 ? min1 : min2)) {
		e = d;
		d = p / q;
		}
	    else {
		d = xm;
		e = d;
		}
	    }
	else {
	    d = xm;
	    e = d;
	    }
	a = b;
	fa = fb;
	if (fabs(d) > tol1)
	    b += d;
	else
	    b += (xm > 0.0 ? fabs(tol1) : - fabs(tol1));
	fb = my_f(b,params);
	}
    abort();
    return 0.0;
    }

struct theta_params { double alpha; };
static double theta_function(double x, void * params) {
    struct theta_params *p = params;
    double x2,x4,ixm,ixm2;
    x2 = x*x;
    x4 = x2*x2;
    ixm = 1.0 / (1-x);
    ixm2 = ixm*ixm;
    return p->alpha - x4*x2*x * ixm2;
    }

void pkdSetThetaTable(PKD pkd,double dThetaMin,double dThetaMax) {
    struct theta_params params;
    double dCritRatio;
    double dMass;
    int i;

    if (dThetaMax!=pkd->dCritThetaMax) {
	pkd->dCritThetaMax = dThetaMax;
	pkd->dCritThetaMin = dThetaMin;
	if (pkd->fCritMass==NULL) {
	    pkd->nCritBins = 128;
	    pkd->fCritMass = malloc(sizeof(float) * pkd->nCritBins);
	    assert(pkd->fCritMass!=NULL);
	    pkd->fCritTheta = malloc(sizeof(float) * pkd->nCritBins);
	    assert(pkd->fCritTheta!=NULL);
	    }

	dCritRatio = pow(pkd->dCritThetaMin,7.0) / pow(1-pkd->dCritThetaMin,2.0);
	for( i=0; i< pkd->nCritBins; i++ ) {
	    pkd->fCritMass[i] = dMass = exp( -0.5*(pkd->nCritBins-i-1) );
	    params.alpha = dCritRatio * pow(dMass,-1.0/3.0);
	    pkd->fCritTheta[i] = brent(0.0,0.9999999, theta_function, &params);
	    pkd->fCritTheta[i] = (pkd->fCritTheta[i]-pkd->dCritThetaMin)
		/ (1-pkd->dCritThetaMin)
		* (pkd->dCritThetaMax-pkd->dCritThetaMin) + pkd->dCritThetaMin;
	    }
	}
    }

/*
** Estimate for log(1+x)
*/
static inline float lge(float x) {
    float x2, x4;
    x2 = x*x;
    x4 = x2*x2;
    return x - x2/2.0 + x2*x/3.0 - x4/4.0 + x4*x/5.0 - x4*x2/6.0;
    }

static inline float getTheta(PKD pkd,float fMass) {
    int i,j,m;
    float fInterval;

    if (fMass<=pkd->fCritMass[0]) return pkd->dCritThetaMax;
    if (fMass>=1.0) return pkd->dCritThetaMin;
    assert(fMass>0.0 && fMass <= 1.0);
    assert(fMass>=pkd->fCritMass[0]);

    /* Locate the correct mass bin - binary search is better than ln? */
    i=0;
    j=pkd->nCritBins-1;
    while(i<j) {
	m = (j-i)/2 + i;
	if (pkd->fCritMass[m] > fMass) j = m;
	else i = m+1;
	}
    --i;

    assert(i>=0 && i<pkd->nCritBins);
    assert(pkd->fCritMass[i] <= fMass );
    assert(pkd->fCritMass[i+1] >= fMass );

    fInterval = lge((fMass-pkd->fCritMass[i])/pkd->fCritMass[i])*2.0;
    assert(fInterval>=0.0 && fInterval <= 1.0);
    return fInterval * (pkd->fCritTheta[i+1]-pkd->fCritTheta[i]) + pkd->fCritTheta[i];
    }
#endif

static inline int getCell(PKD pkd,int iCell,int id,float *pcOpen,KDN **pc) {
    KDN *c;
    int nc;
    assert(iCell > 0);
    if (id == pkd->idSelf) {
	*pc = c = pkdTreeNode(pkd,iCell);
	nc = c->pUpper - c->pLower + 1;
	}
    else if (id < 0) {
        *pc = c = pkdTopNode(pkd,iCell);
	nc = 1000000000; /* we never allow pp with this cell */
	}
    else {
	*pc = c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCell,id));
	nc = c->pUpper - c->pLower + 1;
	}
    if (*pcOpen < 0.0f) {
#ifdef USE_DEHNEN_THETA
	*pcOpen = c->bMax/getTheta(pkd,pkdNodeMom(pkd,c)->m/pkdNodeMom(pkd,pkdTreeNode(pkd,ROOT))->m);
#else
	*pcOpen = c->bMax * pkd->fiCritTheta;
#endif
	}
    return nc;
    }

#ifdef USE_SIMD_OPEN
static const struct CONSTS {
    v4 zero;
    v4 one;
    v4 threehalves;
    v4 two;
    v4 fMonopoleThetaFac2;
    } consts = {
        {SIMD_CONST(0.0)},
	{SIMD_CONST(1.0)},
	{SIMD_CONST(1.5)},
	{SIMD_CONST(2.0)},
	{SIMD_CONST(1.6f*1.6f)},
};
static const struct ICONSTS {
    i4 zero;
    i4 one;
    i4 two;
    i4 three;
    i4 four;
    i4 five;
    i4 six;
    i4 seven;
    i4 eight;
    /* no nine */
    i4 ten;
    i4 sixtyfour;
    i4 walk_min_multipole;
    } iconsts = {
        {SIMD_CONST(0)},
	{SIMD_CONST(1)},
	{SIMD_CONST(2)},
	{SIMD_CONST(3)},
	{SIMD_CONST(4)},
	{SIMD_CONST(5)},
	{SIMD_CONST(6)},
	{SIMD_CONST(7)},
	{SIMD_CONST(8)},
	{SIMD_CONST(10)},
	{SIMD_CONST(PKD_GROUP_SIZE)},
	{SIMD_CONST(3)},
};

static union {
    uint32_t u[SIMD_WIDTH];
    v4sf p;
    } const_fabs = {SIMD_CONST(0x7fffffff)};

/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version will also open buckets ("new" criteria)
*/
static void iOpenOutcomeSIMD(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin) {
    v4sf T0,T1,T2,T3,T4,T5,T6,T7,P1,P2,P3,P4;
    v4sf T,xc,yc,zc,dx,dy,dz,d2,diCrit,cOpen,cOpen2,d2Open,mink2,minbnd2,fourh2;
    int i,n,iEnd,nLeft;
    CL_BLK *blk;
    v4sf iOpen,iOpenA,iOpenB;
    pBND kbnd;
    v4sf k_xCenter, k_yCenter, k_zCenter, k_xMax, k_yMax, k_zMax;
    v4sf k_xMinBnd, k_yMinBnd, k_zMinBnd, k_xMaxBnd, k_yMaxBnd, k_zMaxBnd;
    v4sf k_x, k_y, k_z, k_m, k_4h2, k_bMax, k_Open;
    v4i  k_iLower, k_nk;

    assert ( pkdNodeMom(pkd,k)->m > 0.0f );

    diCrit = SIMD_SPLAT(1.0f/dThetaMin);

    pkdNodeBnd(pkd,k,&kbnd);
    k_xMinBnd = SIMD_SPLAT(kbnd.fCenter[0]-kbnd.fMax[0]);
    k_yMinBnd = SIMD_SPLAT(kbnd.fCenter[1]-kbnd.fMax[1]);
    k_zMinBnd = SIMD_SPLAT(kbnd.fCenter[2]-kbnd.fMax[2]);
    k_xMaxBnd = SIMD_SPLAT(kbnd.fCenter[0]+kbnd.fMax[0]);
    k_yMaxBnd = SIMD_SPLAT(kbnd.fCenter[1]+kbnd.fMax[1]);
    k_zMaxBnd = SIMD_SPLAT(kbnd.fCenter[2]+kbnd.fMax[2]);
    k_xCenter = SIMD_SPLAT(kbnd.fCenter[0]);
    k_yCenter = SIMD_SPLAT(kbnd.fCenter[1]);
    k_zCenter = SIMD_SPLAT(kbnd.fCenter[2]);
    k_xMax = SIMD_SPLAT(kbnd.fMax[0]);
    k_yMax = SIMD_SPLAT(kbnd.fMax[1]);
    k_zMax = SIMD_SPLAT(kbnd.fMax[2]);
    k_x = SIMD_SPLAT(k->r[0]);
    k_y = SIMD_SPLAT(k->r[1]);
    k_z = SIMD_SPLAT(k->r[2]);
    k_m = SIMD_SPLAT(pkdNodeMom(pkd,k)->m);
    k_4h2 = SIMD_SPLAT(4.0f*k->fSoft2);
    k_bMax = SIMD_SPLAT(k->bMax);
    k_iLower = SIMD_SPLATI32(k->iLower);
    k_nk = SIMD_SPLATI32(k->pUpper-k->pLower+1);
    k_Open = SIMD_MUL(consts.threehalves.p,SIMD_MUL(k_bMax,diCrit));

    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	iEnd = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	iEnd = (iEnd+SIMD_MASK) >> SIMD_BITS;
	for(i=0; i<iEnd; ++i) {
	    T = SIMD_OR(SIMD_CMP_GT(blk->m.p[i],consts.zero.p),SIMD_CMP_GT(k_4h2,consts.zero.p));
	    fourh2 = SIMD_OR(SIMD_AND(T,blk->fourh2.p[i]),SIMD_ANDNOT(T,consts.one.p));
	    fourh2 = SIMD_DIV(SIMD_MUL(SIMD_MUL(SIMD_ADD(k_m,blk->m.p[i]),k_4h2),fourh2),
		SIMD_ADD(SIMD_MUL(fourh2,k_m),SIMD_MUL(k_4h2,blk->m.p[i])));
	    xc = SIMD_ADD(blk->x.p[i],blk->xOffset.p[i]);
	    yc = SIMD_ADD(blk->y.p[i],blk->yOffset.p[i]);
	    zc = SIMD_ADD(blk->z.p[i],blk->zOffset.p[i]);
	    dx = SIMD_SUB(k_x,xc);
	    dy = SIMD_SUB(k_y,yc);
	    dz = SIMD_SUB(k_z,zc);
	    d2 = SIMD_MADD(dx,dx,SIMD_MADD(dy,dy,SIMD_MUL(dz,dz)));
	    cOpen = blk->cOpen.p[i];
	    cOpen2 = SIMD_MUL(cOpen,cOpen);
	    d2Open = SIMD_MUL(consts.two.p,SIMD_MAX(cOpen,k_Open));
	    d2Open = SIMD_MUL(d2Open,d2Open);
	    dx = SIMD_SUB(SIMD_AND(const_fabs.p,SIMD_SUB(xc,k_xCenter)),k_xMax);
	    dy = SIMD_SUB(SIMD_AND(const_fabs.p,SIMD_SUB(yc,k_yCenter)),k_yMax);
	    dz = SIMD_SUB(SIMD_AND(const_fabs.p,SIMD_SUB(zc,k_zCenter)),k_zMax);
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    dy = SIMD_AND(dy,SIMD_CMP_GT(dy,consts.zero.p));
	    dz = SIMD_AND(dz,SIMD_CMP_GT(dz,consts.zero.p));
	    mink2 = SIMD_MADD(dx,dx,SIMD_MADD(dy,dy,SIMD_MUL(dz,dz)));

	    minbnd2 = consts.zero.p;

	    dx = SIMD_SUB(k_xMinBnd,SIMD_ADD(blk->xCenter.p[i],SIMD_ADD(blk->xOffset.p[i],blk->xMax.p[i])));
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    minbnd2 = SIMD_MADD(dx,dx,minbnd2);
	    dx = SIMD_SUB(SIMD_SUB(SIMD_ADD(blk->xCenter.p[i],blk->xOffset.p[i]),blk->xMax.p[i]),k_xMaxBnd);
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    minbnd2 = SIMD_MADD(dx,dx,minbnd2);

	    dx = SIMD_SUB(k_yMinBnd,SIMD_ADD(blk->yCenter.p[i],SIMD_ADD(blk->yOffset.p[i],blk->yMax.p[i])));
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    minbnd2 = SIMD_MADD(dx,dx,minbnd2);
	    dx = SIMD_SUB(SIMD_SUB(SIMD_ADD(blk->yCenter.p[i],blk->yOffset.p[i]),blk->yMax.p[i]),k_yMaxBnd);
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    minbnd2 = SIMD_MADD(dx,dx,minbnd2);

	    dx = SIMD_SUB(k_zMinBnd,SIMD_ADD(blk->zCenter.p[i],SIMD_ADD(blk->zOffset.p[i],blk->zMax.p[i])));
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    minbnd2 = SIMD_MADD(dx,dx,minbnd2);
	    dx = SIMD_SUB(SIMD_SUB(SIMD_ADD(blk->zCenter.p[i],blk->zOffset.p[i]),blk->zMax.p[i]),k_zMaxBnd);
	    dx = SIMD_AND(dx,SIMD_CMP_GT(dx,consts.zero.p));
	    minbnd2 = SIMD_MADD(dx,dx,minbnd2);

	    T0 = SIMD_CMP_GT(blk->m.p[i],consts.zero.p);
	    T1 = SIMD_AND(SIMD_CMP_GT(d2,d2Open),SIMD_CMP_GT(minbnd2,fourh2));
	    T2 = SIMD_I2F(SIMD_CMP_EQ_EPI32(blk->iLower.p[i],iconsts.zero.p));
	    T3 = SIMD_OR(SIMD_I2F(SIMD_CMP_GT_EPI32(iconsts.walk_min_multipole.p,blk->nc.p[i])),SIMD_CMP_LE(mink2,cOpen2));
	    T4 = SIMD_CMP_GT(minbnd2,fourh2);
	    T5 = SIMD_CMP_GT(mink2,SIMD_MUL(consts.fMonopoleThetaFac2.p,cOpen2));
	    T6 = SIMD_CMP_GT(cOpen,k_Open);
	    /*invert:T7 = SIMD_I2F(SIMD_CMP_EQ_EPI32(k_iLower,iconsts.zero.p));*/
	    T7 = SIMD_I2F(SIMD_CMP_GT_EPI32(k_nk,iconsts.sixtyfour.p));
	    iOpenA = SIMD_OR(SIMD_AND(T2,iconsts.one.pf),SIMD_ANDNOT(T2,iconsts.three.pf));
	    iOpenB = SIMD_OR(SIMD_AND(T3,iOpenA),SIMD_ANDNOT(T3,
		    SIMD_OR(SIMD_AND(T4,iconsts.four.pf),SIMD_ANDNOT(T4,
			    SIMD_OR(SIMD_AND(T5,iconsts.five.pf),SIMD_ANDNOT(T5,iOpenA))))));
#ifdef NEW_OPENING_CRIT
	    P1 = SIMD_OR(SIMD_AND(T2,iconsts.two.pf),SIMD_ANDNOT(T2,iconsts.three.pf));
#else
	    P1 = SIMD_OR(SIMD_AND(T2,iOpenB),SIMD_ANDNOT(T2,iconsts.three.pf));
#endif
	    P2 = SIMD_OR(SIMD_AND(T7,iconsts.zero.pf),SIMD_ANDNOT(T7,iOpenB));
	    P3 = SIMD_OR(SIMD_AND(T6,P1),SIMD_ANDNOT(T6,P2));
	    P4 = SIMD_OR(SIMD_AND(T1,iconsts.eight.pf),SIMD_ANDNOT(T1,P3));
	    iOpen = SIMD_OR(SIMD_AND(T0,P4),SIMD_ANDNOT(T0,iconsts.ten.pf));
	    blk->iOpen.pf[i] = iOpen;
	    }
	}
    }

#endif
//#else

/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
*/
static void iOpenOutcomeOldCL(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin) {
    const float fMonopoleThetaFac2 = 1.6f * 1.6f;
    const int walk_min_multipole = 3;

    float dx,dy,dz,mink2,d2,d2Open,xc,yc,zc,fourh2,minbnd2,kOpen,cOpen,diCrit;
    int j,i,nk;
    int iOpen,iOpenA,iOpenB;
    CL_BLK *blk;
    int n, nLeft;
    pBND kbnd;

    pkdNodeBnd(pkd,k,&kbnd);
    nk = k->pUpper - k->pLower + 1;

    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	n = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	for(i=0; i<n; ++i) {
	    if (blk->m.f[i] <= 0) iOpen = 10;  /* ignore this cell */
	    else {
		fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,blk->m.f[i],blk->fourh2.f[i]);
		xc = blk->x.f[i] + blk->xOffset.f[i];
		yc = blk->y.f[i] + blk->yOffset.f[i];
		zc = blk->z.f[i] + blk->zOffset.f[i];
		d2 = pow(k->r[0]-xc,2) + pow(k->r[1]-yc,2) + pow(k->r[2]-zc,2);
		diCrit = 1.0 / dThetaMin;
		kOpen = 1.5*k->bMax*diCrit;
		cOpen = blk->cOpen.f[i];
		d2Open = pow(2.0*fmax(cOpen,kOpen),2);
		dx = fabs(xc - kbnd.fCenter[0]) - kbnd.fMax[0];
		dy = fabs(yc - kbnd.fCenter[1]) - kbnd.fMax[1];
		dz = fabs(zc - kbnd.fCenter[2]) - kbnd.fMax[2];
		mink2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);
		minbnd2 = 0;

		dx = kbnd.fCenter[0] - kbnd.fMax[0] -  blk->xCenter.f[i] - blk->xOffset.f[i] - blk->xMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->xCenter.f[i] + blk->xOffset.f[i] - blk->xMax.f[i] - kbnd.fCenter[0] - kbnd.fMax[0];
		if (dx > 0) minbnd2 += dx*dx;

		dx = kbnd.fCenter[1] - kbnd.fMax[1] - blk->yCenter.f[i] - blk->yOffset.f[i] - blk->yMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->yCenter.f[i] + blk->yOffset.f[i] - blk->yMax.f[i] - kbnd.fCenter[1] - kbnd.fMax[1];
		if (dx > 0) minbnd2 += dx*dx;

		dx = kbnd.fCenter[2] - kbnd.fMax[2] - blk->zCenter.f[i] - blk->zOffset.f[i] - blk->zMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->zCenter.f[i] + blk->zOffset.f[i] - blk->zMax.f[i] - kbnd.fCenter[2] - kbnd.fMax[2];
		if (dx > 0) minbnd2 += dx*dx;

		if (d2 > d2Open && minbnd2 > fourh2) iOpen = 8;
		else {
		    if (blk->iLower.i[i] == 0) iOpenA = 1;
		    else iOpenA = 3;
		    if (blk->nc.i[i] < walk_min_multipole || mink2 <= cOpen*cOpen) iOpenB = iOpenA;
		    else if (minbnd2 > fourh2) iOpenB = 4;
		    else if (mink2 > fMonopoleThetaFac2*cOpen*cOpen) iOpenB = 5;
		    else iOpenB = iOpenA;
		    if (cOpen > kOpen) {
			if (blk->iLower.i[i]) iOpen = 3;
			else iOpen = iOpenB;
			}
		    else {
/*		    if (k->iLower) iOpen = 0;*/
			if (nk>PKD_GROUP_SIZE) iOpen = 0;
			else iOpen = iOpenB;
			}
		    }
		}
	    blk->iOpen.i[i] = iOpen;
	    }
	}
    }

/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version has been changed by adding the ability to open buckets.
*/
static void iOpenOutcomeNewCL(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin) {
    const float fMonopoleThetaFac2 = 1.6f * 1.6f;
    const int walk_min_multipole = 3;

    float dx,dy,dz,mink2,d2,d2Open,xc,yc,zc,fourh2,minbnd2,kOpen,cOpen,diCrit;
    int j,i,nk;
    int iOpen,iOpenA,iOpenB;
    CL_BLK *blk;
    int n, nLeft;
    pBND kbnd;

    pkdNodeBnd(pkd,k,&kbnd);
    nk = k->pUpper - k->pLower + 1;

    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	n = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	for(i=0; i<n; ++i) {
	    if (blk->m.f[i] <= 0) iOpen = 10;  /* ignore this cell */
	    else {
		fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,blk->m.f[i],blk->fourh2.f[i]);
		xc = blk->x.f[i] + blk->xOffset.f[i];
		yc = blk->y.f[i] + blk->yOffset.f[i];
		zc = blk->z.f[i] + blk->zOffset.f[i];
		d2 = pow(k->r[0]-xc,2) + pow(k->r[1]-yc,2) + pow(k->r[2]-zc,2);
		diCrit = 1.0f / dThetaMin;
		kOpen = 1.5f * k->bMax * diCrit;
		cOpen = blk->cOpen.f[i];
		d2Open = pow(2.0f * fmax(cOpen,kOpen),2);
		dx = fabs(xc - kbnd.fCenter[0]) - kbnd.fMax[0];
		dy = fabs(yc - kbnd.fCenter[1]) - kbnd.fMax[1];
		dz = fabs(zc - kbnd.fCenter[2]) - kbnd.fMax[2];
		mink2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);
		minbnd2 = 0;

		dx = kbnd.fCenter[0] - kbnd.fMax[0] -  blk->xCenter.f[i] - blk->xOffset.f[i] - blk->xMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->xCenter.f[i] + blk->xOffset.f[i] - blk->xMax.f[i] - kbnd.fCenter[0] - kbnd.fMax[0];
		if (dx > 0) minbnd2 += dx*dx;

		dx = kbnd.fCenter[1] - kbnd.fMax[1] - blk->yCenter.f[i] - blk->yOffset.f[i] - blk->yMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->yCenter.f[i] + blk->yOffset.f[i] - blk->yMax.f[i] - kbnd.fCenter[1] - kbnd.fMax[1];
		if (dx > 0) minbnd2 += dx*dx;

		dx = kbnd.fCenter[2] - kbnd.fMax[2] - blk->zCenter.f[i] - blk->zOffset.f[i] - blk->zMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->zCenter.f[i] + blk->zOffset.f[i] - blk->zMax.f[i] - kbnd.fCenter[2] - kbnd.fMax[2];
		if (dx > 0) minbnd2 += dx*dx;

		if (d2 > d2Open && minbnd2 > fourh2) iOpen = 8;
		else if (cOpen > kOpen) {
		    if (blk->iLower.i[i]) iOpen = 3;
		    else iOpen = 2;
		    }
		else {
		    if (blk->iLower.i[i] == 0) iOpenA = 1;
		    else iOpenA = 3;
		    if (blk->nc.i[i] < walk_min_multipole || mink2 <= cOpen*cOpen) iOpenB = iOpenA;
		    else if (minbnd2 > fourh2) iOpenB = 4;
		    else if (mink2 > fMonopoleThetaFac2*cOpen*cOpen) iOpenB = 5;
		    else iOpenB = iOpenA;
		    if (nk>PKD_GROUP_SIZE) iOpen = 0;
		    else iOpen = iOpenB;
		    }
		}
#ifdef USE_SIMD_OPEN
	    /* This function is only called with SIMD if we are checking the two. Print the differences. */
	    if (blk->iOpen.i[i] != iOpen) {
		printf("SIMD=%d, iOpen=%d, d2=%.8g > d2Open=%.8g, minbnd2=%.8g > fourh2=%.8g, kOpen=%.8g, cOpen=%.8g\n",
		    blk->iOpen.i[i], iOpen, d2, d2Open, minbnd2, fourh2,
		    kOpen, blk->cOpen.f[i]);
		}
#else
	    blk->iOpen.i[i] = iOpen;
#endif
	    }
	}
    }

/*
** This doesn't work.
*/
static void iOpenOutcomeExperemental(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin) {
    const float fMonopoleThetaFac2 = 1.6f * 1.6f;
    const int walk_min_multipole = 3;

    float dx,dy,dz,mink2,minc2,d2,d2Open,xc,yc,zc,fourh2,minbnd2,kOpen,cOpen,diCrit;
    int j,i,nk;
    int iOpen,iOpenA,iOpenB, T1;
    CL_BLK *blk;
    int n, nLeft;
    pBND kbnd;

    pkdNodeBnd(pkd,k,&kbnd);
    nk = k->pUpper - k->pLower + 1;

    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	n = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	for(i=0; i<n; ++i) {
	    T1 = (blk->iCell.i[i] < 0);
	    if (blk->m.f[i] <= 0) iOpen = 10;  /* ignore this cell */
	    else {
		fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,blk->m.f[i],blk->fourh2.f[i]);
		diCrit = 1.0 / dThetaMin;
		kOpen = 1.5*k->bMax*diCrit;
		if (T1) {
		    dx = fabs(k->r[0] - blk->xCenter.f[i] - blk->xOffset.f[i]) - blk->xMax.f[i];
		    dy = fabs(k->r[1] - blk->yCenter.f[i] - blk->yOffset.f[i]) - blk->yMax.f[i];
		    dz = fabs(k->r[2] - blk->zCenter.f[i] - blk->zOffset.f[i]) - blk->zMax.f[i];
		    minc2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);	
		    if (k->iLower == 0) iOpenA = 1;
		    else iOpenA = 0;
		    if (minc2 <= kOpen*kOpen) iOpen = iOpenA;  /* this cell simply stays on the checklist */
		    else if (minc2 > fourh2) iOpen = 6;  /* Add particles of C to local expansion */
		    else if (minc2 > fMonopoleThetaFac2*kOpen*kOpen) iOpen = 7;
		    else iOpen = iOpenA;
		    }
		else {
		    xc = blk->x.f[i] + blk->xOffset.f[i];
		    yc = blk->y.f[i] + blk->yOffset.f[i];
		    zc = blk->z.f[i] + blk->zOffset.f[i];
		    d2 = pow(k->r[0]-xc,2) + pow(k->r[1]-yc,2) + pow(k->r[2]-zc,2);
		    cOpen = blk->cOpen.f[i];
		    d2Open = pow(2.0*fmax(cOpen,kOpen),2);
		    dx = fabs(xc - kbnd.fCenter[0]) - kbnd.fMax[0];
		    dy = fabs(yc - kbnd.fCenter[1]) - kbnd.fMax[1];
		    dz = fabs(zc - kbnd.fCenter[2]) - kbnd.fMax[2];
		    mink2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);
		    minbnd2 = 0;

		    dx = kbnd.fCenter[0] - kbnd.fMax[0] -  blk->xCenter.f[i] - blk->xOffset.f[i] - blk->xMax.f[i];
		    if (dx > 0) minbnd2 += dx*dx;
		    dx = blk->xCenter.f[i] + blk->xOffset.f[i] - blk->xMax.f[i] - kbnd.fCenter[0] - kbnd.fMax[0];
		    if (dx > 0) minbnd2 += dx*dx;

		    dx = kbnd.fCenter[1] - kbnd.fMax[1] - blk->yCenter.f[i] - blk->yOffset.f[i] - blk->yMax.f[i];
		    if (dx > 0) minbnd2 += dx*dx;
		    dx = blk->yCenter.f[i] + blk->yOffset.f[i] - blk->yMax.f[i] - kbnd.fCenter[1] - kbnd.fMax[1];
		    if (dx > 0) minbnd2 += dx*dx;

		    dx = kbnd.fCenter[2] - kbnd.fMax[2] - blk->zCenter.f[i] - blk->zOffset.f[i] - blk->zMax.f[i];
		    if (dx > 0) minbnd2 += dx*dx;
		    dx = blk->zCenter.f[i] + blk->zOffset.f[i] - blk->zMax.f[i] - kbnd.fCenter[2] - kbnd.fMax[2];
		    if (dx > 0) minbnd2 += dx*dx;

		    if (blk->iLower.i[i] == 0) iOpenA = 2; /* add this cell as an "opened" bucket */
		    else iOpenA = 3; /* open this cell and add its children to the checklist */
		    /*
		    ** iOpenB handles all cases where a C-C interaction is NOT acceptable.
		    */
		    if (cOpen > kOpen) iOpenB = iOpenA;
		    else if (k->iLower != 0) iOpenB = 0; /* keep this cell on the checklist */
		    else if (mink2 <= cOpen*cOpen) iOpenB = iOpenA;
		    else if (minbnd2 > fourh2) iOpenB = 4; /* add this cell to the P-C list */
		    else if (mink2 > fMonopoleThetaFac2*cOpen*cOpen) iOpenB = 5; /* use this cell as a softened monopole */
		    else iOpenB = iOpenA;
		    //else if (blk->iLower.i[i] == 0) iOpenB = 1; /* open the bucket P-P */
		    //else iOpenB = 3; /* open this cell and add its children to the checklist */

		    if (d2 <= d2Open) iOpen = iOpenB;
		    else if (minbnd2 > fourh2) iOpen = 8;
		    else if (d2 > fMonopoleThetaFac2*d2Open) iOpen = 9; /* it is absolutely critical to include this case for large softening */
		    else iOpen = iOpenB;
		    }
		}
	    blk->iOpen.i[i] = iOpen;
	    }
	}
    }
//#endif


#if 0
/*
** This implements the original pkdgrav Barnes Hut opening criterion.
*/
static inline int iOpenOutcomeBarnesHut(PKD pkd,KDN *k,CELT *check,KDN **pc,float dThetaMin) {
    const int walk_min_multipole = 3;
    double dMin,dMax,min2,max2,d2,fourh2,cOpen2;
    KDN *c;
    int iCell,nc;
    int iOpenA;
    int j;
    pBND kbnd;

    assert(check->iCell > 0);
    iCell = check->iCell;
    if (check->id == pkd->idSelf) {
	c = pkdTreeNode(pkd,iCell);
	nc = c->pUpper - c->pLower + 1;
	}
    else if (check->id < 0) {
        c = pkdTopNode(pkd,iCell);
	assert(c->iLower != 0);
	nc = walk_min_multipole-1; /* we never allow pp with this cell */
	}
    else {
	c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCell,check->id));
	nc = c->pUpper - c->pLower + 1;
	}

    pkdNodeBnd(pkd, k, &kbnd);
    *pc = c;

    cOpen2 = c->bMax * c->bMax / (dThetaMin*dThetaMin);

    if (c->iLower == 0) iOpenA = 1;
    else iOpenA = 3;

    if (pkdNodeMom(pkd,c)->m <= 0) return(10);  /* ignore this cell */
    else if (k->iLower) {
	/*
	** If this cell is not a bucket calculate the min distance
	** from the center of mass of the check cell to the bounding
	** box of the cell we are calculating gravity on. We also
	** calculate the size of the ball (from checkcell CM) which
	** just contains the this cell (radius^2 given by max2).
	*/
	min2 = 0;
	max2 = 0;
	for (j=0;j<3;++j) {
	    dMin = fabs(c->r[j] + check->rOffset[j] - kbnd.fCenter[j]);
	    dMax = dMin + kbnd.fMax[j];
	    dMin -= kbnd.fMax[j];
	    if (dMin > 0) min2 += dMin*dMin;
	    max2 += dMax*dMax;
	    }
	if (max2 <= cOpen2) return(iOpenA); /* open it for all particles of c */
	else if (min2 > cOpen2) {
	    /*
	    ** For symmetrized softening we need to calculate the distance between
	    ** center of mass between the two cells.
	    */
	    d2 = 0;
	    for (j=0;j<3;++j) {
		d2 += pow(c->r[j] + check->rOffset[j] - k->r[j],2);
		}
	    fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,pkdNodeMom(pkd,c)->m,4*c->fSoft2);
	    if (d2 > fourh2) {
		if (nc >= walk_min_multipole) return(4);  /* accept multipole */
		else return(iOpenA);  /* open the cell for performance reasons */
		}
	    else return(0);   /* in all other cases we can't decide until we get down to a bucket */
	    }
	else return(0);
	}
    else {
	/*
	** If this cell is a bucket we have to either open the checkcell
	** and get the particles, or accept the multipole. For this
	** reason we only need to calculate min2.
	*/
	min2 = 0;
	for (j=0;j<3;++j) {
	    dMin = fabs(c->r[j] + check->rOffset[j] - kbnd.fCenter[j]);
	    dMin -= kbnd.fMax[j];
	    if (dMin > 0) min2 += dMin*dMin;
	    }
	/*
	** By default we open the cell!
	*/
	if (min2 > cOpen2) {
	    /*
	    ** For symmetrized softening we need to calculate the distance between
	    ** center of mass between the two cells.
	    */
	    d2 = 0;
	    for (j=0;j<3;++j) {
		d2 += pow(c->r[j] + check->rOffset[j] - k->r[j],2);
		}
	    fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,pkdNodeMom(pkd,c)->m,4*c->fSoft2);
	    if (d2 > fourh2) {
		if (nc >= walk_min_multipole) return(4);  /* accept multipole */
		else return(iOpenA);  /* open the cell for performance reasons */
		}
	    else return(5); /* means we treat this cell as a softened monopole */
	    }
	else return(iOpenA);  /* open the cell */
	}
    }
#endif

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,
		int bVeryActive,double dThetaMin,double dThetaMax,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p;
    KDN *k,*c,*kFind;
    FMOMR *momc,*momk;
    FMOMR monoPole;
    FLOCR L;
    double fWeight = 0.0;
    double dShiftFlop;
    double dRhoFac;
    const double *v, *a;
    double dOffset[3];
    double xParent,yParent,zParent;
    double cx,cy,cz,d2c;
    double d2,fourh2;
    double dx[3],dir,dir2;
    float fOffset[3];
    float bMaxParent;
    float dirLsum,normLsum,adotai,maga;
    float tax,tay,taz;
    float fMass,fSoft;
    uint64_t iOrder;
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxInitCheck;
    int id,iCell,iSib,iLower,iCheckCell,iCheckLower,iCellDescend;
    int i,j,jTile,pi,pj,nActive,nTotActive;
    float cOpen,kOpen;
    pBND cbnd,kbnd;
    static const float  fZero3[] = {0.0f,0.0f,0.0f};
    static const double dZero3[] = {0.0,0.0,0.0};
    int nc,nk;
    ILPTILE tile;
    ILCTILE ctile;
    CL clTemp;
    CLTILE cltile;
#ifdef USE_SIMD_MOMR
    int ig,iv;
#endif
    SMX smx;
    SMF smf;
    double tempI;
    double dEwFlop = 0.0;

#ifdef PROFILE_GRAVWALK
    VTResume();
#endif

    monoPole.m = 0.0f;
    monoPole.xx = monoPole.yy = monoPole.xy = monoPole.xz = monoPole.yz = 0.0f;
    monoPole.xxx = monoPole.xyy = monoPole.xxy = monoPole.yyy = monoPole.xxz = monoPole.yyz = monoPole.xyz = 0.0f;
    monoPole.xxxx = monoPole.xyyy = monoPole.xxxy = monoPole.yyyy = monoPole.xxxz =
	monoPole.yyyz = monoPole.xxyy = monoPole.xxyz = monoPole.xyyz = 0.0f;

    /*
    ** If necessary, calculate the theta interpolation tables.
    */
#ifdef USE_DEHNEN_THETA
    pkdSetThetaTable(pkd,dThetaMin,dThetaMax);
#else
    pkd->fiCritTheta = 1.0f / dThetaMin;
#endif

    a = dZero3;
    v = dZero3;
    assert(pkd->oNodeMom);
    if (pkd->param.bGravStep) {
	assert(pkd->oNodeAcceleration);
	if (pkd->param.iTimeStepCrit == 1) {
	    assert(pkd->oNodeVelocity);
	    assert(pkd->oVelocity);
	}
    }

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    ** Otherwise we check that the ROOT has active particles!
    */
    if (bVeryActive) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkdTreeNode(pkd,VAROOT)->pUpper - pkdTreeNode(pkd,VAROOT)->pLower + 1);
	}
    else if (!pkdIsCellActive(pkdTreeNode(pkd,ROOT),uRungLo,uRungHi)) return 0;
    /*
    ** Initially we set our cell pointer to
    ** point to the top tree.
    */
    nTotActive = 0;
    ilpClear(pkd->ilp);
    ilcClear(pkd->ilc);
    clClear(pkd->cl);
    /*
    ** Setup smooth for calculating local densities when a particle has too few P-P interactions.
    */
    if (pkd->param.bGravStep) {
	smInitializeRO(&smx,pkd,&smf,pkd->param.nPartRhoLoc,nReps?1:0,SMX_DENSITY_F1);
	smSmoothInitialize(smx);
	/* No particles are inactive for density calculation */
	for (pi=0;pi<pkd->nLocal;++pi) {
	    p = pkdParticle(pkd,pi);
	    smx->ea[pi].bInactive = 0;
	    }
	}
    else smx = NULL;

    iStack = -1;
    /*
    ** Clear local expansion and the timestepping sums.
    */
    momClearFlocr(&L);
    dirLsum = 0;
    normLsum = 0;
    /*
    ** Precalculate RhoFac if required.
    */
    if (pkd->param.bGravStep) {
	double a = csmTime2Exp(pkd->param.csm,dTime);
	dRhoFac = 1.0/(a*a*a);
	}
    else dRhoFac = 0.0;
    /*
    ** First we add any replicas of the entire box
    ** to the Checklist.
    */
    for (ix=-nReps;ix<=nReps;++ix) {
	fOffset[0] = ix*pkd->fPeriod[0];
	for (iy=-nReps;iy<=nReps;++iy) {
	    fOffset[1] = iy*pkd->fPeriod[1];
	    for (iz=-nReps;iz<=nReps;++iz) {
		fOffset[2] = iz*pkd->fPeriod[2];
		bRep = ix || iy || iz;
		if (bRep || bVeryActive) {
		    /* 
		    ** Use leaf of the top tree and NOT the root of the local tree here.
		    */
		    cOpen = -1.0f;
		    id = -1;
		    iLower = pkdTopNode(pkd,ROOT)->iLower;
		    if (!iLower) iLower = ROOT;  /* something other than zero for openening crit - iLower usually can't be == ROOT */
		    nc = getCell(pkd,ROOT,id,&cOpen,&c);
		    pkdNodeBnd(pkd,c,&cbnd);
		    clAppend(pkd->cl,ROOT,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
		    }
		if (bRep && bVeryActive) {
		    /*
		    ** Add the images of the very active tree to the checklist.
		    */
		    cOpen = -1.0f;
		    id = mdlSelf(pkd->mdl);
		    iLower = pkdTopNode(pkd,VAROOT)->iLower;
		    if (!iLower) iLower = ROOT;  /* something other than zero for openening crit - iLower usually can't be == ROOT */
		    nc = getCell(pkd,VAROOT,id,&cOpen,&c);
		    pkdNodeBnd(pkd,c,&cbnd);
		    clAppend(pkd->cl,VAROOT,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
		    }
		}
	    }
	}
    if (!bVeryActive) {
	/*
	** This adds all siblings of a chain leading from the local tree leaf in the top
	** tree up to the ROOT of the top tree.
	*/
	for (j=0;j<3;++j) fOffset[j] = 0.0f;
	iCell = pkd->iTopRoot;
	iSib = SIBLING(iCell);
	while (iSib) {
	    cOpen = -1.0f;
	    id = -1;
	    iLower = pkdTopNode(pkd,iSib)->iLower;
	    if (!iLower) iLower = ROOT;  /* something other than zero for openening crit - iLower usually can't be == ROOT */
	    nc = getCell(pkd,iSib,id,&cOpen,&c);
	    pkdNodeBnd(pkd,c,&cbnd);
	    clAppend(pkd->cl,iSib,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
	    iCell = pkdTopNode(pkd,iCell)->iParent;
	    iSib = SIBLING(iCell);
	    }
	}
    /*
    ** We are now going to work on the local tree.
    ** Make iCell point to the root of the tree again.
    */
    if (bVeryActive) k = pkdTreeNode(pkd,iCell = VAROOT);
    else k = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
	/*
	** Find the next active particle that will be encountered in the walk algorithm below
	** in order to set a good initial center for the P-P interaction list.
	*/
	kFind = k;
	while (kFind->iLower) {
	    kFind = pkdTreeNode(pkd,iCellDescend = kFind->iLower);
	    if (!kFind->nActive) {
		/*
		** Move onto processing the sibling.
		*/
		kFind = pkdTreeNode(pkd,++iCellDescend);
	    }
	}
	for (pj=kFind->pLower;pj<=kFind->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    if (!pkdIsActive(pkd,p)) continue;
	    cx = p->r[0];
	    cy = p->r[1];
	    cz = p->r[2];
	    goto found_it;
	}
	assert(0); /* We didn't find an active particle */
    found_it:
	d2c = (cx - pkd->ilp->cx)*(cx - pkd->ilp->cx) + (cy - pkd->ilp->cy)*(cy - pkd->ilp->cy) +
	      (cz - pkd->ilp->cz)*(cz - pkd->ilp->cz);
	if ( d2c > 1e-5 ) {
	    /*
	    ** Correct all remaining PP interactions to this new center.
	    */
	    ILP_LOOP( pkd->ilp, tile ) {
		ILP_BLK *blk = tile->blk;
		int nLeft, n, prt;
		for( nLeft=tile->lstTile.nBlocks; nLeft >= 0; --nLeft,blk++ ) {
		    n = (nLeft ? pkd->ilp->lst.nPerBlock : tile->lstTile.nInLast);
		    for (prt=0; prt<n; ++prt) {
			blk->dx.f[prt] += cx - pkd->ilp->cx;
			blk->dy.f[prt] += cy - pkd->ilp->cy;
			blk->dz.f[prt] += cz - pkd->ilp->cz;
			}
		    }
		}
	    pkd->ilp->cx = cx;
	    pkd->ilp->cy = cy;
	    pkd->ilp->cz = cz;
	    /*
	    ** Correct all remaining PC interactions to this new center.
	    */
	    ILC_LOOP( pkd->ilc, ctile ) {
		ILC_BLK *blk = ctile->blk;
		int nLeft, n, prt;
		for( nLeft=ctile->lstTile.nBlocks; nLeft >= 0; --nLeft,blk++ ) {
		    n = (nLeft ? pkd->ilc->lst.nPerBlock : ctile->lstTile.nInLast);
		    for (prt=0; prt<n; ++prt) {
			blk->dx.f[prt] += cx - pkd->ilc->cx;
			blk->dy.f[prt] += cy - pkd->ilc->cy;
			blk->dz.f[prt] += cz - pkd->ilc->cz;
			}
		    }
		}
	    pkd->ilc->cx = cx;
	    pkd->ilc->cy = cy;
	    pkd->ilc->cz = cz;
	    }
	while (1) {
	    /*
	    ** Process the Checklist.
	    */
	    tempI = *pdFlop;
	    tempI += dEwFlop;
	    if (pkd->param.bGravStep) {
		a = pkdNodeAccel(pkd,k);
		maga = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	    }
	    /*
	    ** For cells which will remain on the checklist for a further deeper level of 
	    ** level of the tree we will be using the stack cl (even if no push happens later).
	    */
	    clClear(pkd->S[iStack+1].cl);
	    do {
		CL_LOOP(pkd->cl,cltile) {
#ifdef LOCAL_EXPANSION
#ifdef USE_SIMD_OPEN
		    iOpenOutcomeSIMD(pkd,k,pkd->cl,cltile,dThetaMin);
		    /*Verify:iOpenOutcomeNewCL(pkd,k,pkd->cl,cltile,dThetaMin);*/
#else
#ifdef NEW_OPENING_CRIT
		    iOpenOutcomeNewCL(pkd,k,pkd->cl,cltile,dThetaMin);
#else
		    iOpenOutcomeOldCL(pkd,k,pkd->cl,cltile,dThetaMin);
#endif
#endif
#else
		    assert(NULL);
#endif
		    }

		clClear(pkd->clNew);
		CL_LOOP(pkd->cl,cltile) {
		    CL_BLK *blk = cltile->blk;
		    int nLeft;
		    for(nLeft=cltile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
			int n = nLeft ? pkd->cl->lst.nPerBlock : cltile->lstTile.nInLast;
			for (jTile=0;jTile<n;++jTile) {
			    switch (blk->iOpen.i[jTile]) {
			    case 0:
				/*
				** This checkcell stays on the checklist.
				*/
				clAppendItem(pkd->S[iStack+1].cl,blk,jTile);
				break;
			    case 1:
				/*
				** This checkcell's particles are added to the P-P list.
				*/
				iCheckCell = blk->iCell.i[jTile];
				if (iCheckCell < 0) {
				    id = blk->id.i[jTile];
				    pj = -1 - iCheckCell;
				    assert(id >= 0);
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
				    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
				    iOrder = p->iOrder;
				    ilpAppend(pkd->ilp,
					p->r[0] + blk->xOffset.f[jTile],
					p->r[1] + blk->yOffset.f[jTile],
					p->r[2] + blk->zOffset.f[jTile],
					blk->m.f[jTile], blk->fourh2.f[jTile],
					p->iOrder, v[0], v[1], v[2]);
				    if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
				    }
				else {
				    id = blk->id.i[jTile];
				    assert(id >= 0);
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				    for (pj=c->pLower;pj<=c->pUpper;++pj) {
					if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
					else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
					fMass = pkdMass(pkd,p);
					fSoft = pkdSoft(pkd,p);
					if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
					ilpAppend(pkd->ilp,
					    p->r[0] + blk->xOffset.f[jTile],
					    p->r[1] + blk->yOffset.f[jTile],
					    p->r[2] + blk->zOffset.f[jTile],
					    fMass, 4*fSoft*fSoft,
					    p->iOrder, v[0], v[1], v[2]);
					if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
					}
				    if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				    }
				break;
			    case 2:
				/*
				** Now I am trying to open a bucket, which means I add each particle back on the
				** checklist with a cell size of zero.
				*/
				iCheckCell = blk->iCell.i[jTile];
				assert(iCheckCell>=0);
				id = blk->id.i[jTile];
				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];
				assert(id >= 0);
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				for (pj=c->pLower;pj<=c->pUpper;++pj) {
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
				    fMass = pkdMass(pkd,p);
				    fSoft = pkdSoft(pkd,p);
				    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
				    clAppend(pkd->clNew,-1 - pj,id,0,1,0.0,fMass,4.0f*fSoft*fSoft,
					p->r,    /* center of mass */
					fOffset, /* fOffset */
					p->r,    /* center of box */
					fZero3); /* size of box */
				    if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
				    }
				if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				break;
			    case 3:
				/*
				** Open the cell.
				** Here we ASSUME that the children of
				** c are all in sequential memory order!
				** (the new tree build assures this)
				** (also true for the top tree)
				** We could do a prefetch here for non-local
				** cells.
				*/
				iCheckCell = blk->iCell.i[jTile];
				assert(iCheckCell>=0);
				iCheckLower = blk->iLower.i[jTile];
				assert(iCheckLower > 0);
				id = blk->id.i[jTile];
				if (iCheckLower == ROOT) {
				    /* We must progress to the children of this local tree root cell. */
				    assert(id < 0);
				    id = pkdTopNode(pkd,iCheckCell)->pLower;
				    assert(id >= 0);
				    if (id == pkd->idSelf) {
					c = pkdTreeNode(pkd,ROOT);
					iCheckLower = c->iLower;
					}
				    else {
					c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,ROOT,id));
					iCheckLower = c->iLower;
					mdlRelease(pkd->mdl,CID_CELL,c);
					}
				    if (!iCheckLower) {
					/*
					** The ROOT of a local tree is actually a bucket! An irritating case...
					** This is a rare case though, and as such we simply check the local ROOT bucket once more.
					*/
					nc = c->pUpper - c->pLower + 1;
					pkdNodeBnd(pkd,c,&cbnd);
					fOffset[0] = blk->xOffset.f[jTile];
					fOffset[1] = blk->yOffset.f[jTile];
					fOffset[2] = blk->zOffset.f[jTile];
					clAppend(pkd->clNew,ROOT,id,0,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
					break; /* finished, don't add children */
					}
				    }			    
				cOpen = -1.0f;
				nc = getCell(pkd,iCheckLower,id,&cOpen,&c);
				pkdNodeBnd(pkd,c,&cbnd);
				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];
				iLower = c->iLower;
				if (id == -1 && !iLower) iLower = ROOT;  /* something other than zero for openening crit - iLower usually can't be == ROOT */
				clAppend(pkd->clNew,iCheckLower,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
				if (id >= 0 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				/*
				** Also add the sibling of check->iLower.
				*/
				++iCheckLower;
				cOpen = -1.0f;
				nc = getCell(pkd,iCheckLower,id,&cOpen,&c);
				pkdNodeBnd(pkd,c,&cbnd);
				iLower = c->iLower;
				if (id == -1 && !iLower) iLower = ROOT;  /* something other than zero for openening crit - iLower usually can't be == ROOT */
				clAppend(pkd->clNew,iCheckLower,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
				if (id >= 0 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				break;
			    case 4:
				/*
				** Accept multipole!
				** Interact += Moment(c);
				*/
				iCheckCell = blk->iCell.i[jTile];
				/* Add a particle as a monopole */
				if (iCheckCell < 0) {
				    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) {
					id = blk->id.i[jTile];
					assert(id >= 0);
					pj = -1 - iCheckCell;
					if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
					else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
					v = pkdVel(pkd,p);
					}
				    monoPole.m = blk->m.f[jTile];
				    ilcAppend(pkd->ilc,
					blk->x.f[jTile] + blk->xOffset.f[jTile],
					blk->y.f[jTile] + blk->yOffset.f[jTile],
					blk->z.f[jTile] + blk->zOffset.f[jTile],
					&monoPole,0.0,v[0],v[1],v[2]);
				    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1 && id != pkd->idSelf)
					mdlRelease(pkd->mdl,CID_PARTICLE,p);
				    }
				else {
				    id = blk->id.i[jTile];
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else if (id == -1) c = pkdTopNode(pkd,iCheckCell);
				    else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				    /*
				    ** Center of mass velocity is used by the planets code to get higher derivatives of the 
				    ** acceleration and could be used to drift cell moments as an approximation.
				    */
				    if (pkd->oNodeVelocity) v = pkdNodeVel(pkd,c);
				    ilcAppend(pkd->ilc,
					blk->x.f[jTile] + blk->xOffset.f[jTile],
					blk->y.f[jTile] + blk->yOffset.f[jTile],
					blk->z.f[jTile] + blk->zOffset.f[jTile],
					pkdNodeMom(pkd,c),c->bMax,
					v[0],v[1],v[2]);
				    if (id != -1 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				    }
				break;
			    case 5:
				/*
				** We accept this multipole from the opening criterion, but it is a softened
				** interaction, so we need to treat is as a softened monopole by putting it
				** on the particle interaction list.
				*/
				iCheckCell = blk->iCell.i[jTile];
				if (iCheckCell<0) {
				    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) {
					id = blk->id.i[jTile];
					pj = -1 - iCheckCell;
					assert(id >= 0);
					if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
					else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
					v = pkdVel(pkd,p);
					}
				    ilpAppend(pkd->ilp,
					blk->x.f[jTile] + blk->xOffset.f[jTile],
					blk->y.f[jTile] + blk->yOffset.f[jTile],
					blk->z.f[jTile] + blk->zOffset.f[jTile],
					blk->m.f[jTile], blk->fourh2.f[jTile],
					-1,v[0], v[1], v[2]);
				    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1 && id != pkd->idSelf)
					mdlRelease(pkd->mdl,CID_PARTICLE,p);
				    }
				else {
				    id = blk->id.i[jTile];
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else if (id == -1) c = pkdTopNode(pkd,iCheckCell);
				    else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				    if (pkd->oNodeVelocity) v = pkdNodeVel(pkd,c);
				    ilpAppend(pkd->ilp,
					blk->x.f[jTile] + blk->xOffset.f[jTile],
					blk->y.f[jTile] + blk->yOffset.f[jTile],
					blk->z.f[jTile] + blk->zOffset.f[jTile],
					blk->m.f[jTile], blk->fourh2.f[jTile],
					-1, /* set iOrder to negative value for time step criterion */
					v[0], v[1], v[2]);
				    if (id != -1 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				    }
				break;
			    case 6:
				/*
				** This is accepting an "opened" bucket's particles as monopoles for the
				** local expansion.
				*/
				assert(0);
				iCheckCell = -blk->iCell.i[jTile];
				id = blk->id.i[jTile];
				dOffset[0] = blk->xOffset.f[jTile];
				dOffset[1] = blk->yOffset.f[jTile];
				dOffset[2] = blk->zOffset.f[jTile];
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				for (pj=c->pLower;pj<=c->pUpper;++pj) {
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
				    if (pkdMass(pkd,p) <= 0.0) continue;
				    /*
				    ** Monopole Local expansion accepted!
				    */
				    d2 = 0;
				    for (j=0;j<3;++j) {
					dx[j] = k->r[j] - (p->r[j] + dOffset[j]);
					d2 += dx[j]*dx[j];
					}
				    dir = 1.0/sqrt(d2);
				    *pdFlop += momFlocrAddMono5(&L,k->bMax,pkdMass(pkd,p),dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);

				    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
				    if (adotai > 0) {
					adotai /= maga;
					dirLsum += dir*adotai*adotai;
					normLsum += adotai*adotai;
					}
				    if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
				    }
				if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				break;
			    case 7:
				/*
				** This is accepting an "opened" bucket's particles with this (k) cell as a softened monopole.
				** This is the inverse of accepting a cell as a softened monopole, here we calculate the first 
				** order local expansion of each softened particle of the checkcell.
				*/
				assert(0);
				iCheckCell = -blk->iCell.i[jTile];
				id = blk->id.i[jTile];
				dOffset[0] = blk->xOffset.f[jTile];
				dOffset[1] = blk->yOffset.f[jTile];
				dOffset[2] = blk->zOffset.f[jTile];
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				for (pj=c->pLower;pj<=c->pUpper;++pj) {
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlAquire(pkd->mdl,CID_PARTICLE,pj,id));
				    fMass = pkdMass(pkd,p);
				    if (fMass <= 0.0) continue;
				    fSoft = pkdSoft(pkd,p);
				    momk = pkdNodeMom(pkd,k);
				    /*
				    ** Monopole Local expansion accepted!
				    */
				    d2 = 0;
				    for (j=0;j<3;++j) {
					dx[j] = k->r[j] - (p->r[j] + dOffset[j]);
					d2 += dx[j]*dx[j];
					}
				    dir = 1.0/sqrt(d2);
				    fourh2 = softmassweight(fMass,4*fSoft*fSoft,momk->m,4*k->fSoft2);
				    if (d2 > fourh2) {
					dir = 1.0/sqrt(d2);
					dir2 = dir*dir*dir;
					}
				    else {
					/*
					** This uses the Dehnen K1 kernel function now, it's fast!
					*/
					dir = 1.0/sqrt(fourh2);
					dir2 = dir*dir;
					d2 *= dir2;
					dir2 *= dir;
					d2 = 1 - d2;
					dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
					dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
					}
				    dir2 *= fMass;
				    tax = -dx[0]*dir2;
				    tay = -dx[1]*dir2;
				    taz = -dx[2]*dir2;
				    L.m -= fMass*dir;
				    L.x -= tax;
				    L.y -= tay;
				    L.z -= taz;
				    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
				    if (adotai > 0) {
					adotai /= maga;
					dirLsum += dir*adotai*adotai;
					normLsum += adotai*adotai;
					}
				    if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
				    }
				if (id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				break;
			    case 8:
				/*
				** Local expansion accepted!
				*/
#if 1
				iCheckCell = blk->iCell.i[jTile];
				if (iCheckCell<0) {
				    fOffset[0] = blk->xOffset.f[jTile];
				    fOffset[1] = blk->yOffset.f[jTile];
				    fOffset[2] = blk->zOffset.f[jTile];
				    dx[0] = k->r[0] - (blk->x.f[jTile] + blk->xOffset.f[jTile]);
				    dx[1] = k->r[1] - (blk->y.f[jTile] + blk->yOffset.f[jTile]);
				    dx[2] = k->r[2] - (blk->z.f[jTile] + blk->zOffset.f[jTile]);
				    d2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
				    dir = 1.0/sqrt(d2);
				    monoPole.m = blk->m.f[jTile];

				    *pdFlop += momFlocrAddFmomr5cm(&L,k->bMax,&monoPole,0.0,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);

				    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
				    if (adotai > 0) {
					adotai /= maga;
					dirLsum += dir*adotai*adotai;
					normLsum += adotai*adotai;
					}
				    }
				else {
				    id = blk->id.i[jTile];
				    dOffset[0] = blk->xOffset.f[jTile];
				    dOffset[1] = blk->yOffset.f[jTile];
				    dOffset[2] = blk->zOffset.f[jTile];
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else if (id == -1) c = pkdTopNode(pkd,iCheckCell);
				    else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				    d2 = 0;
				    for (j=0;j<3;++j) {
					dx[j] = k->r[j] - (c->r[j] + dOffset[j]);
					d2 += dx[j]*dx[j];
					}
				    dir = 1.0/sqrt(d2);

				    *pdFlop += momFlocrAddFmomr5cm(&L,k->bMax,pkdNodeMom(pkd,c),c->bMax,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);

				    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
				    if (adotai > 0) {
					adotai /= maga;
					dirLsum += dir*adotai*adotai;
					normLsum += adotai*adotai;
					}
				    if (id != -1 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				    }
#endif
				break;
			    case 9:
				/*
				** Here we compute the local expansion due to a single monopole term which could be softened.
				*/
				assert(0);
				iCheckCell = blk->iCell.i[jTile];
				id = blk->id.i[jTile];
				dOffset[0] = blk->xOffset.f[jTile];
				dOffset[1] = blk->yOffset.f[jTile];
				dOffset[2] = blk->zOffset.f[jTile];
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else if (id == -1) c = pkdTopNode(pkd,iCheckCell);
				else c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iCheckCell,id));
				momk = pkdNodeMom(pkd,k);
				momc = pkdNodeMom(pkd,c);
				d2 = 0;
				for (j=0;j<3;++j) {
				    dx[j] = k->r[j] - (c->r[j] + dOffset[j]);
				    d2 += dx[j]*dx[j];
				    }
				dir = 1.0/sqrt(d2);
				fourh2 = softmassweight(momc->m,4*c->fSoft2,momk->m,4*k->fSoft2);
				if (d2 > fourh2) {
				    dir = 1.0/sqrt(d2);
				    dir2 = dir*dir*dir;
				    }
				else {
				    /*
				    ** This uses the Dehnen K1 kernel function now, it's fast!
				    */
				    dir = 1.0/sqrt(fourh2);
				    dir2 = dir*dir;
				    d2 *= dir2;
				    dir2 *= dir;
				    d2 = 1 - d2;
				    dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
				    dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
				    }
				dir2 *= momc->m;
				tax = -dx[0]*dir2;
				tay = -dx[1]*dir2;
				taz = -dx[2]*dir2;
				L.m -= momc->m*dir;
				L.x -= tax;
				L.y -= tay;
				L.z -= taz;
				adotai = a[0]*tax + a[1]*tay + a[2]*taz;
				if (adotai > 0) {
				    adotai /= maga;
				    dirLsum += dir*adotai*adotai;
				    normLsum += adotai*adotai;
				    }
				if (id != -1 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				break;
			    case 10:
				/*
				** This checkcell is removed from the checklist since it has zero/negative mass.
				*/
				break;		
				}
			    }
			}
		    } /* end of CL_LOOP */
		clTemp = pkd->cl;
		pkd->cl = pkd->clNew;
		pkd->clNew = clTemp;
		} while (clCount(pkd->cl));
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
//	    if (!k->iLower) break;
	    if ((k->pUpper-k->pLower+1)<=PKD_GROUP_SIZE) break;
	    xParent = k->r[0];
	    yParent = k->r[1];
	    zParent = k->r[2];
	    bMaxParent = k->bMax;
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    kOpen = -1.0f;
	    iCell = k->iLower;
	    nk = getCell(pkd,iCell,pkd->idSelf,&kOpen,&k);
	    pkdNodeBnd(pkd,k,&kbnd);
	    /*
	    ** Check iCell is active. We eventually want to just to a
	    ** rung check here when we start using tree repair, but
	    ** for now this is just as good.
	    */
	    if (k->nActive) {
		/*
		** iCell is active, continue processing it.
		** Put the sibling onto the checklist.
		*/
		cOpen = -1.0f;
		iSib = iCell+1;
		nc = getCell(pkd,iSib,pkd->idSelf,&cOpen,&c);
		if (c->nActive) {
		    /*
		    ** Sibling is active so we need to clone the checklist!
		    */
		    clClone(pkd->cl,pkd->S[iStack+1].cl);
		    }
		else {
		    /*
		    ** Otherwise we can simple grab it.
		    */
		    clTemp = pkd->cl;
		    pkd->cl = pkd->S[iStack+1].cl;
		    pkd->S[iStack+1].cl = clTemp;
		    }
		pkdNodeBnd(pkd,c,&cbnd);
		clAppend(pkd->cl,iSib,pkd->idSelf,c->iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd.fCenter,cbnd.fMax);
		/*
		** Test whether the sibling is active as well.
		** If not we don't push it onto the stack, but we
		** have to be careful to not pop the stack when we
		** hit the sibling. See the goto InactiveAscend below
		** for how this is done.
		*/
		if (c->nActive) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    /*
		    ** At this point, the ILP is normally empty if you never do P-P except when reaching a bucket.
		    ** Softened multi-poles are also an exception.
		    */
		    ilpCheckPt(pkd->ilp,&pkd->S[iStack].PartChkPt);
		    ilcCheckPt(pkd->ilc,&pkd->S[iStack].CellChkPt);
		    /*
		    ** Note here we already have the correct elements in S[iStack] (iStack+1 was used previously), just need to add one.
		    */
		    clAppend(pkd->S[iStack].cl,iCell,pkd->idSelf,k->iLower,nk,kOpen,pkdNodeMom(pkd,k)->m,4.0f*k->fSoft2,k->r,fOffset,kbnd.fCenter,kbnd.fMax);
		    pkd->S[iStack].L = L;
		    pkd->S[iStack].dirLsum = dirLsum;
		    pkd->S[iStack].normLsum = normLsum;
		    dShiftFlop = momShiftFlocr(&pkd->S[iStack].L,bMaxParent,
					      c->r[0] - xParent,
					      c->r[1] - yParent,
					      c->r[2] - zParent);
		    momRescaleFlocr(&pkd->S[iStack].L,c->bMax,bMaxParent);
		    pkd->S[iStack].fWeight = (*pdFlop-tempI) + dShiftFlop;
		    pkd->S[iStack].fWeight += dEwFlop;
		    }
		}
	    else {
		/*
		** Skip iCell, but add it to the Checklist.
		** No need to push anything onto the stack here.
		** We can simply grab it the checklist from the Stack.
		*/
		clTemp = pkd->cl;
		pkd->cl = pkd->S[iStack+1].cl;
		pkd->S[iStack+1].cl = clTemp;
		clAppend(pkd->cl,iCell,pkd->idSelf,k->iLower,nk,kOpen,pkdNodeMom(pkd,k)->m,4.0f*k->fSoft2,k->r,fOffset,kbnd.fCenter,kbnd.fMax);
		/*
		** Move onto processing the sibling.
		*/
		k = pkdTreeNode(pkd,++iCell);
		}
	    *pdFlop += momShiftFlocr(&L,bMaxParent,k->r[0] - xParent,
				    k->r[1] - yParent,
				    k->r[2] - zParent);
	    momRescaleFlocr(&L,k->bMax,bMaxParent);
	    }
	/*
	** Now the interaction list should be complete and the
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
//	nActive = pkdGravInteract(pkd,uRungLo,uRungHi,k,&L,pkd->ilp,pkd->ilc,
//				  dirLsum,normLsum,bEwald,pdFlop,&dEwFlop,dRhoFac,
//				  smx, &smf);
	/*
	** Update the limit for a shift of the center here based on the opening radius of this
	** cell (the one we just evaluated).
	*/
	fWeight += (*pdFlop-tempI);
	fWeight += dEwFlop;
	if (nActive) {
	    fWeight /= nActive;
	    /*
	    ** Here we used to set the weights of particles based on the work done, but now we just assume that
	    ** all particles cost the same in domain decomposition, so we really don't need to set anything here.
	    */
	    *pdPartSum += nActive*ilpCount(pkd->ilp);
	    *pdCellSum += nActive*ilcCount(pkd->ilc);
	    nTotActive += nActive;
	    }

	while (iCell & 1) {
	InactiveAscend:
	    k = pkdTreeNode(pkd,iCell = k->iParent);
	    if (!iCell) {
		/*
		** Make sure stack is empty.
		*/
		assert(iStack == -1);
#ifdef PROFILE_GRAVWALK
		VTPause();
#endif
		*pdFlop += dEwFlop;   /* Finally add the ewald score to get a proper float count */
		if (smx) {
		    smSmoothFinish(smx);
		    smFinish(smx,&smf);
		    }
		return(nTotActive);
		}
	    }
	k = pkdTreeNode(pkd,++iCell);
	if (!k->nActive) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	ilpRestore(pkd->ilp,&pkd->S[iStack].PartChkPt);
	ilcRestore(pkd->ilc,&pkd->S[iStack].CellChkPt);
	/*
	** Grab the checklist from the stack.
	*/
	clTemp = pkd->cl;
	assert(clCount(pkd->cl) == 0);
	pkd->cl = pkd->S[iStack].cl;
	pkd->S[iStack].cl = clTemp;
	L = pkd->S[iStack].L;
	dirLsum = pkd->S[iStack].dirLsum;
	normLsum = pkd->S[iStack].normLsum;
	fWeight = pkd->S[iStack].fWeight;
	tempI = *pdFlop;
	tempI += dEwFlop;
	--iStack;
	}
    }


int pkdRungWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,void *pParams,void *doFunc(PKD pkd,PARTICLE *p,void *pParams)) {
    KDN *k;
    PARTICLE *p;
    int iCell,pj;
    int nTotActive = 0;

    k = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
	if (k->uMaxRung >= uRungLo || k->uMinRung <= uRungHi ) {
	    /*
	    ** Descend.
	    */
	    if (!k->iLower) {
		/*
		** Process bucket!
		*/
		for (pj=k->pLower;pj<=k->pUpper;++pj) {
		    p = pkdParticle(pkd,pj);
		    if (pkdIsRungRange(p,uRungLo,uRungHi)) {
			if (doFunc) doFunc(pkd,p,pParams);
			++nTotActive;
			}
		    }
		}
	    else k = pkdTreeNode(pkd,iCell = k->iLower);
	    }
	while (iCell & 1) {
	    k = pkdTreeNode(pkd,iCell = k->iParent);
	    if (!iCell) return(nTotActive);
	    }
	k = pkdTreeNode(pkd,++iCell);
	}
    }
