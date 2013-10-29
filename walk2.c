#ifdef HAVE_CONFIG_H
#include "config.h"
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
    vfloat zero;
    vfloat one;
    vfloat threehalves;
    vfloat two;
    } consts = {
        {SIMD_CONST(0.0)},
	{SIMD_CONST(1.0)},
	{SIMD_CONST(1.5)},
	{SIMD_CONST(2.0)},
};
static const struct ICONSTS {
    vint zero;
    vint one;
    vint two;
    vint three;
    vint four;
    vint five;
    vint six;
    vint seven;
    vint eight;
    /* no nine */
    vint ten;
    vint walk_min_multipole;
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
	{SIMD_CONST(3)},
};

static union {
    uint32_t u[SIMD_WIDTH];
    v_sf p;
    } const_fabs = {SIMD_CONST(0x7fffffff)};

/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version will also open buckets ("new" criteria)
*/
static void iOpenOutcomeSIMD(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin, int nGroup) {
    v_sf T0,T1,T2,T3,T4,T6,T7,P1,P2,P3,P4;
    v_sf T,xc,yc,zc,dx,dy,dz,d2,diCrit,cOpen,cOpen2,d2Open,mink2,minbnd2,fourh2;
    int i,n,iEnd,nLeft;
    CL_BLK *blk;
    v_sf iOpen,iOpenA,iOpenB;
    const BND *kbnd;
    v_sf k_xCenter, k_yCenter, k_zCenter, k_xMax, k_yMax, k_zMax;
    v_sf k_xMinBnd, k_yMinBnd, k_zMinBnd, k_xMaxBnd, k_yMaxBnd, k_zMaxBnd;
    v_sf k_x, k_y, k_z, k_m, k_bMax, k_Open;
    v_i  k_nk;
    vint k_nGroup = {SIMD_CONST(nGroup)};

    assert ( pkdNodeMom(pkd,k)->m > 0.0f );

    diCrit = SIMD_SPLAT(1.0f/dThetaMin);

    kbnd = pkdNodeBnd(pkd,k);
    k_xMinBnd = SIMD_SPLAT(kbnd->fCenter[0]-kbnd->fMax[0]);
    k_yMinBnd = SIMD_SPLAT(kbnd->fCenter[1]-kbnd->fMax[1]);
    k_zMinBnd = SIMD_SPLAT(kbnd->fCenter[2]-kbnd->fMax[2]);
    k_xMaxBnd = SIMD_SPLAT(kbnd->fCenter[0]+kbnd->fMax[0]);
    k_yMaxBnd = SIMD_SPLAT(kbnd->fCenter[1]+kbnd->fMax[1]);
    k_zMaxBnd = SIMD_SPLAT(kbnd->fCenter[2]+kbnd->fMax[2]);
    k_xCenter = SIMD_SPLAT(kbnd->fCenter[0]);
    k_yCenter = SIMD_SPLAT(kbnd->fCenter[1]);
    k_zCenter = SIMD_SPLAT(kbnd->fCenter[2]);
    k_xMax = SIMD_SPLAT(kbnd->fMax[0]);
    k_yMax = SIMD_SPLAT(kbnd->fMax[1]);
    k_zMax = SIMD_SPLAT(kbnd->fMax[2]);
    k_x = SIMD_SPLAT(k->r[0]);
    k_y = SIMD_SPLAT(k->r[1]);
    k_z = SIMD_SPLAT(k->r[2]);
    k_m = SIMD_SPLAT(pkdNodeMom(pkd,k)->m);
    k_bMax = SIMD_SPLAT(k->bMax);
    k_nk = SIMD_SPLATI32(k->pUpper-k->pLower+1);
    k_Open = SIMD_MUL(consts.threehalves.p,SIMD_MUL(k_bMax,diCrit));

    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	iEnd = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	iEnd = (iEnd+SIMD_MASK) >> SIMD_BITS;
	for(i=0; i<iEnd; ++i) {
	    fourh2 = blk->fourh2.p[i];
	    xc = SIMD_ADD(blk->x.p[i],blk->xOffset.p[i]);
	    yc = SIMD_ADD(blk->y.p[i],blk->yOffset.p[i]);
	    zc = SIMD_ADD(blk->z.p[i],blk->zOffset.p[i]);
	    dx = SIMD_SUB(k_x,xc);
	    dy = SIMD_SUB(k_y,yc);
	    dz = SIMD_SUB(k_z,zc);
	    d2 = SIMD_MADD(dx,dx,SIMD_MADD(dy,dy,SIMD_MUL(dz,dz)));
	    cOpen = blk->cOpen.p[i];
	    cOpen2 = SIMD_MUL(cOpen,cOpen);
	    d2Open = SIMD_ADD(cOpen,k_Open);
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


	    T6 = SIMD_CMP_GT(cOpen,k_Open);
	    T7 = SIMD_I2F(SIMD_CMP_GT_EPI32(k_nk,k_nGroup.p));
	    iOpenA = SIMD_OR(SIMD_AND(T2,iconsts.one.pf),SIMD_ANDNOT(T2,iconsts.three.pf));
	    iOpenB = SIMD_OR(SIMD_AND(T3,iOpenA),SIMD_ANDNOT(T3,
		    SIMD_OR(SIMD_AND(T4,iconsts.four.pf),SIMD_ANDNOT(T4,iOpenA))));
	    P1 = SIMD_OR(SIMD_AND(T2,iconsts.two.pf),SIMD_ANDNOT(T2,iconsts.three.pf));
	    P2 = SIMD_OR(SIMD_AND(T7,iconsts.zero.pf),SIMD_ANDNOT(T7,iOpenB));
	    P3 = SIMD_OR(SIMD_AND(T6,P1),SIMD_ANDNOT(T6,P2));
	    P4 = SIMD_OR(SIMD_AND(T1,iconsts.eight.pf),SIMD_ANDNOT(T1,P3));
	    iOpen = SIMD_OR(SIMD_AND(T0,P4),SIMD_ANDNOT(T0,iconsts.ten.pf));
	    blk->iOpen.pf[i] = iOpen;
	    }
	}
    }

#endif
/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version has been changed by adding the ability to open buckets.
*/
static void iOpenOutcomeCL(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin,int nGroup) {
    const int walk_min_multipole = 3;
    float dx,dy,dz,mink2,d2,d2Open,xc,yc,zc,fourh2,minbnd2,kOpen,cOpen,diCrit;
    int j,i,nk;
    int iOpen,iOpenA,iOpenB;
    CL_BLK *blk;
    int n, nLeft;
    const BND *kbnd;

    kbnd = pkdNodeBnd(pkd,k);
    nk = k->pUpper - k->pLower + 1;

    diCrit = 1.0f / dThetaMin;
    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	n = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	for(i=0; i<n; ++i) {
	    if (blk->m.f[i] <= 0) iOpen = 10;  /* ignore this cell */
	    else {
		fourh2 = blk->fourh2.f[i];
		xc = blk->x.f[i] + blk->xOffset.f[i];
		yc = blk->y.f[i] + blk->yOffset.f[i];
		zc = blk->z.f[i] + blk->zOffset.f[i];
		d2 = pow(k->r[0]-xc,2) + pow(k->r[1]-yc,2) + pow(k->r[2]-zc,2);
		kOpen = 1.5f * k->bMax * diCrit;
		cOpen = blk->cOpen.f[i];
		d2Open = pow(cOpen+kOpen,2);
		dx = fabs(xc - kbnd->fCenter[0]) - kbnd->fMax[0];
		dy = fabs(yc - kbnd->fCenter[1]) - kbnd->fMax[1];
		dz = fabs(zc - kbnd->fCenter[2]) - kbnd->fMax[2];
		mink2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);
		minbnd2 = 0;

		dx = kbnd->fCenter[0] - kbnd->fMax[0] -  blk->xCenter.f[i] - blk->xOffset.f[i] - blk->xMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->xCenter.f[i] + blk->xOffset.f[i] - blk->xMax.f[i] - kbnd->fCenter[0] - kbnd->fMax[0];
		if (dx > 0) minbnd2 += dx*dx;

		dx = kbnd->fCenter[1] - kbnd->fMax[1] - blk->yCenter.f[i] - blk->yOffset.f[i] - blk->yMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->yCenter.f[i] + blk->yOffset.f[i] - blk->yMax.f[i] - kbnd->fCenter[1] - kbnd->fMax[1];
		if (dx > 0) minbnd2 += dx*dx;

		dx = kbnd->fCenter[2] - kbnd->fMax[2] - blk->zCenter.f[i] - blk->zOffset.f[i] - blk->zMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->zCenter.f[i] + blk->zOffset.f[i] - blk->zMax.f[i] - kbnd->fCenter[2] - kbnd->fMax[2];
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
		    else iOpenB = iOpenA;
		    if (nk>nGroup) iOpen = 0;
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
** Returns total number of active particles for which gravity was calculated.
*/
static int processCheckList(PKD pkd, SMX smx, SMF smf, int iRoot, int iVARoot,
    uint8_t uRungLo,uint8_t uRungHi, double dRhoFac, int bEwald, int nGroup, double dThetaMin,
    int bGravStep, double *pdFlop, double *pdPartSum,double *pdCellSum) {
    KDN *k,*c,*kFind;
    int id,iCell,iSib,iLower,iCheckCell,iCheckLower,iCellDescend;
    PARTICLE *p;
    FMOMR *momc,*momk;
    FMOMR monoPole;
    LOCR L;
    double cx,cy,cz,d2c;
    double fWeight = 0.0;
    double dShiftFlop;
    const double *v, *a;
    double dOffset[3];
    double xParent,yParent,zParent;
    double d2,fourh2;
    double dx[3],dir,dir2;
    double tax,tay,taz;
    float fOffset[3];
    float bMaxParent;
    float dirLsum,normLsum,adotai,maga;
    float fMass,fSoft;
    uint64_t iOrder;
    int iStack;
    int j,jTile,pi,pj,nActive,nTotActive;
    float cOpen,kOpen;
    const BND *cbnd,*kbnd;
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
    double tempI;
    double dEwFlop = 0.0;

    iStack = -1;

    /*
    ** Clear monopole sentinel and local expansion and the timestepping sums.
    */
    momClearFmomr(&monoPole);
    momClearLocr(&L);
    dirLsum = 0;
    normLsum = 0;
    nTotActive = 0;

    a = dZero3;
    v = dZero3;


    /*
    ** We are now going to work on the local tree.
    ** Make iCell point to the root of the tree again.
    */
    if (iVARoot) k = pkdTreeNode(pkd,iCell = iVARoot);
    else k = pkdTreeNode(pkd,iCell = iRoot);

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
	    if (bGravStep) {
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
		    iOpenOutcomeSIMD(pkd,k,pkd->cl,cltile,dThetaMin,nGroup);
		    /*Verify:iOpenOutcomeNewCL(pkd,k,pkd->cl,cltile,dThetaMin);*/
#else
		    iOpenOutcomeCL(pkd,k,pkd->cl,cltile,dThetaMin,nGroup);
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
				    if (bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
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
					if (bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
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
				    if (bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
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
				if (iCheckLower == iRoot) {
				    /* We must progress to the children of this local tree root cell. */
				    assert(id < 0);
				    id = pkdTopNode(pkd,iCheckCell)->pLower;
				    assert(id >= 0);
				    if (id == pkd->idSelf) {
					c = pkdTreeNode(pkd,iRoot);
					iCheckLower = c->iLower;
					}
				    else {
					c = CAST(KDN *,mdlAquire(pkd->mdl,CID_CELL,iRoot,id));
					iCheckLower = c->iLower;
					mdlRelease(pkd->mdl,CID_CELL,c);
					}
				    if (!iCheckLower) {
					/*
					** The iRoot of a local tree is actually a bucket! An irritating case...
					** This is a rare case though, and as such we simply check the local iRoot bucket once more.
					*/
					nc = getCell(pkd,iRoot,id,&cOpen,&c);
					cbnd = pkdNodeBnd(pkd,c);
					fOffset[0] = blk->xOffset.f[jTile];
					fOffset[1] = blk->yOffset.f[jTile];
					fOffset[2] = blk->zOffset.f[jTile];
					clAppend(pkd->clNew,iRoot,id,0,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
					break; /* finished, don't add children */
					}
				    }			    
				cOpen = -1.0f;
				nc = getCell(pkd,iCheckLower,id,&cOpen,&c);
				cbnd = pkdNodeBnd(pkd,c);
				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];
				iLower = c->iLower;
				if (id == -1 && !iLower) iLower = iRoot;  /* something other than zero for openening crit - iLower usually can't be == iRoot */
				clAppend(pkd->clNew,iCheckLower,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
				if (id >= 0 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				/*
				** Also add the sibling of check->iLower.
				*/
				++iCheckLower;
				cOpen = -1.0f;
				nc = getCell(pkd,iCheckLower,id,&cOpen,&c);
				cbnd = pkdNodeBnd(pkd,c);
				iLower = c->iLower;
				if (id == -1 && !iLower) iLower = iRoot;  /* something other than zero for openening crit - iLower usually can't be == iRoot */
				clAppend(pkd->clNew,iCheckLower,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
				if (id >= 0 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				break;
			    case 4:
				/*
				** Accept multipole!
				** Interact += Moment(c);
				*/
				iCheckCell = blk->iCell.i[jTile];
				assert(iCheckCell>=0);
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
				break;
			    case 8:
				/*
				** Local expansion accepted!
				*/
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
				    /* monoPole.m = blk->m.f[jTile];*/
				    /* *pdFlop += momLocrAddFmomr5cm(&L,&monoPole,0.0,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);*/
				    *pdFlop += momLocrAddMono5(&L,blk->m.f[jTile],dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
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

				    if (pkd->param.bCenterOfMassExpand) { 
					*pdFlop += momLocrAddFmomr5cm(&L,pkdNodeMom(pkd,c),c->bMax,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
					}
				    else {
					*pdFlop += momLocrAddFmomr5(&L,pkdNodeMom(pkd,c),c->bMax,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
					}
				    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
				    if (adotai > 0) {
					adotai /= maga;
					dirLsum += dir*adotai*adotai;
					normLsum += adotai*adotai;
					}
				    if (id != -1 && id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
				    }
				break;
			    case 10:
				/*
				** This checkcell is removed from the checklist since it has zero/negative mass.
				*/
				break;		
			    default:
				assert(0);
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
	    if ((k->pUpper-k->pLower+1)<=nGroup) break;
	    xParent = k->r[0];
	    yParent = k->r[1];
	    zParent = k->r[2];
	    bMaxParent = k->bMax;
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    kOpen = -1.0f;
	    iCell = k->iLower;
	    nk = getCell(pkd,iCell,pkd->idSelf,&kOpen,&k);
	    kbnd = pkdNodeBnd(pkd,k);
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
		cbnd = pkdNodeBnd(pkd,c);
		clAppend(pkd->cl,iSib,pkd->idSelf,c->iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
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
		    clAppend(pkd->S[iStack].cl,iCell,pkd->idSelf,k->iLower,nk,kOpen,pkdNodeMom(pkd,k)->m,4.0f*k->fSoft2,k->r,fOffset,kbnd->fCenter,kbnd->fMax);
		    pkd->S[iStack].L = L;
		    pkd->S[iStack].dirLsum = dirLsum;
		    pkd->S[iStack].normLsum = normLsum;
		    dShiftFlop = momShiftLocr(&pkd->S[iStack].L,
					      c->r[0] - xParent,
					      c->r[1] - yParent,
					      c->r[2] - zParent);
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
		clAppend(pkd->cl,iCell,pkd->idSelf,k->iLower,nk,kOpen,pkdNodeMom(pkd,k)->m,4.0f*k->fSoft2,k->r,fOffset,kbnd->fCenter,kbnd->fMax);
		/*
		** Move onto processing the sibling.
		*/
		k = pkdTreeNode(pkd,++iCell);
		}
	    *pdFlop += momShiftLocr(&L,k->r[0] - xParent,
				    k->r[1] - yParent,
				    k->r[2] - zParent);
	    }
	/*
	** Now the interaction list should be complete and the
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
	if (!pkd->param.bNoGrav) {
	    nActive = pkdGravInteract(pkd,uRungLo,uRungHi,k,&L,pkd->ilp,pkd->ilc,
		dirLsum,normLsum,bEwald,bGravStep,nGroup,pdFlop,&dEwFlop,dRhoFac,
		smx, &smf);
	    }
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

static void initGravWalk(PKD pkd,double dTime,double dThetaMin,double dThetaMax,int bPeriodic,int bGravStep,
    SMX *smx, SMF *smf, double *dRhoFac) {
    int pi;
    PARTICLE *p;

    /*
    ** If necessary, calculate the theta interpolation tables.
    */
#ifdef USE_DEHNEN_THETA
    pkdSetThetaTable(pkd,dThetaMin,dThetaMax);
#else
    pkd->fiCritTheta = 1.0f / dThetaMin;
#endif

    assert(pkd->oNodeMom);
    if (bGravStep) {
	assert(pkd->oNodeAcceleration);
	if (pkd->param.iTimeStepCrit == 1) {
	    assert(pkd->oNodeVelocity);
	    assert(pkd->oVelocity);
	    }
	}

    /*
    ** Setup smooth for calculating local densities when a particle has too few P-P interactions.
    */
    if (bGravStep) {
	smInitializeRO(smx,pkd,smf,pkd->param.nPartRhoLoc,bPeriodic,SMX_DENSITY_F1);
	smSmoothInitialize(*smx);
	/* No particles are inactive for density calculation */
	for (pi=0;pi<pkd->nLocal;++pi) {
	    p = pkdParticle(pkd,pi);
	    (*smx)->ea[pi].bInactive = 0;
	    }
	}
    else (*smx) = NULL;

    /*
    ** Precalculate RhoFac if required.
    */
    if (bGravStep) {
	double a = csmTime2Exp(pkd->param.csm,dTime);
	*dRhoFac = 1.0/(a*a*a);
	}
    else *dRhoFac = 0.0;


    }

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalkHop(PKD pkd,double dTime,int nGroup, double dThetaMin,double dThetaMax,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p;
    KDN *c;
    int id,iRoot,iRootSelf;
    float fOffset[3];
    int nActive;
    int i,j,gid;
    float cOpen;
    const BND *cbnd;
    int nc;
    double dRhoFac;
    SMX smx;
    SMF smf;

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd), pkdLocal(pkd));
    initGravWalk(pkd,dTime,dThetaMin,dThetaMax,0,0,&smx,&smf,&dRhoFac);
    nActive = 0;
    for(gid=1; gid<pkd->nGroups; ++gid) {
	if (!pkd->hopGroups[gid].bNeedGrav) continue;
	pkd->hopGroups[gid].bNeedGrav = 0;
	ilpClear(pkd->ilp);
	ilcClear(pkd->ilc);
	clClear(pkd->cl);
	iRootSelf = pkd->hopRootIndex[gid];
	for (i=pkd->hopRootIndex[gid]; i<pkd->hopRootIndex[gid+1]; ++i) {
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    cOpen = -1.0f;
	    id = pkd->hopRoots[i].iPid;
	    iRoot = pkd->hopRoots[i].iIndex;
	    assert(iRoot>0);
	    nc = getCell(pkd,iRoot,id,&cOpen,&c);
	    cbnd = pkdNodeBnd(pkd,c);
	    clAppend(pkd->cl,iRoot,id,c->iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
	    }
	assert(pkd->hopRoots[iRootSelf].iPid==pkd->idSelf);
	nActive += processCheckList(pkd, smx, smf, pkd->hopRoots[iRootSelf].iIndex, 0, 0, MAX_RUNG, dRhoFac, 0, nGroup, dThetaMin, 0, pdFlop, pdPartSum, pdCellSum);
	}

    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    return nActive;
}


/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,int nGroup, int iRoot, int iVARoot,
		double dThetaMin,double dThetaMax,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p;
    KDN *c;
    int id,iCell,iSib,iLower;
    float fOffset[3];
    int ix,iy,iz,bRep;
    int pi;
    int j;
    float cOpen;
    const BND *cbnd;
    int nc;
    double dRhoFac;
    SMX smx;
    SMF smf;

    initGravWalk(pkd,dTime,dThetaMin,dThetaMax,nReps?1:0,pkd->param.bGravStep,&smx,&smf,&dRhoFac);

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    ** Otherwise we check that the iRoot has active particles!
    */
    if (iVARoot) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkdTreeNode(pkd,iVARoot)->pUpper - pkdTreeNode(pkd,iVARoot)->pLower + 1);
	}
    else if (!pkdIsCellActive(pkdTreeNode(pkd,iRoot),uRungLo,uRungHi)) return 0;
    /*
    ** Initially we set our cell pointer to
    ** point to the top tree.
    */
    ilpClear(pkd->ilp);
    ilcClear(pkd->ilc);
    clClear(pkd->cl);

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
		if (bRep || iVARoot) {
		    /* 
		    ** Use leaf of the top tree and NOT the root of the local tree here.
		    */
		    cOpen = -1.0f;
		    id = -1;
		    iLower = pkdTopNode(pkd,iRoot)->iLower;
		    if (!iLower) iLower = iRoot;  /* something other than zero for openening crit - iLower usually can't be == iRoot */
		    nc = getCell(pkd,iRoot,id,&cOpen,&c);
		    cbnd = pkdNodeBnd(pkd,c);
		    clAppend(pkd->cl,iRoot,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
		    }
		if (bRep && iVARoot) {
		    /*
		    ** Add the images of the very active tree to the checklist.
		    */
		    cOpen = -1.0f;
		    id = mdlSelf(pkd->mdl);
		    iLower = pkdTopNode(pkd,iVARoot)->iLower;
		    if (!iLower) iLower = iRoot;  /* something other than zero for openening crit - iLower usually can't be == iRoot */
		    nc = getCell(pkd,iVARoot,id,&cOpen,&c);
		    cbnd = pkdNodeBnd(pkd,c);
		    clAppend(pkd->cl,iVARoot,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
		    }
		}
	    }
	}
    if (!iVARoot) {
	/*
	** This adds all siblings of a chain leading from the local tree leaf in the top
	** tree up to the iRoot of the top tree.
	*/
	for (j=0;j<3;++j) fOffset[j] = 0.0f;
	iCell = pkd->iTopRoot;
	iSib = SIBLING(iCell);
	while (iSib) {
	    cOpen = -1.0f;
	    id = -1;
	    iLower = pkdTopNode(pkd,iSib)->iLower;
	    if (!iLower) iLower = iRoot;  /* something other than zero for openening crit - iLower usually can't be == iRoot */
	    nc = getCell(pkd,iSib,id,&cOpen,&c);
	    cbnd = pkdNodeBnd(pkd,c);
	    clAppend(pkd->cl,iSib,id,iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
	    iCell = pkdTopNode(pkd,iCell)->iParent;
	    iSib = SIBLING(iCell);
	    }
	}

    return processCheckList(pkd, smx, smf, iRoot, iVARoot, uRungLo, uRungHi, dRhoFac, bEwald, nGroup, dThetaMin, pkd->param.bGravStep, pdFlop, pdPartSum, pdCellSum);
    }

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalkGroups(PKD pkd,double dTime,int nGroup, double dThetaMin,double dThetaMax,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p;
    KDN *c;
    int id,iRoot;
    float fOffset[3];
    int ix,iy,iz,bRep;
    int pi;
    int i,j,k;
    float cOpen;
    const BND *cbnd;
    int nc;
    double dRhoFac;
    SMX smx;
    SMF smf;

    initGravWalk(pkd,dTime,dThetaMin,dThetaMax,0,0,&smx,&smf,&dRhoFac);


    /*
    ** Initially we set our cell pointer to
    ** point to the top tree.
    */

    struct psGroup *gd = pkd->psGroupTable.pGroup;
    int nActive=0;
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	ilpClear(pkd->ilp);
	ilcClear(pkd->ilc);
	clClear(pkd->cl);
#if 1
	for (k=1; k < gd[i].nTreeRoots; k++)
	{
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    cOpen = -1.0f;
	    id = gd[i].treeRoots[k].iPid;
	    iRoot = gd[i].treeRoots[k].iLocalRootId;
	    nc = getCell(pkd,iRoot,id,&cOpen,&c);
	    cbnd = pkdNodeBnd(pkd,c);
	    clAppend(pkd->cl,iRoot,id,c->iLower,nc,cOpen,pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
        }
#endif
	nActive += processCheckList(pkd, smx, smf, gd[i].treeRoots[0].iLocalRootId, 0, 0, MAX_RUNG, dRhoFac, 0, nGroup, dThetaMin, 0, pdFlop, pdPartSum, pdCellSum);
    }
    return nActive;
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
