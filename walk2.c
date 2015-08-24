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

static inline int getCell(PKD pkd,int iCell,int id,float *pcOpen,KDN **pc) {
    KDN *c;
    int nc;
    assert(iCell > 0);
    assert(id >= 0);
    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCell);
    else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCell,id));
    *pc = c;
    if (c->bRemote) nc = 1000000000; /* we never allow pp with this cell */
    else nc = c->pUpper - c->pLower + 1;
    *pcOpen = c->bMax * pkd->fiCritTheta;
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
    v_sf xc,yc,zc,dx,dy,dz,d2,diCrit,cOpen,cOpen2,d2Open,mink2,minbnd2,fourh2;
    int i,iEnd,nLeft;
    CL_BLK *blk;
    v_sf iOpen,iOpenA,iOpenB;
    const BND *kbnd;
    v_sf k_xCenter, k_yCenter, k_zCenter, k_xMax, k_yMax, k_zMax;
    v_sf k_xMinBnd, k_yMinBnd, k_zMinBnd, k_xMaxBnd, k_yMaxBnd, k_zMaxBnd;
    v_sf k_x, k_y, k_z, k_bMax, k_Open;
    v_i  k_nk;
    vint k_nGroup = {SIMD_CONST(nGroup)};

    assert ( pkdNodeMom(pkd,k)->m > 0.0f );

    diCrit = SIMD_SPLAT(1.0f/dThetaMin);

    kbnd = pkdNodeBnd(pkd,k);
    k_xMinBnd = SIMD_SPLAT((float)(kbnd->fCenter[0]-kbnd->fMax[0]));
    k_yMinBnd = SIMD_SPLAT((float)(kbnd->fCenter[1]-kbnd->fMax[1]));
    k_zMinBnd = SIMD_SPLAT((float)(kbnd->fCenter[2]-kbnd->fMax[2]));
    k_xMaxBnd = SIMD_SPLAT((float)(kbnd->fCenter[0]+kbnd->fMax[0]));
    k_yMaxBnd = SIMD_SPLAT((float)(kbnd->fCenter[1]+kbnd->fMax[1]));
    k_zMaxBnd = SIMD_SPLAT((float)(kbnd->fCenter[2]+kbnd->fMax[2]));
    k_xCenter = SIMD_SPLAT((float)(kbnd->fCenter[0]));
    k_yCenter = SIMD_SPLAT((float)(kbnd->fCenter[1]));
    k_zCenter = SIMD_SPLAT((float)(kbnd->fCenter[2]));
    k_xMax = SIMD_SPLAT((float)(kbnd->fMax[0]));
    k_yMax = SIMD_SPLAT((float)(kbnd->fMax[1]));
    k_zMax = SIMD_SPLAT((float)(kbnd->fMax[2]));
    k_x = SIMD_SPLAT((float)(k->r[0]));
    k_y = SIMD_SPLAT((float)(k->r[1]));
    k_z = SIMD_SPLAT((float)(k->r[2]));
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
 	    iOpenA = SIMD_SELECT(iconsts.three.pf,iconsts.one.pf,T2);
	    iOpenB = SIMD_SELECT(SIMD_SELECT(iOpenA,iconsts.four.pf,T4),iOpenA,T3);
	    P1 = SIMD_SELECT(iconsts.three.pf,iconsts.two.pf,T2);
	    P2 = SIMD_SELECT(iOpenB,iconsts.zero.pf,T7);
	    P3 = SIMD_SELECT(P2,P1,T6);
	    P4 = SIMD_SELECT(P3,iconsts.eight.pf,T1);
	    iOpen = SIMD_SELECT(iconsts.ten.pf,P4,T0);
	    blk->iOpen.pf[i] = iOpen;
	    }
	}
    }
#else
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
    int i,nk;
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
#endif

static void addChild(PKD pkd, CL cl, int iChild, int id, float *fOffset) {
    int idLower, iLower, idUpper, iUpper;
    float cOpen;
    KDN *c;
    int nc = getCell(pkd,iChild,id,&cOpen,&c);
    const BND *cbnd = pkdNodeBnd(pkd,c);

    iLower = c->iLower;
    idLower = id;
    if (c->bRemote) { idUpper = c->pUpper; iUpper = c->pLower; }
    else            { idUpper = id;        iUpper = iLower+1;  }
    clAppend(cl,iChild,idLower,iLower,idUpper,iUpper,nc,cOpen,
	pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
    }
/*
** Returns total number of active particles for which gravity was calculated.
*/
static int processCheckList(PKD pkd, SMX smx, SMF smf, int iRoot, int iVARoot,
    uint8_t uRungLo,uint8_t uRungHi, double dRhoFac, int bEwald, int nGroup,
    double dThetaMin, int bGravStep, double *pdFlop, double *pdPartSum,double *pdCellSum) {
    KDN *k,*c,*kFind;
    int id,idUpper,iCell,iSib,iLower,iUpper,iCheckCell,iCheckLower,iCellDescend;
    PARTICLE *p;
    FMOMR monoPole;
    LOCR L;
    double cx,cy,cz,d2c;
    double fWeight = 0.0;
    double dShiftFlop;
    const vel_t *v;
    const float *a;
    double dOffset[3];
    double xParent,yParent,zParent;
    double d2;
    double dx[3],dir;
    double tax,tay,taz;
    float fOffset[3];
    float dirLsum,normLsum,adotai,maga;
    float fMass,fSoft;
    int iStack;
    int j,jTile,pj,nActive,nTotActive;
    float cOpen,kOpen;
    const BND *cbnd,*kbnd;
    static const float  fZero3[] = {0.0f,0.0f,0.0f};
    static const vel_t vZero3[] = {0.0,0.0,0.0};
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

    pkdGravStartEwald(pkd);

    iStack = -1;

    /*
    ** Clear monopole sentinel and local expansion and the timestepping sums.
    */
    momClearFmomr(&monoPole);
    momClearLocr(&L);
    dirLsum = 0;
    normLsum = 0;
    nTotActive = 0;

    a = fZero3;
    v = vZero3;


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
	    if (kFind->uMinRung>uRungHi || kFind->uMaxRung < uRungLo) {
		/*
		** Move onto processing the sibling.
		*/
		kFind = pkdTreeNode(pkd,++iCellDescend);
	    }
	}
	for (pj=kFind->pLower;pj<=kFind->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    if (!pkdIsActive(pkd,p)) continue;
	    cx = pkdPos(p->r,0);
	    cy = pkdPos(p->r,1);
	    cz = pkdPos(p->r,2);
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
			blk->dx.f[prt] += (float)(cx - pkd->ilp->cx);
			blk->dy.f[prt] += (float)(cy - pkd->ilp->cy);
			blk->dz.f[prt] += (float)(cz - pkd->ilp->cz);
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
			blk->dx.f[prt] += (float)(cx - pkd->ilc->cx);
			blk->dy.f[prt] += (float)(cy - pkd->ilc->cy);
			blk->dz.f[prt] += (float)(cz - pkd->ilc->cz);
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
		maga = sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	    }
	    /*
	    ** For cells which will remain on the checklist for a further deeper level of 
	    ** level of the tree we will be using the stack cl (even if no push happens later).
	    */
	    clClear(pkd->S[iStack+1].cl);
	    do {
		CL_LOOP(pkd->cl,cltile) {
#ifdef USE_SIMD_OPEN
		    iOpenOutcomeSIMD(pkd,k,pkd->cl,cltile,dThetaMin,nGroup);
		    /*Verify:iOpenOutcomeNewCL(pkd,k,pkd->cl,cltile,dThetaMin);*/
#else
		    iOpenOutcomeCL(pkd,k,pkd->cl,cltile,dThetaMin,nGroup);
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
				    id = blk->idLower.i[jTile];
				    pj = -1 - iCheckCell;
				    assert(id >= 0);
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,CID_PARTICLE,pj,id));
				    if (bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
				    ilpAppend(pkd->ilp,
					pkdPos(p->r,0) + blk->xOffset.f[jTile],
					pkdPos(p->r,1) + blk->yOffset.f[jTile],
					pkdPos(p->r,2) + blk->zOffset.f[jTile],
					blk->m.f[jTile], blk->fourh2.f[jTile],
					p->iOrder, v[0], v[1], v[2]);
				    }
				else {
				    id = blk->idLower.i[jTile];
				    assert(id >= 0);
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCheckCell,id));
				    for (pj=c->pLower;pj<=c->pUpper;++pj) {
					if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
					else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,CID_PARTICLE,pj,id));
					fMass = pkdMass(pkd,p);
					fSoft = pkdSoft(pkd,p);
					if (bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
					ilpAppend(pkd->ilp,
					    pkdPos(p->r,0) + blk->xOffset.f[jTile],
					    pkdPos(p->r,1) + blk->yOffset.f[jTile],
					    pkdPos(p->r,2) + blk->zOffset.f[jTile],
					    fMass, 4*fSoft*fSoft,
					    p->iOrder, v[0], v[1], v[2]);
					}
				    }
				break;
			    case 2:
				/*
				** Now I am trying to open a bucket, which means I add each particle back on the
				** checklist with a cell size of zero.
				*/
				iCheckCell = blk->iCell.i[jTile];
				assert(iCheckCell>=0);
				id = blk->idLower.i[jTile];
				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];
				assert(id >= 0);
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCheckCell,id));
				for (pj=c->pLower;pj<=c->pUpper;++pj) {
				    double r[3];
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,CID_PARTICLE,pj,id));
				    for (j=0;j<3;++j) r[j] = pkdPos(p->r,j);
				    fMass = pkdMass(pkd,p);
				    fSoft = pkdSoft(pkd,p);
				    if (bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
				    clAppend(pkd->clNew,-1 - pj,id,0,0,0,1,0.0,fMass,4.0f*fSoft*fSoft,
					r,       /* center of mass */
					fOffset, /* fOffset */
					r,       /* center of box */
					fZero3); /* size of box */
				    }
				break;
			    case 3:
				/*
				** Open the cell.
				** We could do a prefetch here for non-local
				** cells.
				*/
				iCheckCell = blk->iCell.i[jTile];                 assert(iCheckCell >= 0);
				iCheckLower = blk->iLower.i[jTile];               assert(iCheckLower > 0);

				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];

				addChild(pkd,pkd->clNew,blk->iLower.i[jTile],blk->idLower.i[jTile],fOffset);
				addChild(pkd,pkd->clNew,blk->iUpper.i[jTile],blk->idUpper.i[jTile],fOffset);

				break;
			    case 4:
				/*
				** Accept multipole!
				** Interact += Moment(c);
				*/
				iCheckCell = blk->iCell.i[jTile];
				assert(iCheckCell>=0);
				id = blk->idLower.i[jTile];
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCheckCell,id));
				/*
				** Center of mass velocity is used by the planets code to get higher derivatives of the 
				** acceleration and could be used to drift cell moments as an approximation.
				*/
				if (pkd->oNodeVelocity) v = pkdNodeVel(pkd,c);
				ilcAppend(pkd->ilc,
				    blk->x.f[jTile] + blk->xOffset.f[jTile],
				    blk->y.f[jTile] + blk->yOffset.f[jTile],
				    blk->z.f[jTile] + blk->zOffset.f[jTile],
				    pkdNodeMom(pkd,c),c->bMax,v[0],v[1],v[2]);
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
				    if (bGravStep) {
					adotai = a[0]*tax + a[1]*tay + a[2]*taz;
					if (adotai > 0) {
					    adotai /= maga;
					    dirLsum += dir*adotai*adotai;
					    normLsum += adotai*adotai;
					    }
					}
				    }
				else {
				    id = blk->idLower.i[jTile];
				    dOffset[0] = blk->xOffset.f[jTile];
				    dOffset[1] = blk->yOffset.f[jTile];
				    dOffset[2] = blk->zOffset.f[jTile];
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCheckCell,id));
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
				    if (bGravStep) {
					adotai = a[0]*tax + a[1]*tay + a[2]*taz;
					if (adotai > 0) {
					    adotai /= maga;
					    dirLsum += dir*adotai*adotai;
					    normLsum += adotai*adotai;
					    }
					}
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
		assert(pkd->cl!=NULL);
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
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    iCell = k->iLower;
	    nk = getCell(pkd,iCell,pkd->idSelf,&kOpen,&k);
	    kbnd = pkdNodeBnd(pkd,k);
	    /*
	    ** Check iCell is active. We eventually want to just to a
	    ** rung check here when we start using tree repair, but
	    ** for now this is just as good.
	    */
	    if (k->uMinRung<=uRungHi && k->uMaxRung >= uRungLo) {
		/*
		** iCell is active, continue processing it.
		** Put the sibling onto the checklist.
		*/
		iSib = iCell+1;
		nc = getCell(pkd,iSib,pkd->idSelf,&cOpen,&c);
		if (c->uMinRung<=uRungHi && c->uMaxRung >= uRungLo) {
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
		    assert(pkd->cl!=NULL);
		    pkd->S[iStack+1].cl = clTemp;
		    }
		/*
		** Test whether the sibling is active as well.
		** If not we don't push it onto the stack, but we
		** have to be careful to not pop the stack when we
		** hit the sibling.
		*/
		if (c->uMinRung<=uRungHi && c->uMaxRung >= uRungLo) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    pkd->S[iStack].iNodeIndex = iSib;
		    /*
		    ** At this point, the ILP is normally empty if you never do P-P except when reaching a bucket.
		    ** Softened multi-poles are also an exception.
		    */
		    ilpCheckPt(pkd->ilp,&pkd->S[iStack].PartChkPt);
		    ilcCheckPt(pkd->ilc,&pkd->S[iStack].CellChkPt);
		    /*
		    ** Note here we already have the correct elements in S[iStack] (iStack+1 was used previously), just need to add one.
		    */
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
		assert(pkd->cl!=NULL);
		pkd->S[iStack+1].cl = clTemp;
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
	else nActive = 0;
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
	/* Get the next cell to process from the stack */
	if (iStack == -1) goto doneCheckList;
	k = pkdTreeNode(pkd,iCell = pkd->S[iStack].iNodeIndex);
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
	assert(pkd->cl!=NULL);
	pkd->S[iStack].cl = clTemp;
	L = pkd->S[iStack].L;
	dirLsum = pkd->S[iStack].dirLsum;
	normLsum = pkd->S[iStack].normLsum;
	fWeight = pkd->S[iStack].fWeight;
	tempI = *pdFlop;
	tempI += dEwFlop;
	--iStack;
	}
doneCheckList:
    pkdGravFinishEwald(pkd);
    *pdFlop += dEwFlop;   /* Finally add the ewald score to get a proper float count */
    return(nTotActive);
    }

static void doneGravWalk(PKD pkd,SMX smx,SMF *smf) {
    if (smx) {
	smSmoothFinish(smx);
	smFinish(smx,smf);
	}
    }


static void initGravWalk(PKD pkd,double dTime,double dThetaMin,int bPeriodic,int bGravStep,
    SMX *smx, SMF *smf, double *dRhoFac) {
    int pi;

    pkd->dEnergyU = 0.0;
    pkd->fiCritTheta = 1.0f / dThetaMin;

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
int pkdGravWalkHop(PKD pkd,double dTime,int nGroup, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum) {
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
    initGravWalk(pkd,dTime,dThetaMin,0,0,&smx,&smf,&dRhoFac);
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
	    id = pkd->hopRoots[i].iPid;
	    iRoot = pkd->hopRoots[i].iIndex;
	    assert(iRoot>0);
	    addChild(pkd,pkd->cl,iRoot,id,fOffset);
	    }
	assert(pkd->hopRoots[iRootSelf].iPid==pkd->idSelf);
	nActive += processCheckList(pkd, smx, smf, pkd->hopRoots[iRootSelf].iIndex, 0, 0, MAX_RUNG,
	    dRhoFac, 0, nGroup, dThetaMin, 0, pdFlop, pdPartSum, pdCellSum);
	}
    doneGravWalk(pkd,smx,&smf);
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    return nActive;
}


/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,int nGroup,
    int iRoot1, int iRoot2,int iVARoot,
    double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    KDN *c;
    int id,iLower;
    int iLocalRoot;
    float fOffset[3];
    int ix,iy,iz,bRep;
    float cOpen;
    const BND *cbnd;
    int nc;
    int nActive = 0;
    double dRhoFac;
    SMX smx;
    SMF smf;

    initGravWalk(pkd,dTime,dThetaMin,nReps?1:0,pkd->param.bGravStep,&smx,&smf,&dRhoFac);

    /*
    ** Walk tree 1 against trees 1 (and optionally 2).
    ** Go to the local tree of iRoot1. This is the tree we will walk. We walk against
    ** the tree at iRoot1, and optionally the tree at iRoot2.
    */
    c = pkdTreeNode(pkd,iLocalRoot=iRoot1);
    while(c->bRemote) c = pkdTreeNode(pkd,iLocalRoot = c->iLower);
    /* Check that the iRoot has active particles! */
    if (pkdIsCellActive(c,uRungLo,uRungHi)) {
	/*
	** Initially we set our cell pointer to
	** point to the top tree.
	*/
	ilpClear(pkd->ilp);
	ilcClear(pkd->ilc);
	clClear(pkd->cl);

	/*
	** Add all replicas of the entire box to the Checklist.
	** We add at least one box (0,0,0). The root cell is alway on processor 0.
	*/
	for (ix=-nReps;ix<=nReps;++ix) {
	    fOffset[0] = ix*pkd->fPeriod[0];
	    for (iy=-nReps;iy<=nReps;++iy) {
		fOffset[1] = iy*pkd->fPeriod[1];
		for (iz=-nReps;iz<=nReps;++iz) {
		    fOffset[2] = iz*pkd->fPeriod[2];
		    bRep = ix || iy || iz;
		    /* 
		    ** Use leaf of the top tree and NOT the root of the local tree here.
		    */
		    id = 0;
		    addChild(pkd,pkd->cl,iRoot1,id,fOffset);
		    if (iRoot2>0) addChild(pkd,pkd->cl,iRoot2,id,fOffset);
		    }
		}
	    }
	nActive += processCheckList(pkd, smx, smf, iLocalRoot, iVARoot, uRungLo, uRungHi,
	    dRhoFac, bEwald, nGroup, dThetaMin, pkd->param.bGravStep, pdFlop, pdPartSum, pdCellSum);
	}
    /*
    ** Walk tree 2 against tree 1.
    */
    if (iRoot2>0) {
	c = pkdTreeNode(pkd,iLocalRoot=iRoot2);
	while(c->bRemote) c = pkdTreeNode(pkd,iLocalRoot = c->iLower);
	/* Check that the iRoot has active particles! */
	if (!pkdIsCellActive(pkdTreeNode(pkd,iLocalRoot),uRungLo,uRungHi)) return 0;
	ilpClear(pkd->ilp);
	ilcClear(pkd->ilc);
	clClear(pkd->cl);
	/*
	** Add all replicas of the entire box to the Checklist.
	** We add at least one box (0,0,0). The root cell is alway on processor 0.
	*/
	for (ix=-nReps;ix<=nReps;++ix) {
	    fOffset[0] = ix*pkd->fPeriod[0];
	    for (iy=-nReps;iy<=nReps;++iy) {
		fOffset[1] = iy*pkd->fPeriod[1];
		for (iz=-nReps;iz<=nReps;++iz) {
		    fOffset[2] = iz*pkd->fPeriod[2];
		    bRep = ix || iy || iz;
		    /* 
		    ** Use leaf of the top tree and NOT the root of the local tree here.
		    */
		    id = 0;
		    addChild(pkd,pkd->cl,iRoot1,id,fOffset);
		    }
		}
	    }
	nActive += processCheckList(pkd, smx, smf, iLocalRoot, iVARoot, uRungLo, uRungHi,
	    dRhoFac, 0, nGroup, dThetaMin, pkd->param.bGravStep, pdFlop, pdPartSum, pdCellSum);
	}
    doneGravWalk(pkd,smx,&smf);
    return nActive;
    }

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalkGroups(PKD pkd,double dTime,int nGroup, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    KDN *c;
    int id,iRoot;
    float fOffset[3];
    int i,j,k;
    float cOpen;
    const BND *cbnd;
    int nc;
    double dRhoFac;
    SMX smx;
    SMF smf;

    initGravWalk(pkd,dTime,dThetaMin,0,0,&smx,&smf,&dRhoFac);
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
	    id = gd[i].treeRoots[k].iPid;
	    iRoot = gd[i].treeRoots[k].iLocalRootId;
	    addChild(pkd,pkd->cl,iRoot,id,fOffset);
        }
#endif
	nActive += processCheckList(pkd, smx, smf, gd[i].treeRoots[0].iLocalRootId, 0, 0, MAX_RUNG,
	    dRhoFac, 0, nGroup, dThetaMin, 0, pdFlop, pdPartSum, pdCellSum);
    }
    doneGravWalk(pkd,smx,&smf);
    return nActive;
}
