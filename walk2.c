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
#include "moments.h"
#include "vmoments.h"
#include "cl.h"
#include "cudautil.h"
#include "cudapppc.h"
#include "cudaewald.h"
#include "clutil.h"

static inline int getCell(PKD pkd,int iCache,int iCell,int id,float *pcOpen,KDN **pc) {
    KDN *c;
    int nc;
    assert(iCell > 0);
    assert(id >= 0);
    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCell);
    else c = CAST(KDN *,mdlFetch(pkd->mdl,iCache,iCell,id));
    *pc = c;
    if (c->bRemote|c->bTopTree) nc = 1000000000; /* we never allow pp with this cell */
    else nc = c->pUpper - c->pLower + 1;
    *pcOpen = c->bMax * pkd->fiCritTheta;
    return nc;
    }


#ifdef USE_SIMD_OPEN
void iOpenOutcomeSIMD(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin,SPHOptions *SPHoptions);
#endif
#if 1
/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version has been changed by adding the ability to open buckets.
*/
static void iOpenOutcomeCL(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin,SPHOptions *SPHoptions) {
    const int walk_min_multipole = 3;
    float dx,dy,dz,mink2,d2,d2Open,xc,yc,zc,fourh2,minbnd2,kOpen,cOpen,diCrit,distk2,distc2;
    int i;
    int iOpen,iOpenA,iOpenB;
    CL_BLK *blk;
    int n, nLeft;
    BND kbnd;
    double k_r[3];

    distk2 = HUGE_VALF;
    distc2 = HUGE_VALF;

    kbnd = pkdNodeGetBnd(pkd,k);
    pkdNodeGetPos(pkd,k,k_r);

    diCrit = 1.0f / dThetaMin;
    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	n = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	for(i=0; i<n; ++i) {
	    if (blk->m.f[i] <= 0) iOpen = 10;  /* ignore this cell */
	    else {
		fourh2 = blk->fourh2.f[i];

        if (SPHoptions->doDensity || SPHoptions->doSPHForces) {
            distk2 = 0.0f;
            dx = kbnd.fCenter[0] -  blk->xCenter.f[i] - blk->xOffset.f[i] - blk->xMax.f[i];
            if (dx > 0) distk2 += dx*dx;
            dx = blk->xCenter.f[i] + blk->xOffset.f[i] - blk->xMax.f[i] - kbnd.fCenter[0];
            if (dx > 0) distk2 += dx*dx;

            dx = kbnd.fCenter[1] - blk->yCenter.f[i] - blk->yOffset.f[i] - blk->yMax.f[i];
            if (dx > 0) distk2 += dx*dx;
            dx = blk->yCenter.f[i] + blk->yOffset.f[i] - blk->yMax.f[i] - kbnd.fCenter[1];
            if (dx > 0) distk2 += dx*dx;

            dx = kbnd.fCenter[2] - blk->zCenter.f[i] - blk->zOffset.f[i] - blk->zMax.f[i];
            if (dx > 0) distk2 += dx*dx;
            dx = blk->zCenter.f[i] + blk->zOffset.f[i] - blk->zMax.f[i] - kbnd.fCenter[2];
            if (dx > 0) distk2 += dx*dx;
        }
        if (SPHoptions->doSPHForces) {
            distc2 = 0.0f;
            dx = kbnd.fCenter[0] - kbnd.fMax[0] -  blk->xCenter.f[i] - blk->xOffset.f[i];
            if (dx > 0) distc2 += dx*dx;
            dx = blk->xCenter.f[i] + blk->xOffset.f[i] - kbnd.fCenter[0] - kbnd.fMax[0];
            if (dx > 0) distc2 += dx*dx;

            dx = kbnd.fCenter[1] - kbnd.fMax[1] - blk->yCenter.f[i] - blk->yOffset.f[i];
            if (dx > 0) distc2 += dx*dx;
            dx = blk->yCenter.f[i] + blk->yOffset.f[i] - kbnd.fCenter[1] - kbnd.fMax[1];
            if (dx > 0) distc2 += dx*dx;

            dx = kbnd.fCenter[2] - kbnd.fMax[2] - blk->zCenter.f[i] - blk->zOffset.f[i];
            if (dx > 0) distc2 += dx*dx;
            dx = blk->zCenter.f[i] + blk->zOffset.f[i] - kbnd.fCenter[2] - kbnd.fMax[2];
            if (dx > 0) distc2 += dx*dx;
        }

		xc = blk->x.f[i] + blk->xOffset.f[i];
		yc = blk->y.f[i] + blk->yOffset.f[i];
		zc = blk->z.f[i] + blk->zOffset.f[i];
		d2 = pow(k_r[0]-xc,2) + pow(k_r[1]-yc,2) + pow(k_r[2]-zc,2);
		kOpen = 1.5f * k->bMax * diCrit;
		cOpen = blk->cOpen.f[i];
		d2Open = pow(cOpen+kOpen,2);
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

		if (d2 > d2Open && minbnd2 > fourh2 && distk2 > k->fBoBr2 && distc2 > blk->fBoBr2.f[i]) iOpen = 8;
		else if (cOpen > kOpen) {
		    if (blk->iLower.i[i]) iOpen = 3;
		    else iOpen = 2;
		    }
		else {
		    if (blk->iLower.i[i] == 0) iOpenA = 1;
		    else iOpenA = 3;
		    if (blk->nc.i[i] < walk_min_multipole || mink2 <= cOpen*cOpen) iOpenB = iOpenA;
		    else if (minbnd2 > fourh2 && distk2 > k->fBoBr2 && distc2 > blk->fBoBr2.f[i]) iOpenB = 4;
		    else iOpenB = iOpenA;
		    if (!k->bGroup) iOpen = 0;
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
    double dFlop = COST_FLOP_OPEN*(tile->lstTile.nBlocks*CL_PART_PER_BLK  + tile->lstTile.nInLast);
    pkd->dFlop += dFlop;
    pkd->dFlopSingleCPU += dFlop;
    }
#endif

static void addChild(PKD pkd, int iCache, CL cl, int iChild, int id, float *fOffset) {
    int idLower, iLower, idUpper, iUpper;
    float cOpen;
    KDN *c;
    double c_r[3];
    int nc = getCell(pkd,iCache,iChild,id,&cOpen,&c);
    pkdNodeGetPos(pkd,c,c_r);
    BND cbnd = pkdNodeGetBnd(pkd,c);
    pkdGetChildCells(c,id,idLower,iLower,idUpper,iUpper);
    clAppend(cl,iCache,id,iChild,idLower,iLower,idUpper,iUpper,nc,cOpen,
	pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c_r,fOffset,cbnd.fCenter,cbnd.fMax,c->fBoBr2);
    }
/*
** Returns total number of active particles for which gravity was calculated.
*/
static int processCheckList(PKD pkd, SMX smx, SMF smf, int iRoot, int iRoot2, 
    struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
    double dTime,int bEwald,
    double dThetaMin, double *pdFlop, double *pdPartSum,double *pdCellSum,SPHOptions *SPHoptions) {
    KDN *k,*c,*kFind;
    int id,idUpper,iCell,iSib,iLower,iUpper,iCheckCell,iCheckLower,iCellDescend;
    PARTICLE *p;
    FMOMR monoPole;
    LOCR L;
    FLOCR Lf;
    double cx,cy,cz,d2c;
    double dShiftFlop;
    const vel_t *v;
    const float *a;
    double r[3], k_r[3], c_r[3];
    double dOffset[3];
    double xParent,yParent,zParent;
    double d2;
    double dx[3],dir;
    double tax,tay,taz;
    float fOffset[3];
    float dirLsum,normLsum,adotai;
    float fMass,fSoft;
    int iStack;
    int j,jTile,pj,nActive,nTotActive;
    float cOpen,kOpen;
    static const float  fZero3[] = {0.0f,0.0f,0.0f};
    static const vel_t vZero3[] = {0.0,0.0,0.0};
    ILPTILE tile;
    ILCTILE ctile;
    CL clTemp;
    CLTILE cltile;
    int iCidPart, iCidCell;
    int bReferenceFound;
    uint64_t iOrder;
#ifdef USE_SIMD_FMM
    float fzero[3] = {0,0,0};
#else
    float imaga;
#endif
    double dFlop;

    pkd->dFlop = 0.0; /* Flops are accumulated here! */
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
    k = pkdTreeNode(pkd,iCell = iRoot);
    pkdNodeGetPos(pkd,k,k_r);
    while (1) {
#ifdef ILP_ILC_CAN_BE_NON_EMPTY
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
	    pkdGetPos3(pkd,p,cx,cy,cz);
	    goto found_it;
	}
	printf("%d: TREE ERROR\n", pkd->idSelf);
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
	bReferenceFound = 1;
#else
	bReferenceFound = 0;
#endif
	while (1) {
	    /*
	    ** Process the Checklist.
	    */
#ifdef USE_SIMD_FMM
	    if (ts->bGravStep) a = pkdNodeAccel(pkd,k);
	    else a = fzero;
#else
	    if (ts->bGravStep) {
		a = pkdNodeAccel(pkd,k);
		imaga = 1.0 / sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	    }
#endif
	    /*
	    ** For cells which will remain on the checklist for a further deeper level of 
	    ** level of the tree we will be using the stack cl (even if no push happens later).
	    */
	    clClear(pkd->S[iStack+1].cl);
#ifdef USE_SIMD_FMM
	    ilcClear(pkd->ill);
#endif
	    do {
		CL_LOOP(pkd->cl,cltile) {
#ifdef USE_SIMD_OPEN
		    iOpenOutcomeSIMD(pkd,k,pkd->cl,cltile,dThetaMin,SPHoptions);
		    /*Verify:iOpenOutcomeNewCL(pkd,k,pkd->cl,cltile,dThetaMin);*/
#else
		    iOpenOutcomeCL(pkd,k,pkd->cl,cltile,dThetaMin,SPHoptions);
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
				id = blk->idCell.i[jTile];
				if (iCheckCell < 0) {
				    pj = -1 - iCheckCell;
				    assert(id >= 0);
				    iCidPart = blk->iCache.i[jTile]==CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,iCidPart,pj,id));
				    if (ts->bGravStep && ts->iTimeStepCrit == 1) v = pkdVel(pkd,p);
				    pkdGetPos1(pkd,p,r);
				    if (!bReferenceFound) {
					bReferenceFound=1;
					pkd->ilp->cx=r[0]; pkd->ilp->cy=r[1]; pkd->ilp->cz=r[2];
					pkd->ilc->cx=r[0]; pkd->ilc->cy=r[1]; pkd->ilc->cz=r[2];
					}
				    iOrder = pkd->bNoParticleOrder ? 0 : p->iOrder;
                    v = pkdVel(pkd,p);
                    NEWSPHFIELDS *pNewSph = pkdNewSph(pkd,p);
                    float dtPredDrift = getDtPredDrift(kick,p->bMarked,ts->uRungLo,p->uRung);
                    float Omega = pNewSph->Omega;                     /* should be the Omega field of the sph fields, nyi */
                    float P = 0.0f;                         /* should be calculated by the EOS, nyi */
                    float cs = 0.0f;                        /* should be calculated by the EOS, nyi */
                    const float *ap = pkdAccel(pkd,p);
                    if (SPHoptions->doSPHForces) {
                        P = EOSPCofRhoU(pkdDensity(pkd,p),pNewSph->u + dtPredDrift * pNewSph->uDot,&cs,SPHoptions);
                    }
				    ilpAppend(pkd->ilp,
					r[0] + blk->xOffset.f[jTile],
					r[1] + blk->yOffset.f[jTile],
					r[2] + blk->zOffset.f[jTile],
					blk->m.f[jTile], blk->fourh2.f[jTile],
					iOrder, v[0] + dtPredDrift * ap[0], v[1] + dtPredDrift * ap[1], v[2] + dtPredDrift * ap[2],
                    pkdBall(pkd,p), Omega, pkdDensity(pkd,p), P, cs, pkdSpecies(pkd,p));
				    }
				else {
				    assert(id >= 0);
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else {
					c = CAST(KDN *,mdlFetch(pkd->mdl,blk->iCache.i[jTile],iCheckCell,id));
					}
				    iCidPart = blk->iCache.i[jTile]==CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
				    if (!bReferenceFound) {
					bReferenceFound=1;
					if (id == pkd->idSelf) p = pkdParticle(pkd,c->pLower);
					else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,iCidPart,c->pLower,id));
					pkdGetPos1(pkd,p,r);
					pkd->ilp->cx=r[0]; pkd->ilp->cy=r[1]; pkd->ilp->cz=r[2];
					pkd->ilc->cx=r[0]; pkd->ilc->cy=r[1]; pkd->ilc->cz=r[2];
					}
				    for (pj=c->pLower;pj<=c->pUpper;++pj) {
					if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
					else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,iCidPart,pj,id));
					fMass = pkdMass(pkd,p);
					fSoft = pkdSoft(pkd,p);
					if (ts->bGravStep && ts->iTimeStepCrit == 1) v = pkdVel(pkd,p);
					pkdGetPos1(pkd,p,r);
					iOrder = pkd->bNoParticleOrder ? 0 : p->iOrder;
                    v = pkdVel(pkd,p);
                    NEWSPHFIELDS *pNewSph = pkdNewSph(pkd,p);
                    float dtPredDrift = getDtPredDrift(kick,p->bMarked,ts->uRungLo,p->uRung);
                    float Omega = pNewSph->Omega;                 /* should be the Omega field of the sph fields, nyi */
                    float P = 0.0f;                     /* should be calculated by the EOS, nyi */
                    float cs = 0.0f;                    /* should be calculated by the EOS, nyi */
                    const float *ap = pkdAccel(pkd,p);
                    if (SPHoptions->doSPHForces) {
                        P = EOSPCofRhoU(pkdDensity(pkd,p),pNewSph->u + dtPredDrift * pNewSph->uDot,&cs,SPHoptions);
                    }
					ilpAppend(pkd->ilp,
					    r[0] + blk->xOffset.f[jTile],
					    r[1] + blk->yOffset.f[jTile],
					    r[2] + blk->zOffset.f[jTile],
					    fMass, 4*fSoft*fSoft,
					    iOrder, v[0] + dtPredDrift * ap[0], v[1] + dtPredDrift * ap[1], v[2] + dtPredDrift * ap[2],
                        pkdBall(pkd,p), Omega, pkdDensity(pkd,p), P, cs, pkdSpecies(pkd,p));
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
				id = blk->idCell.i[jTile];
				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];
				assert(id >= 0);
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else {
				    c = CAST(KDN *,mdlFetch(pkd->mdl,blk->iCache.i[jTile],iCheckCell,id));
				    }
				iCidPart = blk->iCache.i[jTile]==CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
				for (pj=c->pLower;pj<=c->pUpper;++pj) {
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,iCidPart,pj,id));
				    pkdGetPos1(pkd,p,r);
				    fMass = pkdMass(pkd,p);
				    fSoft = pkdSoft(pkd,p);
				    if (ts->bGravStep && ts->iTimeStepCrit == 1) v = pkdVel(pkd,p);
				    clAppend(pkd->clNew,blk->iCache.i[jTile],id,-1 - pj,0,0,0,0,1,0.0,fMass,4.0f*fSoft*fSoft,
					r,       /* center of mass */
					fOffset, /* fOffset */
					r,       /* center of box */
					fZero3,  /* size of box */
                    blk->fBoBr2.f[jTile]);
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

				addChild(pkd,blk->iCache.i[jTile],pkd->clNew,blk->iLower.i[jTile],blk->idLower.i[jTile],fOffset);
				addChild(pkd,blk->iCache.i[jTile],pkd->clNew,blk->iUpper.i[jTile],blk->idUpper.i[jTile],fOffset);

				break;
			    case 4:
				/*
				** Accept multipole!
				** Interact += Moment(c);
				*/
                if (SPHoptions->doGravity) {
				iCheckCell = blk->iCell.i[jTile];
				assert(iCheckCell>=0);
				id = blk->idCell.i[jTile];
				if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				else {
				    c = CAST(KDN *,mdlFetch(pkd->mdl,blk->iCache.i[jTile],iCheckCell,id));
				    }
				r[0] = blk->x.f[jTile] + blk->xOffset.f[jTile];
				r[1] = blk->y.f[jTile] + blk->yOffset.f[jTile];
				r[2] = blk->z.f[jTile] + blk->zOffset.f[jTile];
				if (!bReferenceFound) {
				    bReferenceFound=1;
				    pkd->ilp->cx=r[0]; pkd->ilp->cy=r[1]; pkd->ilp->cz=r[2];
				    pkd->ilc->cx=r[0]; pkd->ilc->cy=r[1]; pkd->ilc->cz=r[2];
				    }
				ilcAppend(pkd->ilc,r[0],r[1],r[2],pkdNodeMom(pkd,c),c->bMax);
                }
				break;
			    case 8:
				/*
				** Local expansion accepted!
				*/
                if (SPHoptions->doGravity) {
				iCheckCell = blk->iCell.i[jTile];
				if (iCheckCell<0) {
				    fOffset[0] = blk->xOffset.f[jTile];
				    fOffset[1] = blk->yOffset.f[jTile];
				    fOffset[2] = blk->zOffset.f[jTile];
				    dx[0] = k_r[0] - (blk->x.f[jTile] + blk->xOffset.f[jTile]);
				    dx[1] = k_r[1] - (blk->y.f[jTile] + blk->yOffset.f[jTile]);
				    dx[2] = k_r[2] - (blk->z.f[jTile] + blk->zOffset.f[jTile]);
				    d2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
				    dir = 1.0/sqrt(d2);
#ifdef USE_SIMD_FMM
				    monoPole.m = blk->m.f[jTile];
				    ilcAppendFloat(pkd->ill,dx[0],dx[1],dx[2],&monoPole,1.0);
#else
				    /* monoPole.m = blk->m.f[jTile];*/
				    /* *pdFlop += momLocrAddFmomr5cm(&L,&monoPole,0.0,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);*/
				    dFlop = momLocrAddMono5(&L,blk->m.f[jTile],dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
				    *pdFlop += dFlop;
				    pkd->dFlopDoubleCPU += dFlop;
				    if (ts->bGravStep) {
					adotai = a[0]*tax + a[1]*tay + a[2]*taz;
					if (adotai > 0) {
					    adotai *= imaga;
					    dirLsum += dir*adotai*adotai;
					    normLsum += adotai*adotai;
					    }
					}
#endif
				    }
				else {
				    id = blk->idCell.i[jTile];
				    dOffset[0] = blk->xOffset.f[jTile];
				    dOffset[1] = blk->yOffset.f[jTile];
				    dOffset[2] = blk->zOffset.f[jTile];
				    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCheckCell);
				    else {
					c = CAST(KDN *,mdlFetch(pkd->mdl,blk->iCache.i[jTile],iCheckCell,id));
					}
				    pkdNodeGetPos(pkd,c,c_r);
#ifdef USE_SIMD_FMM
				    for (j=0;j<3;++j) dx[j] = k_r[j] - (c_r[j] + dOffset[j]);
				    ilcAppendFloat(pkd->ill,dx[0],dx[1],dx[2],pkdNodeMom(pkd,c),c->bMax);
#else
				    d2 = 0;
				    for (j=0;j<3;++j) {
					dx[j] = k_r[j] - (c_r[j] + dOffset[j]);
					d2 += dx[j]*dx[j];
					}
				    dir = 1.0/sqrt(d2);
				    dFlop = momLocrAddFmomr5cm(&L,pkdNodeMom(pkd,c),c->bMax,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
				    *pdFlop += dFlop;
				    pkd->dFlopDoubleCPU += dFlop;
				    if (ts->bGravStep) {
					adotai = a[0]*tax + a[1]*tay + a[2]*taz;
					if (adotai > 0) {
					    adotai *= imaga;
					    dirLsum += dir*adotai*adotai;
					    normLsum += adotai*adotai;
					    }
					}
#endif
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
	    ** Now calculate the local expansion.
	    */
#ifdef USE_SIMD_FMM
	    // Need to do something here with dirLsum and normLsum for GravStep
	    // Need to get the scaling factor correct here
	    if (ilcCount(pkd->ill)) {
		float v = k->bMax;
		dFlop = momFlocrSetVFmomr5cm(&Lf,v,pkd->ill,a,&dirLsum,&normLsum);
		*pdFlop += dFlop;
		pkd->dFlopSingleCPU += dFlop;
		momLocrAddFlocr(&L,&Lf,v);
		}
#endif
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (k->bGroup) break; /* A bucket is ALWAYS a group */
	    xParent = k_r[0];
	    yParent = k_r[1];
	    zParent = k_r[2];
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    iCell = k->iLower;
	    getCell(pkd,-1,iCell,pkd->idSelf,&kOpen,&k);
	    pkdNodeGetPos(pkd,k,k_r);
	    /*
	    ** Check iCell is active. We eventually want to just to a
	    ** rung check here when we start using tree repair, but
	    ** for now this is just as good.
	    */
	    if (k->uMinRung<=ts->uRungHi && k->uMaxRung >= ts->uRungLo) {
		/*
		** iCell is active, continue processing it.
		** Put the sibling onto the checklist.
		*/
		iSib = iCell+1;
		getCell(pkd,-1,iSib,pkd->idSelf,&cOpen,&c);
		if (c->uMinRung<=ts->uRungHi && c->uMaxRung >= ts->uRungLo) {
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
		if (c->uMinRung<=ts->uRungHi && c->uMaxRung >= ts->uRungLo) {
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
		    pkdNodeGetPos(pkd,c,c_r);
		    dShiftFlop = momShiftLocr(&pkd->S[iStack].L,
					      c_r[0] - xParent,
					      c_r[1] - yParent,
					      c_r[2] - zParent);
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
		pkdNodeGetPos(pkd,k,k_r);
		}
	    dFlop = momShiftLocr(&L,k_r[0] - xParent,
				    k_r[1] - yParent,
				    k_r[2] - zParent);
	    *pdFlop += dFlop;
	    pkd->dFlopDoubleCPU += dFlop;
	    }
	/*
	** Now the interaction list should be complete and the
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
	nActive = pkdGravInteract(pkd,kick,lc,ts,
	    k,&L,pkd->ilp,pkd->ilc,dirLsum,normLsum,bEwald,pdFlop,
	    smx, &smf, iRoot, iRoot2, SPHoptions);
	/*
	** Update the limit for a shift of the center here based on the opening radius of this
	** cell (the one we just evaluated).
	*/
	if (nActive) {
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
	pkdNodeGetPos(pkd,k,k_r);
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
	--iStack;
	}
doneCheckList:
#ifdef USE_CUDA
    CudaClientFlush(pkd->cudaClient);
#endif
    mdlCompleteAllWork(pkd->mdl);
    *pdFlop += pkd->dFlop; /* Accumulate work flops (notably Ewald) */
    return(nTotActive);
    }

static void doneGravWalk(PKD pkd,SMX smx,SMF *smf) {
    if (smx) {
	smSmoothFinish(smx);
	smFinish(smx,smf);
	}
    }

static void initGravWalk(PKD pkd,double dTime,double dThetaMin,int bPeriodic,int bGravStep,int nPartRhoLoc,int iTimeStepCrit,
    SMX *smx, SMF *smf) {
    int pi;

    pkd->dEnergyU = 0.0;
    pkd->dEnergyT = 0.0;
    pkd->dEnergyW = 0.0;
    pkd->dEnergyF[0] = pkd->dEnergyF[1] = pkd->dEnergyF[2] = 0.0;
    pkd->dEnergyL[0] = pkd->dEnergyL[1] = pkd->dEnergyL[2] = 0.0;

    pkd->fiCritTheta = 1.0f / dThetaMin;

    assert(pkd->oNodeMom);
    if (bGravStep) {
	assert(pkd->oNodeAcceleration);
	if (iTimeStepCrit == 1) {
	    assert(pkd->oNodeVelocity);
	    assert(pkd->oFieldOffset[oVelocity]);
	    }
	}

    /*
    ** Setup smooth for calculating local densities when a particle has too few P-P interactions.
    */
    if (bGravStep) {
	smInitializeRO(smx,pkd,smf,nPartRhoLoc,bPeriodic,SMX_DENSITY_F1);
	smSmoothInitialize(*smx);
	}
    else (*smx) = NULL;
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
    SMX smx;
    SMF smf;

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd), pkdLocal(pkd));
    initGravWalk(pkd,dTime,dThetaMin,0,0,0,0,&smx,&smf);
    nActive = 0;
    for(gid=1; gid<pkd->nGroups; ++gid) {
	if (!pkd->hopGroups[gid].bNeedGrav) continue;
	pkd->hopGroups[gid].bNeedGrav = 0;
	ilpClear(pkd->ilp);
	ilcClear(pkd->ilc);
	clClear(pkd->cl);
	iRootSelf = pkd->hopGroups[gid].iAllRoots;
	for (i=iRootSelf; i<iRootSelf + pkd->hopGroups[gid].nRemote+1; ++i) {
	    for (j=0;j<3;++j) fOffset[j] = 0.0f;
	    id = pkd->hopRoots[i].iPid;
	    iRoot = pkd->hopRoots[i].iIndex;
	    assert(iRoot>0);
	    addChild(pkd,CID_CELL,pkd->cl,iRoot,id,fOffset);
	    }
	assert(pkd->hopRoots[iRootSelf].iPid==pkd->idSelf);
	// nActive += processCheckList(pkd, smx, smf, pkd->hopRoots[iRootSelf].iIndex, 0, 0, MAX_RUNG,
	//     NULL,NULL,1.0,dTime,
	//     0, dThetaMin, 0, 0, pdFlop, pdPartSum, pdCellSum);
	}
    doneGravWalk(pkd,smx,&smf);
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    return nActive;
}


/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
    double dTime,int nReps,int bEwald,int nGroup,
    int iLocalRoot1, int iLocalRoot2,int iVARoot,
    double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum,SPHOptions *SPHoptions) {
    int id;
    float fOffset[3];
    int ix,iy,iz,bRep;
    float cOpen;
    const BND *cbnd;
    int nc;
    int nActive = 0;
    SMX smx;
    SMF smf;
    int iTop1, iTop2;

    initGravWalk(pkd,dTime,dThetaMin,nReps?1:0,ts->bGravStep,ts->nPartRhoLoc,ts->iTimeStepCrit,&smx,&smf);

    iTop1 = pkd->iTopTree[iLocalRoot1];
    iTop2 = pkd->iTopTree[iLocalRoot2];
    id = pkd->idSelf;

    /*
    ** Walk tree 1 against trees 1 (and optionally 2) if there are active particles
    */
    KDN *k = pkdTreeNode(pkd,iLocalRoot1);
    if (k->pLower<=k->pUpper && pkdIsCellActive(k,ts->uRungLo,ts->uRungHi)) {
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
		    addChild(pkd,CID_CELL,pkd->cl,iTop1,id,fOffset);
#ifndef SINGLE_CACHES
		    if (iLocalRoot2>0) addChild(pkd,CID_CELL2,pkd->cl,iTop2,id,fOffset);
#else
		    if (iLocalRoot2>0) addChild(pkd,CID_CELL,pkd->cl,iTop2,id,fOffset);
#endif
		    }
		}
	    }
	nActive += processCheckList(pkd, smx, smf, iLocalRoot1, iLocalRoot2, kick,lc,ts,
	    dTime,bEwald, dThetaMin, pdFlop, pdPartSum, pdCellSum, SPHoptions);
	}
#if 0
    /*
    ** Walk tree 2 against tree 1.
    */
    if (iLocalRoot2>0) {
	/* Check that the iRoot has active particles! */
	if (!pkdIsCellActive(pkdTreeNode(pkd,iLocalRoot2),uRungLo,uRungHi)) return 0;
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
		    addChild(pkd,CID_CELL,pkd->cl,iTop1,id,fOffset);
		    }
		}
	    }
	nActive += processCheckList(pkd, smx, smf, iLocalRoot2, iLocalRoot1,
	    kick,lc,ts,dTime,0, dThetaMin, pdFlop, pdPartSum, pdCellSum, SPHoptions);
	}
#endif
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
    SMX smx;
    SMF smf;

    initGravWalk(pkd,dTime,dThetaMin,0,0,0,0,&smx,&smf);
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
	    addChild(pkd,CID_CELL,pkd->cl,iRoot,id,fOffset);
        }
#endif
	// nActive += processCheckList(pkd, smx, smf, gd[i].treeRoots[0].iLocalRootId, 0, 0, MAX_RUNG,
	//     NULL,NULL,1.0,dTime,
	//     0, dThetaMin, 0, 0, pdFlop, pdPartSum, pdCellSum);
    }
    doneGravWalk(pkd,smx,&smf);
    return nActive;
}
