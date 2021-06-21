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
#define MPICH_SKIP_MPICXX
#include "simd.h"
#include "pkd.h"

/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version will also open buckets ("new" criteria)
*/
extern "C"
void iOpenOutcomeSIMD(PKD pkd,KDN *k,CL cl,CLTILE tile,float dThetaMin,SPHOptions SPHoptions) {
    const float walk_min_multipole = 3;
    fmask T0,T1,T2,T3,T4,T6,T7;
    i32v P1,P2,P3,P4;
    fvec xc,yc,zc,dx,dy,dz,d2,diCrit,cOpen,cOpen2,d2Open,mink2,minbnd2,fourh2;
    int i,iEnd,nLeft;
    CL_BLK *blk;
    i32v iOpen,iOpenA,iOpenB;
    BND kbnd;
    fvec k_xCenter, k_yCenter, k_zCenter, k_xMax, k_yMax, k_zMax;
    fvec k_xMinBnd, k_yMinBnd, k_zMinBnd, k_xMaxBnd, k_yMaxBnd, k_zMaxBnd;
    fvec k_x, k_y, k_z, k_bMax, k_Open;
    fvec k_fBoBr2, blk_fBoBr2;
    fmask k_notgrp;
    double k_r[3];
    pkdNodeGetPos(pkd,k,k_r);

    assert ( pkdNodeMom(pkd,k)->m > 0.0f );

    diCrit = 1.0f/dThetaMin;

    kbnd = pkdNodeGetBnd(pkd,k);
    k_xMinBnd = kbnd.fCenter[0]-kbnd.fMax[0];
    k_yMinBnd = kbnd.fCenter[1]-kbnd.fMax[1];
    k_zMinBnd = kbnd.fCenter[2]-kbnd.fMax[2];
    k_xMaxBnd = kbnd.fCenter[0]+kbnd.fMax[0];
    k_yMaxBnd = kbnd.fCenter[1]+kbnd.fMax[1];
    k_zMaxBnd = kbnd.fCenter[2]+kbnd.fMax[2];
    k_xCenter = kbnd.fCenter[0];
    k_yCenter = kbnd.fCenter[1];
    k_zCenter = kbnd.fCenter[2];
    k_xMax = kbnd.fMax[0];
    k_yMax = kbnd.fMax[1];
    k_zMax = kbnd.fMax[2];
    k_x = k_r[0];
    k_y = k_r[1];
    k_z = k_r[2];
    k_bMax = k->bMax;
    k_notgrp = cvt_fvec(i32v(k->bGroup)) == 0.0;
    k_Open = 1.5f*k_bMax*diCrit;
    k_fBoBr2 = k->fBoBr2;

    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	iEnd = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	iEnd = (iEnd+fvec::mask()) >> SIMD_BITS;
	for(i=0; i<iEnd; ++i) {
	    fourh2 = blk->fourh2.p[i];

        if (SPHoptions.doDensity || SPHoptions.doSPHforces) {
            T0 = fourh2 > k_fBoBr2;
            fourh2 = mask_mov(k_fBoBr2,T0,fourh2);
        }
        if (SPHoptions.doSPHforces) {
            blk_fBoBr2 = blk->fBoBr2.p[i];
            T0 = fourh2 > blk_fBoBr2;
            fourh2 = mask_mov(blk_fBoBr2,T0,fourh2);
        }

	    xc = fvec(blk->x.p[i]) + fvec(blk->xOffset.p[i]);
	    yc = fvec(blk->y.p[i]) + fvec(blk->yOffset.p[i]);
	    zc = fvec(blk->z.p[i]) + fvec(blk->zOffset.p[i]);
	    dx = k_x - xc;
	    dy = k_y - yc;
	    dz = k_z - zc;
	    d2 = dx*dx + dy*dy + dz*dz;
	    cOpen = blk->cOpen.p[i];
	    cOpen2 = cOpen*cOpen;
	    d2Open = cOpen + k_Open;
	    d2Open = d2Open*d2Open;

	    dx = abs(xc-k_xCenter) - k_xMax;
	    dy = abs(yc-k_yCenter) - k_yMax;
	    dz = abs(zc-k_zCenter) - k_zMax;

	    dx = maskz_mov(dx>0,dx);
	    dy = maskz_mov(dy>0,dy);
	    dz = maskz_mov(dz>0,dz);
	    mink2 = dx*dx + dy*dy + dz*dz;
	    minbnd2 = 0.0f;

	    dx = k_xMinBnd - fvec(blk->xCenter.p[i]) - fvec(blk->xOffset.p[i]) - fvec(blk->xMax.p[i]);
	    minbnd2 += maskz_mov(dx>0,dx*dx);
	    dx = fvec(blk->xCenter.p[i]) + fvec(blk->xOffset.p[i]) - fvec(blk->xMax.p[i]) - k_xMaxBnd;
	    minbnd2 += maskz_mov(dx>0,dx*dx);

	    dx = k_yMinBnd - fvec(blk->yCenter.p[i]) - fvec(blk->yOffset.p[i]) - fvec(blk->yMax.p[i]);
	    minbnd2 += maskz_mov(dx>0,dx*dx);
	    dx = fvec(blk->yCenter.p[i]) + fvec(blk->yOffset.p[i]) - fvec(blk->yMax.p[i]) - k_yMaxBnd;
	    minbnd2 += maskz_mov(dx>0,dx*dx);

	    dx = k_zMinBnd - fvec(blk->zCenter.p[i]) - fvec(blk->zOffset.p[i]) - fvec(blk->zMax.p[i]);
	    minbnd2 += maskz_mov(dx>0,dx*dx);
	    dx = fvec(blk->zCenter.p[i]) + fvec(blk->zOffset.p[i]) - fvec(blk->zMax.p[i]) - k_zMaxBnd;
	    minbnd2 += maskz_mov(dx>0,dx*dx);

	    T0 = fvec(blk->m.p[i]) > fvec(0.0f);
	    T1 = (d2>d2Open) & (minbnd2>fourh2);
	    T2 = cvt_fvec(i32v(blk->iLower.p[i])) == 0.0;
	    T3 = (walk_min_multipole > cvt_fvec(i32v(blk->nc.p[i]))) | (mink2<=cOpen2);
	    T4 = minbnd2 > fourh2;
	    T6 = cOpen > k_Open;
	    T7 = k_notgrp;
 	    iOpenA = mask_mov(i32v(3),T2,i32v(1));
	    iOpenB = mask_mov(mask_mov(iOpenA,T4,i32v(4)),T3,iOpenA);
	    P1 = mask_mov(i32v(3),T2,i32v(2));
	    P2 = mask_mov(iOpenB,T7,i32v(0));
	    P3 = mask_mov(P2,T6,P1);
	    P4 = mask_mov(P3,T1,i32v(8));
	    iOpen = mask_mov(i32v(10),T0,P4);
	    blk->iOpen.p[i] = iOpen;
	    }
	}
    double dFlop = COST_FLOP_OPEN*(tile->lstTile.nBlocks*CL_PART_PER_BLK  + tile->lstTile.nInLast);
    pkd->dFlop += dFlop;
    pkd->dFlopSingleCPU += dFlop;
    }
