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
#include "core/simd.h"
#include "pkd.h"
#include "opening.h"
#include "../SPH/SPHOptions.h"

/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
**
** This version will also open buckets ("new" criteria)
*/

void iOpenOutcomeBlock(PKD pkd,treeStore::NodePointer k,int n,clBlock &blk,float dThetaMin,SPHOptions *SPHoptions) {
    const float walk_min_multipole = 3;
    fmask T0,T1,T2,T3,T4,T6,T7;
    i32v P1,P2,P3,P4;
    fvec xc,yc,zc,dx,dy,dz,d2,diCrit,cOpen,cOpen2,d2Open,mink2,minbnd2,fourh2;
#if SPHBALLOFBALLS
    fvec distk2,distc2;
#endif
#if SPHBOXOFBALLS
    fvec box1xMin, box1xMax, box1yMin, box1yMax, box1zMin, box1zMax, box2xMin, box2xMax, box2yMin, box2yMax, box2zMin, box2zMax;
#endif
    int i;
    i32v iOpen,iOpenA,iOpenB;
    fvec k_xCenter, k_yCenter, k_zCenter, k_xMax, k_yMax, k_zMax;
    fvec k_xMinBnd, k_yMinBnd, k_zMinBnd, k_xMaxBnd, k_yMaxBnd, k_zMaxBnd;
    fvec k_x, k_y, k_z, k_bMax, k_Open;
    fmask k_notgrp;
    auto k_r = k->position();

    assert ( k->mass() > 0.0f );

    diCrit = 1.0f/dThetaMin;

    auto kbnd = k->bound();
    k_xMinBnd = kbnd.lower(0);
    k_yMinBnd = kbnd.lower(1);
    k_zMinBnd = kbnd.lower(2);
    k_xMaxBnd = kbnd.upper(0);
    k_yMaxBnd = kbnd.upper(1);
    k_zMaxBnd = kbnd.upper(2);
    k_xCenter = kbnd.center(0);
    k_yCenter = kbnd.center(1);
    k_zCenter = kbnd.center(2);
    k_xMax = kbnd.apothem(0);
    k_yMax = kbnd.apothem(1);
    k_zMax = kbnd.apothem(2);
    k_x = k_r[0];
    k_y = k_r[1];
    k_z = k_r[2];
    k_bMax = k->bMax();
    k_notgrp = cvt_fvec(i32v(k->is_group())) == 0.0;
    k_Open = 1.5f*k_bMax*diCrit;
    const auto &SPHbob = k->have_BOB() ? k->BOB() : SPHBOB();
#if SPHBALLOFBALLS
    fvec k_fBoBr2 = SPHbob.fBoBr * SPHbob.fBoBr;
    fvec blk_fBoBr2 = 0.0f;
    distk2 = HUGE_VALF;
    distc2 = HUGE_VALF;
#endif

    fmask intersect1, intersect2;
    intersect1 = 0;
    intersect2 = 0;

    n = (n + fvec::mask()) / fvec::width(); // Now number of blocks
    for (i=0; i<n; ++i) {
        fourh2 = blk.fourh2.v[i];

        if (pkd->tree.present(KDN_FIELD::oNodeBOB)) {
            if (SPHoptions->doDensity || SPHoptions->doDensityCorrection || SPHoptions->doSPHForces || SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) {
#if SPHBALLOFBALLS
                distk2 = 0.0f;
                dx = SPHbob.fBoBCenter[0] - blk.xCenter.v[i] - blk.xOffset.v[i] - blk.xMax.v[i];
                distk2 += maskz_mov(dx>0,dx*dx);
                dx = blk.xCenter.v[i] + blk.xOffset.v[i] - blk.xMax.v[i] - SPHbob.fBoBCenter[0];
                distk2 += maskz_mov(dx>0,dx*dx);

                dx = SPHbob.fBoBCenter[1] - blk.yCenter.v[i] - blk.yOffset.v[i] - blk.yMax.v[i];
                distk2 += maskz_mov(dx>0,dx*dx);
                dx = blk.yCenter.v[i] + blk.yOffset.v[i] - blk.yMax.v[i] - SPHbob.fBoBCenter[1];
                distk2 += maskz_mov(dx>0,dx*dx);

                dx = SPHbob.fBoBCenter[2] - blk.zCenter.v[i] - blk.zOffset.v[i] - blk.zMax.v[i];
                distk2 += maskz_mov(dx>0,dx*dx);
                dx = blk.zCenter.v[i] + blk.zOffset.v[i] - blk.zMax.v[i] - SPHbob.fBoBCenter[2];
                distk2 += maskz_mov(dx>0,dx*dx);
                intersect1 = distk2 < k_fBoBr2;
#endif
#if SPHBOXOFBALLS
                box1xMin = SPHbob.fBoBMin[0];
                box1xMax = SPHbob.fBoBMax[0];
                box1yMin = SPHbob.fBoBMin[1];
                box1yMax = SPHbob.fBoBMax[1];
                box1zMin = SPHbob.fBoBMin[2];
                box1zMax = SPHbob.fBoBMax[2];
                box2xMin = blk.xCenter.v[i] + blk.xOffset.v[i] - blk.xMax.v[i];
                box2xMax = blk.xCenter.v[i] + blk.xOffset.v[i] + blk.xMax.v[i];
                box2yMin = blk.yCenter.v[i] + blk.yOffset.v[i] - blk.yMax.v[i];
                box2yMax = blk.yCenter.v[i] + blk.yOffset.v[i] + blk.yMax.v[i];
                box2zMin = blk.zCenter.v[i] + blk.zOffset.v[i] - blk.zMax.v[i];
                box2zMax = blk.zCenter.v[i] + blk.zOffset.v[i] + blk.zMax.v[i];
                intersect1 = (box1xMin < box2xMax) & (box2xMin < box1xMax) & (box1yMin < box2yMax) & (box2yMin < box1yMax) & (box1zMin < box2zMax) & (box2zMin < box1zMax);
#endif
            }
            if (SPHoptions->doSPHForces || SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) {
#if SPHBALLOFBALLS
                blk_fBoBr2 = blk.fBoBr.v[i] * blk.fBoBr.v[i];
                distc2 = 0.0f;
                dx = k_xMinBnd - blk.fBoBxCenter.v[i] - blk.xOffset.v[i];
                distc2 += maskz_mov(dx>0,dx*dx);
                dx = blk.fBoBxCenter.v[i] + blk.xOffset.v[i] - k_xMaxBnd;
                distc2 += maskz_mov(dx>0,dx*dx);

                dx = k_yMinBnd - blk.fBoByCenter.v[i] - blk.yOffset.v[i];
                distc2 += maskz_mov(dx>0,dx*dx);
                dx = blk.fBoByCenter.v[i] + blk.yOffset.v[i] - k_yMaxBnd;
                distc2 += maskz_mov(dx>0,dx*dx);

                dx = k_zMinBnd - blk.fBoBzCenter.v[i] - blk.zOffset.v[i];
                distc2 += maskz_mov(dx>0,dx*dx);
                dx = blk.fBoBzCenter.v[i] + blk.zOffset.v[i] - k_zMaxBnd;
                distc2 += maskz_mov(dx>0,dx*dx);
                intersect2 = distc2 < blk_fBoBr2;
#endif
#if SPHBOXOFBALLS
                box1xMin = k_xMinBnd;
                box1xMax = k_xMaxBnd;
                box1yMin = k_yMinBnd;
                box1yMax = k_yMaxBnd;
                box1zMin = k_zMinBnd;
                box1zMax = k_zMaxBnd;
                box2xMin = blk.fBoBxMin.v[i] + blk.xOffset.v[i];
                box2xMax = blk.fBoBxMax.v[i] + blk.xOffset.v[i];
                box2yMin = blk.fBoByMin.v[i] + blk.yOffset.v[i];
                box2yMax = blk.fBoByMax.v[i] + blk.yOffset.v[i];
                box2zMin = blk.fBoBzMin.v[i] + blk.zOffset.v[i];
                box2zMax = blk.fBoBzMax.v[i] + blk.zOffset.v[i];
                intersect2 = (box1xMin < box2xMax) & (box2xMin < box1xMax) & (box1yMin < box2yMax) & (box2yMin < box1yMax) & (box1zMin < box2zMax) & (box2zMin < box1zMax);
#endif
            }
        }
        xc = blk.x.v[i] + blk.xOffset.v[i];
        yc = blk.y.v[i] + blk.yOffset.v[i];
        zc = blk.z.v[i] + blk.zOffset.v[i];
        dx = k_x - xc;
        dy = k_y - yc;
        dz = k_z - zc;
        d2 = dx*dx + dy*dy + dz*dz;
        cOpen = blk.cOpen.v[i];
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

        dx = k_xMinBnd - blk.xCenter.v[i] - blk.xOffset.v[i] - blk.xMax.v[i];
        minbnd2 += maskz_mov(dx>0,dx*dx);
        dx = blk.xCenter.v[i] + blk.xOffset.v[i] - blk.xMax.v[i] - k_xMaxBnd;
        minbnd2 += maskz_mov(dx>0,dx*dx);

        dx = k_yMinBnd - blk.yCenter.v[i] - blk.yOffset.v[i] - blk.yMax.v[i];
        minbnd2 += maskz_mov(dx>0,dx*dx);
        dx = blk.yCenter.v[i] + blk.yOffset.v[i] - blk.yMax.v[i] - k_yMaxBnd;
        minbnd2 += maskz_mov(dx>0,dx*dx);

        dx = k_zMinBnd - blk.zCenter.v[i] - blk.zOffset.v[i] - blk.zMax.v[i];
        minbnd2 += maskz_mov(dx>0,dx*dx);
        dx = blk.zCenter.v[i] + blk.zOffset.v[i] - blk.zMax.v[i] - k_zMaxBnd;
        minbnd2 += maskz_mov(dx>0,dx*dx);

        T0 = blk.m.v[i] > 0.0f;
        T1 = (d2>d2Open) & (minbnd2>fourh2) & ~intersect1 & ~intersect2;
        T2 = cvt_fvec(blk.iLower.v[i]) == 0.0;
        T3 = (walk_min_multipole > cvt_fvec(blk.nc.v[i])) | (mink2<=cOpen2);
        T4 = (minbnd2 > fourh2) & ~intersect1 & ~intersect2;
        T6 = cOpen > k_Open;
        T7 = k_notgrp;
        iOpenA = mask_mov(i32v(3),T2,i32v(1));
        iOpenB = mask_mov(mask_mov(iOpenA,T4,i32v(4)),T3,iOpenA);
        P1 = mask_mov(i32v(3),T2,i32v(2));
        P2 = mask_mov(iOpenB,T7,i32v(0));
        P3 = mask_mov(P2,T6,P1);
        P4 = mask_mov(P3,T1,i32v(8));
        iOpen = mask_mov(i32v(10),T0,P4);
        blk.iOpen.v[i] = iOpen;
    }

    double dFlop = COST_FLOP_OPEN*n;
    pkd->dFlop += dFlop;
    pkd->dFlopSingleCPU += dFlop;
}

void iOpenOutcomeSIMD(PKD pkd,treeStore::NodePointer k,clTile &tile,float dThetaMin,SPHOptions *SPHoptions) {
    // Do the full blocks, followed by the last (probably not full) block
    auto nBlocks = tile.count() / tile.width;
    for (auto iBlock=0; iBlock<nBlocks; ++iBlock) {
        iOpenOutcomeBlock(pkd,k,tile.width,tile[iBlock],dThetaMin,SPHoptions);
    }
    iOpenOutcomeBlock(pkd,k,tile.count() - nBlocks*tile.width,tile[nBlocks],dThetaMin,SPHoptions);
}
