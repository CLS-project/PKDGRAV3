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
#ifndef xxUSE_SIMD_PC
#define MPICH_SKIP_MPICXX
#include "core/simd.h"
#include "pkd.h"
#include "pc.h"

void pkdGravEvalPC(PINFOIN *pPart, int nBlocks, int nInLast, ILC_BLK *blk,  PINFOOUT *pOut ) {
    fvec x,y,z;
    fvec adotai;
    fvec tx,ty,tz;
    fvec xx,xy,xz,yy,yz,zz;
    fvec xxx,xxz,yyy,yyz,xxy,xyy,xyz;

    int j, nLeft;

    fvec fx = pPart->r[0];
    fvec fy = pPart->r[1];
    fvec fz = pPart->r[2];
    fvec pSmooth2 = pPart->fSmooth2;
    fvec Pax = pPart->a[0];
    fvec Pay = pPart->a[1];
    fvec Paz = pPart->a[2];

    fvec ax,ay,az,fPot,dirsum,normsum;
    float *a =pPart->a;
    float a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    fvec imaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;

    /* Pad the last value if necessary */
    for ( j = nInLast; j&fvec::mask(); j++) {
        blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
        blk[nBlocks].m.f[j] = 0.0f;
        blk[nBlocks].u.f[j] = 0.0f;
    }

    ax = 0.0;
    ay = 0.0;
    az = 0.0;
    fPot= 0.0;
    dirsum = 0.0;
    normsum = 0.0;
    //int nIntr = nBlocks * ILP_PART_PER_BLK + nInLast;
    for ( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
        int n = (nLeft ? ILP_PART_PER_BLK : nInLast + fvec::mask()) >> SIMD_BITS;
        for (j=0; j<n; ++j) {
            fvec Idx = blk->dx.p[j];
            fvec Idy = blk->dy.p[j];
            fvec Idz = blk->dz.p[j];
            fvec Ixxxx = blk->xxxx.p[j];
            fvec Ixxxy = blk->xxxy.p[j];
            fvec Ixxxz = blk->xxxz.p[j];
            fvec Ixxyz = blk->xxyz.p[j];
            fvec Ixxyy = blk->xxyy.p[j];
            fvec Iyyyz = blk->yyyz.p[j];
            fvec Ixyyz = blk->xyyz.p[j];
            fvec Ixyyy = blk->xyyy.p[j];
            fvec Iyyyy = blk->yyyy.p[j];
            fvec Ixxx = blk->xxx.p[j];
            fvec Ixyy = blk->xyy.p[j];
            fvec Ixxy = blk->xxy.p[j];
            fvec Iyyy = blk->yyy.p[j];
            fvec Ixxz = blk->xxz.p[j];
            fvec Iyyz = blk->yyz.p[j];
            fvec Ixyz = blk->xyz.p[j];
            fvec Ixx = blk->xx.p[j];
            fvec Ixy = blk->xy.p[j];
            fvec Ixz = blk->xz.p[j];
            fvec Iyy = blk->yy.p[j];
            fvec Iyz = blk->yz.p[j];
#ifdef USE_DIAPOLE
            fvec Ix = blk->x.p[j];
            fvec Iy = blk->y.p[j];
            fvec Iz = blk->z.p[j];
#endif
            fvec Im = blk->m.p[j];
            fvec Iu = blk->u.p[j];

            auto result = EvalPC<fvec,fmask,true>(
                              fx, fy, fz, pSmooth2,
                              Idx, Idy, Idz, Im, Iu,
                              Ixxxx, Ixxxy, Ixxxz, Ixxyz, Ixxyy, Iyyyz, Ixyyz, Ixyyy, Iyyyy,
                              Ixxx, Ixyy, Ixxy, Iyyy, Ixxz, Iyyz, Ixyz, Ixx, Ixy, Ixz, Iyy, Iyz,
#ifdef USE_DIAPOLE
                              Ix, Iy, Iz,
#endif
                              Pax, Pay, Paz,imaga);

            dirsum += result.ir;
            normsum += result.norm;
            fPot += result.pot;
            ax += result.ax;
            ay += result.ay;
            az += result.az;
        }
    }

    pOut->a[0] += hadd(ax);
    pOut->a[1] += hadd(ay);
    pOut->a[2] += hadd(az);
    pOut->fPot += hadd(fPot);
    pOut->dirsum += hadd(dirsum);
    pOut->normsum += hadd(normsum);
}
#endif/*USE_SIMD_PC*/
