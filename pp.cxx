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
#ifndef xxxUSE_SIMD_PP
#define MPICH_SKIP_MPICXX
#include "simd.h"
#include "pkd.h"
#include "pp.h"

extern "C"
void pkdGravEvalPP(PINFOIN *pPart, int nBlocks, int nInLast, ILP_BLK *blk,  PINFOOUT *pOut ) {
    fvec t1, t2, t3, pot;
    fvec pax, pay, paz, pfx, pfy, pfz;
    fvec piax, piay, piaz;
    fvec ppot /*,pmass,p4soft2*/;
    fvec pimaga,psmooth2,pirsum,pnorms;

    float ax,ay,az,fPot,dirsum,normsum;
    float dimaga;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    float fsmooth2 = pPart->fSmooth2;
    float *a = pPart->a;
    int nLeft, j;

    dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0) {
	dimaga = 1.0f/sqrtf(dimaga);
	}
    pimaga = dimaga;


    /*
    ** This is a little trick to speed up the calculation. By setting
    ** unused entries in the list to have a zero mass, the resulting
    ** forces are zero. Setting the distance to a large value avoids
    ** softening the non-existent forces which is slightly faster.
    */
    for( j = nInLast; j&fvec::mask(); j++) {
	blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
	blk[nBlocks].m.f[j] = 0.0f;
	blk[nBlocks].fourh2.f[j] = 1e-18f;
	}

    piax    = a[0];
    piay    = a[1];
    piaz    = a[2];
    pfx     = fx;
    pfy     = fy;
    pfz     = fz;
    psmooth2= fsmooth2;

    pax = pay = paz = ppot = pirsum = pnorms = 0.0;

    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = (nLeft ? ILP_PART_PER_BLK : nInLast+fvec::mask()) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    fvec Idx = blk->dx.p[j];
	    fvec Idy = blk->dy.p[j];
	    fvec Idz = blk->dz.p[j];
	    fvec Im = blk->m.p[j];
	    fvec fourh2 = blk->fourh2.p[j];
	    fvec pir,norm;
	    EvalPP<fvec,fmask,true>(pfx,pfy,pfz,psmooth2,Idx,Idy,Idz,fourh2,Im,t1,t2,t3,pot,
	    	piax,piay,piaz,pimaga,pir,norm);
	    pirsum += pir;
	    pnorms += norm;
	    ppot += pot;
	    pax += t1;
	    pay += t2;
	    paz += t3;
	    }
	}
    ax = hadd(pax);
    ay = hadd(pay);
    az = hadd(paz);
    fPot = hadd(ppot);
    dirsum = hadd(pirsum);
    normsum = hadd(pnorms);

    pOut->a[0] += ax;
    pOut->a[1] += ay;
    pOut->a[2] += az;
    pOut->fPot += fPot;
    pOut->dirsum += dirsum;
    pOut->normsum += normsum;
}

extern "C"
void pkdDensityEval(PINFOIN *pPart, int nBlocks, int nInLast, ILP_BLK *blk,  PINFOOUT *pOut, SPHOptions *SPHoptions ) {
    fvec t1, t2, t3, t4, t5;
    fvec parho, padrhodfball, panden, padndendfball, pfx, pfy, pfz, pfBall;
    fvec pnSmooth;

    float arho, adrhodfball, anden, adndendfball;
    float anSmooth;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    float fBall = pPart->fBall;
    int nLeft, j;

    /*
    ** This is a little trick to speed up the calculation. By setting
    ** unused entries in the list to have a zero mass, the resulting
    ** forces are zero. Setting the distance to a large value avoids
    ** softening the non-existent forces which is slightly faster.
    */
    for( j = nInLast; j&fvec::mask(); j++) {
	blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
	blk[nBlocks].m.f[j] = 0.0f;
	}

    pfx     = fx;
    pfy     = fy;
    pfz     = fz;
    pfBall  = fBall;

    parho = padrhodfball = panden = padndendfball = 0.0f;
    pnSmooth = 0.0f;

    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = (nLeft ? ILP_PART_PER_BLK : nInLast+fvec::mask()) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    fvec Idx = blk->dx.p[j];
	    fvec Idy = blk->dy.p[j];
	    fvec Idz = blk->dz.p[j];
	    fvec Im = blk->m.p[j];
	    EvalDensity<fvec,fmask>(pfx,pfy,pfz,Idx,Idy,Idz,Im,pfBall,t1,t2,t3,t4,t5,SPHoptions);
	    parho += t1;
	    padrhodfball += t2;
        panden += t3;
        padndendfball += t4;
        pnSmooth += t5;
	    }
	}
    arho = hadd(parho);
    adrhodfball = hadd(padrhodfball);
    anden = hadd(panden);
    adndendfball = hadd(padndendfball);
    anSmooth = hadd(pnSmooth);

    pOut->rho += arho;
    pOut->drhodfball += adrhodfball;
    pOut->nden += anden;
    pOut->dndendfball += adndendfball;
    pOut->nSmooth += anSmooth;
}

extern "C"
void pkdSPHForcesEval(PINFOIN *pPart, int nBlocks, int nInLast, ILP_BLK *blk,  PINFOOUT *pOut, SPHOptions *SPHoptions ) {
    fvec t1, t2, t3, t4, t5, t6;
    fvec PfBall, POmega, Pdx, Pdy, Pdz, Pvx, Pvy, Pvz, Prho, PP, Pc;
    fvec IfBall, IOmega, Idx, Idy, Idz, Ivx, Ivy, Ivz, Irho, IP, Ic, Im;
    i32v Pspecies, Ispecies;
    fvec puDot, pax, pay, paz, pdivv, pdtEst;
    float auDot, aax, aay, aaz, adivv, adtEst;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    float fBall = pPart->fBall;
    float fOmega = pPart->Omega;
    float fvx = pPart->v[0];
    float fvy = pPart->v[1];
    float fvz = pPart->v[2];
    float frho = pPart->rho;
    float fP = pPart->P;
    float fc = pPart->c;
    int32_t ispecies = pPart->species;
    int nLeft, j;

    /*
    ** This is a little trick to speed up the calculation. By setting
    ** unused entries in the list to have a zero mass, the resulting
    ** forces are zero. Setting the distance to a large value avoids
    ** softening the non-existent forces which is slightly faster.
    */
    for( j = nInLast; j&fvec::mask(); j++) {
	blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
    blk[nBlocks].fBall.f[j] = blk[nBlocks].Omega.f[j] = blk[nBlocks].rho.f[j] = 1e18f;
	blk[nBlocks].m.f[j] = 0.0f;
    blk[nBlocks].P.f[j] = 0.0f;
	}

    Pdx     = fx;
    Pdy     = fy;
    Pdz     = fz;
    PfBall  = fBall;
    POmega  = fOmega;
    Pvx     = fvx;
    Pvy     = fvy;
    Pvz     = fvz;
    Prho    = frho;
    PP      = fP;
    Pc      = fc;
    Pspecies= ispecies;

    puDot = pax = pay = paz = pdivv = 0.0f;
    pdtEst = HUGE_VALF;

    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = (nLeft ? ILP_PART_PER_BLK : nInLast+fvec::mask()) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    fvec Idx = blk->dx.p[j];
	    fvec Idy = blk->dy.p[j];
	    fvec Idz = blk->dz.p[j];
	    fvec Im = blk->m.p[j];
        fvec IfBall = blk->fBall.p[j];
        fvec IOmega = blk->Omega.p[j];
        fvec Ivx = blk->vx.p[j];
        fvec Ivy = blk->vy.p[j];
        fvec Ivz = blk->vz.p[j];
        fvec Irho = blk->rho.p[j];
        fvec IP = blk->P.p[j];
        fvec Ic = blk->c.p[j];
        i32v Ispecies = blk->species.p[j];

	    EvalSPHForces<fvec,fmask,i32v>(Pdx,Pdy,Pdz,PfBall,POmega,Pvx,Pvy,Pvz,Prho,PP,Pc,Pspecies,
            Idx,Idy,Idz,Im,IfBall,IOmega,Ivx,Ivy,Ivz,Irho,IP,Ic,Ispecies,
            t1,t2,t3,t4,t5,t6,
            SPHoptions);
            puDot += t1;
            pax += t2;
            pay += t3;
            paz += t4;
            pdivv += t5;
            pdtEst = min(pdtEst,t6);
	    }
	}
    auDot = hadd(puDot);
    aax = hadd(pax);
    aay = hadd(pay);
    aaz = hadd(paz);
    adivv = hadd(pdivv);
    adtEst = HUGE_VALF; 
    // This should be a horizontal minimum for an fvec, resulting in a float containing the smallest float in the fvec
    for (int k = 0; k < pdtEst.width(); k++) {
        adtEst = fmin(adtEst,pdtEst[k]);
    }
    assert(adtEst > 0);

    pOut->uDot += auDot;
    pOut->a[0] += aax;
    pOut->a[1] += aay;
    pOut->a[2] += aaz;
    pOut->divv += adivv;
    pOut->dtEst = fmin(pOut->dtEst,adtEst);
}
#endif/*USE_SIMD_PP*/
