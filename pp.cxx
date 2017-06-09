#ifdef HAVE_CONFIG_H
#include "config.h"
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
#endif/*USE_SIMD_PP*/
