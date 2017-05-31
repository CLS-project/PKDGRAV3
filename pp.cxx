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
    float tax,tay,taz;
    float dimaga;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    float fsmooth2 = pPart->fSmooth2;
    float *a = pPart->a;
    int nSoft, nIntr;
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

    nSoft = 0;
    nIntr = nBlocks * ILP_PART_PER_BLK + nInLast;
    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = (nLeft ? ILP_PART_PER_BLK : nInLast+fvec::mask()) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    fvec Idx = blk->dx.p[j];
	    fvec Idy = blk->dy.p[j];
	    fvec Idz = blk->dz.p[j];
	    fvec Im = blk->m.p[j];
//Needed?   fvec fourh2 = max(fvec(blk->fourh2.p[j]),minSoftening); /* There is always a self interaction */
	    fvec fourh2 = blk->fourh2.p[j];
	    fvec pir,norm;
#if 1
	    EvalPP<fvec,fmask,true>(pfx,pfy,pfz,psmooth2,Idx,Idy,Idz,fourh2,Im,t1,t2,t3,pot,
	    	piax,piay,piaz,pimaga,pir,norm);
#else
	    static const float minSoftening = 1e-18f;
	    fvec dx = Idx + pfx;
	    fvec dy = Idy + pfy;
	    fvec dz = Idz + pfz;
	    fvec d2 = dx*dx + dy*dy + dz*dz;
	    fmask vcmp = d2 < fourh2;
	    fvec td2 = max(minSoftening,max(d2,fourh2));
	    pir = rsqrt(td2);
	    fvec pir2 = pir * pir;
	    td2 = pir2 * d2; /* for SOFTENED */
	    pir2 = pir2 * pir;

	    /* pir and pir2 are valid now for both softened and unsoftened particles */
	    /* Now we apply the fix to softened particles only */
	    if (!testz(vcmp)) {
		++nSoft;
		td2 = maskz_mov(vcmp,fvec(1.0f - td2));
		pir *= 1.0f + td2*(0.5f + td2*(3.0f/8.0f + td2*(45.0f/32.0f)));
		pir2 *= 1.0f + td2*(1.5f + td2*(135.0f/16.0f));
		}
	    pir2 *= -Im;
	    t1 = dx * pir2;
	    t2 = dy * pir2;
	    t3 = dz * pir2;
	    pot = -Im*pir;

	    /* Time stepping criteria stuff */
	    fvec padotai = piaz*t3 + piay*t2 + piax*t1;
	    vcmp = (padotai>0.0f) & (d2>=psmooth2);
	    maskz_mov(vcmp,padotai);
	    padotai= padotai * pimaga;
	    norm = padotai * padotai;
	    pir *= norm;
#endif
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
