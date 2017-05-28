#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef USE_SIMD_PP
#define MPICH_SKIP_MPICXX
#include "simd.h"
#include "pkd.h"

static const float minSoftening = 1e-18f;

extern "C"
void pkdGravEvalPP(PINFOIN *pPart, int nBlocks, int nInLast, ILP_BLK *blk,  PINFOOUT *pOut ) {

    fvec t1, t2, t3, pd2;
    fvec pax, pay, paz, pfx, pfy, pfz, pdx, pdy, pdz;
    fvec piax, piay, piaz;
    fvec ppot /*,pmass,p4soft2*/;
    fvec padotai,pimaga,psmooth2,pirsum,pnorms;

    float ax,ay,az,fPot,dirsum,normsum;
    float tax,tay,taz;
    float dimaga;

    float fx = pPart->r[0];
    float fy = pPart->r[1];
    float fz = pPart->r[2];
    /*float fMass = pPart->.fMass;*/
    /*float fSoft = pPart->.fSoft;*/
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
    for( j = nInLast; j&SIMD_MASK; j++) {
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
    /*pmass   = fMass;*/
    /*p4soft2 = 4.0*fSoft*fSoft;*/
    psmooth2= fsmooth2;

    pax.zero();
    pay.zero();
    paz.zero();
    ppot.zero();
    pirsum.zero();
    pnorms.zero();

    nSoft = 0;
    nIntr = nBlocks * ILP_PART_PER_BLK + nInLast;
    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = ((nLeft ? ILP_PART_PER_BLK : nInLast) + SIMD_MASK) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    fvec pfourh2, td2, pir, pir2;
	    fmask vcmp;

	    pdx = blk->dx.p[j] + pfx;
	    pdy = blk->dy.p[j] + pfy;
	    pdz = blk->dz.p[j] + pfz;
	    pd2 = pdz*pdz + pdy*pdy + pdx*pdx;

	    pfourh2 = max(fvec(blk->fourh2.p[j]),minSoftening); /* There is always a self interaction */
	    vcmp = pd2 < pfourh2;
	    td2 = max(pd2,pfourh2);

	    td2 = max(minSoftening,td2);
	    pir = rsqrt(td2);
	    pir2 = pir * pir;
	    td2 = pir2 * pd2; /* for SOFTENED */
	    pir2 = pir2 * pir;

	    /* pir and pir2 are valid now for both softened and unsoftened particles */
	    /* Now we apply the fix to softened particles only */
	    if (!testz(vcmp)) {
		++nSoft;
		td2 = maskz_mov(vcmp,1.0f - td2);
		pir *= 1.0f + td2*(0.5f + td2*(3.0f/8.0f + td2*(45.0f/32.0f)));
		pir2 *= 1.0f + td2*(1.5f + td2*(135.0f/16.0f));
		}
	    pir2 *= -fvec(blk->m.p[j]);
	    t1 = pdx * pir2;
	    t2 = pdy * pir2;
	    t3 = pdz * pir2;
	    ppot -= fvec(blk->m.p[j])*pir;

#if 0
	    /* Time stepping criteria stuff */
	    padotai = SIMD_MADD(piaz,t3,SIMD_MADD(piay,t2,SIMD_MUL(piax,t1)));
	    vcmp = SIMD_AND(SIMD_CMP_GT(padotai,consts.zero.p),SIMD_CMP_GE(pd2,psmooth2));
	    padotai= SIMD_AND(padotai,vcmp);
	    padotai= SIMD_MUL(padotai,pimaga);
	    td2 = SIMD_MUL(padotai,padotai);
	    pirsum = SIMD_MADD(pir,td2,pirsum);
	    pnorms = SIMD_ADD(pnorms,td2);
#endif
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
