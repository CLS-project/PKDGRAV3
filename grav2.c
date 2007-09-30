#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include "pkd.h"
#include "moments.h"
#include "meval.h"
#include "qeval.h"
#include "ewald.h"
#include "grav.h"

#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    } 

void HEAPrholocal(int n, int k, RHOLOCAL ra[]) {
    int l,j,ir,i;
    RHOLOCAL rra;

    l=(n>>1)+1;
    ir=n;

    while (1) {
        if (l > 1) {
            rra = ra[--l-1];
            }
        else {
            rra = ra[ir-1];
            ra[ir-1] = ra[0];
            if (--ir == n-k) {
                ra[0] = rra;
                return;
                }
            }
        i = l;
        j = (l<<1);
        while (j <= ir) {
            if (j < ir && ra[j-1].d2 > ra[j].d2) ++j; /* k smallest elements at the end */
            if (rra.d2 > ra[j-1].d2) {
                ra[i-1] = ra[j-1];
                j += (i=j);
                }
            else j = ir+1;
            }
        ra[i-1] = rra;
        }
    }

/*
** This version of grav.c does all the operations inline, including 
** v_sqrt's and such.
** Returns nActive.
*/
int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *pLoc,ILP *ilp,int nPart,ILC *ilc,int nCell,double dirLsum,double normLsum,
		    int bEwald,int bEwaldKicking,double *pdFlop,double *pdEwFlop) {
    PARTICLE *p = pkd->pStore;
    KDN *pkdn = pBucket;
    const double onethird = 1.0/3.0;
    momFloat ax,ay,az,fPot;
    double x,y,z,d2,dir,dir2;
    momFloat adotai,maga,dimaga,dirsum,normsum;
    momFloat tax,tay,taz,tmon;
    double rholoc,dirDTS,d2DTS,dsmooth2;
#ifndef USE_SIMD_MOMR
    double g2,g3,g4;
    double xx,xy,xz,yy,yz,zz;
    double xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    double tx,ty,tz;
#else
    int nCellILC;
#endif
    double fourh2;
#ifdef SOFTSQUARE
    double ptwoh2;
#endif
    int i,j,k,nN,nSP,nSoft,nActive; 

#ifdef USE_SIMD_MOMR
    nCellILC = nCell;
    momPadSIMDMomr( &nCellILC, ilc );
#endif
    /*
    ** dynamical time-stepping stuff
    */
    RHOLOCAL *rholocal;
    nN = nPart+pkdn->pUpper-pkdn->pLower; /* total number of neighbouring particles (without particle itself) */
    rholocal = malloc(nN*sizeof(RHOLOCAL));
    assert(rholocal != NULL);
    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	if (!pkdIsActive(pkd,&p[i])) continue;
	++nActive;
	p[i].dtGrav = 0.0;
	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
        rholoc = 0;
        dsmooth2 = 0;
	tx = p[i].a[0] + p[i].ae[0];
	ty = p[i].a[1] + p[i].ae[1];
	tz = p[i].a[2] + p[i].ae[2];
	dimaga = tx*tx + ty*ty + tz*tz;
	if (dimaga > 0) {
	    dimaga = 1.0/sqrt(dimaga);
	    }
	dirsum = dirLsum;
	normsum = normLsum;
	if (pLoc) {
	    /*
	    ** Evaluate local expansion.
	    */
	    x = p[i].r[0] - pkdn->r[0];
	    y = p[i].r[1] - pkdn->r[1];
	    z = p[i].r[2] - pkdn->r[2];
	    momEvalLocr(pLoc,x,y,z,&fPot,&ax,&ay,&az);
	}

#ifdef USE_SIMD_MOMR
	momEvalSIMDMomr( nCellILC, ilc, p[i].r, p[i].a,
			 &ax, &ay, &az, &fPot, &dirsum, &normsum );
#else
	for (j=0;j<nCell;++j) {
	    x = p[i].r[0] - ilc[j].x;
	    y = p[i].r[1] - ilc[j].y;
	    z = p[i].r[2] - ilc[j].z;
	    d2 = x*x + y*y + z*z;
	    SQRT1(d2,dir);
            dirDTS = dir;
	    dir2 = dir*dir;
	    g2 = 3*dir*dir2*dir2;
	    g3 = 5*g2*dir2;
	    g4 = 7*g3*dir2;
	    /*
	    ** Calculate the funky distance terms.
	    */
	    xx = 0.5*x*x;
	    xy = x*y;
	    xz = x*z;
	    yy = 0.5*y*y;
	    yz = y*z;
	    zz = 0.5*z*z;
	    xxx = x*(onethird*xx - zz);
	    xxz = z*(xx - onethird*zz);
	    yyy = y*(onethird*yy - zz);
	    yyz = z*(yy - onethird*zz);
	    xx -= zz;
	    yy -= zz;
	    xxy = y*xx;
	    xyy = x*yy;
	    xyz = xy*z;
	    /*
	    ** Now calculate the interaction up to Hexadecapole order.
	    */
	    tx = g4*(ilc[j].mom.xxxx*xxx + ilc[j].mom.xyyy*yyy + ilc[j].mom.xxxy*xxy + ilc[j].mom.xxxz*xxz + ilc[j].mom.xxyy*xyy + ilc[j].mom.xxyz*xyz + ilc[j].mom.xyyz*yyz);
	    ty = g4*(ilc[j].mom.xyyy*xyy + ilc[j].mom.xxxy*xxx + ilc[j].mom.yyyy*yyy + ilc[j].mom.yyyz*yyz + ilc[j].mom.xxyy*xxy + ilc[j].mom.xxyz*xxz + ilc[j].mom.xyyz*xyz);
	    tz = g4*(-ilc[j].mom.xxxx*xxz - (ilc[j].mom.xyyy + ilc[j].mom.xxxy)*xyz - ilc[j].mom.yyyy*yyz + ilc[j].mom.xxxz*xxx + ilc[j].mom.yyyz*yyy - ilc[j].mom.xxyy*(xxz + yyz) + ilc[j].mom.xxyz*xxy + ilc[j].mom.xyyz*xyy);
	    g4 = 0.25*(tx*x + ty*y + tz*z);
	    xxx = g3*(ilc[j].mom.xxx*xx + ilc[j].mom.xyy*yy + ilc[j].mom.xxy*xy + ilc[j].mom.xxz*xz + ilc[j].mom.xyz*yz);
	    xxy = g3*(ilc[j].mom.xyy*xy + ilc[j].mom.xxy*xx + ilc[j].mom.yyy*yy + ilc[j].mom.yyz*yz + ilc[j].mom.xyz*xz);
	    xxz = g3*(-(ilc[j].mom.xxx + ilc[j].mom.xyy)*xz - (ilc[j].mom.xxy + ilc[j].mom.yyy)*yz + ilc[j].mom.xxz*xx + ilc[j].mom.yyz*yy + ilc[j].mom.xyz*xy);
	    g3 = onethird*(xxx*x + xxy*y + xxz*z);
	    xx = g2*(ilc[j].mom.xx*x + ilc[j].mom.xy*y + ilc[j].mom.xz*z);
	    xy = g2*(ilc[j].mom.yy*y + ilc[j].mom.xy*x + ilc[j].mom.yz*z);
	    xz = g2*(-(ilc[j].mom.xx + ilc[j].mom.yy)*z + ilc[j].mom.xz*x + ilc[j].mom.yz*y);
	    g2 = 0.5*(xx*x + xy*y + xz*z);
	    tmon = ilc[j].mom.m*dir;
	    dir2 *= tmon + 5*g2 + 7*g3 + 9*g4;
	    fPot -= tmon + g2 + g3 + g4;
	    tax = xx + xxx + tx - x*dir2;
	    tay = xy + xxy + ty - y*dir2;
	    taz = xz + xxz + tz - z*dir2;
	    adotai = (p[i].a[0]+p[i].ae[0])*tax + (p[i].a[1]+p[i].ae[1])*tay + (p[i].a[2]+p[i].ae[2])*taz; 
            if (adotai > 0) {
		adotai *= dimaga;
                dirsum += dirDTS*adotai*adotai;
                normsum += adotai*adotai;
                }
	    ax += tax;
	    ay += tay;
	    az += taz;
	    } /* end of cell list gravity loop */
#endif		
	mdlCacheCheck(pkd->mdl);
	/*
	** Calculate local density and kernel smoothing length for dynamical time-stepping
	*/
	if(pkd->param.bGravStep) {
	    /*
	    ** Add particles to array rholocal first
	    */
	    for (j=0;j<nPart;++j) {
		x = p[i].r[0] - ilp[j].x;
		y = p[i].r[1] - ilp[j].y;
		z = p[i].r[2] - ilp[j].z;
		d2 = x*x + y*y + z*z;
		rholocal[j].m = ilp[j].m;	
		rholocal[j].d2 = d2;
		}
	    /*
	    ** Add bucket particles to the array rholocal as well!
	    ** Not including yourself!!
	    */
	    k = nPart;
	    for (j=pkdn->pLower;j<=pkdn->pUpper;++j) {
		if(p[i].iOrder == p[j].iOrder) continue;
		x = p[i].r[0] - p[j].r[0];
		y = p[i].r[1] - p[j].r[1];
		z = p[i].r[2] - p[j].r[2];
		d2 = x*x + y*y + z*z;
		rholocal[k].m = p[j].fMass;
		rholocal[k].d2 = d2;
		k += 1;
		}
	    /*
	    ** Calculate local density only in the case of more than 1 neighbouring particle!
	    */
	    if (nN > 1) {
		nSP = (nN < pkd->param.nPartRhoLoc)?nN:pkd->param.nPartRhoLoc;
		HEAPrholocal(nN,nSP,rholocal);
		dsmooth2 = rholocal[nN-nSP].d2;
		SQRT1(dsmooth2,dir);
		dir *= dir;
		for (j=(nN - nSP);j<nN;++j) {
		    d2 = rholocal[j].d2*dir;
		    d2 = (1-d2);
		    rholoc += d2*rholocal[j].m;
		    }
		rholoc = 1.875*M_1_PI*rholoc*dir*sqrt(dir); /* F1 Kernel (15/8) */
		}
	    assert(rholoc >= 0);
	    }
#ifdef USE_SIMD_PP
	PPInteractSIMD( nPart,ilp,p[i].r,p[i].a,p[i].fMass,p[i].fSoft,
			&ax, &ay, &az, &fPot, &dirsum, &normsum );
#else
	for (j=0;j<nPart;++j) {
	    x = p[i].r[0] - ilp[j].x;
	    y = p[i].r[1] - ilp[j].y;
	    z = p[i].r[2] - ilp[j].z;
	    d2 = x*x + y*y + z*z;
	    d2DTS = d2;
	    fourh2 = softmassweight(p[i].fMass,4*p[i].fSoft*p[i].fSoft,ilp[j].m,ilp[j].fourh2);
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
	    }
	    else {
		/*
		** This uses the Dehnen K1 kernel function now, it's fast!
		*/
		SQRT1(fourh2,dir);
		dir2 = dir*dir;
		d2 *= dir2;
		dir2 *= dir;
		d2 = 1 - d2;
		dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
		dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
		++nSoft;
		}
	    dir2 *= ilp[j].m;
	    tax = -x*dir2;
	    tay = -y*dir2;
	    taz = -z*dir2;
	    adotai = (p[i].a[0]+p[i].ae[0])*tax + (p[i].a[1]+p[i].ae[1])*tay + (p[i].a[2]+p[i].ae[2])*taz;
	    if (adotai > 0 && d2DTS >= dsmooth2) {
		adotai *= dimaga;
		dirsum += dir*adotai*adotai;
		normsum += adotai*adotai;
		}
	    fPot -= ilp[j].m*dir;
	    ax += tax;
	    ay += tay;
	    az += taz;
	    } /* end of particle list gravity loop */
#endif
        /*
        ** Finally set new acceleration and potential.
        ** Note that after this point we cannot use the new timestepping criterion since we
        ** overwrite the acceleration.
        */
	p[i].fPot = fPot;
	p[i].a[0] = ax;
	p[i].a[1] = ay;
	p[i].a[2] = az;	
	/*
	** Now finally calculate the Ewald correction for this particle, if it is
	** required.
	*/
	if (bEwald) {
	    *pdEwFlop += pkdParticleEwald(pkd,&p[i],bEwaldKicking);
	}
	/*
	** Set value for time-step, note that we have the current ewald acceleration
	** in this now as well!
	*/
	if (pkd->param.bGravStep) {
	    /*
	    ** Use new acceleration here!
	    */
	    tx = p[i].a[0] + p[i].ae[0];
	    ty = p[i].a[1] + p[i].ae[1];
	    tz = p[i].a[2] + p[i].ae[2];
	    maga = sqrt(tx*tx + ty*ty + tz*tz);
	    p[i].dtGrav = maga*dirsum/normsum + pkd->param.dPreFacRhoLoc*rholoc;
	    p[i].fDensity = pkd->param.dPreFacRhoLoc*rholoc;
	    }

/* 	if (p[i].iOrder < 100) { */
/* 	    printf("PID %d a0 %g a1 %g a2 %g\n",p[i].iOrder,p[i].a[0],p[i].a[1],p[i].a[2]); */
/* 	    } */

	} /* end of i-loop cells & particles */
    /*
    ** Free time-step lists
    */ 
    free(rholocal);
    *pdFlop += nActive*(nPart*40 + nCell*200) + nSoft*15;
    return(nActive);
    }
