#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stddef.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include "pkd.h"
#include "moments.h"
#include "grav.h"

#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    } 

double softmassweight(double m1,double fourh12,double m2,double fourh22){
    return((m1+m2)*(fourh12*fourh22)/(fourh22*m1+fourh12*m2));
    }

/*
** This version of grav.c does all the operations inline, including 
** v_sqrt's and such.
** Returns nActive.
*/
int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *pLoc,ILP *ilp,int nPart,ILC *ilc,int nCell,ILPB *ilpb,int nPartBucket,double *pdFlop) {
    PARTICLE *p = pkd->pStore;
    PARTICLE *pi,*pj;
    KDN *pkdn = pBucket;
    const double onethird = 1.0/3.0;
    momFloat ax,ay,az,fPot;
    double x,y,z,d2,dir,dir2;
    double rhosum,maisum;
    momFloat magai,adotai;
    momFloat tax,tay,taz,tmon;
#ifndef USE_SIMD
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
    int i,j,k,l,na,nia,nSoft,nActive;
#ifdef USE_SIMD_MOMR
    nCellILC = nCell;
    momPadSIMDMomr( &nCellILC, ilc );
#endif
    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	if (!TYPEQueryACTIVE(&p[i])) continue;
	++nActive;
	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
	rhosum = 0;
	maisum = 0;
	/*
	** Evaluate local expansion.
	*/
	x = p[i].r[0] - pkdn->r[0];
	y = p[i].r[1] - pkdn->r[1];
	z = p[i].r[2] - pkdn->r[2];
	momEvalLocr(pLoc,x,y,z,&fPot,&ax,&ay,&az);

#ifdef USE_SIMD_MOMR
	momEvalSIMDMomr( nCellILC, ilc, p[i].r, &ax, &ay, &az, &fPot );
#else
	for (j=0;j<nCell;++j) {
	    x = p[i].r[0] - ilc[j].x;
	    y = p[i].r[1] - ilc[j].y;
	    z = p[i].r[2] - ilc[j].z;
	    d2 = x*x + y*y + z*z;
	    SQRT1(d2,dir);
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
	    adotai = p[i].a[0]*tax + p[i].a[1]*tay + p[i].a[2]*taz; 
	    if (adotai >= 0) {
	      rhosum += adotai*dir;
	      maisum += sqrt(tax*tax + tay*tay + taz*taz);
	    }
	    ax += tax;
	    ay += tay;
	    az += taz;
	    } /* end of cell list gravity loop */
#endif		
	mdlCacheCheck(pkd->mdl);
	
#ifdef SOFTSQUARE
	ptwoh2 = 2*p[i].fSoft*p[i].fSoft;
#endif
	for (j=0;j<nPart;++j) {
	    x = p[i].r[0] - ilp[j].x;
	    y = p[i].r[1] - ilp[j].y;
	    z = p[i].r[2] - ilp[j].z;
	    d2 = x*x + y*y + z*z;
#ifdef SOFTSQUARE
	    fourh2 = ptwoh2 + ilp[j].twoh2;
#endif
#ifdef SOFTLINEAR
	    fourh2 = p[i].fSoft + ilp[j].h;
	    fourh2 *= fourh2;
#endif		       
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE) 
	    /* fourh2 = ilp[j].fourh2; old softening */
	    fourh2 = softmassweight(p[i].fMass,4*p[i].fSoft*p[i].fSoft,ilp[j].m,ilp[j].fourh2);
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		magai = dir*dir;
	    }
	    else {
		/*
		** This uses the Dehnen K1 kernel function now, it's fast!
		*/
		magai = sqrt(d2);
		SQRT1(fourh2,dir);
		dir2 = dir*dir;
		d2 *= dir2;
		dir2 *= dir;
		d2 = 1 - d2;
		dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
		dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
		magai *= dir2;
		++nSoft;
		}
	    dir2 *= ilp[j].m;
	    tax = -x*dir2;
	    tay = -y*dir2;
	    taz = -z*dir2;

	    adotai = p[i].a[0]*tax + p[i].a[1]*tay + p[i].a[2]*taz;
	    if (adotai >= 0) {
	      rhosum += adotai*dir;
	      maisum += ilp[j].m*magai;
	    }

	    fPot -= ilp[j].m*dir;
	    ax += tax;
	    ay += tay;
	    az += taz;
	    } /* end of particle list gravity loop */
	/*
	** Set the density value using the new timestepping formula.
	*/
	if (pkd->param.bGravStep && rhosum > 0) {
	    p[i].dtGrav = rhosum/maisum;
	    }

        /*
        ** Finally set new acceleration and potential.
        ** Note that after this point we cannot use the new timestepping criteri
on since we
        ** overwrite the acceleration.
        */
	p[i].fPot = fPot;
	p[i].a[0] = ax;
	p[i].a[1] = ay;
	p[i].a[2] = az;	
	} /* end of i-loop cells & particles */
    *pdFlop += nActive*(nPart*40 + nCell*200) + nSoft*15;
    return(nActive);
    }
