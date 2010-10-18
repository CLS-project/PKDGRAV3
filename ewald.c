#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include "ewald.h"
#include "pkd.h"
#include "meval.h"
#include "qeval.h"
#include "moments.h"

int pkdParticleEwald(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,PARTICLE *p) {
    MOMC mom = pkd->momRoot;
    float *pPot;
    double fPot,ax,ay,az;
    double dx,dy,dz,x,y,z,r2,dir,dir2,a,alphan,L;
    double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
    double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;
    double g0,g1,g2,g3,g4,g5;
    double onethird = 1.0/3.0;
    double hdotx,s,c;
    float *pa;
    int i,ix,iy,iz,bInHole,bInHolex,bInHolexy;
    int nFlop;
    int nLoop = 0;

    assert(pkd->oAcceleration); /* Validate memory model */
    assert(pkd->oPotential); /* Validate memory model */

    if (!pkdIsDstActive(p,uRungLo,uRungHi)) return 0;
    pa = pkdAccel(pkd,p);
    pPot = pkdPot(pkd,p);

    L = pkd->fPeriod[0];
    fPot = mom.m*pkd->ew.k1;
    ax = 0.0;
    ay = 0.0;
    az = 0.0;
    dx = p->r[0] - pkdTopNode(pkd,ROOT)->r[0];
    dy = p->r[1] - pkdTopNode(pkd,ROOT)->r[1];
    dz = p->r[2] - pkdTopNode(pkd,ROOT)->r[2];
    for (ix=-pkd->ew.nEwReps;ix<=pkd->ew.nEwReps;++ix) {
	bInHolex = (ix >= -pkd->ew.nReps && ix <= pkd->ew.nReps);
	x = dx + ix*L;
	for (iy=-pkd->ew.nEwReps;iy<=pkd->ew.nEwReps;++iy) {
	    bInHolexy = (bInHolex && iy >= -pkd->ew.nReps && iy <= pkd->ew.nReps);
	    y = dy + iy*L;
	    for (iz=-pkd->ew.nEwReps;iz<=pkd->ew.nEwReps;++iz) {
		bInHole = (bInHolexy && iz >= -pkd->ew.nReps && iz <= pkd->ew.nReps);
		/*
		** Scoring for Ewald inner stuff = (+,*)
		**		Visible ops 		= (104,161)
		**		sqrt, 1/sqrt est. 	= (6,11)
		**     division            = (6,11)  same as sqrt.
		**		exp est.			= (6,11)  same as sqrt.
		**		erf/erfc est.		= (12,22) twice a sqrt.
		**		Total			= (128,205) = 333
		**     Old scoring				    = 447
		*/
		z = dz + iz*L;
		r2 = x*x + y*y + z*z;
		if (r2 > pkd->ew.fEwCut2 && !bInHole) continue;
		if (r2 < pkd->ew.fInner2) {
		    /*
		     * For small r, series expand about
		     * the origin to avoid errors caused
		     * by cancellation of large terms.
		     */
		    alphan = pkd->ew.ka;
		    r2 *= pkd->ew.alpha2;
		    g0 = alphan*((1.0/3.0)*r2 - 1.0);
		    alphan *= 2*pkd->ew.alpha2;
		    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
		    }
		else {
		    dir = 1/sqrt(r2);
		    dir2 = dir*dir;
		    a = exp(-r2*pkd->ew.alpha2);
		    a *= pkd->ew.ka*dir2;
		    if (bInHole) {
			g0 = -erf(pkd->ew.alpha/dir);
			}
		    else {
			g0 = erfc(pkd->ew.alpha/dir);
			}
		    g0 *= dir;
		    g1 = g0*dir2 + a;
		    alphan = 2*pkd->ew.alpha2;
		    g2 = 3*g1*dir2 + alphan*a;
		    alphan *= 2*pkd->ew.alpha2;
		    g3 = 5*g2*dir2 + alphan*a;
		    alphan *= 2*pkd->ew.alpha2;
		    g4 = 7*g3*dir2 + alphan*a;
		    alphan *= 2*pkd->ew.alpha2;
		    g5 = 9*g4*dir2 + alphan*a;
		    }
		xx = 0.5*x*x;
		xxx = onethird*xx*x;
		xxy = xx*y;
		xxz = xx*z;
		yy = 0.5*y*y;
		yyy = onethird*yy*y;
		xyy = yy*x;
		yyz = yy*z;
		zz = 0.5*z*z;
		zzz = onethird*zz*z;
		xzz = zz*x;
		yzz = zz*y;
		xy = x*y;
		xyz = xy*z;
		xz = x*z;
		yz = y*z;
		Q2mirx = mom.xx*x + mom.xy*y + mom.xz*z;
		Q2miry = mom.xy*x + mom.yy*y + mom.yz*z;
		Q2mirz = mom.xz*x + mom.yz*y + mom.zz*z;
		Q3mirx = mom.xxx*xx + mom.xxy*xy + mom.xxz*xz + mom.xyy*yy + mom.xyz*yz + mom.xzz*zz;
		Q3miry = mom.xxy*xx + mom.xyy*xy + mom.xyz*xz + mom.yyy*yy + mom.yyz*yz + mom.yzz*zz;
		Q3mirz = mom.xxz*xx + mom.xyz*xy + mom.xzz*xz + mom.yyz*yy + mom.yzz*yz + mom.zzz*zz;
		Q4mirx = mom.xxxx*xxx + mom.xxxy*xxy + mom.xxxz*xxz + mom.xxyy*xyy + mom.xxyz*xyz +
			 mom.xxzz*xzz + mom.xyyy*yyy + mom.xyyz*yyz + mom.xyzz*yzz + mom.xzzz*zzz;
		Q4miry = mom.xxxy*xxx + mom.xxyy*xxy + mom.xxyz*xxz + mom.xyyy*xyy + mom.xyyz*xyz +
			 mom.xyzz*xzz + mom.yyyy*yyy + mom.yyyz*yyz + mom.yyzz*yzz + mom.yzzz*zzz;
		Q4mirz = mom.xxxz*xxx + mom.xxyz*xxy + mom.xxzz*xxz + mom.xyyz*xyy + mom.xyzz*xyz +
			 mom.xzzz*xzz + mom.yyyz*yyy + mom.yyzz*yyz + mom.yzzz*yzz + mom.zzzz*zzz;
		Q4x = pkd->ew.Q4xx*x + pkd->ew.Q4xy*y + pkd->ew.Q4xz*z;
		Q4y = pkd->ew.Q4xy*x + pkd->ew.Q4yy*y + pkd->ew.Q4yz*z;
		Q4z = pkd->ew.Q4xz*x + pkd->ew.Q4yz*y + pkd->ew.Q4zz*z;
		Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (pkd->ew.Q3x*x + pkd->ew.Q3y*y + pkd->ew.Q3z*z) + pkd->ew.Q4;
		Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
		Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
		Qta = g1*mom.m - g2*pkd->ew.Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
		fPot -= g0*mom.m - g1*pkd->ew.Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
		ax += g2*(Q2mirx - pkd->ew.Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
		ay += g2*(Q2miry - pkd->ew.Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
		az += g2*(Q2mirz - pkd->ew.Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
		++nLoop;
		}
	    }
	}
    /*
    ** Scoring for the h-loop (+,*)
    ** 	Without trig = (10,14)
    **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
    **		Total        = (22,36)
    **					 = 58
    */
    for (i=0;i<pkd->ew.nEwhLoop;++i) {
	hdotx = pkd->ew.ewt[i].hx*dx + pkd->ew.ewt[i].hy*dy + pkd->ew.ewt[i].hz*dz;
	c = cos(hdotx);
	s = sin(hdotx);
	fPot += pkd->ew.ewt[i].hCfac*c + pkd->ew.ewt[i].hSfac*s;
	ax += pkd->ew.ewt[i].hx*(pkd->ew.ewt[i].hCfac*s - pkd->ew.ewt[i].hSfac*c);
	ay += pkd->ew.ewt[i].hy*(pkd->ew.ewt[i].hCfac*s - pkd->ew.ewt[i].hSfac*c);
	az += pkd->ew.ewt[i].hz*(pkd->ew.ewt[i].hCfac*s - pkd->ew.ewt[i].hSfac*c);
	}
    *pPot += fPot;
    pa[0] += ax;
    pa[1] += ay;
    pa[2] += az;
    nFlop = nLoop*447 + pkd->ew.nEwhLoop*58;
    return(nFlop);
    }

void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut) {
    MOMC mom = pkd->momRoot;
    int i,hReps,hx,hy,hz,h2;
    double k4,L;
    double gam[6],mfacc,mfacs;
    double ax,ay,az;
    const int iOrder = 4;

    /*
    ** Set up traces of the complete multipole moments.
    */
    pkd->ew.Q4xx = 0.5*(mom.xxxx + mom.xxyy + mom.xxzz);
    pkd->ew.Q4xy = 0.5*(mom.xxxy + mom.xyyy + mom.xyzz);
    pkd->ew.Q4xz = 0.5*(mom.xxxz + mom.xyyz + mom.xzzz);
    pkd->ew.Q4yy = 0.5*(mom.xxyy + mom.yyyy + mom.yyzz);
    pkd->ew.Q4yz = 0.5*(mom.xxyz + mom.yyyz + mom.yzzz);
    pkd->ew.Q4zz = 0.5*(mom.xxzz + mom.yyzz + mom.zzzz);
    pkd->ew.Q4 = 0.25*(pkd->ew.Q4xx + pkd->ew.Q4yy + pkd->ew.Q4zz);
    pkd->ew.Q3x = 0.5*(mom.xxx + mom.xyy + mom.xzz);
    pkd->ew.Q3y = 0.5*(mom.xxy + mom.yyy + mom.yzz);
    pkd->ew.Q3z = 0.5*(mom.xxz + mom.yyz + mom.zzz);
    pkd->ew.Q2 = 0.5*(mom.xx + mom.yy + mom.zz);
    pkd->ew.nReps = nReps;
    pkd->ew.nEwReps = ceil(fEwCut);
    L = pkd->fPeriod[0];
    pkd->ew.fEwCut2 = fEwCut*fEwCut*L*L;
    pkd->ew.fInner2 = 1.2e-3*L*L;
    pkd->ew.nEwReps = pkd->ew.nEwReps > nReps ? pkd->ew.nEwReps : nReps;
    pkd->ew.alpha = 2.0/L;
    pkd->ew.alpha2 = pkd->ew.alpha*pkd->ew.alpha;
    pkd->ew.k1 = M_PI/(pkd->ew.alpha2*L*L*L);
    pkd->ew.ka = 2.0*pkd->ew.alpha/sqrt(M_PI);
    /*
    ** Now setup stuff for the h-loop.
    */
    hReps = ceil(fhCut);
    k4 = M_PI*M_PI/(pkd->ew.alpha*pkd->ew.alpha*L*L);
    i = 0;
    for (hx=-hReps;hx<=hReps;++hx) {
	for (hy=-hReps;hy<=hReps;++hy) {
	    for (hz=-hReps;hz<=hReps;++hz) {
		h2 = hx*hx + hy*hy + hz*hz;
		if (h2 == 0) continue;
		if (h2 > fhCut*fhCut) continue;
		if (i == pkd->ew.nMaxEwhLoop) {
		    pkd->ew.nMaxEwhLoop *= 2;
		    pkd->ew.ewt = realloc(pkd->ew.ewt,pkd->ew.nMaxEwhLoop*sizeof(EWT));
		    assert(pkd->ew.ewt != NULL);
		    }
		gam[0] = exp(-k4*h2)/(M_PI*h2*L);
		gam[1] = 2*M_PI/L*gam[0];
		gam[2] = -2*M_PI/L*gam[1];
		gam[3] = 2*M_PI/L*gam[2];
		gam[4] = -2*M_PI/L*gam[3];
		gam[5] = 2*M_PI/L*gam[4];
		gam[1] = 0.0;
		gam[3] = 0.0;
		gam[5] = 0.0;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		mfacc = 0.0;
		QEVAL(iOrder,pkd->momRoot,gam,hx,hy,hz,ax,ay,az,mfacc);
		gam[0] = exp(-k4*h2)/(M_PI*h2*L);
		gam[1] = 2*M_PI/L*gam[0];
		gam[2] = -2*M_PI/L*gam[1];
		gam[3] = 2*M_PI/L*gam[2];
		gam[4] = -2*M_PI/L*gam[3];
		gam[5] = 2*M_PI/L*gam[4];
		gam[0] = 0.0;
		gam[2] = 0.0;
		gam[4] = 0.0;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		mfacs = 0.0;
		QEVAL(iOrder,pkd->momRoot,gam,hx,hy,hz,ax,ay,az,mfacs);
		pkd->ew.ewt[i].hx = 2*M_PI/L*hx;
		pkd->ew.ewt[i].hy = 2*M_PI/L*hy;
		pkd->ew.ewt[i].hz = 2*M_PI/L*hz;
		pkd->ew.ewt[i].hCfac = mfacc;
		pkd->ew.ewt[i].hSfac = mfacs;
		++i;
		}
	    }
	}
    pkd->ew.nEwhLoop = i;
    }


