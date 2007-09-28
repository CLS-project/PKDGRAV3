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
#include "grav.h"

#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    } 

#define ECCFACMAX 10000

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
		    int bEwald,int nReps,double fEwCut,double *pdFlop,double *pdEwFlop) {
    PARTICLE *p = pkd->pStore;
    PARTICLE *pi,*pj;
    KDN *pkdn = pBucket;
    const double onethird = 1.0/3.0;
    double ax,ay,az,fPot;
    double x,y,z,d2,dir,dir2,g2,g3,g4;
    double xx,xy,xz,yy,yz,zz;
    double xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    double tx,ty,tz;
    double fourh2;
    double rholoc,rhopmax,rhopmaxlocal,dirDTS,d2DTS,dsmooth2;   
    momFloat tax,tay,taz,adotai,maga,dirsum,normsum;
    /*
    ** Ewald Variables
    */
    MOMC mom = pkd->momRoot;
    double dx,dy,dz,r2,g0,g1,g5,xzz,yzz,zzz;
    double Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
    double L,fEwCut2,fInner2,alpha,alpha2,alphan,k1,ka,a;
    double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;
    double hdotx,s,c;
    int ix,iy,iz,nEwReps,bInHole,bInHolex,bInHolexy,nLoop;
#ifdef HERMITE
    double adx,ady,adz;
    double dir5;
    double rv, a3;
    double vx,vy,vz,v2;
#endif
    double summ;
#ifdef PLANETS
    double dSunMass = pkd->dSunMass;
#ifdef SYMBA
    double drmin, a1, a2;
#endif
#endif
#ifdef SOFTSQUARE
    double ptwoh2;
#endif
    int i,j,k,l,nN,nSP,na,nia,nSoft,nActive;

    /*
    ** Ewald Stuff (if needed)
    */
    if (bEwald) {
	/*
	 ** Set up traces of the complete multipole moments.
	 */
	Q4xx = 0.5*(mom.xxxx + mom.xxyy + mom.xxzz);
	Q4xy = 0.5*(mom.xxxy + mom.xyyy + mom.xyzz);
	Q4xz = 0.5*(mom.xxxz + mom.xyyz + mom.xzzz);
	Q4yy = 0.5*(mom.xxyy + mom.yyyy + mom.yyzz);
	Q4yz = 0.5*(mom.xxyz + mom.yyyz + mom.yzzz);
	Q4zz = 0.5*(mom.xxzz + mom.yyzz + mom.zzzz);
	Q4 = 0.25*(Q4xx + Q4yy + Q4zz);
	Q3x = 0.5*(mom.xxx + mom.xyy + mom.xzz);
	Q3y = 0.5*(mom.xxy + mom.yyy + mom.yzz);
	Q3z = 0.5*(mom.xxz + mom.yyz + mom.zzz);
	Q2 = 0.5*(mom.xx + mom.yy + mom.zz);	

	nEwReps = ceil(fEwCut);
	L = pkd->fPeriod[0];
	fEwCut2 = fEwCut*fEwCut*L*L;
	fInner2 = 3.0e-3*L*L;
	nEwReps = nEwReps > nReps ? nEwReps : nReps;
	alpha = 2.0/L;
	alpha2 = alpha*alpha;
	k1 = M_PI/(alpha2*L*L*L);
	ka = 2.0*alpha/sqrt(M_PI);
    }
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
#ifdef HERMITE
	adx = 0;
	ady = 0;
	adz = 0;
#endif
#ifdef SYMBA
	p[i].drmin = 10000.0;
	p[i].n_VA = 0; /* number of particles within 3 hill radius*/
#endif
	rhopmax = 0;
	rholoc = 0;
	dsmooth2 = 0;
	maga = sqrt(p[i].a[0]*p[i].a[0] + p[i].a[1]*p[i].a[1] + p[i].a[2]*p[i].a[2]);
	dirsum = dirLsum;
	normsum = normLsum;
	
	for (j=0;j<nCell;++j) {
	    x = p[i].r[0] - ilc[j].x;
	    y = p[i].r[1] - ilc[j].y;
	    z = p[i].r[2] - ilc[j].z;
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE	   
	    vx = p[i].v[0] - ilc[j].vx;
	    vy = p[i].v[1] - ilc[j].vy;
	    vz = p[i].v[2] - ilc[j].vz;
	    rv = x*vx + y*vy + z*vz;
#endif	    	   
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
	    dir *= ilc[j].mom.m;
#ifdef HERMITE
            /* contribution to ad only from monopole */
	    adx -= vx*dir2*dir;
	    ady -= vy*dir2*dir;
	    adz -= vz*dir2*dir;
	    dir5  = 3.0*rv*dir2*dir2*dir;
	    adx += x*dir5;
	    ady += y*dir5;
	    adz += z*dir5;
#endif
	    dir2 *= dir + 5*g2 + 7*g3 + 9*g4;
	    fPot -= dir + g2 + g3 + g4;
	    tax = xx + xxx + tx - x*dir2;
	    tay = xy + xxy + ty - y*dir2;
	    taz = xz + xxz + tz - z*dir2;
	    adotai = p[i].a[0]*tax + p[i].a[1]*tay + p[i].a[2]*taz; 
	    if (adotai > 0) {
		adotai /= maga;
		dirsum += dirDTS*adotai*adotai;
		normsum += adotai*adotai;
		}
	    ax += tax;
	    ay += tay;
	    az += taz;
	    } /* end of cell list gravity loop */
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
#ifdef SOFTSQUARE
	ptwoh2 = 2*p[i].fSoft*p[i].fSoft;
#endif
	for (j=0;j<nPart;++j) {
	    x = p[i].r[0] - ilp[j].x;
	    y = p[i].r[1] - ilp[j].y;
	    z = p[i].r[2] - ilp[j].z;
	    d2 = x*x + y*y + z*z;
	    d2DTS = d2;
#ifdef HERMITE
	    if(pkd->param.bHermite || pkd->param.iTimeStepCrit > 0){	 
		vx = p[i].v[0] - ilp[j].vx;
		vy = p[i].v[1] - ilp[j].vy;
		vz = p[i].v[2] - ilp[j].vz;  
		v2 = vx*vx + vy*vy + vz*vz;   
		rv = x*vx + y*vy + z*vz; 
		}
#endif
#ifdef SYMBA
	    a1 = p[i].r[0]*p[i].r[0]+p[i].r[1]*p[i].r[1]+p[i].r[2]*p[i].r[2];
	    a2 = ilp[j].x*ilp[j].x + ilp[j].y*ilp[j].y + ilp[j].z*ilp[j].z;
	    a1 = 0.5*(sqrt(a1) + sqrt(a2));
	    a1 *= cbrt((p[i].fMass + ilp[j].m)/(3.0*dSunMass));
	    a2 = sqrt(d2)/a1;
	    p[i].drmin = (a2 < p[i].drmin)?a2:p[i].drmin;
	    if(a2 < 3.0){
		p[i].iOrder_VA[p[i].n_VA] = ilp[j].iOrder;
		p[i].hill_VA[p[i].n_VA] = a1;
		p[i].n_VA++;
		assert(p[i].iOrder != ilp[j].iOrder);
		if(a2 > rsym2){
		    a1 = symfac*(1.0-a2/3.0); /*symfac is defined in pkd.h */ 
		    a1 *= a1*(2.0*a1 -3.0);
		    a1 += 1.0;
		    }else{ 
			a1 = 0.0;
			/* d3 = a2*a2;
			   d3 *= (9.0 -5.0*d3*d3)/(108.0*a1*a1*a1);*/
			}
		}
#endif 
#ifdef SOFTSQUARE
	    fourh2 = ptwoh2 + ilp[j].twoh2;
#endif
#ifdef SOFTLINEAR
	    fourh2 = p[i].fSoft + ilp[j].h;
	    fourh2 *= fourh2;
#endif		       
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE) 
#ifdef SOFTENING_NOT_MASS_WEIGHTED
	    fourh2 = ilp[j].fourh2;
#else
	    fourh2 = softmassweight(p[i].fMass,4*p[i].fSoft*p[i].fSoft,ilp[j].m,ilp[j].fourh2);
#endif
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {		
#ifdef PLANETS 
		if(pkd->param.bCollision){ /* have to use soft-linear */
		    pkd->iCollisionflag = 1; /*this is for veryactivehermite */
		    p[i].iColflag = 1;
		    p[i].iOrderCol = ilp[j].iOrder;		       
		    p[i].dtCol = 1.0*p[i].iOrgIdx;	
		    printf("dr = %e, dr = %e, pi = %i, pj = %i, ilp  \n",
			   sqrt(d2),sqrt(fourh2),p[i].iOrgIdx,ilp[j].iOrder);
		    }
#endif
		/*
		** This uses the Dehnen K1 kernel function now, it's fast!
		*/	
		SQRT1(fourh2,dir);
		dir2 = dir*dir;
		d2 *= dir2;
		dir2 *= dir;
		d2 = 1.0 - d2;
		dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
		dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
		++nSoft;
		}
	    /* GravStep 
	       iTimeStepCrit = 1: MZ eccentricity correction
	       2: Normal
	       3: Planet
	    */	   
	    if(pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 && 
	       (ilp[j].iOrder < pkd->param.nPartColl || p[i].iOrder < pkd->param.nPartColl)) {
		
		summ = p[i].fMass+ilp[j].m;
		rhopmaxlocal = summ*dir2;
#ifdef HERMITE	 
		if((pkd->param.iTimeStepCrit == 1 || pkd->param.iTimeStepCrit == 3) && ilp[j].m > 0){
		    a3 = p[i].r[0]*p[i].r[0]+p[i].r[1]*p[i].r[1]+p[i].r[2]*p[i].r[2];
		    a3 = a3*sqrt(a3);
		    rhopmaxlocal = pkdRho(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
					  vx,vy,vz,rv,v2,a3,p[i].iOrder,ilp[j].iOrder);
		    }
#endif
		rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;	
		}

	    dir2 *= ilp[j].m;
#ifdef SYMBA
	    if(a2 < 3.0) dir2 *= a1; 
#endif
	    tax = -x*dir2;
	    tay = -y*dir2;
	    taz = -z*dir2;
	    adotai = p[i].a[0]*tax + p[i].a[1]*tay + p[i].a[2]*taz;
	    if (adotai > 0 && d2DTS >= dsmooth2) {
		adotai /= maga;
		dirsum += dir*adotai*adotai;
		normsum += adotai*adotai;
		}
	    fPot -= ilp[j].m*dir;
	    ax += tax;
	    ay += tay;
	    az += taz;
#ifdef HERMITE
	    adx -= vx*dir2;
	    ady -= vy*dir2;
	    adz -= vz*dir2;
	    dir5  = 3.0*rv*dir2*dir*dir;
	    adx += x*dir5;
	    ady += y*dir5;
	    adz += z*dir5;
#endif
	    } /* end of particle list gravity loop */
	/*
	** Finally set new acceleration and potential.
	** Note that after this point we cannot use the new timestepping criterion since we
	** overwrite the acceleration.
	*/
	p[i].fPot = fPot;
	p[i].a[0] = ax;
	p[i].a[1] = ay;
	p[i].a[2] = az;
#ifdef HERMITE
	p[i].ad[0] = adx;
	p[i].ad[1] = ady;
	p[i].ad[2] = adz;
#endif
	/*
	** Now finally calculate the Ewald correction for this particle, if it is
	** required.
	*/
	if (bEwald) {
	    fPot = mom.m*k1;
	    ax = 0.0;
	    ay = 0.0;
	    az = 0.0;
	    dx = p[i].r[0] - pkd->kdTop[ROOT].r[0];
	    dy = p[i].r[1] - pkd->kdTop[ROOT].r[1];
	    dz = p[i].r[2] - pkd->kdTop[ROOT].r[2];
	    for (ix=-nEwReps;ix<=nEwReps;++ix) {
		bInHolex = (ix >= -nReps && ix <= nReps);
		x = dx + ix*L;
		for(iy=-nEwReps;iy<=nEwReps;++iy) {
		    bInHolexy = (bInHolex && iy >= -nReps && iy <= nReps);
		    y = dy + iy*L;
		    for(iz=-nEwReps;iz<=nEwReps;++iz) {
			bInHole = (bInHolexy && iz >= -nReps && iz <= nReps);
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
			if (r2 > fEwCut2 && !bInHole) continue;
			if (r2 < fInner2) {
			    /*
			     * For small r, series expand about
			     * the origin to avoid errors caused
			     * by cancellation of large terms.
			     */
			    alphan = ka;
			    r2 *= alpha2;
			    g0 = alphan*((1.0/3.0)*r2 - 1.0);
			    alphan *= 2*alpha2;
			    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
			    alphan *= 2*alpha2;
			    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
			    alphan *= 2*alpha2;
			    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
			    alphan *= 2*alpha2;
			    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
			    alphan *= 2*alpha2;
			    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
			}
			else {
			    dir = 1/sqrt(r2);
			    dir2 = dir*dir;
			    a = exp(-r2*alpha2);
			    a *= ka*dir2;
			    if (bInHole) {
				g0 = -erf(alpha/dir);
			    }
			    else {
				g0 = erfc(alpha/dir);
			    }
			    g0 *= dir;
			    g1 = g0*dir2 + a;
			    alphan = 2*alpha2;
			    g2 = 3*g1*dir2 + alphan*a;
			    alphan *= 2*alpha2;
			    g3 = 5*g2*dir2 + alphan*a;
			    alphan *= 2*alpha2;
			    g4 = 7*g3*dir2 + alphan*a;
			    alphan *= 2*alpha2;
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
			Q4x = Q4xx*x + Q4xy*y + Q4xz*z;
			Q4y = Q4xy*x + Q4yy*y + Q4yz*z;
			Q4z = Q4xz*x + Q4yz*y + Q4zz*z;
			Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (Q3x*x + Q3y*y + Q3z*z) + Q4;
			Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
			Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
			Qta = g1*mom.m - g2*Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
			fPot -= g0*mom.m - g1*Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
			ax += g2*(Q2mirx - Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
			ay += g2*(Q2miry - Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
			az += g2*(Q2mirz - Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
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
	    for (j=0;j<pkd->nEwhLoop;++j) {
		hdotx = pkd->ewt[j].hx*dx + pkd->ewt[j].hy*dy + pkd->ewt[j].hz*dz;
		c = cos(hdotx);
		s = sin(hdotx);
		fPot += pkd->ewt[j].hCfac*c + pkd->ewt[j].hSfac*s;
		ax += pkd->ewt[j].hx*(pkd->ewt[j].hCfac*s - pkd->ewt[j].hSfac*c);
		ay += pkd->ewt[j].hy*(pkd->ewt[j].hCfac*s - pkd->ewt[j].hSfac*c);
		az += pkd->ewt[j].hz*(pkd->ewt[j].hCfac*s - pkd->ewt[j].hSfac*c);
	    }
	    p[i].fPot += fPot;
	    p[i].ae[0] = ax;
	    p[i].ae[1] = ay;
	    p[i].ae[2] = az;
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
	    if(pkd->param.iTimeStepCrit > 0) {
		p[i].dtGrav = (rhopmax > p[i].dtGrav)?rhopmax:p[i].dtGrav;
		}
	    }
	} /* end of i-loop over particles in the bucket */
    /*
    ** Free time-step lists
    */ 
    free(rholocal);
    /*
    ** Perform intra-bucket interactions
    */
    na = 0;
    nia = 0;
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	if (pkdIsActive(pkd,&p[i])) pkd->piActive[na++] = &p[i];
	else pkd->piInactive[nia++] = &p[i];
    }
    /*
    ** Active-Active Interactions.
    */
    for (i=0;i<na-1;++i) {
	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
#ifdef HERMITE
	adx = 0;
	ady = 0;
	adz = 0;
#endif
	pi = pkd->piActive[i];
#ifdef SOFTSQUARE
	ptwoh2 = 2*pi->fSoft*pi->fSoft;
#endif
	for (j=i+1;j<na;++j) {
	    pj = pkd->piActive[j];
	    x = pi->r[0] - pj->r[0];
	    y = pi->r[1] - pj->r[1];
	    z = pi->r[2] - pj->r[2];
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE
	    if(pkd->param.bHermite || pkd->param.iTimeStepCrit > 0){
		vx = pi->v[0] - pj->v[0];
		vy = pi->v[1] - pj->v[1];
		vz = pi->v[2] - pj->v[2];
		v2 = vx*vx + vy*vy + vz*vz;
		rv = x*vx+ y*vy + z*vz;
		}
#endif    
#ifdef SYMBA
	    a1 = pi->r[0]*pi->r[0]+pi->r[1]*pi->r[1]+pi->r[2]*pi->r[2];
	    a2 = pj->r[0]*pj->r[0]+pj->r[1]*pj->r[1]+pj->r[2]*pj->r[2]; 
	    a1 = 0.5*(sqrt(a1) + sqrt(a2));
	    a1 *= cbrt((pi->fMass + pj->fMass)/(3.0*dSunMass));
	    a2 = sqrt(d2)/a1;
	    pi->drmin = (a2 < pi->drmin)?a2:pi->drmin; 
	    pj->drmin = (a2 < pj->drmin)?a2:pj->drmin; 
	    if(a2 < 3.0){
		/* printf("Active-Active iOrder %d jOrder %d \n", 
		   pi->iOrder,pj->iOrder);*/
		pi->iOrder_VA[pi->n_VA] = pj->iOrder;
		pj->iOrder_VA[pj->n_VA] = pi->iOrder;
		pi->hill_VA[pi->n_VA] = a1;
		pj->hill_VA[pj->n_VA] = a1;
		pi->n_VA++;
		pj->n_VA++;
		assert(pj->iOrder!=pi->iOrder);
		if(a2 > rsym2){
		    a1 = symfac*(1.0-a2/3.0);
		    a1 *= a1*(2.0*a1 -3.0);
		    a1 += 1.0;
		}else{ 
		    a1 = 0.0;		 
		}
	    }
#endif   
#ifdef SOFTSQUARE
	    fourh2 = ptwoh2 + 2*pj->fSoft*pj->fSoft;
#endif
#ifdef SOFTLINEAR
	    fourh2 = pi->fSoft + pj->fSoft;
	    fourh2*= fourh2;
#endif		       
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE) 
#ifdef SOFTENING_NOT_MASS_WEIGHTED
	    fourh2 = 4*pj->fSoft*pj->fSoft;
#else
	    fourh2 = softmassweight(pi->fMass,4*pi->fSoft*pi->fSoft,pj->fMass,4*pj->fSoft*pj->fSoft);
#endif
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {	
#ifdef PLANETS 
	 if(pkd->param.bCollision){	
	     pkd->iCollisionflag = 1; /*this is for veryactivehermite */
	     pi->iColflag = 1;
	     pi->iOrderCol = pj->iOrder;	    
	     pi->dtCol = 1.0*pi->iOrgIdx;
	     printf("dr = %e, r1+r2 = %e, pi = %i, pj = %i active-active \n",
		    sqrt(d2),sqrt(fourh2),pi->iOrgIdx,pj->iOrgIdx);       	 
	 }
#endif	 
	        /*
		** This uses the Dehnen K1 kernel function now, it's fast!
		** (no more lookup tables)
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

	    /* GravStep */	
	    if(pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 && 
	       (pj->iOrder < pkd->param.nPartColl || pi->iOrder < pkd->param.nPartColl)) {
	      
		summ = pi->fMass+pj->fMass; 
		rhopmaxlocal = summ*dir2;
#ifdef HERMITE	  	      
		if((pkd->param.iTimeStepCrit == 1 || pkd->param.iTimeStepCrit == 3) && pj->fMass> 0){
		    a3 = pi->r[0]*pi->r[0]+pi->r[1]*pi->r[1]+pi->r[2]*pi->r[2];		
		    a3 = a3*sqrt(a3);
		    rhopmaxlocal = pkdRho(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
					  vx,vy,vz,rv,v2,a3,pi->iOrder,pj->iOrder);				
		    }
#endif	    
		pi->dtGrav = (rhopmaxlocal > pi->dtGrav)?rhopmaxlocal:pi->dtGrav;
		pj->dtGrav = (rhopmaxlocal > pj->dtGrav)?rhopmaxlocal:pj->dtGrav; 
		}

#ifdef SYMBA	  
	    if(a2 < 3.0) dir2 *= a1;
#endif
	    fPot += dir*pj->fMass;
	    ax += x*dir2*pj->fMass;
	    ay += y*dir2*pj->fMass;
	    az += z*dir2*pj->fMass;
	    pj->fPot -= dir*pi->fMass;
	    pj->a[0] += x*dir2*pi->fMass;
	    pj->a[1] += y*dir2*pi->fMass;
	    pj->a[2] += z*dir2*pi->fMass;
#ifdef HERMITE	    
	    dir5 = 3.0*rv*dir2*dir*dir;
	    adx += (vx*dir2-x*dir5)*pj->fMass;
	    ady += (vy*dir2-y*dir5)*pj->fMass;
	    adz += (vz*dir2-z*dir5)*pj->fMass;
	    pj->ad[0] += (vx*dir2-x*dir5)*pi->fMass;
	    pj->ad[1] += (vy*dir2-y*dir5)*pi->fMass;
	    pj->ad[2] += (vz*dir2-z*dir5)*pi->fMass;
#endif

	    } /* end of j-loop */
	pi->fPot -= fPot;
	pi->a[0] -= ax;
	pi->a[1] -= ay;
	pi->a[2] -= az;
#ifdef HERMITE
	pi->ad[0] -= adx;
	pi->ad[1] -= ady;
	pi->ad[2] -= adz;
#endif
	} /* end of i-loop active-active */
    /*
    ** Active-Inactive interactions.
    */
    for (i=0;i<na;++i) {
	fPot = 0;
	ax = 0;
	ay = 0;
	az = 0;
#ifdef HERMITE
	adx = 0;
	ady = 0;
	adz = 0;
#endif
	pi = pkd->piActive[i];
#ifdef SOFTSQUARE
	ptwoh2 = 2*pi->fSoft*pi->fSoft;
#endif
	for (j=0;j<nia;++j) {
	    pj = pkd->piInactive[j];
	    x = pi->r[0] - pj->r[0];
	    y = pi->r[1] - pj->r[1];
	    z = pi->r[2] - pj->r[2];
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE
	    if(pkd->param.bHermite || pkd->param.iTimeStepCrit > 0){
	      vx = pi->v[0] - pj->v[0];
	      vy = pi->v[1] - pj->v[1];
	      vz = pi->v[2] - pj->v[2];    
	      rv = vx*x + vy*y + vz*z;
	      v2 = vx*vx + vy*vy + vz*vz;  
	    }
#endif	   
#ifdef SYMBA 
	    assert(0); /* all particles should be active for Symba*/
#endif 
#ifdef SOFTSQUARE
	    fourh2 = ptwoh2 + 2*pj->fSoft*pj->fSoft;
#endif
#ifdef SOFTLINEAR
	    fourh2 = pi->fSoft + pj->fSoft;
	    fourh2*= fourh2;
#endif		       
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE) 
#ifdef SOFTENING_NOT_MASS_WEIGHTED
	    fourh2 = 4*pj->fSoft*pj->fSoft;
#else
	    fourh2 = softmassweight(pi->fMass,4*pi->fSoft*pi->fSoft,pj->fMass,4*pj->fSoft*pj->fSoft);
#endif
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {
		/*
		** This uses the Dehnen K1 kernel function now, it's fast!
		** (no more lookup tables)
		*/
#ifdef PLANETS 
		if(pkd->param.bCollision){
		    pkd->iCollisionflag = 1; /*this is for veryactivehermite */
		    pi->iColflag = 1;
		    pi->iOrderCol = pj->iOrder;		   
		    pi->dtCol = 1.0*pi->iOrgIdx;	
		    printf("dr = %e, r1+r2 = %e, pi = %i, pj = %i, active-inactive \n",sqrt(d2),sqrt(fourh2),pi->iOrgIdx,pj->iOrgIdx);
		}
#endif
	
		SQRT1(fourh2,dir);
		dir2 = dir*dir;
		d2 *= dir2;
		dir2 *= dir;
		d2 = 1 - d2;
		dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
		dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
		++nSoft;
	    }

	    /* GravStep */	
	    if(pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 && 
	       (pj->iOrder < pkd->param.nPartColl || pi->iOrder < pkd->param.nPartColl)) {
	      
		summ = pi->fMass+pj->fMass; 
		rhopmaxlocal = summ*dir2;
	      
#ifdef HERMITE	      
		if((pkd->param.iTimeStepCrit == 1 || pkd->param.iTimeStepCrit == 3) && pj->fMass> 0){
		    a3 = pi->r[0]*pi->r[0]+pi->r[1]*pi->r[1]+pi->r[2]*pi->r[2];		
		    a3 = a3*sqrt(a3);
		    rhopmaxlocal = pkdRho(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
					  vx,vy,vz,rv,v2,a3,pi->iOrder,pj->iOrder);			
		    }
#endif	      
		pi->dtGrav = (rhopmaxlocal > pi->dtGrav)?rhopmaxlocal:pi->dtGrav;
		}
#ifdef SYMBA	  
	    if(a2 < 3.0) dir2 *= a1;
#endif
	    fPot += dir*pj->fMass;
	    ax += x*dir2*pj->fMass;
	    ay += y*dir2*pj->fMass;
	    az += z*dir2*pj->fMass;
#ifdef HERMITE	    
	    dir5 = 3.0*rv*dir2*dir*dir;
	    adx += (vx*dir2-x*dir5)*pj->fMass;
	    ady += (vy*dir2-y*dir5)*pj->fMass;
	    adz += (vz*dir2-z*dir5)*pj->fMass;
#endif

	} /* end of j-loop */
	pi->fPot -= fPot;
	pi->a[0] -= ax;
	pi->a[1] -= ay;
	pi->a[2] -= az;
#ifdef HERMITE
	pi->ad[0] -= adx;
	pi->ad[1] -= ady;
	pi->ad[2] -= adz;
#endif	
	} /* end of i-loop active-inactive */
    *pdEwFlop += nLoop*447 + nActive*pkd->nEwhLoop*58;
    *pdFlop += nActive*(nPart*40 + nCell*200) + nSoft*15;
    return(nActive);
}
    
#ifdef HERMITE
double pkdRho(PKD pkd, double rhopmaxlocal,double summ, double sumr, double *dir2, 
	      double *dir, double x, double y, double z, double vx, double vy, double vz, 
	      double rv, double v2, double a3, int iOrder,int jOrder){
  
  if(pkd->param.iTimeStepCrit == 1){	
    double Etot, L2, ecc, eccfac;
    /* Etot and L are normalized by the reduced mass */
    Etot = 0.5*v2 - summ*(*dir);
    L2 = (y*vz - z*vy)*(y*vz - z*vy) + (z*vx - x*vz)*(z*vx - x*vz) + 
      (x*vy - y*vx)*(x*vy - y*vx);	      
    ecc = 1+2*Etot*L2/summ/summ;
    ecc = (ecc <= 0)?0:sqrt(ecc);	
    eccfac = (2*ecc + 1)/fabs(1-ecc);
    eccfac = (eccfac > ECCFACMAX)?ECCFACMAX:eccfac;
    if(eccfac > 1.0) rhopmaxlocal *= eccfac; 
  }
#ifdef PLANETS 
  if(pkd->param.iTimeStepCrit == 3){		  
    double hill, vd, vesc2;
    double dr, dr3, dr5;
    double fhill = 3.0;
    double rhof;
    
    hill = a3*summ/3.0;		 
    if(fhill*fhill*fhill*(*dir2)*hill > 1.0){
      hill = pow(hill, 1.0/3.0);
      
      rhopmaxlocal = 3.0*(summ*(*dir2)+ 0.5*sumr*rv*rv*(*dir2)*(*dir2)/(*dir)); 
 
      vesc2 = 2.0*summ/sumr;     
      rhof = (rv*(*dir))*(rv*(*dir))/vesc2;
      if (rhof < 3.0) rhof = 3.0;

      rhopmaxlocal *= sqrt((*dir)*hill);
      a3 /= rhof;
      if(rhopmaxlocal*a3 > 1.0){
	dr = rv*(*dir)*pkd->param.dEta/sqrt(rhopmaxlocal);
      }else{
	dr = rv*(*dir)*pkd->param.dEta*sqrt(a3);
      }
      
      dr = 1.0/(*dir)+0.5*dr;
      
      if(dr < 0.0){
	vesc2 = 2.0*summ/sumr;
	printf("pi %d, pj %d, dr0/rhill %e, dr1/rhill %e, v/vesc %e, \n",
	       iOrder,jOrder,1.0/((*dir)*hill), dr/hill, sqrt(v2/vesc2)); 
	dr = fhill*hill;
	*dir = 1.0/dr;
	*dir2 = (*dir)*(*dir)*(*dir);

      }
      
      /* if mutual distance < fhill Hill radius */  
      if(dr < fhill*hill){
	dr3 = dr*dr*dr*sqrt(dr/hill);
	dr5 = dr3*dr*dr;	   			 
	rhopmaxlocal = 3.0*(summ/dr3 + 0.5*sumr*rv*rv/dr5);
	  /*rhopmaxlocal = rhof*summ/dr3;*/
	rhopmaxlocal = (rhopmaxlocal > 1.0/a3)?rhopmaxlocal:(1.0/a3);		   
      }
    }
  }
#endif
  return(rhopmaxlocal);
}

#endif
