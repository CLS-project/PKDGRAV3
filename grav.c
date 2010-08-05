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
#include "ewald.h"

const char *grav_c_module_id = "$Id$";
const char *grav_h_module_id = GRAV_H_MODULE_ID;

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
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,KDN *pBucket,LOCR *pLoc,ILP *ilp,int nPart,ILC *ilc,int nCell,
		    double dirLsum,double normLsum,int bEwald,double *pdFlop,double *pdEwFlop,double dRhoFac) {
    PARTICLE *p;
    float *a;
    double *vi,*vj;
    PARTICLE *pi,*pj;
    KDN *pkdn = pBucket;
    pBND kbnd;
    const double onethird = 1.0/3.0;
    float *pPot, fPot;
    double ax,ay,az;
    double x,y,z,d2,dir,dir2,g2,g3,g4;
    double xx,xy,xz,yy,yz,zz;
    double xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    double tx,ty,tz;
    double fourh2;
    double dtGrav,dT;
    double rholoc,rhopmax,rhopmaxlocal,dirDTS,d2DTS,dsmooth2,size2;
    momFloat tax,tay,taz,adotai,maga,dimaga,dirsum,normsum;
    double vx,vy,vz;
#ifdef HERMITE
    double v2;
    double adx,ady,adz;
    double dir5;
    double rv, a3;
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
    int i,j,k,nN,nSP,nTN,na,nia,nSoft,nActive;

    kbnd = pkdNodeBnd(pkd, pkdn);

    /*
    ** dynamical time-stepping stuff
    */
    RHOLOCAL *rholocal;
    nN = nPart+pkdn->pUpper-pkdn->pLower; /* Total number of particles on ilp list (without particle itself) Attention: not all are true neigbours! */
    rholocal = malloc(nN*sizeof(RHOLOCAL));
    assert(rholocal != NULL);
    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;
    for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
	p = pkdParticle(pkd,i);
	a = pkdAccel(pkd,p);
	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;
	++nActive;
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
	p->drmin = 10000.0;
	p->n_VA = 0; /* number of particles within 3 hill radius*/
#endif
	rhopmax = 0;
	dsmooth2 = 0;
	tx = a[0];
	ty = a[1];
	tz = a[2];
	maga = tx*tx + ty*ty + tz*tz;
	dimaga = maga;
	if (dimaga > 0) {
	    dimaga = 1.0/sqrt(dimaga);
	    }
	dirsum = dirLsum;
	normsum = normLsum;

	for (j=0;j<nCell;++j) {
	    x = p->r[0] - ilc[j].x;
	    y = p->r[1] - ilc[j].y;
	    z = p->r[2] - ilc[j].z;
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE
	    vx = p->v[0] - ilc[j].vx;
	    vy = p->v[1] - ilc[j].vy;
	    vz = p->v[2] - ilc[j].vz;
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
	    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
	    if (adotai > 0) {
		adotai *= dimaga;
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
	if (pkd->param.bGravStep) {
	    /*
	    ** Add particles to array rholocal first
	    */
	    rholoc = 0;
	    nTN = 0; /* Number of true neighbours, i.e. particles that fulfil the locality criterion defined by bnd.size */
	    size2 = pkd->param.dFacExcludePart * *kbnd.size; 
	    size2 *= size2; /* size2 estimate from tree.c that also works if there is only 1 particle in bucket */
	    for (j=0;j<nPart;++j) {
		x = p->r[0] - ilp[j].x;
		y = p->r[1] - ilp[j].y;
		z = p->r[2] - ilp[j].z;
		d2 = x*x + y*y + z*z;
		if (d2 > size2) continue;
		rholocal[nTN].m = ilp[j].m;
		rholocal[nTN].d2 = d2;
		nTN++;
		}
	    /*
            ** Add bucket particles to the array rholocal as well!
            ** Not including yourself!!
            */
            for (j=pkdn->pLower;j<=pkdn->pUpper;++j) {
		pj = pkdParticle(pkd,j);
		if (p->iOrder == pj->iOrder) continue;
		x = p->r[0] - pj->r[0];
		y = p->r[1] - pj->r[1];
		z = p->r[2] - pj->r[2];
		d2 = x*x + y*y + z*z;
		rholocal[nTN].m = pkdMass(pkd,pj);
		rholocal[nTN].d2 = d2;
                nTN++;
                }
	    /*
	    ** Calculate local density only in the case of more than 1 neighbouring particle!
	    */
	    if (nTN > 1) {
		nSP = (nTN < pkd->param.nPartRhoLoc)?nTN:pkd->param.nPartRhoLoc;
		HEAPrholocal(nTN,nSP,rholocal);
		dsmooth2 = rholocal[nTN-nSP].d2;
		SQRT1(dsmooth2,dir);
		dir2 = dir * dir;
		for (j=(nTN-nSP);j<nTN;++j) {
		    d2 = rholocal[j].d2*dir2;
		    d2 = (1-d2);
		    if (d2 < 0) d2 = 0.0;
		    rholoc += d2*rholocal[j].m;
		    }
		rholoc = 1.875*M_1_PI*rholoc*dir2*dir; /* F1 Kernel (15/8) */
		}
	    assert(rholoc >= 0);
	    }
#ifdef SOFTSQUARE
	ptwoh2 = 2*pkdSoft(pkd,p)*pkdSoft(pkd,p);
#endif
	for (j=0;j<nPart;++j) {
	    x = p->r[0] - ilp[j].x;
	    y = p->r[1] - ilp[j].y;
	    z = p->r[2] - ilp[j].z;
	    d2 = x*x + y*y + z*z;
	    d2DTS = d2;
#ifdef HERMITE
	    if (pkd->param.bHermite || pkd->param.iTimeStepCrit > 0) {
		vx = p->v[0] - ilp[j].vx;
		vy = p->v[1] - ilp[j].vy;
		vz = p->v[2] - ilp[j].vz;
		v2 = vx*vx + vy*vy + vz*vz;
		rv = x*vx + y*vy + z*vz;
		}
#endif
#ifdef SYMBA
	    a1 = p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2];
	    a2 = ilp[j].x*ilp[j].x + ilp[j].y*ilp[j].y + ilp[j].z*ilp[j].z;
	    a1 = 0.5*(sqrt(a1) + sqrt(a2));
	    a1 *= cbrt((pkdMass(pkd,p) + ilp[j].m)/(3.0*dSunMass));
	    a2 = sqrt(d2)/a1;
	    p->drmin = (a2 < p->drmin)?a2:p->drmin;
	    if (a2 < 3.0) {
		p->iOrder_VA[p->n_VA] = ilp[j].iOrder;
		p->hill_VA[p->n_VA] = a1;
		p->n_VA++;
		assert(p->iOrder != ilp[j].iOrder);
		if (a2 > rsym2) {
		    a1 = symfac*(1.0-a2/3.0); /*symfac is defined in pkd.h */
		    a1 *= a1*(2.0*a1 -3.0);
		    a1 += 1.0;
		    }
		else {
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
	    fourh2 = pkdSoft(pkd,p) + ilp[j].h;
	    fourh2 *= fourh2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
#ifdef SOFTENING_NOT_MASS_WEIGHTED
	    fourh2 = ilp[j].fourh2;
#else
	    fourh2 = softmassweight(pkdMass(pkd,p),4*pkdSoft(pkd,p)*pkdSoft(pkd,p),ilp[j].m,ilp[j].fourh2);
#endif
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {
#ifdef PLANETS
		if (pkd->param.bCollision) { /* have to use soft-linear */
		    pkd->iCollisionflag = 1; /*this is for veryactivehermite */
		    p->iColflag = 1;
		    p->iOrderCol = ilp[j].iOrder;
		    p->dtCol = 1.0*p->iOrgIdx;
		    printf("dr = %e, dr = %e, pi = %i, pj = %i, ilp  \n",
			   sqrt(d2),sqrt(fourh2),p->iOrgIdx,ilp[j].iOrder);
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
	    /*
	    ** GravStep if iTimeStepCrit =
	    ** 0: Mean field regime for dynamical time (normal/standard setting)
	    ** 1: Gravitational scattering regime for dynamical time with eccentricity correction
	    ** 2: Normal
	    ** 3: Planet
	    */
	    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 && 
		(p->iOrder < pkd->param.nPartColl || ilp[j].iOrder < pkd->param.nPartColl)) {
		summ = pkdMass(pkd,p)+ilp[j].m;
		rhopmaxlocal = summ*dir2;
		/*
		** Gravitational scattering regime (iTimeStepCrit=1)
		*/
		if (pkd->param.iTimeStepCrit == 1) {
		    vi = pkdVel(pkd,p);
		    vx = vi[0] - ilp[j].vx;
		    vy = vi[1] - ilp[j].vy;
		    vz = vi[2] - ilp[j].vz;
		    rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,x,y,z,vx,vy,vz,pkd->param.dEccFacMax);
		    }
#ifdef HERMITE
		if (pkd->param.iTimeStepCrit == 3 && ilp[j].m > 0) {
		    a3 = p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2];
		    a3 = a3*sqrt(a3);
		    rhopmaxlocal = pkdRho3(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
					  vx,vy,vz,rv,v2,a3,p->iOrder,ilp[j].iOrder);
		    }
#endif
		rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;
		}

	    dir2 *= ilp[j].m;
#ifdef SYMBA
	    if (a2 < 3.0) dir2 *= a1;
#endif
	    tax = -x*dir2;
	    tay = -y*dir2;
	    taz = -z*dir2;
	    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
	    if (adotai > 0 && d2DTS >= dsmooth2) {
		adotai *= dimaga;
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
	pPot = pkdPot(pkd,p);
	*pPot = fPot;
	a[0] = ax;
	a[1] = ay;
	a[2] = az;
#ifdef HERMITE
	p->ad[0] = adx;
	p->ad[1] = ady;
	p->ad[2] = adz;
#endif
	/*
	** Now finally calculate the Ewald correction for this particle, if it is
	** required.
	*/
	if (bEwald) {
	    *pdEwFlop += pkdParticleEwald(pkd,uRungLo,uRungHi,p);
	    }

	/*
	** Set value for time-step, note that we have the current ewald acceleration
	** in this now as well!
	*/
	if (pkd->param.bGravStep) {
	    /*
	    ** If this is the first time through, the accelerations will have 
	    ** all been zero resulting in zero for normsum (and nan for dtGrav).
	    ** We repeat this process again, so dtGrav will be correct.
	    */
	    if (normsum > 0.0) {
		/*
		** Use new acceleration here!
		*/
		tx = a[0];
		ty = a[1];
		tz = a[2];
		maga = sqrt(tx*tx + ty*ty + tz*tz);
		dtGrav = maga*dirsum/normsum;
		}
	    else dtGrav = 0.0;
	    dtGrav += pkd->param.dPreFacRhoLoc*rholoc;
	    if (pkd->param.iTimeStepCrit > 0) {
		dtGrav = (rhopmax > dtGrav)?rhopmax:dtGrav;
		}
	    if (dtGrav > 0.0) {
		dT = pkd->param.dEta/sqrt(dtGrav*dRhoFac);
		p->uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
		}
	    else p->uNewRung = 0;
	    p->fDensity = rholoc;
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
	p = pkdParticle(pkd,i);
	if (pkdIsDstActive(p,uRungLo,uRungHi)) pkd->piActive[na++] = p;
	else pkd->piInactive[nia++] = p;
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
	ptwoh2 = 2*pkdSoft(pkd,pi)*pkdSoft(pkd,pi);
#endif
	for (j=i+1;j<na;++j) {
	    pj = pkd->piActive[j];
	    x = pi->r[0] - pj->r[0];
	    y = pi->r[1] - pj->r[1];
	    z = pi->r[2] - pj->r[2];
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE
	    if (pkd->param.bHermite || pkd->param.iTimeStepCrit > 0) {
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
	    a1 *= cbrt((pkdMass(pkd,pi) + pkdMass(pkd,pj))/(3.0*dSunMass));
	    a2 = sqrt(d2)/a1;
	    pi->drmin = (a2 < pi->drmin)?a2:pi->drmin;
	    pj->drmin = (a2 < pj->drmin)?a2:pj->drmin;
	    if (a2 < 3.0) {
		/* printf("Active-Active iOrder %d jOrder %d \n",
		   pi->iOrder,pj->iOrder);*/
		pi->iOrder_VA[pi->n_VA] = pj->iOrder;
		pj->iOrder_VA[pj->n_VA] = pi->iOrder;
		pi->hill_VA[pi->n_VA] = a1;
		pj->hill_VA[pj->n_VA] = a1;
		pi->n_VA++;
		pj->n_VA++;
		assert(pj->iOrder!=pi->iOrder);
		if (a2 > rsym2) {
		    a1 = symfac*(1.0-a2/3.0);
		    a1 *= a1*(2.0*a1 -3.0);
		    a1 += 1.0;
		    }
		else {
		    a1 = 0.0;
		    }
		}
#endif
#ifdef SOFTSQUARE
	    fourh2 = ptwoh2 + 2*pkdSoft(pkd,pj)*pkdSoft(pkd,pj);
#endif
#ifdef SOFTLINEAR
	    fourh2 = pkdSoft(pkd,pi) + pkdSoft(pkd,pj);
	    fourh2*= fourh2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
#ifdef SOFTENING_NOT_MASS_WEIGHTED
	    fourh2 = 4*pkdSoft(pkd,pj)*pkdSoft(pkd,pj);
#else
	    fourh2 = softmassweight(pkdMass(pkd,pi),4*pkdSoft(pkd,pi)*pkdSoft(pkd,pi),pkdMass(pkd,pj),4*pkdSoft(pkd,pj)*pkdSoft(pkd,pj));
#endif
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {
#ifdef PLANETS
		if (pkd->param.bCollision) {
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
	    /* 
	    ** GravStep
	    */
	    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 &&
		(pi->iOrder < pkd->param.nPartColl || pj->iOrder < pkd->param.nPartColl)) {
		summ = pkdMass(pkd,pi)+pkdMass(pkd,pj);
		rhopmaxlocal = summ*dir2;
		/*
		** Gravitational scattering regime (iTimeStepCrit=1)
		*/
		if (pkd->param.iTimeStepCrit == 1) {
		    vi = pkdVel(pkd,pi);
		    vj = pkdVel(pkd,pj);
		    vx = vi[0] - vj[0];
		    vy = vi[1] - vj[1];
		    vz = vi[2] - vj[2];
		    rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,x,y,z,vx,vy,vz,pkd->param.dEccFacMax);
		    }
#ifdef HERMITE
		if (pkd->param.iTimeStepCrit == 3 && ilp[j].m > 0) {
		    a3 = p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2];
		    a3 = a3*sqrt(a3);
		    rhopmaxlocal = pkdRho3(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
					  vx,vy,vz,rv,v2,a3,p->iOrder,ilp[j].iOrder);
		    }
#endif
		if ( rhopmaxlocal > 0.0 ) {
		    dT = pkd->param.dEta/sqrt(rhopmaxlocal*dRhoFac);
		    uint8_t uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
		    if ( uNewRung > pi->uNewRung ) pi->uNewRung = uNewRung;
		    if ( uNewRung > pj->uNewRung ) pj->uNewRung = uNewRung;
		    }
		}

#ifdef SYMBA
	    if (a2 < 3.0) dir2 *= a1;
#endif
	    fPot += dir*pkdMass(pkd,pj);
	    ax += x*dir2*pkdMass(pkd,pj);
	    ay += y*dir2*pkdMass(pkd,pj);
	    az += z*dir2*pkdMass(pkd,pj);
	    pPot = pkdPot(pkd,pj);
	    *pPot -= dir*pkdMass(pkd,pi);
	    a = pkdAccel(pkd,pj);
	    a[0] += x*dir2*pkdMass(pkd,pi);
	    a[1] += y*dir2*pkdMass(pkd,pi);
	    a[2] += z*dir2*pkdMass(pkd,pi);
#ifdef HERMITE
	    dir5 = 3.0*rv*dir2*dir*dir;
	    adx += (vx*dir2-x*dir5)*pkdMass(pkd,pj);
	    ady += (vy*dir2-y*dir5)*pkdMass(pkd,pj);
	    adz += (vz*dir2-z*dir5)*pkdMass(pkd,pj);
	    pj->ad[0] += (vx*dir2-x*dir5)*pkdMass(pkd,pi);
	    pj->ad[1] += (vy*dir2-y*dir5)*pkdMass(pkd,pi);
	    pj->ad[2] += (vz*dir2-z*dir5)*pkdMass(pkd,pi);
#endif

	    } /* end of j-loop */
	pPot = pkdPot(pkd,pi);
	*pPot -= fPot;
	a = pkdAccel(pkd,pi);
	a[0] -= ax;
	a[1] -= ay;
	a[2] -= az;
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
	ptwoh2 = 2*pkdSoft(pkd,pi)*pkdSoft(pkd,pi);
#endif
	for (j=0;j<nia;++j) {
	    pj = pkd->piInactive[j];
	    x = pi->r[0] - pj->r[0];
	    y = pi->r[1] - pj->r[1];
	    z = pi->r[2] - pj->r[2];
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE
	    if (pkd->param.bHermite || pkd->param.iTimeStepCrit > 0) {
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
	    fourh2 = ptwoh2 + 2*pkdSoft(pkd,pj)*pkdSoft(pkd,pj);
#endif
#ifdef SOFTLINEAR
	    fourh2 = pkdSoft(pkd,pi) + pkdSoft(pkd,pj);
	    fourh2*= fourh2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
#ifdef SOFTENING_NOT_MASS_WEIGHTED
	    fourh2 = 4*pkdSoft(pkd,pj)*pkdSoft(pkd,pj);
#else
	    fourh2 = softmassweight(pkdMass(pkd,pi),4*pkdSoft(pkd,pi)*pkdSoft(pkd,pi),pkdMass(pkd,pj),4*pkdSoft(pkd,pj)*pkdSoft(pkd,pj));
#endif
#endif
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		}
	    else {
#ifdef PLANETS
		if (pkd->param.bCollision) {
		    pkd->iCollisionflag = 1; /*this is for veryactivehermite */
		    pi->iColflag = 1;
		    pi->iOrderCol = pj->iOrder;
		    pi->dtCol = 1.0*pi->iOrgIdx;
		    printf("dr = %e, r1+r2 = %e, pi = %i, pj = %i, active-inactive \n",sqrt(d2),sqrt(fourh2),pi->iOrgIdx,pj->iOrgIdx);
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
	    /* 
	    ** GravStep
	    */
	    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 &&
		(pi->iOrder < pkd->param.nPartColl || pj->iOrder < pkd->param.nPartColl)) {
		summ = pkdMass(pkd,pi)+pkdMass(pkd,pj);
		rhopmaxlocal = summ*dir2;
		/*
		** Gravitational scattering regime (iTimeStepCrit=1)
		*/
		if (pkd->param.iTimeStepCrit == 1) {
		    vi = pkdVel(pkd,pi);
		    vj = pkdVel(pkd,pj);
		    vx = vi[0] - vj[0];
		    vy = vi[1] - vj[1];
		    vz = vi[2] - vj[2];
		    rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,x,y,z,vx,vy,vz,pkd->param.dEccFacMax);
		    }
#ifdef HERMITE
		if (pkd->param.iTimeStepCrit == 3 && ilp[j].m > 0) {
		    a3 = p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2];
		    a3 = a3*sqrt(a3);
		    rhopmaxlocal = pkdRho3(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
					  vx,vy,vz,rv,v2,a3,p->iOrder,ilp[j].iOrder);
		    }
#endif
		if ( rhopmaxlocal > 0.0 ) {
		    dT = pkd->param.dEta/sqrt(rhopmaxlocal*dRhoFac);
		    uint8_t uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
		    if ( uNewRung > pi->uNewRung ) pi->uNewRung = uNewRung;
		    }
		}
#ifdef SYMBA
	    if (a2 < 3.0) dir2 *= a1;
#endif
	    fPot += dir*pkdMass(pkd,pj);
	    ax += x*dir2*pkdMass(pkd,pj);
	    ay += y*dir2*pkdMass(pkd,pj);
	    az += z*dir2*pkdMass(pkd,pj);
#ifdef HERMITE
	    dir5 = 3.0*rv*dir2*dir*dir;
	    adx += (vx*dir2-x*dir5)*pkdMass(pkd,pj);
	    ady += (vy*dir2-y*dir5)*pkdMass(pkd,pj);
	    adz += (vz*dir2-z*dir5)*pkdMass(pkd,pj);
#endif

	    } /* end of j-loop */
	pPot = pkdPot(pkd,pi);
	*pPot -= fPot;
	a = pkdAccel(pkd,pi);
	a[0] -= ax;
	a[1] -= ay;
	a[2] -= az;
#ifdef HERMITE
	pi->ad[0] -= adx;
	pi->ad[1] -= ady;
	pi->ad[2] -= adz;
#endif
	} /* end of i-loop active-inactive */
    *pdFlop += nActive*(nPart*40 + nCell*200) + nSoft*15;
    return(nActive);
    }

/*
** Gravitational scattering regime (iTimeStepCrit=1)
*/
double pkdRho1(double rhopmaxlocal, double summ, double dir, double x, double y, double z, double vx, double vy, double vz, double EccFacMax) {

    double Etot, L2, ecc, eccfac, v2;
    /*
    ** Etot and L are normalized by the reduced mass 
    */
    v2 = vx*vx + vy*vy + vz*vz;
    Etot = 0.5*v2 - summ*dir;
    L2 = (y*vz - z*vy)*(y*vz - z*vy) + (z*vx - x*vz)*(z*vx - x*vz) + (x*vy - y*vx)*(x*vy - y*vx);
    ecc = 1+2*Etot*L2/(summ*summ);
    ecc = (ecc <= 0)?0:sqrt(ecc);
    eccfac = (1 + 2*ecc)/fabs(1-ecc);
    eccfac = (eccfac > EccFacMax)?EccFacMax:eccfac;
    if (eccfac > 1.0) rhopmaxlocal *= eccfac;
    return rhopmaxlocal;
    }
/*
** Planet (iTimeStepCrit=3)
*/
#ifdef HERMITE
#ifdef PLANETS
double pkdRho3(PKD pkd, double rhopmaxlocal,double summ, double sumr, double *dir2,
	      double *dir, double x, double y, double z, double vx, double vy, double vz,
	      double rv, double v2, double a3, int iOrder,int jOrder) {

    double hill, vd, vesc2;
    double dr, dr3, dr5;
    double fhill = 3.0;
    double rhof;
    
    hill = a3*summ/3.0;
    if (fhill*fhill*fhill*(*dir2)*hill > 1.0) {
	hill = pow(hill, 1.0/3.0);
	
	rhopmaxlocal = 3.0*(summ*(*dir2)+ 0.5*sumr*rv*rv*(*dir2)*(*dir2)/(*dir));
	
	vesc2 = 2.0*summ/sumr;
	rhof = (rv*(*dir))*(rv*(*dir))/vesc2;
	if (rhof < 3.0) rhof = 3.0;
	
	rhopmaxlocal *= sqrt((*dir)*hill);
	a3 /= rhof;
	if (rhopmaxlocal*a3 > 1.0) {
	    dr = rv*(*dir)*pkd->param.dEta/sqrt(rhopmaxlocal);
	    }
	else {
	    dr = rv*(*dir)*pkd->param.dEta*sqrt(a3);
	    }
	
	dr = 1.0/(*dir)+0.5*dr;
	
	if (dr < 0.0) {
	    vesc2 = 2.0*summ/sumr;
	    printf("pi %d, pj %d, dr0/rhill %e, dr1/rhill %e, v/vesc %e, \n",
		   iOrder,jOrder,1.0/((*dir)*hill), dr/hill, sqrt(v2/vesc2));
	    dr = fhill*hill;
	    *dir = 1.0/dr;
	    *dir2 = (*dir)*(*dir)*(*dir);
	    
	    }
	
	/* if mutual distance < fhill Hill radius */
	if (dr < fhill*hill) {
	    dr3 = dr*dr*dr*sqrt(dr/hill);
	    dr5 = dr3*dr*dr;
	    rhopmaxlocal = 3.0*(summ/dr3 + 0.5*sumr*rv*rv/dr5);
	    /*rhopmaxlocal = rhof*summ/dr3;*/
	    rhopmaxlocal = (rhopmaxlocal > 1.0/a3)?rhopmaxlocal:(1.0/a3);
	    }
	}
    return(rhopmaxlocal);
    }
#endif
#endif

