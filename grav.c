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

#define CAcceptAngle 0.75

#define ECCFACMAX 10000

#define NMAXPLD 16

void HEAPheapstruct(int n, int k, HEAPSTRUCT ra[]) {
    int l,j,ir,i;
    HEAPSTRUCT rra;

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
	    if (j < ir && ra[j-1].rhoenc < ra[j].rhoenc) ++j; /* k largest elements at the end */
	    if (rra.rhoenc < ra[j-1].rhoenc) {
		ra[i-1] = ra[j-1];
		j += (i=j);
		}
	    else j = ir+1;
	    }
	ra[i-1] = rra;
	}
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
int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *pLoc,ILP *ilp,int nPart,ILC *ilc,int nCell,ILPB *ilpb,int nPartBucket,double *pdFlop) {
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
    double rhocadd,rhocaddlocal,rholoc,rhopmax,rhopmaxlocal,costheta;
    double vx,vy,vz,v2;
    momFloat tax,tay,taz,magai,adotai;
    double rhosum,maisum;
#ifdef HERMITE
    double adx,ady,adz;
    double dir5;
#endif
    double rv, a3, summ;
#ifdef SOFTSQUARE
    double ptwoh2;
#endif
    int i,j,k,l,nN,nC,nSP,nSC,nSCmin,nSPB,nSPBmin,na,nia,nSoft,nActive;
    /*
    ** time-step lists
    */
    RHOENC *rhoenc;
    RHOLOCAL *rholocal;	
    HEAPSTRUCT *heapstruct;
    int *jmax;
    /*
    ** Determine particle numbers and sizes of time-step lists
    */
    nC = nCell+nPartBucket;                     /* total number of cells */
    nN = nPart+pkdn->pUpper-pkdn->pLower;       /* total number of neighbouring particles (without particle itself) */
    nSC  = (int) (0.005*nCell);                 /* estimate of number of cells to check for maximum */
    nSPB = (int) (0.005*nPartBucket);           /* estimate of number of particle-buckets to check for maximum */
    nSCmin  = (nCell < 2)?nCell:2;              /* minimum number of cells to check for maximum */
    nSPBmin = (nPartBucket < 1)?nPartBucket:1;  /* minimum number of particle-buckets to check for maximum */
    nSC  = (nSC < nSCmin)?nSCmin:nSC;           /* final number of cells to check for maximum */
    nSPB = (nSPB < nSPBmin)?nSPBmin:nSPB;       /* final number of particle-buckets to check for maximum */
    /*
    ** Allocate time-step lists
    */
    rhoenc = malloc(nC*sizeof(RHOENC));
    assert(rhoenc != NULL);
    rholocal = malloc(nN*sizeof(RHOLOCAL));
    assert(rholocal != NULL);
    heapstruct = malloc(nC*sizeof(HEAPSTRUCT));
    assert(heapstruct != NULL);
    jmax = malloc((nSC+nSPB)*sizeof(int));
    assert(jmax != NULL);
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
	rhocadd = 0;
	rholoc = 0;
	rhopmax = 0;
	rhosum = 0;
	maisum = 0;
	
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
	    rhoenc[j].index = j;
	    rhoenc[j].x = x;
	    rhoenc[j].y = y;
	    rhoenc[j].z = z;
	    rhoenc[j].dir = dir;
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
	    rhoenc[j].rhoenc = dir2;
	    adotai = p[i].a[0]*tax + p[i].a[1]*tay + p[i].a[2]*taz; 
	    if (adotai >= 0) {
	      rhosum += adotai*rhoenc[j].dir;
	      maisum += sqrt(tax*tax + tay*tay + taz*taz);
	    }
	    ax += tax;
	    ay += tay;
	    az += taz;
	    } /* end of cell list gravity loop */
	
	if(pkd->param.bGravStep && pkd->param.iTimeStepCrit >= 0) {
	    /*
	    ** Add particle-bucket stuff
	    */
	    for(j=nCell,k=0;j<nC;++j,++k) {
		rhoenc[j].index = j;
		rhoenc[j].x = p[i].r[0] - ilpb[k].x;
		rhoenc[j].y = p[i].r[1] - ilpb[k].y;
		rhoenc[j].z = p[i].r[2] - ilpb[k].z;
		d2 = rhoenc[j].x*rhoenc[j].x + rhoenc[j].y*rhoenc[j].y + rhoenc[j].z*rhoenc[j].z; 
#ifdef SOFTSQUARE
		ptwoh2 = 2*p[i].fSoft*p[i].fSoft;
		fourh2 = ptwoh2 + ilpb[k].twoh2;
#endif
#ifdef SOFTLINEAR
		fourh2 = p[i].fSoft + ilpb[k].h;
		fourh2 *= fourh2;
#endif		       
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE) 
#ifdef SOFTENING_NOT_MASS_WEIGHTED
		fourh2 = ilpb[k].fourh2;
#else
		fourh2 = softmassweight(p[i].fMass,4*p[i].fSoft*p[i].fSoft,ilpb[k].m,ilpb[k].fourh2);
#endif
#endif
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
		    }

		rhoenc[j].dir = dir;
		dir2 *=  ilpb[k].m;
		rhoenc[j].rhoenc = dir2;
	    }
	    /*
	    ** Determine the nSC maximum cells in the cell list
	    */
	    if (nCell > 0) {
	      for (j=0;j<nCell;++j) {
		heapstruct[j].index = rhoenc[j].index;
		heapstruct[j].rhoenc = rhoenc[j].rhoenc;
	      }
	      HEAPheapstruct(nCell,nSC,heapstruct);
	      for (j=0;j<nSC;++j) {
		jmax[j] = heapstruct[nCell-1-j].index;
	      }
	    }
	    /*
	    ** Determine the nSPB maximum cells in the particle-bucket list
	    */
	    if (nPartBucket > 0) {
	      for (j=0;j<nPartBucket;++j) {
		heapstruct[j].index = rhoenc[nCell+j].index;
		heapstruct[j].rhoenc = rhoenc[nCell+j].rhoenc;
	      }
	      HEAPheapstruct(nPartBucket,nSPB,heapstruct);
	      for (j=0;j<nSPB;++j) {
		jmax[nSC+j] = heapstruct[nPartBucket-1-j].index;
	      }
	    }
	    /*
	    ** Determine the maximum of the rhocadd values
	    */
	    if (nC > 0) {
	      for (j=0;j<(nSC+nSPB);++j) {
		l = jmax[j];
		rhocaddlocal = 0;
		for (k=0;k<nC;++k) {
		  costheta = (rhoenc[l].x*rhoenc[k].x + rhoenc[l].y*rhoenc[k].y + rhoenc[l].z*rhoenc[k].z)*rhoenc[l].dir*rhoenc[k].dir;
		  if (costheta > CAcceptAngle && 2*rhoenc[k].dir > rhoenc[l].dir) rhocaddlocal += rhoenc[k].rhoenc;
		}
		rhocadd = (rhocaddlocal > rhocadd)?rhocaddlocal:rhocadd;
	      }
	    }
	    assert(rhocadd >= 0);
	    } /* end of cell & particle-bucket list time-step loop */
	
	mdlCacheCheck(pkd->mdl);
	
#ifdef SOFTSQUARE
	ptwoh2 = 2*p[i].fSoft*p[i].fSoft;
#endif
	for (j=0;j<nPart;++j) {
	    x = p[i].r[0] - ilp[j].x;
	    y = p[i].r[1] - ilp[j].y;
	    z = p[i].r[2] - ilp[j].z;
	    d2 = x*x + y*y + z*z;
#ifdef HERMITE
	    if(pkd->param.bHermite || pkd->param.iTimeStepCrit > 0){	 
	      vx = p[i].v[0] - ilp[j].vx;
	      vy = p[i].v[1] - ilp[j].vy;
	      vz = p[i].v[2] - ilp[j].vz;  
	      v2 = vx*vx + vy*vy + vz*vz;   
	      rv = x*vx + y*vy + z*vz; 
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
	    rholocal[j].m = ilp[j].m;	
	    rholocal[j].d2 = d2;
	    if (d2 > fourh2) {
		SQRT1(d2,dir);
		dir2 = dir*dir*dir;
		magai = dir*dir;
		}
	    else {		
#ifdef PLANETS 
		    if(pkd->param.bCollision){
		      p[i].iColflag = 1;
		      p[i].iOrderCol = ilp[j].iOrder;
		      p[i].dtCol = 1.0*p[i].iOrgIdx;	
		      printf("dr = %e, dr = %e, pi = %i, pj = %i, ilp  \n",sqrt(d2),sqrt(fourh2),p[i].iOrgIdx,ilp[j].iOrder);
		    }
#endif
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
	    /* GravStep 
	       iTimeStepCrit = 1: MZ eccentricity correction
	                       2: Normal
			       3: Planet
	    */	   
	    if(pkd->param.bGravStep && pkd->param.iTimeStepCrit > 0 && 
	       (ilp[j].iOrder < pkd->param.nPColl || p[i].iOrder < pkd->param.nPColl)) {

	      summ = p[i].fMass+ilp[j].m;
	      rhopmaxlocal = summ*dir2;

	      if((pkd->param.iTimeStepCrit == 1 || pkd->param.iTimeStepCrit == 3) && ilp[j].m > 0){
		a3 = p[i].r[0]*p[i].r[0]+p[i].r[1]*p[i].r[1]+p[i].r[2]*p[i].r[2];
		a3 = a3*sqrt(a3);
		rhopmaxlocal = pkdRho(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
				      vx,vy,vz,rv,v2,a3,p[i].iOrder,ilp[j].iOrder);
	      }
		rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;	
	    }

	    dir2 *= ilp[j].m;
	    tax = -x*dir2;
	    tay = -y*dir2;
	    taz = -z*dir2;
#if 0
	    adotai = p[i].a[0]*tax + p[i].a[1]*tay + p[i].a[2]*taz;
	    if (adotai >= 0) {
	      rhosum += adotai*dir;
	      maisum += ilp[j].m*magai;
	    }
#endif
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

	if(pkd->param.bGravStep && pkd->param.iTimeStepCrit >= 0) {
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
	    ** Add local density from particle list to mean field value.
	    ** Calculate local density only in the case of more than 1 neighbouring particle! 
	    */
	    if (nN > 1) {
		nSP = (nN < NMAXPLD)?nN:NMAXPLD;
		HEAPrholocal(nN,nSP,rholocal);
		SQRT1(rholocal[nN-nSP].d2,dir);
		dir *= dir;
		for (j=(nN - nSP);j<nN;++j) {
		    d2 = rholocal[j].d2*dir;
		    d2 =(1-d2)*(1-d2)*(1-d2)*(1-d2)*(1-d2);
		    rholoc += d2*rholocal[j].m; 
		    }
		rholoc = 9009/1024*M_1_PI*rholoc*dir*sqrt(dir); /* F5 Kernel */
		}
	    assert(rholoc >= 0);
	    rhocadd += rholoc;
	    p[i].dtGrav = (rhocadd > p[i].dtGrav)?rhocadd:p[i].dtGrav;
	    /*
	    ** Mean field value from cells including local density is now set.
	    ** Check if particle-particle interaction with nPColl particles (e.g. black holes) dominate!
	    */
	    if(pkd->param.iTimeStepCrit > 0) {
		p[i].dtGrav = (rhopmax > p[i].dtGrav)?rhopmax:p[i].dtGrav;
		}
	    } /* end of particle-particle list time-step loop */
	/*
	** Set the density value using the new timestepping formula.
	*/
	if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == -1 && rhosum > 0) {
	    p[i].dtGrav = rhosum/maisum;
	    }
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
	} /* end of i-loop over particles in the bucket */
    /*
    ** Free time-step lists
    */ 
    free(rhoenc);
    free(rholocal);
    free(heapstruct);
    free(jmax);
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
	       (pj->iOrder < pkd->param.nPColl || pi->iOrder < pkd->param.nPColl)) {
	      
	      summ = pi->fMass+pj->fMass; 
	      rhopmaxlocal = summ*dir2;
	      
	      if((pkd->param.iTimeStepCrit == 1 || pkd->param.iTimeStepCrit == 3) && pj->fMass> 0){
		a3 = pi->r[0]*pi->r[0]+pi->r[1]*pi->r[1]+pi->r[2]*pi->r[2];		
		a3 = a3*sqrt(a3);
		rhopmaxlocal = pkdRho(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
				      vx,vy,vz,rv,v2,a3,pi->iOrder,pj->iOrder);				
	      }
	    
	      pi->dtGrav = (rhopmaxlocal > pi->dtGrav)?rhopmaxlocal:pi->dtGrav;
	      pj->dtGrav = (rhopmaxlocal > pj->dtGrav)?rhopmaxlocal:pj->dtGrav; 
	    }

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
	       (pj->iOrder < pkd->param.nPColl || pi->iOrder < pkd->param.nPColl)) {
	      
	      summ = pi->fMass+pj->fMass; 
	      rhopmaxlocal = summ*dir2;
	      
	      if((pkd->param.iTimeStepCrit == 1 || pkd->param.iTimeStepCrit == 3) && pj->fMass> 0){
		a3 = pi->r[0]*pi->r[0]+pi->r[1]*pi->r[1]+pi->r[2]*pi->r[2];		
		a3 = a3*sqrt(a3);
		rhopmaxlocal = pkdRho(pkd,rhopmaxlocal,summ,sqrt(fourh2),&dir2,&dir,x,y,z,
				      vx,vy,vz,rv,v2,a3,pi->iOrder,pj->iOrder);			
  	
	      }
	      
	      pi->dtGrav = (rhopmaxlocal > pi->dtGrav)?rhopmaxlocal:pi->dtGrav;
	    }

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
    *pdFlop += nActive*(nPart*40 + nCell*200) + nSoft*15;
    return(nActive);
}
    
    
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

