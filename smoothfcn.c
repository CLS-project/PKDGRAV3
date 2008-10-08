#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *smoothfcn_module_id = "$Id$";

#include <math.h>
#include <assert.h>
#include "smoothfcn.h"

#ifdef M43D
/* M43D Creates a 3D kernel by convolution of 3D tophats the way M4(1D) is made in 1D */
#define BALL2(a) ((a)->fBall*(a)->fBall)
#define KERNEL(ak,ar2) { \
		ak = sqrt(ar2); \
		if (ar2 < 1.0) ak = 6.*0.25/350./3. *(1360+ar2*(-2880 \
			 +ar2*(3528+ak*(-1890+ak*(-240+ak*(270-6*ar2)))))); \
		else if (ar2 < 4.0) ak = 6.*0.25/350./3. *(7040-1152/ak+ak*(-10080+ak*(2880+ak*(4200 \
			 +ak*(-3528+ak*(630+ak*(240+ak*(-90+2*ar2)))))))); \
		else ak = 0.0;\
		}
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) adk = 6.*0.25/350./3. * (-2880*2 \
			 +ar2*(3528*4+ adk*(-1890*5 + adk*(-240*6+ adk*(270*7-6*9*ar2))))); \
		else if (ar2 < 4.0) adk = 6.*0.25/350./3. *((1152/ar2-10080)/adk+(2880*2+adk*(4200*3 \
			 +adk*(-3528*4+adk*(630*5+adk*(240*6 +adk*(-90*7+2*9*ar2))))))); \
		else adk = 0.0;\
		}

#else
#ifdef HSHRINK
/* HSHRINK M4 Kernel uses an effective h of (pi/6)^(1/3) times h for nSmooth neighbours */
#define dSHRINKFACTOR        0.805995977
#define BALL2(a) ((a)->fBall*(a)->fBall*(dSHRINKFACTOR*dSHRINKFACTOR))
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
		else ak = 0.0; \
		}
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else if (ar2 < 4.0) { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
		else adk = 0.0; \
		}

#else
/* Standard M_4 Kernel */
#define BALL2(a) ((a)->fBall*(a)->fBall)
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
		else ak = 0.0;\
		}
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else if (ar2 < 4.0) { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
		else adk = 0.0;\
		}
#endif
#endif

void NullSmooth(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    }

void initDensity(void *vpkd, void *p) {
    ((PARTICLE *)p)->fDensity = 0.0;
    }

void combDensity(void *vpkd, void *p1,void *p2) {
    ((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
    }

void Density(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    float ih2,r2,rs,fDensity,fMass;
    int i;

    ih2 = 4.0/BALL2(p);
    fDensity = 0.0;
    for (i=0;i<nSmooth;++i) {
	fMass = pkdMass(pkd,nnList[i].pPart);
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	fDensity += rs*fMass;
	}
    p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity;
    }

void DensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    float fNorm,ih2,r2,rs,fMassQ,fMassP;
    int i;
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	fMassQ = pkdMass(pkd,q);
	p->fDensity += rs*fMassQ;
	q->fDensity += rs*fMassP;
	}
    }


void initMeanVel(void *vpkd, void *pvoid) {
    PKD pkd = (PKD)vpkd;
    assert(pkd);
    PARTICLE *p = pvoid;
    VELSMOOTH *pvel = pkdField(p,pkd->oVelSmooth);
    int j;
    for(j=0;j<3;++j) pvel->vmean[j] = 0.0;
    }

void combMeanVel(void *vpkd, void *p1void,void *p2void) {
    PKD pkd = (PKD)vpkd;
    assert(pkd);
    PARTICLE *p1 = p1void;
    PARTICLE *p2 = p2void;
    VELSMOOTH *p1vel = pkdField(p1,pkd->oVelSmooth);
    VELSMOOTH *p2vel = pkdField(p2,pkd->oVelSmooth);
    int j;

    for (j=0;j<3;++j) p1vel->vmean[j] += p2vel->vmean[j];
    }

void MeanVel(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    double *qv;
    VELSMOOTH *pvel;
    float v[3];
    float ih2,r2,rs,fMass;
    int i,j;

    pvel = pkdField(p,pkd->oVelSmooth);
    ih2 = 4.0/BALL2(p);
    for (j=0;j<3;++j) v[j] = 0.0;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	q = nnList[i].pPart;
	fMass = pkdMass(pkd,q);
	qv = pkdVel(pkd,q);
	for (j=0;j<3;++j) v[j] += rs*fMass/q->fDensity*qv[j];
	}
    for (j=0;j<3;++j) pvel->vmean[j] = M_1_PI*sqrt(ih2)*ih2*v[j];
    }

void MeanVelSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    VELSMOOTH *pvel, *qvel;
    double *pv, *qv;
    float fNorm,ih2,r2,rs,fMassQ,fMassP;
    int i,j;

    pv = pkdVel(pkd,p);
    pvel = pkdField(p,pkd->oVelSmooth);
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	qv = pkdVel(pkd,q);
	qvel = pkdField(q,pkd->oVelSmooth);
	fMassQ = pkdMass(pkd,q);
	for (j=0;j<3;++j) {
	    pvel->vmean[j] += rs*fMassQ/q->fDensity*qv[j];
	    qvel->vmean[j] += rs*fMassP/p->fDensity*pv[j];
	    }
	}
    }



void initDivv(void *vpkd, void *pvoid) {
    PKD pkd = (PKD)vpkd;
    assert(pkd);
    PARTICLE *p = pvoid;
    VELSMOOTH *pvel = pkdField(p,pkd->oVelSmooth);
    pvel->divv = 0.0;
    }

void combDivv(void *vpkd, void *p1void,void *p2void) {
    PKD pkd = (PKD)vpkd;
    assert(pkd);
    PARTICLE *p1 = p1void;
    PARTICLE *p2 = p2void;
    VELSMOOTH *p1vel = pkdField(p1,pkd->oVelSmooth);
    VELSMOOTH *p2vel = pkdField(p2,pkd->oVelSmooth);
    p1vel->divv += p2vel->divv;
    }

void Divv(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    double *qv;
    double *pv;
    VELSMOOTH *pvel;
    float fNorm,ih2,r2,rs,fMass,dvdotdr;
    int i;

    pvel = pkdField(p,pkd->oVelSmooth);
    pv = pkdVel(pkd,p);
    ih2 = 4.0/BALL2(p);
    fNorm = M_1_PI*ih2*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	DKERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	fMass = pkdMass(pkd,q);
	qv = pkdVel(pkd,q);
	dvdotdr = (qv[0] - pv[0])*nnList[i].dx + (qv[1] - pv[1])*nnList[i].dy + (qv[2] - pv[2])*nnList[i].dz;
	pvel->divv += rs*fMass/q->fDensity*dvdotdr;
	}
    }

void DivvSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    VELSMOOTH *pvel, *qvel;
    double *pv, *qv;
    float fNorm,ih2,r2,rs,fMassQ,fMassP,dvdotdr;
    int i;

    pvel = pkdField(p,pkd->oVelSmooth);
    pv = pkdVel(pkd,p);
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*ih2*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	DKERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	qv = pkdVel(pkd,q);
	qvel = pkdField(q,pkd->oVelSmooth);
	fMassQ = pkdMass(pkd,q);
	dvdotdr = (qv[0] - pv[0])*nnList[i].dx + (qv[1] - pv[1])*nnList[i].dy + (qv[2] - pv[2])*nnList[i].dz;
	pvel->divv += rs*fMassQ/q->fDensity*dvdotdr;
	qvel->divv += rs*fMassP/p->fDensity*dvdotdr;
	}
    }



void initVelDisp2(void *vpkd, void *pvoid) {
    PKD pkd = (PKD)vpkd;
    assert(pkd);
    PARTICLE *p = pvoid;
    VELSMOOTH *pvel = pkdField(p,pkd->oVelSmooth);
    pvel->veldisp2 = 0.0;
    }

void combVelDisp2(void *vpkd, void *p1void,void *p2void) {
    PKD pkd = (PKD)vpkd;
    assert(pkd);
    PARTICLE *p1 = p1void;
    PARTICLE *p2 = p2void;
    VELSMOOTH *p1vel = pkdField(p1,pkd->oVelSmooth);
    VELSMOOTH *p2vel = pkdField(p2,pkd->oVelSmooth);
    p1vel->veldisp2 += p2vel->veldisp2;
    }

void VelDisp2(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    double *qv;
    double *pv;
    VELSMOOTH *pvel;
    float fNorm,ih2,r2,rs,fMass,tv,tv2;
    int i;

    pvel = pkdField(p,pkd->oVelSmooth);
    pv = pkdVel(pkd,p);
    ih2 = 4.0/BALL2(p);
    fNorm = M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	fMass = pkdMass(pkd,q);
	qv = pkdVel(pkd,q);

	tv = qv[0] - pvel->vmean[0] - pvel->divv*nnList[i].dx;
	tv2 = tv*tv;
	tv = qv[1] - pvel->vmean[1] - pvel->divv*nnList[i].dy;
	tv2 += tv*tv;
	tv = qv[2] - pvel->vmean[2] - pvel->divv*nnList[i].dz;
	tv2 += tv*tv;

	pvel->veldisp2 += rs*fMass/q->fDensity*tv2;
	}
    }

void VelDisp2Sym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    VELSMOOTH *pvel, *qvel;
    double *pv, *qv;
    float fNorm,ih2,r2,rs,fMassQ,fMassP,tv,tv2;
    int i;

    pvel = pkdField(p,pkd->oVelSmooth);
    pv = pkdVel(pkd,p);
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	qv = pkdVel(pkd,q);
	qvel = pkdField(q,pkd->oVelSmooth);
	fMassQ = pkdMass(pkd,q);

	tv = qv[0] - pvel->vmean[0] - pvel->divv*nnList[i].dx;
	tv2 = tv*tv;
	tv = qv[1] - pvel->vmean[1] - pvel->divv*nnList[i].dy;
	tv2 += tv*tv;
	tv = qv[2] - pvel->vmean[2] - pvel->divv*nnList[i].dz;
	tv2 += tv*tv;

	pvel->veldisp2 += rs*fMassQ/q->fDensity*tv2;
	qvel->veldisp2 += rs*fMassP/p->fDensity*tv2;
	}
    }



void initGroupMerge(void *vpkd, void *g) {
    }
void combGroupMerge(void *vpkd, void *g1, void *g2) {
    FOFGD * gd1 = (FOFGD *)g1;
    FOFGD * gd2 = (FOFGD *)g2;

    /* accept second (remote) group number, if its smaller */
    if (gd1->iGlobalId < gd2->iGlobalId)
	gd2->iGlobalId = gd1->iGlobalId;
    else gd1->iGlobalId = gd2->iGlobalId;
    /* if someone says that this not my group, accept it */
    gd1->bMyGroup *= gd2->bMyGroup;
    gd2->bMyGroup *= gd1->bMyGroup;
    }
void initGroupBins(void *vpkd, void *b) {
    FOFBIN * gb1 = (FOFBIN *)b;

    gb1->nMembers = 0;
    gb1->fMassInBin = 0.0;
    gb1->fMassEnclosed = 0.0;
    gb1->v2[0] = 0.0;
    gb1->v2[1] = 0.0;
    gb1->v2[2] = 0.0;
    gb1->L[0] = 0.0;
    gb1->L[1] = 0.0;
    gb1->L[2] = 0.0;
    }
void combGroupBins(void *vpkd, void *b1, void *b2) {
    FOFBIN * gb1 = (FOFBIN *)b1;
    FOFBIN * gb2 = (FOFBIN *)b2;

    /* add entries */
    gb1->nMembers += gb2->nMembers;
    gb1->fMassInBin += gb2->fMassInBin;
    gb1->v2[0] += gb2->v2[0];
    gb1->v2[1] += gb2->v2[1];
    gb1->v2[2] += gb2->v2[2];
    gb1->L[0]  += gb2->L[0];
    gb1->L[1]  += gb2->L[1];
    gb1->L[2]  += gb2->L[2];
    }

void AddRelaxation(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    FLOAT  fNorm,ih2,vSigma2,L,fRel,e,pmax,pmin,vMean[3],feps;
    FLOAT beta,gamma,rho;
    int i,j;
    PARTICLE *q;
    double *v;
    double *pRelax;

    pRelax = pkdField(p,pkd->oRelaxation);

    beta = 10.0; /* pmax=beta*l, large beta parameter reduces eps dependence  */
    gamma = 0.17; /* gamma scales overall relaxation time */
    feps = 1.0; /* cuttoff for pmin at feps*epsilon or p_90 deg */

    for (j=0;j<3;++j) {
	vMean[j] = 0.0;
	}
    ih2 = 4.0/BALL2(p);
    vSigma2 = 0.0;
    rho = 0.0;
    fNorm = M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].pPart;
	v = pkdVel(pkd,q);
	rho += pkdMass(pkd,q);
	for (j=0;j<3;++j) {
	    vSigma2 += v[j]*v[j];
	    vMean[j] += v[j];
	    }
	}
    vSigma2 /= nSmooth;
    rho /= 4.18*pow(p->fBall,3.0);
    for (j=0;j<3;++j) {
	vMean[j] /= nSmooth;
	vSigma2 -= vMean[j]*vMean[j];
	}
    vSigma2 = 0.33333*vSigma2; /* now its the one dimensional vel.dispersion squared */
    pmin = 2.0*pkdMass(pkd,p)/(6.0*vSigma2);  /* pmin in comoving units */
    e = feps*pkdSoft(pkd,p);
    if (pmin < e) pmin = e; /* pmin is the bigger one of epsilon and p90 */
    pmax = beta*pow(pkdMass(pkd,p)/rho, 0.333333); /* pmax=beta*mean interp. separation */
    if (pmax < pmin || vSigma2 < 0) return;  /* no relaxation if pmin > pmax */
    L = pmax/pmin;
    fRel = 0.5*( log(1+L*L) - L*L/(1+L*L) );
    fRel *= pkdMass(pkd,p)*rho*pow(vSigma2, -1.5)/gamma;
    *pRelax += fRel * smf->dDeltaT;
    }

#ifdef SYMBA
void DrmininDrift(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    int i,j,k,ial;
    PARTICLE *q;
    double dDelta = smf->pkd->param.dDelta;
    double dSunMass = smf->dSunMass;
    double fMass1 = p->fMass;
    double piOrder = p->iOrder;
    double a1,a2;
    double x0,y0,z0,x1,y1,z1,vx0,vy0,vz0,vx1,vy1,vz1;
    double dx0,dy0,dz0,dvx0,dvy0,dvz0,dvx1,dvy1,dvz1;
    double dr0,dr1,dv0,dv1;
    double a,b,c,d;
    double tmin,drmin,hill;
    /*double a1b, a2b, hillb;*/
    /*
    **  calculates drmin during drift using third order intepolation
    ** dr0, dr1: relative distance before and after drift
    ** dv0, dv1: relative velocity before and after drift
    */

    p->drmin2 = 1000.0;
    x0 = p->rb[0];
    y0 = p->rb[1];
    z0 = p->rb[2];
    vx0 = p->vb[0];
    vy0 = p->vb[1];
    vz0 = p->vb[2];
    x1 = p->r[0];
    y1 = p->r[1];
    z1 = p->r[2];
    vx1 = p->v[0];
    vy1 = p->v[1];
    vz1 = p->v[2];
    a1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    /* a1b = sqrt(x0*x0 + y0*y0 + z0*z0);*/

    for (i=0;i<nSmooth;++i) {
	q = nnList[i].pPart;
	if (piOrder==(q->iOrder)) goto enstep;
	a2 = sqrt(q->r[0]*q->r[0]+q->r[1]*q->r[1]+q->r[2]*q->r[2]);
	/*a2b = sqrt(q->rb[0]*q->rb[0]+q->rb[1]*q->rb[1]+q->rb[2]*q->rb[2]);*/

	hill = 0.5*cbrt((fMass1 + q->fMass)/(3.0*dSunMass));
	/* hillb = hill*(a1b+a2b);*/
	hill *= (a1+a2);

	dr1 = sqrt(nnList[i].dx*nnList[i].dx + nnList[i].dy*nnList[i].dy
		   + nnList[i].dz*nnList[i].dz);

	dx0 = x0 - q->rb[0];
	dy0 = y0 - q->rb[1];
	dz0 = z0 - q->rb[2];
	dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

	if (dr1 < 3.0*hill) {
	    p->drmin2 = dr1/hill;
	    /* printf("dr1 iOrder %d jOrder %d dr0 %e\n",
	    p->iOrder,q->iOrder,dr0/hillb);*/
	    /* check if particle q is already in the list*/
	    ial = 0;
	    if (p->n_VA > 0) {
		for (k=0;k<p->n_VA;k++) {
		    if (p->iOrder_VA[k] == q->iOrder) {
			ial = 1;
			}
		    }
		}

	    if (ial==0) {
		/*if(dr0 <= 3.0*hillb)printf("dr0 %e, drmin%e \n",dr0/hillb,p->drmin);*/
		p->iOrder_VA[p->n_VA] = q->iOrder;
		p->hill_VA[p->n_VA] = hill;
		p->n_VA++;
		assert(piOrder!=(q->iOrder));
		}
	    goto enstep;
	    }

	dvx0 = vx0 - q->vb[0];
	dvy0 = vy0 - q->vb[1];
	dvz0 = vz0 - q->vb[2];
	dvx1 = vx1 - q->v[0];
	dvy1 = vy1 - q->v[1];
	dvz1 = vz1 - q->v[2];

	dv0 = sqrt(dvx0*dx0 + dvy0*dy0 + dvz0*dz0)/dr0;
	dv1 = (dvx1*nnList[i].dx + dvy1*nnList[i].dy + dvz1*nnList[i].dz)/dr1;

	a = 2.0*(dr0-dr1) + (dv0 + dv1)*dDelta;
	b = 3.0*(dr1-dr0) - (2.0*dv0 + dv1)*dDelta;
	c = dv0*dDelta;
	d = 4.0*b*b -12.0*a*c; /* BB-4AC: A=3a, B=2b C=c */

	if (d < 0.0) {
	    /*printf("iOrder %d, 2, drmin2 %e \n",p->iOrder, p->drmin2);*/
	    goto enstep; /* no encounter */
	    }
	else {
	    tmin = 0.5*(-b+sqrt(d))/a;
	    if (tmin < 0.0 || tmin > 1.0) {
		/*printf("iOrder %d, 3, drmin2 %e \n",p->iOrder, p->drmin2);*/
		goto enstep; /* no encounter */
		}
	    else {
		drmin = ((a*tmin + b)*tmin + c)*tmin + dr0;
		if (drmin < 3.0*hill) {
		    p->drmin2 = drmin/hill;
		    /* check if particle q is already in the list*/
		    ial = 0;
		    if (p->n_VA > 0) {
			for (k=0;k<p->n_VA;k++) {
			    if (p->iOrder_VA[k] == q->iOrder) {
				ial = 1;
				}
			    }
			}

		    if (ial==0) {
			/*if(dr0 <= 3.0*hillb)printf("dr0 %e, drmin%e \n",dr0/hillb,p->drmin);*/
			p->iOrder_VA[p->n_VA] = q->iOrder;
			p->hill_VA[p->n_VA] = hill;
			p->n_VA++;
			assert(piOrder!=(q->iOrder));
			}
		    goto enstep;
		    }
		}
	    }
    enstep:
	continue;
	}
    }


#endif
