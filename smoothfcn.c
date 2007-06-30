#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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
		else ak = 6.*0.25/350./3. *(7040-1152/ak+ak*(-10080+ak*(2880+ak*(4200 \
	                 +ak*(-3528+ak*(630+ak*(240+ak*(-90+2*ar2)))))))); \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) adk = 6.*0.25/350./3. * (-2880*2 \
	                 +ar2*(3528*4+ adk*(-1890*5 + adk*(-240*6+ adk*(270*7-6*9*ar2))))); \
		else adk = 6.*0.25/350./3. *((1152/ar2-10080)/adk+(2880*2+adk*(4200*3 \
	                 +adk*(-3528*4+adk*(630*5+adk*(240*6 +adk*(-90*7+2*9*ar2))))))); \
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
		else ak = 0.25*ak*ak*ak; \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
                }
#endif
#endif

void NullSmooth(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    }

void initDensity(void *p)
    {
    ((PARTICLE *)p)->fDensity = 0.0;
    }

void combDensity(void *p1,void *p2)
    {
    ((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
    }

void Density(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    FLOAT ih2,r2,rs,fDensity;
    int i;

    ih2 = 4.0/BALL2(p);
    fDensity = 0.0;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	fDensity += rs*nnList[i].pPart->fMass;
	}
    p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
    }

void DensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    PARTICLE *q;
    FLOAT fNorm,ih2,r2,rs;
    int i;

    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	r2 = nnList[i].fDist2*ih2;
	KERNEL(rs,r2);
	rs *= fNorm;
	q = nnList[i].pPart;
	p->fDensity += rs*q->fMass;
	q->fDensity += rs*p->fMass;
	}
    }

void initParticleMarkDensity(void *p)
    {
    ((PARTICLE *)p)->fDensity = 0.0;
    TYPESet((PARTICLE *) p,TYPE_DensZeroed);
    }

void initMarkDensity(void *p)
    {
    ((PARTICLE *)p)->fDensity = 0.0;
    }

void combMarkDensity(void *p1,void *p2)
    {
    if (TYPETest((PARTICLE *) p1,TYPE_DensZeroed)) 
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
    else if (TYPETest((PARTICLE *) p2,TYPE_DensZeroed)) {
	((PARTICLE *)p1)->fDensity = ((PARTICLE *)p2)->fDensity;
	}
    ((PARTICLE *)p1)->iActive |= ((PARTICLE *)p2)->iActive;
    }

void MarkDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    assert(0);
    }

void MarkDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    PARTICLE *q;
    FLOAT fNorm,ih2,r2,rs;
    int i;
    unsigned int qiActive;

    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    if (TYPETest(p,TYPE_ACTIVE)) {
	TYPESet(p,TYPE_NbrOfACTIVE);
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    q = nnList[i].pPart;
	    p->fDensity += rs*q->fMass;
	    if (TYPETest(q,TYPE_DensZeroed)) 
		q->fDensity += rs*p->fMass;
	    else {
		q->fDensity = rs*p->fMass;
		TYPESet(q, TYPE_DensZeroed);
		}
	    TYPESet(q,TYPE_NbrOfACTIVE);
	    }
	} 
    else {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    q = nnList[i].pPart;
	    if (TYPETest(p,TYPE_DensZeroed)) 
		p->fDensity += rs*q->fMass;
	    else {
		p->fDensity = rs*q->fMass;
		TYPESet(p,TYPE_DensZeroed);
		}
	    if (TYPETest(q,TYPE_DensZeroed)) 
		q->fDensity += rs*p->fMass;
	    else {
		q->fDensity = rs*p->fMass;
		TYPESet(q, TYPE_DensZeroed);
		}
	    qiActive |= q->iActive;
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    }

void initParticleMarkIIDensity(void *p)
    {
    if (TYPEFilter((PARTICLE *) p,TYPE_DensACTIVE|TYPE_DensZeroed,
		   TYPE_DensACTIVE)) {
	((PARTICLE *)p)->fDensity = 0.0;
	TYPESet((PARTICLE *)p,TYPE_DensZeroed);
	}
    }

void initMarkIIDensity(void *p)
    {
    ((PARTICLE *) p)->fDensity = 0.0;
    }

void combMarkIIDensity(void *p1,void *p2)
    {
    if (TYPETest((PARTICLE *) p1,TYPE_DensACTIVE)) {
	if (TYPETest((PARTICLE *) p1,TYPE_DensZeroed)) 
	    ((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	else if (TYPETest((PARTICLE *) p2,TYPE_DensZeroed)) {
	    ((PARTICLE *)p1)->fDensity = ((PARTICLE *)p2)->fDensity;
	    }
	}
    ((PARTICLE *)p1)->iActive |= ((PARTICLE *)p2)->iActive;
    }

void MarkIIDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    assert(0);
    }

void MarkIIDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    PARTICLE *q;
    FLOAT fNorm,ih2,r2,rs;
    int i;
    unsigned int qiActive;

    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    if (TYPETest(p,TYPE_DensACTIVE)) {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    qiActive |= q->iActive;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    p->fDensity += rs*q->fMass;
	    if (TYPETest(q,TYPE_DensACTIVE)) {
		if (TYPETest(q,TYPE_DensZeroed)) 
		    q->fDensity += rs*p->fMass;
		else {
		    q->fDensity = rs*p->fMass;
		    TYPESet(q,TYPE_DensZeroed);
		    }
		}
	    if (TYPETest(p,TYPE_ACTIVE)) TYPESet(q,TYPE_NbrOfACTIVE);
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    else if (TYPETest(p,TYPE_ACTIVE)) {
	TYPESet( p,TYPE_NbrOfACTIVE);
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    TYPESet(q,TYPE_NbrOfACTIVE);
	    if (!TYPETest(q,TYPE_DensACTIVE)) continue;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    if (TYPETest(q,TYPE_DensZeroed)) 
		q->fDensity += rs*p->fMass;
	    else {
		q->fDensity = rs*p->fMass;
		TYPESet(q,TYPE_DensZeroed);
		}
	    }
	}
    else {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    qiActive |= q->iActive;
	    if (!TYPETest(q,TYPE_DensACTIVE)) continue;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    if (TYPETest(q,TYPE_DensZeroed)) 
		q->fDensity += rs*p->fMass;
	    else {
		q->fDensity = rs*p->fMass;
		TYPESet(q,TYPE_DensZeroed);
		}
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    }

void initMark(void *p)
    {
    }

void combMark(void *p1,void *p2)
    {
    ((PARTICLE *)p1)->iActive |= ((PARTICLE *)p2)->iActive;
    }


void initGroupMerge(void *g)
    {
    /* nothing to do here */
    }
void combGroupMerge(void *g1, void *g2)
    {
    FOFGD * gd1 = (FOFGD *)g1;
    FOFGD * gd2 = (FOFGD *)g2;

    /* accept second (remote) group number, if its smaller */	
    if(gd1->iGlobalId < gd2->iGlobalId)       
	gd2->iGlobalId = gd1->iGlobalId;
    else gd1->iGlobalId = gd2->iGlobalId;
    /* if someone says its not my group, accept it */	
    gd1->bMyGroup *= gd2->bMyGroup;
    gd2->bMyGroup *= gd1->bMyGroup;
    }
void initGroupBins(void *b)
    {
    /* nothing to do here */
    }
void combGroupBins(void *b1, void *b2)
    {
    FOFBIN * gb1 = (FOFBIN *)b1;
    FOFBIN * gb2 = (FOFBIN *)b2;

    /* add entries */	
    gb1->nMembers += gb2->nMembers;
    gb1->fMassInBin += gb2->fMassInBin;	
    gb1->fMassEnclosed += gb2->fMassEnclosed;	
    gb1->v2[0] += gb2->v2[0];
    gb1->v2[1] += gb2->v2[1];
    gb1->v2[2] += gb2->v2[2];
    gb1->L[0]  += gb2->L[0];
    gb1->L[1]  += gb2->L[1];
    gb1->L[2]  += gb2->L[2];
    }

#ifdef RELAXATION 
void AddRelaxation(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    FLOAT  fNorm,ih2,vSigma2,L,fRel,e,pmax,pmin,vMean[3],feps;
    FLOAT beta,gamma,rho;
    int i,j;
    PARTICLE *q;

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
	rho += q->fMass;
	for (j=0;j<3;++j) { 
	    vSigma2 += q->v[j]*q->v[j];
	    vMean[j] += q->v[j];
	    } 
	}	
    vSigma2 /= nSmooth;
    rho /= 4.18*pow(p->fBall,3.0);
    for (j=0;j<3;++j) {
	vMean[j] /= nSmooth;
	vSigma2 -= vMean[j]*vMean[j];
   	}
    vSigma2 = 0.33333*vSigma2; /* now its the one dimensional vel.dispersion squared */
    pmin = 2.0*p->fMass/(6.0*vSigma2);  /* pmin in comoving units */
    e = feps*p->fSoft;
    if(pmin < e) pmin = e;  /* pmin is the bigger one of epsilon and p90 */
    pmax = beta*pow(p->fMass/rho, 0.333333); /* pmax=beta*mean interp. separation */
    if(pmax < pmin || vSigma2 < 0) return;   /* no relaxation if pmin > pmax */
    L = pmax/pmin;
    fRel = 0.5*( log(1+L*L) - L*L/(1+L*L) );
    fRel *= p->fMass*rho*pow(vSigma2, -1.5)/gamma;
    p->fRelax += fRel * smf->dDeltaT;	
    }
#endif  /* RELAXATION */

#ifdef SYMBA
void DrmininDrift(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf){
    int i,j;
    PARTICLE *q;
    double dDelta = smf->pkd->param.dDelta;
    double dSunMass = smf->dSunMass;
    double fMass1 = p->fMass;
    double a1,a2;
    double x0,y0,z0,x1,y1,z1,vx0,vy0,vz0,vx1,vy1,vz1;
    double dx0,dy0,dz0,dvx0,dvy0,dvz0,dvx1,dvy1,dvz1;
    double dr0,dr1,dv0,dv1;
    double a,b,c,d;
    double tmin,drmin,hill;
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
  
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].pPart;
	if(p==q)continue; 		
	a2 = sqrt(q->r[0]*q->r[0]+q->r[1]*q->r[1]+q->r[2]*q->r[2]);
	hill = 0.5*(a1+a2)*pow((fMass1 + q->fMass)/(3.0*dSunMass),1.0/3.0);
 
	dr1 = sqrt(nnList[i].dx*nnList[i].dx + nnList[i].dy*nnList[i].dy 
		   + nnList[i].dz*nnList[i].dz);
	
	if(dr1 < 3.0*hill){
	    p->drmin2 = dr1/hill;
	    return;
	}
	    
	dx0 = x0 - q->rb[0];
	dy0 = y0 - q->rb[1];
	dz0 = z0 - q->rb[2];
	dvx0 = vx0 - q->vb[0];
	dvy0 = vy0 - q->vb[1];
	dvz0 = vz0 - q->vb[2];
	dvx1 = vx1 - q->v[0];
	dvy1 = vy1 - q->v[1];
	dvz1 = vz1 - q->v[2];
	
	dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0); 
	dv0 = sqrt(dvx0*dx0 + dvy0*dy0 + dvz0*dz0)/dr0; 
	dv1 = (dvx1*nnList[i].dx + dvy1*nnList[i].dy + dvz1*nnList[i].dz)/dr1; 

	a = 2.0*(dr0-dr1) + (dv0 + dv1)*dDelta;
	b = 3.0*(dr1-dr0) - (2.0*dv0 + dv1)*dDelta; 
	c = dv0*dDelta;
	d = 4.0*b*b -12.0*a*c; /* BB-4AC: A=3a, B=2b C=c */  

	if(d < 0.0){ 
	    return; /* no encounter */
	}else{
	    tmin = 0.5*(-b+sqrt(d))/a;
	    if(tmin < 0.0 || tmin > 1.0){
		return; /* no encounter */
	    }else{
		drmin = ((a*tmin + b)*tmin + c)*tmin + dr0;
		if(drmin < 3.0*hill){
		    p->drmin2 = drmin/hill;		
		    return;
		}
		
	    }
	}   	
    }   
}


#endif 
