#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef PLANETS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "collision.h"

#define BOUNCE_OK 0
#define NEAR_MISS 1


void
pkdNextCollision(PKD pkd,double *dt,int *iOrder1,int *iOrder2)
{
	/*
	 ** Returns time and iOrder of particles for earliest predicted
	 ** collision on this processor. Argument initialization must
	 ** occur in calling routine.
         ** (!!Note temporary we are using iOrgIdx in stead of dtCol)
	 */

	PARTICLE *p;
	int i;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (!p->iColflag) continue; /* skip particles wihout collision flag */
		if (!pkdIsActive(pkd,p)) continue; /* skip over inactive particles */
		if (p->iOrder < 0) continue; /* skip over deleted particles */				
		if (p->dtCol < *dt) {
			*dt = p->dtCol;
			*iOrder1 = p->iOrder;
			*iOrder2 = p->iOrderCol;		
			}
		}
	}

void
pkdGetColliderInfo(PKD pkd,int iOrder,COLLIDER *c)
{
	/*
	 ** Returns collider info for particle with matching iOrder.
	 */

	PARTICLE *p;
	int i,j,k;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (p->iOrder == iOrder) {		 
			c->id.iPid = pkd->idSelf;		
			c->id.iOrder = iOrder;
			c->id.iIndex = i;
			c->id.iOrgIdx = p->iOrgIdx;			
			c->fMass = p->fMass;
			c->fRadius = p->fSoft;		
			for (j=0;j<3;j++) {
				c->r[j] = p->r[j];
				c->v[j] = p->v[j];
				c->w[j] = p->w[j];
#ifdef HERMITE
				c->a[j] = p->a[j];
				c->ad[j] = p->ad[j];
#endif
				}
#ifdef SYMBA
			c->drmin = p->drmin;
			c->n_VA = p->n_VA;
			for(k=0;k<p->n_VA;k++){
			    c->iOrder_VA[k]=p->iOrder_VA[k];
			    c->i_VA[k]=p->i_VA[k];
			    c->hill_VA[k]=p->hill_VA[k];
			    }
#endif
			c->iColor = p->iColor;
			c->dt = p->dt;
			c->iRung = p->iRung;		
			return;
			}
		}
	}

void PutColliderInfo(const COLLIDER *c,int iOrder2,PARTICLE *p,double dt)
{
	/*
	 ** Stores collider info in particle structure (except id & color).
	 ** Also, dt is stored in dtPrevCol for inelastic collapse checks.
	 **
	 ** NOTE: Because colliding particles have their positions traced back to
	 ** the start of the step using their NEW velocities, it is possible for
	 ** a neighbour to no longer lie inside a previous collider's search ball.
	 ** To fix this, the ball radius is expanded by a conservative estimate of
	 ** the displacement amount (using the "Manhattan metric").
	 */

	int i,k;
	double r;

	p->fMass = c->fMass;
	p->fSoft = c->fRadius;
	r = fabs(c->r[0] - p->r[0]) +
		fabs(c->r[1] - p->r[1]) +
		fabs(c->r[2] - p->r[2]);
	for (i=0;i<3;i++) {
	  p->r[i] = c->r[i];
	  p->v[i] = c->v[i];
	  p->w[i] = c->w[i];
#ifdef HERMITE
	  p->a[i] = c->a[i];
	  p->ad[i] = c->ad[i];
#endif
	}
		
#ifdef SYMBA
	p->drmin = c->drmin;
	p->n_VA = c->n_VA;
	for(k=0;k<c->n_VA;k++){
	    p->iOrder_VA[k]=c->iOrder_VA[k];
	    p->i_VA[k]=c->i_VA[k];
	    p->hill_VA[k]=c->hill_VA[k];
	}
#endif
	p->iRung = c->iRung;
	/*p->fBall2 += 2*sqrt(p->fBall2)*r + r*r;
	p->dtPrevCol = dt;
	p->iPrevCol = iOrder2;  stored to avoid false collisions */
	p->dtCol = DBL_MAX; /* nessesary in pkdNextcollision */  
	p->iColflag = 0;  /* here we reset collision flag */
	}

void
pkdMerge(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,double dDensity,
		 COLLIDER **cOut,int *pnOut)
{
	/*
	 ** Merges colliders into new centre-of-mass particle (cOut[0])
	 ** with radius determined by the supplied fixed bulk density.
	 ** If dDensity is zero, the radius of the merged particle is
	 ** instead set by the bulk density of the largest collider.
	 */

	const double dDenFac = 4.0/3*M_PI;
	int iOrderTemp;
	
	COLLIDER *c;
	FLOAT com_pos[3],com_vel[3],rc1[3],rc2[3],vc1[3],vc2[3],ang_mom[3];
	FLOAT m1,m2,m,r1,r2,r,i1,i2,i;
	int k,j;
	FLOAT com_a[3],com_ad[3];

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	if (dDensity)
		r = pow(m/(dDenFac*dDensity),1.0/3);
	else
		r = pow(r1*r1*r1 + r2*r2*r2,1.0/3); /* conserves volume *//*DEBUG used to be r = (m2 > m1 ? r2*pow(m/m2,1.0/3) : r1*pow(m/m1,1.0/3));*/
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	i = 0.4*m*r*r;

	for (k=0;k<3;k++) {
		com_pos[k] = (m1*c1->r[k] + m2*c2->r[k])/m;
		rc1[k] = c1->r[k] - com_pos[k];
		rc2[k] = c2->r[k] - com_pos[k];
		com_vel[k] = (m1*c1->v[k] + m2*c2->v[k])/m;
		vc1[k] = c1->v[k] - com_vel[k];
		vc2[k] = c2->v[k] - com_vel[k];
#ifdef HERMITE
		com_a[k] = (m1*c1->a[k] + m2*c2->a[k])/m;
		com_ad[k] = (m1*c1->ad[k] + m2*c2->ad[k])/m;
#endif
	}

	ang_mom[0] = m1*(rc1[1]*vc1[2] - rc1[2]*vc1[1]) + i1*c1->w[0] +
		         m2*(rc2[1]*vc2[2] - rc2[2]*vc2[1]) + i2*c2->w[0];
	ang_mom[1] = m1*(rc1[2]*vc1[0] - rc1[0]*vc1[2]) + i1*c1->w[1] +
				 m2*(rc2[2]*vc2[0] - rc2[0]*vc2[2]) + i2*c2->w[1];
	ang_mom[2] = m1*(rc1[0]*vc1[1] - rc1[1]*vc1[0]) + i1*c1->w[2] +
				 m2*(rc2[0]*vc2[1] - rc2[1]*vc2[0]) + i2*c2->w[2];

	*pnOut = 1;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut != NULL);

	c = &(*cOut)[0];

	/* Note: id info set in pkdDoCollision(), used only for log */

	c->fMass = m;
	c->fRadius = r;

	for (k=0;k<3;k++) {
		c->r[k] = com_pos[k];
		c->v[k] = com_vel[k];
		c->w[k] = ang_mom[k]/i;
#ifdef HERMITE
		c->a[k] = com_a[k];
		c->ad[k] = com_ad[k];
#endif
		}
#ifdef SYMBA
	c->drmin = c1->drmin;
	c->n_VA = 0;
	
	if(c1->n_VA > 1){
	    /* copy neiboring list of c1 to that of c excluding particle c2*/
	    for(k=0;k<c1->n_VA;k++){
		if(c1->iOrder_VA[k] != c2->id.iOrder){
		    c->iOrder_VA[c->n_VA] = c1->iOrder_VA[k];
		    c->i_VA[c->n_VA] = c1->i_VA[k];
		    c->hill_VA[c->n_VA] = c1->hill_VA[k];
		    /* Since mass has changed, this is not exactly correct!
		       But probably does not cause any problem */ 	 
		    c->n_VA++;
		}
	    }
	}
	
	if(c2->n_VA > 1){
	    for(k=0;k<c2->n_VA;k++){
		iOrderTemp = c2->iOrder_VA[k];
		/* check if neigboring particles for c2 are already included 
		   in the list for c */ 
		if(iOrderTemp == c1->id.iOrder)goto skip_VA;
		for(j=0;j<c->n_VA;j++){
		    if(iOrderTemp == c->iOrder_VA[j])goto skip_VA;
		}			
		c->iOrder_VA[c->n_VA] = iOrderTemp;
		c->i_VA[c->n_VA] = c2->i_VA[k];
		c->hill_VA[c->n_VA] = c2->hill_VA[k];
		c->n_VA++;
	    skip_VA: 
		continue;
	    }
	}
	
#endif
	
	/* Set merger's timestep to iRung of largest mass. */
	/* XXX there is a bug in changing timesteps during a collision 
	   but this makes it less bothersome. */
	/*c->iRung = (c2->fMass > c1->fMass ? c2->iRung : c1->iRung);*/
	c->iRung = (c2->iRung > c1->iRung ? c2->iRung : c1->iRung);
	
}

int
pkdBounce(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
		  double dEpsN,double dEpsT,int bFixCollapse,
		  COLLIDER **cOut,int *pnOut)
{
	/* Bounces colliders, preserving particle order */

	COLLIDER *co1,*co2;
	FLOAT n[3],s1[3],s2[3],v[3],s[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m1,m2,m,r1,r2,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;

/*DEBUG verbose EpsN, EpsT output
	(void) printf("e_n = %g e_t = %g\n",dEpsN,dEpsT);
*/

	*pnOut = 2;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut != NULL);

	(*cOut)[0] = *c1;
	(*cOut)[1] = *c2;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	mu = m1*m2/m;
	alpha = 2.5*(1/m1 + 1/m2);
	beta = 1/(1 + alpha*mu);

	a = 0;
	for (i=0;i<3;i++) {
		n[i] = (c2->r[i] - c1->r[i]);
		a += n[i]*n[i];
		}
	a = 1/sqrt(a);
	for (i=0;i<3;i++)
		n[i] *= a;

	s1[0] = r1*(c1->w[1]*n[2] - c1->w[2]*n[1]);
	s1[1] = r1*(c1->w[2]*n[0] - c1->w[0]*n[2]);
	s1[2] = r1*(c1->w[0]*n[1] - c1->w[1]*n[0]);

	s2[0] = r2*(c2->w[2]*n[1] - c2->w[1]*n[2]);
	s2[1] = r2*(c2->w[0]*n[2] - c2->w[2]*n[0]);
	s2[2] = r2*(c2->w[1]*n[0] - c2->w[0]*n[1]);

	for (i=0;i<3;i++) {
		v[i] = c2->v[i] - c1->v[i];
		s[i] = s2[i] - s1[i];
		}

	for (i=0;i<3;i++)
		u[i] = v[i] + s[i];

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
	if (a >= 0) {
#if (INTERNAL_WARNINGS)
		static int bGiveWarning = 1;
		if (bGiveWarning) {
			(void) fprintf(stderr,"WARNING: %i & %i -- near miss? (a = %g)\n",
						   c1->id.iOrder,c2->id.iOrder,a);
#if (INTERNAL_WARNINGS_ONCE)
			bGiveWarning = 0;
#endif
			}
#endif /* INTERNAL WARNINGS */

		if (bFixCollapse)
			return NEAR_MISS; /* particles remain unchanged */
		else
			assert(0); /* near miss not allowed */
		}
	for (i=0;i<3;i++) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

	a = (1 + dEpsN);
	b = beta*(1 - dEpsT);
	for (i=0;i<3;i++)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

	a =   m2/m;
	b = - m1/m;
	c = r1/i1;
	d = r2/i2;

	co1 = &(*cOut)[0];
	co2 = &(*cOut)[1];

	for (i=0;i<3;i++) {
		co1->v[i] += a*p[i];
		co2->v[i] += b*p[i];
		co1->w[i] += c*q[i];
		co2->w[i] += d*q[i];
		}

/*DEBUG -- dT check
	{
	double dT;
	dT =
		- mu*(v[0]*p[0] + v[1]*p[1] + v[2]*p[2]) +
			0.5*mu*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) +
				(r1*c1->w[0] + r2*c2->w[0])*q[0] +
					(r1*c1->w[1] + r2*c2->w[1])*q[1] +
						(r1*c1->w[2] + r2*c2->w[2])*q[2] +
							0.5*alpha*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
	(void) printf("COLLISION %i & %i e_n=%f e_t=%f dT %e\n",
				  c1->id.iOrder,c2->id.iOrder,dEpsN,dEpsT,dT);
	}
*/

	return BOUNCE_OK;
	}

void
pkdFrag(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
		COLLIDER **cOut,int *pnOut)
{
	/*
	 ** Fragments colliders into at most MAX_NUM_FRAG pieces.
	 ** Not implemented yet. Development notes:
	 ** - be sure new particles are ACTIVE
	 ** - may need to assert(*pnOut <= MAX_NUM_FRAG)
	 ** - remember to set id info for logging purposes
	 ** - and decide on starting iRung values
	 ** - and assign colors
	 */
	}

void
pkdDoCollision(PKD pkd,double dt,const COLLIDER *pc1,const COLLIDER *pc2,
			   int bPeriodic,const COLLISION_PARAMS *CP,int *piOutcome,
			   double *dT,COLLIDER *cOut,int *pnOut)
{
	/*
	 ** Processes collision by advancing particle coordinates to impact
	 ** time, determining the collision outcome, and moving the resulting
	 ** particles back to their (new) start-of-step positions. Diagnostic
	 ** info is returned if the first collider resides on the processor
	 ** (the second collider may not be local, in which case this routine
	 ** is called twice in parallel, once for each particle). Local copies
	 ** of the collider info are used because the particle coordinates need
	 ** to be advanced to the moment of collision. With periodic boundary
	 ** conditions, the second collider may be offset by one period as well.
	 ** The offset is removed before calling PutColliderInfo().
	 */

	COLLIDER c1,c2,*c;
	double v2,ve2,r2,dImpactEnergy;
	int bDiagInfo,iOutcome,i,j,n; /*DEBUG would bReportInfo be better?*/

/*DEBUG verbose collision output
	(void) printf("COLLISION %i & %i (dt = %.16e)\n",
				  pc1->id.iOrder,pc2->id.iOrder,dt);
*/

	/* Get local copies of collider info */

	c1 = *pc1; /* struct copy */
	c2 = *pc2;

	bDiagInfo = (c1.id.iPid == pkd->idSelf);

	if (bDiagInfo && dT) *dT = 0;

	/* Advance coordinates to impact time (no need now)

	for (i=0;i<3;i++) {
		c1.r[i] += c1.v[i]*dt;
		c2.r[i] += c2.v[i]*dt;
	}*/

	/* Determine collision outcome */

	v2 = 0.0;
        r2 = 0.0;
	for (i=0;i<3;i++){
		v2 += (c2.v[i] - c1.v[i])*(c2.v[i] - c1.v[i]);
		r2 += (c2.r[i] - c1.r[i])*(c2.r[i] - c1.r[i]);
	}
	ve2 = 2*(c1.fMass + c2.fMass)/(c1.fRadius + c2.fRadius);

	dImpactEnergy = 0.5*c1.fMass*c2.fMass/(c1.fMass + c2.fMass)*v2;

	iOutcome = MISS;

	if (CP->iOutcomes == MERGE ||
		((CP->iOutcomes & MERGE) &&
		 v2 <= CP->dBounceLimit*CP->dBounceLimit*ve2)) {
		iOutcome = MERGE;
		pkdMerge(pkd,&c1,&c2,CP->dDensity,&c,&n);
		assert(n == 1);
		if (CP->iOutcomes & (BOUNCE | FRAG)) { /* check if spinning too fast */
			double w2max,w2=0;
			w2max = c->fMass/(c->fRadius*c->fRadius*c->fRadius);
			for (i=0;i<3;i++)
				w2 += c->w[i]*c->w[i];
			if (w2 > w2max) {
				int rv;
				free((void *) c);
				iOutcome = BOUNCE;
				rv = pkdBounce(pkd,&c1,&c2,CP->dEpsN,CP->dEpsT,CP->bFixCollapse,&c,&n);
				assert(rv == BOUNCE_OK);
				assert(n == 2);
				}
			}
		}
	else if (CP->iOutcomes & BOUNCE) {
		double dEpsN=1,dEpsT=1;
		iOutcome = BOUNCE;

		        dEpsN = CP->dEpsN;
			dEpsT = CP->dEpsT;
	
		if (pkdBounce(pkd,&c1,&c2,dEpsN,dEpsT,CP->bFixCollapse,&c,&n) == NEAR_MISS)
			iOutcome = MISS;
		assert(n == 2);
		}
	else if (CP->iOutcomes & FRAG) {
		iOutcome = FRAG;
		pkdFrag(pkd,&c1,&c2,&c,&n);
		assert(n <= MAX_NUM_FRAG);
		}

	if (!CP->bFixCollapse)
		assert(iOutcome != MISS); /* SOMETHING has to happen... */

	if (bDiagInfo) *piOutcome = iOutcome;

	assert(n > 0);

	if (bDiagInfo) {
		if (dT) {
			double moi;
			*dT = 0; /* redundant */
			for (j=0;j<n;j++) {
				moi = 0.4*c[j].fMass*c[j].fRadius*c[j].fRadius;
				for (i=0;i<3;i++)
					*dT += c[j].fMass*(c[j].v[i]*c[j].v[i]) +
						moi*(c[j].w[i]*c[j].w[i]);
				}
			moi = 0.4*c1.fMass*c1.fRadius*c1.fRadius;
			for (i=0;i<3;i++)
				*dT -= c1.fMass*(c1.v[i]*c1.v[i]) + moi*(c1.w[i]*c1.w[i]);
			moi = 0.4*c2.fMass*c2.fRadius*c2.fRadius;
			for (i=0;i<3;i++)
				*dT -= c2.fMass*(c2.v[i]*c2.v[i]) + moi*(c2.w[i]*c2.w[i]);
			*dT *= 0.5;
			if(iOutcome = MERGE) /* add potential change*/
                        {
			  *dT += c1.fMass*c2.fMass/sqrt(r2);
                        }
			}
		/*DEBUG (void) printf("Compare: dT = %e\n",*dT); */
		for (i=0;i<n;i++)
			cOut[i] = c[i];
		*pnOut = n;
		}

	/* Trace particles back to start of step (modify this part later)

	for (i=0;i<n;i++)
		for (j=0;j<3;j++)
			c[i].r[j] -= c[i].v[j]*dt; */

	/* Handle output cases */

	if (n == 1) { /* merge */
		PARTICLE_ID *pMrg=NULL,*pDel=NULL,*pOth=NULL;
		if (c1.id.iPid == pkd->idSelf) { /* local particle */
			/*
			 ** Keep this particle if it has larger mass, or, if both
			 ** particles are the same mass, keep this particle if it
			 ** has a lower original index, or, if both particles have
			 ** same original index, keep this particle if it has a
			 ** lower processor number, or, if both particles are on
			 ** same processor, keep particle with lower iOrder.
			 */
			if (c1.fMass > c2.fMass ||
				(c1.fMass == c2.fMass &&
				 (c1.id.iOrgIdx < c2.id.iOrgIdx ||
				  (c1.id.iOrgIdx == c2.id.iOrgIdx &&
				   (c1.id.iPid < c2.id.iPid ||
					(c1.id.iPid == c2.id.iPid &&
					 c1.id.iOrder < c2.id.iOrder))))))
				pMrg = &c1.id; /* new centre-of-mass particle */
			else
				pDel = &c1.id; /* this particle is removed */
			pOth = &c2.id;
			}
		if (c2.id.iPid == pkd->idSelf) { /* same thing for other particle */
			if (c2.fMass > c1.fMass ||
				(c2.fMass == c1.fMass &&
				 (c2.id.iOrgIdx < c1.id.iOrgIdx ||
				  (c2.id.iOrgIdx == c1.id.iOrgIdx &&
				   (c2.id.iPid < c1.id.iPid ||
					(c2.id.iPid == c1.id.iPid &&
					 c2.id.iOrder < c1.id.iOrder))))))
				pMrg = &c2.id;
			else
				pDel = &c2.id;
			pOth = &c1.id;
			}
		/* 
		 ** Store or delete particle as appropriate, and make
		 ** sure master knows ID of surviving particle for
		 ** tracking purposes.
		 */
		if (pMrg) {
			PutColliderInfo(&c[0],INT_MAX,&pkd->pStore[pMrg->iIndex],dt);
			if (bDiagInfo) cOut[0].id = *pMrg; /* struct copy */
			}
		if (pDel) {
			pkdDeleteParticle(pkd,&pkd->pStore[pDel->iIndex]);
			if (bDiagInfo && !pMrg) cOut[0].id = *pOth; /* may need merger info */
			}
		}
	else if (n == 2) { /* bounce or mass transfer */
		if (c1.id.iPid == pkd->idSelf)
			PutColliderInfo(&c[0],c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
		if (c2.id.iPid == pkd->idSelf) {			
			PutColliderInfo(&c[1],c1.id.iOrder,&pkd->pStore[c2.id.iIndex],dt);
			}
		}
	else { /* fragmentation */
		assert(0); /* not implemented yet */
		/* note in this case new iOrders start at pkd->nMaxOrderDark */
		}

	/* Free resources */

	free((void *) c);
	}

void
pkdDoCollisionVeryActive(PKD pkd,double dTime)
{
	COLLIDER c1out,c2out,cOut;
	COLLIDER *c1,*c2,*c;
	double sec;
	unsigned int nCol=0,nMis=0,nMrg=0,nBnc=0,nFrg=0;		 
        int bPeriodic = pkd->param.bPeriodic;
	int iOrder1, iOrder2, piOutcome, pnOut;
	double dt, dT;      
      	
	do{	
		  dt = DBL_MAX;
		  pkdNextCollision(pkd, &dt, &iOrder1, &iOrder2);
		  	
		  /* process the collision */
		  if (COLLISION(dt)) {
		    /*printf("%i,%i\n",iOrder1,iOrder2);*/
		        assert(iOrder1 >= 0);
		        assert(iOrder2 >= 0);		
		    
		        pkdGetColliderInfo(pkd,iOrder1,&c1out);		
			c1 = &c1out; /* struct copy */				
			assert(c1->id.iOrder == iOrder1);
		
			pkdGetColliderInfo(pkd,iOrder2,&c2out);
			c2 = &c2out; /* struct copy */
			assert(c2->id.iOrder == iOrder2);

		  	pkdDoCollision(pkd, dt, c1, c2, bPeriodic,
				       &pkd->param.CP, &piOutcome, &dT, &cOut, &pnOut);
				
			pkd->dDeltaEcoll += dT; /*account for kinetic energy loss + (potential)*/
			++nCol;
			switch (piOutcome) {
			case MISS:
				++nMis;
				--nCol;
				break;
			case MERGE:
				++nMrg;
				break;
			case BOUNCE:
				++nBnc;
				break;
			case FRAG:
				++nFrg;
				break;
			default:
				assert(0); /* unknown outcome */
				}

		switch (pkd->param.iCollLogOption) { /* log collision if requested */
		case COLL_LOG_NONE:
			break;
		case COLL_LOG_VERBOSE:
			{
			FILE *fp;
			int i;

			fp = fopen(pkd->param.achCollLog,"a");
			assert(fp != NULL);

			/* for (i=0;i<3;i++) {
				c1->r[i] += c1->v[i]*dt;
				c2->r[i] += c2->v[i]*dt;
				} */

			fprintf(fp,"%s-%s COLLISION:T=%e\n",
					_pkdParticleLabel(pkd,c1->iColor),
					_pkdParticleLabel(pkd,c2->iColor),dTime);

			fprintf(fp,"***1:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
					c1->id.iPid,c1->id.iOrder,c1->id.iIndex,c1->id.iOrgIdx,
					c1->fMass,c1->fRadius,c1->dt,c1->iRung,
					c1->r[0],c1->r[1],c1->r[2],
					c1->v[0],c1->v[1],c1->v[2],
					c1->w[0],c1->w[1],c1->w[2]);

			fprintf(fp,"***2:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
					c2->id.iPid,c2->id.iOrder,c2->id.iIndex,c2->id.iOrgIdx,
					c2->fMass,c2->fRadius,c2->dt,c2->iRung,
					c2->r[0],c2->r[1],c2->r[2],
					c2->v[0],c2->v[1],c2->v[2],
					c2->w[0],c2->w[1],c2->w[2]);
			fprintf(fp,"***OUTCOME=%s dT=%e\n",
					piOutcome == MISS ? "MISS" :
					piOutcome == MERGE ? "MERGE" :
					piOutcome == BOUNCE ? "BOUNCE" :
					piOutcome == FRAG ? "FRAG" : "UNKNOWN",dT);
			for (i=0;i<(pnOut < MAX_NUM_FRAG ? pnOut : MAX_NUM_FRAG);i++) {
				c = &((&cOut)[i]);

				fprintf(fp,"***out%i:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",i,
					c->id.iPid,c->id.iOrder,c->id.iIndex,c->id.iOrgIdx,
					c->fMass,c->fRadius,c->iRung,
					c->r[0],c->r[1],c->r[2],
					c->v[0],c->v[1],c->v[2],
					c->w[0],c->w[1],c->w[2]);
					}
			fclose(fp);
			break;
			}
		case COLL_LOG_TERSE:
			{
			/*
			 ** FORMAT: For each event, time (double), collider 1 iOrgIdx
			 ** (int), collider 2 iOrgIdx (int), number of post-collision
			 ** particles (int), iOrgIdx for each of these (n * int).
			 */

			FILE *fp;
			XDR xdrs;		       
			int i;

			if (piOutcome != MERGE && piOutcome != FRAG)
			break; /* only care when particle indices change */
			       fp = fopen(pkd->param.achCollLog,"a");
			       assert(fp != NULL);
			       xdrstdio_create(&xdrs,fp,XDR_ENCODE);
			    
				(void) xdr_double(&xdrs,&dTime);
				/* MERGE =1, BOUNCE =2*/
				(void) xdr_int(&xdrs,&piOutcome); 

                                /* info for c1*/
				(void) xdr_int(&xdrs,&c1->iColor);
				(void) xdr_int(&xdrs,&c1->id.iOrgIdx);			
				(void) xdr_double(&xdrs,&c1->fMass);
				(void) xdr_double(&xdrs,&c1->fRadius);
				for (i=0;i<N_DIM;i++)
				 (void)xdr_double(&xdrs,&c1->r[i]);
			        for (i=0;i<N_DIM;i++)
				 (void)xdr_double(&xdrs,&c1->v[i]);
		                for (i=0;i<N_DIM;i++)
				 (void)xdr_double(&xdrs,&c1->w[i]);

			        /* info for c2*/
				(void) xdr_int(&xdrs,&c2->iColor);
				(void) xdr_int(&xdrs,&c2->id.iOrgIdx);			
				(void) xdr_double(&xdrs,&c2->fMass);
				(void) xdr_double(&xdrs,&c2->fRadius);
				for (i=0;i<N_DIM;i++)
				 (void)xdr_double(&xdrs,&c2->r[i]);
			        for (i=0;i<N_DIM;i++)
				 (void)xdr_double(&xdrs,&c2->v[i]);
		                for (i=0;i<N_DIM;i++)
				 (void)xdr_double(&xdrs,&c2->w[i]);
			  				
				xdr_destroy(&xdrs);
				(void) fclose(fp);
				break;
			}
			default:
				assert(0); /* invalid collision log option */
				} /* logging */
			} /* if collision */
		} while (COLLISION(dt));
	
	/* Adddelete is taken care of in master level.*/	

	if (pkd->param.bVStep) {
		double dsec = dTime - sec;
		printf("%i collision%s: %i miss%s, %i merger%s, %i bounce%s, %i frag%s\n",
			   nCol,nCol==1?"":"s",nMis,nMis==1?"":"es",nMrg,nMrg==1?"":"s",
			   nBnc,nBnc==1?"":"s",nFrg,nFrg==1?"":"s");
		printf("Collision search completed, time = %g sec\n",dsec);
		}
	pkd->iCollisionflag = 0;
	}



static char *
_pkdParticleLabel(PKD pkd,int iColor)
{
	switch (iColor) {
	case SUN:
		return "SUN";
	case JUPITER:
		return "JUPITER";
	case SATURN:
		return "SATURN";
	case URANUS:
		return "URANUS";
	case NEPTUNE:
		return "NEPTUNE";
	case PLANETESIMAL:
		return "PLANETESIMAL";
	default:
		return "UNKNOWN";
		}
	}


void pkdGetVariableVeryActive(PKD pkd, double *dDeltaEcoll)
{
  *dDeltaEcoll = pkd->dDeltaEcoll;
  pkd->dDeltaEcoll = 0.0;
}


void pkdCheckHelioDist(PKD pkd,double *dT,double *dSM){
    int i,j,k,n;
    double a2;
    double rsun = 0.1; /* solar radius */
    double resc = 100.0; /* escape distance */
    PARTICLE *p = pkd->pStore;  
    n = pkd->nLocal;
    
    *dT = 0.0;
    *dSM = 0.0;

    for(i=0;i<n;++i) {
	if(p[i].iOrder < 0) continue; 
	a2 = (p[i].r[0]*p[i].r[0] + p[i].r[1]*p[i].r[1] + p[i].r[2]*p[i].r[2]);

	if(a2 < rsun*rsun || a2 > resc*resc){
	    a2 = sqrt(a2);
	    printf("particle %d is deleted with heliocentric distance %e",
		   p[i].iOrder,a2);
	    double moi;
	    /* kinetic and rotational energy */
	    moi = 0.4*p[i].fMass*p[i].fSoft*p[i].fSoft;
	    for (k=0;k<3;k++){
		*dT -=  0.5*(p[i].fMass*(p[i].v[k]*p[i].v[k]) +
		    moi*(p[i].w[k]*p[i].w[k]));
	    }		      	   
	    /* add potential change*/	    
	    *dT += p[i].fMass*pkd->dSunMass/a2;
	    if(a2 < rsun){
	    *dSM += p[i].fMass;
	    }else if (a2 > resc){
		for(j=0;j<n;++j) {
		    if(i==j)continue;
		    a2=0.0;
		    for (k=0;k<3;k++){
			a2 += (p[i].r[k] - p[j].r[k])*(p[i].r[k] - p[j].r[k]);
		    }
		    *dT += p[i].fMass*p[j].fMass/sqrt(a2);  
		}
	    }
	    printf(" dE = %e \n",*dT);
	    pkdDeleteParticle(pkd,&p[i]);
	}       
    }
}


#endif /* PLANETS */
