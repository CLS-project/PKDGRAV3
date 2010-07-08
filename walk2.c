#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef PROFILE_GRAVWALK
#include "VtuneApi.h"
#endif
#ifdef __SSE__
#include <xmmintrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <string.h>
#include "pkd.h"
#include "walk.h"
#include "grav.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "moments.h"
#ifdef TIME_WALK_WORK
#include <sys/time.h>
#endif

const char *walk2_module_id = "$Id$";
const char *walk_h_module_id = WALK_H_MODULE_ID;


#define WALK_MINMULTIPOLE	3


#ifdef TIME_WALK_WORK
typedef struct {
    double fTimer;
    struct timeval tv1;
    } TIMER;

void startTimer(TIMER *t) {
    struct timezone tz;
    gettimeofday(&t->tv1, &tz);
    }

void clearTimer(TIMER *t) {
    t->fTimer = 0.0;
    startTimer(t);
    }

float stopTimer(TIMER *t) {
    struct timezone tz;
    struct timeval tv2;

    gettimeofday(&tv2, &tz);
    tv2.tv_sec -= t->tv1.tv_sec;
    if ( tv2.tv_usec < t->tv1.tv_usec ) {
	tv2.tv_sec--;
	tv2.tv_usec += 1000000;
	}
    tv2.tv_usec -= t->tv1.tv_usec;
    t->fTimer += tv2.tv_sec + (float)tv2.tv_usec/1000000.0;
    return t->fTimer;
    }

float getTimer(TIMER *t) {
    struct timezone tz;
    struct timeval tv2;

    gettimeofday(&tv2, &tz);
    tv2.tv_sec -= t->tv1.tv_sec;
    if ( tv2.tv_usec < t->tv1.tv_usec ) {
	tv2.tv_sec--;
	tv2.tv_usec += 1000000;
	}
    tv2.tv_usec -= t->tv1.tv_usec;
    return t->fTimer + tv2.tv_sec + (float)tv2.tv_usec/1000000.0;
    }
#endif


/*
** This is the new momentum conserving version of the opening criterion.
*/
static inline int iOpenOutcome(PKD pkd,KDN *k,CELT *check,KDN **pc) {
    const double fMonopoleThetaFac = 1.5;
    double dx,dy,dz,minc2,mink2,d2,d2Open,xc,yc,zc,fourh2;
    KDN *c;
    int T1,iCell;
    int iOpen,iOpenA,iOpenB;
        
    if (T1 = (check->iCell < 0)) iCell = -check->iCell;
    else iCell = check->iCell;
    if (check->id == pkd->idSelf) {
	c = pkdTreeNode(pkd,iCell);
    }
    else if (check->id < 0) {
        c = pkdTopNode(pkd,iCell);
	assert(c->iLower != 0);
    }
    else {
	c = mdlAquire(pkd->mdl,CID_CELL,iCell,check->id);
    }

    pBND cbnd = pkdNodeBnd(pkd, c);
    pBND kbnd = pkdNodeBnd(pkd, k);


    *pc = c;
    if (pkdNodeMom(pkd,c)->m <= 0) iOpen = 10;  /* ignore this cell */
    else {
	dx = fabs(k->r[0] - cbnd.fCenter[0] - check->rOffset[0]) - cbnd.fMax[0];
	dy = fabs(k->r[1] - cbnd.fCenter[1] - check->rOffset[1]) - cbnd.fMax[1];
	dz = fabs(k->r[2] - cbnd.fCenter[2] - check->rOffset[2]) - cbnd.fMax[2];
	minc2 = (dx>0)?dx*dx:0 + (dy>0)?dy*dy:0 + (dz>0)?dz*dz:0;
	fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,pkdNodeMom(pkd,c)->m,4*c->fSoft2);
	if (T1) {
	    if (k->iLower == 0) iOpenA = 1; /* add this cell as an "opened" bucket */
	    else iOpenA = 0; /* open this cell and add its children to the checklist */
	    if (minc2 <= k->fOpen*k->fOpen) iOpen = iOpenA;  /* this cell simply stays on the checklist */
	    else if (minc2 > fourh2) iOpen = 6;  /* Add particles of C to local expansion */
	    else if (minc2 > pow(fMonopoleThetaFac*k->fOpen,2)) iOpen = 7;
	    else iOpen = iOpenA;
	} 
	else {
	    if (c->iLower == 0) iOpenA = 2; /* add this cell as an "opened" bucket */
	    else iOpenA = 3; /* open this cell and add its children to the checklist */

	    xc = c->r[0] + check->rOffset[0];
	    yc = c->r[1] + check->rOffset[1];
	    zc = c->r[2] + check->rOffset[2];
	    d2 = pow(k->r[0]-xc,2) + pow(k->r[1]-yc,2) + pow(k->r[2]-zc,2);
	    d2Open = pow(c->fOpen + k->fOpen,2);
	    dx = fabs(xc - kbnd.fCenter[0]) - kbnd.fMax[0];
	    dy = fabs(yc - kbnd.fCenter[1]) - kbnd.fMax[1];
	    dz = fabs(zc - kbnd.fCenter[2]) - kbnd.fMax[2];
	    mink2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);

	    if (c->fOpen > k->fOpen) iOpenB = iOpenA;
	    else if (k->iLower != 0) iOpenB = 0; /* keep this cell on the checklist */
	    else if (mink2 <= c->fOpen*c->fOpen) iOpenB = iOpenA;
	    else if (mink2 > fourh2) iOpenB = 4; /* add this cell to the P-C list */
	    else if (mink2 > pow(fMonopoleThetaFac*c->fOpen,2)) iOpenB = 5; /* use this cell as a softened monopole */
	    else iOpenB = iOpenA;

	    if (d2 <= d2Open) iOpen = iOpenB;
	    else if (minc2 > fourh2 && mink2 > fourh2) iOpen = 8;
	    else if (d2 > fMonopoleThetaFac*fMonopoleThetaFac*d2Open) iOpen = 9; /* it is absolutely critical to include this case for large softening */
	    else iOpen = iOpenB;
	}
    }
    return(iOpen);
}


/*
** This implements the original pkdgrav2m opening criterion, which has been
** well tested, gives good force accuracy, but may not be the most efficient
** and also doesn't explicitly conserve momentum.
*/
static inline int iOpenOutcomeOld(PKD pkd,KDN *k,CELT *check,KDN **pc) {
    const double fMonopoleThetaFac = 1.5;
    const int walk_min_multipole = 3;
    double dx,dy,dz,mink2,d2,d2Open,xc,yc,zc,fourh2;
    KDN *c;
    int iCell,nc;
    int iOpen,iOpenA,iOpenB;
        
    assert(check->iCell > 0);
    iCell = check->iCell;
    if (check->id == pkd->idSelf) {
	c = pkdTreeNode(pkd,iCell);
	nc = c->pUpper - c->pLower + 1;
    }
    else if (check->id < 0) {
        c = pkdTopNode(pkd,iCell);
	assert(c->iLower != 0);
	nc = walk_min_multipole; /* we never allow pp with this cell */
    }
    else {
	c = mdlAquire(pkd->mdl,CID_CELL,iCell,check->id);
	nc = c->pUpper - c->pLower + 1;
    }

    pBND kbnd = pkdNodeBnd(pkd, k);

    *pc = c;
    if (pkdNodeMom(pkd,c)->m <= 0) iOpen = 10;  /* ignore this cell */
    else {
	fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,pkdNodeMom(pkd,c)->m,4*c->fSoft2);
	xc = c->r[0] + check->rOffset[0];
	yc = c->r[1] + check->rOffset[1];
	zc = c->r[2] + check->rOffset[2];
	d2 = pow(k->r[0]-xc,2) + pow(k->r[1]-yc,2) + pow(k->r[2]-zc,2);
	d2Open = pow(2.0*fmax(c->fOpen,k->fOpen),2);
	dx = fabs(xc - kbnd.fCenter[0]) - kbnd.fMax[0];
	dy = fabs(yc - kbnd.fCenter[1]) - kbnd.fMax[1];
	dz = fabs(zc - kbnd.fCenter[2]) - kbnd.fMax[2];
	mink2 = ((dx>0)?dx*dx:0) + ((dy>0)?dy*dy:0) + ((dz>0)?dz*dz:0);
	if (d2 > d2Open && d2 > fourh2) iOpen = 8;
	else {
	    if (c->iLower == 0) iOpenA = 1;
	    else iOpenA = 3;
	    if (nc < walk_min_multipole || mink2 <= c->fOpen*c->fOpen) iOpenB = iOpenA;
	    else if (mink2 > fourh2) return(4);
	    else if (mink2 > pow(fMonopoleThetaFac*c->fOpen,2)) return(5);
	    else iOpenB = iOpenA;
	    if (c->fOpen > k->fOpen) {
		if (c->iLower) iOpen = 3;
		else iOpen = iOpenB;
	    }
	    else {
		if (k->iLower) iOpen = 0;
		else iOpen = iOpenB;
	    }
	}
    }
    return(iOpen);
}

/*
** This implements the original pkdgrav Barnes Hut opening criterion.
*/
static inline int iOpenOutcomeBarnesHut(PKD pkd,KDN *k,CELT *check,KDN **pc) {
  const double fMonopoleThetaFac = 1.5;    /* Note: here we don't use this */
    const int walk_min_multipole = 3;
    double dMin,dMax,min2,max2,d2,fourh2;
    KDN *c;
    int iCell,nc;
    int iOpenA;
    int j;
        
    assert(check->iCell > 0);
    iCell = check->iCell;
    if (check->id == pkd->idSelf) {
	c = pkdTreeNode(pkd,iCell);
	nc = c->pUpper - c->pLower + 1;
    }
    else if (check->id < 0) {
        c = pkdTopNode(pkd,iCell);
	assert(c->iLower != 0);
	nc = walk_min_multipole-1; /* we never allow pp with this cell */
    }
    else {
	c = mdlAquire(pkd->mdl,CID_CELL,iCell,check->id);
	nc = c->pUpper - c->pLower + 1;
    }

    pBND kbnd = pkdNodeBnd(pkd, k);

    *pc = c;

    if (c->iLower == 0) iOpenA = 1;
    else iOpenA = 3;

    if (pkdNodeMom(pkd,c)->m <= 0) return(10);  /* ignore this cell */
    else if (k->iLower) {
	/*
	** If this cell is not a bucket calculate the min distance
	** from the center of mass of the check cell to the bounding
	** box of the cell we are calculating gravity on. We also
	** calculate the size of the ball (from checkcell CM) which
	** just contains the this cell (radius^2 given by max2).
	*/
	min2 = 0;
	max2 = 0;
	for (j=0;j<3;++j) {
	  dMin = fabs(c->r[j] + check->rOffset[j] - kbnd.fCenter[j]);
	  dMax = dMin + kbnd.fMax[j];
	  dMin -= kbnd.fMax[j];
	  if (dMin > 0) min2 += dMin*dMin;
	  max2 += dMax*dMax;
	}
	if (max2 <= c->fOpen*c->fOpen) return(iOpenA); /* open it for all particles of c */
	else if (min2 > c->fOpen*c->fOpen) {
	  /*
	  ** For symmetrized softening we need to calculate the distance between
	  ** center of mass between the two cells.
	  */
	  d2 = 0;
	  for (j=0;j<3;++j) {
	    d2 += pow(c->r[j] + check->rOffset[j] - k->r[j],2);
	  }
	  fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,pkdNodeMom(pkd,c)->m,4*c->fSoft2);
	  if (d2 > fourh2) {
	    if (nc >= walk_min_multipole) return(4);  /* accept multipole */
 	    else return(iOpenA);  /* open the cell for performance reasons */
	  }
	  else return(0);   /* in all other cases we can't decide until we get down to a bucket */
	}
	else return(0);
    }
    else {
      /*
      ** If this cell is a bucket we have to either open the checkcell
      ** and get the particles, or accept the multipole. For this
      ** reason we only need to calculate min2.
      */
      min2 = 0;
      for (j=0;j<3;++j) {
	dMin = fabs(c->r[j] + check->rOffset[j] - kbnd.fCenter[j]);
	dMin -= kbnd.fMax[j];
	if (dMin > 0) min2 += dMin*dMin;
      }
      /*
      ** By default we open the cell!
      */
      if (min2 > c->fOpen*c->fOpen) {
	/*
	** For symmetrized softening we need to calculate the distance between
	** center of mass between the two cells.
	*/
	d2 = 0;
	for (j=0;j<3;++j) {
	    d2 += pow(c->r[j] + check->rOffset[j] - k->r[j],2);
	}
	fourh2 = softmassweight(pkdNodeMom(pkd,k)->m,4*k->fSoft2,pkdNodeMom(pkd,c)->m,4*c->fSoft2);
	if (d2 > fourh2) {
	  if (nc >= walk_min_multipole) return(4);  /* accept multipole */
	  else return(iOpenA);  /* open the cell for performance reasons */
	}
	else return(5); /* means we treat this cell as a softened monopole */
      }
      else return(iOpenA);  /* open the cell */
    }
}

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,
		int bVeryActive,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p;
    KDN *k,*c,*kFind,*kSib;
    MOMR *momc,*momk;
    LOCR L;
    double dirLsum,normLsum,adotai,maga;
    double tax,tay,taz;
    double fWeight = 0.0;
    double dShiftFlop;
    double dRhoFac;
    double *v, *a, zero[3];
    FLOAT d2,fourh2;
    FLOAT fMass,fSoft;
    FLOAT rOffset[3];
    FLOAT xParent,yParent,zParent;
    FLOAT dx[3],dir,dir2;
    FLOAT cx,cy,cz,d2c;
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxInitCheck,nCheck;
    int iCell,iSib,iCheckCell,iCellDescend;
    int i,ii,j,pj,nActive,nTotActive;
    int iOpen;
    ILPTILE tile;
    ILCTILE ctile;
#ifdef USE_SIMD_MOMR
    int ig,iv;
#endif

    double tempI;

    double dEwFlop = 0.0;

#ifdef PROFILE_GRAVWALK
    VTResume();
#endif

    for (j=0;j<3;++j) zero[j] = 0.0;
    a = &zero[0];
    v = &zero[0];
    assert(pkd->oNodeMom);
    if (pkd->param.bGravStep) {
	assert(pkd->oNodeAcceleration);
	if (pkd->param.iTimeStepCrit == 1) {
	    assert(pkd->oNodeVelocity);
	    assert(pkd->oVelocity);
	}
    }

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    ** Otherwise we check that the ROOT has active particles!
    */
    if (bVeryActive) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkdTreeNode(pkd,VAROOT)->pUpper - pkdTreeNode(pkd,VAROOT)->pLower + 1);
	}
    else if (!pkdIsCellActive(pkdTreeNode(pkd,ROOT),uRungLo,uRungHi)) return 0;
    /*
    ** Initially we set our cell pointer to
    ** point to the top tree.
    */
    nTotActive = 0;
    ilpClear(pkd->ilp);
    ilcClear(pkd->ilc);

    /*
    ** Allocate Checklist.
    */
    nMaxInitCheck = 2*nReps+1;
    nMaxInitCheck = nMaxInitCheck*nMaxInitCheck*nMaxInitCheck;	/* all replicas */
    iCell = pkd->iTopRoot;
    while ((iCell = pkdTopNode(pkd,iCell)->iParent)) ++nMaxInitCheck; /* all top tree siblings */
    assert(nMaxInitCheck < pkd->nMaxCheck);  /* we should definitely have enough to cover us here! */
    nCheck = 0;
    iStack = -1;
    /*
    ** Clear local expansion and the timestepping sums.
    */
    momClearLocr(&L);
    dirLsum = 0;
    normLsum = 0;
    /*
    ** Precalculate RhoFac if required.
    */
    if (pkd->param.bGravStep) {
	double a = csmTime2Exp(pkd->param.csm,dTime);
	dRhoFac = 1.0/(a*a*a);
	}
    else dRhoFac = 0.0;
    /*
    ** First we add any replicas of the entire box
    ** to the Checklist.
    */
    for (ix=-nReps;ix<=nReps;++ix) {
	rOffset[0] = ix*pkd->fPeriod[0];
	for (iy=-nReps;iy<=nReps;++iy) {
	    rOffset[1] = iy*pkd->fPeriod[1];
	    for (iz=-nReps;iz<=nReps;++iz) {
		rOffset[2] = iz*pkd->fPeriod[2];
		bRep = ix || iy || iz;
		if (bRep || bVeryActive) {
		    pkd->Check[nCheck].iCell = ROOT;
		    /* If leaf of top tree, use root of
		       local tree.
		    */
		    if (pkdTopNode(pkd,ROOT)->iLower) {
			pkd->Check[nCheck].id = -1;
			}
		    else {
			pkd->Check[nCheck].id = pkdTopNode(pkd,ROOT)->pLower;
			}
		    for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = rOffset[j];
		    ++nCheck;
		    }
		if (bRep && bVeryActive) {
		    /*
		    ** Add the images of the very active tree to the checklist.
		    */
		    pkd->Check[nCheck].iCell = VAROOT;
		    pkd->Check[nCheck].id = mdlSelf(pkd->mdl);
		    for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = rOffset[j];
		    ++nCheck;
		    }
		}
	    }
	}
    if (!bVeryActive) {
	/*
	** This adds all siblings of a chain leading from the local tree leaf in the top
	** tree up to the ROOT of the top tree.
	*/
	iCell = pkd->iTopRoot;
	iSib = SIBLING(iCell);
	while (iSib) {
	    if (pkdTopNode(pkd,iSib)->iLower) {
		pkd->Check[nCheck].iCell = iSib;
		pkd->Check[nCheck].id = -1;
		}
	    else {
		/* If leaf of top tree, use root of local tree */
		pkd->Check[nCheck].iCell = ROOT;
		pkd->Check[nCheck].id = pkdTopNode(pkd,iSib)->pLower;
		}
	    for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
	    ++nCheck;
	    iCell = pkdTopNode(pkd,iCell)->iParent;
	    iSib = SIBLING(iCell);
	    }
	}
    /*
    ** Initialize the PP interaction list center, just in case it has not been done yet!
    */
    pkd->ilp->cx = 0;
    pkd->ilp->cy = 0;
    pkd->ilp->cz = 0;
    pkd->ilp->d2cmax = 0;

    /*
    ** We are now going to work on the local tree.
    ** Make iCell point to the root of the tree again.
    */
    if (bVeryActive) k = pkdTreeNode(pkd,iCell = VAROOT);
    else k = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
	/*
	** Find the next active particle that will be encountered in the walk algorithm below
	** in order to set a good initial center for the P-P interaction list.
	*/
	kFind = k;
	while (kFind->iLower) {
	    kFind = pkdTreeNode(pkd,iCellDescend = kFind->iLower);
	    if (!kFind->nActive) {
		/*
		** Move onto processing the sibling.
		*/
		kFind = pkdTreeNode(pkd,++iCellDescend);
	    }
	}
	for (pj=kFind->pLower;pj<=kFind->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    if (!pkdIsActive(pkd,p)) continue;
	    cx = p->r[0];
	    cy = p->r[1];
	    cz = p->r[2];
	    goto found_it;
	}
	assert(0); /* We didn't find an active particle */
    found_it:
	d2c = (cx - pkd->ilp->cx)*(cx - pkd->ilp->cx) + (cy - pkd->ilp->cy)*(cy - pkd->ilp->cy) +
	      (cz - pkd->ilp->cz)*(cz - pkd->ilp->cz);
	/*if (d2c > pkd->ilp->d2cmax) {*/
	if ( d2c > 1e-5) {
/*	    printf("%d:Shift of center too large for the coming interactions! old:(%.10g,%.10g,%.10g) new:(%.10g,%.10g,%.10g)\n",
  mdlSelf(pkd->mdl),pkd->ilp->cx,pkd->ilp->cy,pkd->ilp->cz,cx,cy,cz); */
	    /*
	    ** Correct all remaining PP interactions to this new center.
	    */
	    ILP_LOOP( pkd->ilp, tile ) {
		for ( j=0; j<tile->nPart; ++j ) {
		    tile->d.dx.f[j] += cx - pkd->ilp->cx;
		    tile->d.dy.f[j] += cy - pkd->ilp->cy;
		    tile->d.dz.f[j] += cz - pkd->ilp->cz;
		    }
		}
	    pkd->ilp->cx = cx;
	    pkd->ilp->cy = cy;
	    pkd->ilp->cz = cz;
	    pkd->ilp->d2cmax = k->fOpen*k->fOpen;
	    }
	while (1) {
	    /*
	    ** Process the Checklist.
	    */
	    tempI = *pdFlop;
	    tempI += dEwFlop;
	    ii = 0;
	    if (pkd->param.bGravStep) {
		a = pkdNodeAccel(pkd,k);
		maga = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	    }
/*
	    printf("\nCELL:%d ",iCell);
*/
	    for (i=0;i<nCheck;++i) {
		iOpen = iOpenOutcomeOld(pkd,k,&pkd->Check[i],&c);
/*
		printf("%1d",iOpen);
*/
		switch (iOpen) {
		case 0:
		    /*
		    ** This checkcell stays on the checklist.
		    */
		    pkd->Check[ii++] = pkd->Check[i];
		    break;
		case 1:
		    /*
		    ** This checkcell's particles are added to the P-P list.
		    */
		    for (pj=c->pLower;pj<=c->pUpper;++pj) {
			if (pkd->Check[i].id == pkd->idSelf) p = pkdParticle(pkd,pj);
			else p = mdlAquire(pkd->mdl,CID_PARTICLE,pj,pkd->Check[i].id);
			fMass = pkdMass(pkd,p);
			fSoft = pkdSoft(pkd,p);
			if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdVel(pkd,p);
			ilpAppend(pkd->ilp,
				  p->r[0] + pkd->Check[i].rOffset[0],
				  p->r[1] + pkd->Check[i].rOffset[1],
				  p->r[2] + pkd->Check[i].rOffset[2],
				  fMass, 4*fSoft*fSoft,
				  p->iOrder, v[0], v[1], v[2]);
			if (pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
		    }
		    break;
		case 2:
		    /*
		    ** Now I am trying to open a bucket, which means I keep this cell as an "opened bucket"
		    ** on the checklist. This is marked by a negative cell id. Could prefetch all the 
		    ** bucket's particles at this point if we want.
		    */
		    if (nCheck + 1 > pkd->nMaxCheck) {
			pkd->nMaxCheck += 1000;
			pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			assert(pkd->Check != NULL);
			for (ism=0;ism<pkd->nMaxStack;++ism) {
			    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->S[ism].Check != NULL);
			}
			printf("Case 2: CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		    }
		    pkd->Check[nCheck] = pkd->Check[i];
		    pkd->Check[nCheck].iCell = -pkd->Check[nCheck].iCell;
		    nCheck += 1;
		    break;
		case 3:
		    /*
		    ** Open the cell.
		    ** Here we ASSUME that the children of
		    ** c are all in sequential memory order!
		    ** (the new tree build assures this)
		    ** (also true for the top tree)
		    ** We could do a prefetch here for non-local
		    ** cells.
		    */
		    if (nCheck + 2 > pkd->nMaxCheck) {
			pkd->nMaxCheck += 1000;
			pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			assert(pkd->Check != NULL);
			for (ism=0;ism<pkd->nMaxStack;++ism) {
			    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->S[ism].Check != NULL);
			}
			printf("Case 3: CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		    }
		    iCheckCell = c->iLower;
		    pkd->Check[nCheck] = pkd->Check[i];
		    pkd->Check[nCheck+1] = pkd->Check[i];
		    /*
		    ** If we are opening a leaf of the top tree
		    ** we need to correctly set the processor id.
		    ** (this is a bit tricky)
		    */
		    if (pkd->Check[i].id < 0) {
			if (pkdTopNode(pkd,iCheckCell)->pLower >= 0) {
			    pkd->Check[nCheck].iCell = ROOT;
			    pkd->Check[nCheck].id = pkdTopNode(pkd,iCheckCell)->pLower;
			}
			else {
			    pkd->Check[nCheck].iCell = iCheckCell;
			    assert(pkd->Check[nCheck].id == -1);
			}
			if (pkdTopNode(pkd,iCheckCell+1)->pLower >= 0) {
			    pkd->Check[nCheck+1].iCell = ROOT;
			    pkd->Check[nCheck+1].id = pkdTopNode(pkd,iCheckCell+1)->pLower;
			}
			else {
			    pkd->Check[nCheck+1].iCell = iCheckCell+1;
			    assert(pkd->Check[nCheck+1].id == -1);
			}
		    }
		    else {
			pkd->Check[nCheck].iCell = iCheckCell;
			pkd->Check[nCheck+1].iCell = iCheckCell+1;
			assert(pkd->Check[nCheck].id == pkd->Check[i].id);
			assert(pkd->Check[nCheck+1].id == pkd->Check[i].id);
		    }
		    nCheck += 2;
		    break;
		case 4:
		    /*
		    ** Accept multipole!
		    ** Interact += Moment(c);
		    */
		    ctile = pkd->ilc->tile;
		    if (ctile->nCell == ctile->nMaxCell) {
			ctile = ilcExtend(pkd->ilc);
		    }
		    j = ctile->nCell;
		    momc = pkdNodeMom(pkd,c);
		    ctile->d[j].x.f = c->r[0] + pkd->Check[i].rOffset[0];
		    ctile->d[j].y.f = c->r[1] + pkd->Check[i].rOffset[1];
		    ctile->d[j].z.f = c->r[2] + pkd->Check[i].rOffset[2];
		    ctile->d[j].m.f = momc->m;
		    ctile->d[j].xx.f = momc->xx;
		    ctile->d[j].yy.f = momc->yy;
		    ctile->d[j].xy.f = momc->xy;
		    ctile->d[j].xz.f = momc->xz;
		    ctile->d[j].yz.f = momc->yz;
		    ctile->d[j].xxx.f = momc->xxx;
		    ctile->d[j].xyy.f = momc->xyy;
		    ctile->d[j].xxy.f = momc->xxy;
		    ctile->d[j].yyy.f = momc->yyy;
		    ctile->d[j].xxz.f = momc->xxz;
		    ctile->d[j].yyz.f = momc->yyz;
		    ctile->d[j].xyz.f = momc->xyz;
		    ctile->d[j].xxxx.f = momc->xxxx;
		    ctile->d[j].xyyy.f = momc->xyyy;
		    ctile->d[j].xxxy.f = momc->xxxy;
		    ctile->d[j].yyyy.f = momc->yyyy;
		    ctile->d[j].xxxz.f = momc->xxxz;
		    ctile->d[j].yyyz.f = momc->yyyz;
		    ctile->d[j].xxyy.f = momc->xxyy;
		    ctile->d[j].xxyz.f = momc->xxyz;
		    ctile->d[j].xyyz.f = momc->xyyz;
		    ++ctile->nCell;
		    break;
		case 5:
		    /*
		    ** We accept this multipole from the opening criterion, but it is a softened
		    ** interaction, so we need to treat is as a softened monopole by putting it
		    ** on the particle interaction list.
		    */
		    if (pkd->param.bGravStep && pkd->param.iTimeStepCrit == 1) v = pkdNodeVel(pkd,c);
		    ilpAppend(pkd->ilp,
			      c->r[0] + pkd->Check[i].rOffset[0],
			      c->r[1] + pkd->Check[i].rOffset[1],
			      c->r[2] + pkd->Check[i].rOffset[2],
			      pkdNodeMom(pkd,c)->m, 4*c->fSoft2,
			      -1, /* set iOrder to negative value for time step criterion */
			      v[0], v[1], v[2]);
		    break;
		case 6:
		    /*
		    ** This is accepting an "opened" bucket's particles as monopoles for the
		    ** local expansion.
		    */
		    for (pj=c->pLower;pj<=c->pUpper;++pj) {
			/*
			** NOTE: We can only use iClass here if the class table was
			** merged during a parallel read with GetClasses/SetClasses.
			** There is no problem id a serial read was performed.
			*/
			if (pkd->Check[i].id == pkd->idSelf) p = pkdParticle(pkd,pj);
			else p = mdlAquire(pkd->mdl,CID_PARTICLE,pj,pkd->Check[i].id);
			if (pkdMass(pkd,p) <= 0.0) continue;
			/*
			** Monopole Local expansion accepted!
			*/
			d2 = 0;
			for (j=0;j<3;++j) {
			    dx[j] = k->r[j] - (p->r[j] + pkd->Check[i].rOffset[j]);
			    d2 += dx[j]*dx[j];
			}
			dir = 1.0/sqrt(d2);
			*pdFlop += momLocrAddMono5(&L,pkdMass(pkd,p),dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
			adotai = a[0]*tax + a[1]*tay + a[2]*taz;
			if (adotai > 0) {
			    adotai /= maga;
			    dirLsum += dir*adotai*adotai;
			    normLsum += adotai*adotai;
			}
			if (pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
		    }
		    break;
		case 7:
		    /*
		    ** This is accepting an "opened" bucket's particles with this (k) cell as a softened monopole.
		    ** This is the inverse of accepting a cell as a softened monopole, here we calculate the first 
		    ** order local expansion of each softened particle of the checkcell.
		    */
		    for (pj=c->pLower;pj<=c->pUpper;++pj) {
			/*
			** NOTE: We can only use iClass here if the class table was
			** merged during a parallel read with GetClasses/SetClasses.
			** There is no problem id a serial read was performed.
			*/
			if (pkd->Check[i].id == pkd->idSelf) p = pkdParticle(pkd,pj);
			else p = mdlAquire(pkd->mdl,CID_PARTICLE,pj,pkd->Check[i].id);
			fMass = pkdMass(pkd,p);
			if (fMass <= 0.0) continue;
			fSoft = pkdSoft(pkd,p);
			momk = pkdNodeMom(pkd,k);
			/*
			** Monopole Local expansion accepted!
			*/
			d2 = 0;
			for (j=0;j<3;++j) {
			    dx[j] = k->r[j] - (p->r[j] + pkd->Check[i].rOffset[j]);
			    d2 += dx[j]*dx[j];
			}
			dir = 1.0/sqrt(d2);
			fourh2 = softmassweight(fMass,4*fSoft*fSoft,momk->m,4*k->fSoft2);
			if (d2 > fourh2) {
			    dir = 1.0/sqrt(d2);
			    dir2 = dir*dir*dir;
			}
			else {
			    /*
			    ** This uses the Dehnen K1 kernel function now, it's fast!
			    */
			    dir = 1.0/sqrt(fourh2);
			    dir2 = dir*dir;
			    d2 *= dir2;
			    dir2 *= dir;
			    d2 = 1 - d2;
			    dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
			    dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
			}
			dir2 *= fMass;
			tax = -dx[0]*dir2;
			tay = -dx[1]*dir2;
			taz = -dx[2]*dir2;
			L.m -= fMass*dir;
			L.x -= tax;
			L.y -= tay;
			L.z -= taz;
			adotai = a[0]*tax + a[1]*tay + a[2]*taz;
			if (adotai > 0) {
			    adotai /= maga;
			    dirLsum += dir*adotai*adotai;
			    normLsum += adotai*adotai;
			}
			if (pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
		    }
		    break;
		case 8:
		    /*
		    ** Local expansion accepted!
		    */
		    d2 = 0;
		    for (j=0;j<3;++j) {
			dx[j] = k->r[j] - (c->r[j] + pkd->Check[i].rOffset[j]);
			d2 += dx[j]*dx[j];
		    }
		    dir = 1.0/sqrt(d2);
		    *pdFlop += momLocrAddMomr5(&L,pkdNodeMom(pkd,c),dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
		    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
		    if (adotai > 0) {
			adotai /= maga;
			dirLsum += dir*adotai*adotai;
			normLsum += adotai*adotai;
		    }
		    break;
		case 9:
		    /*
		    ** Here we compute the local expansion due to a single monopole term which could be softened.
		    */
		    momk = pkdNodeMom(pkd,k);
		    momc = pkdNodeMom(pkd,c);
		    d2 = 0;
		    for (j=0;j<3;++j) {
			dx[j] = k->r[j] - (c->r[j] + pkd->Check[i].rOffset[j]);
			d2 += dx[j]*dx[j];
		    }
		    dir = 1.0/sqrt(d2);
		    fourh2 = softmassweight(momc->m,4*c->fSoft2,momk->m,4*k->fSoft2);
		    if (d2 > fourh2) {
			dir = 1.0/sqrt(d2);
			dir2 = dir*dir*dir;
		    }
		    else {
			/*
			** This uses the Dehnen K1 kernel function now, it's fast!
			*/
			dir = 1.0/sqrt(fourh2);
			dir2 = dir*dir;
			d2 *= dir2;
			dir2 *= dir;
			d2 = 1 - d2;
			dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
			dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
		    }
		    dir2 *= momc->m;
		    tax = -dx[0]*dir2;
		    tay = -dx[1]*dir2;
		    taz = -dx[2]*dir2;
		    L.m -= momc->m*dir;
		    L.x -= tax;
		    L.y -= tay;
		    L.z -= taz;
		    adotai = a[0]*tax + a[1]*tay + a[2]*taz;
		    if (adotai > 0) {
			adotai /= maga;
			dirLsum += dir*adotai*adotai;
			normLsum += adotai*adotai;
		    }
		    break;
		case 10:
		    /*
		    ** This checkcell is removed from the checklist since it has zero/negative mass.
		    */
		    break;		
		}
		if (pkd->Check[i].id >= 0 && pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
	    }
	    nCheck = ii;
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (!k->iLower) break;
	    xParent = k->r[0];
	    yParent = k->r[1];
	    zParent = k->r[2];
	    k = pkdTreeNode(pkd,iCell = k->iLower);
	    /*
	    ** Make sure all the check lists are long enough to handle 1 more cell.
	    */
	    if (nCheck == pkd->nMaxCheck) {
		pkd->nMaxCheck += 1000;
		pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
		assert(pkd->Check != NULL);
		for (ism=0;ism<pkd->nMaxStack;++ism) {
		    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
		    assert(pkd->S[ism].Check != NULL);
		    }
		printf("F CPU:%d increased check list size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		}
	    /*
	    ** Check iCell is active. We eventually want to just to a
	    ** rung check here when we start using tree repair, but
	    ** for now this is just as good.
	    */
	    if (k->nActive) {
		/*
		** iCell is active, continue processing it.
		** Put the sibling onto the checklist.
		*/
		pkd->Check[nCheck].iCell = iCell + 1;
		pkd->Check[nCheck].id = pkd->idSelf;
		for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
		++nCheck;
		/*
		** Test whether the sibling is active as well.
		** If not we don't push it onto the stack, but we
		** have to be careful to not pop the stack when we
		** hit the sibling. See the goto InactiveAscend below
		** for how this is done.
		*/
		kSib = pkdTreeNode(pkd,iCell+1);
		if (kSib->nActive) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    ilpCheckPt(pkd->ilp,&pkd->S[iStack].PartChkPt);
		    ilcCheckPt(pkd->ilc,&pkd->S[iStack].CellChkPt);
		    pkd->S[iStack].nCheck = nCheck;
		    /*
		    ** Maybe use a memcpy here!
		    ** for (i=0;i<nCheck;++i) pkd->S[iStack].Check[i] = pkd->Check[i];
		    */
		    memcpy(pkd->S[iStack].Check,pkd->Check,nCheck*sizeof(CELT));
		    pkd->S[iStack].Check[nCheck-1].iCell = iCell;
		    pkd->S[iStack].L = L;
		    pkd->S[iStack].dirLsum = dirLsum;
		    pkd->S[iStack].normLsum = normLsum;
		    dShiftFlop = momShiftLocr(&pkd->S[iStack].L,
					      kSib->r[0] - xParent,
					      kSib->r[1] - yParent,
					      kSib->r[2] - zParent);
		    pkd->S[iStack].fWeight = (*pdFlop-tempI) + dShiftFlop;
		    pkd->S[iStack].fWeight += dEwFlop;
		    }
		}
	    else {
		/*
		** Skip iCell, but add it to the Checklist.
		** No need to push anything onto the stack here.
		*/
		pkd->Check[nCheck].iCell = iCell;
		pkd->Check[nCheck].id = pkd->idSelf;
		for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
		++nCheck;
		/*
		** Move onto processing the sibling.
		*/
		k = pkdTreeNode(pkd,++iCell);
		}
	    *pdFlop += momShiftLocr(&L,k->r[0] - xParent,
				    k->r[1] - yParent,
				    k->r[2] - zParent);
	    }
	/*
	** Now the interaction list should be complete and the
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
	assert(nCheck == 0);

	/*
	** Now calculate gravity on this bucket!
	*/
	nActive = pkdGravInteract(pkd,uRungLo,uRungHi,k,&L,pkd->ilp,pkd->ilc,
				  dirLsum,normLsum,bEwald,pdFlop,&dEwFlop,dRhoFac);
	/*
	** Update the limit for a shift of the center here based on the opening radius of this
	** cell (the one we just evaluated).
	*/
	pkd->ilp->d2cmax = k->fOpen*k->fOpen;
	fWeight += (*pdFlop-tempI);
	fWeight += dEwFlop;
	if (nActive) {
	    fWeight /= nActive;
	    /*
	    ** Here we used to set the weights of particles based on the work done, but now we just assume that
	    ** all particles cost the same in domain decomposition, so we really don't need to set anything here.
	    */
	    *pdPartSum += nActive*ilpCount(pkd->ilp);
	    *pdCellSum += nActive*ilcCount(pkd->ilc);
	    nTotActive += nActive;
	    }

	while (iCell & 1) {
	InactiveAscend:
	    k = pkdTreeNode(pkd,iCell = k->iParent);
	    if (!iCell) {
		/*
		** Make sure stack is empty.
		*/
		assert(iStack == -1);
#ifdef PROFILE_GRAVWALK
		VTPause();
#endif
		*pdFlop += dEwFlop;   /* Finally add the ewald score to get a proper float count */

		return(nTotActive);
		}
	    }
	k = pkdTreeNode(pkd,++iCell);
	if (!k->nActive) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	ilpRestore(pkd->ilp,&pkd->S[iStack].PartChkPt);
	ilcRestore(pkd->ilc,&pkd->S[iStack].CellChkPt);
	nCheck = pkd->S[iStack].nCheck;
	/*
	** Use a memcpy here. This is where we would win with lists since we could simply take the
	** pointer to this list. We would have to link the old checklist into the freelist.
	** for (i=0;i<nCheck;++i) Check[i] = S[iStack].Check[i];
	*/
	memcpy(pkd->Check,pkd->S[iStack].Check,nCheck*sizeof(CELT));
	L = pkd->S[iStack].L;
	dirLsum = pkd->S[iStack].dirLsum;
	normLsum = pkd->S[iStack].normLsum;
	fWeight = pkd->S[iStack].fWeight;
	tempI = *pdFlop;
	tempI += dEwFlop;
	--iStack;
	}
    printf("\n");
    }

