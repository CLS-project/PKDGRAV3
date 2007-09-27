#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "walk.h"
#include "grav.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "ewald.h"

typedef struct CheckElt {
    int iCell;
    int id;
    FLOAT rOffset[3];
    } CELT;

typedef struct CheckStack {
    int nPart;
    int nCell;
    int nPartBucket;
    int nCheck;
    CELT *Check;
    } CSTACK;

#define WALK_MINMULTIPOLE	3

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,double dTime,int nReps,int bEwald,int bEwaldKick,int bVeryActive,double fEwCut,
		double *pdFlop,double *pdPartSum,double *pdCellSum)
    {
    PARTICLE *p = pkd->pStore;
    PARTICLE *pRemote;
    KDN *c;
    KDN *pkdc;
    CSTACK *S;
    CELT *Check;
    ILP *ilp;
    ILC *ilc;
    ILPB *ilpb;
    double dDriftFac;
    double fWeight;
    double tempI;
    double dEwaldFlop;
    FLOAT dMin,dMax,min2,max2,d2,h2;
    FLOAT rCheck[3];
    FLOAT rOffset[3];
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxCheck,nCheck;
    int iCell,iSib,iCheckCell;
    int i,ii,j,n,id,pj,nActive,nTotActive;
    int iOpen;
    int nPart,nMaxPart;
    int nCell,nMaxCell;
    int nPartBucket,nMaxPartBucket;
    double dSyncDelta;

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    */
    if (bVeryActive) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkd->kdNodes[VAROOT].pUpper - pkd->kdNodes[VAROOT].pLower + 1);
	}	
    else if (!pkdIsCellActive(pkd,&pkd->kdNodes[ROOT])) return 0;
    /*
    ** SyncDelta is used to compare the current time to the time of a cell to decide
    ** if they are synchronous. We need some sort of sensible minimum difference
    ** and we can use the maximum rung for this purpose.
    */
    dSyncDelta = pkd->param.dDelta*pow(2.0,-(pkd->param.iMaxRung+2.0));
    /*
    ** Initially we set our cell pointer to 
    ** point to the top tree.
    */
    c = pkd->kdTop;
    /*
    ** Allocate initial interaction lists.
    */
    nTotActive = 0;
    nPart = 0;
    nMaxPart = 500;
    ilp = malloc(nMaxPart*sizeof(ILP));
    assert(ilp != NULL);
    nCell = 0;
    nMaxCell = 500;
    ilc = malloc(nMaxCell*sizeof(ILC));
    assert(ilc != NULL);
    nPartBucket = 0;
    nMaxPartBucket = 500;
    ilpb = malloc(nMaxPartBucket*sizeof(ILPB));
    assert(ilpb != NULL);
    /*
    ** Allocate Checklist.
    */
    nMaxCheck = 2*nReps+1;
    nMaxCheck = nMaxCheck*nMaxCheck*nMaxCheck;	/* all replicas */
    iCell = pkd->iTopRoot;
    while ((iCell = c[iCell].iParent)) ++nMaxCheck; /* all top tree siblings */
    nMaxCheck = (nMaxCheck>1000)?nMaxCheck:1000;
    Check = malloc(nMaxCheck*sizeof(CELT));
    assert(Check != NULL);
    nCheck = 0;
    /*
    ** Allocate the stack.
    */
    S = malloc(pkd->nMaxDepth*sizeof(CSTACK));
    assert(S != NULL);
    iStack = -1;
    for (ism=0;ism<pkd->nMaxDepth;++ism) {
	S[ism].Check = malloc(nMaxCheck*sizeof(CELT));
	assert(S[ism].Check != NULL);
	}
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
		    Check[nCheck].iCell = ROOT;
		    /* If leaf of top tree, use root of
		       local tree.
		    */
		    if (c[ROOT].iLower) {
			Check[nCheck].id = -1;
			}
		    else {
			Check[nCheck].id = c[ROOT].pLower;
			}				    
		    for (j=0;j<3;++j) Check[nCheck].rOffset[j] = rOffset[j];
		    ++nCheck;
		    }
		if (bRep && bVeryActive) {
		    /*
		    ** Add the images of the very active tree to the checklist.
		    */
		    Check[nCheck].iCell = VAROOT;
		    Check[nCheck].id = mdlSelf(pkd->mdl);
		    for (j=0;j<3;++j) Check[nCheck].rOffset[j] = rOffset[j];
		    ++nCheck;
		    }
		}
	    }
	}
    if (!bVeryActive) {
	iCell = pkd->iTopRoot;
	iSib = SIBLING(iCell);
	while (iSib) {
	    if (c[iSib].iLower) {
		Check[nCheck].iCell = iSib;
		Check[nCheck].id = -1;
		}
	    else {
		/* If leaf of top tree, use root of local tree */
		Check[nCheck].iCell = ROOT;
		Check[nCheck].id = c[iSib].pLower;
		}
	    for (j=0;j<3;++j) Check[nCheck].rOffset[j] = 0.0;
	    ++nCheck;
	    iCell = c[iCell].iParent;
	    iSib = SIBLING(iCell);
	    }
	}
    /*
    ** Now switch the cell pointer to point to 
    ** the local tree.
    */
    c = pkd->kdNodes;
    if (bVeryActive) iCell = VAROOT;
    else iCell = ROOT;
    while (1) {
	while (1) {
	    /*
	    ** Process the Checklist.
	    */
	    ii = 0;
	    for (i=0;i<nCheck;++i) {
		id = Check[i].id;
		if (id == pkd->idSelf) {
		    pkdc = &pkd->kdNodes[Check[i].iCell];
		    n = pkdc->pUpper - pkdc->pLower + 1;
		    }
		else if (id < 0) {
		    pkdc = &pkd->kdTop[Check[i].iCell];
		    assert(pkdc->iLower != 0);
		    n = WALK_MINMULTIPOLE;  /* See check below */
		    }
		else {
		    pkdc = mdlAquire(pkd->mdl,CID_CELL,Check[i].iCell,id);
		    n = pkdc->pUpper - pkdc->pLower + 1;
		    }
		/*
		** If the cell is not time synchronous, then work out a drift factor
		** for this cell.
		*/
		if (fabs(pkdc->dTimeStamp-dTime) > dSyncDelta) {
		  /*
		  ** We need to account for cosmological drift factor here!
		  */
		  if (pkd->param.csm->bComove && pkd->param.bCannonical) {
		    /*
		    ** This might get called quite a bit in this code. Better might
		    ** be to store a dDriftFac within the CheckElt structure, thereby
		    ** reducing the number of calls to csmComoveDriftFac. Otherwise
		    ** we may need to speed this function up.
		    */
		    dDriftFac = csmComoveDriftFac(pkd->param.csm,pkdc->dTimeStamp,dTime - pkdc->dTimeStamp);
		  }
		  else {
		    dDriftFac = dTime - pkdc->dTimeStamp;
		  }
		  for (j=0;j<3;++j) rCheck[j] = pkdc->r[j] + 
		    dDriftFac*pkdc->v[j] + Check[i].rOffset[j];
		}
		else {
		  dDriftFac = 0.0;
		  for (j=0;j<3;++j) rCheck[j] = pkdc->r[j] + Check[i].rOffset[j];
		}
		if (c[iCell].iLower) {
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
			dMin = fabs(rCheck[j] - c[iCell].bnd.fCenter[j]);
			dMax = dMin + c[iCell].bnd.fMax[j];
			dMin -= c[iCell].bnd.fMax[j];
			if (dMin > 0) min2 += dMin*dMin;
			max2 += dMax*dMax;
			}
		    /*
		    ** By default we just keep this checkcell on the 
		    ** checklist.
		    */
		    if (max2 <= pkdc->fOpen2 || n < WALK_MINMULTIPOLE) iOpen = 1;
		    else if (min2 > pkdc->fOpen2) {
#ifdef SOFTLINEAR
			/*
			** For symmetrized softening we need to calculate the distance between
			** center of mass between the two cells.
			*/
			d2 = 0;
			for (j=0;j<3;++j) {
			    d2 += (rCheck[j] - c[iCell].r[j])*(rCheck[j] - c[iCell].r[j]);
			    }
			h2 = sqrt(pkdc->fSoft2) + sqrt(c[iCell].fSoft2);
			h2 *= h2;
			if (d2 > h2) iOpen = -1;
			else iOpen = 0;   /* in all other cases we can't decide until we get down to a bucket */
#endif
#ifdef SOFTSQUARE
			/*
			** For symmetrized softening we need to calculate the distance between
			** center of mass between the two cells.
			*/
			d2 = 0;
			for (j=0;j<3;++j) {
			    d2 += (rCheck[j] - c[iCell].r[j])*(rCheck[j] - c[iCell].r[j]);
			    }
			h2 = 2*(pkdc->fSoft2 + c[iCell].fSoft2);
			if (d2 > h2) iOpen = -1;
			else iOpen = 0;   /* in all other cases we can't decide until we get down to a bucket */
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
			if (min2 > 4*pkdc->fSoft2) iOpen = -1;
			else if (max2 < 4*pkdc->fSoft2) iOpen = -2;  /* means we accept this cell as a softened monopole */
			else iOpen = 0;  /* keep checking this cell */
#endif
			}
		    else iOpen = 0;
		    }
		else {
		    /*
		    ** If this cell is a bucket we have to either open the checkcell
		    ** and get the particles, or accept the multipole. For this
		    ** reason we only need to calculate min2.
		    */
		    min2 = 0;
		    for (j=0;j<3;++j) {
			dMin = fabs(rCheck[j] - c[iCell].bnd.fCenter[j]);
			dMin -= c[iCell].bnd.fMax[j];
			if (dMin > 0) min2 += dMin*dMin;
			}
		    /*
		    ** By default we open the cell!
		    */
		    if (min2 > pkdc->fOpen2 && n >= WALK_MINMULTIPOLE) {
#ifdef SOFTLINEAR
			/*
			** For symmetrized softening we need to calculate the distance between
			** center of mass between the two cells.
			*/
			d2 = 0;
			for (j=0;j<3;++j) {
			    d2 += (rCheck[j] - c[iCell].r[j])*(rCheck[j] - c[iCell].r[j]);
			    }
			h2 = sqrt(pkdc->fSoft2) + sqrt(c[iCell].fSoft2);
			h2 *= h2;
			if (d2 > h2) iOpen = -1;
			else iOpen = -2; /* means we treat this cell as a softened monopole */
#endif
#ifdef SOFTSQUARE
			/*
			** For symmetrized softening we need to calculate the distance between
			** center of mass between the two cells.
			*/
			d2 = 0;
			for (j=0;j<3;++j) {
			    d2 += (rCheck[j] - c[iCell].r[j])*(rCheck[j] - c[iCell].r[j]);
			    }
			h2 = 2*(pkdc->fSoft2 + c[iCell].fSoft2);
			if (d2 > h2) iOpen = -1;
			else iOpen = -2; /* means we treat this cell as a softened monopole */
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
			if (min2 > 4*pkdc->fSoft2) iOpen = -1;
			else iOpen = -2; /* means we treat this cell as a softened monopole */
#endif
			}
		    else iOpen = 1;
		    }
/*
  printf("   i:%6d iCheck:%6d id:%2d iOpen:%2d\n",i,Check[i].iCell,id,iOpen);
*/
		if (iOpen > 0) {
		    /*
		    ** Contained! (or intersected in the case of reaching the bucket)
		    */
		    iCheckCell = pkdc->iLower;
		    if (iCheckCell) {
			/*
			** Open the cell.
			** Here we ASSUME that the children of 
			** pkdc are all in sequential memory order!
			** (the new tree build assures this)
			** (also true for the top tree)
			** We could do a prefetch here for non-local
			** cells.
			*/
			if (nCheck + 2 > nMaxCheck) {
			    nMaxCheck += 1000;
			    Check = realloc(Check,nMaxCheck*sizeof(CELT));
			    assert(Check != NULL);
			    for (ism=0;ism<pkd->nMaxDepth;++ism) {
				S[ism].Check = realloc(S[ism].Check,nMaxCheck*sizeof(CELT));
				assert(S[ism].Check != NULL);
				}
			    }
			Check[nCheck] = Check[i];
			Check[nCheck+1] = Check[i];
			/*
			** If we are opening a leaf of the top tree
			** we need to correctly set the processor id.
			** (this is a bit tricky)
			*/
			if (id < 0 ) {
			    if(pkd->kdTop[iCheckCell].pLower >= 0) {
				Check[nCheck].iCell = ROOT;
				Check[nCheck].id = pkd->kdTop[iCheckCell].pLower;
				}
			    else {
				Check[nCheck].iCell = iCheckCell;
				}
			    if(pkd->kdTop[iCheckCell+1].pLower >= 0) {
				Check[nCheck+1].iCell = ROOT;
				Check[nCheck+1].id = pkd->kdTop[iCheckCell+1].pLower;
				}
			    else {
				Check[nCheck+1].iCell = iCheckCell+1;
				}

			    }
			else {
			    Check[nCheck].iCell = iCheckCell;
			    Check[nCheck+1].iCell = iCheckCell+1;
			    }
			nCheck += 2;
			}
		    else {
			/*
			** Now I am trying to open a bucket, which means I place particles on the ilp
			** interaction list. I also have to make sure that I place the
			** particles in time synchronous positions.
			*/
			if (id == pkd->idSelf) {
			    /*
			    ** Local Bucket Interaction.
			    ** Interact += Pacticles(pkdc);
			    */
			    if (nPart + n > nMaxPart) {
				nMaxPart += 500 + n;
				ilp = realloc(ilp,nMaxPart*sizeof(ILP));
				assert(ilp != NULL);	
				}
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				ilp[nPart].iOrder = p[pj].iOrder;
				ilp[nPart].m = p[pj].fMass;
				ilp[nPart].x = p[pj].r[0] + dDriftFac*p[pj].v[0] + Check[i].rOffset[0];
				ilp[nPart].y = p[pj].r[1] + dDriftFac*p[pj].v[1] + Check[i].rOffset[1];
				ilp[nPart].z = p[pj].r[2] + dDriftFac*p[pj].v[2] + Check[i].rOffset[2];
				ilp[nPart].vx = p[pj].v[0]; 
				ilp[nPart].vy = p[pj].v[1];
				ilp[nPart].vz = p[pj].v[2];
#ifdef SOFTLINEAR
				ilp[nPart].h = p[pj].fSoft;
#endif
#ifdef SOFTSQUARE
				ilp[nPart].twoh2 = 2*p[pj].fSoft*p[pj].fSoft;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
				ilp[nPart].fourh2 = 4*p[pj].fSoft*p[pj].fSoft;
#endif
				++nPart;
				}
			    }
			else {
			    /*
			    ** Remote Bucket Interaction.
			    ** Interact += Pacticles(pkdc);
			    */
			    if (nPart + n > nMaxPart) {
				nMaxPart += 500 + n;
				ilp = realloc(ilp,nMaxPart*sizeof(ILP));
				assert(ilp != NULL);	
				}
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				pRemote = mdlAquire(pkd->mdl,CID_PARTICLE,pj,id);
				ilp[nPart].iOrder = pRemote->iOrder;
				ilp[nPart].m = pRemote->fMass;
				ilp[nPart].x = pRemote->r[0] + dDriftFac*pRemote->v[0] + Check[i].rOffset[0];
				ilp[nPart].y = pRemote->r[1] + dDriftFac*pRemote->v[1] + Check[i].rOffset[1];
				ilp[nPart].z = pRemote->r[2] + dDriftFac*pRemote->v[2] + Check[i].rOffset[2];
				ilp[nPart].vx = pRemote->v[0]; 
				ilp[nPart].vy = pRemote->v[1];
				ilp[nPart].vz = pRemote->v[2];
#ifdef SOFTLINEAR
				ilp[nPart].h = pRemote->fSoft;
#endif
#ifdef SOFTSQUARE
				ilp[nPart].twoh2 = 2*pRemote->fSoft*pRemote->fSoft;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
				ilp[nPart].fourh2 = 4*pRemote->fSoft*pRemote->fSoft;
#endif
				++nPart;
				mdlRelease(pkd->mdl,CID_PARTICLE,pRemote);
				}
			    }
			/*
			** ...and we need to consider it in the timestepping part so place this cell
			** onto the particle-bucket list. This list is not used for force evaluation 
			** though.
			*/
			if (nPartBucket == nMaxPartBucket) {
			    nMaxPartBucket += 500;
			    ilpb = realloc(ilpb,nMaxPartBucket*sizeof(ILPB));
			    assert(ilpb != NULL);
			    }
			ilpb[nPartBucket].x = rCheck[0];
			ilpb[nPartBucket].y = rCheck[1];
			ilpb[nPartBucket].z = rCheck[2];
			ilpb[nPartBucket].m = pkdc->mom.m;   /* we really only need the mass here */
#ifdef SOFTLINEAR
		        ilpb[nPartBucket].h = sqrt(pkdc->fSoft2);
#endif
#ifdef SOFTSQUARE
		        ilpb[nPartBucket].twoh2 = 2*pkdc->fSoft2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
		        ilpb[nPartBucket].fourh2 = 4*pkdc->fSoft2;
#endif
			++nPartBucket;
			}  /* end of opening a bucket */
		    }
		else if (iOpen == -1) {
		    /*
		    ** No intersection, accept multipole!
		    ** Interact += Moment(pkdc);
		    */
		    if (nCell == nMaxCell) {
			nMaxCell += 500;
			ilc = realloc(ilc,nMaxCell*sizeof(ILC));
			assert(ilc != NULL);
			}
		    ilc[nCell].x = rCheck[0];
		    ilc[nCell].y = rCheck[1];
		    ilc[nCell].z = rCheck[2];
		    ilc[nCell].mom = pkdc->mom;
#ifdef HERMITE
		    ilc[nCell].vx = pkdc->v[0];
		    ilc[nCell].vy = pkdc->v[1];
		    ilc[nCell].vz = pkdc->v[2];
#endif
		    ++nCell;
		    }
		else if (iOpen == -2) {
		    /*
		    ** We accept this multipole from the opening criterion, but it is a softened
		    ** interaction, so we need to treat is as a softened monopole by putting it
		    ** on the particle interaction list.
		    */
		    if (nPart == nMaxPart) {
			nMaxPart += 500;
			ilp = realloc(ilp,nMaxPart*sizeof(ILP));
			assert(ilp != NULL);	
			}
		    ilp[nPart].iOrder = -1; /* set iOrder to negative value for time step criterion */
		    ilp[nPart].m = pkdc->mom.m;
		    ilp[nPart].x = rCheck[0];
		    ilp[nPart].y = rCheck[1];
		    ilp[nPart].z = rCheck[2];
		    ilp[nPart].vx = pkdc->v[0];
		    ilp[nPart].vy = pkdc->v[1];
		    ilp[nPart].vz = pkdc->v[2];
#ifdef SOFTLINEAR
		    ilp[nPart].h = sqrt(pkdc->fSoft2);
#endif
#ifdef SOFTSQUARE
		    ilp[nPart].twoh2 = 2*pkdc->fSoft2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
		    ilp[nPart].fourh2 = 4*pkdc->fSoft2;
#endif
		    ++nPart;
		    /*
		    ** ...and we need to consider it in the timestepping part so place this cell
		    ** onto the particle-bucket list. This list is not used for force evaluation 
		    ** though.
		    */
		    if (nPartBucket == nMaxPartBucket) {
			nMaxPartBucket += 500;
			ilpb = realloc(ilpb,nMaxPartBucket*sizeof(ILPB));
			assert(ilpb != NULL);
			}
		    ilpb[nPartBucket].x = rCheck[0];
		    ilpb[nPartBucket].y = rCheck[1];
		    ilpb[nPartBucket].z = rCheck[2];
		    ilpb[nPartBucket].m = pkdc->mom.m;   /* we really only need the mass here */
#ifdef SOFTLINEAR
		    ilpb[nPartBucket].h = sqrt(pkdc->fSoft2);
#endif
#ifdef SOFTSQUARE
		    ilpb[nPartBucket].twoh2 = 2*pkdc->fSoft2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
		    ilpb[nPartBucket].fourh2 = 4*pkdc->fSoft2;
#endif
		    ++nPartBucket;
		    }
		else {
		    Check[ii++] = Check[i];
		    }
		if (id >= 0 && id != pkd->idSelf) {
		    mdlRelease(pkd->mdl,CID_CELL,pkdc);
		    }
		}
	    nCheck = ii;
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (!c[iCell].iLower) break;
	    iCell = c[iCell].iLower;
	    /*
	    ** Now add the siblings of iCell to the 
	    ** Checklist.
	    */
	    if (nCheck == nMaxCheck) {
		nMaxCheck += 1000;
		Check = realloc(Check,nMaxCheck*sizeof(CELT));
		assert(Check != NULL);
		for (ism=0;ism<pkd->nMaxDepth;++ism) {
		    S[ism].Check = realloc(S[ism].Check,nMaxCheck*sizeof(CELT));
		    assert(S[ism].Check != NULL);
		    }
		}
	    /*
	    ** Check iCell is active. We eventually want to just to a 
	    ** rung check here when we start using tree repair, but 
	    ** for now this is just as good.
	    */
	    if (pkdIsCellActive(pkd,c+iCell)) {
		/*
		** iCell is active, continue processing it.
		** Put the sibling onto the checklist.
		*/
		Check[nCheck].iCell = iCell + 1;
		Check[nCheck].id = pkd->idSelf;
		for (j=0;j<3;++j) Check[nCheck].rOffset[j] = 0.0;
		++nCheck;
		/*
		** Test whether the sibling is active as well.
		** If not we don't push it onto the stack, but we
		** have to be careful to not pop the stack when we
		** hit the sibling. See the goto InactiveAscend below
		** for how this is done.
		*/
		if (pkdIsCellActive(pkd,c+iCell+1)) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    S[iStack].nPart = nPart;
		    S[iStack].nCell = nCell;
		    S[iStack].nPartBucket = nPartBucket;
		    S[iStack].nCheck = nCheck;
		    for (i=0;i<nCheck;++i) S[iStack].Check[i] = Check[i];
		    S[iStack].Check[nCheck-1].iCell = iCell;
		    }
		}
	    else {
		/*
		** Skip iCell, but add it to the Checklist.
		** No need to push anything onto the stack here.
		*/
		Check[nCheck].iCell = iCell;
		Check[nCheck].id = pkd->idSelf;
		for (j=0;j<3;++j) Check[nCheck].rOffset[j] = 0.0;
		++nCheck;
		/*
		** Move onto processing the sibling.
		*/
		++iCell;
		}
	    }
	/*
	** Now the interaction list should be complete and the 
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
	assert(nCheck == 0);
	/*
	** We no longer add *this bucket to any interaction list, this is now done with an 
	** N(N-1)/2 loop in pkdBucketInteract().
	*/
	pkdc = &c[iCell];
	/*
	** Now calculate gravity on this bucket!
	*/
	tempI = *pdFlop;
	nActive = pkdGravInteract(pkd,pkdc,NULL,ilp,nPart,ilc,nCell,ilpb,nPartBucket,0,0,pdFlop);
	/*
	** Note that if Ewald is being performed we need to factor this
	** constant cost into the load balancing weights.
	**
	** CAREFUL: pkdBucketEwald MUST follow pkdGravInteract since we want to use the 
	** particle's old acceleration in pkdGravInteract!
	*/
	if (bEwald) {
	    dEwaldFlop = pkdBucketEwald(pkd,&pkd->kdNodes[iCell],nReps,fEwCut,bEwaldKick);
	    if (!bEwaldKick) {
		*pdFlop += dEwaldFlop;
	    }
	}
	if (nActive) {
	    fWeight = (*pdFlop-tempI)/nActive;
	    pkdBucketWeight(pkd,iCell,fWeight);
	    *pdPartSum += nActive*nPart;
	    *pdCellSum += nActive*nCell;
	    nTotActive += nActive;
	    }

	while (iCell & 1) {
	InactiveAscend:
	    iCell = c[iCell].iParent;
	    if (!iCell) {
		/*
		** Make sure stack is empty and free its storage.
		*/
		assert(iStack == -1);
		for (ism=0;ism<pkd->nMaxDepth;++ism) {
		    free(S[ism].Check);
		    }
		free(S);
		/*
		** Free checklist storage.
		*/
		free(Check);
		/*
		** Free interaction lists.
		*/
		free(ilp);
		free(ilc);
		free(ilpb);
		return(nTotActive);
		}
	    }
	++iCell;
	if (!pkdIsCellActive(pkd,c+iCell)) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	nPart = S[iStack].nPart;
	nCell = S[iStack].nCell;
	nPartBucket = S[iStack].nPartBucket;
	nCheck = S[iStack].nCheck;
	for (i=0;i<nCheck;++i) Check[i] = S[iStack].Check[i];
	--iStack;
	}
    }

