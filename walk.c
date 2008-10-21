#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *walk_module_id = "$Id$";

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

#define WALK_MINMULTIPOLE	4

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,
		int bVeryActive,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p;
    PARTICLE *pRemote;
    KDN *c;
    KDN *pkdc;
    double dDriftFac;
    double fWeight;
    double tempI;
    double dEwFlop;
    double dRhoFac;
    FLOAT dMin,dMax,min2,max2;
#if defined(SOFTLINEAR) || defined(SOFTSQUARE)
    FLOAT d2,h2;
#endif
    FLOAT rCheck[3];
    FLOAT rOffset[3];
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxInitCheck,nCheck;
    int iCell,iSib,iCheckCell;
    int i,ii,j,n,id,pj,nActive,nTotActive;
    int iOpen;
    int nPart;
    int nCell;
    double dSyncDelta;
    double *v;

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    */
    if (bVeryActive) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkd->kdNodes[VAROOT].pUpper - pkd->kdNodes[VAROOT].pLower + 1);
	}
    else if (!pkdIsCellActive(&pkd->kdNodes[ROOT],uRungLo,uRungHi)) return 0;
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
    nCell = 0;
    nMaxInitCheck = 2*nReps+1;
    nMaxInitCheck = nMaxInitCheck*nMaxInitCheck*nMaxInitCheck;	/* all replicas */
    iCell = pkd->iTopRoot;
    while ((iCell = c[iCell].iParent)) ++nMaxInitCheck; /* all top tree siblings */
    assert(nMaxInitCheck < pkd->nMaxCheck);
    nCheck = 0;
    iStack = -1;
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
		    if (c[ROOT].iLower) {
			pkd->Check[nCheck].id = -1;
			}
		    else {
			pkd->Check[nCheck].id = c[ROOT].pLower;
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
	iCell = pkd->iTopRoot;
	iSib = SIBLING(iCell);
	while (iSib) {
	    if (c[iSib].iLower) {
		pkd->Check[nCheck].iCell = iSib;
		pkd->Check[nCheck].id = -1;
		}
	    else {
		/* If leaf of top tree, use root of local tree */
		pkd->Check[nCheck].iCell = ROOT;
		pkd->Check[nCheck].id = c[iSib].pLower;
		}
	    for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
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
		id = pkd->Check[i].id;
		if (id == pkd->idSelf) {
		    pkdc = &pkd->kdNodes[pkd->Check[i].iCell];
		    n = pkdc->pUpper - pkdc->pLower + 1;
		    }
		else if (id < 0) {
		    pkdc = &pkd->kdTop[pkd->Check[i].iCell];
		    assert(pkdc->iLower != 0);
		    n = WALK_MINMULTIPOLE;  /* See check below */
		    }
		else {
		    pkdc = mdlAquire(pkd->mdl,CID_CELL,pkd->Check[i].iCell,id);
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
                    if (pkd->param.csm->bComove) {
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
			dDriftFac*pkdc->v[j] + pkd->Check[i].rOffset[j];
		    }
		else {
		    dDriftFac = 0.0;
		    for (j=0;j<3;++j) rCheck[j] = pkdc->r[j] + pkd->Check[i].rOffset[j];
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
		  printf("   i:%6d iCheck:%6d id:%2d iOpen:%2d\n",i,pkd->Check[i].iCell,id,iOpen);
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
			if (nCheck + 2 > pkd->nMaxCheck) {
			    pkd->nMaxCheck += 1000;
			    pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->Check != NULL);
			    for (ism=0;ism<pkd->nMaxStack;++ism) {
				pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
				assert(pkd->S[ism].Check != NULL);
				}
			    }
			pkd->Check[nCheck] = pkd->Check[i];
			pkd->Check[nCheck+1] = pkd->Check[i];
			/*
			** If we are opening a leaf of the top tree
			** we need to correctly set the processor id.
			** (this is a bit tricky)
			*/
			if (id < 0 ) {
			    if (pkd->kdTop[iCheckCell].pLower >= 0) {
				pkd->Check[nCheck].iCell = ROOT;
				pkd->Check[nCheck].id = pkd->kdTop[iCheckCell].pLower;
				}
			    else {
				pkd->Check[nCheck].iCell = iCheckCell;
				}
			    if (pkd->kdTop[iCheckCell+1].pLower >= 0) {
				pkd->Check[nCheck+1].iCell = ROOT;
				pkd->Check[nCheck+1].id = pkd->kdTop[iCheckCell+1].pLower;
				}
			    else {
				pkd->Check[nCheck+1].iCell = iCheckCell+1;
				}

			    }
			else {
			    pkd->Check[nCheck].iCell = iCheckCell;
			    pkd->Check[nCheck+1].iCell = iCheckCell+1;
			    }
			nCheck += 2;
			}
		    else {
			/*
			** Now I am trying to open a bucket, which means I place particles on the pkd->ilp
			** interaction list. I also have to make sure that I place the
			** particles in time synchronous positions.
			*/
			if (id == pkd->idSelf) {
			    /*
			    ** Local Bucket Interaction.
			    ** Interact += Pacticles(pkdc);
			    */
			    if (nPart + n > pkd->nMaxPart) {
				pkd->nMaxPart += 500 + n;
				pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
				assert(pkd->ilp != NULL);
				}
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				p = pkdParticle(pkd,pj);
				v = pkdVel(pkd,p);
				pkd->ilp[nPart].iOrder = p->iOrder;
				pkd->ilp[nPart].m = pkdMass(pkd,p);
				pkd->ilp[nPart].x = p->r[0] + dDriftFac*v[0] + pkd->Check[i].rOffset[0];
				pkd->ilp[nPart].y = p->r[1] + dDriftFac*v[1] + pkd->Check[i].rOffset[1];
				pkd->ilp[nPart].z = p->r[2] + dDriftFac*v[2] + pkd->Check[i].rOffset[2];
				pkd->ilp[nPart].vx = v[0];
				pkd->ilp[nPart].vy = v[1];
				pkd->ilp[nPart].vz = v[2];
#ifdef SOFTLINEAR
				pkd->ilp[nPart].h = p->fSoft;
#endif
#ifdef SOFTSQUARE
				pkd->ilp[nPart].twoh2 = 2*p->fSoft*p->fSoft;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
				pkd->ilp[nPart].fourh2 = 4*pkdSoft(pkd,p)*pkdSoft(pkd,p);
#endif
				++nPart;
				}
			    }
			else {
			    /*
			    ** Remote Bucket Interaction.
			    ** Interact += Pacticles(pkdc);
			    */
			    if (nPart + n > pkd->nMaxPart) {
				pkd->nMaxPart += 500 + n;
				pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
				assert(pkd->ilp != NULL);
				}
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				pRemote = mdlAquire(pkd->mdl,CID_PARTICLE,pj,id);
				v = pkdVel(pkd,pRemote);
				pkd->ilp[nPart].iOrder = pRemote->iOrder;
				pkd->ilp[nPart].m = pkdMass(pkd,pRemote);
				pkd->ilp[nPart].x = pRemote->r[0] + dDriftFac*v[0] + pkd->Check[i].rOffset[0];
				pkd->ilp[nPart].y = pRemote->r[1] + dDriftFac*v[1] + pkd->Check[i].rOffset[1];
				pkd->ilp[nPart].z = pRemote->r[2] + dDriftFac*v[2] + pkd->Check[i].rOffset[2];
				pkd->ilp[nPart].vx = v[0];
				pkd->ilp[nPart].vy = v[1];
				pkd->ilp[nPart].vz = v[2];
#ifdef SOFTLINEAR
				pkd->ilp[nPart].h = pRemote->fSoft;
#endif
#ifdef SOFTSQUARE
				pkd->ilp[nPart].twoh2 = 2*pRemote->fSoft*pRemote->fSoft;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
				pkd->ilp[nPart].fourh2 = 4*pkdSoft(pkd,pRemote)*pkdSoft(pkd,pRemote);
#endif
				++nPart;
				mdlRelease(pkd->mdl,CID_PARTICLE,pRemote);
				}
			    }
			}  /* end of opening a bucket */
		    }
		else if (iOpen == -1) {
		    /*
		    ** No intersection, accept multipole!
		    ** Interact += Moment(pkdc);
		    */
		    if (nCell == pkd->nMaxCell) {
			pkd->nMaxCell += 500;
			pkd->ilc = realloc(pkd->ilc,pkd->nMaxCell*sizeof(ILC));
			assert(pkd->ilc != NULL);
			}
		    pkd->ilc[nCell].x = rCheck[0];
		    pkd->ilc[nCell].y = rCheck[1];
		    pkd->ilc[nCell].z = rCheck[2];
		    pkd->ilc[nCell].mom = pkdc->mom;
#ifdef HERMITE
		    pkd->ilc[nCell].vx = pkdc->v[0];
		    pkd->ilc[nCell].vy = pkdc->v[1];
		    pkd->ilc[nCell].vz = pkdc->v[2];
#endif
		    ++nCell;
		    }
		else if (iOpen == -2) {
		    /*
		    ** We accept this multipole from the opening criterion, but it is a softened
		    ** interaction, so we need to treat is as a softened monopole by putting it
		    ** on the particle interaction list.
		    */
		    if (nPart == pkd->nMaxPart) {
			pkd->nMaxPart += 500;
			pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
			assert(pkd->ilp != NULL);
			}
		    pkd->ilp[nPart].iOrder = -1; /* set iOrder to negative value for time step criterion */
		    pkd->ilp[nPart].m = pkdc->mom.m;
		    pkd->ilp[nPart].x = rCheck[0];
		    pkd->ilp[nPart].y = rCheck[1];
		    pkd->ilp[nPart].z = rCheck[2];
		    pkd->ilp[nPart].vx = pkdc->v[0];
		    pkd->ilp[nPart].vy = pkdc->v[1];
		    pkd->ilp[nPart].vz = pkdc->v[2];
#ifdef SOFTLINEAR
		    pkd->ilp[nPart].h = sqrt(pkdc->fSoft2);
#endif
#ifdef SOFTSQUARE
		    pkd->ilp[nPart].twoh2 = 2*pkdc->fSoft2;
#endif
#if !defined(SOFTLINEAR) && !defined(SOFTSQUARE)
		    pkd->ilp[nPart].fourh2 = 4*pkdc->fSoft2;
#endif
		    ++nPart;
		    }
		else {
		    pkd->Check[ii++] = pkd->Check[i];
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
	    if (nCheck == pkd->nMaxCheck) {
		pkd->nMaxCheck += 1000;
		pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
		assert(pkd->Check != NULL);
		for (ism=0;ism<pkd->nMaxStack;++ism) {
		    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
		    assert(pkd->S[ism].Check != NULL);
		    }
		}
	    /*
	    ** Check iCell is active. We eventually want to just to a
	    ** rung check here when we start using tree repair, but
	    ** for now this is just as good.
	    */
	    if (pkdIsCellActive(c+iCell,uRungLo,uRungHi)) {
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
		if (pkdIsCellActive(c+iCell+1,uRungLo,uRungHi)) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    pkd->S[iStack].nPart = nPart;
		    pkd->S[iStack].nCell = nCell;
		    pkd->S[iStack].nCheck = nCheck;
		    for (i=0;i<nCheck;++i) pkd->S[iStack].Check[i] = pkd->Check[i];
		    pkd->S[iStack].Check[nCheck-1].iCell = iCell;
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
	tempI += dEwFlop;
	nActive = pkdGravInteract(pkd,uRungLo,uRungHi,pkdc,NULL,pkd->ilp,nPart,pkd->ilc,nCell,
				  0.0,0.0,bEwald,pdFlop,&dEwFlop,dRhoFac);

	if (nActive) {
	    fWeight = (*pdFlop-tempI+dEwFlop)/nActive;
	    /*
	    ** Here we used to set the weight of each particle based on the work done, but now we just
	    ** assume that all active particles cost the same.
	    */
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
		*pdFlop += dEwFlop;   /* Finally add the ewald score to get a proper float count */
		return(nTotActive);
		}
	    }
	++iCell;
	if (!pkdIsCellActive(c+iCell,uRungLo,uRungHi)) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	nPart = pkd->S[iStack].nPart;
	nCell = pkd->S[iStack].nCell;
	nCheck = pkd->S[iStack].nCheck;
	for (i=0;i<nCheck;++i) pkd->Check[i] = pkd->S[iStack].Check[i];
	--iStack;
	}
    }

