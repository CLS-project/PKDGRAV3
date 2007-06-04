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
#include "pkd.h"
#include "walk.h"
#include "grav.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "ewald.h"
#include "moments.h"
#ifdef TIME_WALK_WORK
#include <sys/time.h>
#endif


typedef struct CheckElt {
    int iCell;
    int id;
    FLOAT rOffset[3];
    FLOAT fOpen;
    } CELT;

typedef struct CheckStack {
    int nPart;
    int nCell;
    int nPartBucket;
    int nCheck;
    CELT *Check;
    LOCR L;
    double fWeight;
    } CSTACK;

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

double dMonopoleThetaFac = 1.5;

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,double dTime,int nReps,int bEwald,int bVeryActive,double fEwCut,
		double *pdFlop,double *pdPartSum,double *pdCellSum)
    {
    PARTICLE *p = pkd->pStore;
    PARTICLE *pRemote;
    KDN *c;
    KDN *pkdc, *next_pkdc;
    CSTACK *S;
    CELT *Check;
    LOCR L;
    ILP *ilp;
    ILC *ilc;
    double fWeight = 0.0;
    FLOAT dMin,dMax,min2,max2,d2,fourh2;
    double dDriftFac;
    FLOAT rCheck[3];
    FLOAT rOffset[3];
    FLOAT xParent,yParent,zParent;
    FLOAT dx[3],dir;
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxCheck,nCheck;
    int iCell,iSib,iCheckCell;
    int i,ii,j,n,id,pj,nActive,nTotActive;
    int iOpen;
    int nPart,nMaxPart;
    int nCell,nMaxCell;
#ifdef USE_SIMD_MOMR
    int ig,iv;
#endif
#ifdef USE_SIMD_LOCR
    v4sf vdir;
    float  sdir;
    int nGlam,nMaxGlam;
    GLAM *ilglam;
#ifdef GLAM_STATS
    int nTotGlam=0, nAvgGlam=0;
#endif
#else
    momFloat t1, t2, t3r, t4r;
#endif
#ifdef TIME_WALK_WORK
    TIMER tv;
#else
    double tempI;
#endif
    double dSyncDelta;

#ifdef PROFILE_GRAVWALK
    VTResume();
#endif

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    */
    if (bVeryActive) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkd->kdNodes[VAROOT].pUpper - pkd->kdNodes[VAROOT].pLower + 1);
	}	
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
#ifdef USE_SIMD_MOMR
    nMaxCell = 1000;
    ilc = SIMD_malloc(nMaxCell/4*sizeof(ILC));
#else
    nMaxCell = 500;
    ilc = malloc(nMaxCell*sizeof(ILC));
#endif
    assert(ilc != NULL);
#ifdef USE_SIMD_LOCR
    nMaxGlam = 1000;
    ilglam = SIMD_malloc(nMaxGlam*sizeof(GLAM));
    assert(ilglam != 0);
    nGlam = 0;
#endif

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
    ** Clear local expansion.
    */
    momClearLocr(&L);
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

#ifdef TIME_WALK_WORK
	    clearTimer(&tv);
#else
	    tempI = *pdFlop;
#endif
	    ii = 0;

	    for (i=0;i<nCheck;++i) {
		id = Check[i].id;

#ifdef __SSE__strange
		    /* Does nothing */
		    _mm_prefetch( (char *)(c+iCell)+0,_MM_HINT_T0 );
		    _mm_prefetch( (char *)(c+iCell)+64,_MM_HINT_T0 );
		    _mm_prefetch( (char *)(c+iCell)+128,_MM_HINT_T0 );
#endif
		if (id == pkd->idSelf) {
		    pkdc = pkd->kdNodes + Check[i].iCell;
		    n = pkdc->pUpper - pkdc->pLower + 1;
		    }
		else if (id < 0) {
		    pkdc = pkd->kdTop + Check[i].iCell;
		    assert(pkdc->iLower != 0);
		    n = WALK_MINMULTIPOLE;  /* See check below */
		    }
		else {
#ifdef TIME_WALK_WORK
		    stopTimer(&tv);
#endif
		    pkdc = mdlAquire(pkd->mdl,CID_CELL,Check[i].iCell,id);
#ifdef TIME_WALK_WORK
		    startTimer(&tv);
#endif
		    n = pkdc->pUpper - pkdc->pLower + 1;
		    }
#if 0
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
#else
		dDriftFac = 0.0;
		for (j=0;j<3;++j) rCheck[j] = pkdc->r[j] + Check[i].rOffset[j];
#endif
		d2 = 0;
		for (j=0;j<3;++j) {
		    d2 += (rCheck[j] - c[iCell].r[j])*(rCheck[j] - c[iCell].r[j]);
		}
		iOpen = 0;
		if (d2 > (c[iCell].fOpen + pkdc->fOpen)*(c[iCell].fOpen + pkdc->fOpen)) {
		    /*
		    ** Accept local expansion, but check softening.
		    */
		    fourh2 = softmassweight(c[iCell].mom.m,4*c[iCell].fSoft2,pkdc->mom.m,4*pkdc->fSoft2);
		    if (d2 > fourh2) {
			/*
			** Local expansion accepted!
			** Add to the GLAM list to be evaluated later.
			*/
#ifdef USE_SIMD_LOCR
			if (nGlam == nMaxGlam) {
			    nMaxGlam += 1000;
			    ilglam = SIMD_realloc(ilglam,
						  (nMaxGlam-1000)*sizeof(GLAM),
						  nMaxGlam*sizeof(GLAM));
			    assert(ilglam != 0);
			}
#ifdef __SSE__
			vdir = _mm_rsqrt_ss(_mm_set_ss(d2));
			/* Better: sdir = _mm_cvtss_f32(vdir); */
			_mm_store_ss(&sdir,vdir);
			sdir *= ((3.0 - sdir * sdir * (float)d2) * 0.5);
#else
			sdir = 1.0/sqrt(d2);
#endif
			ilglam[nGlam].q = pkdc->mom;
			ilglam[nGlam].dir = sdir;
			ilglam[nGlam].g0 = -sdir;
			ilglam[nGlam].t1 = -sdir;
			ilglam[nGlam].t2 = -3*sdir;
			ilglam[nGlam].t3r = -5;
			ilglam[nGlam].t4r = -7;
			ilglam[nGlam].x = dx[0];
			ilglam[nGlam].y = dx[1];
			ilglam[nGlam].z = dx[2];
			ilglam[nGlam].zero = 0.0;
			++nGlam;
#else
#if 1
			dir = 1.0/sqrt(d2);
			t1 = -dir;
			t2 = -3*dir;
			t3r = -5;
			t4r = -7;
			momGenLocrAddMomr(&L,&pkdc->mom,dir,-dir,
					  t1,t2,t3r,t4r,dx[0],dx[1],dx[2]);
#else
			momLocrAddMomr(&L,&pkdc->mom,dir,dx[0],dx[1],dx[2]);
#endif
#endif
		    }
		    else {
			/*
			** We want to test if it can be used as a softened monopole.
			** Now we calculate minimum distance from the cm of the 
			** checkcell to the edge of the current cell's bounding box.
			*/
			min2 = 0;	
			for (j=0;j<3;++j) {
			    dMin = fabs(rCheck[j] - c[iCell].bnd.fCenter[j]);
			    dMin -= c[iCell].bnd.fMax[j];
			    if (dMin > 0) min2 += dMin*dMin;
			}
			if (min2 > pkdc->fOpen*dMonopoleThetaFac*pkdc->fOpen*dMonopoleThetaFac) {
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
#ifndef USE_SIMD
			    ilp[nPart].iOrder = -1; /* set iOrder to negative value for time step criterion */
#endif
			    ilp[nPart].m = pkdc->mom.m;
			    ilp[nPart].x = rCheck[0];
			    ilp[nPart].y = rCheck[1];
			    ilp[nPart].z = rCheck[2];
#ifndef USE_SIMD
			    ilp[nPart].vx = pkdc->v[0];
			    ilp[nPart].vy = pkdc->v[1];
			    ilp[nPart].vz = pkdc->v[2];
#endif
			    ilp[nPart].fourh2 = 4*pkdc->fSoft2;
			    ++nPart;
			}
			else {
			    /*
			    ** Otherwise we must open this checkcell.
			    */
			    iOpen = 1;
			}
		    }
		} /* end of basic accept local expansion */
		else {
		    /*
		    ** If the checklist has the larger fOpen then Open it, otherwise keep it on 
		    ** the checklist (open the current cell eventually).
		    */
		    if (pkdc->fOpen > c[iCell].fOpen) {
			if (pkdc->iLower) {
			    iOpen = 1;
			}
			else {
			    /*
			    ** We can't open pkdc.  At this point we would prefer to accept
			    ** the cell as a multipole expansion, but if that isn't possible,
			    ** we open it.  If we were to always open it, that leads to
			    ** degenerate cases where a large bucket on the checklist will
			    ** result in most particles being added to the interaction list.
			    ** This happens in cosmological simulations where the outer regions
			    ** are binned and have a few very large particles.
			    */
			    if (n >= WALK_MINMULTIPOLE) {
				min2 = 0;	
				for (j=0;j<3;++j) {
				    dMin = fabs(rCheck[j] - c[iCell].bnd.fCenter[j]);
				    dMin -= c[iCell].bnd.fMax[j];
				    if (dMin > 0) min2 += dMin*dMin;
				}
				if (min2 > pkdc->fOpen*pkdc->fOpen) {
				    fourh2 = softmassweight(c[iCell].mom.m,4*c[iCell].fSoft2,pkdc->mom.m,4*pkdc->fSoft2);
				    if (min2 > fourh2) {
					/*
					** The multipole is also ok as far as softening goes, so accept it.
					*/
					iOpen = -2;
				    }
				    else if (min2 > pkdc->fOpen*dMonopoleThetaFac*pkdc->fOpen*dMonopoleThetaFac) {
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
#ifndef USE_SIMD
					ilp[nPart].iOrder = -1; /* set iOrder to negative value for time step criterion */
#endif
					ilp[nPart].m = pkdc->mom.m;
					ilp[nPart].x = rCheck[0];
					ilp[nPart].y = rCheck[1];
					ilp[nPart].z = rCheck[2];
#ifndef USE_SIMD
					ilp[nPart].vx = pkdc->v[0];
					ilp[nPart].vy = pkdc->v[1];
					ilp[nPart].vz = pkdc->v[2];
#endif
					ilp[nPart].fourh2 = 4*pkdc->fSoft2;
					++nPart;
				    }
				    else {
					/*
					** Unfortunately the multipole does not meet the softening criteria for 
					** an unsoftened hexadecapole nor for a softened monopole. We must open it.
					*/
					iOpen = 1;
				    }
				}
				else {
				    /*
				    ** This bucket cannot be accepted as a multipole.
				    ** Opening will produce particles on the particle interaction list.
				    */
				    iOpen = 1;
				}
			    }
			    else {
				/*
				** This bucket did not have enough particle to make it worth accepting as a
				** multipole, since it is faster to simply add P-P interactions at this stage.
				*/
				iOpen = 1;
			    }
			}
		    }
		    /*
		    ** From here we know that this c[iCell] has the larger opening radius.
		    */
		    else if (c[iCell].iLower) {
			Check[ii++] = Check[i];
		    }
		    else {
			/*
			** In this case we cannot open the current cell despite it having the
			** larger opening radius. We must try to use pkdc as a P-C interaction in 
			** this case, otherwise we have the danger of opening too many cells.
			*/
			if (n >= WALK_MINMULTIPOLE) {
			    min2 = 0;	
			    for (j=0;j<3;++j) {
				dMin = fabs(rCheck[j] - c[iCell].bnd.fCenter[j]);
				dMin -= c[iCell].bnd.fMax[j];
				if (dMin > 0) min2 += dMin*dMin;
			    }
			    if (min2 > pkdc->fOpen*pkdc->fOpen) {
				fourh2 = softmassweight(c[iCell].mom.m,4*c[iCell].fSoft2,pkdc->mom.m,4*pkdc->fSoft2);
				if (min2 > fourh2) {
				    /*
				    ** The multipole is also ok as far as softening goes, so accept it.
				    */
				    iOpen = -2;
				}
				else if (min2 > pkdc->fOpen*dMonopoleThetaFac*pkdc->fOpen*dMonopoleThetaFac) {
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
#ifndef USE_SIMD
				    ilp[nPart].iOrder = -1; /* set iOrder to negative value for time step criterion */
#endif
				    ilp[nPart].m = pkdc->mom.m;
				    ilp[nPart].x = rCheck[0];
				    ilp[nPart].y = rCheck[1];
				    ilp[nPart].z = rCheck[2];
#ifndef USE_SIMD
				    ilp[nPart].vx = pkdc->v[0];
				    ilp[nPart].vy = pkdc->v[1];
				    ilp[nPart].vz = pkdc->v[2];
#endif
				    ilp[nPart].fourh2 = 4*pkdc->fSoft2;
				    ++nPart;
				}
				else {
				    /*
				    ** Unfortunately the multipole does not meet the softening criteria for 
				    ** an unsoftened hexadecapole nor for a softened monopole. We must open it.
				    */
				    iOpen = 1;
				}
			    }
			    else {
				/*
				** We don't accept the particle-cell interaction, open the check cell!
				*/
				iOpen = 1;
			    }
			}
			else {
			    iOpen = 1;
			}
		    }
		}
		if (iOpen > 0) {
		    /*
		    ** Here we go through the opening of a checkcell!
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
			** interaction list.
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
#ifndef USE_SIMD
				ilp[nPart].iOrder = p[pj].iOrder;
#endif
				ilp[nPart].m = p[pj].fMass;
				ilp[nPart].x = p[pj].r[0] + dDriftFac*p[pj].v[0] + Check[i].rOffset[0];
				ilp[nPart].y = p[pj].r[1] + dDriftFac*p[pj].v[1] + Check[i].rOffset[1];
				ilp[nPart].z = p[pj].r[2] + dDriftFac*p[pj].v[2] + Check[i].rOffset[2];
#ifndef USE_SIMD
				ilp[nPart].vx = p[pj].v[0]; 
				ilp[nPart].vy = p[pj].v[1];
				ilp[nPart].vz = p[pj].v[2];
#endif
				ilp[nPart].fourh2 = 4*p[pj].fSoft*p[pj].fSoft;
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
#ifdef TIME_WALK_WORK
			    stopTimer(&tv);
#endif
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				pRemote = mdlAquire(pkd->mdl,CID_PARTICLE,pj,id);
#ifndef USE_SIMD
				ilp[nPart].iOrder = pRemote->iOrder;
#endif
				ilp[nPart].m = pRemote->fMass;
				ilp[nPart].x = pRemote->r[0] + dDriftFac*pRemote->v[0] + Check[i].rOffset[0];
				ilp[nPart].y = pRemote->r[1] + dDriftFac*pRemote->v[1] + Check[i].rOffset[1];
				ilp[nPart].z = pRemote->r[2] + dDriftFac*pRemote->v[2] + Check[i].rOffset[2];
#ifndef USE_SIMD
				ilp[nPart].vx = pRemote->v[0]; 
				ilp[nPart].vy = pRemote->v[1];
				ilp[nPart].vz = pRemote->v[2];
#endif
				ilp[nPart].fourh2 = 4*pRemote->fSoft*pRemote->fSoft;
				++nPart;
				mdlRelease(pkd->mdl,CID_PARTICLE,pRemote);
			    }
#ifdef TIME_WALK_WORK
			    startTimer(&tv);
#endif
			}  /* end of opening a bucket */
		    }
		}
		else if (iOpen == -2) {
		    /*
		    ** No intersection, accept multipole!
		    ** Interact += Moment(pkdc);
		    */
		    if (nCell == nMaxCell) {
			nMaxCell += 500;
#ifdef USE_SIMD_MOMR
			ilc = SIMD_realloc(ilc, (nMaxCell-500)/4*sizeof(ILC),
					   nMaxCell/4*sizeof(ILC));
#else
			ilc = realloc(ilc,nMaxCell*sizeof(ILC));
#endif
			assert(ilc != NULL);
			}
#ifdef USE_SIMD_MOMR
		    ig = nCell >> 2;
		    iv = nCell&3;
		    
		    ilc[ig].x[iv] = rCheck[0];
		    ilc[ig].y[iv] = rCheck[1];
		    ilc[ig].z[iv] = rCheck[2];
		    ilc[ig].m.f[iv] = pkdc->mom.m;
		    ilc[ig].xx.f[iv] = pkdc->mom.xx;
		    ilc[ig].yy.f[iv] = pkdc->mom.yy;
		    ilc[ig].xy.f[iv] = pkdc->mom.xy;
		    ilc[ig].xz.f[iv] = pkdc->mom.xz;
		    ilc[ig].yz.f[iv] = pkdc->mom.yz;
		    ilc[ig].xxx.f[iv] = pkdc->mom.xxx;
		    ilc[ig].xyy.f[iv] = pkdc->mom.xyy;
		    ilc[ig].xxy.f[iv] = pkdc->mom.xxy;
		    ilc[ig].yyy.f[iv] = pkdc->mom.yyy;
		    ilc[ig].xxz.f[iv] = pkdc->mom.xxz;
		    ilc[ig].yyz.f[iv] = pkdc->mom.yyz;
		    ilc[ig].xyz.f[iv] = pkdc->mom.xyz;
		    ilc[ig].xxxx.f[iv] = pkdc->mom.xxxx;
		    ilc[ig].xyyy.f[iv] = pkdc->mom.xyyy;
		    ilc[ig].xxxy.f[iv] = pkdc->mom.xxxy;
		    ilc[ig].yyyy.f[iv] = pkdc->mom.yyyy;
		    ilc[ig].xxxz.f[iv] = pkdc->mom.xxxz;
		    ilc[ig].yyyz.f[iv] = pkdc->mom.yyyz;
		    ilc[ig].xxyy.f[iv] = pkdc->mom.xxyy;
		    ilc[ig].xxyz.f[iv] = pkdc->mom.xxyz;
		    ilc[ig].xyyz.f[iv] = pkdc->mom.xyyz;
#else
		    ilc[nCell].x = rCheck[0];
		    ilc[nCell].y = rCheck[1];
		    ilc[nCell].z = rCheck[2];
		    ilc[nCell].mom = pkdc->mom;
#endif
		    ++nCell;
		    }
		if (id >= 0 && id != pkd->idSelf) {
		    mdlRelease(pkd->mdl,CID_CELL,pkdc);
		    }
		}
	    nCheck = ii;
	    /*
	    ** Evaluate the GLAM list here.
	    */
#ifdef USE_SIMD_LOCR
	    *pdFlop += momGenLocrAddSIMDMomr(&L,nGlam,ilglam,0,0,0,0,0);
#ifdef GLAM_STATS
	    nAvgGlam += nGlam;
	    nTotGlam++;
#endif
	    nGlam = 0;
#endif
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (!c[iCell].iLower) break;
	    xParent = c[iCell].r[0];
	    yParent = c[iCell].r[1];
	    zParent = c[iCell].r[2];
	    iCell = c[iCell].iLower;
	    /*
	    ** Make sure all the check lists are long enough to handle 1 more cell.
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
	    if (c[iCell].iActive & TYPE_ACTIVE) {
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
		if (c[iCell+1].iActive & TYPE_ACTIVE) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.		    */
		    ++iStack;
		    S[iStack].nPart = nPart;
		    S[iStack].nCell = nCell;
/*
		    S[iStack].nPartBucket = nPartBucket;
*/
		    S[iStack].nCheck = nCheck;
		    for (i=0;i<nCheck;++i) S[iStack].Check[i] = Check[i];
		    S[iStack].Check[nCheck-1].iCell = iCell;
		    S[iStack].L = L;
		    momShiftLocr(&S[iStack].L,
				 c[iCell+1].r[0] - xParent,
				 c[iCell+1].r[1] - yParent,
				 c[iCell+1].r[2] - zParent);
#ifdef TIME_WALK_WORK
		    S[iStack].fWeight = getTimer(&tv);
#else
		    S[iStack].fWeight = (*pdFlop-tempI);
#endif
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
	    momShiftLocr(&L,c[iCell].r[0] - xParent,
			 c[iCell].r[1] - yParent,
			 c[iCell].r[2] - zParent);
	    }
	/*
	** Now the interaction list should be complete and the 
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
#ifdef USE_SIMD_LOCR
	assert(nGlam == 0);
#endif
	assert(nCheck == 0);
	/*
	** Add Bucket self interactions.
	*/
	pkdc = &c[iCell];
	/*
	** Local Self-Bucket Interaction.
	** Interact += Pacticles(pkdc);
	*/
	n = pkdc->pUpper - pkdc->pLower + 1;
	if (nPart + n > nMaxPart) {
	    nMaxPart += 500 + n;
	    ilp = realloc(ilp,nMaxPart*sizeof(ILP));
	    assert(ilp != NULL);	
	}
	for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
#ifndef USE_SIMD
	    ilp[nPart].iOrder = p[pj].iOrder;
#endif
	    ilp[nPart].m = p[pj].fMass;
	    /*
	    ** We will assume that all the particles in my bucket are at the same time here so 
	    ** we will not have a drift factor to worry about.
	    */
	    ilp[nPart].x = p[pj].r[0];
	    ilp[nPart].y = p[pj].r[1];
	    ilp[nPart].z = p[pj].r[2];
#ifndef USE_SIMD
	    ilp[nPart].vx = p[pj].v[0]; 
	    ilp[nPart].vy = p[pj].v[1];
	    ilp[nPart].vz = p[pj].v[2];
#endif
	    ilp[nPart].fourh2 = 4*p[pj].fSoft*p[pj].fSoft;
	    ++nPart;
	}

	/*
	** Now calculate gravity on this bucket!
	*/

	nActive = pkdGravInteract(pkd,pkdc,&L,ilp,nPart,ilc,nCell,NULL,0,pdFlop);
	/*
	** Note that if Ewald is being performed we need to factor this
	** constant cost into the load balancing weights.
	*/
	if (bEwald) {
	    *pdFlop += pkdBucketEwald(pkd,&pkd->kdNodes[iCell],nReps,fEwCut,4);
	}

#ifdef TIME_WALK_WORK
	fWeight += getTimer(&tv);
#else
	fWeight += (*pdFlop-tempI);
#endif
	if (nActive) {
#ifdef TIME_WALK_WORK
	    fWeight /= nActive;
#else
	    fWeight = (*pdFlop-tempI)/nActive;
#endif
	    pkdBucketWeight(pkd,iCell,fWeight);
/*
  printf("%6d nPart:%5d nCell:%5d\n",iCell,nPart,nCell);
*/
	    *pdPartSum += nActive*nPart;
	    *pdCellSum += nActive*nCell;
	    nTotActive += nActive;
	    }

	while (iCell & 1) {
	InactiveAscend:
	    iCell = c[iCell].iParent;
	    if (!iCell) {
#if defined(USE_SIMD) && defined(GLAM_STATS)
		printf( "%d: nCalls=%d, AvgGlam=%f\n",
			mdlSelf(pkd->mdl),
			nTotGlam, (float)nAvgGlam/nTotGlam);
#endif
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
#ifdef USE_SIMD
		SIMD_free(ilc);
#else
		free(ilc);
#endif
#ifdef PROFILE_GRAVWALK
    VTPause();
#endif
		return(nTotActive);
		}
	    }
	++iCell;
	if (!(c[iCell].iActive & TYPE_ACTIVE)) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	nPart = S[iStack].nPart;
	nCell = S[iStack].nCell;
/*
	nPartBucket = S[iStack].nPartBucket;
*/
	nCheck = S[iStack].nCheck;
	for (i=0;i<nCheck;++i) Check[i] = S[iStack].Check[i];
	L = S[iStack].L;
#ifdef TIME_WALK_WORK
	fWeight = S[iStack].fWeight;
	clearTimer(&tv);
#else
	fWeight = S[iStack].fWeight;
	tempI = *pdFlop;
#endif
	--iStack;
	}
    }

