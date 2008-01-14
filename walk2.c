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
int pkdGravWalk(PKD pkd,double dTime,int nReps,int bEwald,int bVeryActive,
		double *pdFlop,double *pdPartSum,double *pdCellSum) {
    PARTICLE *p = pkd->pStore;
    PARTICLE *pRemote;
    KDN *c;
    KDN *pkdc;
    LOCR L;
    double dirLsum,normLsum,adotai,maga;
    double tax,tay,taz;
    double fWeight = 0.0;
    double dEwaldFlop;
    double dShiftFlop;
    FLOAT dMin,min2,d2,fourh2;
    FLOAT rCheck[3];
    FLOAT rOffset[3];
    FLOAT xParent,yParent,zParent;
    FLOAT dx[3],dir;
    FLOAT cx,cy,cz,d2c;
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxInitCheck,nCheck;
    int iCell,iSib,iCheckCell,iCellDescend;
    int i,ii,j,n,id,pj,nActive,nTotActive;
    int iOpen;
    ILPTILE tile;
    int nCell;
#ifdef USE_SIMD_MOMR
    int ig,iv;
#endif
#ifdef TIME_WALK_WORK
    TIMER tv;
#else
    double tempM;
    double tempI;
#endif
    double dEwFlop = 0.0;

#ifdef PROFILE_GRAVWALK
    VTResume();
#endif

    /*
    ** If we are doing the very active gravity then check that there is a very active tree!
    ** Otherwise we check that the ROOT has active particles!
    */
    if (bVeryActive) {
	assert(pkd->nVeryActive != 0);
	assert(pkd->nVeryActive == pkd->kdNodes[VAROOT].pUpper - pkd->kdNodes[VAROOT].pLower + 1);
	}
    else if (!pkdIsCellActive(pkd,&pkd->kdNodes[ROOT])) return 0;
    /*
    ** Initially we set our cell pointer to 
    ** point to the top tree.
    */
    c = pkd->kdTop;
    nTotActive = 0;
    ilpClear(pkd->ilp);
    nCell = 0;
    /*
    ** Allocate Checklist.
    */
    nMaxInitCheck = 2*nReps+1;
    nMaxInitCheck = nMaxInitCheck*nMaxInitCheck*nMaxInitCheck;	/* all replicas */
    iCell = pkd->iTopRoot;
    while ((iCell = c[iCell].iParent)) ++nMaxInitCheck; /* all top tree siblings */
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
    ** Initialize the PP interaction list center, just in case it has not been done yet!
    */
    pkd->ilp->cx = 0;
    pkd->ilp->cy = 0;
    pkd->ilp->cz = 0;	    
    pkd->ilp->d2cmax = 0;
    /*
    ** Now switch the cell pointer to point to 
    ** the local tree.
    */
    c = pkd->kdNodes;
    /*
    ** Make iCell point to the root of the tree again.
    */
    if (bVeryActive) iCell = VAROOT;
    else iCell = ROOT;
    while (1) {

	/*
	** Find the next active particle that will be encountered in the walk algorithm below
	** in order to set a good initial center for the P-P interaction list.
	*/
	iCellDescend = iCell;
	while (c[iCellDescend].iLower) {
	    iCellDescend = c[iCellDescend].iLower;
	    if (!pkdIsCellActive(pkd,&c[iCellDescend])) {
		/*
		** Move onto processing the sibling.
		*/
		++iCellDescend;
		}
	    }
	for (pj=c[iCellDescend].pLower;pj<=c[iCellDescend].pUpper;++pj) {
	    if (!pkdIsActive(pkd,&p[pj])) continue;
	    cx = p[pj].r[0];
	    cy = p[pj].r[1];
	    cz = p[pj].r[2];
	    break;
	    }
	assert(pj <= c[iCell].pUpper);  /* otherwise we did not come to an active particle */
	d2c = (cx - pkd->ilp->cx)*(cx - pkd->ilp->cx) + (cy - pkd->ilp->cy)*(cy - pkd->ilp->cy) + 
	    (cz - pkd->ilp->cz)*(cz - pkd->ilp->cz);
	//if (d2c > pkd->ilp->d2cmax) {
	if ( d2c > 1e-5) {
//	    printf("%d:Shift of center too large for the coming interactions! old:(%.10g,%.10g,%.10g) new:(%.10g,%.10g,%.10g)\n",
//		   mdlSelf(pkd->mdl),pkd->ilp->cx,pkd->ilp->cy,pkd->ilp->cz,cx,cy,cz);
	    /*
	    ** Correct all remaining PP interactions to this new center.
	    */
	    for( tile=pkd->ilp->first; tile!=pkd->ilp->tile->next; tile=tile->next ) {
		for( j=0; j<tile->nPart; ++j ) {
		    tile->dx.f[j] += cx - pkd->ilp->cx;
		    tile->dy.f[j] += cy - pkd->ilp->cy;
		    tile->dz.f[j] += cz - pkd->ilp->cz;
		}
	    }
	    pkd->ilp->cx = cx;
	    pkd->ilp->cy = cy;
	    pkd->ilp->cz = cz;	    
	    pkd->ilp->d2cmax = c[iCellDescend].fOpen*c[iCellDescend].fOpen;
	}

	while (1) {
	    /*
	    ** Process the Checklist.
	    */

#ifdef TIME_WALK_WORK
	    clearTimer(&tv);
#else
	    tempI = *pdFlop;
	    tempI += dEwFlop;
	    tempM = 0.001*(pkd->mdl->cache[CID_PARTICLE].nMiss + pkd->mdl->cache[CID_CELL].nMiss);
#endif
	    ii = 0;

	    for (i=0;i<nCheck;++i) {
		id = pkd->Check[i].id;

#ifdef __SSE__strange
		    /* Does nothing */
		    _mm_prefetch( (char *)(c+iCell)+0,_MM_HINT_T0 );
		    _mm_prefetch( (char *)(c+iCell)+64,_MM_HINT_T0 );
		    _mm_prefetch( (char *)(c+iCell)+128,_MM_HINT_T0 );
#endif
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
#ifdef TIME_WALK_WORK
		    stopTimer(&tv);
#endif
		    pkdc = mdlAquire(pkd->mdl,CID_CELL,pkd->Check[i].iCell,id);
#ifdef TIME_WALK_WORK
		    startTimer(&tv);
#endif
		    n = pkdc->pUpper - pkdc->pLower + 1;
		    }
		for (j=0;j<3;++j) rCheck[j] = pkdc->r[j] + pkd->Check[i].rOffset[j];
		d2 = 0;
		for (j=0;j<3;++j) {
		    dx[j] = c[iCell].r[j] - rCheck[j];
		    d2 += dx[j]*dx[j];
		}
		fourh2 = softmassweight(c[iCell].mom.m,4*c[iCell].fSoft2,pkdc->mom.m,4*pkdc->fSoft2);
		iOpen = 0;
		if (d2 > (c[iCell].fOpen + pkdc->fOpen)*(c[iCell].fOpen + pkdc->fOpen)) {
		    /*
		    ** Accept local expansion, but check softening.
		    */
		    if (d2 > fourh2) iOpen = -1;
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
			if (min2 > pkdc->fOpen*dMonopoleThetaFac*pkdc->fOpen*dMonopoleThetaFac) iOpen = -3;
			/*
			** Otherwise we must open one of the two cells. We defer this decision to the tests
			** which follow. In this case we still have iOpen  = 0.
			*/
		    }
		}
		if (iOpen == 0) {
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
			    ** We can't open pkdc (it is a bucket).  At this point we would prefer to accept
			    ** the bucket as a multipole expansion, but if that isn't possible use its particles.
			    */
			    if (n >= WALK_MINMULTIPOLE) {
				min2 = 0;	
				for (j=0;j<3;++j) {
				    dMin = fabs(rCheck[j] - c[iCell].bnd.fCenter[j]);
				    dMin -= c[iCell].bnd.fMax[j];
				    if (dMin > 0) min2 += dMin*dMin;
				}
				if (min2 > pkdc->fOpen*pkdc->fOpen) {
				    if (min2 > fourh2) {
					/*
					** The multipole is also ok as far as softening goes, so accept it.
					*/
					iOpen = -2;
				    }
				    else if (min2 > pkdc->fOpen*dMonopoleThetaFac*pkdc->fOpen*dMonopoleThetaFac) {
					iOpen = -3;
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
				** This bucket did not have enough particles to make it worth accepting as a
				** multipole, since it is faster to simply add P-P interactions at this stage.
				*/
				iOpen = 1;
			    }
			}
		    }
		    /*
		    ** From here we know that this c[iCell] has the larger opening radius.
		    */
		    else if (!c[iCell].iLower) {
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
				if (min2 > fourh2) {
				    /*
				    ** The multipole is also ok as far as softening goes, so accept it.
				    */
				    iOpen = -2;
				}
				else if (min2 > pkdc->fOpen*dMonopoleThetaFac*pkdc->fOpen*dMonopoleThetaFac) {
				    iOpen = -3;
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
		if (iOpen == 0) {
		    pkd->Check[ii++] = pkd->Check[i];
		    }
		else if (iOpen > 0) {
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
			if (nCheck + 2 > pkd->nMaxCheck) {
			    pkd->nMaxCheck += 1000;
			    pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->Check != NULL);
			    for (ism=0;ism<pkd->nMaxStack;++ism) {
				pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
				assert(pkd->S[ism].Check != NULL);
				}
			    printf("A CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
			    }
			pkd->Check[nCheck] = pkd->Check[i];
			pkd->Check[nCheck+1] = pkd->Check[i];
			/*
			** If we are opening a leaf of the top tree
			** we need to correctly set the processor id.
			** (this is a bit tricky)
			*/
			if (id < 0) {
			    if (pkd->kdTop[iCheckCell].pLower >= 0) {
				pkd->Check[nCheck].iCell = ROOT;
				pkd->Check[nCheck].id = pkd->kdTop[iCheckCell].pLower;
				}
			    else {
				pkd->Check[nCheck].iCell = iCheckCell;
				assert(pkd->Check[nCheck].id == -1);
				}
			    if (pkd->kdTop[iCheckCell+1].pLower >= 0) {
				pkd->Check[nCheck+1].iCell = ROOT;
				pkd->Check[nCheck+1].id = pkd->kdTop[iCheckCell+1].pLower;
				}
			    else {
				pkd->Check[nCheck+1].iCell = iCheckCell+1;
				assert(pkd->Check[nCheck+1].id == -1);
				}

			    }
			else {
			    pkd->Check[nCheck].iCell = iCheckCell;
			    pkd->Check[nCheck+1].iCell = iCheckCell+1;
			    assert(pkd->Check[nCheck].id == id);
			    assert(pkd->Check[nCheck+1].id == id);
			    }
			nCheck += 2;
			}
		    else {
			/*
			** Now I am trying to open a bucket, which means I place particles on the ilp
			** interaction list.
			*/
			tile = pkd->ilp->tile;
			if (id == pkd->idSelf) {
			    /*
			    ** Local Bucket Interaction.
			    ** Interact += Pacticles(pkdc);
			    */
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				if ( tile->nPart == tile->nMaxPart )
				    tile = ilpExtend(pkd->ilp);
				tile->dx.f[tile->nPart] = pkd->ilp->cx - (p[pj].r[0] + pkd->Check[i].rOffset[0]);
				tile->dy.f[tile->nPart] = pkd->ilp->cy - (p[pj].r[1] + pkd->Check[i].rOffset[1]);
				tile->dz.f[tile->nPart] = pkd->ilp->cz - (p[pj].r[2] + pkd->Check[i].rOffset[2]);
				tile->m.f[tile->nPart] = p[pj].fMass;
				tile->fourh2.f[tile->nPart] = 4*p[pj].fSoft*p[pj].fSoft;
#ifdef HERMITE
				tile->vx.f[tile->nPart] = p[pj].v[0];
				tile->vy.f[tile->nPart] = p[pj].v[1];
				tile->vz.f[tile->nPart] = p[pj].v[2];
#endif
#if defined(SYMBA) || defined(PLANETS)
				tile->iOrder.i[tile->nPart] = p[pj].iOrder;
#endif
				++tile->nPart;
			    }
			}
			else {
			    /*
			    ** Remote Bucket Interaction.
			    ** Interact += Pacticles(pkdc);
			    */
			    //ADDBACK:printf("C CPU:%d increased particle list size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxPart);
#ifdef TIME_WALK_WORK
			    stopTimer(&tv);
#endif
			    for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
				pRemote = mdlAquire(pkd->mdl,CID_PARTICLE,pj,id);
				if ( tile->nPart == tile->nMaxPart )
				    tile = ilpExtend(pkd->ilp);

				tile->dx.f[tile->nPart] = pkd->ilp->cx - (pRemote->r[0] + pkd->Check[i].rOffset[0]);
				tile->dy.f[tile->nPart] = pkd->ilp->cy - (pRemote->r[1] + pkd->Check[i].rOffset[1]);
				tile->dz.f[tile->nPart] = pkd->ilp->cz - (pRemote->r[2] + pkd->Check[i].rOffset[2]);
				tile->m.f[tile->nPart] = pRemote->fMass;
				tile->fourh2.f[tile->nPart] = 4*pRemote->fSoft*pRemote->fSoft;
#if defined(SYMBA) || defined(PLANETS)
				tile->iOrder.i[tile->nPart] = pRemote->iOrder;
#endif
#ifdef HERMITE
				tile->vx.f[tile->nPart] = pRemote->v[0]; 
				tile->vy.f[tile->nPart] = pRemote->v[1]; 
				tile->vz.f[tile->nPart] = pRemote->v[2]; 
#endif
				++tile->nPart;
				mdlRelease(pkd->mdl,CID_PARTICLE,pRemote);
			    }
#ifdef TIME_WALK_WORK
			    startTimer(&tv);
#endif
			}  /* end of opening a bucket */
		    }
		}
		else if (iOpen == -1) {
		    /*
		    ** Local expansion accepted!
		    ** Add to the GLAM list to be evaluated later.
		    */
		    dir = 1.0/sqrt(d2);
		    *pdFlop += momLocrAddMomr5(&L,&pkdc->mom,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
		    adotai = pkdc->a[0]*(-tax) + pkdc->a[1]*(-tay) + pkdc->a[2]*(-taz); /* temporary hack to get it right */
		    if (adotai > 0) {
			maga = sqrt(pkdc->a[0]*pkdc->a[0] + pkdc->a[1]*pkdc->a[1] + pkdc->a[2]*pkdc->a[2]);
			adotai /= maga;
			dirLsum += dir*adotai*adotai;
			normLsum += adotai*adotai;
			}
		}
		else if (iOpen == -2) {
		    /*
		    ** No intersection, accept multipole!
		    ** Interact += Moment(pkdc);
		    */
		    if (nCell == pkd->nMaxCell) {
			pkd->nMaxCell += 500;
#ifdef USE_SIMD_MOMR
			pkd->ilc = SIMD_realloc(pkd->ilc, (pkd->nMaxCell-500)/4*sizeof(ILC),
					   pkd->nMaxCell/4*sizeof(ILC));
#else
			pkd->ilc = realloc(pkd->ilc,pkd->nMaxCell*sizeof(ILC));
#endif
			assert(pkd->ilc != NULL);
			printf("D CPU:%d increased cell list size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCell);
			}
#ifdef USE_SIMD_MOMR
		    ig = nCell >> 2;
		    iv = nCell&3;
		    
		    pkd->ilc[ig].x[iv] = rCheck[0];
		    pkd->ilc[ig].y[iv] = rCheck[1];
		    pkd->ilc[ig].z[iv] = rCheck[2];
		    pkd->ilc[ig].m.f[iv] = pkdc->mom.m;
		    pkd->ilc[ig].xx.f[iv] = pkdc->mom.xx;
		    pkd->ilc[ig].yy.f[iv] = pkdc->mom.yy;
		    pkd->ilc[ig].xy.f[iv] = pkdc->mom.xy;
		    pkd->ilc[ig].xz.f[iv] = pkdc->mom.xz;
		    pkd->ilc[ig].yz.f[iv] = pkdc->mom.yz;
		    pkd->ilc[ig].xxx.f[iv] = pkdc->mom.xxx;
		    pkd->ilc[ig].xyy.f[iv] = pkdc->mom.xyy;
		    pkd->ilc[ig].xxy.f[iv] = pkdc->mom.xxy;
		    pkd->ilc[ig].yyy.f[iv] = pkdc->mom.yyy;
		    pkd->ilc[ig].xxz.f[iv] = pkdc->mom.xxz;
		    pkd->ilc[ig].yyz.f[iv] = pkdc->mom.yyz;
		    pkd->ilc[ig].xyz.f[iv] = pkdc->mom.xyz;
		    pkd->ilc[ig].xxxx.f[iv] = pkdc->mom.xxxx;
		    pkd->ilc[ig].xyyy.f[iv] = pkdc->mom.xyyy;
		    pkd->ilc[ig].xxxy.f[iv] = pkdc->mom.xxxy;
		    pkd->ilc[ig].yyyy.f[iv] = pkdc->mom.yyyy;
		    pkd->ilc[ig].xxxz.f[iv] = pkdc->mom.xxxz;
		    pkd->ilc[ig].yyyz.f[iv] = pkdc->mom.yyyz;
		    pkd->ilc[ig].xxyy.f[iv] = pkdc->mom.xxyy;
		    pkd->ilc[ig].xxyz.f[iv] = pkdc->mom.xxyz;
		    pkd->ilc[ig].xyyz.f[iv] = pkdc->mom.xyyz;
#else
		    pkd->ilc[nCell].x = rCheck[0];
		    pkd->ilc[nCell].y = rCheck[1];
		    pkd->ilc[nCell].z = rCheck[2];
		    pkd->ilc[nCell].mom = pkdc->mom;
#endif
		    ++nCell;
		    }
		else if (iOpen == -3) {
		    /*
		    ** We accept this multipole from the opening criterion, but it is a softened
		    ** interaction, so we need to treat is as a softened monopole by putting it
		    ** on the particle interaction list.
		    */
			//ADDBACK: printf("E CPU:%d increased particle list size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxPart);

		    tile = pkd->ilp->tile;
		    if ( tile->nPart == tile->nMaxPart )
			tile = ilpExtend(pkd->ilp);
		    tile->m.f[tile->nPart] = pkdc->mom.m;
		    tile->dx.f[tile->nPart] = pkd->ilp->cx - rCheck[0];
		    tile->dy.f[tile->nPart] = pkd->ilp->cy - rCheck[1];
		    tile->dz.f[tile->nPart] = pkd->ilp->cz - rCheck[2];
		    tile->fourh2.f[tile->nPart] = 4*pkdc->fSoft2;
#if defined(SYMBA) || defined(PLANETS)
		    tile->iOrder.i[tile->nPart] = -1; /* set iOrder to negative value for time step criterion */
#endif
#ifdef HERMITE
		    tile->vx.f[tile->nPart] = pkdc->v[0];
		    tile->vy.f[tile->nPart] = pkdc->v[1];
		    tile->vz.f[tile->nPart] = pkdc->v[2];
#endif
		    ++tile->nPart;
		}
		else {
		    mdlassert(pkd->mdl,iOpen >= -3 && iOpen <= 1);
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
	    xParent = c[iCell].r[0];
	    yParent = c[iCell].r[1];
	    zParent = c[iCell].r[2];
	    iCell = c[iCell].iLower;
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
	    if (pkdIsCellActive(pkd,&c[iCell])) {
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
		if (pkdIsCellActive(pkd,&c[iCell+1])) {
		    /*
		    ** Sibling is active as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.		    
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    ilpCheckPt(pkd->ilp,&pkd->S[iStack].PartChkPt);
		    pkd->S[iStack].nCell = nCell;
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
					      c[iCell+1].r[0] - xParent,
					      c[iCell+1].r[1] - yParent,
					      c[iCell+1].r[2] - zParent);
#ifdef TIME_WALK_WORK
		    pkd->S[iStack].fWeight = getTimer(&tv);
#elif defined(COUNT_MISSES)
		    pkd->S[iStack].fWeight = 0.001*(pkd->mdl->cache[CID_PARTICLE].nMiss + pkd->mdl->cache[CID_CELL].nMiss) - tempM;
#else
		    pkd->S[iStack].fWeight = (*pdFlop-tempI) + dShiftFlop;
		    pkd->S[iStack].fWeight += dEwFlop;
#endif
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
	    *pdFlop += momShiftLocr(&L,c[iCell].r[0] - xParent,
				    c[iCell].r[1] - yParent,
				    c[iCell].r[2] - zParent);
	    }
	/*
	** Now the interaction list should be complete and the 
	** Checklist should be empty! Calculate gravity on this
	** Bucket!
	*/
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
	//ADDBACK:printf("G CPU:%d increased particle list size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxPart);
	for (pj=pkdc->pLower;pj<=pkdc->pUpper;++pj) {
	    tile = pkd->ilp->tile;
	    if ( tile->nPart == tile->nMaxPart )
		tile = ilpExtend(pkd->ilp);
	    tile->m.f[tile->nPart] = p[pj].fMass;
	    /*
	    ** We will assume that all the particles in my bucket are at the same time here so 
	    ** we will not have a drift factor to worry about.
	    */
	    tile->dx.f[tile->nPart] = pkd->ilp->cx - p[pj].r[0];
	    tile->dy.f[tile->nPart] = pkd->ilp->cy - p[pj].r[1];
	    tile->dz.f[tile->nPart] = pkd->ilp->cz - p[pj].r[2];
	    tile->fourh2.f[tile->nPart] = 4*p[pj].fSoft*p[pj].fSoft;
#if defined(SYMBA) || defined(PLANETS)
	    tile->iOrder.i[tile->nPart] = p[pj].iOrder;
#endif
#ifdef HERMITE
	    tile->vx.f[tile->nPart] = p[pj].v[0];
	    tile->vy.f[tile->nPart] = p[pj].v[1];
	    tile->vz.f[tile->nPart] = p[pj].v[2];
#endif
	    ++tile->nPart;
	}

	/*
	** Now calculate gravity on this bucket!
	*/
	nActive = pkdGravInteract(pkd,pkdc,&L,pkd->ilp,pkd->ilc,nCell,dirLsum,normLsum,
				  bEwald,pdFlop,&dEwFlop);
	/*
	** Update the limit for a shift of the center here based on the opening radius of this 
	** cell (the one we just evaluated).
	*/
	pkd->ilp->d2cmax = c[iCell].fOpen*c[iCell].fOpen;

#ifdef TIME_WALK_WORK
	fWeight += getTimer(&tv);
#elif defined(COUNT_MISSES)
	fWeight += 0.001*(pkd->mdl->cache[CID_PARTICLE].nMiss + pkd->mdl->cache[CID_CELL].nMiss) - tempM;
#else
	fWeight += (*pdFlop-tempI);
	fWeight += dEwFlop;
#endif
	if (nActive) {
	    fWeight /= nActive;
#if !defined(COUNT_MISSES) && !defined(TIME_WALK_WORK)
	    /*
	    ** The simplest thing to do is to set the weight for all particles to 1.0 which 
	    ** seems to work the best.
	    */
	    fWeight = 1.0;
#endif
	    pkdBucketWeight(pkd,iCell,/*fWeight*/1.0);
/*
  printf("%6d nPart:%5d nCell:%5d\n",iCell,nPart,nCell);
*/
	    *pdPartSum += nActive*ilpCount(pkd->ilp);
	    *pdCellSum += nActive*nCell;
	    nTotActive += nActive;
	    }

	while (iCell & 1) {
	InactiveAscend:
	    iCell = c[iCell].iParent;
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
	++iCell;
	if (!pkdIsCellActive(pkd,&c[iCell])) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	ilpRestore(pkd->ilp,&pkd->S[iStack].PartChkPt);
	nCell = pkd->S[iStack].nCell;
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
#ifdef TIME_WALK_WORK
	fWeight = pkd->S[iStack].fWeight;
	clearTimer(&tv);
#else
	fWeight = pkd->S[iStack].fWeight;
	tempI = *pdFlop;
	tempI += dEwFlop;
	tempM = 0.001*(pkd->mdl->cache[CID_PARTICLE].nMiss + pkd->mdl->cache[CID_CELL].nMiss);
#endif
	--iStack;
	}
    }

