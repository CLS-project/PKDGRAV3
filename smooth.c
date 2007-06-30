#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "smoothfcn.h"

int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bGasOnly,
		 int bPeriodic,int bSymmetric,int iSmoothType, 
		 int eParticleTypes,
		 double dfBall2OverSoft2 ) {
    SMX smx;
    void (*initParticle)(void *) = NULL;
    void (*init)(void *) = NULL;
    void (*comb)(void *,void *) = NULL;
    int pi;
    int nTree;
    int iTopDepth;

    smx = malloc(sizeof(struct smContext));
    assert(smx != NULL);
    smx->pkd = pkd;
    if (smf != NULL) smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->bGasOnly = bGasOnly;
    smx->bPeriodic = bPeriodic;
    smx->eParticleTypes = eParticleTypes;

    switch (iSmoothType) {
    case SMX_NULL:
	smx->fcnSmooth = NullSmooth;
	initParticle = NULL; /* Original Particle */
	init = NULL; /* Cached copies */
	comb = NULL;
	smx->fcnPost = NULL;
	break;
    case SMX_DENSITY:
	smx->fcnSmooth = bSymmetric?DensitySym:Density;
	initParticle = initDensity; /* Original Particle */
	init = initDensity; /* Cached copies */
	comb = combDensity;
	smx->fcnPost = NULL;
	break;
    case SMX_MARKDENSITY:
	smx->fcnSmooth = bSymmetric?MarkDensitySym:MarkDensity;
	initParticle = initParticleMarkDensity; /* Original Particle */
	init = initMarkDensity; /* Cached copies */
	comb = combMarkDensity;
	smx->fcnPost = NULL;
	break;
    case SMX_MARKIIDENSITY:
	smx->fcnSmooth = bSymmetric?MarkIIDensitySym:MarkIIDensity;
	initParticle = initParticleMarkIIDensity; /* Original Particle */
	init = initMarkIIDensity; /* Cached copies */
	comb = combMarkIIDensity;
	smx->fcnPost = NULL;
	break;
    case SMX_MARK:
	smx->fcnSmooth = NULL;
	initParticle = NULL;
	init = initMark;
	comb = combMark;
	smx->fcnPost = NULL;
	break;
    case SMX_FOF:
	assert(bSymmetric == 0);
	smx->fcnSmooth = NULL;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
#ifdef RELAXATION
    case SMX_RELAXATION:
	assert(bSymmetric == 0);
	smx->fcnSmooth = AddRelaxation;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
#endif /* RELAXATION */
#ifdef SYMBA
    case SMX_SYMBA:
	assert(bSymmetric == 0);
	smx->fcnSmooth = DrmininDrift;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
#endif /* SYMBA */

    default:
	assert(0);
	}
    /*
    ** Initialize the ACTIVE particles in the tree.
    ** There are other particles in the tree -- just not active.
    */
    nTree = pkd->kdNodes[ROOT].pUpper + 1;
    if (initParticle != NULL) {
	for (pi=0;pi<nTree;++pi) {
	    if (TYPETest(&(pkd->pStore[pi]),smx->eParticleTypes)) {
		initParticle(&pkd->pStore[pi]);
	    }
	}
    }
    /*
    ** Start particle caching space (cell cache is already active).
    */
    if (bSymmetric) {
	mdlCOcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
		   nTree,init,comb);
	}
    else {
	mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
		   nTree);
	}
    /*
    ** Allocate Nearest-Neighbor List.
    ** And also the associated list bRemote flags.
    */
    smx->nnListSize = 0;
    smx->nnListMax = NNLIST_INCREMENT;
    smx->nnList = malloc(smx->nnListMax*sizeof(NN));
    assert(smx->nnList != NULL);
    smx->nnbRemote = malloc(smx->nnListMax*sizeof(int));
    assert(smx->nnbRemote != NULL);
    /*
    ** Allocate priority queue.
    */
    smx->pq = malloc(nSmooth*sizeof(PQ));
    assert(smx->pq != NULL);
    PQ_INIT(smx->pq,nSmooth);
    /*
    ** Allocate special stacks for searching.
    ** There is a mistake here, since I use these stacks for the remote trees as well.
    ** This can be easily fixed, but a hack for now.
    */
    smx->S = malloc(1024*sizeof(int));
    assert(smx->S != NULL);
    smx->Smin = malloc(1024*sizeof(FLOAT));
    assert(smx->Smin != NULL);
    /*
    ** Allocate special stacks for searching within the top tree.
    ** Calculate the number of levels in the top tree. 
    */
    iTopDepth = 1+(int)ceil(log((double)smx->pkd->nThreads)/log(2.0));
    smx->ST = malloc(iTopDepth*sizeof(int));
    assert(smx->ST != NULL);
    smx->SminT = malloc(iTopDepth*sizeof(FLOAT));
    assert(smx->SminT != NULL);
    *psmx = smx;	
    return(1);
    }


void smFinish(SMX smx,SMF *smf)
    {
    PKD pkd = smx->pkd;
    int pi;
    char achOut[128];

    /*
     * Output statistics.
     */
    sprintf(achOut, "Cell Accesses: %g\n",
	    mdlNumAccess(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Min ratio: %g\n",
	    mdlMinRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "Particle Accesses: %g\n",
	    mdlNumAccess(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Min ratio: %g\n",
	    mdlMinRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    /*
    ** Stop particle caching space.
    */
    mdlFinishCache(smx->pkd->mdl,CID_PARTICLE);
    /*
    ** Now do any post calculations, these ususlly involve some sort of
    ** normalizations of the smoothed quantities, usually division by
    ** the local density! Do NOT put kernel normalizations in here as
    ** these do not depend purely on local properties in the case of
    ** "Gather-Scatter" kernel.
    */
    if (smx->fcnPost != NULL) {
	for (pi=0;pi<pkd->nLocal;++pi) {
	    if (TYPETest(&(pkd->pStore[pi]),smx->eParticleTypes)) {
		smx->fcnPost(&pkd->pStore[pi],smf);
		}
	    }
	}
    /*
    ** Free up context storage.
    */
    free(smx->S);
    free(smx->Smin);
    free(smx->ST);
    free(smx->SminT);
    free(smx->pq);
    free(smx->nnList);
    free(smx->nnbRemote);
    free(smx);
    }


/*
** This function performs a local nearest neighbor search.
** Note that this function cannot be called for a periodic
** replica of the local domain, that can be done with the
** pqSearchRemote function setting id == idSelf.
*/
PQ *pqSearchLocal(SMX smx,PARTICLE *pi,FLOAT r[3],int *pbDone) 
    {
    PARTICLE *p = smx->pkd->pStore;
    KDN *c = smx->pkd->kdNodes;
    PQ *pq;
    FLOAT dx,dy,dz,dMin,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int i,j,n,pj,pWant,pEnd,iCell,iParent;
    int sp = 0;
    int sm = 0;

    *pbDone = 1;	/* assume that we will complete the search */
    assert(smx->nQueue == 0);
    pq = smx->pq;
    /*
    ** Decide where the first containment test needs to
    ** be performed. If no particle is specfied then we
    ** don't perform containment tests except at the 
    ** root, so that the pbDone flag can be correctly
    ** set.
    */
    if (pi) {
	iCell = pi->iBucket;
	while (iCell != ROOT) {
#ifdef GASOLINE
	    if (smx->bGasOnly) n = c[iCell].nGas;
	    else n = c[iCell].pUpper - c[iCell].pLower + 1;
#else
	    n = c[iCell].pUpper - c[iCell].pLower + 1;
#endif
	    if (n < smx->nSmooth) iCell = c[iCell].iParent;
	    else break;
	    }
	S[sp] = iCell;
	iCell = pi->iBucket;
	}
    else {
	iCell = ROOT;
	S[sp] = iCell;
	}
    /*
    ** Start of PRIOQ Loading loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
	    MINDIST(c[iCell].bnd,r,min1);
	    ++iCell;
	    MINDIST(c[iCell].bnd,r,min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		--iCell;
		}
	    else {
		Smin[sm++] = min1;
		}
#ifdef GASOLINE
	    if (smx->bGasOnly && c[iCell].nGas == 0) goto LoadNotContained;
#endif
	    }
	pWant = c[iCell].pLower + smx->nSmooth - smx->nQueue - 1;
#ifdef GASOLINE
	if (smx->bGasOnly) pEnd = c[iCell].pLower + c[iCell].nGas - 1;
	else pEnd = c[iCell].pUpper;
#else
	pEnd = c[iCell].pUpper;
#endif
	if (pWant > pEnd) {
	    for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
		dx = r[0] - p[pj].r[0];
		dy = r[1] - p[pj].r[1];
		dz = r[2] - p[pj].r[2];
		pq[smx->nQueue].pPart = &p[pj];
		pq[smx->nQueue].fDist2 = dx*dx + dy*dy + dz*dz;
		++smx->nQueue;
		}
	    }
	else {
	    for (pj=c[iCell].pLower;pj<=pWant;++pj) {
		dx = r[0] - p[pj].r[0];
		dy = r[1] - p[pj].r[1];
		dz = r[2] - p[pj].r[2];
		pq[smx->nQueue].pPart = &p[pj];
		pq[smx->nQueue].fDist2 = dx*dx + dy*dy + dz*dz;
		++smx->nQueue;
		}
	    PQ_BUILD(pq,smx->nSmooth,pq);
	    for (;pj<=pEnd;++pj) {
		dx = r[0] - p[pj].r[0];
		dy = r[1] - p[pj].r[1];
		dz = r[2] - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < pq->fDist2) {
		    pq->pPart = &p[pj];
		    pq->fDist2 = fDist2;
		    PQ_REPLACE(pq);
		    }
		}
	    goto NoIntersect;  /* done loading phase */
	    }
#ifdef GASOLINE
    LoadNoIntersect:
#endif
	while (iCell == S[sp]) {
	    if (!sp) {
		*pbDone = 0;
		/*
		** Set dx,dy and dz and bRemote before leaving this
		** function.
		*/
		for (i=0;i<smx->nQueue;++i) {
		    pq[i].dx = r[0] - pq[i].pPart->r[0];
		    pq[i].dy = r[1] - pq[i].pPart->r[1];
		    pq[i].dz = r[2] - pq[i].pPart->r[2];
		    pq[i].bRemote = 0;
		    }
		return NULL;		/* EXIT, could not load enough particles! */
		}
	    --sp;
	    iCell = c[iCell].iParent;
	    }
#ifdef GASOLINE
    LoadNotContained:
#endif
	iCell ^= 1;
	if (sm) --sm;
#ifdef GASOLINE
	if (smx->bGasOnly && c[iCell].nGas == 0) {
	    iCell = c[iCell].iParent;
	    goto LoadNoIntersect;   			
	    }
#endif
	S[++sp] = iCell;
	}
    /*
    ** Start of PRIOQ searching loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
#ifdef GASOLINE
	    if (smx->bGasOnly) {
		if (c[iCell].nGas == 0) min1 = pq->fDist2;
		else {
		    MINDIST(c[iCell].bnd,r,min1);
		    }
		++iCell;
		if (c[iCell].nGas == 0) min2 = pq->fDist2;
		else {
		    MINDIST(c[iCell].bnd,r,min2);
		    }
		}
	    else {
		MINDIST(c[iCell].bnd,r,min1);
		++iCell;
		MINDIST(c[iCell].bnd,r,min2);
		}
#else
	    MINDIST(c[iCell].bnd,r,min1);
	    ++iCell;
	    MINDIST(c[iCell].bnd,r,min2);
#endif
	    if (min1 < min2) {
		Smin[sm++] = min2;
		--iCell;
		if (min1 >= pq->fDist2) goto NotContained;
		}
	    else {
		Smin[sm++] = min1;
		if (min2 >= pq->fDist2) goto NotContained;
		}
	    }
#ifdef GASOLINE
	if (smx->bGasOnly) pEnd = c[iCell].pLower + c[iCell].nGas - 1;
	else pEnd = c[iCell].pUpper;
#else
	pEnd = c[iCell].pUpper;
#endif
	for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
	    dx = r[0] - p[pj].r[0];
	    dy = r[1] - p[pj].r[1];
	    dz = r[2] - p[pj].r[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < pq->fDist2) {
		pq->pPart = &p[pj];
		pq->fDist2 = fDist2;
		PQ_REPLACE(pq);
		}
	    }
    NoIntersect:
	while (iCell == S[sp]) {
	    if (sp) {
		--sp;
		iCell = c[iCell].iParent;
		}
	    else {
		/*
		** Containment Test! 
		*/
		for (j=0;j<3;++j) {
		    dMin = c[iCell].bnd.fMax[j] - 
			fabs(c[iCell].bnd.fCenter[j] - r[j]);
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = c[iCell].iParent;
			if (!iParent) {
			    *pbDone = 0;		/* EXIT, not contained! */
			    break;
			    }
			S[sp] = iParent;
			goto NotContained;
			}
		    }
		/*
		** Set dx,dy and dz and bRemote before leaving this
		** function.
		*/
		for (i=0;i<smx->nQueue;++i) {
		    smx->pq[i].dx = r[0] - smx->pq[i].pPart->r[0];
		    smx->pq[i].dy = r[1] - smx->pq[i].pPart->r[1];
		    smx->pq[i].dz = r[2] - smx->pq[i].pPart->r[2];
		    smx->pq[i].bRemote = 0;
		    }
		return pq;
		}
	    }
    NotContained:
	iCell ^= 1;		
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
#ifdef GASOLINE
	    if (smx->bGasOnly && c[iCell].nGas == 0) {
		iCell = c[iCell].iParent;
		goto NoIntersect;
		}
#endif
	    MINDIST(c[iCell].bnd,r,min2);
	    }
	if (min2 >= pq->fDist2) {
	    iCell = c[iCell].iParent;
	    goto NoIntersect;
	    }
	S[++sp] = iCell;
	}
    }



PQ *pqSearchRemote(SMX smx,PQ *pq,int id,FLOAT r[3]) 
    {
    MDL mdl = smx->pkd->mdl;
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p;
    KDN *pkdn,*pkdu;
    FLOAT dx,dy,dz,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int pj,pWant,pEnd,iCell;
    int sp = 0;
    int sm = 0;
    int idSelf = smx->pkd->idSelf;

    iCell = ROOT;
    S[sp] = iCell;
    if (id == idSelf) pkdn = &c[iCell];
    else pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    if (smx->nQueue == smx->nSmooth) goto StartSearch;
    /*
    ** Start of PRIOQ Loading loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (pkdn->iLower) {
	    iCell = pkdn->iLower;
	    if (id == idSelf) pkdn = &c[iCell];
	    else {
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		}
	    MINDIST(pkdn->bnd,r,min1);
	    ++iCell;
	    if (id == idSelf) pkdu = &c[iCell];
	    else pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
	    MINDIST(pkdu->bnd,r,min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		if (id != idSelf) mdlRelease(mdl,CID_CELL,pkdu);
		--iCell;
		}
	    else {
		Smin[sm++] = min1;
		if (id != idSelf) mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = pkdu;
		}
#ifdef GASOLINE
	    if (smx->bGasOnly && pkdn->nGas == 0) goto LoadNotContained;
#endif
	    }
	pWant = pkdn->pLower + smx->nSmooth - smx->nQueue - 1;
#ifdef GASOLINE
	if (smx->bGasOnly) pEnd = pkdn->pLower + pkdn->nGas - 1;
	else pEnd = pkdn->pUpper;
#else
	pEnd = pkdn->pUpper;
#endif
	if (pWant > pEnd) {
	    for (pj=pkdn->pLower;pj<=pEnd;++pj) {
		if (id == idSelf) {
		    p = &smx->pkd->pStore[pj];
		    pq[smx->nQueue].bRemote = 0;
		    }
		else {
		    p = mdlAquire(mdl,CID_PARTICLE,pj,id);
		    pq[smx->nQueue].bRemote = 1;
		    }
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		pq[smx->nQueue].pPart = p;
		pq[smx->nQueue].fDist2 = dx*dx + dy*dy + dz*dz;
		pq[smx->nQueue].dx = dx;
		pq[smx->nQueue].dy = dy;
		pq[smx->nQueue].dz = dz;
		++smx->nQueue;
		}
	    }
	else {
	    for (pj=pkdn->pLower;pj<=pWant;++pj) {
		if (id == idSelf) {
		    p = &smx->pkd->pStore[pj];
		    pq[smx->nQueue].bRemote = 0;
		    }
		else {
		    p = mdlAquire(mdl,CID_PARTICLE,pj,id);
		    pq[smx->nQueue].bRemote = 1;
		    }
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		pq[smx->nQueue].pPart = p;
		pq[smx->nQueue].fDist2 = dx*dx + dy*dy + dz*dz;
		pq[smx->nQueue].dx = dx;
		pq[smx->nQueue].dy = dy;
		pq[smx->nQueue].dz = dz;
		++smx->nQueue;
		}
	    PQ_BUILD(pq,smx->nSmooth,pq);
	    for (;pj<=pEnd;++pj) {
		if (id == idSelf) p = &smx->pkd->pStore[pj];
		else p = mdlAquire(mdl,CID_PARTICLE,pj,id);
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < pq->fDist2) {
		    if (pq->bRemote) mdlRelease(mdl,CID_PARTICLE,pq->pPart);
		    if (id == idSelf) pq->bRemote = 0;
		    else pq->bRemote = 1;
		    pq->pPart = p;
		    pq->fDist2 = fDist2;
		    pq->dx = dx;
		    pq->dy = dy;
		    pq->dz = dz;
		    PQ_REPLACE(pq);
		    }
		else if (id != idSelf) mdlRelease(mdl,CID_PARTICLE,p);
		}
	    goto NoIntersect;  /* done loading phase */
	    }
#ifdef GASOLINE
    LoadNoIntersect:
#endif
	while (iCell == S[sp]) {
	    if (!sp) {
		if (id != idSelf) mdlRelease(mdl,CID_CELL,pkdn);
		return NULL;		/* EXIT, could not load enough particles! */
		}
	    --sp;
	    iCell = pkdn->iParent;
	    if (id == idSelf) pkdn = &c[iCell];
	    else {
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		}
	    }
#ifdef GASOLINE
    LoadNotContained:
#endif
	iCell ^= 1;
	if (id == idSelf) pkdn = &c[iCell];
	else {
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    }
	if (sm) --sm;
#ifdef GASOLINE
	if (smx->bGasOnly && pkdn->nGas == 0) {
	    iCell = pkdn->iParent;
	    if (id == idSelf) pkdn = &c[iCell];
	    else {
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		}
	    goto LoadNoIntersect;   			
	    }
#endif
	S[++sp] = iCell;
	}
    StartSearch:
    /*
    ** Start of PRIOQ searching loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (pkdn->iLower) {
	    iCell = pkdn->iLower;
	    if (id == idSelf) pkdn = &c[iCell];
	    else {
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		}
#ifdef GASOLINE
	    if (smx->bGasOnly) {
		if (pkdn->nGas == 0) min1 = pq->fDist2;
		else {
		    MINDIST(pkdn->bnd,r,min1);
		    }
		++iCell;
		if (id == idSelf) pkdu = &c[iCell];
		else pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
		if (pkdu->nGas == 0) min2 = pq->fDist2;
		else {
		    MINDIST(pkdu->bnd,r,min2);
		    }
		}
	    else {
		MINDIST(pkdn->bnd,r,min1);
		++iCell;
		if (id == idSelf) pkdu = &c[iCell];
		else pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
		MINDIST(pkdu->bnd,r,min2);
		}
#else
	    MINDIST(pkdn->bnd,r,min1);
	    ++iCell;
	    if (id == idSelf) pkdu = &c[iCell];
	    else pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
	    MINDIST(pkdu->bnd,r,min2);
#endif
	    if (min1 < min2) {
		Smin[sm++] = min2;
		--iCell;
		if (id != idSelf) mdlRelease(mdl,CID_CELL,pkdu);
		if (min1 >= pq->fDist2) goto NotContained;
		}
	    else {
		Smin[sm++] = min1;
		if (id != idSelf) mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = pkdu;
		if (min2 >= pq->fDist2) goto NotContained;
		}
	    }
#ifdef GASOLINE
	if (smx->bGasOnly) pEnd = pkdn->pLower + pkdn->nGas - 1;
	else pEnd = pkdn->pUpper;
#else
	pEnd = pkdn->pUpper;
#endif
	for (pj=pkdn->pLower;pj<=pEnd;++pj) {
	    if (id == idSelf) p = &smx->pkd->pStore[pj];
	    else p = mdlAquire(mdl,CID_PARTICLE,pj,id);
	    dx = r[0] - p->r[0];
	    dy = r[1] - p->r[1];
	    dz = r[2] - p->r[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < pq->fDist2) {
		if (pq->bRemote) mdlRelease(mdl,CID_PARTICLE,pq->pPart);
		if (id == idSelf) pq->bRemote = 0;
		else pq->bRemote = 1;
		pq->pPart = p;
		pq->fDist2 = fDist2;
		pq->dx = dx;
		pq->dy = dy;
		pq->dz = dz;
		PQ_REPLACE(pq);
		}
	    else if (id != idSelf) mdlRelease(mdl,CID_PARTICLE,p);
	    }
    NoIntersect:
	while (iCell == S[sp]) {
	    if (!sp) {
		if (id != idSelf) mdlRelease(mdl,CID_CELL,pkdn);
		return pq;
		}
	    --sp;
	    iCell = pkdn->iParent;
	    if (id == idSelf) pkdn = &c[iCell];
	    else {
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		}
	    }
    NotContained:
	iCell ^= 1;		
	if (id == idSelf) pkdn = &c[iCell];
	else {
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    }
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
#ifdef GASOLINE
	    if (smx->bGasOnly && pkdn->nGas == 0) {
		iCell = pkdn->iParent;
		if (id == idSelf) pkdn = &c[iCell];
		else {
		    mdlRelease(mdl,CID_CELL,pkdn);
		    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		    }
		goto NoIntersect;
		}
#endif
	    MINDIST(pkdn->bnd,r,min2);
	    }
	if (min2 >= pq->fDist2) {
	    iCell = pkdn->iParent;
	    if (id == idSelf) pkdn = &c[iCell];
	    else {
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
		}
	    goto NoIntersect;
	    }
	S[++sp] = iCell;
	}
    }


PQ *pqSearch(SMX smx,PQ *pq,PARTICLE *pi,FLOAT r[3],int bReplica,int *pbDone) {
    KDN *c = smx->pkd->kdTop;
    int idSelf = smx->pkd->idSelf;
    FLOAT *Smin = smx->SminT;
    int *S = smx->ST;
    FLOAT dMin,min1,min2;
    int j,iCell,id,iParent;
    int sp = 0;
    int sm = 0;

    *pbDone = 0;
    if (bReplica) iCell = ROOT;
    else {
	iCell = smx->pkd->iTopRoot;
	assert(c[iCell].pLower == idSelf);
	}
    if (iCell != ROOT) S[sp] = c[iCell].iParent;
    else S[sp] = iCell;	
    if (smx->nQueue == smx->nSmooth) goto StartSearch;
    /*
    ** Start of PRIOQ Loading loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
	    MINDIST(c[iCell].bnd,r,min1);
	    ++iCell;
	    MINDIST(c[iCell].bnd,r,min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		--iCell;
		}
	    else {
		Smin[sm++] = min1;
		}
#ifdef GASOLINE
	    if (smx->bGasOnly && c[iCell].nGas == 0) goto LoadNotContained;
#endif
	    }
	id = c[iCell].pLower;	/* this is the thread id in LTT */
	if (bReplica || id != idSelf) {
	    pq = pqSearchRemote(smx,pq,id,r);
	    }
	else {
	    pq = pqSearchLocal(smx,pi,r,pbDone);
	    if (*pbDone) return pq;	/* early exit */
	    }
	if (smx->nQueue == smx->nSmooth) goto NoIntersect;  /* done loading phase */
#ifdef GASOLINE
    LoadNoIntersect:
#endif
	while (iCell == S[sp]) {
	    if (!sp) {
		return NULL;		/* EXIT, could not load enough particles! */
		}
	    --sp;
	    iCell = c[iCell].iParent;
	    }
#ifdef GASOLINE
    LoadNotContained:
#endif
	iCell ^= 1;
	if (sm) --sm;
#ifdef GASOLINE
	if (smx->bGasOnly && c[iCell].nGas == 0) {
	    iCell = c[iCell].iParent;
	    goto LoadNoIntersect;   			
	    }
#endif
	S[++sp] = iCell;
	}
    /*
    ** Start of PRIOQ searching loop.
    */
    StartSearch:
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
#ifdef GASOLINE
	    if (smx->bGasOnly) {
		if (c[iCell].nGas == 0) min1 = pq->fDist2;
		else {
		    MINDIST(c[iCell].bnd,r,min1);
		    }
		++iCell;
		if (c[iCell].nGas == 0) min2 = pq->fDist2;
		else {
		    MINDIST(c[iCell].bnd,r,min2);
		    }
		}
	    else {
		MINDIST(c[iCell].bnd,r,min1);
		++iCell;
		MINDIST(c[iCell].bnd,r,min2);
		}
#else
	    MINDIST(c[iCell].bnd,r,min1);
	    ++iCell;
	    MINDIST(c[iCell].bnd,r,min2);
#endif
	    if (min1 < min2) {
		Smin[sm++] = min2;
		--iCell;
		if (min1 >= pq->fDist2) goto NotContained;
		}
	    else {
		Smin[sm++] = min1;
		if (min2 >= pq->fDist2) goto NotContained;
		}
	    }
	id = c[iCell].pLower;	/* this is the thread id in LTT */
	pq = pqSearchRemote(smx,pq,id,r);
    NoIntersect:
	while (iCell == S[sp]) {
	    if (sp) {
		--sp;
		iCell = c[iCell].iParent;
		}
	    else if (!bReplica) {
		/*
		** Containment Test! 
		*/
		for (j=0;j<3;++j) {
		    dMin = c[iCell].bnd.fMax[j] - 
			fabs(c[iCell].bnd.fCenter[j] - r[j]);
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = c[iCell].iParent;
			if (!iParent) {
			    *pbDone = 0;
			    return pq;
			    }
			S[sp] = iParent;
			goto NotContained;
			}
		    }
		*pbDone = 1;
		return pq;
		}
	    else return pq;
	    }
    NotContained:
	iCell ^= 1;		
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
#ifdef GASOLINE
	    if (smx->bGasOnly && c[iCell].nGas == 0) {
		iCell = c[iCell].iParent;
		goto NoIntersect;
		}
#endif
	    MINDIST(c[iCell].bnd,r,min2);
	    }
	if (min2 >= pq->fDist2) {
	    iCell = c[iCell].iParent;
	    goto NoIntersect;
	    }
	S[++sp] = iCell;
	}
    }


void smSmooth(SMX smx,SMF *smf)
    {
    PKD pkd = smx->pkd;
    PARTICLE *p = pkd->pStore;
    PQ *pq;
    FLOAT r[3],fBall;
    int iStart[3],iEnd[3];
    int pi,i,j,bDone;
    int ix,iy,iz;
   
    for (pi=0;pi<pkd->nLocal;++pi) {
	if (!TYPETest(&(p[pi]),smx->eParticleTypes)) continue;
	pq = NULL;
	smx->nQueue = 0;
	pq = pqSearch(smx,pq,&p[pi],p[pi].r,0,&bDone);
	/*
	** Search in replica boxes if it is required.
	*/
	if (!bDone && smx->bPeriodic) {
	    /*
	    ** Note for implementing SLIDING PATCH, the offsets for particles are
	    ** negative here, reflecting the relative +ve offset of the simulation
	    ** volume.
	    */
	    fBall = sqrt(pq->fDist2);
	    for (j=0;j<3;++j) {
		iStart[j] = floor((p[pi].r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((p[pi].r[j] + fBall)/pkd->fPeriod[j] + 0.5);
		}
	    for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = p[pi].r[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		    r[1] = p[pi].r[1] - iy*pkd->fPeriod[1];
		    for (iz=iStart[2];iz<=iEnd[2];++iz) {
			r[2] = p[pi].r[2] - iz*pkd->fPeriod[2];
			if (ix || iy || iz) {
			    pq = pqSearch(smx,pq,&p[pi],r,1,&bDone);
			    }
			}
		    }	
		}
	    }
	/*
	** Should we ever get tripped by this assert, it means that 
	** for some reason there are less particles in the box (and 
	** replicas) of the desired type than are requested by 
	** nSmooth! We die on this condition at the moment, but maybe
	** there are sensible cases to be dealt with here.
	*/
	assert(pq != NULL);
	/*
	** Maybe should give a warning if the search radius is not conained
	** within the replica volume.
	*/
	
	p[pi].fBall = sqrt(pq->fDist2);
	for (i=0;i<smx->nSmooth;++i) {
	    smx->nnList[i].pPart = smx->pq[i].pPart;
	    smx->nnList[i].fDist2 = smx->pq[i].fDist2;			
	    smx->nnList[i].dx = smx->pq[i].dx;
	    smx->nnList[i].dy = smx->pq[i].dy;
	    smx->nnList[i].dz = smx->pq[i].dz;
	    }
	
	/*
	** Apply smooth funtion to the neighbor list.
	*/
	smx->fcnSmooth(&p[pi],smx->nSmooth,smx->nnList,smf);
	/*
	** Call mdlCacheCheck to make sure we are making progress!
	*/
	mdlCacheCheck(pkd->mdl);
	/*
	** Release aquired pointers.
	*/
	for (i=0;i<smx->nSmooth;++i) {
	    if (smx->pq[i].bRemote) {
		mdlRelease(pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
		}
	    }
	}
    }


void smGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[3])
    {
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,nCnt,pEnd;
    int idSelf;

    idSelf = smx->pkd->idSelf;
    nCnt = smx->nnListSize;
    iCell = ROOT;
    while (1) {
#ifdef GASOLINE
	if (smx->bGasOnly && c[iCell].nGas == 0) {
	    goto NoIntersect;
	    }
#endif
	MINDIST(c[iCell].bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
	    S[sp++] = iCell+1;
	    continue;
	    }
	else {
#ifdef GASOLINE
	    if (smx->bGasOnly) pEnd = c[iCell].pLower + c[iCell].nGas - 1;
	    else pEnd = c[iCell].pUpper;
#else
	    pEnd = c[iCell].pUpper;
#endif
	    for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
		dx = r[0] - p[pj].r[0];
		dy = r[1] - p[pj].r[1];
		dz = r[2] - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    if(nCnt >= smx->nnListMax) {
			smx->nnListMax += NNLIST_INCREMENT;
			smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
			assert(smx->nnList != NULL);
			smx->nnbRemote = realloc(smx->nnbRemote,smx->nnListMax*sizeof(int));
			assert(smx->nnbRemote != NULL);
			}
		    smx->nnList[nCnt].fDist2 = fDist2;
		    smx->nnList[nCnt].dx = dx;
		    smx->nnList[nCnt].dy = dy;
		    smx->nnList[nCnt].dz = dz;
		    smx->nnList[nCnt].pPart = &p[pj];
		    smx->nnList[nCnt].iPid = idSelf;
		    smx->nnList[nCnt].iIndex = pj;
		    smx->nnbRemote[nCnt] = 0;
		    ++nCnt;
		    }
		}
	    }
    NoIntersect:
	if (sp) iCell = S[--sp];
	else break;
	}
    smx->nnListSize = nCnt;
    }


void smGatherRemote(SMX smx,FLOAT fBall2,FLOAT r[3],int id)
    {
    MDL mdl = smx->pkd->mdl;
    KDN *pkdn;
    PARTICLE *pp;
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int pj,nCnt,pEnd;
    int iCell;

    assert(id != smx->pkd->idSelf);
    nCnt = smx->nnListSize;
    iCell = ROOT;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    while (1) {
#ifdef GASOLINE
	if (smx->bGasOnly && pkdn->nGas == 0) {
	    goto NoIntersect;
	    }
#endif
	MINDIST(pkdn->bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (pkdn->iLower) {
	    iCell = pkdn->iLower;
	    S[sp++] = iCell+1;
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    continue;
	    }
	else {
#ifdef GASOLINE
	    if (smx->bGasOnly) pEnd = pkdn->pLower + pkdn->nGas - 1;
	    else pEnd = pkdn->pUpper;
#else
	    pEnd = pkdn->pUpper;
#endif
	    for (pj=pkdn->pLower;pj<=pEnd;++pj) {
		pp = mdlAquire(mdl,CID_PARTICLE,pj,id);
		dx = r[0] - pp->r[0];
		dy = r[1] - pp->r[1];
		dz = r[2] - pp->r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    if(nCnt >= smx->nnListMax) {
			smx->nnListMax += NNLIST_INCREMENT;
			smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
			assert(smx->nnList != NULL);
			smx->nnbRemote = realloc(smx->nnbRemote,smx->nnListMax*sizeof(int));
			assert(smx->nnbRemote != NULL);
			}
		    smx->nnList[nCnt].fDist2 = fDist2;
		    smx->nnList[nCnt].dx = dx;
		    smx->nnList[nCnt].dy = dy;
		    smx->nnList[nCnt].dz = dz;
		    smx->nnList[nCnt].pPart = pp;
		    smx->nnList[nCnt].iPid = id;
		    smx->nnList[nCnt].iIndex = pj;
		    smx->nnbRemote[nCnt] = 1;
		    ++nCnt;
		    }
		}
	    }
    NoIntersect:
	if (sp) {
	    iCell = S[--sp];
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    }
	else break;
	}
    mdlRelease(mdl,CID_CELL,pkdn);
    smx->nnListSize = nCnt;
    }


void smGather(SMX smx,FLOAT fBall2,FLOAT r[3])
    {
    KDN *c = smx->pkd->kdTop;
    int idSelf = smx->pkd->idSelf;
    int *S = smx->ST;
    FLOAT min2;
    int iCell,id;
    int sp = 0;

    iCell = ROOT;
    while (1) {
#ifdef GASOLINE
	if (smx->bGasOnly && c[iCell].nGas == 0) {
	    goto NoIntersect;
	    }
#endif
	MINDIST(c[iCell].bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
	    S[sp++] = iCell+1;
	    continue;
	    }
	else {
	    id = c[iCell].pLower; /* this is the thread id in LTT */
	    if (id != idSelf) {
		smGatherRemote(smx,fBall2,r,id);
		}
	    else {
		smGatherLocal(smx,fBall2,r);
		}
	    }
    NoIntersect:
	if (sp) iCell = S[--sp];
	else break;
	}
    }


void smReSmooth(SMX smx,SMF *smf)
    {
    PKD pkd = smx->pkd;
    PARTICLE *p = pkd->pStore;
    FLOAT r[3],fBall;
    int iStart[3],iEnd[3];
    int pi,i,j;
    int ix,iy,iz;

    for (pi=0;pi<pkd->nLocal;++pi) {
	if (!TYPETest(&(p[pi]),smx->eParticleTypes)) continue;
	smx->nnListSize = 0;
	/*
	** Note for implementing SLIDING PATCH, the offsets for particles are
	** negative here, reflecting the relative +ve offset of the simulation
	** volume.
	*/
	fBall = p[pi].fBall;
	if (smx->bPeriodic) {
	    for (j=0;j<3;++j) {
		iStart[j] = floor((p[pi].r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((p[pi].r[j] + fBall)/pkd->fPeriod[j] + 0.5);
		}
	    for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = p[pi].r[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		    r[1] = p[pi].r[1] - iy*pkd->fPeriod[1];
		    for (iz=iStart[2];iz<=iEnd[2];++iz) {
			r[2] = p[pi].r[2] - iz*pkd->fPeriod[2];
			smGather(smx,fBall*fBall,r);
			}
		    }
		}
	    }
	else {
	    smGather(smx,fBall*fBall,p[pi].r);
	    }
	/*
	** Apply smooth funtion to the neighbor list.
	*/
	smx->fcnSmooth(&p[pi],smx->nnListSize,smx->nnList,smf);
	/*
	** Release aquired pointers.
	*/
	for (i=0;i<smx->nnListSize;++i) {
	    if (smx->nnbRemote[i]) {
		mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
		}
	    }
	}
    }


typedef struct CheckElt {
    int iCell;
    int id;
    FLOAT rOffset[3];
    } CELT;

typedef struct CheckStack {
    CELT *Check;
    int nCheck;
    } CSTACK;


#define WALK_MINMULTIPOLE	3

#ifdef GASOLINE
/*
** This algorithm is much like the one used in walk.c for gravity.
** It depends on having the bound-of-balls in the tree structure and
** becuase this is not present in the pure gravity case this is only 
** a meaningful function for GASOLINE and maybe COLLISIONS.
*/   
void smReSmoothWalk(SMX smx,SMF *smf)
    {
    PKD pkd = smx->pkd;
    PARTICLE *p = smx->pkd->pStore;
    PARTICLE *pRemote;
    KDN *c;
    KDN *pkdc;
    CSTACK *S;
    int iStack,ism;
    CELT *Check;
    FLOAT dMin,min2,fBall2;
    FLOAT rOffset[3];
    int ix,iy,iz,bRep;
    int nMaxCheck,nCheck;
    int iCell,iSib,iCheckCell;
    int i,ii,j,n,id,pi,pj,nTotActive;
    int nPart,nMaxPart;
    ILP *ilp;
    FLOAT x,y,z,d2;
    int nAccept,nFail;

    int nReps = 1;

    /*
    ** Initially we set our cell pointer to 
    ** point to the top tree.
    */
    c = pkd->kdTop;
    /*
    ** Allocate initial interaction lists.
    */
    nFail = 0;
    nTotActive = 0;
    nPart = 0;
    nMaxPart = 500;
    ilp = malloc(nMaxPart*sizeof(ILP));
    assert(ilp != NULL);
    /*
    ** Allocate Checklist.
    */
    nMaxCheck = 2*nReps+1;
    nMaxCheck = nMaxCheck*nMaxCheck*nMaxCheck;	/* all replicas */
    iCell = pkd->iTopRoot;
    while (iCell = c[iCell].iParent) ++nMaxCheck; /* all top tree siblings */
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
#ifdef SLIDING_PATCH
	    rOffset[1] += SHEAR(ix,pkd->fPeriod[0],pkd->fPeriod[1],
				pkd->dTime,pkd->dOrbFreq);
#endif
	    for (iz=-nReps;iz<=nReps;++iz) {
		rOffset[2] = iz*pkd->fPeriod[2];
		bRep = ix || iy || iz;
		if (!bRep) continue;
		Check[nCheck].iCell = ROOT;
		Check[nCheck].id = -1;
		for (j=0;j<3;++j) Check[nCheck].rOffset[j] = rOffset[j];
		++nCheck;
		}
	    }
	}
    iCell = pkd->iTopRoot;
    iSib = SIBLING(iCell);
    while (iSib) {
	Check[nCheck].iCell = iSib;
	Check[nCheck].id = -1;
	for (j=0;j<3;++j) Check[nCheck].rOffset[j] = 0.0;
	++nCheck;
	iCell = c[iCell].iParent;
	iSib = SIBLING(iCell);
	}
    /*
    ** Now switch the cell pointer to point to 
    ** the local tree.
    */
    c = pkd->kdNodes;
    iCell = ROOT;
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
		    }
		else if (id < 0) {
		    pkdc = &pkd->kdTop[Check[i].iCell];
		    }
		else {
		    pkdc = mdlAquire(pkd->mdl,CID_CELL,Check[i].iCell,id);
		    }
		if (c[iCell].iLower) {
		    for (j=0;j<3;++j) {
			dMin = fabs(pkdc->bnd.fCenter[j] + Check[i].rOffset[j] - c[iCell].bndBall.fCenter[j]);
			dMin -= c[iCell].bndBall.fMax[j] + pkdc->bnd.fMax[j];
			if (dMin > 0) goto SkipCheck;
			}
		    for (j=0;j<3;++j) {
			dMin = fabs(pkdc->bnd.fCenter[j] + Check[i].rOffset[j] - c[iCell].bndBall.fCenter[j]);
			dMin -= c[iCell].bndBall.fMax[j] - pkdc->bnd.fMax[j];
			if (dMin > 0) goto SkipOpen;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			dMin = fabs(pkdc->bnd.fCenter[j] + Check[i].rOffset[j] - c[iCell].bndBall.fCenter[j]);
			dMin -= c[iCell].bndBall.fMax[j] + pkdc->bnd.fMax[j];
			if (dMin > 0) goto SkipCheck;
			}
		    }
		if (pkdc->iLower) {
		    iCheckCell = pkdc->iLower;
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
		    Check[nCheck].iCell = iCheckCell;
		    Check[nCheck+1] = Check[i];
		    Check[nCheck+1].iCell = iCheckCell+1;
		    /*
		    ** If we are opening a leaf of the top tree
		    ** we need to correctly set the processor id.
		    ** (this is a bit tricky)
		    */
		    if (id < 0 && pkdc->pLower >= 0) {
			Check[nCheck].id = pkdc->pLower;
			Check[nCheck+1].id = pkdc->pLower;
			}
		    nCheck += 2;
		    goto SkipCheck;		/* don't want to save this one */
		    }
	    SkipOpen:
		Check[ii++] = Check[i];	/* save this one */
	    SkipCheck:
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
	    Check[nCheck].iCell = iCell + 1;
	    Check[nCheck].id = pkd->idSelf;
	    for (j=0;j<3;++j) Check[nCheck].rOffset[j] = 0.0;
	    ++nCheck;
	    /*
	    ** Push Checklist for the sibling onto the stack
	    ** before proceeding deeper in the tree.
	    */
	    ++iStack;
	    S[iStack].nCheck = nCheck;
	    for (i=0;i<nCheck;++i) S[iStack].Check[i] = Check[i];
	    S[iStack].Check[nCheck-1].iCell = iCell;
	    }
	/*
	** Now the interaction list should be complete.
	** Now test each check cell for containment in the ball of each 
	** particle.
	*/
#if 0
	nPart = 0;
	for (pj=c[iCell].pLower;pj<=c[iCell].pUpper;++pj) {
	    fBall2 = p[pj].fBall*p[pj].fBall;
	    nAccept = 0;
	    for (i=0;i<nCheck;++i) {
		id = Check[i].id;
		if (id == pkd->idSelf) {
		    pkdc = &pkd->kdNodes[Check[i].iCell];
		    }
		else if (id < 0) {
		    pkdc = &pkd->kdTop[Check[i].iCell];
		    }
		else {
		    pkdc = mdlAquire(pkd->mdl,CID_CELL,Check[i].iCell,id);
		    }
		for (j=0;j<3;++j) {
		    dMin = fabs(pkdc->bnd.fCenter[j] + Check[i].rOffset[j] - p[pj].r[j]);
		    dMin -= pkdc->bnd.fMax[j] + p[pj].fBall;
		    if (dMin > 0) goto SkipCell;
		    }
		min2 = 0;
		for (j=0;j<3;++j) {
		    dMin = fabs(pkdc->bnd.fCenter[j] + Check[i].rOffset[j] - p[pj].r[j]);
		    dMin -= pkdc->bnd.fMax[j];
		    if (dMin > 0) min2 += dMin*dMin;
		    }
		if (min2 > fBall2) goto SkipCell;

		nPart += pkdc->pUpper - pkdc->pLower + 1;

		/*
		** Now open this bucket and calculate real distances.
		*/
		if (id == pkd->idSelf) {
		    for (pi=pkdc->pLower;pi<=pkdc->pUpper;++pi) {						
			x = p[pi].r[0] - p[pj].r[0] + Check[i].rOffset[0];
			y = p[pi].r[1] - p[pj].r[1] + Check[i].rOffset[1];
			z = p[pi].r[2] - p[pj].r[2] + Check[i].rOffset[2];
			d2 = x*x + y*y + z*z;
			if (d2 < fBall2) ++nAccept;
			}
		    }
		else {
		    for (pi=pkdc->pLower;pi<=pkdc->pUpper;++pi) {
			pRemote = mdlAquire(pkd->mdl,CID_PARTICLE,pi,id);
			x = pRemote->r[0] - p[pj].r[0] + Check[i].rOffset[0];
			y = pRemote->r[1] - p[pj].r[1] + Check[i].rOffset[1];
			z = pRemote->r[2] - p[pj].r[2] + Check[i].rOffset[2];
			d2 = x*x + y*y + z*z;
			if (d2 < fBall2) ++nAccept;
			mdlRelease(pkd->mdl,CID_PARTICLE,pRemote);
			}
		    }
	    SkipCell:
		if (id >= 0 && id != pkd->idSelf) {
		    mdlRelease(pkd->mdl,CID_CELL,pkdc);
		    }
		}
	    if (nAccept < 32) ++nFail;
	    }

	*pdPartSum += nPart;

	n = c[iCell].pUpper - c[iCell].pLower + 1;
	nTotActive += n;

#endif
	free(Check);
	while (iCell & 1) {
	    iCell = c[iCell].iParent;
	    if (!iCell) {
		/*
		** Make sure stack is empty and free its storage.
		*/
		assert(iStack == -1);
		free(S);
		/*
		** Free checklist storage.
		*/
		free(Check);
		}
	    }
	++iCell;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	nCheck = S[iStack].nCheck;
	for (i=0;i<nCheck;++i) Check[i] = S[iStack].Check[i];
	--iStack;
	}

    }


void smBallScatterLocal(SMX smx,FLOAT r[3],int iMarkType)
    {
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    FLOAT dMin,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,j,pj,pEnd;

    assert(smx->bGasOnly);
    iCell = ROOT;
    while (1) {
	if (c[iCell].nGas == 0) {
	    goto NoIntersect;
	    }
	for (j=0;j<3;++j) {
	    dMin = fabs(c[iCell].bndBall.fCenter[j] - r[j]);
	    if (dMin > c[iCell].bndBall.fMax[j]) goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
	    S[sp++] = iCell+1;
	    continue;
	    }
	else {
	    pEnd = c[iCell].pLower + c[iCell].nGas - 1;
	    for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
		dx = r[0] - p[pj].r[0];
		dy = r[1] - p[pj].r[1];
		dz = r[2] - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < p[pj].fBallMax*p[pj].fBallMax) {
		    TYPESet(&p[pj],iMarkType);
		    }
		}
	    }
    NoIntersect:
	if (sp) iCell = S[--sp];
	else break;
	}
    }


void smBallScatterRemote(SMX smx,FLOAT r[3],int iMarkType,int id)
    {
    MDL mdl = smx->pkd->mdl;
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p = smx->pkd->pStore;
    KDN *pkdn;
    PARTICLE *pp;
    int *S = smx->S;
    FLOAT dMin,dx,dy,dz,fDist2;
    int sp = 0;
    int iCell,j,pj,pEnd;

    assert(smx->bGasOnly);
    assert(id != smx->pkd->idSelf);
    iCell = ROOT;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    while (1) {
	if (pkdn->nGas == 0) {
	    goto NoIntersect;
	    }
	for (j=0;j<3;++j) {
	    dMin = fabs(pkdn->bndBall.fCenter[j] - r[j]);
	    if (dMin > pkdn->bndBall.fMax[j]) {
		goto NoIntersect;
		}
	    }
	/*
	** We have an intersection to test.
	*/
	if (pkdn->iLower) {
	    iCell = pkdn->iLower;
	    S[sp++] = iCell+1;
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    continue;
	    }
	else {
	    pEnd = pkdn->pLower + pkdn->nGas - 1;
	    for (pj=pkdn->pLower;pj<=pEnd;++pj) {
		pp = mdlAquire(mdl,CID_PARTICLE,pj,id);
		dx = r[0] - pp->r[0];
		dy = r[1] - pp->r[1];
		dz = r[2] - pp->r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < pp->fBallMax*pp->fBallMax) {
		    TYPESet(pp,iMarkType);						
		    }
		mdlRelease(mdl,CID_PARTICLE,pp);
		}
	    }
    NoIntersect:
	if (sp) {
	    iCell = S[--sp];
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    }
	else break;
	}
    mdlRelease(mdl,CID_CELL,pkdn);
    }


void smBallScatter(SMX smx,FLOAT r[3],int iMarkType)
    {
    KDN *c = smx->pkd->kdTop;
    int *S = smx->ST;
    FLOAT dMin;
    int idSelf = smx->pkd->idSelf;
    int j,iCell,id,iParent;
    int sp = 0;

    assert(smx->bGasOnly);
    iCell = ROOT;
    while (1) {
	if (c[iCell].nGas == 0) {
	    goto NoIntersect;
	    }
	for (j=0;j<3;++j) {
	    dMin = fabs(c[iCell].bndBall.fCenter[j] - r[j]);
	    if (dMin > c[iCell].bndBall.fMax[j]) {
		goto NoIntersect;
		}
	    }
	/*
	** We have an intersection to test.
	*/
	if (c[iCell].iLower) {
	    iCell = c[iCell].iLower;
	    S[sp++] = iCell+1;
	    continue;
	    }
	else {
	    id = c[iCell].pLower; /* this is the thread id in LTT */
	    if (id != idSelf) {
		smBallScatterRemote(smx,r,iMarkType,id);
		}
	    else {
		smBallScatterLocal(smx,r,iMarkType);
		}
	    }
    NoIntersect:
	if (sp) iCell = S[--sp];
	else break;
	}
    }


void smMarkSmooth(SMX smx,int iMarkType)
    {
    PKD pkd = smx->pkd;
    PARTICLE *p = pkd->pStore;
    PQ *pq;
    BND *pbnd;
    FLOAT r[3];
    FLOAT min[3],max[3];
    int iStart[3],iEnd[3];
    int pi,i,j;
    int ix,iy,iz;

    if (smx->bPeriodic) {
	pbnd = &smx->pkd->kdTop[ROOT].bndBall;
	for (j=0;j<3;++j) {
	    min[j] = pbnd->fCenter[j] + pbnd->fMax[j] - 0.5*pkd->fPeriod[j];
	    if (min[j] < 0) min[j] = 0;
	    max[j] = pbnd->fCenter[j] - pbnd->fMax[j] + 0.5*pkd->fPeriod[j];
	    if (max[j] < 0) max[j] = 0;
	    }
	}
    for (pi=0;pi<pkd->nLocal;++pi) {
	if (!TYPETest(&(p[pi]),smx->eParticleTypes)) continue;
	/*
	** Note for implementing SLIDING PATCH, the offsets for particles are
	** negative here, reflecting the relative +ve offset of the simulation
	** volume.
	*/
	if (smx->bPeriodic) {
	    for (j=0;j<3;++j) {
		iStart[j] = floor((p[pi].r[j] - min[j])/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((p[pi].r[j] + max[j])/pkd->fPeriod[j] + 0.5);
		}
	    for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = p[pi].r[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		    r[1] = p[pi].r[1] - iy*pkd->fPeriod[1];
		    for (iz=iStart[2];iz<=iEnd[2];++iz) {
			r[2] = p[pi].r[2] - iz*pkd->fPeriod[2];
			smBallScatter(smx,r,iMarkType);
			}
		    }
		}
	    }
	else {
	    smBallScatter(smx,p[pi].r,iMarkType);
	    }
	}
    }
#endif /* GASOLINE */

void smFof(SMX smx,int nFOFsDone,SMF *smf)
    {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p = smx->pkd->pStore;
    FOFRM* rm;
    FOFRM* rmTmp;
    FOFPG* protoGroup;		
    FOFBIN *bin;

    int pi,pn,pnn,nCnt,i,j,k,nRmCnt, iRmIndex,nRmListSize,iListIndex,nListSize;
    int nFifo, iHead, iTail, iMaxGroups, iGroup;
    int *Fifo, *lista, *listb, *eClass;
    int iStart[3],iEnd[3];
    int ix,iy,iz;

    FLOAT r[3],l[3],relpos[3],lx,ly,lz,fBall,fBall2Max,rho;
    int nTree,cnt,tmp;
    cnt = 0;

    if (smx->bPeriodic) {
	lx = pkd->fPeriod[0];
	ly = pkd->fPeriod[1];
	lz = pkd->fPeriod[2];
	}
    else {
	lx = FLOAT_MAXVAL;
	ly = FLOAT_MAXVAL;
	lz = FLOAT_MAXVAL;
	}
    l[0] = lx;
    l[1] = ly;
    l[2] = lz;
	
    iHead = 0;
    iTail = 0;
    iMaxGroups = pkd->kdNodes[ROOT].pUpper + 1;
    nRmListSize = iMaxGroups;
    nListSize = iMaxGroups;
    tmp = pkd->nDark+pkd->nGas+pkd->nStar;

    iGroup = 0; 
    nTree = pkd->kdNodes[ROOT].pUpper + 1;
    nFifo = pkdLocal(pkd);
    Fifo = (int *)malloc(nFifo*sizeof(int));
	
    protoGroup = (FOFPG *)malloc(iMaxGroups*sizeof(FOFPG));
    /*used something smaller than FOFGD here to reduce memory usage*/
    assert(protoGroup != NULL);

    rm = (FOFRM *)malloc(nRmListSize*sizeof(FOFRM));
    assert(rm != NULL);

    lista = (int *)malloc(nListSize*sizeof(int));
    assert(lista != NULL);
	
    listb = (int *)malloc(nListSize*sizeof(int));
    assert(listb != NULL);

    iRmIndex = 0;
    iListIndex = 0;
    pkd->nGroups = 0;
    pkd->nMaxRm = 0;
    fBall2Max = 0.0;
    if( nFOFsDone > 0){
	mdlROcache(mdl,CID_BIN,pkd->groupBin,sizeof(FOFBIN),pkd->nBins);
	if(pkd->idSelf != 0){
	    for(i=0; i< pkd->nBins; i++){
		bin = mdlAquire(mdl,CID_BIN,i,0);
		mdlRelease(mdl,CID_BIN,bin);
		pkd->groupBin[i] = *bin;	
		}
	    }	
	mdlFinishCache(mdl,CID_BIN);
	for (pn=0;pn<nTree;pn++) {
	    if( p[pn].pBin >= 0 ){
		for(j = 0; j < 3; j++)	{
		    relpos[j] = corrPos(pkd->groupBin[p[pn].pBin].com[j], p[pn].r[j], l[j]) 
			- pkd->groupBin[p[pn].pBin].com[j];
		    }
		rho = pkd->groupBin[p[pn].pBin].fDensity;
		if(rho > p[pn].fDensity) rho =  p[pn].fDensity;
		p[pn].fBall = pow(p[pn].fMass/(rho*smf->fContrast),2.0/3.0);
		/* set velocity linking lenght in case of a phase space FOF */
		p[pn].fBallv2 = pkd->groupBin[p[pn].pBin].v2[0]+ 
		    pkd->groupBin[p[pn].pBin].v2[1]+ pkd->groupBin[p[pn].pBin].v2[2];
		p[pn].fBallv2 *= 2.0;
		p[pn].fBallv2 *= pow(smf->fContrast,-2.0/3.0);
		if(p[pn].fBall > smf->dTau2*pow(p[pn].fMass/ smf->fContrast,2.0/3.0) )
		    p[pn].fBall = smf->dTau2*pow(p[pn].fMass/ smf->fContrast,2.0/3.0);
		} else {
		    p[pn].fBall = 0.0;
		    }
	    if(p[pn].fBall > fBall2Max)fBall2Max = p[pn].fBall;
	    p[pn].pGroup = 0;
	    }	
	} else {
	    for (pn=0;pn<nTree;pn++) {
		p[pn].pGroup = 0;
		p[pn].fBallv2 = -1.0;
		if(smf->bTauAbs) {
		    p[pn].fBall = smf->dTau2;
		    } else {
			p[pn].fBall = smf->dTau2*pow(p[pn].fMass,0.6666);
			}
		if(p[pn].fBall > fBall2Max)fBall2Max = p[pn].fBall;
		if(p[pn].fBall > fBall2Max)fBall2Max = p[pn].fBall;
		}	
	    }
    /* Have to restart particle chache, since we will need 
     * the updated p[pn].fBall now */
    mdlFinishCache(mdl,CID_PARTICLE);
    mdlROcache(mdl,CID_PARTICLE,p,sizeof(PARTICLE),nTree);	
    /* Starting FOF search now... */
    for (pn=0;pn<nTree;pn++) {
	if (p[pn].pGroup ) continue;
	iGroup++;
	protoGroup[iGroup].nMembers = 0;
	protoGroup[iGroup].iId = iGroup;		
	nRmCnt = 0;	
	/*
	** Mark particle and add it to the do-fifo
	*/
	p[pn].pGroup = iGroup;
	Fifo[iTail++] = pn;
	if(iTail == nFifo) iTail = 0;
	while(iHead != iTail) {
	    pi = Fifo[iHead++];
	    if(iHead == nFifo) iHead=0;
	    /*
	    ** Do a Ball Gather at the radius p[pi].fBall
	    */
	    smx->nnListSize =0;
	    fBall = sqrt(p[pi].fBall);
	    if (smx->bPeriodic) {
		for (j=0;j<3;++j) {
		    iStart[j] = floor((p[pi].r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		    iEnd[j] = floor((p[pi].r[j] + fBall)/pkd->fPeriod[j] + 0.5);
		    }
		for (ix=iStart[0];ix<=iEnd[0];++ix) {
		    r[0] = p[pi].r[0] - ix*pkd->fPeriod[0];
		    for (iy=iStart[1];iy<=iEnd[1];++iy) {
			r[1] = p[pi].r[1] - iy*pkd->fPeriod[1];
			for (iz=iStart[2];iz<=iEnd[2];++iz) {
			    r[2] = p[pi].r[2] - iz*pkd->fPeriod[2];
			    smGather(smx,p[pi].fBall,r);
			    }
			}
		    }
		}
	    else {
		smGather(smx,p[pi].fBall,p[pi].r);
		}
	    nCnt = smx->nnListSize;
	    for(pnn=0;pnn<nCnt;++pnn ){
		if(smx->nnbRemote[pnn] == 0){
		    /* Do not add particles that are allready in this group*/
		    if (smx->nnList[pnn].pPart->pGroup == iGroup ) continue;
		    /* Check phase space distance */
		    if (p->fBallv2 > 0.0){
			if(phase_dist(p[pi],*smx->nnList[pnn].pPart) > 1.0) continue;
			}
		    /* Add to equivalence relations if this particle is allready in another group */
		    if (smx->nnList[pnn].pPart->pGroup > 0){
			if(lista[iListIndex-1] == iGroup && listb[iListIndex-1] == smx->nnList[pnn].pPart->pGroup) continue;
			iListIndex++;
			if(iListIndex >= nListSize){
			    nListSize = 2 * nListSize;
			    lista = (int *) realloc(lista,nListSize*sizeof(int));
			    assert(lista != NULL);
			    listb = (int *) realloc(listb,nListSize*sizeof(int));
			    assert(listb != NULL);
			    }
			lista[iListIndex] = iGroup;
			listb[iListIndex] = smx->nnList[pnn].pPart->pGroup;
			continue;
			}
		    /* 
		    **  Mark particle and add it to the do-fifo
		    */
		    smx->nnList[pnn].pPart->pGroup = iGroup;
		    Fifo[iTail++] = smx->nnList[pnn].iIndex;
		    if(iTail == nFifo) iTail = 0;
		    } else {	 /* Nonlocal neighbors: */
			/* Make remote member linking symmetric by using smaller linking lenght if different: */
			if (p[pi].fBallv2 > 0.0) { /* Check phase space distance */
			    if(phase_dist(p[pi],*smx->nnList[pnn].pPart) > 1.0 || phase_dist(*smx->nnList[pnn].pPart,p[pi]) < 1.0) continue;
			    } else { /* real space distance */
				if(smx->nnList[pnn].fDist2 > smx->nnList[pnn].pPart->fBall) continue;
				}
 				
			/* Add to RM list if new */ 
			rm[iRmIndex].iIndex = smx->nnList[pnn].iIndex ;
			rm[iRmIndex].iPid = smx->nnList[pnn].iPid;						 				  
			if(bsearch (rm+iRmIndex,rm+iRmIndex-nRmCnt,nRmCnt, sizeof(FOFRM),CmpRMs) == NULL ){ 
			    nRmCnt++;	
			    iRmIndex++;
			    qsort(rm+iRmIndex-nRmCnt,nRmCnt,sizeof(FOFRM),CmpRMs);
			    }
			}					
		}
	    }	
	/* FIFO done for this group, add remote Members to the group data before doing next group: */
	protoGroup[iGroup].iFirstRm = iRmIndex - nRmCnt;
	protoGroup[iGroup].nRemoteMembers = nRmCnt;
	if( nRmCnt > pkd->nMaxRm ) pkd->nMaxRm = nRmCnt;
	}
	
    free(Fifo); 
    if(iListIndex > 1){
	eClass = (int *)malloc((iGroup+1)*sizeof(int));
	assert(eClass != NULL);
	for(k=1; k<=iGroup; k++) eClass[k] = k;
	for(i=1; i<=iListIndex; i++) {
	    j=lista[i];
	    while(eClass[j] != j) j=eClass[j];
	    k=listb[i];
	    while(eClass[k] != k) k=eClass[k];
	    if(j != k) eClass[j]=k;
	    }
	for(j=1;j<=iGroup;j++){
	    while(eClass[j] != eClass[eClass[j]]) eClass[j]=eClass[eClass[j]];	
	    }
	/* Make smallest id the eClass id: */
	for(j=1;j<=iGroup;j++){
	    if(j < eClass[j]){
		for(k=j+1; k<=iGroup; k++)if( eClass[k] == eClass[j]) eClass[k] = j;
		eClass[j] = j;
		}
	    }
	/* Update particle group ids: */ 
	for (pn=0; pn<nTree ; pn++) {
	    p[pn].pGroup = eClass[p[pn].pGroup]; 
	    }	
	/* Update RM list: */ 
	rmTmp = (FOFRM *)malloc(pkd->nMaxRm*sizeof(FOFRM)); 
	assert(rmTmp != NULL); 
	for (i=iGroup; i > 0 ;i--) { 
	    if(eClass[i] != i && protoGroup[i].nRemoteMembers > 0) { 
			
		/*copy these RM to tmp */
		for (k=0;k < protoGroup[i].nRemoteMembers ;k++){  
		    rmTmp[k].iPid = rm[protoGroup[i].iFirstRm+k].iPid;	
		    rmTmp[k].iIndex = rm[protoGroup[i].iFirstRm+k].iIndex;	
		    }
		/* move all from protoGroup[eClass[i] + 1].iFirstRm */
		/* to protoGroup[i].iFirstRM by protoGroup[i].nRemoteMembers upward */
		for (k = protoGroup[i].iFirstRm-1; k>= protoGroup[eClass[i]+ 1].iFirstRm;k--){	 
		    rm[k + protoGroup[i].nRemoteMembers].iPid = rm[k].iPid;	
		    rm[k + protoGroup[i].nRemoteMembers].iIndex  = rm[k].iIndex;	
		    } 
		/* copy tmp to protoGroup[eClass[i] +1].iFirstRM */ 
		for (k=0;k < protoGroup[i].nRemoteMembers ;k++){ 
		    rm[protoGroup[eClass[i] + 1].iFirstRm+k] = rmTmp[k] ;				 
		    }	 
		/* update iFirstRm where neccesary */
		for (j= eClass[i] + 1 ; j <= i ;j++) {
		    protoGroup[j].iFirstRm += protoGroup[i].nRemoteMembers; 
		    }
		protoGroup[eClass[i]].nRemoteMembers += protoGroup[i].nRemoteMembers; 
		protoGroup[i].nRemoteMembers = 0; 
		} 
	    }	 
	free(rmTmp); 
	free(eClass); 
	}
    /*
    ** Now we can already reject small groups if they are local 
    */
    for (pn=0; pn<nTree ; pn++){
	++(protoGroup[p[pn].pGroup].nMembers);
	} 	
    /*
    ** Create a remapping and give unique local Ids !
    */
    iMaxGroups = iGroup;
    iGroup= 1+ pkd->idSelf;
    pkd->nGroups = 0;
    protoGroup[0].iId = tmp;
    for (i=1;i <= iMaxGroups ;i++) {
	protoGroup[i].iId = iGroup;
	if(protoGroup[i].nMembers < smf->nMinMembers && protoGroup[i].nRemoteMembers == 0 ) {
	    protoGroup[i].iId = tmp;
	    }	
	else {
	    iGroup += pkd->nThreads;
	    ++pkd->nGroups;
	    }		
	}	 
    /*
    ** Update the particle groups ids.
    */
    for (pi=0; pi<nTree ; pi++) {
	p[pi].pGroup = protoGroup[p[pi].pGroup].iId;
	}
    /*
    ** Allocate memory for group data
    */
    if(pkd->nGroups>0){
	pkd->groupData = (FOFGD *) malloc(pkd->nGroups*sizeof(FOFGD)); 
	assert(pkd->groupData != NULL);
	}
    k=1;
    for (i=0 ; i < pkd->nGroups;i++){
	while(protoGroup[k].iId == tmp) k++; 
	pkd->groupData[i].iGlobalId = protoGroup[k].iId;
	pkd->groupData[i].iLocalId = protoGroup[k].iId;
	pkd->groupData[i].nLocal = protoGroup[k].nMembers;
	pkd->groupData[i].iFirstRm = protoGroup[k].iFirstRm;
	pkd->groupData[i].nRemoteMembers = protoGroup[k].nRemoteMembers;
	k++;
	pkd->groupData[i].bMyGroup = 1;
	pkd->groupData[i].fMass = 0.0;
	pkd->groupData[i].fStarMass = 0.0;
	pkd->groupData[i].fGasMass = 0.0;
	pkd->groupData[i].fAvgDens = 0.0;
	pkd->groupData[i].fVelDisp = 0.0;
	for(j=0;j<3;j++){		
	    pkd->groupData[i].fVelSigma2[j] = 0.0;
	    pkd->groupData[i].r[j] = 0.0;
	    pkd->groupData[i].v[j] = 0.0;
	    pkd->groupData[i].rmax[j] = -l[j];
	    pkd->groupData[i].rmin[j] = l[j];
	    }
	pkd->groupData[i].fDeltaR2 = 0.0;
	pkd->groupData[i].potmin = -1.0;
	pkd->groupData[i].vcircMax = 0.0;
	pkd->groupData[i].rvcircMax = 0.0;
	pkd->groupData[i].rvir = 0.0;
	pkd->groupData[i].Mvir = 0.0;
	pkd->groupData[i].rhoBG = 0.0;
	pkd->groupData[i].lambda = 0.0;
	if(pkd->groupData[i].nRemoteMembers == 0) cnt++;
	}
	
    rm  = (FOFRM *) realloc(rm,(iRmIndex)*sizeof(FOFRM)); 
    pkd->remoteMember = rm; 
    pkd->nRm = iRmIndex;
    /*
    ** Calculate local group properties
    */
    for (pi=0;pi<nTree ;++pi) {
	if(p[pi].pGroup != tmp){	
	    i = (p[pi].pGroup - 1 - pkd->idSelf)/pkd->nThreads;
	    if(TYPETest(&p[pi],TYPE_GAS) ) 
		pkd->groupData[i].fGasMass += p[pi].fMass;
	    if(TYPETest(&p[pi],TYPE_STAR) ) 
		pkd->groupData[i].fStarMass += p[pi].fMass;
	    pkd->groupData[i].fVelDisp += (p[pi].v[0]*p[pi].v[0]+p[pi].v[1]*p[pi].v[1]+p[pi].v[2]*p[pi].v[2])*p[pi].fMass;
	    for(j=0;j<3;j++){
		pkd->groupData[i].fVelSigma2[j] += ( p[pi].v[j]*p[pi].v[j] )*p[pi].fMass;
		if(pkd->groupData[i].fMass != 0.0) 
		    r[j] = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, p[pi].r[j], l[j]);
		else  r[j] = p[pi].r[j]; 
		pkd->groupData[i].r[j] += r[j]*p[pi].fMass;
		pkd->groupData[i].fDeltaR2 +=  r[j]* r[j]*p[pi].fMass;
		if(r[j] > pkd->groupData[i].rmax[j]) pkd->groupData[i].rmax[j] = r[j];
		if(r[j] < pkd->groupData[i].rmin[j]) pkd->groupData[i].rmin[j] = r[j];
		pkd->groupData[i].v[j] += p[pi].v[j]*p[pi].fMass;
		}
	    pkd->groupData[i].fMass += p[pi].fMass;
	    /* use absolute values of particle potential, sign of pot has changed in pkdgrav2*/
	    if(fabs(p[pi].fPot) > pkd->groupData[i].potmin){
		pkd->groupData[i].potmin = fabs(p[pi].fPot);
		for(j=0;j<3;j++)pkd->groupData[i].rpotmin[j]=r[j];
		}
	    if(nFOFsDone > 0){
		rho = p[pn].fMass/(pow(p[pn].fBall ,3.0/2.0)*smf->fContrast);
		if(rho > pkd->groupData[i].rhoBG)pkd->groupData[i].rhoBG = rho;
		} else { 
		    pkd->groupData[i].rhoBG = 1.0;
		    }
	    }
	}
    for (i=0; i<pkd->nGroups;++i){
	pkd->groupData[i].fRadius = 0.0;
	pkd->groupData[i].fRadius += pkd->groupData[i].rmax[0]-pkd->groupData[i].rmin[0];
	pkd->groupData[i].fRadius += pkd->groupData[i].rmax[1]-pkd->groupData[i].rmin[1];
	pkd->groupData[i].fRadius += pkd->groupData[i].rmax[2]-pkd->groupData[i].rmin[2];
	pkd->groupData[i].fRadius /= 6.0;
	if(pkd->groupData[i].nLocal > 1 ){
	    pkd->groupData[i].fAvgDens = pkd->groupData[i].fMass*0.238732414/pkd->groupData[i].fRadius/pkd->groupData[i].fRadius/pkd->groupData[i].fRadius; /*assuming spherical*/		
	    } else
		pkd->groupData[i].fAvgDens = 0.0;
	}	
    }  

FLOAT corrPos(FLOAT com, FLOAT r,FLOAT l){
    if(com > 0.2*l && r< -0.2*l) return r + l;
    else if(com < - 0.2*l && r > 0.2*l) return r - l;
    else return r;
    }
int CmpParticleGroupIds(const void *v1,const void *v2)
    {
    PARTICLE *p1 = (PARTICLE *)v1;
    PARTICLE *p2 = (PARTICLE *)v2;
    return(p1->pGroup - p2->pGroup);
    }
FLOAT phase_dist(PARTICLE pa,PARTICLE pb)
    { 
    int i;
    FLOAT dx;

    dx=0.0;
    for(i=0; i<3; i++) dx = pow(pa.r[i]-pb.r[i],2);
    dx /= pa.fBall;
    return dx;
    /* UNCOMMENT this to use phase space FOF for the adaptive FOF runs:  --J.D.--
    ** FLOAT dv=0.0; 
    ** for(i=0; i<3; i++) dv = pow(pa.v[i]-pb.v[i],2); 
    ** dv /= pa.fBallv2; 
    ** return dx+dv; */ 
    }
int CmpRMs(const void *v1,const void *v2)
    {
    FOFRM *rm1 = (FOFRM *)v1;
    FOFRM *rm2 = (FOFRM *)v2;
    if(rm1->iPid != rm1->iPid) return(rm1->iPid - rm2->iPid);
    else return(rm1->iIndex - rm2->iIndex);
    }
int CmpProtoGroups(const void *v1,const void *v2)
    { 
    FOFPG *g1 = (FOFPG *)v1;
    FOFPG *g2 = (FOFPG *)v2;
    return(g1->iId - g2->iId);
    }
int CmpGroups(const void *v1,const void *v2)
    {
    FOFGD *g1 = (FOFGD *)v1;
    FOFGD *g2 = (FOFGD *)v2;
    return(g1->nTotal - g2->nTotal);
    }
int smGroupMerge(SMF *smf,int bPeriodic)
    {
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p = smf->pkd->pStore;
    PARTICLE *pPart;
    FLOAT l[3], r,min,max,corr;
    int pi,id,i,j,k,index,listSize, sgListSize, lsgListSize;
    int nLSubGroups,nSubGroups,nMyGroups;
    int iHead, iTail, nFifo,tmp, nTree;

    FOFGD *sG;
    FOFRM *rmFifo;
    FOFRM *remoteRM;
    FOFRM rm;
    FOFGD **subGroup; /* Array of group data pointers */
    FOFGD **lSubGroup; /* Array of group data pointers */

    if (bPeriodic) {
	for(j=0;j<3;j++) l[j] = pkd->fPeriod[j];
	}
    else {
	for(j=0;j<3;j++) l[j] = FLOAT_MAXVAL;
	}

    tmp = pkd->nDark+pkd->nGas+pkd->nStar;
    nTree = pkd->kdNodes[ROOT].pUpper + 1;
    nFifo = 30*pkd->nMaxRm + 1;
    sgListSize = 10*pkd->nThreads;
    lsgListSize = 10*pkd->nThreads;
	
    subGroup = (FOFGD **)malloc(sgListSize*sizeof(FOFGD *));
    assert(subGroup != NULL);

    lSubGroup = (FOFGD **)malloc(lsgListSize*sizeof(FOFGD *));
    assert(lSubGroup != NULL);
    rmFifo = (FOFRM *)malloc(nFifo*sizeof(FOFRM));	
    assert(rmFifo != NULL);	
    /*		
    ** Start RO particle cache.
    */
    mdlROcache(mdl,CID_PARTICLE,p,sizeof(PARTICLE), nTree);
    /*
    ** Start CO group data cache.
    */
    mdlCOcache(mdl,CID_GROUP,pkd->groupData,sizeof(FOFGD), pkd->nGroups,initGroupMerge,combGroupMerge);
    /*			
    ** Start RO remote member cache.
    */
    mdlROcache(mdl,CID_RM,pkd->remoteMember,sizeof(FOFRM),pkd->nRm);
	
    for (i=0; i < pkd->nGroups ;i++){
	nSubGroups = 0;
	nLSubGroups = 0;
	iHead = 0;
	iTail = 0;
	pkd->groupData[i].nTotal = pkd->groupData[i].nLocal;   
	if(pkd->groupData[i].bMyGroup == 0)	goto NextGroup;
	if(pkd->groupData[i].nRemoteMembers && pkd->groupData[i].bMyGroup){		
			
	    /* Add all remote members to the Fifo: */		
	    for (j=pkd->groupData[i].iFirstRm; j < pkd->groupData[i].iFirstRm + pkd->groupData[i].nRemoteMembers ;j++)
		rmFifo[iTail++] = pkd->remoteMember[j];
	    while(iHead != iTail){
		rm = rmFifo[iHead++];
		if(iHead == nFifo) iHead = 0;
		if(rm.iPid == pkd->idSelf){
		    pPart = pkd->pStore + rm.iIndex;
		    /* Local: Have I got this group already? If RM not in a group, ignore it */
		    if(pPart->pGroup == pkd->groupData[i].iLocalId || pPart->pGroup == tmp
		       || pPart->pGroup == 0 ) goto NextRemoteMember;	
		    for (k=0; k < nLSubGroups ;k++){
			if(lSubGroup[k]->iLocalId == pPart->pGroup){
			    goto NextRemoteMember;
			    }
			}
		    /* Local: New subgroup found, add to list: */	
		    sG = pkd->groupData + (pPart->pGroup - 1 - pkd->idSelf)/pkd->nThreads;
		    if(nLSubGroups >= lsgListSize){
			lsgListSize *= 1.5;
			lSubGroup = (FOFGD **)realloc(lSubGroup, lsgListSize*sizeof(FOFGD *));
			assert(lSubGroup != NULL);						
			}
		    lSubGroup[nLSubGroups++] = sG;
		    if(sG->fAvgDens > pkd->groupData[i].fAvgDens) {
			pkd->groupData[i].bMyGroup = 0;
			goto NextGroup;	
			}
		    if(sG->iLocalId < pkd->groupData[i].iGlobalId)pkd->groupData[i].iGlobalId = sG->iLocalId; 
		    pkd->groupData[i].nTotal += sG->nLocal;
		    /* Add all its remote members to the Fifo: */		
		    for (j=sG->iFirstRm; j< sG->iFirstRm + sG->nRemoteMembers ;j++){					
			rmFifo[iTail++] = pkd->remoteMember[j];
			if(iTail == nFifo) iTail = 0;
			}
		    } else {
			pPart = mdlAquire(mdl,CID_PARTICLE,rm.iIndex,rm.iPid);					
			mdlRelease(mdl,CID_PARTICLE,pPart);
			/* Remote: ignore if not in a group */	
			if(pPart->pGroup == tmp){
			    goto NextRemoteMember;
			    }
			/* Remote: Have I got this group already? */	
			for (k=0; k < nSubGroups ;++k){
			    if(pPart->pGroup == subGroup[k]->iLocalId){
				goto NextRemoteMember;						
				}
			    }	
			/* Remote: New subgroup found, add to list: */	
			index = (pPart->pGroup - 1 - rm.iPid)/pkd->nThreads ;
			sG = mdlAquire(mdl,CID_GROUP,index,rm.iPid);
			if(nSubGroups >= sgListSize){
			    sgListSize *= 1.5;
			    subGroup = (FOFGD **)realloc(subGroup, sgListSize*sizeof(FOFGD *));
			    assert(subGroup != NULL);						
			    }
			subGroup[nSubGroups++] = sG;
			if(sG->fAvgDens > pkd->groupData[i].fAvgDens) {
			    pkd->groupData[i].bMyGroup = 0;
			    goto NextGroup;
			    }
			if(sG->iLocalId < pkd->groupData[i].iGlobalId)pkd->groupData[i].iGlobalId = sG->iLocalId; 
			pkd->groupData[i].nTotal += sG->nLocal;
			/* Add all its remote members to the Fifo: */		
			for (j=sG->iFirstRm; j < sG->iFirstRm + sG->nRemoteMembers ;j++){
			    remoteRM = mdlAquire(mdl,CID_RM,j,rm.iPid);		
			    rmFifo[iTail++] = *remoteRM;
			    if(iTail == nFifo) iTail = 0;
			    mdlRelease(mdl,CID_RM,remoteRM);
			    }
			}
	    NextRemoteMember:
		if(0){
		    }	
		}
	    if(pkd->groupData[i].nTotal < smf->nMinMembers){ 
		/* 
		** Nonlocal group too small:
		*/
		pkd->groupData[i].iGlobalId = 0;
		pkd->groupData[i].bMyGroup = 0;
		for (k=0;k < nSubGroups;++k) {
		    subGroup[k]->iGlobalId = 0;
		    subGroup[k]->bMyGroup = 0;
		    }	
		for (k=0;k < nLSubGroups;++k) {
		    lSubGroup[k]->iGlobalId = 0;
		    lSubGroup[k]->bMyGroup = 0;
		    }				
		} else {		
		    /*
		    ** Nonlocal group big enough: calculate properties
		    */
		    for (k=0;k < nSubGroups;++k) {
			sG = subGroup[k];
			sG->iGlobalId = pkd->groupData[i].iGlobalId;			 
			sG->bMyGroup = 0;
			if( sG->potmin > pkd->groupData[i].potmin){
			    pkd->groupData[i].potmin = sG->potmin;
			    for(j=0;j<3;j++)pkd->groupData[i].rpotmin[j]=
						corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->rpotmin[j], l[j]);
			    }
			for(j=0;j<3;j++){
			    r = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->r[j]/sG->fMass, l[j]);
			    pkd->groupData[i].r[j] += r*sG->fMass;
			    max = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->rmax[j], l[j]);
			    if(max > pkd->groupData[i].rmax[j]) pkd->groupData[i].rmax[j] = max;
			    min = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->rmin[j], l[j]);
			    if(min < pkd->groupData[i].rmin[j]) pkd->groupData[i].rmin[j] = min;
			    pkd->groupData[i].fVelSigma2[j] += sG->fVelSigma2[j];
			    pkd->groupData[i].v[j] += sG->v[j];
			    corr = r - sG->r[j]/sG->fMass;
			    pkd->groupData[i].fDeltaR2 += 2.0*corr*sG->r[j] + sG->fMass*corr*corr;
			    }
			pkd->groupData[i].fStarMass += sG->fStarMass;
			pkd->groupData[i].fGasMass += sG->fGasMass;
			pkd->groupData[i].fVelDisp += sG->fVelDisp;
			pkd->groupData[i].fDeltaR2 += sG->fDeltaR2;
			pkd->groupData[i].fMass += sG->fMass;				
			}	
		    for(k=0;k < nLSubGroups;++k) {

			sG = lSubGroup[k];
			sG->iGlobalId = pkd->groupData[i].iGlobalId;			
			sG->bMyGroup = 0;
			if( sG->potmin > pkd->groupData[i].potmin){
			    pkd->groupData[i].potmin = sG->potmin;
			    for(j=0;j<3;j++)pkd->groupData[i].rpotmin[j]=
						corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->rpotmin[j], l[j]);
			    }
			for(j=0;j<3;j++){
			    r = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->r[j]/sG->fMass, l[j]);
			    pkd->groupData[i].r[j] += r*sG->fMass;
			    max = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->rmax[j], l[j]);
			    if(max > pkd->groupData[i].rmax[j]) pkd->groupData[i].rmax[j] = max;
			    min = corrPos(pkd->groupData[i].r[j]/pkd->groupData[i].fMass, sG->rmin[j], l[j]);
			    if(min < pkd->groupData[i].rmin[j]) pkd->groupData[i].rmin[j] = min;
			    pkd->groupData[i].fVelSigma2[j] += sG->fVelSigma2[j];
			    pkd->groupData[i].v[j] += sG->v[j];
			    corr = r - sG->r[j]/sG->fMass;
			    pkd->groupData[i].fDeltaR2 += 2.0*corr*sG->r[j] + sG->fMass*corr*corr;
			    }
			pkd->groupData[i].fStarMass += sG->fStarMass;
			pkd->groupData[i].fGasMass += sG->fGasMass;
			pkd->groupData[i].fVelDisp += sG->fVelDisp;
			pkd->groupData[i].fDeltaR2 += sG->fDeltaR2;
			pkd->groupData[i].fMass += sG->fMass;
			}	
		    }
	NextGroup:
	    /*
	    ** Release non-local pointers.
	    */
	    for (k=0;k < nSubGroups;k++) {
		mdlRelease(mdl,CID_GROUP,subGroup[k]);
		} 
	    } 
	}
    mdlFinishCache(mdl,CID_PARTICLE);			
    mdlFinishCache(mdl,CID_GROUP);
    mdlFinishCache(mdl,CID_RM);
    free(subGroup);
    free(lSubGroup);
    free(rmFifo);
    free(pkd->remoteMember);
    pkd->nRm = 0;
    /*  Update the groups ids of the local particles */	
    for (pi=0;pi<nTree ;pi++){		
	index = (p[pi].pGroup - 1 - pkd->idSelf)/pkd->nThreads ;
	if(index >= 0 && index < pkd->nGroups )
	    p[pi].pGroup = pkd->groupData[index].iGlobalId;
	else p[pi].pGroup = 0;
	}		
    /*	
    ** Move real groups to low memory and normalize their properties.
    */	
    nMyGroups=0;
    for (i=0; i< pkd->nGroups;i++){
	if(pkd->groupData[i].bMyGroup && pkd->groupData[i].iGlobalId != 0){
	    for(j=0;j<3;j++){
		pkd->groupData[i].r[j] /= pkd->groupData[i].fMass;
		if(pkd->groupData[i].r[j] > 0.5*l[j]) pkd->groupData[i].r[j] -= l[j];
		if(pkd->groupData[i].r[j] < -0.5*l[j]) pkd->groupData[i].r[j] += l[j];	 
		pkd->groupData[i].v[j] /= pkd->groupData[i].fMass;
		pkd->groupData[i].fVelSigma2[j] = 
		    pkd->groupData[i].fVelSigma2[j]/pkd->groupData[i].fMass - pkd->groupData[i].v[j]*pkd->groupData[i].v[j];
		pkd->groupData[i].fVelSigma2[j] = pow(pkd->groupData[i].fVelSigma2[j], 0.5);
		}	
	    pkd->groupData[i].fRadius = 0.0;
	    pkd->groupData[i].fRadius += pkd->groupData[i].rmax[0]-pkd->groupData[i].rmin[0];
	    pkd->groupData[i].fRadius += pkd->groupData[i].rmax[1]-pkd->groupData[i].rmin[1];
	    pkd->groupData[i].fRadius += pkd->groupData[i].rmax[2]-pkd->groupData[i].rmin[2];
	    pkd->groupData[i].fRadius /= 6.0;
	    pkd->groupData[i].fVelDisp = 
		pkd->groupData[i].fVelDisp/pkd->groupData[i].fMass - pkd->groupData[i].v[0]*pkd->groupData[i].v[0]
		- pkd->groupData[i].v[1]*pkd->groupData[i].v[1]- pkd->groupData[i].v[2]*pkd->groupData[i].v[2];
	    pkd->groupData[i].fVelDisp = pow(pkd->groupData[i].fVelDisp, 0.5);
	    pkd->groupData[i].fDeltaR2 /= pkd->groupData[i].fMass; 
	    pkd->groupData[i].fDeltaR2 -= pkd->groupData[i].r[0]*pkd->groupData[i].r[0]
		+ pkd->groupData[i].r[1]*pkd->groupData[i].r[1] + pkd->groupData[i].r[2]*pkd->groupData[i].r[2];
	    pkd->groupData[i].fDeltaR2 = pow(pkd->groupData[i].fDeltaR2, 0.5);
	    pkd->groupData[i].fAvgDens = 
		pkd->groupData[i].fMass*0.238732414/pkd->groupData[i].fRadius/pkd->groupData[i].fRadius/pkd->groupData[i].fRadius; /*assuming spherical halos*/
	    pkd->groupData[nMyGroups] = pkd->groupData[i];
	    nMyGroups++;	
	    }
	}	
    if(nMyGroups == pkd->nGroups){
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,(nMyGroups+1)*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	}
    pkd->groupData[nMyGroups].bMyGroup = 0;
    /*
    ** Start RO group data cache and master reads and saves all the group data.
    */
    mdlROcache(mdl,CID_GROUP,pkd->groupData,sizeof(FOFGD), nMyGroups + 1);
    if(pkd->idSelf == 0) {
	listSize = pkd->nThreads*(pkd->nGroups+1);
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,listSize*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	for(id=1; id < pkd->nThreads; id++){
	    index = 0;
	    while(1){
		sG = mdlAquire(mdl,CID_GROUP,index,id);
		mdlRelease(mdl,CID_GROUP,sG);
		if(sG->bMyGroup != 0){
		    if(nMyGroups >= listSize-1){
			listSize *= 1.5; 
			pkd->groupData = (FOFGD *) realloc(pkd->groupData, listSize*sizeof(FOFGD)); 
			assert(pkd->groupData != NULL);
			}
		    pkd->groupData[nMyGroups] = *sG;
		    index++;
		    nMyGroups++;
		    } else {
			break;
			}
		}
	    }
	pkd->groupData[nMyGroups].bMyGroup = 0;
	}
    mdlFinishCache(mdl,CID_GROUP);
    if(pkd->idSelf != 0){
	nMyGroups = 0;
	} else {
	    /*
	    ** Master orders groupData according to the number of members
	    */	
	    qsort(pkd->groupData,nMyGroups,sizeof(FOFGD), CmpGroups);
	    }
    pkd->nGroups = nMyGroups;
    return nMyGroups;
    }

int smGroupProfiles(SMX smx, SMF *smf,int bPeriodic, int nTotalGroups,int bLogBins, int nFOFsDone)
    {
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p = smf->pkd->pStore;
    double dx2;
    FLOAT l[3],L[3],r[3],relvel[3],com[3],V,Vprev,vcirc,vcircMax,rvcircMax,M,R,binFactor;
    FLOAT rvir,Mvir,Delta,fBall;
    int pn,i,j,k,iBin,nBins,maxId,nTree,index,nCnt,pnn;
    int* iGroupIndex;
    int iStart[3],iEnd[3];
    int ix,iy,iz;
    FOFGD *gdp;
    FOFBIN *bin;
    Delta = smf->Delta; /* density contrast over critial density within rvir  */
    binFactor = smf->binFactor; /* to assure that the bins always reach out to the tidal/virial radius */
    M=0.0;
    R=0.0;
    if (bPeriodic) {
	for(j=0;j<3;j++) l[j] = pkd->fPeriod[j];
	}
    else {
	for(j=0;j<3;j++) l[j] = FLOAT_MAXVAL;
	}
    nTree = pkd->kdNodes[ROOT].pUpper + 1;
    /*
    ** Start RO group data cache and read all if you are not master
    */
    if(pkd->idSelf != 0){
	pkd->groupData = (FOFGD *) malloc(nTotalGroups*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	}
    mdlROcache(mdl,CID_GROUP,pkd->groupData,sizeof(FOFGD),pkd->nGroups);
    if(pkd->idSelf != 0){
	for(i=0; i< nTotalGroups; i++){
	    gdp = mdlAquire(mdl,CID_GROUP,i,0);
	    mdlRelease(mdl,CID_GROUP,gdp);
	    pkd->groupData[i] = *gdp;			
	    }
	}	
    mdlFinishCache(mdl,CID_GROUP);
 	 
    /*
    ** Calculate nBins and allocate bin memory.
    */
    maxId = 0;
    nBins = 0;
    for(i=0; i< nTotalGroups; i++){
	if( pkd->groupData[i].nTotal >= smf->nMinProfile ){ /* use RM field for bin pointers now...*/ 
	    pkd->groupData[i].nRemoteMembers = smf->nBins; 
	    }	else pkd->groupData[i].nRemoteMembers = 0;
	pkd->groupData[i].iFirstRm = nBins;
	nBins += pkd->groupData[i].nRemoteMembers;
	if(pkd->groupData[i].iGlobalId > maxId ) maxId = pkd->groupData[i].iGlobalId;
	}
    pkd->groupBin = (FOFBIN *) malloc(nBins*sizeof(FOFBIN));
    assert(pkd->groupBin != NULL);
    /*
    ** Create a map group id -> index of group in the pkd->groupData array.
    */
    iGroupIndex = (int *) malloc( (maxId+1)*sizeof(int));
    assert(iGroupIndex != NULL);
    for(i=0; i< maxId+1; i++)iGroupIndex[i]= -1;
    for(i=0; i< nTotalGroups; i++){
	iGroupIndex[pkd->groupData[i].iGlobalId] = i; 
	}
    /*			
    ** Initalize bin array
    */	
    iBin = 0;
    for (i=0; i< nTotalGroups; i++) {
	for(j=0; j < pkd->groupData[i].nRemoteMembers; j++){
	    assert(iBin < nBins);
	    pkd->groupBin[iBin].iId = pkd->groupData[i].iGlobalId;
	    pkd->groupBin[iBin].nMembers = 0;
	    if(bLogBins){/* logarithmic bins */
		dx2 = pkd->groupData[i].nRemoteMembers-(j+1);
		pkd->groupBin[iBin].fRadius = pkd->groupData[i].fRadius*binFactor* 
		    pow(1e-4, dx2/pkd->groupData[i].nRemoteMembers);
		} else { /* linear bins */
		    dx2 = j+1;
		    pkd->groupBin[iBin].fRadius = pkd->groupData[i].fRadius*binFactor*  
			dx2/pkd->groupData[i].nRemoteMembers;  
		    } 
	    pkd->groupBin[iBin].fMassInBin = 0.0;
	    pkd->groupBin[iBin].fMassEnclosed = 0.0;
	    pkd->groupBin[iBin].v2[0] = 0.0;
	    pkd->groupBin[iBin].v2[1] = 0.0;
	    pkd->groupBin[iBin].v2[2] = 0.0;
	    pkd->groupBin[iBin].L[0] = 0.0;
	    pkd->groupBin[iBin].L[1] = 0.0;
	    pkd->groupBin[iBin].L[2] = 0.0;	
	    pkd->groupBin[iBin].a = 0.0;	
	    pkd->groupBin[iBin].b = 0.0;	
	    pkd->groupBin[iBin].c = 0.0;	
	    pkd->groupBin[iBin].phi = 0.0;	
	    pkd->groupBin[iBin].theta = 0.0;	
	    pkd->groupBin[iBin].psi = 0.0;	
	    iBin++;
	    } 
	}
    /*			
    ** Add my particles to the correspondig bin
    */	
    for(pn=0;pn<nTree;pn++) {
	p[pn].pBin = -1;
	}
    for (index=0; index< nTotalGroups; index++){ 
	if(pkd->groupData[index].nRemoteMembers > 0){
	    k = pkd->groupData[index].iFirstRm + pkd->groupData[index].nRemoteMembers -1;
	    for(j = 0; j < 3; j++){
		if(nFOFsDone >= smf->bUsePotmin){	
		    com[j] =  pkd->groupData[index].r[j];	
		    } else {		
			com[j] = pkd->groupData[index].rpotmin[j];						
			}	
		}
	    smx->nnListSize = 0;
	    fBall =pkd->groupBin[k].fRadius;
	    if (bPeriodic) {
		for (j=0;j<3;++j) {
		    iStart[j] = floor((com[j] - fBall)/pkd->fPeriod[j] + 0.5);
		    iEnd[j] = floor((com[j] + fBall)/pkd->fPeriod[j] + 0.5);
		    }
		for (ix=iStart[0];ix<=iEnd[0];++ix) {
		    r[0] = com[0] - ix*pkd->fPeriod[0];
		    for (iy=iStart[1];iy<=iEnd[1];++iy) {
			r[1] = com[1] - iy*pkd->fPeriod[1];
			for (iz=iStart[2];iz<=iEnd[2];++iz) {
			    r[2] = com[2] - iz*pkd->fPeriod[2];
			    smGatherLocal(smx,fBall*fBall,r);
			    }
			}
		    }
		}
	    else {
		smGatherLocal(smx,fBall*fBall,com);
		}
	    nCnt = smx->nnListSize;
	    for(pnn=0;pnn<nCnt;++pnn ){
		dx2=0.0;
		for(j = 0; j < 3; j++){
		    r[j] = corrPos(com[j], smx->nnList[pnn].pPart->r[j], l[j])	- com[j];
		    relvel[j] = smx->nnList[pnn].pPart->v[j] - pkd->groupData[index].v[j];
		    dx2 += pow(r[j],2);
		    } 
		iBin = pkd->groupData[index].iFirstRm;
		k = 0;
		while(dx2 > pkd->groupBin[iBin+k].fRadius*pkd->groupBin[iBin+k].fRadius){
		    k++;
		    if( k == pkd->groupData[index].nRemoteMembers) goto nextParticle;
		    }
		assert(iBin+k < nBins);
		pkd->groupBin[iBin+k].nMembers++;
		pkd->groupBin[iBin+k].fMassInBin += smx->nnList[pnn].pPart->fMass;
		pkd->groupBin[iBin+k].fMassEnclosed += smx->nnList[pnn].pPart->fMass;
		pkd->groupBin[iBin+k].v2[0] += smx->nnList[pnn].pPart->fMass*pow(relvel[0],2.0);
		pkd->groupBin[iBin+k].v2[1] += smx->nnList[pnn].pPart->fMass*pow(relvel[1],2.0);
		pkd->groupBin[iBin+k].v2[2] += smx->nnList[pnn].pPart->fMass*pow(relvel[2],2.0);
		pkd->groupBin[iBin+k].com[0] = com[0];
		pkd->groupBin[iBin+k].com[1] = com[1];
		pkd->groupBin[iBin+k].com[2] = com[2];
		pkd->groupBin[iBin+k].L[0] += smx->nnList[pnn].pPart->fMass*(r[1]*relvel[2] - r[2]*relvel[1]);
		pkd->groupBin[iBin+k].L[1] += smx->nnList[pnn].pPart->fMass*(r[2]*relvel[0] - r[0]*relvel[2]);
		pkd->groupBin[iBin+k].L[2] += smx->nnList[pnn].pPart->fMass*(r[0]*relvel[1] - r[1]*relvel[0]);	
		smx->nnList[pnn].pPart->pBin = iBin+k;
		while(k < pkd->groupData[index].nRemoteMembers-1){
		    k++;
		    pkd->groupBin[iBin+k].fMassEnclosed += smx->nnList[pnn].pPart->fMass;
		    }
	    nextParticle:
		;
		}		
	    }
	}
    /*			
    ** Start CO group profiles cache.
    */	 
    mdlCOcache(mdl,CID_BIN,pkd->groupBin,sizeof(FOFBIN),nBins,initGroupBins,combGroupBins);
    if(pkd->idSelf != 0) { 
	for(i=0; i< nBins; i++){ 
	    bin = mdlAquire(mdl,CID_BIN,i,0);	   
	    *bin = pkd->groupBin[i];  	
	    mdlRelease(mdl,CID_BIN,bin);      	
	    }
	}
    mdlFinishCache(mdl,CID_BIN);
    /*			
    ** Caculate densities, vcirc max, rvir:
    */	
    if(pkd->idSelf == 0){
	Vprev = 0.0;
	vcircMax = 0.0;
	rvcircMax = 0.0;
	rvir = 0.0;
	Mvir = 0.0;
	for(k=0;k<3;k++) L[k] = 0.0;
	for (i=0; i< nBins; i++) {
	    if( i > 0 && pkd->groupBin[i].iId != j ){
		Vprev = 0.0;
		pkd->groupData[iGroupIndex[j]].vcircMax = vcircMax;
		pkd->groupData[iGroupIndex[j]].rvcircMax = rvcircMax;
		pkd->groupData[iGroupIndex[j]].rvir = rvir;
		pkd->groupData[iGroupIndex[j]].Mvir = Mvir;
		pkd->groupData[iGroupIndex[j]].lambda = 
		    pow((L[0]*L[0]+L[1]*L[1]+L[2]*L[2])/2.0,0.5)/(M*pow(M*R,0.5) );
		vcircMax = 0.0;
		rvcircMax = 0.0;
		rvir = 0.0;
		Mvir = 0.0;
		for(k=0;k<3;k++) L[k] = 0.0;
		}
	    j = pkd->groupBin[i].iId;
	    V = 4.1887902048*pow(pkd->groupBin[i].fRadius,3.0);
	    pkd->groupBin[i].fDensity = pkd->groupBin[i].fMassInBin/(V-Vprev);
	    for(k=0;k<3;k++)pkd->groupBin[i].v2[k] /= pkd->groupBin[i].fMassInBin;
	    Vprev = V;
	    vcirc = pow(pkd->groupBin[i].fMassEnclosed/pkd->groupBin[i].fRadius,0.5);
	    if(pkd->groupBin[i].fMassEnclosed/V > Delta){
		rvir = pkd->groupBin[i].fRadius;
		Mvir = pkd->groupBin[i].fMassEnclosed;
		}
	    if( vcirc >  vcircMax 
		&& vcirc > 5.0*2.046*pow(pkd->groupData[iGroupIndex[j]].rhoBG,0.5)*pkd->groupBin[i].fRadius){ 
		if(rvir < pkd->groupBin[i].fRadius){
		    }else{
			vcircMax = vcirc;
			rvcircMax = pkd->groupBin[i].fRadius;
			}
		}	
	    if(! (pkd->groupBin[i].fRadius > rvcircMax) ){
		for(k=0;k<3;k++) L[k] += pkd->groupBin[i].L[k];
		M = pkd->groupBin[i].fMassEnclosed;
		R = pkd->groupBin[i].fRadius;
		}
	    }
	pkd->groupData[iGroupIndex[j]].vcircMax = vcircMax;
	pkd->groupData[iGroupIndex[j]].rvcircMax = rvcircMax;
	pkd->groupData[iGroupIndex[j]].rvir = rvir;
	pkd->groupData[iGroupIndex[j]].Mvir = Mvir;
	pkd->groupData[iGroupIndex[j]].lambda = 
	    pow((L[0]*L[0]+L[1]*L[1]+L[2]*L[2])/2.0,0.5)/(M*pow(M*R,0.5) );
	}
    pkd->nBins =  nBins;
    return nBins;
    }	
