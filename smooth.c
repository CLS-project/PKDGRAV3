#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *smooth_module_id = "$Id$";

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
#include "rbtree.h"
#include <sys/stat.h>

int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType) {
    SMX smx;
    void (*initParticle)(void *,void *) = NULL;
    void (*init)(void *,void *) = NULL;
    void (*comb)(void *,void *,void *) = NULL;
    int pi;
    int nTree;
    int iTopDepth;

    smx = malloc(sizeof(struct smContext));
    assert(smx != NULL);
    smx->pkd = pkd;
    if (smf != NULL) smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->bPeriodic = bPeriodic;

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
    case SMX_DENDVDX:
	assert( pkd->oSph ); /* Validate memory model */
	smx->fcnSmooth = DenDVDX;
	initParticle = NULL; /* Original Particle */
	init = initDenDVDX; /* Cached copies */
	comb = combDenDVDX;
	smx->fcnPost = NULL;
	break;
    case SMX_SPHFORCES:
	assert( pkd->oSph ); /* Validate memory model */
	assert( pkd->oAcceleration ); /* Validate memory model */
	smx->fcnSmooth = SphForces;
	initParticle = initSphForcesParticle; /* Original Particle */
	init = initSphForces; /* Cached copies */
	comb = combSphForces;
	smx->fcnPost = NULL;
	break;
    case SMX_MEANVEL:
	assert( pkd->oVelSmooth); /* Validate memory model */
	smx->fcnSmooth = bSymmetric?MeanVelSym:MeanVel;
	initParticle = initMeanVel; /* Original Particle */
	init = initMeanVel; /* Cached copies */
	comb = combMeanVel;
	smx->fcnPost = NULL;
	break;
    case SMX_DIVV:
	assert( pkd->oVelSmooth); /* Validate memory model */
	smx->fcnSmooth = bSymmetric?DivvSym:Divv;
	initParticle = initDivv; /* Original Particle */
	init = initDivv; /* Cached copies */
	comb = combDivv;
	smx->fcnPost = NULL;
	break;
    case SMX_VELDISP2:
	assert( pkd->oVelSmooth); /* Validate memory model */
	smx->fcnSmooth = bSymmetric?VelDisp2Sym:VelDisp2;
	initParticle = initVelDisp2; /* Original Particle */
	init = initVelDisp2; /* Cached copies */
	comb = combVelDisp2;
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
    case SMX_RELAXATION:
	assert( pkd->oRelaxation); /* Validate memory model */
	assert(bSymmetric == 0);
	smx->fcnSmooth = AddRelaxation;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
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
	    PARTICLE *p = pkdParticle(pkd,pi);
	    /*if (TYPETest(p,smx->eParticleTypes))*/
	    if (pkdIsActive(pkd,p)) initParticle(pkd,p);
	    }
	}
    /*
    ** Start particle caching space (cell cache is already active).
    */
    if (bSymmetric) {
	mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,
		   pkdParticleBase(pkd),pkdParticleSize(pkd),
		   nTree,pkd,init,comb);
	}
    else {
	mdlROcache(pkd->mdl,CID_PARTICLE,NULL,
		   pkdParticleBase(pkd),pkdParticleSize(pkd),
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


void smFinish(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
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
	    p = pkdParticle(pkd,pi);
	    if ( pkdIsSrcActive(p,0,MAX_RUNG) && pkdIsDstActive(p,0,MAX_RUNG) )
		smx->fcnPost(p,smf);
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
PQ *pqSearchLocal(SMX smx,FLOAT r[3],int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p;
    PQ *pq;
    FLOAT dx,dy,dz,dMin,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int i,j,pj,pWant,pEnd,iCell,iParent;
    int sp = 0;
    int sm = 0;

    *pbDone = 1;	/* assume that we will complete the search */
    assert(smx->nQueue == 0);
    pq = smx->pq;
    /*
    ** We don't perform containment tests except at the
    ** root, so that the pbDone flag can be correctly
    ** set.
    */
    iCell = ROOT;
    S[sp] = iCell;
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
	    }
	pWant = c[iCell].pLower + smx->nSmooth - smx->nQueue - 1;
	pEnd = c[iCell].pUpper;
	if (pWant > pEnd) {
	    for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		pq[smx->nQueue].pPart = p;
		pq[smx->nQueue].fDist2 = dx*dx + dy*dy + dz*dz;
		++smx->nQueue;
		}
	    }
	else {
	    for (pj=c[iCell].pLower;pj<=pWant;++pj) {
		p = pkdParticle(pkd,pj);
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		pq[smx->nQueue].pPart = p;
		pq[smx->nQueue].fDist2 = dx*dx + dy*dy + dz*dz;
		++smx->nQueue;
		}
	    PQ_BUILD(pq,smx->nSmooth,pq);
	    for (;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < pq->fDist2) {
		    pq->pPart = p;
		    pq->fDist2 = fDist2;
		    PQ_REPLACE(pq);
		    }
		}
	    goto NoIntersect;  /* done loading phase */
	    }
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
	iCell ^= 1;
	if (sm) --sm;
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
	    MINDIST(c[iCell].bnd,r,min1);
	    ++iCell;
	    MINDIST(c[iCell].bnd,r,min2);
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
	pEnd = c[iCell].pUpper;
	for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
	    p = pkdParticle(pkd,pj);
	    if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	    dx = r[0] - p->r[0];
	    dy = r[1] - p->r[1];
	    dz = r[2] - p->r[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < pq->fDist2) {
		pq->pPart = p;
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
	    MINDIST(c[iCell].bnd,r,min2);
	    }
	if (min2 >= pq->fDist2) {
	    iCell = c[iCell].iParent;
	    goto NoIntersect;
	    }
	S[++sp] = iCell;
	}
    }



PQ *pqSearchRemote(SMX smx,PQ *pq,int id,FLOAT r[3]) {
    PKD pkd = smx->pkd;
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
	    }
	pWant = pkdn->pLower + smx->nSmooth - smx->nQueue - 1;
	pEnd = pkdn->pUpper;
	if (pWant > pEnd) {
	    for (pj=pkdn->pLower;pj<=pEnd;++pj) {
		if (id == idSelf) {
		    p = pkdParticle(pkd,pj);
		    pq[smx->nQueue].bRemote = 0;
		    }
		else {
		    p = mdlAquire(mdl,CID_PARTICLE,pj,id);
		    pq[smx->nQueue].bRemote = 1;
		    }
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) {
		    if (id != idSelf) mdlRelease(mdl,CID_PARTICLE,p);
		    continue;
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
		    p = pkdParticle(pkd,pj);
		    pq[smx->nQueue].bRemote = 0;
		    }
		else {
		    p = mdlAquire(mdl,CID_PARTICLE,pj,id);
		    pq[smx->nQueue].bRemote = 1;
		    }
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) {
		    if (id != idSelf) mdlRelease(mdl,CID_PARTICLE,p);
		    continue;
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
		if (id == idSelf) p = pkdParticle(pkd,pj);
		else p = mdlAquire(mdl,CID_PARTICLE,pj,id);
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) {
		    if (id != idSelf) mdlRelease(mdl,CID_PARTICLE,p);
		    continue;
		    }
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
	iCell ^= 1;
	if (id == idSelf) pkdn = &c[iCell];
	else {
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    }
	if (sm) --sm;
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
	    MINDIST(pkdn->bnd,r,min1);
	    ++iCell;
	    if (id == idSelf) pkdu = &c[iCell];
	    else pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
	    MINDIST(pkdu->bnd,r,min2);
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
	pEnd = pkdn->pUpper;
	for (pj=pkdn->pLower;pj<=pEnd;++pj) {
	    if (id == idSelf) p = pkdParticle(pkd,pj);
	    else p = mdlAquire(mdl,CID_PARTICLE,pj,id);
	    if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) {
		if (id != idSelf) mdlRelease(mdl,CID_PARTICLE,p);
		continue;
		}
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


PQ *pqSearch(SMX smx,PQ *pq,FLOAT r[3],int bReplica,int *pbDone) {
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
	    }
	id = c[iCell].pLower;	/* this is the thread id in LTT */
	if (bReplica || id != idSelf) {
	    pq = pqSearchRemote(smx,pq,id,r);
	    }
	else {
	    pq = pqSearchLocal(smx,r,pbDone);
	    if (*pbDone) return pq;	/* early exit */
	    }
	if (smx->nQueue == smx->nSmooth) goto NoIntersect;  /* done loading phase */
	while (iCell == S[sp]) {
	    if (!sp) {
		return NULL;		/* EXIT, could not load enough particles! */
		}
	    --sp;
	    iCell = c[iCell].iParent;
	    }
	iCell ^= 1;
	if (sm) --sm;
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
	    MINDIST(c[iCell].bnd,r,min1);
	    ++iCell;
	    MINDIST(c[iCell].bnd,r,min2);
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
	    MINDIST(c[iCell].bnd,r,min2);
	    }
	if (min2 >= pq->fDist2) {
	    iCell = c[iCell].iParent;
	    goto NoIntersect;
	    }
	S[++sp] = iCell;
	}
    }


void smSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    PQ *pq;
    FLOAT r[3],fBall;
    int iStart[3],iEnd[3];
    int pi,i,j,bDone;
    int ix,iy,iz;

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	pq = NULL;
	smx->nQueue = 0;
	pq = pqSearch(smx,pq,p->r,0,&bDone);
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
		iStart[j] = floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5);
		}
	    for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = p->r[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		    r[1] = p->r[1] - iy*pkd->fPeriod[1];
		    for (iz=iStart[2];iz<=iEnd[2];++iz) {
			r[2] = p->r[2] - iz*pkd->fPeriod[2];
			if (ix || iy || iz) {
			    pq = pqSearch(smx,pq,r,1,&bDone);
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

	p->fBall = sqrt(pq->fDist2);
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
	smx->fcnSmooth(p,smx->nSmooth,smx->nnList,smf);
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


void smGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[3]) {
    PKD pkd = smx->pkd;
    KDN *c = smx->pkd->kdNodes;
    PARTICLE *p;
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,nCnt,pEnd;
    int idSelf;

    idSelf = smx->pkd->idSelf;
    nCnt = smx->nnListSize;
    iCell = ROOT;
    while (1) {
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
	    pEnd = c[iCell].pUpper;
	    for (pj=c[iCell].pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    if (nCnt >= smx->nnListMax) {
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
		    smx->nnList[nCnt].pPart = p;
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


void smGatherRemote(SMX smx,FLOAT fBall2,FLOAT r[3],int id) {
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
	    pEnd = pkdn->pUpper;
	    for (pj=pkdn->pLower;pj<=pEnd;++pj) {
		pp = mdlAquire(mdl,CID_PARTICLE,pj,id);
		if ( !pkdIsSrcActive(pp,0,MAX_RUNG) ) {
		    mdlRelease(mdl,CID_PARTICLE,pp);
		    continue;
		    }
		dx = r[0] - pp->r[0];
		dy = r[1] - pp->r[1];
		dz = r[2] - pp->r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    if (nCnt >= smx->nnListMax) {
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
		else mdlRelease(mdl,CID_PARTICLE,pp);
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


void smGather(SMX smx,FLOAT fBall2,FLOAT r[3]) {
    KDN *c = smx->pkd->kdTop;
    int idSelf = smx->pkd->idSelf;
    int *S = smx->ST;
    FLOAT min2;
    int iCell,id;
    int sp = 0;

    iCell = ROOT;
    while (1) {
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

void smReSmoothOne(SMX smx,SMF *smf,void *p,FLOAT *R,FLOAT fBall) {
    PKD pkd = smx->pkd;
    FLOAT r[3];
    int iStart[3],iEnd[3];
    int i,j;
    int ix,iy,iz;

    smx->nnListSize = 0;
    /*
    ** Note for implementing SLIDING PATCH, the offsets for particles are
    ** negative here, reflecting the relative +ve offset of the simulation
    ** volume.
    */
    if (smx->bPeriodic) {
	for (j=0;j<3;++j) {
	    iStart[j] = floor((R[j] - fBall)/pkd->fPeriod[j] + 0.5);
	    iEnd[j] = floor((R[j] + fBall)/pkd->fPeriod[j] + 0.5);
	    }
	for (ix=iStart[0];ix<=iEnd[0];++ix) {
	    r[0] = R[0] - ix*pkd->fPeriod[0];
	    for (iy=iStart[1];iy<=iEnd[1];++iy) {
		r[1] = R[1] - iy*pkd->fPeriod[1];
		for (iz=iStart[2];iz<=iEnd[2];++iz) {
		    r[2] = R[2] - iz*pkd->fPeriod[2];
		    smGather(smx,fBall*fBall,r);
		    }
		}
	    }
	}
    else {
	smGather(smx,fBall*fBall,R);
	}
    /*
    ** Apply smooth funtion to the neighbor list.
    */
    smx->fcnSmooth(p,smx->nnListSize,smx->nnList,smf);
    /*
    ** Release aquired pointers.
    */
    for (i=0;i<smx->nnListSize;++i) {
	if (smx->nnbRemote[i]) {
	    mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
	    }
	}
    }

void smReSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if ( pkdIsDstActive(p,0,MAX_RUNG) && pkdIsSrcActive(p,0,MAX_RUNG) )
	    smReSmoothOne(smx,smf,p,p->r,p->fBall);
	}
    }


FLOAT phase_dist(PKD pkd,double dvTau2,PARTICLE *pa,PARTICLE *pb,double H) {
    int j;
    FLOAT dx,dv,dx2,dv2;
    double *va, *vb;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */

    va = pkdVel(pkd,pa);
    vb = pkdVel(pkd,pb);

    dx2=0.0;
    for (j=0;j<3;++j) {
	dx = pa->r[j] - pb->r[j];
	dx2 += dx*dx;
	}
    dx2 /= pa->fBall;   /* this is actually fBall2! */
    dv2 = 0.0;
    if (dvTau2 > 0) {
      for (j=0;j<3;++j) {
	dv = (va[j] - vb[j]) + H*(pa->r[j] - pb->r[j]);
	dv2 += dv*dv;
        }
      dv2 /= dvTau2;
      }
    return(dx2 + dv2);
    }


FLOAT corrPos(FLOAT com,FLOAT r,FLOAT l) {
    if (com > 0.25*l && r < -0.25*l) return r + l;
    else if (com < -0.25*l && r > 0.25*l) return r - l;
    else return r;
    }


FLOAT PutInBox(FLOAT r,FLOAT l) {
    if (r < -0.5*l) return r + l;
    else if (r > 0.5*l) return r - l;
    else return r;
    }

typedef struct {
    RB_NODE node;
    FOFRM   data;
    } RM_NODE;

typedef struct protoGroup {
    int iId;
    int nMembers;
    RB_TREE treeRemoteMembers;
    } FOFPG;

/*
**  Copy the given tree (recursively) to an array
*/
int copy_rm( FOFRM *rm,RB_NODE *node) {
    int iCount = 0;
    if ( node != NULL ) {
	RM_NODE *rmnode = (RM_NODE *)(node);
	iCount = copy_rm( rm, node->link[0] );
	rm += iCount;
	*rm++ = rmnode->data;
	iCount++;
	iCount += copy_rm( rm, node->link[1] );
	}
    return iCount;
    }

int CmpRMs(void *ctx,const void *v1,const void *v2) {
    FOFRM *rm1 = (FOFRM *)v1;
    FOFRM *rm2 = (FOFRM *)v2;
    if (rm1->iPid != rm2->iPid) return (rm1->iPid - rm2->iPid);
    else return (rm1->iIndex - rm2->iIndex);
    }

void smFof(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p;
    float *pPot;
    int32_t *pBin;
    int32_t *pGroup;
    int32_t *pPartGroup;
    double *v;
    RB_TYPE rm_type;
    FOFRM   rm_data;
    FOFPG* protoGroup;
    FOFBIN *bin;

    int pi,pn,pnn,nCnt,i,j,k;
    int nRmListSize,nRmCnt,iRmIndex;
    int nFifo, iHead, iTail, iMaxGroups, iGroup;
    int *Fifo;
    int iStart[3],iEnd[3];
    int ix,iy,iz;

    FLOAT r[3],l[3],relpos[3],lx,ly,lz,fBall,fBall2Max,rho,fvBall2;
    FLOAT fMass;
    int nTree,tmp;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */
    assert(pkd->oPotential); /* Validate memory model */

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
    tmp = pkd->nDark+pkd->nGas+pkd->nStar;

    nTree = pkd->kdNodes[ROOT].pUpper + 1;
    iMaxGroups = nTree+1;
    nFifo = nTree;
    Fifo = (int *)malloc(nFifo*sizeof(int));
    assert(Fifo != NULL);

    /*used something smaller than FOFGD here to reduce memory usage*/
    protoGroup = (FOFPG *)malloc(iMaxGroups*sizeof(FOFPG));
    assert(protoGroup != NULL);

    /* This is the "empty" group */
    iGroup = 0;
    protoGroup[iGroup].nMembers = 0;
    protoGroup[iGroup].iId = iGroup;
    protoGroup[iGroup].treeRemoteMembers = NULL;

    nRmListSize = 0;
    rb_type_create(&rm_type,sizeof(RM_NODE),0,CmpRMs,0,0);

    pkd->nGroups = 0;
    pkd->nMaxRm = 0;
    fBall2Max = 0.0;

    /* set spatial linking lenght for each particle (mass dependend when bTauAbs=0) */
    for (pn=0;pn<nTree;pn++) {
      p = pkdParticle(pkd,pn);
      fMass = pkdMass(pkd,p);
      pGroup = pkdInt32(p,pkd->oGroup);
      *pGroup = 0; 
      if (smf->bTauAbs) {
	p->fBall = smf->dTau2;
	/* enforce a real space linking length smaller than the mean particle separation at all times :*/
	if (smf->dTau2 > 0.2*pow(fMass,0.6666) )
	  p->fBall = 0.2*pow(fMass,0.6666);
      }
      else {
	p->fBall = smf->dTau2*pow(fMass,0.6666);
      }
      if (p->fBall > fBall2Max) fBall2Max = p->fBall;
    }

    /* the velocity space linking lenght is the same for all particles. Zero means: do regular 3D FOF */
    fvBall2 = smf->dVTau2;
 
   /* Have to restart particle chache, since we will need the updated p->fBall now */
    mdlFinishCache(mdl,CID_PARTICLE);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),nTree);

    /* Starting FOF search now... */
    for (pn=0;pn<nTree;pn++) {
	p = pkdParticle(pkd,pn);
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	pGroup = pkdInt32(p,pkd->oGroup);
	if (*pGroup ) continue;
	iGroup++;
	assert(iGroup < iMaxGroups);
	protoGroup[iGroup].nMembers = 0;
	protoGroup[iGroup].iId = iGroup;
	protoGroup[iGroup].treeRemoteMembers = NULL;
	nRmCnt = 0;
	/*
	** Mark particle and add it to the do-fifo
	*/
	*pGroup = iGroup;
	Fifo[iTail] = pn; iTail++;
	if (iTail == nFifo) iTail = 0;
	while (iHead != iTail) {
	    pi = Fifo[iHead];iHead++;
	    p = pkdParticle(pkd,pi);
	    if (iHead == nFifo) iHead=0;
	    /*
	    ** Do a Ball Gather at the radius p->fBall
	    */
	    smx->nnListSize =0;
	    fBall = sqrt(p->fBall);
	    if (smx->bPeriodic) {
	      for (j=0;j<3;++j) {
		iStart[j] = floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5);
	      }
	      for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = p->r[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		  r[1] = p->r[1] - iy*pkd->fPeriod[1];
		  for (iz=iStart[2];iz<=iEnd[2];++iz) {
		    r[2] = p->r[2] - iz*pkd->fPeriod[2];
		    smGather(smx,p->fBall,r);
		  }
		}
	      }
	    }
	    else {
	      smGather(smx,p->fBall,p->r);
	    }
	    nCnt = smx->nnListSize;
	    for (pnn=0;pnn<nCnt;++pnn ) {
	      if (smx->nnbRemote[pnn] == 0) { /* Local neighbors: */

		/* Do not add particles which already are in a group*/
		pPartGroup = pkdInt32(smx->nnList[pnn].pPart,pkd->oGroup);
		if (*pPartGroup) continue;
		
		if (fvBall2 > 0.0) {
		  /* Check if this particle is close enough in velocity space  */
		  if (phase_dist(pkd,fvBall2,p,smx->nnList[pnn].pPart,smf->H) > 1.0) continue;
		}
		/*
		**  Mark particle and add it to the do-fifo
		*/
		*pPartGroup = iGroup;
		Fifo[iTail] = smx->nnList[pnn].iIndex;iTail++;
		if (iTail == nFifo) iTail = 0;
	      
	      } else {	 /* Nonlocal neighbors: */
		
		/* Make remote member linking symmetric by using smaller linking length if different: */
		if (fvBall2 > 0.0) { /* Check velocity space distance */
		  if (phase_dist(pkd,fvBall2,p,smx->nnList[pnn].pPart,smf->H) > 1.0 ||
		      phase_dist(pkd,fvBall2,smx->nnList[pnn].pPart,p,smf->H) > 1.0) continue;
		}
		else { /* real space distance */
		  if (smx->nnList[pnn].fDist2 > smx->nnList[pnn].pPart->fBall) continue;
		}
		
		/* Add to remote member (RM) list if new */
		rm_data.iIndex = smx->nnList[pnn].iIndex ;
		rm_data.iPid = smx->nnList[pnn].iPid;

		if ( rb_insert(&rm_type,&protoGroup[iGroup].treeRemoteMembers,&rm_data) ) {
		  nRmCnt++;
		  nRmListSize++;
		}
	      }
	    }
	}
	if ( nRmCnt > pkd->nMaxRm ) pkd->nMaxRm = nRmCnt;
    }
    free(Fifo);
    
    /*
    ** Now we can already reject groups which are too small, if they are entirely local
    */
    for (pn=0; pn<nTree ; pn++) {
	p = pkdParticle(pkd,pn);
	pGroup = pkdInt32(p,pkd->oGroup);
	if (*pGroup >=0 && *pGroup < iMaxGroups)
	    ++(protoGroup[*pGroup].nMembers);
	else
	    printf("ERROR: idSelf=%i , p->pGroup=%i too large. iMaxGroups=%i \n",pkd->idSelf,*pGroup,iMaxGroups);
	}
    /*
    ** Create a remapping and give unique local Ids !
    */
    iMaxGroups = iGroup;
    iGroup= 1 + pkd->idSelf;
    pkd->nGroups = 0;
    protoGroup[0].iId = tmp;
    for (i=1;i<=iMaxGroups;i++) {
	protoGroup[i].iId = iGroup;
	if (protoGroup[i].nMembers < smf->nMinMembers && protoGroup[i].treeRemoteMembers == NULL) {
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
    for (pi=0;pi<nTree;pi++) {
	p = pkdParticle(pkd,pi);
	pGroup = pkdInt32(p,pkd->oGroup);
	*pGroup = protoGroup[*pGroup].iId;
	}
    /*
    ** Allocate the remote members array
    */
    pkd->nRm = nRmListSize;
    pkd->remoteMember = mdlMalloc(mdl,(pkd->nRm+1)*sizeof(FOFRM));
    iRmIndex = 0;
    /*
    ** Allocate memory for group data
    */
    pkd->groupData = (FOFGD *) malloc((1+pkd->nGroups)*sizeof(FOFGD));
    assert(pkd->groupData != NULL);
    k=1;
    for (i=0;i<pkd->nGroups;i++) {
	while (protoGroup[k].iId == tmp) k++;
	pkd->groupData[i].iGlobalId = protoGroup[k].iId;
	pkd->groupData[i].iLocalId = protoGroup[k].iId;
	pkd->groupData[i].nLocal = protoGroup[k].nMembers;
	pkd->groupData[i].iFirstRm = iRmIndex;
	pkd->groupData[i].nRemoteMembers = copy_rm(pkd->remoteMember+iRmIndex,protoGroup[k].treeRemoteMembers);
	iRmIndex += pkd->groupData[i].nRemoteMembers;
	rb_free(&rm_type, &protoGroup[k].treeRemoteMembers);
	k++;
	pkd->groupData[i].bMyGroup = 1;
	pkd->groupData[i].fMass = 0.0;
	for (j=0;j<3;j++) {
	    pkd->groupData[i].rcom[j] = 0.0;
	    pkd->groupData[i].r[j] = 0.0;
	    pkd->groupData[i].v[j] = 0.0;
	    }
	pkd->groupData[i].fRMSRadius = 0.0;
	pkd->groupData[i].potordenmax = -1.0;
    }
    free(protoGroup);
    /* Sanity check: the list size should match the number of elements copied */
    assert( iRmIndex == nRmListSize );
    /*
    ** Calculate local group properties
    */
    for (pi=0;pi<nTree;++pi) {
	p = pkdParticle(pkd,pi);
	pGroup = pkdInt32(p,pkd->oGroup);
	fMass = pkdMass(pkd,p);
	v = pkdVel(pkd,p);
	if (*pGroup != tmp) {
	    i = (*pGroup - 1 - pkd->idSelf)/pkd->nThreads;
	    for (j=0;j<3;j++) {
	      if (pkd->groupData[i].fMass > 0.0)
		    r[j] = corrPos(pkd->groupData[i].rcom[j]/pkd->groupData[i].fMass,p->r[j],l[j]);
		else  r[j] = p->r[j];
		pkd->groupData[i].rcom[j] += r[j]*fMass;
		pkd->groupData[i].fRMSRadius +=  r[j]*r[j]*fMass;
		pkd->groupData[i].v[j] += v[j]*fMass;
		}
	    pkd->groupData[i].fMass += fMass;
	    if(smf->iCenterType == 1){ /* maximum of fabs(potential) is stored in case 1 (centered on potential)*/
	      pPot = pkdPot(pkd,p);
	      if ( fabs(*pPot) > pkd->groupData[i].potordenmax) {
		pkd->groupData[i].potordenmax = fabs(*pPot);
		for (j=0;j<3;j++) pkd->groupData[i].r[j] = r[j];
	      }
	    } else {
	      if (p->fDensity > pkd->groupData[i].potordenmax) {
		pkd->groupData[i].potordenmax = p->fDensity;
		for (j=0;j<3;j++) pkd->groupData[i].r[j] = r[j];
	      }
	    }
	}
    }
}


int CmpGroups(const void *v1,const void *v2) {
    FOFGD *g1 = (FOFGD *)v1;
    FOFGD *g2 = (FOFGD *)v2;
    return g1->nTotal - g2->nTotal;
    }

static void mktmpdir( const char *dirname ) {
    struct stat s;
    if ( stat(dirname,&s) == 0 ) {
	if ( S_ISDIR(s.st_mode) )
	    return;
	}
    mkdir( dirname, 0700 );
    }

int smGroupMerge(SMF *smf,int bPeriodic) {
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p;
    PARTICLE *pPart;
    int32_t *pBin, *pGroup;
    int32_t *pPartGroup, iPartGroup;
    FLOAT l[3], r,min,max,corr;
    int pi,id,i,j,k,index,listSize, sgListSize, lsgListSize;
    int nLSubGroups,nSubGroups,nMyGroups;
    int iHead, iTail, nFifo,tmp, nTree;
    FILE * pFile; /* files for parallel output of group ids, links and densities*/
    FILE * lFile;
    FILE * dFile;
    char filename [30];

    FOFGD *sG;
    FOFRM *rmFifo;
    FOFRM *remoteRM;
    FOFRM rm;
    FOFGD **subGroup; /* Array of group data pointers */
    FOFGD **lSubGroup; /* Array of group data pointers */

    if (bPeriodic) {
	for (j=0;j<3;j++) l[j] = pkd->fPeriod[j];
	}
    else {
	for (j=0;j<3;j++) l[j] = FLOAT_MAXVAL;
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
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd), nTree);
    /*
    ** Start CO group data cache.
    */
    /*printf( "Processor %d cache has %d entries\n", mdlSelf(mdl), pkd->nGroups );*/
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD), pkd->nGroups,pkd,initGroupMerge,combGroupMerge);
    /*
    ** Start RO remote member cache.
    */
    mdlROcache(mdl,CID_RM,NULL,pkd->remoteMember,sizeof(FOFRM),pkd->nRm);

    for (i=0; i < pkd->nGroups ;i++) {
	nSubGroups = 0;
	nLSubGroups = 0;
	iHead = 0;
	iTail = 0;
	pkd->groupData[i].nTotal = pkd->groupData[i].nLocal;
	if (pkd->groupData[i].bMyGroup == 0)	goto NextGroup;
	if (pkd->groupData[i].nRemoteMembers && pkd->groupData[i].bMyGroup) {

	    /* Add all remote members to the Fifo: */
	  for (j=pkd->groupData[i].iFirstRm; j < pkd->groupData[i].iFirstRm + pkd->groupData[i].nRemoteMembers ;j++){
		rmFifo[iTail] = pkd->remoteMember[j];
		iTail++;
	  }
	    while (iHead != iTail) {
		rm = rmFifo[iHead];
		iHead++;
		if (iHead == nFifo) iHead = 0;
		if (rm.iPid == pkd->idSelf) {
		    pPart = pkdParticle(pkd,rm.iIndex);
		    pPartGroup = pkdInt32(pPart,pkd->oGroup);
		    /* Local: Have I got this group already? If RM not in a group, ignore it */
		    if (*pPartGroup == pkd->groupData[i].iLocalId || *pPartGroup == tmp
			    || *pPartGroup == 0 ) goto NextRemoteMember;
		    for (k=0; k < nLSubGroups ;k++) {
			if (lSubGroup[k]->iLocalId == *pPartGroup) {
			    goto NextRemoteMember;
			    }
			}
		    /* Local: New subgroup found, add to list: */
		    sG = pkd->groupData + (*pPartGroup - 1 - pkd->idSelf)/pkd->nThreads;
		    if (nLSubGroups >= lsgListSize) {
			lsgListSize *= 1.5;
			lSubGroup = (FOFGD **)realloc(lSubGroup, lsgListSize*sizeof(FOFGD *));
			assert(lSubGroup != NULL);
			}
		    lSubGroup[nLSubGroups++] = sG;
		    if (sG->potordenmax > pkd->groupData[i].potordenmax) {
			pkd->groupData[i].bMyGroup = 0;
			goto NextGroup;
			}
		    if (sG->iLocalId < pkd->groupData[i].iGlobalId)
		      pkd->groupData[i].iGlobalId = sG->iLocalId;
		    pkd->groupData[i].nTotal += sG->nLocal;
		    /* Add all its remote members to the Fifo: */
		    for (j=sG->iFirstRm; j< sG->iFirstRm + sG->nRemoteMembers ;j++) {
			rmFifo[iTail++] = pkd->remoteMember[j];
			if (iTail == nFifo) iTail = 0;
			}
		    }
		else {
		    pPart = mdlAquire(mdl,CID_PARTICLE,rm.iIndex,rm.iPid);
		    iPartGroup = * pkdInt32(pPart,pkd->oGroup);
		    mdlRelease(mdl,CID_PARTICLE,pPart);
		    
		    /* Remote: ignore if not in a group */
		    if (iPartGroup == tmp) {
			goto NextRemoteMember;
			}
		    /* Remote: Have I got this group already? */
		    for (k=0; k < nSubGroups ;++k) {
			if (iPartGroup == subGroup[k]->iLocalId) {
			    goto NextRemoteMember;
			    }
			}
		    /* Remote: New subgroup found, add to list: */
		    index = (iPartGroup - 1 - rm.iPid)/pkd->nThreads ;
		    sG = mdlAquire(mdl,CID_GROUP,index,rm.iPid);
		    
		    if (nSubGroups >= sgListSize) {
			sgListSize *= 1.5;
			subGroup = (FOFGD **)realloc(subGroup, sgListSize*sizeof(FOFGD *));
			assert(subGroup != NULL);
			}
		    subGroup[nSubGroups++] = sG;
		    if (sG->potordenmax > pkd->groupData[i].potordenmax) {
			pkd->groupData[i].bMyGroup = 0;
			goto NextGroup;
			}
		    if (sG->iLocalId < pkd->groupData[i].iGlobalId) 
		      pkd->groupData[i].iGlobalId = sG->iLocalId;
		    pkd->groupData[i].nTotal += sG->nLocal;
		    /* Add all its remote members to the Fifo: */
		    for (j=sG->iFirstRm; j < sG->iFirstRm + sG->nRemoteMembers ;j++) {
			remoteRM = mdlAquire(mdl,CID_RM,j,rm.iPid);
			rmFifo[iTail++] = *remoteRM;
			if (iTail == nFifo) iTail = 0;
			mdlRelease(mdl,CID_RM,remoteRM);
			}
		    }
	    NextRemoteMember:
		;
		}
	    if (pkd->groupData[i].nTotal < smf->nMinMembers) {
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
		}
	    else {
		/*
		** Nonlocal group big enough: calculate properties
		*/
		for (k=0;k<nSubGroups + nLSubGroups;++k) {
		  if(k < nSubGroups) sG = subGroup[k];
		  else sG = lSubGroup[k-nSubGroups];
		  sG->iGlobalId = pkd->groupData[i].iGlobalId;
		  sG->bMyGroup = 0;
		  for (j=0;j<3;j++) {
		    pkd->groupData[i].rcom[j] += sG->rcom[j];
		    pkd->groupData[i].v[j] += sG->v[j];
		  }
		  pkd->groupData[i].fRMSRadius += sG->fRMSRadius;
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
    mdlFree(mdl,pkd->remoteMember);
    pkd->nRm = 0;

    for (pi=0;pi<nTree ;pi++){
      p = pkdParticle(pkd,pi);
      pGroup = pkdInt32(p,pkd->oGroup);
      index = (*pGroup - 1 - pkd->idSelf)/pkd->nThreads ;
      if(index >= 0 && index < pkd->nGroups )
	 *pGroup = pkd->groupData[index].iGlobalId;
      else
	*pGroup = 0;
    }

    /*
    ** Here the embarasssingly parallel output of particle group ids,
    ** local density and group links were written. This was removed now.
    **
    ** Output of particle arrays with local density and group ids have to be added 
    ** elsewhere. Group link output is not needed, links can be generated from the
    ** group id arrays in postprocessing quite easily.
    **
    ** JD -- Feb 24, 2009
    */

    /* Move real groups to low memory and normalize their properties. */
    nMyGroups=0;
    for (i=0; i< pkd->nGroups;i++) {
	if (pkd->groupData[i].bMyGroup && pkd->groupData[i].iGlobalId != 0) {
	    for (j=0;j<3;j++)
	       pkd->groupData[i].rcom[j] /= pkd->groupData[i].fMass;
	    /* 
	    ** Do not calculate fDeltaR2 with the corrected positions!
	    */
	    pkd->groupData[i].fRMSRadius /= pkd->groupData[i].fMass;
	    pkd->groupData[i].fRMSRadius -= pkd->groupData[i].rcom[0]*pkd->groupData[i].rcom[0]
		+ pkd->groupData[i].rcom[1]*pkd->groupData[i].rcom[1]
		+ pkd->groupData[i].rcom[2]*pkd->groupData[i].rcom[2];
	    pkd->groupData[i].fRMSRadius = sqrt(pkd->groupData[i].fRMSRadius);
	    /* 
	    ** Now put all the positions back into the box and normalise the rest
	    */
	    for (j=0;j<3;j++) {
   	        pkd->groupData[i].rcom[j] = PutInBox(pkd->groupData[i].rcom[j],l[j]);
		pkd->groupData[i].r[j] = PutInBox(pkd->groupData[i].r[j],l[j]);
		pkd->groupData[i].v[j] /= pkd->groupData[i].fMass;
	    }
	    pkd->groupData[nMyGroups] = pkd->groupData[i];
	    nMyGroups++;
	}
    }
    if (nMyGroups == pkd->nGroups) {
      pkd->groupData = (FOFGD *) realloc(pkd->groupData,(nMyGroups+1)*sizeof(FOFGD));
      assert(pkd->groupData != NULL);
    }
    pkd->groupData[nMyGroups].bMyGroup = 0;

    /* Below the master node collects all group information, to write it out later.
    **
    ** Ideally this would be replaced by a scheme where every node writes out the data
    ** of his own groups, resp. sends them to his I/O node to do so for him.
    **
    ** But since the amount of data is rather small (e.g. it gave ascii files of around 50 MB for VL-2)
    ** the current serial output of group data is not expected to become a bottleneck anytime soon.
    /*

    /* Start RO group data cache and master reads and saves all the group data. */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD), nMyGroups + 1);
    if (pkd->idSelf == 0) {
	listSize = pkd->nThreads*(pkd->nGroups+1);
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,listSize*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	for (id=1; id < pkd->nThreads; id++) {
	    index = 0;
	    while (1) {
		sG = mdlAquire(mdl,CID_GROUP,index,id);
		mdlRelease(mdl,CID_GROUP,sG);
		if (sG->bMyGroup != 0) {
		    if (nMyGroups >= listSize-1) {
			listSize *= 2.0;
			pkd->groupData = (FOFGD *) realloc(pkd->groupData, listSize*sizeof(FOFGD));
			assert(pkd->groupData != NULL);
			}
		    pkd->groupData[nMyGroups] = *sG;
		    index++;
		    nMyGroups++;
		    }
		else {
		    break;
		    }
		}
	    }
	pkd->groupData[nMyGroups].bMyGroup = 0;
	}
    mdlFinishCache(mdl,CID_GROUP);
    if (pkd->idSelf != 0) {
	nMyGroups = 0;
	}
    else {
	/*
	** Master orders groupData
	*/
	qsort(pkd->groupData,nMyGroups,sizeof(FOFGD), CmpGroups);
	}
    pkd->nGroups = nMyGroups;
    return nMyGroups;
    }

int smGroupProfiles(SMX smx, SMF *smf, int nTotalGroups) {
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p;
    int32_t *pBin, *pPartBin;
    double *v;
    double dx2;
    FLOAT l[3],L[3],r[3],relvel[3],com[3];
    FLOAT rvir,Mvir,fBall,lastbin,fMass,fAvgDens;
    int pn,i,j,k,iBin,nBins,nTree,index,nCnt,pnn;
    int iStart[3],iEnd[3];
    int ix,iy,iz;
    FOFGD *gdp;
    FOFBIN *bin;

    if (nTotalGroups==0) return 0;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */
    if (smx->bPeriodic) {
	for (j=0;j<3;j++) l[j] = pkd->fPeriod[j];
	}
    else {
	for (j=0;j<3;j++) l[j] = FLOAT_MAXVAL;
	}
    nTree = pkd->kdNodes[ROOT].pUpper + 1;
    /*
    ** Start RO group data cache and read all if you are not master
    **
    ** Every node needs to know all group positions and properties to set up the bins. 
    ** Only konwing the properties of ones own groups would not be enough,
    ** because often local particles contribute mass to a bin around
    ** a group which belonged to another node.
    **
    ** This is another reason to keep having the master collect everything, so it can be
    ** distributed to everybody here. 
    */
    if (pkd->idSelf != 0) {
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,nTotalGroups*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	}
    mdlROcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD),pkd->nGroups);
    if (pkd->idSelf != 0) {
	for (i=0; i< nTotalGroups; i++) {
	    gdp = mdlAquire(mdl,CID_GROUP,i,0);
	    mdlRelease(mdl,CID_GROUP,gdp);
	    pkd->groupData[i] = *gdp;
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Allocate memory for the bins.
    */
    nBins = nTotalGroups*smf->nBins;
    if ( pkd->groupBin != NULL ) free(pkd->groupBin);
    pkd->groupBin = (FOFBIN *) malloc( (nBins)*sizeof(FOFBIN) );
    assert(pkd->groupBin != NULL);
    /*
    ** Initalize bin array
    */
    iBin = 0;
    for (i=0; i< nTotalGroups; i++) {
	if (smf->bLogBins == 2) { /* logarithmic bins with same fixed non-comoving size for all groups */
	    lastbin = smf->binFactor/smf->a; /* NB: lastbin is actually the first bin in this case */
	    }
	else {
	  /* estimate virial radius, assuming isothermal shperes */
	  fAvgDens = 0.5*pkd->groupData[i].fMass
	    *0.238732414/pow(pkd->groupData[i].fRMSRadius,3.0); /*using half mass radius and assuming spherical halos*/
	  lastbin = pow(fAvgDens,0.5)*pkd->groupData[i].fRMSRadius*smf->binFactor;
	}
	for (j=0; j < smf->nBins; j++) {
	    pkd->groupBin[iBin].nMembers = 0;
	    if (smf->bLogBins == 1) {/* logarithmic bins */
		dx2 = smf->nBins-(j+1);
		pkd->groupBin[iBin].fRadius = smf->fMinRadius*pow(lastbin/smf->fMinRadius,((float) j)/( (float)smf->nBins-1.0));
		}
	    else if (smf->bLogBins == 2) {/* logarithmic bins with same fixed non-comoving size for all groups */
		pkd->groupBin[iBin].fRadius = lastbin*pow(10.0, 0.1*j);
		}
	    else { /* linear bins */
		dx2 = j+1;
		pkd->groupBin[iBin].fRadius = lastbin*dx2/smf->nBins;
		}
	    pkd->groupBin[iBin].fMassInBin = 0.0;
	    iBin++;
	    }
	}
    /*
    ** Add local particles to their corresponding bins
    */
    for (index=0; index< nTotalGroups; index++) {
      k = index*smf->nBins + smf->nBins-1;
      fBall = pkd->groupBin[k].fRadius;
      
      for (j = 0; j < 3; j++) {
	if(smf->iCenterType == 0)
	  com[j] = pkd->groupData[index].rcom[j];
	else
	  com[j] = pkd->groupData[index].r[j];
      }
      smx->nnListSize = 0;
      if (smx->bPeriodic) {
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
      for (pnn=0;pnn<nCnt;++pnn ) {
	iBin = index*smf->nBins;
	k = 0;
	while (smx->nnList[pnn].fDist2 > pkd->groupBin[iBin+k].fRadius*pkd->groupBin[iBin+k].fRadius) {
	  k++;
	  if ( k == smf->nBins) goto nextParticle;
	}
	fMass = pkdMass(pkd,smx->nnList[pnn].pPart);
	pkd->groupBin[iBin+k].nMembers++;
	pkd->groupBin[iBin+k].fMassInBin += fMass;
      nextParticle:
	;
      }
    }
    /*
    ** Start CO group profiles cache.
    */
    mdlCOcache(mdl,CID_BIN,NULL,pkd->groupBin,sizeof(FOFBIN),nBins,pkd,initGroupBins,combGroupBins);
    if (pkd->idSelf != 0) {
	for (i=0; i< nBins; i++) {
	    if (pkd->groupBin[i].fMassInBin > 0.0) {
		bin = mdlAquire(mdl,CID_BIN,i,0);
		*bin = pkd->groupBin[i];
		mdlRelease(mdl,CID_BIN,bin);
		}
	    }
	}
    mdlFinishCache(mdl,CID_BIN);
    if (pkd->idSelf != 0) {
	free(pkd->groupData);
	free(pkd->groupBin);
	}
    pkd->nBins =  nBins;
    return nBins;
    }
