#include <math.h>
#include "pkd.h"
#include "psd.h"
#include "knn6d.h"
#include "smooth.h"

extern const int primes[];

/*
** Assumes that p does not already occur in the hash table!!!
*/
void psHashAdd(KNN6D  knn,void *p) {
    struct hashElement *t;
    uint32_t i = ((intptr_t)(p))%knn->nHash;
    if (!knn->pHash[i].p) {
	knn->pHash[i].p = p;
    }
    else {
	t = knn->pFreeHash;
	assert(t != NULL);
	knn->pFreeHash = t->coll;
	t->coll = knn->pHash[i].coll;
	knn->pHash[i].coll = t;
	t->p = p;
    }
}

/*
** Assumes that p is definitely in the hash table!!!
*/
void psHashDel(KNN6D  knn,void *p) {
    struct hashElement *t,*tt;
    uint32_t i = ((intptr_t)(p))%knn->nHash;

    if (!knn->pHash[i].coll) {
	/*
	** It has to be the first element.
	*/
	knn->pHash[i].p = NULL;
    }
    else if (knn->pHash[i].p == p) {
	/*
	** It is the first element, but there are others!
	*/
	t = knn->pHash[i].coll;
	knn->pHash[i].coll = t->coll;
	knn->pHash[i].p = t->p;
	t->coll = knn->pFreeHash;
	knn->pFreeHash = t;
    }
    else {
	tt = &knn->pHash[i];
	while (tt->coll->p != p) tt = tt->coll;
	t = tt->coll;
	tt->coll = t->coll; /* unlink */
	t->coll = knn->pFreeHash;
	knn->pFreeHash = t;	
    }
}


int psHashPresent(KNN6D  knn,void *p) {
    struct hashElement *t;
    uint32_t i = ((intptr_t)(p))%knn->nHash;

    if (knn->pHash[i].p == p) return 1;
    t = knn->pHash[i].coll;
    while (t) {
	if (t->p == p) return 1;
	else t = t->coll;
    }
    return 0;
}

int knn6dInitialize(PKD pkd, KNN6D knn, int pqSize,int bPeriodic) {
    int i,j;
    int iTopDepth;

    /*
    ** Allocate Nearest-Neighbor List.
    */
    knn->nnListSize = 0;
    knn->nnListMax = NNLIST_INCREMENT;
    knn->bPeriodic = bPeriodic;

    /*
    ** Allocate priority queue.
    */
    knn->pqSize = pqSize;
    knn->pq = malloc(pqSize*sizeof(*knn->pq));
    assert(knn->pq != NULL);
    PQ6_INIT(knn->pq,pqSize);
    /*
    ** Allocate hash table entries.
    ** The constant here just sets the hash table loading factor, for numbers larger than
    ** the 1000'th prime we end up using the result here as the hash table modulus.
    */
    knn->nHash = (int)floor(pqSize*1.543765241931);
    for (i=0;i<1000;++i) {
	if (primes[i] > knn->nHash) {
	    knn->nHash = primes[i];
	    break;
	}
    }

    knn->pHash = malloc((knn->nHash+pqSize)*sizeof(struct hashElement));
    assert(knn->pHash != NULL);
    for (i=0;i<knn->nHash;++i) {
	knn->pHash[i].p = NULL;
	knn->pHash[i].coll = NULL;
    }
    /*
    ** set up the extra entries that may be needed for collision chains
    */
    knn->pFreeHash = &knn->pHash[i];
    for (;i<(knn->nHash+pqSize-1);++i) {
	knn->pHash[i].p = NULL;
	knn->pHash[i].coll = &knn->pHash[i+1];
    }
    knn->pHash[i].p = NULL;
    knn->pHash[i].coll = NULL;
    /*
    ** Allocate special stacks for searching.
    ** There is a mistake here, since I use these stacks for the remote trees as well.
    ** This can be easily fixed, but a hack for now.
    */
    knn->S = malloc(1024*sizeof(int));
    assert(knn->S != NULL);
    knn->Smin = malloc(1024*sizeof(FLOAT));
    assert(knn->Smin != NULL);
    /*
    ** Allocate special stacks for searching within the top tree.
    ** Calculate the number of levels in the top tree.
    */
    iTopDepth = 1+(int)ceil(log((double)pkd->nThreads)/log(2.0));
    knn->ST = malloc(iTopDepth*sizeof(int));
    assert(knn->ST != NULL);
    knn->SminT = malloc(iTopDepth*sizeof(FLOAT));
    assert(knn->SminT != NULL);
    /*
    ** Set up the sentinel particle with some very far away distance.
    ** This is used to initially load the priority queue and all references
    ** to this particle should end up being replaced in the priority queue
    ** as long as there are pqSize particles set bSrcActive=1.
    */
    for (j=0;j<3;++j) {
	knn->pSentinel.r[j] = HUGE_VAL;
    }
    knn->pSentinel.bSrcActive = 1;
    knn->pSentinel.bDstActive = 0;
    /*
    ** Need to cast the pLite to an array of extra stuff.
    */
    //assert(sizeof(PLITE) >= sizeof(struct smExtraArray));
    //knn->ea = UNION_CAST(pkd->pLite,PLITE *,struct smExtraArray *);
    knn->ea = malloc((pkd->nLocal+1) * sizeof(struct smExtraArray));
    assert(knn->ea != NULL);
    return(1);
}


void knn6dFree(PKD pkd, KNN6D  smx) {
    char achOut[128];

    /*
     * Output statistics.
     */
#if 0
    sprintf(achOut, "Cell Accesses: %g\n",
	    mdlNumAccess(pkd->mdl,CID_CELL));
    mdlDiag(pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(pkd->mdl,CID_CELL));
    mdlDiag(pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(pkd->mdl,CID_CELL));
    mdlDiag(pkd->mdl, achOut);
    sprintf(achOut, "Particle Accesses: %g\n",
	    mdlNumAccess(pkd->mdl,CID_PARTICLE));
    mdlDiag(pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(pkd->mdl,CID_PARTICLE));
    mdlDiag(pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(pkd->mdl,CID_PARTICLE));
    mdlDiag(pkd->mdl, achOut);
#endif

    /*
    ** Free up context storage.
    */
    free(smx->S);
    free(smx->Smin);
    free(smx->ST);
    free(smx->SminT);
    free(smx->pq);
    free(smx->pHash);
    free(smx->ea);
    free(smx);
}

static PQ6 *updateParticleQueue(PKD pkd, KNN6D  smx, KDN *kdn, PQ6 *pq, FLOAT *r, FLOAT *v, FLOAT *rscale, FLOAT *vscale)
{
    int pEnd;
    int pj;
    FLOAT fDist2;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    PARTICLE *p;
    double *pv;
    int idSelf = pkd->idSelf;

    pEnd = kdn->pUpper;
    for (pj=kdn->pLower;pj<=pEnd;++pj) {
/**/    if (smx->ea[pj].bInactive) continue;
/**/    p = pkdParticle(pkd,pj);
	pv = pkdVel(pkd,p);
	dr0 = (p->r[0] - r[0]) * rscale[0];
	dr1 = (p->r[1] - r[1]) * rscale[1];
	dr2 = (p->r[2] - r[2]) * rscale[2];
	dr3 = (pv[0] - v[0]) * vscale[0];
	dr4 = (pv[1] - v[1]) * vscale[1];
	dr5 = (pv[2] - v[2]) * vscale[2];
	fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
	//assert(!isnan(fDist2));
	if (fDist2 < pq->fDist2) {
	    if (pq->iPid == idSelf) {
		smx->ea[pq->iIndex].bInactive = 0;
	    } 
	    else {
		psHashDel(smx,pq->pPart);  //?
		mdlRelease(pkd->mdl,CID_PARTICLE,pq->pPart);
		pq->iPid = idSelf;
	    }
	    pq->pPart = p;
	    pq->fDist2 = fDist2;
	    pq->dr[0] = dr0;
	    pq->dr[1] = dr1;
	    pq->dr[2] = dr2;
	    pq->dr[3] = dr3;
	    pq->dr[4] = dr4;
	    pq->dr[5] = dr5;
	    pq->iIndex = pj;
	    smx->ea[pj].bInactive = 1; /* de-activate a particle that enters the queue */
	    PQ6_REPLACE(pq);
	}
    }

    return pq;
}


static int descendTree(PKD pkd, KDN **kdn0, PQ6 *pq, FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale, int *iCell0, int *S, int *sp0, FLOAT *Smin, int *sm0)
{
    FLOAT min1, min2;
    pBND bnd[2];
    int iCell = *iCell0;
    int sp = *sp0;
    int sm = *sm0;
    int ret = 1;
    KDN *kdn = *kdn0;

    while (kdn->iLower) {
	kdn = pkdTreeNode(pkd,iCell = kdn->iLower); 
	pkdNodeBnd(pkd, kdn, &bnd[0]);
	pkdNodeVBnd(pkd, kdn, &bnd[1]);
	PSMINDIST(bnd,r,v,rscale,vscale, min1);

	kdn = pkdTreeNode(pkd,++iCell);
	pkdNodeBnd(pkd, kdn, &bnd[0]);
	pkdNodeVBnd(pkd, kdn, &bnd[1]);
	PSMINDIST(bnd,r,v,rscale,vscale, min2);

	if (min1 < min2) {
	    Smin[sm++] = min2;
	    kdn = pkdTreeNode(pkd,--iCell);
	    if (min1 >= pq->fDist2) goto NotContained;
	}
	else {
	    Smin[sm++] = min1;
	    if (min2 >= pq->fDist2) goto NotContained;
	}
    }

    ret = 0;

NotContained:

    *sm0 = sm;
    *sp0 = sp;
    *iCell0 = iCell;
    *kdn0 = kdn;

    return ret;
}

/*
** This function performs a local nearest neighbor search.
*/
static PQ6 *_SearchLocal(PKD pkd, KNN6D smx, PQ6 *pq, FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale, int *pbDone) {
    KDN *kdn;
    PARTICLE *p;
    FLOAT dMin,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int j,pj,iCell,iParent;
    int sp = 0;
    int sm = 0;
    int idSelf = pkd->idSelf;
    pBND bnd[2];

    *pbDone = 1;	/* assume that we will complete the search */
    /*
    ** We don't perform containment tests except at the
    ** root, so that the pbDone flag can be correctly
    ** set.
    */
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    S[sp] = iCell;
    /*
    ** Start of PRIOQ searching loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	if (descendTree(pkd, &kdn, pq, r,v, rscale, vscale, &iCell, S, &sp, Smin, &sm))
	    goto NotContained;
	
	pq = updateParticleQueue(pkd, smx, kdn, pq, r,v, rscale, vscale);

    NoIntersect:
	while (iCell == S[sp]) {
	    if (sp) {
		--sp;
		kdn = pkdTreeNode(pkd,iCell = kdn->iParent);
	    }
	    else {
		/*
		** Containment Test!
		*/
		pkdNodeBnd(pkd, kdn, &bnd[0]);
		for (j=0;j<3;++j) {
		    dMin = (bnd[0].fMax[j] - fabs(bnd[0].fCenter[j] - r[j])) * rscale[j];
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = kdn->iParent;
			if (!iParent) {
			    *pbDone = 0;		/* EXIT, not contained! */
			    break;
			}
			S[sp] = iParent;
			goto NotContained;
		    }
		}

		pkdNodeVBnd(pkd, kdn, &bnd[1]);
		for (j=0;j<3;++j) {
		    dMin = (bnd[1].fMax[j] - fabs(bnd[1].fCenter[j] - v[j])) * vscale[j];
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = kdn->iParent;
			if (!iParent) {
			    *pbDone = 0;		/* EXIT, not contained! */
			    break;
			}
			S[sp] = iParent;
			goto NotContained;
		    }
		}

		//fprintf(stderr, "cnt is %i\n", cnt);
		return pq;
	    }
	}
    NotContained:
	kdn = pkdTreeNode(pkd,iCell ^= 1);
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
	    pkdNodeBnd(pkd, kdn, &bnd[0]);
	    pkdNodeVBnd(pkd, kdn, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min2);
	}
	if (min2 >= pq->fDist2) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iParent);
	    goto NoIntersect;
	}
	S[++sp] = iCell;
    }

}



static PQ6 *_SearchRemote(PKD pkd, KNN6D smx, PQ6 *pq, int id,FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale) {
    KDN *kdn;
    PARTICLE *p;
    double *pv;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    FLOAT min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int pj,pEnd,iCell;
    int sp = 0;
    int sm = 0;
    int idSelf = pkd->idSelf;
    pBND bnd[2];
    MDL mdl = pkd->mdl;
    KDN *pkdn,*pkdu;

    assert(id != idSelf);
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    S[sp] = iCell;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    /*
    ** Start of PRIOQ searching loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (pkdn->iLower) {
	    kdn  = pkdTreeNode(pkd,iCell = pkdn->iLower);
/**/	    mdlRelease(mdl,CID_CELL,pkdn);
/**/	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    pkdNodeBnd(pkd,pkdn, &bnd[0]);
	    pkdNodeVBnd(pkd,pkdn, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min1);

	    kdn  = pkdTreeNode(pkd,++iCell);
/**/	    pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
	    pkdNodeBnd(pkd,pkdu, &bnd[0]);
	    pkdNodeVBnd(pkd,pkdu, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		kdn = pkdTreeNode(pkd,--iCell);
/**/		mdlRelease(mdl,CID_CELL,pkdu);
		if (min1 >= pq->fDist2) goto NotContained;
	    }
	    else {
		Smin[sm++] = min1;
/**/		mdlRelease(mdl,CID_CELL,pkdn);
/**/		pkdn = pkdu;
		if (min2 >= pq->fDist2) goto NotContained;
	    }
	}
	pEnd = pkdn->pUpper;
	for (pj=pkdn->pLower;pj<=pEnd;++pj) {
/**/	    p = mdlAquire(mdl,CID_PARTICLE,pj,id);
	    if (!p->bSrcActive || psHashPresent(smx,p)) {
/**/		mdlRelease(mdl,CID_PARTICLE,p);
/**/		continue;
/**/	    }
	    pv = pkdVel(pkd,p);
	    dr0 = -(r[0] - p->r[0]) * rscale[0];
	    dr1 = -(r[1] - p->r[1]) * rscale[1];
	    dr2 = -(r[2] - p->r[2]) * rscale[2];
	    dr3 = -(v[0] -   pv[0]) * vscale[0];
	    dr4 = -(v[1] -   pv[1]) * vscale[1];
	    dr5 = -(v[2] -   pv[2]) * vscale[2];
	    fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
	    if (fDist2 < pq->fDist2) {
		if (pq->iPid == idSelf) {
		    smx->ea[pq->iIndex].bInactive = 0;
		}
		else {
		    psHashDel(smx,pq->pPart);
		    mdlRelease(mdl,CID_PARTICLE,pq->pPart);
		}
		pq->pPart = p;
		pq->fDist2 = fDist2;
		pq->dr[0] = dr0;
		pq->dr[1] = dr1;
		pq->dr[2] = dr2;
		pq->dr[3] = dr3;
		pq->dr[4] = dr4;
		pq->dr[5] = dr5;
		pq->iIndex = pj;
		pq->iPid = id;
		psHashAdd(smx,p);
		PQ6_REPLACE(pq);
	    }
/**/	    else mdlRelease(mdl,CID_PARTICLE,p);
	}
    NoIntersect:
	while (iCell == S[sp]) {
	    if (!sp) {
		mdlRelease(mdl,CID_CELL,pkdn);
		return pq;
	    }
	    --sp;
	    kdn = pkdTreeNode(pkd,iCell = pkdn->iParent);
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	}
    NotContained:
	kdn = pkdTreeNode(pkd,iCell ^= 1);
/**/    mdlRelease(mdl,CID_CELL,pkdn);
/**/    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) {
	    min2 = Smin[--sm];
	    }
	else {
	    pkdNodeBnd(pkd, pkdn, &bnd[0]);
	    pkdNodeVBnd(pkd, pkdn, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min2);
	    }
	if (min2 >= pq->fDist2) {
	    kdn = pkdTreeNode(pkd,iCell = pkdn->iParent);
/**/	    mdlRelease(mdl,CID_CELL,pkdn);
/**/	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    goto NoIntersect;
	    }
	S[++sp] = iCell;
    }
}


PQ6 *_Search(PKD pkd, KNN6D knn, PQ6 *pq, FLOAT *r,FLOAT *v, FLOAT *rscale,FLOAT *vscale, int bReplica,int *pbDone) {
    KDN *kdn;
    int idSelf = pkd->idSelf;
    FLOAT *Smin = knn->SminT;
    int *S = knn->ST;
    FLOAT dMin,min1,min2;
    int j,iCell,id,iParent;
    int sp = 0;
    int sm = 0;
    pBND bnd[2];

    *pbDone = 0;
    if (bReplica) kdn = pkdTopNode(pkd,iCell = ROOT);
    else {
	kdn = pkdTopNode(pkd,iCell = pkd->iTopRoot);
	assert(kdn->pLower == idSelf);
    }
    if (iCell != ROOT) S[sp] = kdn->iParent;
    else S[sp] = iCell;


    while (1) {
	mdlCacheCheck(pkd->mdl);

	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (kdn->iLower) {
	    kdn = pkdTopNode(pkd,iCell = kdn->iLower); 
	    pkdNodeBnd(pkd, kdn, &bnd[0]);
	    pkdNodeVBnd(pkd, kdn, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min1);

	    kdn = pkdTopNode(pkd,++iCell);	     
	    pkdNodeBnd(pkd, kdn, &bnd[0]);
	    pkdNodeVBnd(pkd, kdn, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		kdn = pkdTopNode(pkd,--iCell);
		if (min1 >= pq->fDist2) goto NotContained;
	    }
	    else {
		Smin[sm++] = min1;
		if (min2 >= pq->fDist2) goto NotContained;
	    }
	}
	id = kdn->pLower;	/* this is the thread id in LTT */
	if (id == pkd->idSelf) {
	    pq = _SearchLocal(pkd,knn,pq, r,v,rscale,vscale,pbDone);
	    if (*pbDone) return(pq);
	}
	else {
	    pq = _SearchRemote(pkd,knn,pq, id,r,v,rscale,vscale);
	}
    NoIntersect:
	while (iCell == S[sp]) {
	    if (sp) {
		--sp;
		kdn = pkdTopNode(pkd,iCell = kdn->iParent);
	    }
	    else if (!bReplica) {
		/*
		** Containment Test!
		*/
		pkdNodeBnd(pkd, kdn, &bnd[0]);
		for (j=0;j<3;++j) {
		    dMin = (bnd[0].fMax[j] - fabs(bnd[0].fCenter[j] - r[j])) * rscale[j];
		    //dMin = (bnd[0].fMax[j] - fabs(bnd[0].fCenter[j] - r[j]));
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = kdn->iParent;
			if (!iParent) {
			    *pbDone = 0;
			    return pq;
			}
			S[sp] = iParent;
			goto NotContained;
		    }
		}

#if 1
		pkdNodeVBnd(pkd, kdn, &bnd[1]);
		for (j=0;j<3;++j) {
		    dMin = (bnd[1].fMax[j] - fabs(bnd[1].fCenter[j] - v[j])) * vscale[j];
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = kdn->iParent;
			if (!iParent) {
			    *pbDone = 0;
			    return pq;
			}
			S[sp] = iParent;
			goto NotContained;
		    }
		}
#endif

		*pbDone = 1;
		return pq;
	    }
	    else return pq;
	}
    NotContained:
	kdn = pkdTopNode(pkd,iCell ^= 1);
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
	    pkdNodeBnd(pkd, kdn, &bnd[0]);
	    pkdNodeVBnd(pkd, kdn, &bnd[1]);
	    PSMINDIST(bnd,r,v,rscale,vscale, min2);
	}
	if (min2 >= pq->fDist2) {
	    kdn = pkdTopNode(pkd,iCell = kdn->iParent);
	    goto NoIntersect;
	}
	S[++sp] = iCell;
    }
}

/*
** Find the k nearest neighbors. If first_time is true the priority-queue
** and knn_data is initialized. Subsequent calls with first_time false
** will update the priority-queue.
*/
void knn6d(PKD pkd, KNN6D knn, int pid, float *fBall, int first_time)
{
    int i,j;
    FLOAT *r,*v;
    PARTICLE *p;
    int iStart[3],iEnd[3];
    FLOAT R[3];
    int ix,iy,iz;

    if (first_time)
    {
	/*
	** Initialize the bInactive flags for all local particles.
	*/
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    knn->ea[i].bInactive = (p->bSrcActive)?0:1;
	}
	knn->ea[pkd->nLocal].bInactive = 0;  /* initialize for Sentinel, but this is not really needed */

	/*
	** Initialize the priority queue first.
	*/
	for (i=0;i<knn->pqSize;++i) {
	    knn->pq[i].pPart = &knn->pSentinel;
	    knn->pq[i].iIndex = pkd->nLocal;
	    knn->pq[i].iPid = pkd->idSelf;
	    knn->pq[i].dr[0] = HUGE_VAL;
	    knn->pq[i].dr[1] = HUGE_VAL;
	    knn->pq[i].dr[2] = HUGE_VAL;
	    knn->pq[i].dr[3] = HUGE_VAL;
	    knn->pq[i].dr[4] = HUGE_VAL;
	    knn->pq[i].dr[5] = HUGE_VAL;
	    knn->pq[i].fDist2 = HUGE_VAL;
	}
	for (j=0;j<3;++j) knn->data.rLast[j] = 0.0;
	for (j=0;j<3;++j) knn->data.vLast[j] = 0.0;
    }
    else
    {
	for (i=0;i<knn->pqSize;++i)
	{
	    for (j=0; j < 3; j++) if (knn->data.rscale[j] != 0) knn->pq[i].dr[j]   /= knn->data.rscale[j];
	    for (j=0; j < 3; j++) if (knn->data.vscale[j] != 0) knn->pq[i].dr[3+j] /= knn->data.vscale[j];
	}
    }

    p = pkdParticle(pkd, pid);

    r = p->r;
    v = pkdVel(pkd, p);

    for (j=0;j<3;++j) knn->data.rscale[j] = knn->psm[pid].rscale[j];
    for (j=0;j<3;++j) knn->data.vscale[j] = knn->psm[pid].vscale[j];

    FLOAT *rscale = knn->data.rscale;
    FLOAT *vscale = knn->data.vscale;

    if (!first_time)
    {
	for (i=0;i<knn->pqSize;++i) 
	{
	    knn->pq[i].fDist2 = 0;
	    for (j=0; j < 3; j++) 
	    {
		knn->pq[i].dr[j]  -= r[j] - knn->data.rLast[j];
		knn->pq[i].dr[j]  *= knn->data.rscale[j];
		knn->pq[i].fDist2 += pow(knn->pq[i].dr[j],2);
		knn->pq[i].dr[3+j]  -= v[j] - knn->data.vLast[j];
		knn->pq[i].dr[3+j]  *= knn->data.vscale[j];
		knn->pq[i].fDist2 += pow(knn->pq[i].dr[3+j],2);
	    }
	}
    }

    for (j=0;j<3;++j) knn->data.rLast[j] = r[j];
    for (j=0;j<3;++j) knn->data.vLast[j] = v[j];
    PQ6_BUILD(knn->pq,knn->pqSize,knn->pqTEMP);
    int bDone;
    knn->pqTEMP = _Search(pkd,knn,knn->pqTEMP, r,v,rscale,vscale,0,&bDone);

    /*
    ** Search in replica boxes if it is required.
    */
    if (!bDone && knn->bPeriodic) {
	double fBall = sqrt(knn->pqTEMP->fDist2);
	for (j=0;j<3;++j) {
	    iStart[j] = -1; //d2i(floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5));
	    iEnd[j] = +1; //d2i(floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5));
	    //fprintf(stderr, "%i %i %i  %i %i %i\n", iStart[0], iStart[1], iStart[2],
	//					    iEnd[0], iEnd[1], iEnd[2]);
	    }
	for (ix=iStart[0];ix<=iEnd[0];++ix) {
	    R[0] = r[0] - ix*pkd->fPeriod[0]; //) * rscale[0];
	    for (iy=iStart[1];iy<=iEnd[1];++iy) {
		R[1] = r[1] - iy*pkd->fPeriod[1]; //) * rscale[1];
		for (iz=iStart[2];iz<=iEnd[2];++iz) {
		    R[2] = r[2] - iz*pkd->fPeriod[2]; //) * rscale[2];
		    if (ix || iy || iz) {
			knn->pqTEMP = _Search(pkd,knn,knn->pqTEMP, R,v,rscale,vscale,1,&bDone);
			}
		    }
		}
	    }
	}

    if (fBall != NULL)
	*fBall = sqrt(knn->pqTEMP->fDist2);
}

void knn6dFinish(PKD pkd, KNN6D  knn)
{
    int i;
    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<knn->pqSize;++i) {
	if (knn->pq[i].iPid == pkd->idSelf) {
	    knn->ea[knn->pq[i].iIndex].bInactive = 0;
	}
	else {
	    psHashDel(knn,knn->pq[i].pPart);
	    mdlRelease(pkd->mdl,CID_PARTICLE,knn->pq[i].pPart);
	}
    }
}

static PQ6 *_GatherLocal(PKD pkd, KNN6D smx, PQ6 *pq, float fBall2, FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale) {
    KDN *kdn;
    PARTICLE *p;
    FLOAT min2,fDist2;
    double *pv;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,nCnt,pEnd;
    int idSelf = pkd->idSelf;
    pBND bnd[2];

    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
        pkdNodeBnd(pkd, kdn, &bnd[0]);
        pkdNodeVBnd(pkd, kdn, &bnd[1]);
	PSMINDIST(bnd,r,v,rscale,vscale,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	}
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iLower);
	    S[sp++] = iCell+1;
	    continue;
	}
	else {
	    pEnd = kdn->pUpper;
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
	/**/    if (smx->ea[pj].bInactive) continue;
	/**/    p = pkdParticle(pkd,pj);
		pv = pkdVel(pkd,p);
		dr0 = (p->r[0] - r[0]) * rscale[0];
		dr1 = (p->r[1] - r[1]) * rscale[1];
		dr2 = (p->r[2] - r[2]) * rscale[2];
		dr3 = (pv[0] - v[0]) * vscale[0];
		dr4 = (pv[1] - v[1]) * vscale[1];
		dr5 = (pv[2] - v[2]) * vscale[2];
		fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
		//assert(!isnan(fDist2));
		if (fDist2 <= fBall2 && fDist2 < pq->fDist2) {
		    if (pq->iPid == idSelf) {
			smx->ea[pq->iIndex].bInactive = 0;
		    } 
		    else {
			psHashDel(smx,pq->pPart);  //?
			mdlRelease(pkd->mdl,CID_PARTICLE,pq->pPart);
			pq->iPid = idSelf;
		    }
		    pq->pPart = p;
		    pq->fDist2 = fDist2;
		    pq->dr[0] = dr0;
		    pq->dr[1] = dr1;
		    pq->dr[2] = dr2;
		    pq->dr[3] = dr3;
		    pq->dr[4] = dr4;
		    pq->dr[5] = dr5;
		    pq->iIndex = pj;
		    smx->ea[pj].bInactive = 1; /* de-activate a particle that enters the queue */
		    PQ6_REPLACE(pq);
		}
	    }
	}
    NoIntersect:
	if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp]);
	else break;
    }
    return pq;
}


static PQ6 *_GatherRemote(PKD pkd, KNN6D smx, PQ6 *pq, int id, float fBall2, FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale) {
    MDL mdl = pkd->mdl;
    KDN *pkdn;
    PARTICLE *p;
    FLOAT min2,fDist2;
    double *pv;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    int *S = smx->S;
    int sp = 0;
    int pj,pEnd;
    int iCell;
    pBND bnd[2];
    int idSelf = pkd->idSelf;

    assert(id != pkd->idSelf);
    iCell = ROOT;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    while (1) {
        pkdNodeBnd(pkd, pkdn, &bnd[0]);
        pkdNodeVBnd(pkd, pkdn, &bnd[1]);
	PSMINDIST(bnd,r,v,rscale,vscale,min2);
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
		p = mdlAquire(mdl,CID_PARTICLE,pj,id);
	        if (!p->bSrcActive || psHashPresent(smx,p)) {
		    mdlRelease(mdl,CID_PARTICLE,p);
		    continue;
		}
		pv = pkdVel(pkd,p);
		dr0 = (p->r[0] - r[0]) * rscale[0];
		dr1 = (p->r[1] - r[1]) * rscale[1];
		dr2 = (p->r[2] - r[2]) * rscale[2];
		dr3 = (  pv[0] - v[0]) * vscale[0];
		dr4 = (  pv[1] - v[1]) * vscale[1];
		dr5 = (  pv[2] - v[2]) * vscale[2];
		fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
		if (fDist2 <= fBall2 && fDist2 < pq->fDist2) {
		    if (pq->iPid == idSelf) {
			smx->ea[pq->iIndex].bInactive = 0;
		    }
		    else {
			psHashDel(smx,pq->pPart);
			mdlRelease(mdl,CID_PARTICLE,pq->pPart);
		    }
		    pq->pPart = p;
		    pq->fDist2 = fDist2;
		    pq->dr[0] = dr0;
		    pq->dr[1] = dr1;
		    pq->dr[2] = dr2;
		    pq->dr[3] = dr3;
		    pq->dr[4] = dr4;
		    pq->dr[5] = dr5;
		    pq->iIndex = pj;
		    pq->iPid = id;
		    psHashAdd(smx,p);
		    PQ6_REPLACE(pq);
		}
		else mdlRelease(mdl,CID_PARTICLE,p);
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
    return pq;
}


static PQ6 *_Gather(PKD pkd, KNN6D  knn, PQ6 *pq, float fBall2, FLOAT *r,FLOAT *v, FLOAT *rscale,FLOAT *vscale) {
    KDN *kdn;
    int *S = knn->ST;
    FLOAT min2;
    int iCell,id;
    int sp = 0;
    pBND bnd[2];

    kdn = pkdTopNode(pkd,iCell = ROOT);
    while (1) {
        pkdNodeBnd(pkd, kdn, &bnd[0]);
        pkdNodeVBnd(pkd, kdn, &bnd[1]);
	PSMINDIST(bnd,r,v,rscale,vscale,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	}
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTopNode(pkd,iCell = kdn->iLower);
	    S[sp++] = iCell+1;
	    continue;
	}
	else {
	    id = kdn->pLower; /* this is the thread id in LTT */
	    if (id != pkd->idSelf) {
		pq = _GatherRemote(pkd,knn,pq, id,fBall2, r,v,rscale,vscale);
	    }
	    else {
		pq = _GatherLocal(pkd,knn,pq, fBall2, r,v,rscale,vscale);
	    }
	}
    NoIntersect:
	if (sp) kdn = pkdTopNode(pkd,iCell = S[--sp]);
	else break;
    }

    return pq;
}

void knn6dGather(PKD pkd, KNN6D knn, float fBall, int pid, int first_time) {

    float fBall2 = fBall * fBall;

    int i,j;
    FLOAT *r,*v;
    PARTICLE *p;
    int iStart[3],iEnd[3];
    FLOAT R[3];
    int ix,iy,iz;

    if (first_time)
    {
	/*
	** Initialize the bInactive flags for all local particles.
	*/
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    knn->ea[i].bInactive = (p->bSrcActive)?0:1;
	}
	knn->ea[pkd->nLocal].bInactive = 0;  /* initialize for Sentinel, but this is not really needed */

	/*
	** Initialize the priority queue first.
	*/
	for (i=0;i<knn->pqSize;++i) {
	    knn->pq[i].pPart = &knn->pSentinel;
	    knn->pq[i].iIndex = pkd->nLocal;
	    knn->pq[i].iPid = pkd->idSelf;
	    knn->pq[i].dr[0] = HUGE_VAL;
	    knn->pq[i].dr[1] = HUGE_VAL;
	    knn->pq[i].dr[2] = HUGE_VAL;
	    knn->pq[i].dr[3] = HUGE_VAL;
	    knn->pq[i].dr[4] = HUGE_VAL;
	    knn->pq[i].dr[5] = HUGE_VAL;
	    knn->pq[i].fDist2 = HUGE_VAL;
	}
	for (j=0;j<3;++j) knn->data.rLast[j] = 0.0;
	for (j=0;j<3;++j) knn->data.vLast[j] = 0.0;
    }
    else
    {
	for (i=0;i<knn->pqSize;++i)
	{
	    for (j=0; j < 3; j++) if (knn->data.rscale[j] != 0) knn->pq[i].dr[j]   /= knn->data.rscale[j];
	    for (j=0; j < 3; j++) if (knn->data.vscale[j] != 0) knn->pq[i].dr[3+j] /= knn->data.vscale[j];
	}
    }

    p = pkdParticle(pkd, pid);

    r = p->r;
    v = pkdVel(pkd, p);

    for (j=0;j<3;++j) knn->data.rscale[j] = knn->psm[pid].rscale[j];
    for (j=0;j<3;++j) knn->data.vscale[j] = knn->psm[pid].vscale[j];

    FLOAT *rscale = knn->data.rscale;
    FLOAT *vscale = knn->data.vscale;

    if (!first_time)
    {
	for (i=0;i<knn->pqSize;++i) 
	{
	    knn->pq[i].fDist2 = 0;
	    for (j=0; j < 3; j++) 
	    {
		knn->pq[i].dr[j]  -= r[j] - knn->data.rLast[j];
		knn->pq[i].dr[j]  *= knn->data.rscale[j];
		knn->pq[i].fDist2 += pow(knn->pq[i].dr[j],2);
		knn->pq[i].dr[3+j]  -= v[j] - knn->data.vLast[j];
		knn->pq[i].dr[3+j]  *= knn->data.vscale[j];
		knn->pq[i].fDist2 += pow(knn->pq[i].dr[3+j],2);
	    }
	}
    }

    for (j=0;j<3;++j) knn->data.rLast[j] = r[j];
    for (j=0;j<3;++j) knn->data.vLast[j] = v[j];
    PQ6_BUILD(knn->pq,knn->pqSize,knn->pqTEMP);

    if (knn->bPeriodic)
    {
	for (j=0;j<3;++j) {
	    iStart[j] = -1; //d2i(floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5));
	    iEnd[j] = +1; //d2i(floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5));
	    //fprintf(stderr, "%i %i %i  %i %i %i\n", iStart[0], iStart[1], iStart[2],
	//					    iEnd[0], iEnd[1], iEnd[2]);
	    }
	for (ix=iStart[0];ix<=iEnd[0];++ix) {
	    R[0] = r[0] - ix*pkd->fPeriod[0];
	    for (iy=iStart[1];iy<=iEnd[1];++iy) {
		R[1] = r[1] - iy*pkd->fPeriod[1];
		for (iz=iStart[2];iz<=iEnd[2];++iz) {
		    R[2] = r[2] - iz*pkd->fPeriod[2];
		    knn->pqTEMP = _Gather(pkd,knn,knn->pqTEMP, fBall2,R,v,rscale,vscale);
		    }
		}
	    }
	}
    else
	knn->pqTEMP = _Gather(pkd,knn,knn->pqTEMP, fBall2,r,v,rscale,vscale);

}

