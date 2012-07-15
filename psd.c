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
#include "psd.h"
#include "pkd.h"
#include "rbtree.h"
#include <sys/stat.h>

#define PROP(pkd,p) (-*pkdPot((pkd), (p)))
//#define PROP(pkd,p) (p->fDensity)

#define PSM(i) (psx->psm[i])
//#define PSM(i) (pkdParticle(pkd, i)->psm)

extern const int primes[];

/*
** The Epanechnikov kernel
*/
static double ekernel(double u2) {
    return (0 <= u2 && u2 <= 1) * (1-u2);
}

static double grad_ekernel(double u) {
    return (fabs(u) <= 1) * (-2*u);
}

/* 
** A spline kernel (unused)
*/
static inline float spline_kernel(float u) {
    if (0 < u && u <= 0.5)
	return 1 - (6*u*u) * (1 + u);
    else if (0.5 < u && u <= 1)
	return 2 * pow(1-u, 3);
    return 0;
}

static inline float grad_spline_kernel(float x, float u) {
    if (0 < u && u <= 0.5)
	return (6*x) * (2 - 3*u);
    else if (0.5 < u && u <= 1)
	return -3*x/u * (1 - 2*u + u*u);
    return 0;
}

/*
** The smoothed phase-space density.
*/
float psdDensity(PKD pkd, PARTICLE *p,int nSmooth,NN6 *nnList, FLOAT *rscale, FLOAT *vscale) {
    double ih2,r2,fDensity,fMass;
    int i;

    ih2 = BALL2(p);
    assert(ih2 > 0);
    fDensity = 0.0;

    fDensity = pkdMass(pkd, p) * ekernel(pow(1 * 6 / (1+6.), 2));
    for (i=0;i<nSmooth;++i) {
	if (nnList[i].pPart == p) continue;
	fMass = pkdMass(pkd,nnList[i].pPart);
	assert(fMass > 0);
	r2 = nnList[i].fDist2/ih2;
	fDensity += ekernel(r2) * fMass;
	}

    float V0=1, V1=1;
    for (i=0; i < 3; i++) {
	if (rscale[i] > 0)
	{
	    V0 *= p->fBall;
	    V1 *= rscale[i];
	}
    }

    for (i=0; i < 3; i++) {
	if (vscale[i] > 0)
	{
	    V0 *= p->fBall;
	    V1 *= vscale[i];
	}
    }

    assert(!isnan(fDensity));
    assert(!isinf(fDensity));

    return 0.77403670 * fDensity * V1 / V0;

}

/*
** Compute the normalized density gradient. Also correct the E0 error
** to reduce the noise.
*/
void psdDensityGrad(PKD pkd, PSX psx, int pid, FLOAT *fDensityGrad) {
    int i,j;
    FLOAT *rscale, *vscale;
    PARTICLE *p = pkdParticle(pkd, pid);
    const FLOAT fBall = p->fBall;
    const FLOAT fDensity = p->fDensity;

    rscale = PSM(pid).rscale;
    vscale = PSM(pid).vscale;

    for (j=0; j < 6; j++) fDensityGrad[j] = 0;

    for (i=0;i < psx->nSmooth;++i)
    {
	double r = sqrt(psx->pq[i].fDist2) / fBall;
	if (r > 0)
	{
	    for (j=0; j < 6; j++)
	    {
		double dx = psx->pq[i].dr[j]/fBall;
		double c = (fDensity - psx->pq[i].pPart->fDensity) / fDensity;
		fDensityGrad[j] += (pkdMass(pkd, psx->pq[i].pPart) * c) * grad_ekernel(dx);
	    }
	}
    }

    double L = 0;
    for (i=0; i < 6; i++) L += pow(fDensityGrad[i], 2);
    L = sqrt(L);
    for (i=0; i < 6; i++) fDensityGrad[i] /= L;
}

/*
** Compute the normalized potential gradient. Also correct the E0 error
** to reduce the noise.
*/
void psdPotentialGrad(PKD pkd, PSX psx, int pid, FLOAT *fPotentialGrad) {
    int i,j;
    FLOAT *rscale, *vscale;
    PARTICLE *p = pkdParticle(pkd, pid);
    const FLOAT fBall = p->fBall;
    const FLOAT fPot = *pkdPot(pkd, p);
    const FLOAT fDensity = p->fDensity;

    rscale = PSM(pid).rscale;
    vscale = PSM(pid).vscale;

    for (j=0; j < 6; j++) fPotentialGrad[j] = 0;

    for (i=0;i < psx->nSmooth;++i)
    {
	double r = sqrt(psx->pq[i].fDist2) / fBall;
	if (r > 0)
	{
	    for (j=0; j < 6; j++)
	    {
		double dx = psx->pq[i].dr[j]/fBall;
		double c = (fPot - *pkdPot(pkd, psx->pq[i].pPart)) / fDensity;
		fPotentialGrad[j] += (pkdMass(pkd, psx->pq[i].pPart) * c) * grad_ekernel(dx);
	    }
	}
    }

    double L = 0;
    for (i=0; i < 6; i++) L += pow(fPotentialGrad[i], 2);
    L = sqrt(L);
    for (i=0; i < 6; i++) fPotentialGrad[i] /= L;
}

/*
** Assumes that p does not already occur in the hash table!!!
*/
void psHashAdd(PSX psx,void *p) {
    struct hashElement *t;
    uint32_t i = ((intptr_t)(p))%psx->nHash;
    if (!psx->pHash[i].p) {
	psx->pHash[i].p = p;
    }
    else {
	t = psx->pFreeHash;
	assert(t != NULL);
	psx->pFreeHash = t->coll;
	t->coll = psx->pHash[i].coll;
	psx->pHash[i].coll = t;
	t->p = p;
    }
}

/*
** Assumes that p is definitely in the hash table!!!
*/
void psHashDel(PSX psx,void *p) {
    struct hashElement *t,*tt;
    uint32_t i = ((intptr_t)(p))%psx->nHash;

    if (!psx->pHash[i].coll) {
	/*
	** It has to be the first element.
	*/
	psx->pHash[i].p = NULL;
    }
    else if (psx->pHash[i].p == p) {
	/*
	** It is the first element, but there are others!
	*/
	t = psx->pHash[i].coll;
	psx->pHash[i].coll = t->coll;
	psx->pHash[i].p = t->p;
	t->coll = psx->pFreeHash;
	psx->pFreeHash = t;
    }
    else {
	tt = &psx->pHash[i];
	while (tt->coll->p != p) tt = tt->coll;
	t = tt->coll;
	tt->coll = t->coll; /* unlink */
	t->coll = psx->pFreeHash;
	psx->pFreeHash = t;	
    }
}


int psHashPresent(PSX psx,void *p) {
    struct hashElement *t;
    uint32_t i = ((intptr_t)(p))%psx->nHash;

    if (psx->pHash[i].p == p) return 1;
    t = psx->pHash[i].coll;
    while (t) {
	if (t->p == p) return 1;
	else t = t->coll;
    }
    return 0;
}

int psdInitialize(PSX smx,PKD pkd, PSF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType, int initCache) {
    int i,j;
    int iTopDepth;

    if (smf != NULL) smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->bPeriodic = bPeriodic;

    /*
    ** Allocate Nearest-Neighbor List.
    */
    smx->nnListSize = 0;
    smx->nnListMax = NNLIST_INCREMENT;
    smx->nnList = malloc(smx->nnListMax*sizeof(NN6));
    assert(smx->nnList != NULL);

    /*
    ** Allocate priority queue.
    */
    smx->pq = malloc(nSmooth*sizeof(PQ6));
    assert(smx->pq != NULL);
    PQ6_INIT(smx->pq,nSmooth);
    /*
    ** Allocate hash table entries.
    ** The constant here just sets the hash table loading factor, for numbers larger than
    ** the 1000'th prime we end up using the result here as the hash table modulus.
    */
    smx->nHash = (int)floor(nSmooth*1.543765241931);
    for (i=0;i<1000;++i) {
	if (primes[i] > smx->nHash) {
	    smx->nHash = primes[i];
	    break;
	}
    }

    smx->pHash = malloc((smx->nHash+nSmooth)*sizeof(struct hashElement));
    assert(smx->pHash != NULL);
    for (i=0;i<smx->nHash;++i) {
	smx->pHash[i].p = NULL;
	smx->pHash[i].coll = NULL;
    }
    /*
    ** set up the extra entries that may be needed for collision chains
    */
    smx->pFreeHash = &smx->pHash[i];
    for (;i<(smx->nHash+nSmooth-1);++i) {
	smx->pHash[i].p = NULL;
	smx->pHash[i].coll = &smx->pHash[i+1];
    }
    smx->pHash[i].p = NULL;
    smx->pHash[i].coll = NULL;
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
    /*
    ** Set up the sentinel particle with some very far away distance.
    ** This is used to initially load the priority queue and all references
    ** to this particle should end up being replaced in the priority queue
    ** as long as there are nSmooth particles set bSrcActive=1.
    */
    for (j=0;j<3;++j) {
	smx->pSentinel.r[j] = HUGE_VAL;
    }
    smx->pSentinel.bSrcActive = 1;
    smx->pSentinel.bDstActive = 0;
    /*
    ** Need to cast the pLite to an array of extra stuff.
    */
    //assert(sizeof(PLITE) >= sizeof(struct smExtraArray));
    //smx->ea = UNION_CAST(pkd->pLite,PLITE *,struct smExtraArray *);
    smx->ea = malloc((pkd->nLocal+1) * sizeof(struct smExtraArray));
    assert(smx->ea != NULL);
    return(1);
}


void psdFinish(PSX smx, PSF *smf) {
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
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "Particle Accesses: %g\n",
	    mdlNumAccess(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);

    /*
    ** Free up context storage.
    */
    free(smx->S);
    free(smx->Smin);
    free(smx->ST);
    free(smx->SminT);
    free(smx->pq);
    free(smx->nnList);
    free(smx->pHash);
    free(smx);

    free(smx->psm); smx->psm = NULL;
}


/*
** This function performs a local nearest neighbor search.
*/
PQ6 *pqSearchLocalPsd(PSX smx,PQ6 *pq,FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale, int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    double *pv;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    FLOAT dMin,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int j,pj,pEnd,iCell,iParent;
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
	pEnd = kdn->pUpper;
	for (pj=kdn->pLower;pj<=pEnd;++pj) {
/**/	if (smx->ea[pj].bInactive) continue;
/**/	p = pkdParticle(pkd,pj);
	    pv = pkdVel(pkd,p);
	    dr0 = -(r[0] - p->r[0]) * rscale[0];
	    dr1 = -(r[1] - p->r[1]) * rscale[1];
	    dr2 = -(r[2] - p->r[2]) * rscale[2];
	    dr3 = -(v[0] -   pv[0]) * vscale[0];
	    dr4 = -(v[1] -   pv[1]) * vscale[1];
	    dr5 = -(v[2] -   pv[2]) * vscale[2];
	    fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
	    assert(!isnan(fDist2));
	    //if (p->fDensity < 1e6) continue;
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



PQ6 *pqSearchRemotePsd(PSX smx,PQ6 *pq,int id,FLOAT *r,FLOAT *v, FLOAT *rscale, FLOAT *vscale) {
    PKD pkd = smx->pkd;
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
    MDL mdl = smx->pkd->mdl;
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
/**/	else mdlRelease(mdl,CID_PARTICLE,p);
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


PQ6 *pqSearchPsd(PSX smx,PQ6 *pq,FLOAT *r,FLOAT *v, FLOAT *rscale,FLOAT *vscale, int bReplica,int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    int idSelf = smx->pkd->idSelf;
    FLOAT *Smin = smx->SminT;
    int *S = smx->ST;
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
	    pq = pqSearchLocalPsd(smx,pq,r,v,rscale,vscale,pbDone);
	    if (*pbDone) return(pq);
	}
	else {
	    pq = pqSearchRemotePsd(smx,pq,id,r,v,rscale,vscale);
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
** Used to keep track of the state of things from a previous call to knn()
*/
struct knn_data
{
    FLOAT rscale[3], vscale[3];
    FLOAT rLast[3], vLast[3];
};

/*
** Find the k nearest neighbors. If first_time is true the priority-queue
** and knn_data is initialized. Subsequent calls with first_time false
** will update the priority-queue.
*/
PQ6 *knn(PSX psx, PQ6 *pq, int pid, struct knn_data *knn_data, int first_time)
{
    int i,j;
    PKD pkd = psx->pkd;
    FLOAT *r,*v;
    PARTICLE *p;

    if (first_time)
    {
	/*
	** Initialize the bInactive flags for all local particles.
	*/
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    psx->ea[i].bInactive = (p->bSrcActive)?0:1;
	}
	psx->ea[pkd->nLocal].bInactive = 0;  /* initialize for Sentinel, but this is not really needed */

	/*
	** Initialize the priority queue first.
	*/
	for (i=0;i<psx->nSmooth;++i) {
	    psx->pq[i].pPart = &psx->pSentinel;
	    psx->pq[i].iIndex = pkd->nLocal;
	    psx->pq[i].iPid = pkd->idSelf;
	    psx->pq[i].dr[0] = HUGE_VAL;
	    psx->pq[i].dr[1] = HUGE_VAL;
	    psx->pq[i].dr[2] = HUGE_VAL;
	    psx->pq[i].dr[3] = HUGE_VAL;
	    psx->pq[i].dr[4] = HUGE_VAL;
	    psx->pq[i].dr[5] = HUGE_VAL;
	    psx->pq[i].fDist2 = HUGE_VAL;
	    //psx->pq[i].fDist2 = 0;
	    //for (j=0; j < 6; j++) psx->pq[i].fDist2 += pow(psx->pq[i].dr[j],2);
	}
	for (j=0;j<3;++j) knn_data->rLast[j] = 0.0;
	for (j=0;j<3;++j) knn_data->vLast[j] = 0.0;
    }
    else
    {
	for (i=0;i<psx->nSmooth;++i)
	{
	    for (j=0; j < 3; j++) psx->pq[i].dr[j]   /= knn_data->rscale[j];
	    for (j=0; j < 3; j++) psx->pq[i].dr[3+j] /= knn_data->vscale[j];
	}
    }

    p = pkdParticle(pkd, pid);

    r = p->r;
    v = pkdVel(pkd, p);

    for (j=0;j<3;++j) knn_data->rscale[j] = PSM(pid).rscale[j];
    for (j=0;j<3;++j) knn_data->vscale[j] = PSM(pid).vscale[j];

    if (!first_time)
    {
	for (i=0;i<psx->nSmooth;++i) 
	{
	    psx->pq[i].fDist2 = 0;
	    for (j=0; j < 3; j++) 
	    {
		if (knn_data->rscale[j] == 0)
		{
		    psx->pq[i].dr[j] = 0;
		}
		else
		{
		    psx->pq[i].dr[j]  -= r[j] - knn_data->rLast[j];
		    psx->pq[i].dr[j]  *= knn_data->rscale[j];
		}
		psx->pq[i].fDist2 += pow(psx->pq[i].dr[j],2);
	    }

	    for (j=0; j < 3; j++) 
	    {
		if (knn_data->vscale[j] == 0)
		{
		    psx->pq[i].dr[3+j] = 0;
		}
		else
		{
		    psx->pq[i].dr[3+j]  -= v[j] - knn_data->vLast[j];
		    psx->pq[i].dr[3+j]  *= knn_data->vscale[j];
		}
		psx->pq[i].fDist2 += pow(psx->pq[i].dr[3+j],2);
	    }
	}
    }

    for (j=0;j<3;++j) knn_data->rLast[j] = r[j];
    for (j=0;j<3;++j) knn_data->vLast[j] = v[j];
    PQ6_BUILD(psx->pq,psx->nSmooth,pq);
    int bDone;
    return pqSearchPsd(psx,pq,r,v,knn_data->rscale,knn_data->vscale,0,&bDone);
}


#if 0
void psdSmooth(PSX smx, PSF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;
    int N;
    printf("psdSmooth (file)\n");

    //FILE *fp = fopen("/home/itp/jonathan/zbox3b/MAD/MockHalos/Hernquist/main.den", "r"); assert(fp != NULL);
    //FILE *fp = fopen("subsubhalo_nfw.den6", "r"); assert(fp != NULL);
    FILE *fp = fopen("A-5.den6", "r"); assert(fp != NULL);
    //FILE *fp = fopen("test.den", "r"); assert(fp != NULL);
    fscanf(fp, "%i\n", &N);
    assert(N == pkd->nLocal);
    float *den = malloc(N * sizeof(*den)); assert(den != NULL);
    for (pi=0;pi<pkd->nLocal;++pi) {
	fscanf(fp, "%f\n", &den[pi]);
    }
    fclose(fp);

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	p->fDensity = den[p->iOrder];
    }
    free(den);

#if 1
    //fp = fopen("test.ball", "r"); assert(fp != NULL);
    //fp = fopen("/home/itp/jonathan/zbox3b/MAD/MockHalos/Hernquist/main.ball", "r"); assert(fp != NULL);
    fp = fopen("A-5.ball", "r"); assert(fp != NULL);
    fscanf(fp, "%i\n", &N);
    assert(N == pkd->nLocal);
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	fscanf(fp, "%f\n", &p->fBall);
    }
    fclose(fp);
#endif
    printf("psdSmooth (file) finished\n");
}

#else


/*
** Compute the smoothed phase-space densities for all particles.
*/
void psdSmooth(PSX psx,  PSF *smf) {
    PARTICLE *p;
    PKD pkd = psx->pkd;
    PQ6 *pq = psx->pq;
    FLOAT fBall;
    int pi,i,bDone=0;

    struct knn_data knn_data;

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);

    for (pi=0;pi<pkd->nLocal;++pi) 
    {
	if (pi % 10000 == 0)
	    fprintf(stdout, "\033[%iC%3i\r", pkd->idSelf*4, (int)(100.*(((double)pi) / pkd->nLocal)));
	p = pkdParticle(pkd,pi);
	/* if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue; */

	pq = knn(psx,pq, pi, &knn_data, pi==0);

	/*
	** Search in replica boxes if it is required.
	*/
	if (!bDone && psx->bPeriodic) {
	    fBall = sqrt(pq->fDist2);
	    /* scaling missing here... */
	    assert(0);
#if 0
	    for (j=0;j<3;++j) {
		iStart[j] = floor((r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j]   = floor((r[j] + fBall)/pkd->fPeriod[j] + 0.5);
		iStart[3+j] = floor((v[j] - fBall)/pkd->fPeriod[3+j] + 0.5);
		iEnd[3+j]   = floor((v[j] + fBall)/pkd->fPeriod[3+j] + 0.5);
	    }
	    void pqSearchRecur(int i, int not_center) /* inline recursive function */
	    {
		int j;
		if (i == 6) {
		    if (not_center) pq = pqSearchPsd(psx,pq,r2,rscale,vscale,1,&bDone);
		    } 
		else {

		    for (j=iStart[i]; j <= iEnd[i]; j++) {
			r2[i] = r[i] - j*pkd->fPeriod[i];
			pqSearchRecur(i+1, not_center || j);
			}
		    }
	    }
	    pqSearchRecur(0,0);
#endif
	}
	p->fBall = sqrt(pq->fDist2);
	/*
	** Apply smooth function to the neighbor list.
	*/
	p->fDensity = psdDensity(pkd, p,psx->nSmooth,psx->pq,PSM(pi).rscale,PSM(pi).vscale);
	/*
	** Call mdlCacheCheck to make sure we are making progress!
	*/
	mdlCacheCheck(pkd->mdl);
    }
    fprintf(stdout, "\033[%iC    \r", pkd->idSelf*4);

    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<psx->nSmooth;++i) {
	if (psx->pq[i].iPid == pkd->idSelf) {
	    psx->ea[psx->pq[i].iIndex].bInactive = 0;
	}
	else {
	    psHashDel(psx,psx->pq[i].pPart);
	    mdlRelease(pkd->mdl,CID_PARTICLE,psx->pq[i].pPart);
	}
    }

    mdlFinishCache(pkd->mdl,CID_PARTICLE);

}
#endif

/*
** For each of the neighbors stored in psx->pq, compute the arc length
** between the vector a and the position of the neighbor.
*/
void calc_arclen(PSX psx, double *a, double *arclen)
{
    int i,pj;
    double L = 0;
    double b[6];

    for (pj=0; pj < psx->nSmooth; pj++)
    {
	L = 0;
	for (i=0; i < 6; i++)
	{
	    b[i] = psx->pq[pj].dr[i];
	    L += pow(b[i], 2);
	}
	L = sqrt(L);
	for (i=0; i < 6; i++)
	    b[i] /= L;

	arclen[pj] = 0;
	for (i=0; i < 6; i++)
	    arclen[pj] += a[i] * b[i];
	arclen[pj] = acos(arclen[pj]);
	arclen[pj] *= L;
    }
}

/*
** Form local groups by forming chains of particles that end at local density maxima.
**
** Each particle links to its neighbor that lies best along the density gradient and
** is denser. Best is defined as having the smallest arc length between the gradient
** vector and the relative neighbor position vector.
**
** If a neighbor already belongs to a chain then the current chain attaches to it.
** If no neighbor is denser then the chain terminates.
** A chain may also terminate if a neighbor exists but is on another processor. In
** this case, the terminal particle and its neighbor form a bridge which will be 
** joined at a later stage after the local groups are built.
*/
void psdSmoothLink(PSX psx, PSF *smf) {
    PARTICLE *p;
    PKD pkd = psx->pkd;
    PQ6 *pq = psx->pq;
    int64_t pi,i;
    int64_t pj;
    int64_t idx;
    int32_t *piGroup;

#define TEMP_S_INCREASE 100
    int *C; NEW_STACK(C, TEMP_S_INCREASE);
    int *G; NEW_STACK(G, TEMP_S_INCREASE);

    int sorted_nbrs[psx->nSmooth];
    double arclen[psx->nSmooth];
    double fDensityGrad[6];

    struct bridge_info *B; NEW_STACK(B, TEMP_S_INCREASE);
    struct bridge_info bi;

    int arclen_compar(const void *a0, const void *b0)
    {
	int a = *(int *)a0;
	int b = *(int *)b0;
	if (arclen[a] < arclen[b]) return -1;
	if (arclen[a] > arclen[b]) return +1;
	return 0;
    }

    int iorder_compar(const void *a0, const void *b0)
    {
	int a = *(int *)a0;
	int b = *(int *)b0;
	if (psx->pq[a].pPart->iOrder < psx->pq[b].pPart->iOrder) return -1;
	if (psx->pq[a].pPart->iOrder > psx->pq[b].pPart->iOrder) return +1;
	return 0;
    }

    int group_size_compar(const void *a0, const void *b0)
    {
	PSGD *a = (PSGD *)a0;
	PSGD *b = (PSGD *)b0;
	if (a->nTotal > b->nTotal) return -1;
	if (a->nTotal < b->nTotal) return +1;
	return 0;
    }

#if 0
    float den_max = 0;
    for (pi=0;pi<pkd->nLocal;pi++) 
    {
	if (pkdParticle(pkd, pi)->fDensity > den_max)
	    den_max = pkdParticle(pkd, pi)->fDensity;
    }

#endif
    for (pi=0;pi<pkd->nLocal;pi++) 
    {
	p = pkdParticle(pkd,pi);
	*pkdGroup(pkd,p) = 0;
	//p->fDensity /= den_max;
    }

    pkd->nGroups = 0;
    pkd->nMaxRm = 0;

    psx->nBridges = 0;

    int nGroups = 1;
    int32_t trial_group;
    int first_time = 1;
    struct knn_data knn_data;

#if 0
    char fname[256];
    sprintf(fname, "link.%i", pkd->idSelf);
    FILE *fp = fopen(fname, "w");
#endif

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);

    int nPeaks = 0;
    int nSingles = 0;
    for (idx=0;idx<pkd->nLocal;++idx) 
    {
	pi = idx;

	if (idx % 10000 == 0)
	    fprintf(stdout, "\033[%iC%3i\r", pkd->idSelf*4, (int)(100.*(((double)idx) / pkd->nLocal)));

	PARTICLE *p0 = pkdParticle(pkd,pi);
	if (*pkdGroup(pkd,p0) != 0) continue;

	int chain_len = 0;

	trial_group = -1;

	assert(STACK_EMPTY(C));
	while (trial_group == -1)
	{
	    EXTEND_STACK(C);

	    PUSH(C,pi);
	    chain_len++;

	    p = pkdParticle(pkd,pi);

	    assert(*pkdGroup(pkd,p) != -1);
	    *pkdGroup(pkd,p) = -1;

	    /* Find our neighbors */
	    pq = knn(psx,pq, pi, &knn_data, first_time);
	    first_time = 0;

	    for (pj=0; pj < psx->nSmooth; pj++)
		sorted_nbrs[pj] = pj;

#if 0
	    qsort(sorted_nbrs, psx->nSmooth, sizeof(*sorted_nbrs), iorder_compar);
	    fprintf(fp, "NBR %ld", p->iOrder);
	    for (pj=0; pj < psx->nSmooth; pj++)
		fprintf(fp, " %ld", psx->pq[sorted_nbrs[pj]].pPart->iOrder);
	    fprintf(fp, "\n");
#endif

	    //psdDensityGrad(pkd, psx, pi, fDensityGrad);
	    psdPotentialGrad(pkd, psx, pi, fDensityGrad);
	    calc_arclen(psx, fDensityGrad, arclen);
	    qsort(sorted_nbrs, psx->nSmooth, sizeof(*sorted_nbrs), arclen_compar);

	    int bridge = 0;

	    /* Find the first sorted particle with a density greater than ours */
	    float max_den=PROP(pkd, pkdParticle(pkd, pi));
	    PARTICLE *next_p = NULL;
	    int max_den_j=-1;
	    PQ6 *nbr;
	    for (pj=0; pj < psx->nSmooth; pj++)
	    {
		nbr = psx->pq + sorted_nbrs[pj];

		if (PROP(pkd, nbr->pPart) > max_den)
		{
#if 0
		    max_den = nbr->pPart->fDensity;
		    max_den_j = pj;
		}
	    }

	    while (max_den_j != -1)
	    {
		pj = max_den_j;
		nbr = psx->pq + sorted_nbrs[pj];
		{

#endif
		    bridge = nbr->iPid != pkd->idSelf;

		    if (bridge)
		    {
			bi.iPid   = nbr->iPid;
			bi.iIndex = nbr->iIndex;
		    }
		    else
		    {
			int32_t nbr_grp = *pkdGroup(pkd, nbr->pPart);

			/* This would be false if we formed a loop. We shouldn't be able to form loops. */
			assert(nbr_grp != -1);

			if (nbr_grp > 0) /* We've found a local chain. Take that group. */
			    trial_group = nbr_grp;

			next_p = nbr->pPart;
			pi = nbr->iIndex;
		    }

		    /* fprintf(fp, "LINK %ld %ld\n", p->iOrder, nbr->pPart->iOrder); */

		    if (pkd->oAcceleration)
		    {
			pkdAccel(pkd, p)[1] = nbr->pPart->r[0] - p->r[0];
			pkdAccel(pkd, p)[2] = nbr->pPart->r[1] - p->r[1];
			pkdAccel(pkd, p)[0] = nbr->pPart->r[2] - p->r[2];
		    }

		    break;
		}
	    }

	    /* We didn't find a new particle to link to. We must be at a peak. */
	    if (next_p == NULL)
	    {
		if (!bridge)
		{
		    /*
		    ** If we are at a peak and the chain is only one particle long, look down the gradient
		    ** and find the first particle with the density less than ours that belongs to a group.
		    */
		    if (chain_len == 1)
			nSingles++;
#if 0
		    if (chain_len == 1)
		    {
			float max_den=pkdParticle(pkd, pi)->fDensity;
			for (pj=0; pj < psx->nSmooth; pj++)
			{
			    PARTICLE *nbr = psx->pq[sorted_nbrs[pj]].pPart;

			    if (nbr->fDensity < max_den)
			    {
				bridge = psx->pq[sorted_nbrs[pj]].iPid != pkd->idSelf;

				int32_t nbr_grp = *pkdGroup(pkd, nbr);

				if (bridge)
				{
				    bi.Pid = psx->pq[sorted_nbrs[pj]].iPid;
				    bi.pid = psx->pq[sorted_nbrs[pj]].iIndex;
				    bi.done = 0;
				}
				else
				{
				    // Do we really ever go into this if?
				    if (nbr_grp > 0)
					trial_group = nbr_grp;
				}

				break;
			    }
			}
		    }
#endif

		    /* We still couldn't find something to link to. Just create a new group. */
		    if (trial_group == -1 && !bridge)
		    {
			trial_group = nGroups++;
			EXTEND_STACK(G);
			PUSH(G,pi);
		    }

		    //fprintf(stdout, "PEAK!\n");
		    nPeaks++;
#if 0
		    pkdAccel(pkd, p)[0] = 0;
		    pkdAccel(pkd, p)[1] = 0;
		    pkdAccel(pkd, p)[2] = 0;
#endif
		}


		if (bridge)
		{
		    trial_group = nGroups++;
		    EXTEND_STACK(G);
		    PUSH(G,pi);

		    bi.local_gid = trial_group;
		    bi.remote_gid = -1;
		    bi.done = 0;
		    EXTEND_STACK(B);
		    PUSH(B, bi);
		    psx->nBridges++;
		}
	    }

#if 0
	    pkdAccel(pkd, p)[0] = fDensityGrad[0];
	    pkdAccel(pkd, p)[1] = fDensityGrad[1];
	    pkdAccel(pkd, p)[2] = fDensityGrad[2];
#endif

	    /*
	    ** Call mdlCacheCheck to make sure we are making progress!
	    */
	    mdlCacheCheck(pkd->mdl);
	}

	assert(trial_group != -1);
	assert(trial_group != 0);

	/* Assign particles to the (new) group */
	while (!STACK_EMPTY(C))
	{
	    int pid = POP(C);
	    *pkdGroup(pkd, pkdParticle(pkd,pid)) = trial_group;
	}
    }
    fprintf(stdout, "\033[%iC    \r", pkd->idSelf*4);

#if 0
    fclose(fp);
#endif
#if 0
    fprintf(stdout, "%i] nPeaks is %i\n",   pkd->idSelf, nPeaks);
    fprintf(stdout, "%i] nBridges is %i\n",   pkd->idSelf, psx->nBridges);
    fprintf(stdout, "%i] nSingles is %i\n",  pkd->idSelf,  nSingles);
    fprintf(stdout, "%i] nGroups is %i\n",  pkd->idSelf,  nGroups);
#endif

    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<psx->nSmooth;++i) {
	if (psx->pq[i].iPid == pkd->idSelf) {
	    psx->ea[psx->pq[i].iIndex].bInactive = 0;
	}
	else {
	    psHashDel(psx,psx->pq[i].pPart);
	    mdlRelease(pkd->mdl,CID_PARTICLE,psx->pq[i].pPart);
	}
    }

    mdlFinishCache(pkd->mdl,CID_PARTICLE);

    /*
    ** Get rid of the stack of bridges and store them in the psx context.
    */
    psx->bridges = malloc(psx->nBridges * sizeof(*(psx->bridges)));
    i=0;
    while (!STACK_EMPTY(B)) { psx->bridges[i++] = POP(B); }
    assert(STACK_EMPTY(B));

    /* When will we free this? */
    pkd->nGroups = nGroups;
    pkd->psGroupData = mdlMalloc(pkd->mdl, pkd->nGroups * sizeof(PSGD));

    /*
    ** Create the local group table
    */
    memset(pkd->psGroupData, 0, pkd->nGroups * sizeof(PSGD));
    PSGD *gd = pkd->psGroupData;

    gd[0].iLocalId = 0;
    gd[0].iPid = pkd->idSelf;
    gd[0].fDensity = 0;
    gd[0].nTotal = pkd->nLocal;
    gd[0].bridge = 0;

    for (i=nGroups-1; i > 0; i--)
    {
	int pid = POP(G);
	PARTICLE *p = pkdParticle(pkd,pid);
	assert(*pkdGroup(pkd,p) == i);
	gd[i].iGlobalId = i;
	gd[i].iLocalId = i;
	gd[i].iPid = pkd->idSelf;
	gd[i].fDensity = PROP(pkd, p);
	gd[i].r[0] = p->r[0];
	gd[i].r[1] = p->r[1];
	gd[i].r[2] = p->r[2];
	gd[i].fMass += pkdMass(pkd, p);
	gd[i].bridge = 0;
	gd[i].dup = 0;
    }
    assert(STACK_EMPTY(G));

    for (i=0; i < psx->nBridges; i++)
	gd[psx->bridges[i].local_gid].bridge = 1;

    FREE_STACK(C);
    FREE_STACK(G);
    FREE_STACK(B);
}

/*
** Join groups that span multiple processors. This will be called many times from master.
**
** Group information is propagated down the chains from the processor that owns the group.
** Ownership is defined as the processor where the density peak is located. The group id
** from the owner defines the group id on other processors.
**
** One edge case occurs when a group chain returns to the owners processor (possible many times).
** This creates a duplicate entry in the group table. Such occurences are noted and later
** the particles in duplicated groups will be reassigned to a single group.
**
** Group information within the particles is not considered here except to retrieve the group id
** of the remote group in a bridge.
*/
int psdJoinBridges(PSX psx, PSF *smf) {
    PKD pkd = psx->pkd;
    int done = 1;
    int i,j;
    MDL mdl = psx->pkd->mdl;
    PARTICLE *p, *p_remote;
    struct bridge_info *bi = psx->bridges;
    PSGD *gd = pkd->psGroupData;

    struct store
    {
	int i;
	PSGD gd;
    };

#define TEMP_S_INCREASE 100
    struct store store;
    struct store *S; NEW_STACK(S, TEMP_S_INCREASE);

    assert(pkd->oGroup);

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd,i);
	if (*pkdGroup(pkd, p) == 0)
	{
	    fprintf(stderr, "%i] Particle %ld has group 0\n", pkd->idSelf, i);
	    assert(0);
	}
    }


    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);
    mdlROcache(mdl,CID_GROUP,NULL,pkd->psGroupData,sizeof(PSGD), pkd->nGroups);

    int nBridgesLeft = 0;

    /*
    ** Find once the group id of the remote particle in a bridge.
    */
    for (i=0; i < psx->nBridges; i++)
    {
	if (!bi[i].done && bi[i].remote_gid == -1)
	{
	    assert(bi[i].iPid != pkd->idSelf);
	    p_remote = mdlAquire(mdl,CID_PARTICLE,bi[i].iIndex,bi[i].iPid);
	    bi[i].remote_gid = *pkdGroup(pkd, p_remote);
	    if (bi[i].remote_gid == 0)
	    {
		fprintf(stderr, "%i] %i %i %ld id=%llu\n",
		    pkd->idSelf, bi[i].remote_gid, bi[i].iPid, bi[i].iIndex,
		    p_remote->iOrder); 
	    }
	    assert(bi[i].remote_gid != 0);
	    mdlRelease(mdl,CID_PARTICLE,p_remote);
	}
    }

    /*
    ** Copy the remote group information from higher up the chain. Store this on a stack rather
    ** than overwrite our local data because another processor might be trying to read it.
    */
    for (i=0; i < psx->nBridges; i++)
    {
	if (bi[i].done) continue;

	PSGD *new_gd = mdlAquire(mdl, CID_GROUP, bi[i].remote_gid, bi[i].iPid);
	store.i = i;
	store.gd = *new_gd;
	EXTEND_STACK(S);
	PUSH(S, store);
	if (store.gd.iLocalId == 0)
	{
	    fprintf(stderr, "%i] store.iLocalId=%i  %i %i %ld\n", pkd->idSelf, store.gd.iLocalId, bi[i].remote_gid, bi[i].iPid, bi[i].iIndex); 
	}
	assert(store.gd.iLocalId != 0);
	mdlRelease(mdl,CID_GROUP,new_gd);
    }

    mdlFinishCache(mdl,CID_GROUP);
    mdlFinishCache(mdl,CID_PARTICLE);

    /*
    ** Now that the cache is closed, we can safely update our local group table.
    */
    while (!STACK_EMPTY(S))
    {
	store = POP(S);
	i = store.i;

	assert(!bi[i].done);

	gd = pkd->psGroupData + bi[i].local_gid;

	if(gd->bridge == 0 && store.gd.bridge == 1)
	    assert(0);

	*gd = store.gd;
	gd->dup = (gd->iPid == pkd->idSelf) * gd->iLocalId;

	bi[i].done = !gd->bridge;
	done = done && bi[i].done;

	nBridgesLeft += !bi[i].done;
    }


#if 0
    for (i=0; i < psx->nBridges; i++)
    {
	if (bi[i].done) continue;

	gd = pkd->psGroupData + bi[i].local_gid;

	if(gd->bridge == 0 && store[i].bridge == 1)
	    assert(0);

	*gd = store[i];
	gd->dup = (gd->iPid == pkd->idSelf) * gd->iLocalId;

	bi[i].done = !gd->bridge;
	done = done && bi[i].done;

	nBridgesLeft += !bi[i].done;
    }
#endif

    FREE_STACK(S);

#if 0
    if (nBridgesLeft > 0)
	fprintf(stdout, "%i] %i bridges left\n", pkd->idSelf, nBridgesLeft);
#endif

    return done;
}

/*
** Count unique local groups. Do not count duplicates that will later be removed. 
*/
int psdCountLocalGroups(PSX psx) {
    PKD pkd = psx->pkd;
    int i;
    PSGD *gd = pkd->psGroupData;

    int count = 0;
    for (i=1; i < pkd->nGroups; i++)
    {
	assert(gd[i].iLocalId != 0);
	if (!gd[i].dup && gd[i].iPid == pkd->idSelf)
	    count++;
    }

    return count;
}

static void initGroups(void *vpkd, void *a)
{
    PSGD *g = (PSGD *)a;
    g->fRMSRadius = 0;
    g->fMass = 0;
    g->nTotal = 0;
    g->nLocal = 0;
}

static void combGroups(void *vpkd, void *a, void *b)
{
    PSGD * g1 = (PSGD *)a;
    PSGD * g2 = (PSGD *)b;
    g1->nTotal += g2->nTotal;
    g1->fMass  += g2->fMass;
    if (g2->fRMSRadius > g1->fRMSRadius)
	g1->fRMSRadius = g2->fRMSRadius;
}

/*
** Update the local group table with global ids. The offset into the range of globally
** unique ids comes from master.
*/
void psdUpdateGroups(PSX psx, int offs, int count)
{
    PKD pkd = psx->pkd;
    int64_t i;
    PSGD *gd = pkd->psGroupData;

    int new_local_id=1;

    struct gdstore
    {
	int i;
	PSGD gd;
    };

#define TEMP_S_INCREASE 100
    struct gdstore *G; NEW_STACK(G, TEMP_S_INCREASE);
    struct gdstore gdstore;

    /*
    ** Update local groups with a unique global id
    */
    for (i=1; i < pkd->nGroups; i++)
    {
	if (!gd[i].dup && gd[i].iPid == pkd->idSelf)
	{
	    gd[i].iGlobalId = new_local_id + offs;
	    new_local_id++;
	}
    }
    assert(new_local_id == count+1);

    /*
    ** To avoid hassles later, we fix up particles that belong to 
    ** duplicate group entries.
    */
    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);

	if (gd[gid].dup != 0)
	{
	    *pkdGroup(pkd, p) = gd[gid].dup;
	}
    }

    /*************/
    /*************/

    /*
    ** Now bring over the global ids from remote groups. We use a stack to save the data
    ** so that the cache doesn't see half updated data.
    */

    //fprintf(stderr, "%i] nGroups %i\n", pkd->idSelf, pkd->nGroups);
    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupData,sizeof(PSGD), pkd->nGroups);
    for (i=1; i < pkd->nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf)
	{
	    assert(gd[i].dup == 0);
	    PSGD *new_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[i].iLocalId, gd[i].iPid);
	    assert(new_gd->iPid == gd[i].iPid);
	    assert(new_gd->dup == 0);
	    gdstore.i = i;
	    gdstore.gd = *new_gd;
	    EXTEND_STACK(G);
	    PUSH(G, gdstore);
	    mdlRelease(pkd->mdl,CID_GROUP,new_gd);
	}
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

    while (!STACK_EMPTY(G))
    {
	gdstore = POP(G);
	gd[gdstore.i] = gdstore.gd;
    }

    /*************/
    /*************/

#if 1

    /*
    ** Compute some group quantities like total mass and radius.
    */
    for (i=1; i < pkd->nGroups; i++)
    {
	gd[i].nTotal = 0;
	gd[i].nLocal = 0;
	gd[i].fMass = 0;
	gd[i].fRMSRadius = 0;
    }

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);

	gd[gid].nLocal++;
	gd[gid].nTotal++;
	gd[gid].fMass += pkdMass(pkd,p);

	double d = pow(gd[gid].r[0] - p->r[0],2)
		 + pow(gd[gid].r[1] - p->r[1],2)
		 + pow(gd[gid].r[2] - p->r[2],2);

	d = sqrt(d);

	//fprintf(stderr, "@@ %i/%i %e %e\n", gid, pkd->nGroups, d, gd[gid].fRMSRadius);
	if (d > gd[gid].fRMSRadius)
	    gd[gid].fRMSRadius = d;
    }

    /*
    ** Combine local data across domains.
    */
    mdlCOcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupData,sizeof(PSGD), pkd->nGroups,pkd,initGroups,combGroups);
    //fprintf(stderr, "%i] nGroups %i\n", pkd->idSelf, pkd->nGroups);

    for (i=1; i < pkd->nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf)
	{
	    assert(gd[i].dup == 0);
	    PSGD *remote_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[i].iLocalId, gd[i].iPid);
	    remote_gd->nTotal = gd[i].nTotal;
	    remote_gd->fMass  = gd[i].fMass;
	    if (gd[i].fRMSRadius > remote_gd->fRMSRadius)
		remote_gd->fRMSRadius = gd[i].fRMSRadius;
	    mdlRelease(pkd->mdl,CID_GROUP,remote_gd);
	}
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

#endif

#if 0
    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupData,sizeof(PSGD), pkd->nGroups);

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);
	if (gd[gid].dup != 0)
	    gid = gd[gid].dup;

	double d = pow(gd[gid].r[0] - p->r[0],2) 
		 + pow(gd[gid].r[1] - p->r[1],2)
		 + pow(gd[gid].r[2] - p->r[2],2);

	double fRMSRadius;
	if (gd[gid].iPid == pkd->idSelf)
	{
	    fRMSRadius = gd[gid].fRMSRadius;
	}
	else
	{
	    PSGD *remote_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[gid].iLocalId, gd[gid].iPid);
	    fRMSRadius = remote_gd->fRMSRadius;
	    mdlRelease(pkd->mdl,CID_GROUP,remote_gd);
	}

	p->fDensity = 1 - pow(d / fRMSRadius,2);
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

#endif

    /*
    ** Update particles with their global group id.
    */
    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
#if 0
	if (gd[*pkdGroup(pkd, p)].fDensity / p->fDensity > 200)
	    *pkdGroup(pkd, p) = 0;
	else
#endif
	    *pkdGroup(pkd, p) = gd[*pkdGroup(pkd, p)].iGlobalId;
	//*pkdGroup(pkd, p) = *pkdGroup(pkd, p) % 256;
    }

    FREE_STACK(G);
}

