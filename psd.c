#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_PSD

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "psd.h"
#include "pkd.h"
#include "rbtree.h"
#include <sys/stat.h>

const char *psd_module_id = "$Id$";
const char *psd_h_module_id = PSD_H_MODULE_ID;

extern const int primes[];

void initPsdDensity(void *vpkd, void *p) {
    ((PARTICLE *)p)->fDensity = 0.0;
    }

void combPsdDensity(void *vpkd, void *p1,void *p2) {
    ((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
    }

static inline float ekernel(float u2) {
    return (0 <= u2 && u2 <= 1) * (1-u2);
}

void psdDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf, FLOAT scale[6]) {
    PKD pkd = smf->pkd;
    float ih2,r2,fDensity,fMass;
    int i,j;

    ih2 = BALL2(p); //4.0/BALL2(p);
    fDensity = 0.0;

    for (i=0;i<nSmooth;++i) {
        //assert(nnList[i].fDist2 > 0);
        fMass = pkdMass(pkd,nnList[i].pPart);
        assert(fMass > 0);
        r2 = nnList[i].fDist2/ih2;
        if (r2 == 0) r2 = p->fBall * 6 / (1+6.);
        fDensity += ekernel(r2) * fMass;
        }

    float V0=1, V1=1;
    for (i=0; i < 6; i++) {
        V0 *= p->fBall;
        V1 *= scale[i];
    }

    assert(p->fBall > 0);
    p->fDensity = 0.77403670 * fDensity * V1 / V0;
    if (isnan(p->fDensity))
        fprintf(stderr, "%g %g %g %g \n", p->fBall, fDensity, V1, V0);

    assert(!isnan(p->fDensity));
    assert(!isinf(p->fDensity));

    }

void psdDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    float fNorm,ih2,r2,rs,fMassQ,fMassP;
    int i;
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        fMassQ = pkdMass(pkd,q);
        p->fDensity += rs*fMassQ;
        q->fDensity += rs*fMassP;
        }
    }


int psdInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType) {
    SMX smx;
    void (*initParticle)(void *,void *) = NULL;
    void (*init)(void *,void *) = NULL;
    void (*comb)(void *,void *,void *) = NULL;
    int i,pi,j;
    int nTree;
    int iTopDepth;

    smx = malloc(sizeof(struct smContext));
    assert(smx != NULL);
    smx->pkd = pkd;
    if (smf != NULL) smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->bPeriodic = bPeriodic;
    /*
    ** Initialize the context for compressed nearest neighbor lists.
    */
    smx->lcmp = lcodeInit(pkd->nThreads,pkd->idSelf,pkd->nLocal,nSmooth);

    switch (iSmoothType) {
    case PSD_DENSITY:
        //smx->fcnSmooth = bSymmetric ? psdDensitySym : psdDensity;
        initParticle = initPsdDensity; /* Original Particle */
        init = initPsdDensity; /* Cached copies */
        comb = combPsdDensity;
        smx->fcnPost = NULL;
        break;
#if 0
    case PSD_DENSITY:
        smx->fcnSmooth = bSymmetric ? psdDensitySym : psdDensity;
        initParticle = initPsdDensity; /* Original Particle */
        init = initPsdDensity; /* Cached copies */
        comb = combPsdDensity;
        smx->fcnPost = NULL;
        break;
#endif
    case PSD_FOF:
        assert(bSymmetric == 0);
        smx->fcnSmooth = NULL;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    default:
        assert(0);
    }
    /*
    ** Initialize the ACTIVE particles in the tree.
    ** There are other particles in the tree -- just not active.
    */
    nTree = pkdTreeNode(pkd,ROOT)->pUpper + 1;
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
    */
    smx->nnListSize = 0;
    smx->nnListMax = NNLIST_INCREMENT;
    smx->nnList = malloc(smx->nnListMax*sizeof(NN));
    assert(smx->nnList != NULL);

    /*
    ** Allocate priority queue.
    */
    smx->pq = malloc(nSmooth*sizeof(PQ));
    assert(smx->pq != NULL);
    PQ_INIT(smx->pq,nSmooth);
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
    assert(sizeof(PLITE) >= sizeof(struct smExtraArray));
    smx->ea = UNION_CAST(pkd->pLite,PLITE *,struct smExtraArray *);
    *psmx = smx;
    return(1);
}


void psdFinish(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;
    char achOut[128];

    printf("psdFinish\n");

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
                smx->fcnPost(pkd,p,smf);
        }
    }
    /*
    ** Finish compressed lists.
    */
    lcodeFinish(smx->lcmp);
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

    printf("psdFinish Finished\n");
}


/*
** This function performs a local nearest neighbor search.
*/
PQ *pqSearchLocalPsd(SMX smx,PQ *pq,FLOAT r[6],FLOAT scale[6], int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    double *v;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    FLOAT dMin,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int i,j,pj,pEnd,iCell,iParent;
    int sp = 0;
    int sm = 0;
    int idSelf = pkd->idSelf;
    pBND bnd;

    *pbDone = 1;        /* assume that we will complete the search */
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
            bnd = pkdNodeBnd(pkd, kdn);
            MINDIST6(bnd,r,min1,scale);
            kdn = pkdTreeNode(pkd,++iCell);
            bnd = pkdNodeBnd(pkd, kdn);
            MINDIST6(bnd,r,min2,scale);
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
            if (smx->ea[pj].bInactive) continue;
            p = pkdParticle(pkd,pj);
            v = pkdVel(pkd,p);
            dr0 = (r[0] - p->r[0]) * scale[0];
            dr1 = (r[1] - p->r[1]) * scale[1];
            dr2 = (r[2] - p->r[2]) * scale[2];
            dr3 = (r[3] -    v[0]) * scale[3];
            dr4 = (r[4] -    v[1]) * scale[4];
            dr5 = (r[5] -    v[2]) * scale[5];
            fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
            if (fDist2 < pq->fDist2) {
                if (pq->iPid == idSelf) {
                    smx->ea[pq->iIndex].bInactive = 0;
                } 
                else {
                    //smHashDel(smx,pq->pPart);
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
                PQ_REPLACE(pq);
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
                bnd = pkdNodeBnd(pkd, kdn);
                for (j=0;j<6;++j) {
                    dMin = (bnd.fMax[j] - fabs(bnd.fCenter[j] - r[j])) * scale[j];
                    if (dMin*dMin < pq->fDist2 || dMin < 0) {
                        iParent = kdn->iParent;
                        if (!iParent) {
                            *pbDone = 0;                /* EXIT, not contained! */
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
            bnd = pkdNodeBnd(pkd, kdn);
            MINDIST6(bnd,r,min2,scale);
        }
        if (min2 >= pq->fDist2) {
            kdn = pkdTreeNode(pkd,iCell = kdn->iParent);
            goto NoIntersect;
        }
        S[++sp] = iCell;
    }
}



PQ *pqSearchRemotePsd(SMX smx,PQ *pq,int id,FLOAT r[6], FLOAT scale[6]) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *kdn;
    PARTICLE *p;
    KDN *pkdn,*pkdu;
    double *v;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    FLOAT min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int pj,pEnd,iCell;
    int sp = 0;
    int sm = 0;
    int idSelf = pkd->idSelf;
    pBND bnd;

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
            mdlRelease(mdl,CID_CELL,pkdn);
            pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
            bnd  = pkdNodeBnd(pkd,pkdn);
            MINDIST6(bnd,r,min1,scale);
            kdn  = pkdTreeNode(pkd,++iCell);
            pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
            bnd  = pkdNodeBnd(pkd,pkdu);
            MINDIST6(bnd,r,min2,scale);
            if (min1 < min2) {
                Smin[sm++] = min2;
                kdn = pkdTreeNode(pkd,--iCell);
                mdlRelease(mdl,CID_CELL,pkdu);
                if (min1 >= pq->fDist2) goto NotContained;
            }
            else {
                Smin[sm++] = min1;
                mdlRelease(mdl,CID_CELL,pkdn);
                pkdn = pkdu;
                if (min2 >= pq->fDist2) goto NotContained;
            }
        }
        pEnd = pkdn->pUpper;
        for (pj=pkdn->pLower;pj<=pEnd;++pj) {
            p = mdlAquire(mdl,CID_PARTICLE,pj,id);
            //if (smHashPresent(smx,p)) continue;
            if (!p->bSrcActive) {
                mdlRelease(mdl,CID_PARTICLE,p);
                continue;
            }
            v = pkdVel(pkd,p);
            dr0 = (r[0] - p->r[0]) * scale[0];
            dr1 = (r[1] - p->r[1]) * scale[1];
            dr2 = (r[2] - p->r[2]) * scale[2];
            dr3 = (r[3] -    v[0]) * scale[3];
            dr4 = (r[4] -    v[1]) * scale[4];
            dr5 = (r[5] -    v[2]) * scale[5];
            fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
            if (fDist2 < pq->fDist2) {
                if (pq->iPid == idSelf) {
                    smx->ea[pq->iIndex].bInactive = 0;
                }
                else {
                    //smHashDel(smx,pq->pPart);
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
                //smHashAdd(smx,p);
                PQ_REPLACE(pq);
            }
            else mdlRelease(mdl,CID_PARTICLE,p);
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
        mdlRelease(mdl,CID_CELL,pkdn);
        pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
        /*
        ** Intersection Test. (ball-test)
        */
        if (sm) {
            min2 = Smin[--sm];
            }
        else {
            bnd = pkdNodeBnd(pkd, pkdn);
            MINDIST6(bnd,r,min2,scale);
            }
        if (min2 >= pq->fDist2) {
            kdn = pkdTreeNode(pkd,iCell = pkdn->iParent);
            mdlRelease(mdl,CID_CELL,pkdn);
            pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
            goto NoIntersect;
            }
        S[++sp] = iCell;
    }
}


PQ *pqSearchPsd(SMX smx,PQ *pq,FLOAT r[6],FLOAT scale[6], int bReplica,int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    int idSelf = smx->pkd->idSelf;
    FLOAT *Smin = smx->SminT;
    int *S = smx->ST;
    FLOAT dMin,min1,min2;
    int j,iCell,id,iParent;
    int sp = 0;
    int sm = 0;
    pBND bnd;

    *pbDone = 0;
    if (bReplica) kdn = pkdTopNode(pkd,iCell = ROOT);
    else {
        kdn = pkdTopNode(pkd,iCell = pkd->iTopRoot);
        assert(kdn->pLower == idSelf);
    }
    if (iCell != ROOT) S[sp] = kdn->iParent;
    else S[sp] = iCell;

    while (1) {
        /*
        ** Descend to bucket via the closest cell at each level.
        */
        while (kdn->iLower) {
            kdn = pkdTopNode(pkd,iCell = kdn->iLower); 
            bnd = pkdNodeBnd(pkd, kdn);
            MINDIST6(bnd,r,min1,scale);
            kdn = pkdTopNode(pkd,++iCell);             
            bnd = pkdNodeBnd(pkd, kdn);
            MINDIST6(bnd,r,min2,scale);
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
        id = kdn->pLower;        /* this is the thread id in LTT */
        if (id == pkd->idSelf) {
            pq = pqSearchLocalPsd(smx,pq,r,scale,pbDone);
            if (*pbDone) return(pq);
        }
        else {
            pq = pqSearchRemotePsd(smx,pq,id,r,scale);
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
                bnd = pkdNodeBnd(pkd, kdn);
                for (j=0;j<6;++j) {
                    dMin = (bnd.fMax[j] - fabs(bnd.fCenter[j] - r[j])) * scale[j];
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
            bnd = pkdNodeBnd(pkd, kdn);
            MINDIST6(bnd,r,min2,scale);
        }
        if (min2 >= pq->fDist2) {
            kdn = pkdTopNode(pkd,iCell = kdn->iParent);
            goto NoIntersect;
        }
        S[++sp] = iCell;
    }
}


#if 1
void psdSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;
    int N;
    printf("psdSmooth (file)\n");

    FILE *fp = fopen("test.den", "r"); assert(fp != NULL);
    fscanf(fp, "%i\n", &N);
    assert(N == pkd->nLocal);
    for (pi=0;pi<pkd->nLocal;++pi) {
        p = pkdParticle(pkd,pi);
        fscanf(fp, "%f\n", &p->fDensity);
    }
    fclose(fp);

    fp = fopen("test.ball", "r"); assert(fp != NULL);
    fscanf(fp, "%i\n", &N);
    assert(N == pkd->nLocal);
    for (pi=0;pi<pkd->nLocal;++pi) {
        p = pkdParticle(pkd,pi);
        fscanf(fp, "%f\n", &p->fBall);
    }
    fclose(fp);
    printf("psdSmooth (file) finished\n");
}

#else
void psdSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    PQ *pq = smx->pq;
    FLOAT r[6],r2[6],fBall;
    FLOAT rLast[6];
    FLOAT scale[6];
    int iStart[6],iEnd[6];
    int pi,i,j,bDone;
    int ix,iy,iz;
    int ivx,ivy,ivz;
    double *v;

    printf("psdSmooth\n");

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (pi=0;pi<pkd->nLocal;++pi) {
        p = pkdParticle(pkd,pi);
        smx->ea[pi].bInactive = (p->bSrcActive)?0:1;
    }
    smx->ea[pkd->nLocal].bInactive = 0;  /* initialize for Sentinel, but this is not really needed */
#if 0
    for (pi=0;pi<pkd->nLocal;++pi) {
        //printf("pi=%i\n", pi);
        p = pkdParticle(pkd,pi);
        if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
        psdDensity(p,smx->nSmooth,NULL,smf,scale);
    }
#else
    /*
    ** Initialize the priority queue first.
    */
    for (i=0;i<smx->nSmooth;++i) {
        smx->pq[i].pPart = &smx->pSentinel;
        smx->pq[i].iIndex = pkd->nLocal;
        smx->pq[i].iPid = pkd->idSelf;
        smx->pq[i].dr[0] = smx->pSentinel.r[0];
        smx->pq[i].dr[1] = smx->pSentinel.r[1];
        smx->pq[i].dr[2] = smx->pSentinel.r[2];
        smx->pq[i].dr[3] = HUGE_VAL;
        smx->pq[i].dr[4] = HUGE_VAL;
        smx->pq[i].dr[5] = HUGE_VAL;
        smx->pq[i].fDist2 = HUGE_VAL;
        //smx->pq[i].fDist2 = 0;
        //for (j=0; j < 6; j++) smx->pq[i].fDist2 += pow(smx->pq[i].dr[j],2);
    }
    for (j=0;j<6;++j) rLast[j] = 0.0;

    for (pi=0;pi<pkd->nLocal;++pi) {
        //fprintf(stderr, "pi=%i\n", pi);
        p = pkdParticle(pkd,pi);
        if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
        v = pkdVel(pkd, p);

        for (j=0;j<3;++j) r[j] = p->r[j];
        for (j=3;j<6;++j) r[j] =    v[j-3];

        for (j=0;j<6;++j) 
        {
            if (pkdPsMetric(pkd, p)->scale[j] == 0)
                scale[j] = 1;
            else
                scale[j] = 1.0 / pkdPsMetric(pkd, p)->scale[j];
            assert(!isinf(scale[j]));
            //fprintf(stderr, "%g ", scale[j]);
        }

        /*
        ** Correct distances and rebuild priority queue.
        */
        for (i=0;i<smx->nSmooth;++i) {
            smx->pq[i].fDist2 = 0;
            for (j=0; j < 6; j++) {
                smx->pq[i].dr[j]  += r[j] - rLast[j];
                smx->pq[i].dr[j]  *= scale[j];
                smx->pq[i].fDist2 += pow(smx->pq[i].dr[j],2);
            }
            //assert(smx->pq[i].fDist2 > 0);
        }
        for (j=0;j<6;++j) rLast[j] = r[j];
        PQ_BUILD(smx->pq,smx->nSmooth,pq);


        pq = pqSearchPsd(smx,pq,r,scale,0,&bDone);

        //printf("%i] %g %g %g %g %g %g\n", pi, scale[0], scale[1], scale[2], scale[3], scale[4],scale[5]);

        /*
        ** Search in replica boxes if it is required.
        */
        if (!bDone && smx->bPeriodic) {
            fBall = sqrt(pq->fDist2);
            /* scaling missing here... */
            assert(0);
            for (j=0;j<6;++j) {
                iStart[j] = floor((r[j] - fBall)/pkd->fPeriod[j] + 0.5);
                iEnd[j]   = floor((r[j] + fBall)/pkd->fPeriod[j] + 0.5);
            }
            void pqSearchRecur(int i, int not_center) /* inline recursive function */
            {
                int j;
                if (i == 6) {
                    if (not_center) pq = pqSearchPsd(smx,pq,r2,scale,1,&bDone);
                    } 
                else {

                    for (j=iStart[i]; j <= iEnd[i]; j++) {
                        r2[i] = r[i] - j*pkd->fPeriod[i];
                        pqSearchRecur(i+1, not_center || j);
                        }
                    }
            }
            pqSearchRecur(0,0);
        }
        p->fBall = sqrt(pq->fDist2);
        /*
        ** Apply smooth function to the neighbor list.
        */
        psdDensity(p,smx->nSmooth,smx->pq,smf,scale);
        /*
        ** Call mdlCacheCheck to make sure we are making progress!
        */
        mdlCacheCheck(pkd->mdl);

        for (i=0;i<smx->nSmooth;++i)
            for (j=0; j < 6; j++) 
                smx->pq[i].dr[j] /= scale[j];
    }
#endif
    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<smx->nSmooth;++i) {
        if (smx->pq[i].iPid == pkd->idSelf) {
            smx->ea[smx->pq[i].iIndex].bInactive = 0;
        }
        else {
            //smHashDel(smx,smx->pq[i].pPart);
            mdlRelease(pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
        }
    }

    printf("psdSmooth Finished\n");
}
#endif

void psdGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[6],FLOAT scale[6]) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    FLOAT min2,fDist2;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,nCnt,pEnd;
    int idSelf = pkd->idSelf;
    double *v;
    pBND bnd;

    nCnt = smx->nnListSize;
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
        bnd = pkdNodeBnd(pkd, kdn);
        MINDIST6(bnd,r,min2,scale);
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
                p = pkdParticle(pkd,pj);
                if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
                v = pkdVel(pkd, p);
                dr0 = (r[0] - p->r[0]) * scale[0];
                dr1 = (r[1] - p->r[1]) * scale[1];
                dr2 = (r[2] - p->r[2]) * scale[2];
                dr3 = (r[3] -    v[0]) * scale[3];
                dr4 = (r[4] -    v[1]) * scale[4];
                dr5 = (r[5] -    v[2]) * scale[5];
                fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
                if (fDist2 <= fBall2) {
                    if (nCnt >= smx->nnListMax) {
                        smx->nnListMax += NNLIST_INCREMENT;
                        smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
                        assert(smx->nnList != NULL);
                    }
                    smx->nnList[nCnt].fDist2 = fDist2;
                    smx->nnList[nCnt].dr[0] = dr0;
                    smx->nnList[nCnt].dr[1] = dr1;
                    smx->nnList[nCnt].dr[2] = dr2;
                    smx->nnList[nCnt].dr[3] = dr3;
                    smx->nnList[nCnt].dr[4] = dr4;
                    smx->nnList[nCnt].dr[5] = dr5;
                    smx->nnList[nCnt].pPart = p;
                    smx->nnList[nCnt].iIndex = pj;
                    smx->nnList[nCnt].iPid = idSelf;
                    ++nCnt;
                }
            }
        }
    NoIntersect:
        if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp]);
        else break;
    }
    smx->nnListSize = nCnt;
}


void psdGatherRemote(SMX smx,FLOAT fBall2,FLOAT r[6],FLOAT scale[6], int id) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *pkdn;
    PARTICLE *pp;
    FLOAT min2,fDist2;
    FLOAT dr0,dr1,dr2,dr3,dr4,dr5;
    int *S = smx->S;
    int sp = 0;
    int pj,nCnt,pEnd;
    int iCell;
    double *v;
    pBND bnd;

    assert(id != smx->pkd->idSelf);
    nCnt = smx->nnListSize;
    iCell = ROOT;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    while (1) {
        bnd = pkdNodeBnd(pkd,pkdn);
        MINDIST6(bnd,r,min2,scale);
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
                v = pkdVel(pkd, pp);
                dr0 = (r[0] - pp->r[0]) * scale[0];
                dr1 = (r[1] - pp->r[1]) * scale[1];
                dr2 = (r[2] - pp->r[2]) * scale[2];
                dr3 = (r[3] -     v[0]) * scale[3];
                dr4 = (r[4] -     v[1]) * scale[4];
                dr5 = (r[5] -     v[2]) * scale[5];
                fDist2 = dr0*dr0 + dr1*dr1 + dr2*dr2 + dr3*dr3 + dr4*dr4 + dr5*dr5;
                if (fDist2 <= fBall2) {
                    if (nCnt >= smx->nnListMax) {
                        smx->nnListMax += NNLIST_INCREMENT;
                        smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
                        assert(smx->nnList != NULL);
                    }
                    smx->nnList[nCnt].fDist2 = fDist2;
                    smx->nnList[nCnt].dr[0] = dr0;
                    smx->nnList[nCnt].dr[1] = dr1;
                    smx->nnList[nCnt].dr[2] = dr2;
                    smx->nnList[nCnt].dr[3] = dr3;
                    smx->nnList[nCnt].dr[4] = dr4;
                    smx->nnList[nCnt].dr[5] = dr5;
                    smx->nnList[nCnt].pPart = pp;
                    smx->nnList[nCnt].iIndex = pj;
                    smx->nnList[nCnt].iPid = id;
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

void psdGather(SMX smx,FLOAT fBall2,FLOAT r[6], FLOAT scale[6]) {
    KDN *kdn;
    PKD pkd = smx->pkd;
    int *S = smx->ST;
    FLOAT min2;
    int iCell,id;
    int sp = 0;
    pBND bnd;

    kdn = pkdTopNode(pkd,iCell = ROOT);
    while (1) {
        bnd = pkdNodeBnd(pkd, kdn);
        MINDIST6(bnd,r,min2,scale);
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
                psdGatherRemote(smx,fBall2,r,scale,id);
            }
            else {
                psdGatherLocal(smx,fBall2,r,scale);
            }
        }
    NoIntersect:
        if (sp) kdn = pkdTopNode(pkd,iCell = S[--sp]);
        else break;
    }
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

static int CmpRMs(void *ctx,const void *v1,const void *v2) {
    FOFRM *rm1 = (FOFRM *)v1;
    FOFRM *rm2 = (FOFRM *)v2;
    if (rm1->iPid != rm2->iPid) return (rm1->iPid - rm2->iPid);
    else return (rm1->iIndex - rm2->iIndex);
}

static int CmpGroups(const void *v1,const void *v2) {
    FOFGD *g1 = (FOFGD *)v1;
    FOFGD *g2 = (FOFGD *)v2;
    return g2->nTotal - g1->nTotal;
}

typedef struct
{
    int pid;
    int gid;
} trial_group_t;

void psdFof(SMX smx,SMF *smf) {
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
    int iMaxGroups, iGroup;
    int iStart[6],iEnd[6];
    int ix,iy,iz;
    int idSelf = pkd->idSelf;

    FLOAT r[6],l[6],relpos[6],lx,ly,lz,fBall,rho;
    FLOAT scale[6];
    FLOAT fMass;
    int nTree,tmp;

#define TEMP_S_INCREASE 100
    int *C;		/* this is the stack */
    int *T;
    NEW_STACK(C, TEMP_S_INCREASE);
    NEW_STACK(T, TEMP_S_INCREASE);

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

    tmp = pkd->nDark+pkd->nGas+pkd->nStar;

    nTree = pkdTreeNode(pkd,ROOT)->pUpper + 1;
    iMaxGroups = nTree+1;

    trial_group_t *tg = malloc(nTree*sizeof(*tg));

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

    for (pn=0;pn<nTree;pn++) 
    {
        p = pkdParticle(pkd,pn);
        fMass = pkdMass(pkd,p);
        pGroup = pkdInt32(p,pkd->oGroup);
        *pGroup = 0; 

        //p->fBall = 1.0;
    }

    /* Have to restart particle cache, since we will need the updated p->fBall now */
    mdlFinishCache(mdl,CID_PARTICLE);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),nTree);

    int trial_group = iGroup;

    int active_groups = 0;

    /* Starting FOF search now... */
    for (pn=0;pn<nTree;pn++) 
    {
        p = pkdParticle(pkd,pn);
        if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;

        pGroup = pkdInt32(p,pkd->oGroup);
        if (*pGroup ) continue;


        if (trial_group == iGroup)
        {
            iGroup++;
            active_groups++;
            printf("New Group %i / %i\n", iGroup, active_groups);
        }

        assert(iGroup < iMaxGroups);
        protoGroup[iGroup].nMembers = 0;
        protoGroup[iGroup].iId = iGroup;
        protoGroup[iGroup].treeRemoteMembers = NULL;
        nRmCnt = 0;

        trial_group = iGroup;

        pi = pn;

        assert(STACK_EMPTY(C));
        assert(STACK_EMPTY(T));
        while (pi != -1) 
        {
            fprintf(stderr, "%i\n", pi);
            EXTEND_STACK(C);
            EXTEND_STACK(T);

            PUSH(C,pi);
            PUSH(T,pi);

            p = pkdParticle(pkd,pi);
            pGroup = pkdInt32(p,pkd->oGroup);
            *pGroup = trial_group;


            /*
            ** Do a Ball Gather at the radius p->fBall
            */

            v = pkdVel(pkd, p);

            for (j=0;j<6;++j) 
            {
                if (pkdPsMetric(pkd, p)->scale[j] == 0)
                    scale[j] = 1;
                else
                    scale[j] = 1.0 / pkdPsMetric(pkd, p)->scale[j];
                assert(!isinf(scale[j]));
            }

            for (j=0;j<3;++j) r[j] = p->r[j];
            for (j=3;j<6;++j) r[j] =    v[j-3];

            float fDensity = p->fDensity; 

            smx->nnListSize = 0;
            fBall = sqrt(p->fBall);

            if (smx->bPeriodic) 
            {
                assert(0);
                for (j=0;j<3;++j) 
                {
                    iStart[j] = floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5);
                    iEnd[j] = floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5);
                }
                for (ix=iStart[0];ix<=iEnd[0];++ix) 
                {
                    r[0] = p->r[0] - ix*pkd->fPeriod[0];
                    for (iy=iStart[1];iy<=iEnd[1];++iy) 
                    {
                        r[1] = p->r[1] - iy*pkd->fPeriod[1];
                        for (iz=iStart[2];iz<=iEnd[2];++iz) 
                        {
                            r[2] = p->r[2] - iz*pkd->fPeriod[2];
                            psdGather(smx,p->fBall,r,scale);
                        }
                    }
                }
            }
            else 
            {
                psdGather(smx,p->fBall,r,scale);
            }
            nCnt = smx->nnListSize;

            int nn_compar(const void *a0, const void *b0)
            {
                NN *a = (NN *)a0;
                NN *b = (NN *)b0;
                if (a->fDist2 < b->fDist2) return -1;
                if (a->fDist2 > b->fDist2) return +1;
                return 0;
            }
            qsort(smx->nnList, nCnt, sizeof(*smx->nnList), nn_compar);

            int pnn_max = -1;
            assert(smx->nnList[0].pPart->r[0] == p->r[0]);
            assert(smx->nnList[0].pPart->r[1] == p->r[1]);
            assert(smx->nnList[0].pPart->r[2] == p->r[2]);

            //fprintf(stderr, "%i\n", pi);
            for (pnn=1;pnn<nCnt;++pnn ) 
            {
                //if (pi == smx->nnList[pnn].iIndex) continue;
                if (*pkdInt32(smx->nnList[pnn].pPart,pkd->oGroup) == trial_group) continue;
                //if ((smx->nnList[pnn].pPart->fDensity - fDensity)/fDensity > 0.8) continue;

                if (smx->nnList[pnn].pPart->fDensity > fDensity)
                //if ((smx->nnList[pnn].pPart->fDensity - fDensity)/fDensity > -0.5)
                //if ((int)(10*log10(smx->nnList[pnn].pPart->fDensity)) >= (int)(10*log10(fDensity)))
                {
                    pnn_max = pnn;
                    break;
                }

                //fprintf(stderr, "%i  %f  %f  %f\n", pi, smx->nnList[pnn].pPart->fDensity, fDensity, (smx->nnList[pnn].pPart->fDensity - fDensity)/fDensity);
            }

            if (pnn_max == -1) /* Must be a peak */
            {
                pi = -1;

                if (p->fDensity < 1e9)
                {
                    if (!STACK_EMPTY(T))
                    {
                        POP(T);
                        if (!STACK_EMPTY(T))
                            pi = POP(T);
                    }
                }
                //fprintf(stderr, "peak\n");
            }
            else 
            {
                pPartGroup = pkdInt32(smx->nnList[pnn_max].pPart,pkd->oGroup);
                if (*pPartGroup != 0) /* Found another group; join it */
                {
                    active_groups--;
                    trial_group = *pPartGroup;
                    pi = -1;
                }
                else 
                {
                    pi = smx->nnList[pnn_max].iIndex;
                }
            }
        }

        CLEAR_STACK(T);

        while (!STACK_EMPTY(C))
        {
            pi = POP(C);
            p = pkdParticle(pkd,pi);
            pGroup = pkdInt32(p,pkd->oGroup);
            *pGroup = trial_group;

            tg[pi].pid = pi;
            tg[pi].gid = *pGroup;
        }

        if ( nRmCnt > pkd->nMaxRm ) pkd->nMaxRm = nRmCnt;
    }

    


    /*
    ** Now we can already reject groups which are too small, if they are entirely local
    */
    for (i=0; i<nTree ; i++) {
        if (tg[i].gid >=0 && tg[i].gid < iMaxGroups)
            protoGroup[tg[i].gid].nMembers++;
        else
            printf("ERROR: idSelf=%i , p->pGroup=%i too large. iMaxGroups=%i \n",pkd->idSelf,tg[i].gid,iMaxGroups);
    }

#if 0
    /*
    ** Create a remapping and give unique local Ids !
    */
    iMaxGroups = iGroup;
    iGroup= 1 + idSelf;
    pkd->nGroups = 0;
    protoGroup[0].iId = tmp;
    for (i=1;i<=iMaxGroups;i++) {
        protoGroup[i].iId = iGroup;
        if (protoGroup[i].nMembers < smf->nMinMembers) {
        //if (protoGroup[i].nMembers < smf->nMinMembers && protoGroup[i].treeRemoteMembers == NULL) {
            protoGroup[i].iId = tmp;
        }
        else {
            iGroup += pkd->nThreads;
            ++pkd->nGroups;
        }
    }
#endif
    /*
    ** Update the particle groups ids.
    */
    for (i=0;i<nTree;i++) 
    {
        if (protoGroup[tg[i].gid].nMembers < smf->nMinMembers)
            tg[i].gid = 0;
    }

    int tg_compar(const void *a0, const void *b0)
    {
        trial_group_t *a = (trial_group_t *)a0;
        trial_group_t *b = (trial_group_t *)b0;
        return a->gid - b->gid;
    }
    qsort(tg, nTree, sizeof(*tg), tg_compar);

    trial_group = 0;
    int last_group = 0;
    for (i=0; i < nTree; i++)
    {
        if (tg[i].gid != last_group)
        {
            last_group = tg[i].gid;
            trial_group++;
        }

        *pkdInt32(pkdParticle(pkd,tg[i].pid),pkd->oGroup) = trial_group;

    }


#if 0
    /*
    ** Allocate the remote members array
    */
    pkd->nRm = nRmListSize;
    pkd->remoteMember = mdlMalloc(mdl,(pkd->nRm+1)*sizeof(FOFRM));
    iRmIndex = 0;
    /*
    ** Allocate memory for group data
    */
#if 1
    pkd->groupData = (FOFGD *) malloc((1+pkd->nGroups)*sizeof(FOFGD));
    assert(pkd->groupData != NULL);
    k=1;
    for (i=0;i<pkd->nGroups;i++) {
        while (protoGroup[k].iId == tmp) k++;
        pkd->groupData[i].iGlobalId = protoGroup[k].iId;
        pkd->groupData[i].iLocalId = protoGroup[k].iId;
        pkd->groupData[i].nLocal = protoGroup[k].nMembers;

        pkd->groupData[i].nTotal = pkd->groupData[i].nLocal;

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
    //assert( iRmIndex == nRmListSize );
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

    /*
    ** Master orders groupData
    */
    qsort(pkd->groupData,pkd->nGroups,sizeof(FOFGD), CmpGroups);
#endif
#endif
}

#endif /* USE_PSD */
