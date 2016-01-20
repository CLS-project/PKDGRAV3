#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <assert.h>
#include "pkd.h"
#include "moments.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif


uint32_t pkdDistribTopTree(PKD pkd, uint32_t uRoot, uint32_t nTop, KDN *pTop) {
    int i, iTop;
    KDN *pLocalRoot = pkdTreeNode(pkd,uRoot);

    iTop = pkd->iTopTree[uRoot] = pkdTreeAllocNodes(pkd, nTop);
    for(i=0; i<nTop; ++i) {
	KDN *pNode = pkdNode(pkd,pTop,i);
	KDN *pLocal = pkdTreeNode(pkd,iTop+i);
	pkdCopyNode(pkd,pLocal,pNode);
	assert(pLocal->bTopTree);
	if (!pLocal->bRemote) { /* Fixup the links if necessary */
	    pLocal->iLower += iTop + i;
	    pLocal->pUpper += iTop + i;
	    }
	else if (pLocal->pLower == pkd->idSelf) {
	    pLocal->pLower = pLocalRoot->pLower;
	    pLocal->pUpper = pLocalRoot->pUpper;
	    pLocal->iLower = pLocalRoot->iLower;
	    pLocal->bRemote = 0;
	    pLocal->bTopTree = 0;
	    }
	}
    return iTop;
    }

static KDN *InitializeRootCommon(PKD pkd,uint32_t uRoot) {
    KDN *pRoot = pkdTreeNode(pkd,uRoot);
//    BND *bnd;
//    int j;

    pRoot->bTopTree = 0;
    pRoot->bRemote = 0;
    *pkdNodeBnd(pkd, pRoot) = pkd->bnd;
//    for (j=0;j<3;++j) {
//	bnd->fCenter[j] = pkd->bnd.fCenter[j];
//	bnd->fMax[j] = pkd->bnd.fMax[j];
//	}

    return pRoot;
    }

/*
** Creates a single root at ROOT with only marked particles
*/
void pkdTreeInitMarked(PKD pkd) {
    KDN *pRoot   = InitializeRootCommon(pkd,ROOT);
    PARTICLE *pMarked, *pNot = NULL;
    local_t iLast;
    int i;

    if (mdlCacheStatus(pkd->mdl,CID_CELL)) mdlFinishCache(pkd->mdl,CID_CELL);

    pMarked = pkdParticle(pkd,i=0);
    pNot = pkdParticle(pkd,iLast=pkd->nLocal-1);
    PARTITION(pMarked<pNot,pMarked<=pNot,
	pMarked=pkdParticle(pkd,++i),pNot=pkdParticle(pkd,--iLast),
	pkdSwapParticle(pkd,pMarked,pNot),
	pMarked->bMarked!=0,pNot->bMarked==0);

    pRoot->pLower = 0;
    pRoot->pUpper = iLast;
    }

void pkdDumpTrees(PKD pkd,int bOnlyVA, uint8_t uRungDD) {
    KDN *pRoot    = InitializeRootCommon(pkd,ROOT);

    if (pkd->nNodes == 0) {
	assert(bOnlyVA==0);
	pkd->nNodes = NRESERVED_NODES;
	pkdTreeNode(pkd,ROOT)->iLower = NRESERVED_NODES;
	pkdTreeNode(pkd,FIXROOT)->iLower = NRESERVED_NODES;
	}

    /* We can always close the "active" caches */
    if (mdlCacheStatus(pkd->mdl,CID_CELL)) mdlFinishCache(pkd->mdl,CID_CELL);
// done elsewhere   if (mdlCacheStatus(pkd->mdl,CID_PARTICLE)) mdlFinishCache(pkd->mdl,CID_PARTICLE);

    /* Full Normal tree build. Only ROOT will be used */
    if (uRungDD == IRUNGMAX) {
	pRoot->pLower = 0;
	pRoot->pUpper = pkd->nLocal - 1;
	pkd->nNodes = NRESERVED_NODES;
	pRoot->iLower = 0;
	}
    /* Just rebuilding (active) ROOT. Truncate it. pLower and pUpper are still valid. */
    else if (bOnlyVA) {
	if (pRoot->iLower) {
	    assert(pRoot->iLower >= NRESERVED_NODES);
	    pkd->nNodes = pRoot->iLower; /* This effectively truncates the nodes used by this tree */
	    pRoot->iLower = 0;
	    }
	}

    /* Need to build two trees which is more complicated. */
    else {
	KDN *pRootFix = InitializeRootCommon(pkd,FIXROOT);
	PARTICLE *p, *pVA = NULL;
	local_t nVeryActive = 0;
	local_t iLast;
	int i;

	pkd->fixbnd = *pkdNodeBnd(pkd, pRoot);

#ifndef SINGLE_CACHES
	/* Here we also have to close the "fixed" caches */
	if (mdlCacheStatus(pkd->mdl,CID_CELL2))     mdlFinishCache(pkd->mdl,CID_CELL2);
	if (mdlCacheStatus(pkd->mdl,CID_PARTICLE2)) mdlFinishCache(pkd->mdl,CID_PARTICLE2);
#endif
	p = pkdParticle(pkd,i=0);
	pVA = pkdParticle(pkd,iLast=pkd->nLocal-1);
	PARTITION(p<pVA,p<=pVA,
	    p=pkdParticle(pkd,++i),pVA=pkdParticle(pkd,--iLast),
	    pkdSwapParticle(pkd,p,pVA),
	    p->uRung <= uRungDD,pVA->uRung > uRungDD);

	pkd->nNodes = NRESERVED_NODES;
	pRoot->iLower = 0;
	pRootFix->iLower = 0;

	pRootFix->pLower = 0;
	pRootFix->pUpper = iLast;
	pRoot->pLower = iLast + 1;
	pRoot->pUpper = pkd->nLocal - 1;

//	for(i=pRootFix->pLower; i<=pRootFix->pUpper; ++i) assert(pkdParticle(pkd,i)->uRung<=uRungDD);
//	for(i=pRoot->pLower; i<=pRoot->pUpper; ++i) assert(pkdParticle(pkd,i)->uRung>uRungDD);
	}
    }

#define MIN_SRATIO    0.05

/*
** Partition the particles between pLower and pUpper (inclusive)
** around "Split" in dimension "d". Return the index of the split.
*/
static int PartPart(PKD pkd,int pLower,int pUpper,int d,pos_t Split) {
	PARTICLE *pi = pkdParticle(pkd,pLower);
	PARTICLE *pj = pkdParticle(pkd,pUpper);
	while (pi <= pj) {
	    if (pkdPosRaw(pkd,pi,d) < Split) pi = (PARTICLE *)(((char *)pi) + pkdParticleSize(pkd));
	    else break;
	    }
	while (pi <= pj) {
	    if (Split < pkdPosRaw(pkd,pj,d)) pj = (PARTICLE *)(((char *)pj) - pkdParticleSize(pkd));
	    else break;
	    }
	if (pi < pj) {
	    pkdSwapParticle(pkd,pi,pj);
	    while (1) {
		while (pkdPosRaw(pkd,(pi = (PARTICLE *)(((char *)pi) + pkdParticleSize(pkd))),d) < Split);
		while (Split < pkdPosRaw(pkd,(pj = (PARTICLE *)(((char *)pj) - pkdParticleSize(pkd))),d));
		if (pi < pj) {
		    pkdSwapParticle(pkd,pi,pj);
		    }
		else break;
		}
	    }
	return ((char *)pi - (char *)pkdParticleBase(pkd)) / pkdParticleSize(pkd);
    }


/*
** M is the bucket size.
** This function assumes that the root node is correctly set up (particularly the bounds).
*/
#define TEMP_S_INCREASE 100
void BuildTemp(PKD pkd,int iNode,int M,double dMaxMax) {
    PARTICLE *pi, *pj;
    KDN *pNode = pkdTreeNode(pkd,iNode);
    BND *bnd,*lbnd,*rbnd;
    KDN *pLeft, *pRight;
    double lrMax;
    pos_t Split; /* Integer or double */
    int *S;		/* this is the stack */
    int s,ns;
    int iLeft,iRight;
    int d;
    int i,j;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;

    pNode->iDepth = 0;

    // Single bucket? We are done.
    if (pNode->pUpper - pNode->pLower + 1 <= M) return;

    /*
    ** Allocate stack!
    */
    ns = TEMP_S_INCREASE;
    s = 0;
    S = CAST(int*,malloc(ns*sizeof(int)));
    assert(S != NULL);

    bnd = pkdNodeBnd(pkd,pNode);
    assert( bnd->fMax[0] > 0.0 ||
	    bnd->fMax[1] > 0.0 ||
	    bnd->fMax[2] > 0.0 );
    while (1) {
	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/
	if (bnd->fMax[0] < bnd->fMax[1]) {
	    if (bnd->fMax[1] < bnd->fMax[2]) d = 2;
	    else d = 1;
	    }
	else if (bnd->fMax[0] < bnd->fMax[2]) d = 2;
	else d = 0;
	Split = pkdDblToPos(pkd,bnd->fCenter[d]);
	/*
	** Now start the partitioning of the particles about
	** fSplit on dimension given by d.
	*/
	i = PartPart(pkd,pNode->pLower,pNode->pUpper,d,Split);
	nl = i - pNode->pLower;
	nr = pNode->pUpper - i + 1;
	if (nl > 0 && nr > 0) {
	    /*
	    ** Allocate 2 new tree nodes making sure that we have
	    ** allocated enough storage.
	    */
	    pkdTreeAllocNodePair(pkd,&iLeft,&iRight);
	    pLeft = pkdTreeNode(pkd,iLeft);
	    pLeft->bTopTree = 0;
	    pLeft->bRemote = 0;
	    pLeft->pLower = pNode->pLower;
	    pLeft->pUpper = i-1;
	    pLeft->iDepth = pNode->iDepth+1;
	    pRight = pkdTreeNode(pkd,iRight);
	    assert(iRight & 1);
	    pRight->bTopTree = 0;
	    pRight->bRemote = 0;
	    pRight->pLower = i;
	    pRight->pUpper = pNode->pUpper;
	    pRight->iDepth = pNode->iDepth+1;
	    pNode->iLower = iLeft;

            lbnd = pkdNodeBnd(pkd, pLeft);
            rbnd = pkdNodeBnd(pkd, pRight);

	    /*
	    ** Now deal with the bounds.
	    */
	    lrMax = 0.0;
	    for (j=0;j<3;++j) {
		if (j == d) {
		    rbnd->fMax[j] = lbnd->fMax[j] = 0.5*bnd->fMax[j];
		    lbnd->fCenter[j] = bnd->fCenter[j] - lbnd->fMax[j];
		    rbnd->fCenter[j] = bnd->fCenter[j] + rbnd->fMax[j];
		    }
		else {
		    lbnd->fCenter[j] = bnd->fCenter[j];
		    lbnd->fMax[j] = bnd->fMax[j];
		    rbnd->fCenter[j] = bnd->fCenter[j];
		    rbnd->fMax[j] = bnd->fMax[j];
		    }
		if (rbnd->fMax[j] > lrMax) lrMax = rbnd->fMax[j];
		}
	    /*
	    ** Now figure out which subfile to process next.
	    */
	    if (lrMax >dMaxMax) {
		lc = (nl > 0); /* this condition means the left child is not a bucket */
		rc = (nr > 0);
		}
	    else {
		lc = (nl > M); /* this condition means the left child is not a bucket */
		rc = (nr > M);
		}
	    if (rc && lc) {
		/* Allocate more stack if required */
		if ( s+1 >= ns ) {
		    assert( s+1 == ns );
		    ns += TEMP_S_INCREASE;
		    S = CAST(int *,realloc(S,ns*sizeof(int)));
		    }
		if (nr > nl) {
		    S[s++] = iRight;	/* push tr */
		    iNode = iLeft;		/* process lower subfile */
		    }
		else {
		    S[s++] = iLeft;	/* push tl */
		    iNode = iRight;		/* process upper subfile */
		    }
		}
	    else if (lc) {
		/*
		** Right must be a bucket in this case!
		*/
		iNode = iLeft;   /* process lower subfile */
		pRight->iLower = 0;
		++nBucket;
		}
	    else if (rc) {
		/*
		** Left must be a bucket in this case!
		*/
		iNode = iRight;   /* process upper subfile */
		pLeft->iLower = 0;
		++nBucket;
		}
	    else {
		/*
		** Both are buckets (we need to pop from the stack to get the next subfile.
		*/
		pLeft->iLower = 0;
		++nBucket;
		pRight->iLower = 0;
		++nBucket;
		if (s) iNode = S[--s];		/* pop tn */
		else break;
		}
	    }
	else {
	    ++pNode->iDepth;
	    /*
	    ** No nodes allocated, Change the bounds if needed!
	    */
	    if (d >= 0 && d < 3) bnd->fMax[d] *= 0.5;
	    lrMax = 0.0;
	    for(j=0; j<3; ++j)
		if (bnd->fMax[j] > lrMax) lrMax = bnd->fMax[j];
	    if (nl > 0) {
		if (d >= 0 && d < 3) bnd->fCenter[d] -= bnd->fMax[d];
		lc = (lrMax>dMaxMax || nl > M); /* this condition means the node is not a bucket */
		if (!lc) {
		    pNode->iLower = 0;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    else {
		if (d >= 0 && d < 3) bnd->fCenter[d] += bnd->fMax[d];
		rc = (lrMax>dMaxMax || nr > M);
		if (!rc) {
		    pNode->iLower = 0;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    }
	pNode = pkdTreeNode(pkd,iNode);
        bnd = pkdNodeBnd(pkd, pNode);
	}
    free(S);
    }


/*
** With more than a single tree, we must be careful to make sure
** that they match or the P-P, P-C and hence checklists will explode.
** This routine takes a template tree (usually the large "fixed" tree)
** and contructs a tree with buckets no larger than the corresponding
** buckets in the template tree. Once we reach a bucket in the template
** tree, then we can either create a bucket, or continue to create the
** new tree in the usual way.
*/
void BuildFromTemplate(PKD pkd,int iNode,int M,int iTemplate) {

    // The tree we are building
    KDN *pNode, *pLeft, *pRight;
    BND *bnd, *lbnd, *rbnd; 
    int iLeft, iRight;

    // The "template" tree
    KDN *pTemp, *ptLeft, *ptRight;
    BND tbnd, *ltbnd, *rtbnd;

    pos_t Split;
    int d;
    int i, j, nr, nl;
    int lg, rg;

    double dMaxMax;
    int ns,s;
    struct buildStack {
	BND tbnd;
	KDN *pTemp;
	int iNode;
	} *S;
    // If we have an empty tree root then we are "done"
    pNode = pkdTreeNode(pkd,iNode);
    pNode->iDepth = 0;
    pNode->iLower = 0;
    if (pNode->pUpper - pNode->pLower + 1 == 0) return;

    // Setup stack and push the two root cells of each tree.
    ns = TEMP_S_INCREASE;
    S = (struct buildStack *)malloc(ns*sizeof(struct buildStack));
    assert(S != NULL);
    s = 1;
    S[0].iNode = iNode;
    S[0].pTemp = pkdTreeNode(pkd,iTemplate);
    S[0].tbnd = pkd->fixbnd;

    do {
	// Pop the next cells to process
	assert(s);
	--s;
	pTemp = S[s].pTemp;
	tbnd = S[s].tbnd;
	iNode = S[s].iNode;
	bnd = pkdNodeBnd(pkd,pNode = pkdTreeNode(pkd,iNode));

	// Follow the template tree to wherever it leads.
	while(pTemp->iLower) {
	    assert(pTemp->iDepth >= pNode->iDepth);
	    ltbnd = pkdNodeBnd(pkd, ptLeft=pkdTreeNode(pkd,pTemp->iLower));
	    rtbnd = pkdNodeBnd(pkd, ptRight=pkdTreeNode(pkd,pTemp->iLower+1));

	    // Calculate the dimension and it's split value based on the algorithm above
	    if (tbnd.fMax[0] < tbnd.fMax[1]) {
		if (tbnd.fMax[1] < tbnd.fMax[2]) d = 2;
		else d = 1;
		}
	    else if (tbnd.fMax[0] < tbnd.fMax[2]) d = 2;
	    else d = 0;
	    Split = pkdDblToPos(pkd,tbnd.fCenter[d]);

	    // We must have done a split along this dimension or something is wacked
	    assert(ltbnd->fCenter[d] != tbnd.fCenter[d]);
	    assert(rtbnd->fCenter[d] != tbnd.fCenter[d]);

	    // Check for ghosting on the template tree
	    lg = (ltbnd->fCenter[d] > tbnd.fCenter[d]);
	    rg = (rtbnd->fCenter[d] < tbnd.fCenter[d]);
	    if (pTemp->iDepth == pNode->iDepth) {
		assert(lg==0 && rg==0); // Neither can be ghosted
		}
	    else {
		assert(lg==0 || rg==0); // Obviously both cannot be ghosted
		}

	    // Partition the particles on either side of the split
	    i = PartPart(pkd,pNode->pLower,pNode->pUpper,d,Split);
	    nl = i - pNode->pLower;
	    nr = pNode->pUpper - i + 1;
	    assert(nl>0 || nr>0);

	    if (nl==0) {
		++pNode->iDepth;
		bnd->fMax[d] = 0.5 * (bnd->fCenter[d] + bnd->fMax[d] - tbnd.fCenter[d]);
		assert(bnd->fMax[d] >0.0);
		bnd->fCenter[d] = tbnd.fCenter[d] + bnd->fMax[d];
		if (rg) break; // Just out and build normally
		else {
		    tbnd.fMax[d] *= 0.5;
		    tbnd.fCenter[d] += tbnd.fMax[d];
		    pTemp = ptRight;
		    continue;
		    }
		}
	    else if (nr==0) {
		++pNode->iDepth;
		bnd->fMax[d] = 0.5 * (tbnd.fCenter[d] - bnd->fCenter[d] + bnd->fMax[d]);
		assert(bnd->fMax[d] >0.0);
		bnd->fCenter[d] = tbnd.fCenter[d] - bnd->fMax[d];
		if (lg) break; // Just out and build normally
		else {
		    tbnd.fMax[d] *= 0.5;
		    tbnd.fCenter[d] -= tbnd.fMax[d];
		    pTemp = ptLeft;
		    continue;
		    }
		}

	    // Good. We have particles on both sides, so we need to split this cell
	    // and setup the appropriate bounds.

	    pkdTreeAllocNodePair(pkd,&iLeft,&iRight);

	    pLeft = pkdTreeNode(pkd,iLeft);
	    pLeft->bTopTree = 0;
	    pLeft->bRemote = 0;
	    pLeft->pLower = pNode->pLower;
	    pLeft->pUpper = i-1;
	    pLeft->iLower = 0;
	    pLeft->iDepth = pNode->iDepth + 1;
	    assert(pLeft->pLower <= pLeft->pUpper);

	    pRight = pkdTreeNode(pkd,iRight);
	    assert(iRight & 1);
	    pRight->bTopTree = 0;
	    pRight->bRemote = 0;
	    pRight->pLower = i;
	    pRight->pUpper = pNode->pUpper;
	    pRight->iLower = 0;
	    pRight->iDepth = pNode->iDepth + 1;
	    assert(pRight->pLower <= pRight->pUpper);

	    pNode->iLower = iLeft;

	    lbnd = pkdNodeBnd(pkd, pLeft);
	    rbnd = pkdNodeBnd(pkd, pRight);

	    iNode = -1;
	    pNode = NULL;

	    /*
	    ** Now deal with the bounds.
	    ** tbnd.fCenter[d] is the "split"
	    */
	    for (j=0;j<3;++j) {
		if (j == d) {
		    rbnd->fMax[d] = 0.5 * (bnd->fCenter[d] + bnd->fMax[d] - tbnd.fCenter[d]);
		    rbnd->fCenter[d] = tbnd.fCenter[d] + rbnd->fMax[d];
		    lbnd->fMax[d] = 0.5 * (tbnd.fCenter[d] - bnd->fCenter[d] + bnd->fMax[d]);
		    lbnd->fCenter[d] = tbnd.fCenter[d] - lbnd->fMax[d];
		    assert(rbnd->fMax[d] >0.0);
		    assert(lbnd->fMax[d] >0.0);
		    }
		else {
		    lbnd->fCenter[j] = bnd->fCenter[j];
		    lbnd->fMax[j] = bnd->fMax[j];
		    rbnd->fCenter[j] = bnd->fCenter[j];
		    rbnd->fMax[j] = bnd->fMax[j];
		    }
		}

	    tbnd.fMax[d] *= 0.5; // Shrink template cell
	    dMaxMax = HUGE_VAL;
	    for (j=0;j<3;++j) if (tbnd.fMax[j]<dMaxMax) dMaxMax = tbnd.fMax[j];
	    assert(dMaxMax>0.0);

	    // Left is a ghost cell or bucket
	    if (lg || ptLeft->iLower==0) {
		BuildTemp(pkd,iLeft,M,dMaxMax);
		tbnd.fCenter[d] += tbnd.fMax[d];
		iNode = iRight;
		pTemp = ptRight;
		}
	    else if (rg || ptRight->iLower==0) {
		BuildTemp(pkd,iRight,M,dMaxMax);
		tbnd.fCenter[d] -= tbnd.fMax[d];
		iNode = iLeft;
		pTemp = ptLeft;
		}
	    // push right and continue left
	    else {
		if ( s+1 >= ns ) {
		    assert( s+1 == ns );
		    ns += TEMP_S_INCREASE;
		    S = (struct buildStack *)realloc(S,ns*sizeof(struct buildStack));
		    }
		S[s].tbnd = tbnd;
		S[s].tbnd.fCenter[d] += tbnd.fMax[d];
		S[s].pTemp = ptRight;
		S[s].iNode = iRight;
		++s;

		tbnd.fCenter[d] -= tbnd.fMax[d];
		pTemp = ptLeft;
		iNode = iLeft;
		}
	    pNode = pkdTreeNode(pkd,iNode);
	    bnd = pkdNodeBnd(pkd,pNode);
	    }

	dMaxMax = HUGE_VAL;
	for (j=0;j<3;++j) if (tbnd.fMax[j]<dMaxMax) dMaxMax = tbnd.fMax[j];
	assert(dMaxMax>0.0);

	// Bucket in the template tree: Now just build
	BuildTemp(pkd,iNode,M,dMaxMax);
	} while(s);

    free(S);
    }



static vel_t zeroV[3] = {0.0,0.0,0.0};
static float  zeroF[3] = {0.0,0.0,0.0};

void Create(PKD pkd,int iRoot) {
    int iNode = iRoot;
    PARTICLE *p;
    KDN *pkdn,*pkdl,*pkdu;
    FMOMR mom;
    SPHBNDS *bn;
    BND *bnd;
    FLOAT fSoft,x,y,z,ax,ay,az,ft,d2,d2Max,dih2,bmin,b;
    float *a, m, fMass, fBall;
    vel_t *v, vx, vy, vz;
    int pj,d,nDepth,ism;
    const int nMaxStackIncrease = 1;

    /* If the tree is empty, we just create a sensible moment and we are done. */
    pkdn = pkdTreeNode(pkd,iNode);
    if (pkdn->pLower > pkdn->pUpper) {
	pkdn->bMax = 1.0;
	pkdn->uMinRung = MAX_RUNG;
	pkdn->uMaxRung = 0;
	pkdn->bSrcActive = pkdn->bDstActive = 0;
	if (pkd->oNodeMom) momClearFmomr(pkdNodeMom(pkd,pkdn));
	return;
	}

    nDepth = 1;
    while (1) {
	while ((pkdn=pkdTreeNode(pkd,iNode))->iLower) {
	    pkd->S[nDepth-1].iNodeIndex = iNode;
	    iNode = pkdn->iLower;
	    ++nDepth;
	    /*
	    ** Is this the deepest in the tree so far? We might need to have more stack
	    ** elements for the tree walk!
	    ** nMaxStack == nDepth guarantees that there is at least one deeper
	    ** stack entry available than what is needed to walk the tree.
	    */
	    if (nDepth > pkd->nMaxStack) {
		pkd->S = CAST(CSTACK *,realloc(pkd->S,(pkd->nMaxStack+nMaxStackIncrease)*sizeof(CSTACK)));
		assert(pkd->S != NULL);
		for (ism=pkd->nMaxStack;ism<(pkd->nMaxStack+nMaxStackIncrease);++ism) {
		    clInitialize(&pkd->S[ism].cl,&pkd->clFreeList);
		    assert(pkd->S[ism].cl != NULL);
		    }
		pkd->nMaxStack += nMaxStackIncrease;
		}
	    }

	/*
	** Now calculate all bucket quantities!
	** This includes M,CoM,Moments and special
	** bounds and iMaxRung.
	*/
	pkdn = pkdTreeNode(pkd,iNode);
        bnd = pkdNodeBnd(pkd, pkdn);
	/*
	** Before squeezing the bounds, calculate a minimum b value based on the splitting bounds alone.
	** This gives us a better feel for the "size" of a bucket with only a single particle.
	*/
	MINSIDE(bnd->fMax,bmin);
	/*
	** Now shrink wrap the bucket bounds.
	*/
	pj = pkdn->pLower;
	p = pkdParticle(pkd,pj);
	for (d=0;d<3;++d) {
	    ft = pkdPos(pkd,p,d);
	    bnd->fCenter[d] = ft;
	    bnd->fMax[d] = ft;
	    }
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    for (d=0;d<3;++d) {
		ft = pkdPos(pkd,p,d);
		if (ft < bnd->fCenter[d])
		    bnd->fCenter[d] = ft;
		else if (ft > bnd->fMax[d])
		    bnd->fMax[d] = ft;
		}
	    }
	for (d=0;d<3;++d) {
	    ft = bnd->fCenter[d];
	    bnd->fCenter[d] = 0.5*(bnd->fMax[d] + ft);
	    bnd->fMax[d] = 0.5*(bnd->fMax[d] - ft);
	    }
	pj = pkdn->pLower;
	p = pkdParticle(pkd,pj);
	a = pkd->oAcceleration ? pkdAccel(pkd,p) : zeroF;
	m = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	v = pkd->oVelocity ? pkdVel(pkd,p) : zeroV;
	fMass = m;
	dih2 = fSoft;
	x = m*pkdPos(pkd,p,0);
	y = m*pkdPos(pkd,p,1);
	z = m*pkdPos(pkd,p,2);
	vx = m*v[0];
	vy = m*v[1];
	vz = m*v[2];
	ax = m*a[0];
	ay = m*a[1];
	az = m*a[2];
	pkdn->uMinRung = pkdn->uMaxRung = p->uRung;
	pkdn->bDstActive = p->bDstActive;
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    a = pkd->oAcceleration ? pkdAccel(pkd,p) : zeroF;
	    m = pkdMass(pkd,p);
	    fSoft = pkdSoft(pkd,p);
	    v = pkd->oVelocity ? pkdVel(pkd,p) : zeroV;
	    fMass += m;
	    if (fSoft>dih2) dih2=fSoft;
	    x += m*pkdPos(pkd,p,0);
	    y += m*pkdPos(pkd,p,1);
	    z += m*pkdPos(pkd,p,2);
	    vx += m*v[0];
	    vy += m*v[1];
	    vz += m*v[2];
	    ax += m*a[0];
	    ay += m*a[1];
	    az += m*a[2];
	    if ( p->uRung > pkdn->uMaxRung ) pkdn->uMaxRung = p->uRung;
	    if ( p->uRung < pkdn->uMinRung ) pkdn->uMinRung = p->uRung;
	    if ( p->bDstActive ) pkdn->bDstActive = 1;
	    }
	m = 1.0f / fMass;
	if (pkd->param.bCenterOfMassExpand) {
	    pkdn->r[0] = m*x;
	    pkdn->r[1] = m*y;
	    pkdn->r[2] = m*z;
	    }
	else {
	    /*
	    ** For now set it to the center of the bounding box, but later
	    ** we want the tightest bounding sphere here.
	    */
	    for (d=0;d<3;++d) pkdn->r[d] = bnd->fCenter[d];
	    }
	if (pkd->oNodeVelocity) {
	    vel_t *pVel = pkdNodeVel(pkd,pkdn);
	    pVel[0] = m*vx;
	    pVel[1] = m*vy;
	    pVel[2] = m*vz;
	    }
	if (pkd->oNodeAcceleration) {
	    float *pAcc = pkdNodeAccel(pkd,pkdn);
	    pAcc[0] = m*ax;
	    pAcc[1] = m*ay;
	    pAcc[2] = m*az;
	    }
	pkdn->fSoft2 = dih2*dih2;
	d2Max = 0.0;
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    x = pkdPos(pkd,p,0) - pkdn->r[0];
	    y = pkdPos(pkd,p,1) - pkdn->r[1];
	    z = pkdPos(pkd,p,2) - pkdn->r[2];
	    d2 = x*x + y*y + z*z;
	    /*
	    ** Update bounding ball and softened bounding ball.
	    */
	    d2Max = (d2 > d2Max)?d2:d2Max;
	    }
#ifdef USE_MAXSIDE
        MAXSIDE(bnd->fMax,b);
#else
	b = sqrt(d2Max);
#endif
	if (b==0.0) b = 1.0f; /* FIXME: Single particle. Perhaps momMakeFmomr should be more robust. */
        else if (b < bmin) b = bmin;
	pkdn->bMax = b;
	assert(pkdn->bMax>=0);
	/*
	** Now calculate the reduced multipole moment.
	** Note that we use the cell's openening radius as the scaling factor!
	*/
	if (pkd->oNodeMom) {
	    momClearFmomr(pkdNodeMom(pkd,pkdn));
	    for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		p = pkdParticle(pkd,pj);
		x = pkdPos(pkd,p,0) - pkdn->r[0];
		y = pkdPos(pkd,p,1) - pkdn->r[1];
		z = pkdPos(pkd,p,2) - pkdn->r[2];
		m = pkdMass(pkd,p);
		momMakeFmomr(&mom,m,pkdn->bMax,x,y,z);
		momAddFmomr(pkdNodeMom(pkd,pkdn),&mom);
	    }
	}
	/*
	** Calculate bucket fast gas bounds.
	*/
	if (pkd->oNodeSphBounds) {
	    bn = pkdNodeSphBounds(pkd,pkdn);
	    /*
	    ** Default bounds always makes the cell look infinitely far away, regardless from where.
	    */
	    for (d=0;d<3;++d) bn->A.min[d] = HUGE_VAL;
	    for (d=0;d<3;++d) bn->A.max[d] = -HUGE_VAL;
	    for (d=0;d<3;++d) bn->B.min[d] = HUGE_VAL;
	    for (d=0;d<3;++d) bn->B.max[d] = -HUGE_VAL;
	    for (d=0;d<3;++d) bn->BI.min[d] = HUGE_VAL;
	    for (d=0;d<3;++d) bn->BI.max[d] = -HUGE_VAL;
	    for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		p = pkdParticle(pkd,pj);
		fBall = pkdBall(pkd,p);
		if (pkdIsGas(pkd,p)) {
		    /*
		    ** This first ball bound over all gas particles is only used for remote searching.
		    */
		    for (d=0;d<3;++d) bn->B.min[d] = fmin(bn->B.min[d],pkdPos(pkd,p,d) - (1+pkd->param.ddHonHLimit)*fBall);
		    for (d=0;d<3;++d) bn->B.max[d] = fmax(bn->B.max[d],pkdPos(pkd,p,d) + (1+pkd->param.ddHonHLimit)*fBall);
		    if (pkdIsActive(pkd,p)) {
			for (d=0;d<3;++d) bn->A.min[d] = fmin(bn->A.min[d],pkdPos(pkd,p,d));
			for (d=0;d<3;++d) bn->A.max[d] = fmax(bn->A.max[d],pkdPos(pkd,p,d));
		    }
		    else {
			for (d=0;d<3;++d) bn->BI.min[d] = fmin(bn->BI.min[d],pkdPos(pkd,p,d) - (1+pkd->param.ddHonHLimit)*fBall);
			for (d=0;d<3;++d) bn->BI.max[d] = fmax(bn->BI.max[d],pkdPos(pkd,p,d) + (1+pkd->param.ddHonHLimit)*fBall);
		    }
		}
	    }
	    /*
	    ** Note that minimums can always safely be increased and maximums safely decreased in parallel, even on 
	    ** a shared memory machine, without needing locking since these bounds should always simply be seen as being
	    ** a conservative bound on the particles in the algorithms. This is true AS LONG AS a double precision store
	    ** operation is atomic (i.e., that the individual bytes are not written one-by-one). We take advantage of 
	    ** this fact in the fast gas algorithm where we locally reduce the bounds to exclude particles which have 
	    ** already been completed in the direct neighbor search phase.
	    */
	}
	/*
	** Finished with the bucket, move onto the next one,
	** or to the parent.
	*/
	while ((iNode & 1) || iNode==iRoot ) {
	    if ( --nDepth == 0) return; /* exit point!!! */
	    iNode = pkd->S[nDepth-1].iNodeIndex;
	    /*
	    ** Now combine quantities from each of the children (2) of
	    ** this cell to form the quantities for this cell.
	    ** First find the CoM, just like for the bucket.
	    */
	    pkdn = pkdTreeNode(pkd,iNode);
            bnd = pkdNodeBnd(pkd, pkdn);
	    /*
	    ** Before squeezing the bounds, calculate a minimum b value based on the splitting bounds alone.
	    ** This gives us a better feel for the "size" of a bucket with only a single particle.
	    */
	    MINSIDE(bnd->fMax,bmin);
	    pj = pkdn->pLower;
	    pkdl = pkdTreeNode(pkd,pkdn->iLower);
	    pkdu = pkdTreeNode(pkd,pkdn->iLower + 1);
	    pkdCombineCells1(pkd,pkdn,pkdl,pkdu);
	    if (pkdn->pUpper - pj < NMAX_OPENCALC) {
		assert(pj<=pkdn->pUpper);
		p = pkdParticle(pkd,pj);
		x = pkdPos(pkd,p,0) - pkdn->r[0];
		y = pkdPos(pkd,p,1) - pkdn->r[1];
		z = pkdPos(pkd,p,2) - pkdn->r[2];
		d2Max = x*x + y*y + z*z;
		for (++pj;pj<=pkdn->pUpper;++pj) {
		    p = pkdParticle(pkd,pj);
		    x = pkdPos(pkd,p,0) - pkdn->r[0];
		    y = pkdPos(pkd,p,1) - pkdn->r[1];
		    z = pkdPos(pkd,p,2) - pkdn->r[2];
		    d2 = x*x + y*y + z*z;
		    d2Max = (d2 > d2Max)?d2:d2Max;
		    }
		assert(d2Max>0);
		/*
		** Now determine the opening radius for gravity.
		*/
#ifdef USE_MAXSIDE
		MAXSIDE(bnd->fMax,b);
		if (b < bmin) b = bmin;
		if (d2Max>b) b = d2Max;
		pkdn->bMax = b;
#else
		pkdn->bMax = sqrt(d2Max);
		if (pkdn->bMax < bmin) pkdn->bMax = bmin;
#endif
		assert(pkdn->bMax >= 0);
		}
	    else {
	      CALCOPEN(pkdn,bmin);  /* set bMax */
	    }
	    pkdCombineCells2(pkd,pkdn,pkdl,pkdu);
	    }
	++iNode;
	}
    }


void pkdCombineCells1(PKD pkd,KDN *pkdn,KDN *p1,KDN *p2) {
    FLOAT m1,m2,ifMass;
    int j;
    BND *bnd, *p1bnd, *p2bnd;

    bnd = pkdNodeBnd(pkd, pkdn);
    p1bnd = pkdNodeBnd(pkd, p1);
    p2bnd = pkdNodeBnd(pkd, p2);
    BND_COMBINE(bnd,p1bnd,p2bnd);
    if (pkd->oNodeMom) {
	m1 = pkdNodeMom(pkd,p1)->m;
	m2 = pkdNodeMom(pkd,p2)->m;
	ifMass = 1/(m1 + m2);
	/*
	** In the case where a cell has all its particles source inactive mom.m == 0, which is ok, but we
	** still need a reasonable center in order to define opening balls in the tree code.
	*/
	if ( m1==0.0 || m2 == 0.0 ) {
	    ifMass = 1.0;
	    m1 = m2 = 0.5;
	    }
	}
    else {
	ifMass = 1.0;
	m1 = m2 = 0.5;
	}
    if (pkd->param.bCenterOfMassExpand) {
	for (j=0;j<3;++j) pkdn->r[j] = ifMass*(m1*p1->r[j] + m2*p2->r[j]);
	}
    else {
	for (j=0;j<3;++j) pkdn->r[j] = bnd->fCenter[j];
	}
    if (pkd->oNodeVelocity) {
	for (j=0;j<3;++j)	
	    pkdNodeVel(pkd,pkdn)[j]
		= ifMass*(m1*pkdNodeVel(pkd,p1)[j] + m2*pkdNodeVel(pkd,p2)[j]);
	}
    if (pkd->oNodeAcceleration) {
	for (j=0;j<3;++j)	
	    pkdNodeAccel(pkd,pkdn)[j]
		= ifMass*(m1*pkdNodeAccel(pkd,p1)[j] + m2*pkdNodeAccel(pkd,p2)[j]);
	}
    pkdn->fSoft2 = p1->fSoft2 > p2->fSoft2 ? p1->fSoft2 : p2->fSoft2;
    pkdn->uMinRung = p1->uMinRung < p2->uMinRung ? p1->uMinRung : p2->uMinRung;
    pkdn->uMaxRung = p1->uMaxRung > p2->uMaxRung ? p1->uMaxRung : p2->uMaxRung;
    pkdn->bDstActive = p1->bDstActive || p2->bDstActive;
    }


void pkdCombineCells2(PKD pkd,KDN *pkdn,KDN *p1,KDN *p2) {
    FMOMR mom;
    SPHBNDS *b1,*b2,*bn;
    float x,y,z;
    int j;

    /*
    ** Now calculate the reduced multipole moment.
    ** Shift the multipoles of each of the children
    ** to the CoM of this cell and add them up.
    */
    if (pkd->oNodeMom) {
	*pkdNodeMom(pkd,pkdn) = *pkdNodeMom(pkd,p1);
	x = (float)(p1->r[0] - pkdn->r[0]);
	y = (float)(p1->r[1] - pkdn->r[1]);
	z = (float)(p1->r[2] - pkdn->r[2]);
	momShiftFmomr(pkdNodeMom(pkd,pkdn),p1->bMax,x,y,z);

	momRescaleFmomr(pkdNodeMom(pkd,pkdn),pkdn->bMax,p1->bMax);

	mom = *pkdNodeMom(pkd,p2);
	x = (float)(p2->r[0] - pkdn->r[0]);
	y = (float)(p2->r[1] - pkdn->r[1]);
	z = (float)(p2->r[2] - pkdn->r[2]);
	momShiftFmomr(&mom,p2->bMax,x,y,z);
	momScaledAddFmomr(pkdNodeMom(pkd,pkdn),pkdn->bMax,&mom,p2->bMax);

	}
    /*
    ** Combine the special fast gas ball bounds for SPH.
    */
    if (pkd->oNodeSphBounds) {
	b1 = pkdNodeSphBounds(pkd,p1);
	b2 = pkdNodeSphBounds(pkd,p2);
	bn = pkdNodeSphBounds(pkd,pkdn);
	for (j=0;j<3;++j) bn->A.min[j] = fmin(b1->A.min[j],b2->A.min[j]);
	for (j=0;j<3;++j) bn->A.max[j] = fmax(b1->A.max[j],b2->A.max[j]);
	for (j=0;j<3;++j) bn->B.min[j] = fmin(b1->B.min[j],b2->B.min[j]);
	for (j=0;j<3;++j) bn->B.max[j] = fmax(b1->B.max[j],b2->B.max[j]);
	for (j=0;j<3;++j) bn->BI.min[j] = fmin(b1->BI.min[j],b2->BI.min[j]);
	for (j=0;j<3;++j) bn->BI.max[j] = fmax(b1->BI.max[j],b2->BI.max[j]);
    }
}

void pkdTreeBuild(PKD pkd,int nBucket, uint32_t uRoot,uint32_t uTemp) {
    int i;
#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Tree Build");
    __itt_string_handle* shMySubtask = __itt_string_handle_create("My SubTask");
#endif
#ifdef USE_ITT
     __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    pkdClearTimer(pkd,0);
    pkdStartTimer(pkd,0);

    /*
    ** The KDN at "uRoot" (e.g., ROOT) is already setup (pLower and pUpper are correct)
    ** For more information look a pkdDumpTrees and the Initialize*() routines above.
    */

    if (uTemp==0) BuildTemp(pkd,uRoot,nBucket,HUGE_VAL);
    else  BuildFromTemplate(pkd,uRoot,nBucket,uTemp);
    Create(pkd,uRoot);
    pkdStopTimer(pkd,0);

    if (uRoot == FIXROOT) {
#ifndef SINGLE_CACHES
	mdlROcache(pkd->mdl,CID_CELL2,pkdTreeNodeGetElement,pkd,pkd->iTreeNodeSize,pkd->nNodes);
	mdlROcache(pkd->mdl,CID_PARTICLE2,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
#endif
	}
    else {
	mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,pkd->iTreeNodeSize,pkd->nNodes);
	}

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
    }

void pkdTreeBuildByGroup(PKD pkd, int nBucket) {
    PARTICLE *p;
    KDN *pNode;
    BND *bnd;
    double r[3], dMin[3], dMax[3];
    int i,j,k,n,gid,gid2,iRoot;
    int iTree;

    assert(0); /* pLite is gone -- this code path needs to be tested */
    if (pkd->nNodes > 0) {
	/*
	** Close cell caching space and free up nodes.
	*/
	mdlFinishCache(pkd->mdl,CID_CELL);
	}

    /*
    ** It is only forseen that there are 4 reserved nodes at present 0-NULL, 1-ROOT, 2-UNUSED, 3-VAROOT.
    */
    pkd->nNodes = NRESERVED_NODES;

    if (pkd->hopSavedRoots == 0) {
	/* Sort particle by group, but with group 0 at the end */
	int *iGrpOffset = malloc(sizeof(int)*(pkd->nGroups+1));
	int *iGrpEnd = malloc(sizeof(int)*(pkd->nGroups+1));

	/* Count the number of particles in each group */
	for (i=0; i<=pkd->nGroups;++i) iGrpOffset[i] = 0;
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    gid = *pkdGroup(pkd,p);
	    ++iGrpOffset[gid+1];
	    }
	iGrpOffset[0] = iGrpOffset[1];
	iGrpOffset[1] = 0;
	/* Calculate starting offsets for particles in a group */
	for(i=2; i<=pkd->nGroups;++i) {
	    iGrpOffset[i] += iGrpOffset[i-1];
	    iGrpEnd[i-1] = iGrpOffset[i];
	    }

	/* Now construct the top tree node for each group */
	i = 0;
	for(gid=1; gid<pkd->nGroups;++gid) {
	    pkdTreeAllocRootNode(pkd,&iRoot);
	    pkd->hopGroups[gid].iTreeRoot = iRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    pNode->pLower = i;
	    i = iGrpOffset[gid];
	    pNode->pUpper = i - 1;
	    }

	/* Reorder the particles into group order */
	for(iTree=1;iTree<pkd->nGroups;++iTree) {
	    for(i=iGrpOffset[iTree]; i<iGrpEnd[iTree]; ) {
		p = pkdParticle(pkd,i);
		gid = *pkdGroup(pkd,p);
		if (gid==0) gid = pkd->nGroups;
		if (gid == iTree) ++i;
		else {
		    PARTICLE *p2 = pkdParticle(pkd,iGrpOffset[gid]++);
		    pkdSwapParticle(pkd,p,p2);
		    }
		}
	    }

	free(iGrpOffset);
	free(iGrpEnd);

	/* Calculate the bounds for each group */
	for (i=0;i<pkd->nLocal;) {
	    p = pkdParticle(pkd,i);
	    gid = *pkdGroup(pkd,p);
	    if (gid==0) break;
	    iRoot = pkd->hopGroups[gid].iTreeRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    bnd = pkdNodeBnd(pkd, pNode);
	    
	    pNode->iLower = 0;
	    assert(pNode->pLower == i);

	    for (j=0;j<3;++j) dMin[j] = dMax[j] = pkdPos(pkd,p,j);
	    for(p = pkdParticle(pkd,++i); i<pkd->nLocal && *pkdGroup(pkd,p)==gid; ++i) {
		for (j=0;j<3;++j) r[j] = pkdPos(pkd,p,j);
		pkdMinMax(r,dMin,dMax);
		}
	    for (j=0;j<3;++j) {
		bnd->fCenter[j] = 0.5*(dMin[j] + dMax[j]);
		bnd->fMax[j] = 0.5*(dMax[j] - dMin[j]);
		}
	    assert(pNode->pUpper == i-1);
	    }
	/* We can use this to quickly rebuild the trees */
	pkd->hopSavedRoots = pkd->nNodes;
	}
    else {
	pkd->nNodes = pkd->hopSavedRoots;
	for(gid2=1; gid2<pkd->nGroups;++gid2) {
	    if (!pkd->hopGroups[gid2].bNeedGrav) continue;
	    iRoot = pkd->hopGroups[gid2].iTreeRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    pNode->iLower = 0;
	    bnd = pkdNodeBnd(pkd, pNode);
	    n = pNode->pUpper;
	    for(i=pNode->pLower; i<=n; ) {
		p = pkdParticle(pkd,i);
		gid = *pkdGroup(pkd,p);
		if (gid) {
		    assert(gid==gid2);
		    if (i==pNode->pLower) {
			for (j=0;j<3;++j) dMin[j] = dMax[j] = pkdPos(pkd,p,j);
			}
		    else {
			for (j=0;j<3;++j) r[j] = pkdPos(pkd,p,j);
			pkdMinMax(r,dMin,dMax);
			}
		    ++i;
		    }
		else {
		    PARTICLE *p2 = pkdParticle(pkd,n--);
		    pkdSwapParticle(pkd,p,p2);
		    }
		}
	    assert(k==n+1);
	    pNode->pUpper = n;
	    for (j=0;j<3;++j) {
		bnd->fCenter[j] = 0.5*(dMin[j] + dMax[j]);
		bnd->fMax[j] = 0.5*(dMax[j] - dMin[j]);
		}
	    }
	}


    for(gid=1; gid<pkd->nGroups; ++gid)
	if (pkd->hopGroups[gid].bNeedGrav)
	    BuildTemp(pkd,pkd->hopGroups[gid].iTreeRoot,nBucket,HUGE_VAL);
    for(gid=1; gid<pkd->nGroups; ++gid)
	if (pkd->hopGroups[gid].bNeedGrav)
	    Create(pkd,pkd->hopGroups[gid].iTreeRoot);
    /*
    ** Finally activate a read only cache for remote access.
    */
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,
	pkd->iTreeNodeSize,pkd->nNodes);


    }

/*
** Hopefully we can bypass this step once we figure out how to do the
** Multipole Ewald with reduced multipoles.
*/
void pkdCalcRoot(PKD pkd,uint32_t uRoot,double *com,MOMC *pmom) {
    PARTICLE *p;
    FLOAT xr = com[0];
    FLOAT yr = com[1];
    FLOAT zr = com[2];
    FLOAT x,y,z;
    FLOAT fMass;
    MOMC mc;
    KDN *kdn = pkdTreeNode(pkd,uRoot);
    int i = kdn->pLower;
    if (kdn->pLower > kdn->pUpper) momClearMomc(pmom);
    else {
	p = pkdParticle(pkd,i);
	x = pkdPos(pkd,p,0) - xr;
	y = pkdPos(pkd,p,1) - yr;
	z = pkdPos(pkd,p,2) - zr;
	fMass = pkdMass(pkd,p);
	momMakeMomc(pmom,fMass,x,y,z);
	for (++i;i<=kdn->pUpper;++i) {
	    p = pkdParticle(pkd,i);
	    fMass = pkdMass(pkd,p);
	    x = pkdPos(pkd,p,0) - xr;
	    y = pkdPos(pkd,p,1) - yr;
	    z = pkdPos(pkd,p,2) - zr;
	    momMakeMomc(&mc,fMass,x,y,z);
	    momAddMomc(pmom,&mc);
	    }
	}
    }


void pkdDistribRoot(PKD pkd,double *r,MOMC *pmom) {
    pkd->ew.r[0] = r[0];
    pkd->ew.r[1] = r[1];
    pkd->ew.r[2] = r[2];
    pkd->ew.mom = *pmom;
    }
