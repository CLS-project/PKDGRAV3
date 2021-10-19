/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <assert.h>
#include "pkd.h"
#include "gravity/moments.h"

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif


uint32_t pkdDistribTopTree(PKD pkd, uint32_t uRoot, uint32_t nTop, KDN *pTop, int allocateMemory) {
    int i, iTop;
    KDN *pLocalRoot = pkdTreeNode(pkd,uRoot);

    if (allocateMemory) {
        iTop = pkd->iTopTree[uRoot] = pkdTreeAllocNodes(pkd, nTop);
    } else {
        iTop = pkd->iTopTree[uRoot];
    }
    for(i=0; i<nTop; ++i) {
	KDN *pNode = pkdNode(pkd,pTop,i);
	KDN *pLocal = pkdTreeNode(pkd,iTop+i);
	pkdCopyNode(pkd,pLocal,pNode);
	assert(pLocal->bTopTree);
	pLocal->bGroup = 0;
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
    pRoot->bTopTree = 0;
    pRoot->bRemote = 0;
    pkdNodeSetBnd(pkd, pRoot, &pkd->bnd);
    return pRoot;
    }

/*
** Creates a single root at ROOT with only marked particles
*/
void pkdTreeInitMarked(PKD pkd) {
    KDN *pRoot   = InitializeRootCommon(pkd,ROOT);
    KDN *pRootFix = InitializeRootCommon(pkd,FIXROOT);
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
    pRoot->pUpper = i-1;
    pRootFix->pLower = i;
    pRootFix->pUpper = pkd->nLocal - 1;

    pkd->nNodes = NRESERVED_NODES;
    pRoot->iLower = 0;
    pRoot->bGroup = 1;
    pRootFix->iLower = 0;
    pRootFix->bGroup = 1;
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
	pRoot->bGroup = 1;
	}
    /* Just rebuilding (active) ROOT. Truncate it. pLower and pUpper are still valid. */
    else if (bOnlyVA) {
	if (pRoot->iLower) {
	    assert(pRoot->iLower >= NRESERVED_NODES);
	    pkd->nNodes = pRoot->iLower; /* This effectively truncates the nodes used by this tree */
	    pRoot->iLower = 0;
	    pRoot->bGroup = 1;
	    }
	}

    /* Need to build two trees which is more complicated. */
    else {
	KDN *pRootFix = InitializeRootCommon(pkd,FIXROOT);
	PARTICLE *p, *pVA = NULL;
	local_t iLast;
	int i;

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
	pRoot->bGroup = 1;
	pRootFix->iLower = 0;
	pRootFix->bGroup = 1;

	pRootFix->pLower = 0;
	pRootFix->pUpper = iLast;
	pRoot->pLower = iLast + 1;
	pRoot->pUpper = pkd->nLocal - 1;

//	for(i=pRootFix->pLower; i<=pRootFix->pUpper; ++i) assert(pkdParticle(pkd,i)->uRung<=uRungDD);
//	for(i=pRoot->pLower; i<=pRoot->pUpper; ++i) assert(pkdParticle(pkd,i)->uRung>uRungDD);
	}
    }

#define MIN_SRATIO    0.05


template<class split_t>
static split_t PosRaw(PKD pkd,PARTICLE *p, int d) {
    return reinterpret_cast<split_t *>(pkdField(p,pkd->oFieldOffset[oPosition]))[d];
    }

/*
** Partition the particles between pLower and pUpper (inclusive)
** around "Split" in dimension "d". Return the index of the split.
*/
template<class split_t>
static int PartPart(PKD pkd,int pLower,int pUpper,int d,split_t Split) {
	PARTICLE *pi = pkdParticle(pkd,pLower);
	PARTICLE *pj = pkdParticle(pkd,pUpper);
	while (pi <= pj) {
	    if (PosRaw<split_t>(pkd,pi,d) < Split) pi = (PARTICLE *)(((char *)pi) + pkdParticleSize(pkd));
	    else break;
	    }
	while (pi <= pj) {
	    if (Split < PosRaw<split_t>(pkd,pj,d)) pj = (PARTICLE *)(((char *)pj) - pkdParticleSize(pkd));
	    else break;
	    }
	if (pi < pj) {
	    pkdSwapParticle(pkd,pi,pj);
	    while (1) {
		while (PosRaw<split_t>(pkd,(pi = (PARTICLE *)(((char *)pi) + pkdParticleSize(pkd))),d) < Split);
		while (Split < PosRaw<split_t>(pkd,(pj = (PARTICLE *)(((char *)pj) - pkdParticleSize(pkd))),d));
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
void BuildTemp(PKD pkd,int iNode,int M,int nGroup,double dMaxMax) {
    KDN *pNode = pkdTreeNode(pkd,iNode);
    KDN *pLeft, *pRight;
    double lrMax;
    int *S;		/* this is the stack */
    int s,ns;
    int iLeft,iRight;
    int d;
    int i;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;

    pNode->iDepth = 0;
    pNode->iSplitDim = 3;

    // Single bucket? We are done.
    if (pNode->pUpper - pNode->pLower + 1 <= M) return;

    /*
    ** Allocate stack!
    */
    ns = TEMP_S_INCREASE;
    s = 0;
    S = CAST(int*,malloc(ns*sizeof(int)));
    assert(S != NULL);

    Bound bnd = pkdNodeGetBnd(pkd,pNode);
    while (1) {
	// Begin new stage! Calculate the appropriate fSplit.
	// Pick longest dimension and split it in half.
	// Calculate the new left and right cells. We will use
	// either the left, the right or (usually) both.
	pNode->iSplitDim = d = bnd.maxdim();
	Bound lbnd,rbnd;
	std::tie(lbnd,rbnd) = bnd.split(d);
	lrMax = 0.5*lbnd.maxside();
	// Now start the partitioning of the particles about
	// fSplit on dimension given by d.
	if (pkd->bIntegerPosition) {
	    int32_t Split = pkdDblToIntPos(pkd,bnd.center(d));
	    i = PartPart(pkd,pNode->pLower,pNode->pUpper,d,Split);
	    }
	else i = PartPart(pkd,pNode->pLower,pNode->pUpper,d,bnd.center(d));
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
	    pLeft->iSplitDim = 3;
	    pRight = pkdTreeNode(pkd,iRight);
	    assert(iRight & 1);
	    pRight->bTopTree = 0;
	    pRight->bRemote = 0;
	    pRight->pLower = i;
	    pRight->pUpper = pNode->pUpper;
	    pRight->iDepth = pNode->iDepth+1;
	    pRight->iSplitDim = 3;
	    pNode->iLower = iLeft;
	    pNode->bGroup = pNode->pUpper - pNode->pLower < nGroup;

            pkdNodeSetBnd(pkd, pLeft, &lbnd);
            pkdNodeSetBnd(pkd, pRight, &rbnd);

	    /*
	    ** Now figure out which subfile to process next.
	    */
	    if (lrMax > dMaxMax) {
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
		pRight->bGroup = 1;
		++nBucket;
		}
	    else if (rc) {
		/*
		** Left must be a bucket in this case!
		*/
		iNode = iRight;   /* process upper subfile */
		pLeft->iLower = 0;
		pLeft->bGroup = 1;
		++nBucket;
		}
	    else {
		/*
		** Both are buckets (we need to pop from the stack to get the next subfile.
		*/
		pLeft->iLower = 0;
		pLeft->bGroup = 1;
		++nBucket;
		pRight->iLower = 0;
		pRight->bGroup = 1;
		++nBucket;
		if (s) iNode = S[--s];		/* pop tn */
		else break;
		}
	    }
	else {
	    ++pNode->iDepth;
	    // No nodes allocated: just change the bounds
	    if (nl > 0) {
		pkdNodeSetBnd(pkd, pNode, &lbnd);
		lc = (lrMax>dMaxMax || nl > M); /* this condition means the node is not a bucket */
		if (!lc) {
		    pNode->iLower = 0;
		    pNode->bGroup = 1;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    else {
		pkdNodeSetBnd(pkd, pNode, &rbnd);
		rc = (lrMax>dMaxMax || nr > M);
		if (!rc) {
		    pNode->iLower = 0;
		    pNode->bGroup = 1;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    }
	pNode = pkdTreeNode(pkd,iNode);
        bnd = pkdNodeGetBnd(pkd, pNode);
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
void BuildFromTemplate(PKD pkd,int iNode,int M,int nGroup,int iTemplate) {

    // The tree we are building
    KDN *pNode, *pLeft, *pRight;
    int iLeft, iRight;

    // The "template" tree
    KDN *pTemp, *ptLeft, *ptRight;

    int d;
    int i, nr, nl;
    int ns,s;
    double dMaxMax;
    struct buildStack {
	KDN *pTemp;
	int iNode;
	} *S;
    // If we have an empty tree root then we are "done"
    pNode = pkdTreeNode(pkd,iNode);
    pNode->iDepth = 0;
    pNode->iLower = 0;
    pNode->bGroup = 1;
    if (pNode->pUpper - pNode->pLower + 1 == 0) return;

    // Setup stack and push the two root cells of each tree.
    ns = TEMP_S_INCREASE;
    S = (struct buildStack *)malloc(ns*sizeof(struct buildStack));
    assert(S != NULL);
    s = 1;
    S[0].iNode = iNode;
    S[0].pTemp = pkdTreeNode(pkd,iTemplate);

    dMaxMax = S[0].pTemp->bMax;


    do {
	// Pop the next cells to process
	assert(s);
	--s;
	pTemp = S[s].pTemp;
	iNode = S[s].iNode;

	Bound bnd = pkdNodeGetBnd(pkd,pNode = pkdTreeNode(pkd,iNode));

	// Follow the template tree to wherever it leads.
	while(pTemp->iLower) {
	    d = pTemp->iSplitDim;
	    if (pTemp->iDepth != pNode->iDepth || d>2) break;

	    // Split is between left and right child nodes in the given dimension
	    Bound ltbnd = pkdNodeGetBnd(pkd, ptLeft=pkdTreeNode(pkd,pTemp->iLower));
	    Bound rtbnd = pkdNodeGetBnd(pkd, ptRight=pkdTreeNode(pkd,pTemp->iLower+1));
	    auto dSplit = 0.5 * (ltbnd.upper(d) + rtbnd.lower(d));

	    // Partition the particles on either side of the split
	    if (pkd->bIntegerPosition) {
		int32_t Split = pkdDblToIntPos(pkd,dSplit);
		i = PartPart(pkd,pNode->pLower,pNode->pUpper,d,Split);
		}
	    else i = PartPart(pkd,pNode->pLower,pNode->pUpper,d,dSplit);
	    nl = i - pNode->pLower;
	    nr = pNode->pUpper - i + 1;
	    assert(nl>0 || nr>0);


	    // Calculate bounding regions
	    Bound lbnd,rbnd;
	    std::tie(lbnd,rbnd) = bnd.split(d,dSplit);
	    assert(rbnd.width(d) > 0.0);
	    assert(lbnd.width(d) > 0.0);

	    if (nl==0) { // Particles on the right only
		++pNode->iDepth;
		pkdNodeSetBnd(pkd,pNode, &rbnd);
		pTemp = ptRight;
		continue;
		}
	    else if (nr==0) { // Particles on the left only
		++pNode->iDepth;
		pkdNodeSetBnd(pkd,pNode, &lbnd);
		pTemp = ptLeft;
		continue;
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
	    pLeft->bGroup = 1;
	    pLeft->iDepth = pNode->iDepth + 1;
	    assert(pLeft->pLower <= pLeft->pUpper);

	    pRight = pkdTreeNode(pkd,iRight);
	    assert(iRight & 1);
	    pRight->bTopTree = 0;
	    pRight->bRemote = 0;
	    pRight->pLower = i;
	    pRight->pUpper = pNode->pUpper;
	    pRight->iLower = 0;
	    pRight->bGroup = 1;
	    pRight->iDepth = pNode->iDepth + 1;
	    assert(pRight->pLower <= pRight->pUpper);

	    pNode->iLower = iLeft;
	    pNode->bGroup = (pNode->pUpper-pNode->pLower < nGroup) && pTemp->bGroup;

	    iNode = -1;
	    pNode = NULL;

	    pkdNodeSetBnd(pkd, pLeft, &lbnd);
	    pkdNodeSetBnd(pkd, pRight, &rbnd);

	    // push right and continue left
	    if ( s+1 >= ns ) {
		assert( s+1 == ns );
		ns += TEMP_S_INCREASE;
		S = (struct buildStack *)realloc(S,ns*sizeof(struct buildStack));
		}
	    S[s].pTemp = ptRight;
	    S[s].iNode = iRight;
	    ++s;
	    pTemp = ptLeft;
	    iNode = iLeft;
	    pNode = pkdTreeNode(pkd,iNode);
	    bnd = pkdNodeGetBnd(pkd,pNode);
	    }
	
	// Bucket in the template tree: Now just build, but set a sensible maximum cell size
	BuildTemp(pkd,iNode,M,nGroup,dMaxMax * pow(2.0,-pNode->iDepth/3.0));
	} while(s);

    free(S);
    }



static vel_t zeroV[3] = {0.0,0.0,0.0};
static float  zeroF[3] = {0.0,0.0,0.0};

void Create(PKD pkd,int iRoot,double ddHonHLimit) {
    int iNode = iRoot;
    PARTICLE *p;
    KDN *pkdn,*pkdl,*pkdu;
    FMOMR mom;
    SPHBNDS *bn;
    double kdn_r[3];
    double fSoft,x,y,z,ax,ay,az,ft[3],d2,d2Max,dih2,bmin,b;
    double dx,dy,dz,fBoBr,fBoBxCenter,fBoByCenter,fBoBzCenter;
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
    pkdn->bHasMarked = 0;
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
        Bound bnd = pkdNodeGetBnd(pkd, pkdn);
	/*
	** Before squeezing the bounds, calculate a minimum b value based on the splitting bounds alone.
	** This gives us a better feel for the "size" of a bucket with only a single particle.
	*/
	bmin = bnd.minside();
	/*
	** Now shrink wrap the bucket bounds.
	*/
	pj = pkdn->pLower;
	p = pkdParticle(pkd,pj);
#if defined(__AVX__) && defined(USE_SIMD)
	if (pkd->bIntegerPosition) {
	    __m128i ivmin, ivmax;
	    ivmin = ivmax = pkdGetPosRaw(pkd,p);
	    for (++pj;pj<=pkdn->pUpper;++pj) {
		p = pkdParticle(pkd,pj);
		__m128i v = pkdGetPosRaw(pkd,p);
		ivmin = _mm_min_epi32(ivmin,v);
		ivmax = _mm_max_epi32(ivmax,v);
		}
	    union {
		double d[4];
		__m256d p;
    		} vmin,vmax;
	    vmin.p = _mm256_mul_pd(_mm256_cvtepi32_pd(ivmin),_mm256_set1_pd(1.0/INTEGER_FACTOR) );
	    vmax.p = _mm256_mul_pd(_mm256_cvtepi32_pd(ivmax),_mm256_set1_pd(1.0/INTEGER_FACTOR) );
	    bnd = Bound(vmin.d,vmax.d);
	    }
	else
#endif
	{
	    pkdGetPos1(pkd,p,ft);
	    double dMin[3] {ft[0],ft[1],ft[2]};
	    double dMax[3] {ft[0],ft[1],ft[2]};
	    for (++pj;pj<=pkdn->pUpper;++pj) {
		p = pkdParticle(pkd,pj);
		pkdGetPos1(pkd,p,ft);
		for (d=0;d<3;++d) {
		    dMin[d] = std::min(dMin[d],ft[d]);
		    dMax[d] = std::max(dMax[d],ft[d]);
		}
		}
	    bnd = Bound(dMin,dMax);
	    }
        pkdNodeSetBnd(pkd, pkdn, &bnd);
	pj = pkdn->pLower;
	p = pkdParticle(pkd,pj);
	a = pkd->oFieldOffset[oAcceleration] ? pkdAccel(pkd,p) : zeroF;
	m = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	v = pkd->oFieldOffset[oVelocity] ? pkdVel(pkd,p) : zeroV;
	fMass = m;
	dih2 = fSoft;
	pkdGetPos3(pkd,p,x,y,z);

    /* initialize ball of balls */
    fBoBr = pkdBall(pkd,p);
    fBoBxCenter = x;
    fBoByCenter = y;
    fBoBzCenter = z;
    /* initialize marked flag */
    pkdn->bHasMarked = p->bMarked;

	x *= m;
	y *= m;
	z *= m;
	vx = m*v[0];
	vy = m*v[1];
	vz = m*v[2];
	ax = m*a[0];
	ay = m*a[1];
	az = m*a[2];
	pkdn->uMinRung = pkdn->uMaxRung = p->uRung;
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    a = pkd->oFieldOffset[oAcceleration] ? pkdAccel(pkd,p) : zeroF;
	    m = pkdMass(pkd,p);
	    fSoft = pkdSoft(pkd,p);
	    v = pkd->oFieldOffset[oVelocity] ? pkdVel(pkd,p) : zeroV;
	    fMass += m;
	    if (fSoft>dih2) dih2=fSoft;
	    pkdGetPos1(pkd,p,ft);

        /* calculate ball of balls */
        // dx = bnd.fCenter[0] - ft[0];
        // dy = bnd.fCenter[1] - ft[1];
        // dz = bnd.fCenter[2] - ft[2];
        // fBoBrq = sqrt(dx*dx + dy*dy + dz*dz) + pkdBall(pkd,p);
        // fBoBr = fBoBrq > fBoBr ? fBoBrq : fBoBr;
        blitz::TinyVector<float,3> p1fBoBCenter = blitz::TinyVector<float,3>(fBoBxCenter,fBoByCenter,fBoBzCenter);
        blitz::TinyVector<float,3> p2fBoBCenter = blitz::TinyVector<float,3>(ft[0],ft[1],ft[2]);
        blitz::TinyVector<float,3> difference = p1fBoBCenter - p2fBoBCenter;
        blitz::TinyVector<float,3> direction = difference / sqrtf(blitz::dot(difference,difference));
        blitz::TinyVector<float,3> point1 = p1fBoBCenter + direction * fBoBr;
        blitz::TinyVector<float,3> point2 = p2fBoBCenter - direction * pkdBall(pkd,p);
        blitz::TinyVector<float,3> midpoint = (point1 - point2) / 2.0f + point2;
        fBoBr = sqrtf(blitz::dot((point1 - point2) / 2.0f,(point1 - point2) / 2.0f));
        fBoBxCenter = midpoint[0];
        fBoByCenter = midpoint[1];
        fBoBzCenter = midpoint[2];
        if (p->bMarked) pkdn->bHasMarked = 1;

	    x += m*ft[0];
	    y += m*ft[1];
	    z += m*ft[2];
	    vx += m*v[0];
	    vy += m*v[1];
	    vz += m*v[2];
	    ax += m*a[0];
	    ay += m*a[1];
	    az += m*a[2];
	    if ( p->uRung > pkdn->uMaxRung ) pkdn->uMaxRung = p->uRung;
	    if ( p->uRung < pkdn->uMinRung ) pkdn->uMinRung = p->uRung;
	    }
	m = 1.0f / fMass;
	kdn_r[0] = m*x;
	kdn_r[1] = m*y;
	kdn_r[2] = m*z;
	pkdNodeSetPos1(pkd,pkdn,kdn_r);
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
    pkdn->fBoBr2 = fBoBr*fBoBr;
    pkdn->fBoBxCenter = fBoBxCenter;
    pkdn->fBoByCenter = fBoByCenter;
    pkdn->fBoBzCenter = fBoBzCenter;
	d2Max = 0.0;
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
#if defined(__AVX__) && defined(USE_SIMD)
	    if (pkd->bIntegerPosition) {
		__m256d v = _mm256_sub_pd(pkdGetPos(pkd,p),_mm256_setr_pd(kdn_r[0],kdn_r[1],kdn_r[2],0.0));
		v = _mm256_mul_pd(v,v);
		__m128d t0 = _mm256_extractf128_pd(v,0);
		__m128d t2 = _mm256_extractf128_pd(v,1);
		__m128d t1 = _mm_unpackhi_pd(t0,t0);
		t0 = _mm_add_sd(t0,t2);
		t0 = _mm_add_sd(t0,t1);
		d2Max = _mm_cvtsd_f64(_mm_max_sd(t0,_mm_set_sd(d2Max)));
	    }
	    else
#endif
	    {
		pkdGetPos1(pkd,p,ft);
		x = ft[0] - kdn_r[0];
		y = ft[1] - kdn_r[1];
		z = ft[2] - kdn_r[2];
		d2 = x*x + y*y + z*z;
		/*
		** Update bounding ball and softened bounding ball.
		*/
		d2Max = (d2 > d2Max)?d2:d2Max;
		}
	    }
#ifdef USE_MAXSIDE
	b = bnd.maxside();
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
		pkdGetPos1(pkd,p,ft);
		x = ft[0] - kdn_r[0];
		y = ft[1] - kdn_r[1];
		z = ft[2] - kdn_r[2];
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
		    double r[3];
		    pkdGetPos1(pkd,p,r);
		    /*
		    ** This first ball bound over all gas particles is only used for remote searching.
		    */
		    for (d=0;d<3;++d) bn->B.min[d] = fmin(bn->B.min[d],r[d] - (1+ddHonHLimit)*fBall);
		    for (d=0;d<3;++d) bn->B.max[d] = fmax(bn->B.max[d],r[d] + (1+ddHonHLimit)*fBall);
		    if (pkdIsActive(pkd,p)) {
			for (d=0;d<3;++d) bn->A.min[d] = fmin(bn->A.min[d],r[d]);
			for (d=0;d<3;++d) bn->A.max[d] = fmax(bn->A.max[d],r[d]);
		    }
		    else {
			for (d=0;d<3;++d) bn->BI.min[d] = fmin(bn->BI.min[d],r[d] - (1+ddHonHLimit)*fBall);
			for (d=0;d<3;++d) bn->BI.max[d] = fmax(bn->BI.max[d],r[d] + (1+ddHonHLimit)*fBall);
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
            bnd = pkdNodeGetBnd(pkd, pkdn);
	    /*
	    ** Before squeezing the bounds, calculate a minimum b value based on the splitting bounds alone.
	    ** This gives us a better feel for the "size" of a bucket with only a single particle.
	    */
	    bmin = bnd.minside();
	    pj = pkdn->pLower;
	    pkdl = pkdTreeNode(pkd,pkdn->iLower);
	    pkdu = pkdTreeNode(pkd,pkdn->iLower + 1);
	    pkdCombineCells1(pkd,pkdn,pkdl,pkdu);
	    pkdNodeGetPos(pkd,pkdn,kdn_r);
	    if (pkdn->pUpper - pj < NMAX_OPENCALC) {
		assert(pj<=pkdn->pUpper);
		d2Max = 0;
		for (;pj<=pkdn->pUpper;++pj) {
		    p = pkdParticle(pkd,pj);
#if defined(__AVX__) && defined(USE_SIMD)
		    if (pkd->bIntegerPosition) {
			__m256d v = _mm256_sub_pd(pkdGetPos(pkd,p),_mm256_setr_pd(kdn_r[0],kdn_r[1],kdn_r[2],0.0));
			v = _mm256_mul_pd(v,v);
			__m128d t0 = _mm256_extractf128_pd(v,0);
			__m128d t2 = _mm256_extractf128_pd(v,1);
			__m128d t1 = _mm_unpackhi_pd(t0,t0);
			t0 = _mm_add_sd(t0,t2);
			t0 = _mm_add_sd(t0,t1);
			d2Max = _mm_cvtsd_f64(_mm_max_sd(t0,_mm_set_sd(d2Max)));
			}
		    else
#endif
		    {
			pkdGetPos3(pkd,p,x,y,z);
			x -= kdn_r[0];
			y -= kdn_r[1];
			z -= kdn_r[2];
			d2 = x*x + y*y + z*z;
			d2Max = (d2 > d2Max)?d2:d2Max;
			}
		    }
		assert(d2Max>0);
		/*
		** Now determine the opening radius for gravity.
		*/
#ifdef USE_MAXSIDE
		b = bnd.maxside();
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
    double m1,m2,ifMass;
    int j;
    double kdn_r[3];

    Bound p1bnd = pkdNodeGetBnd(pkd, p1);
    Bound p2bnd = pkdNodeGetBnd(pkd, p2);
    Bound bnd = p1bnd.combine(p2bnd);
    pkdNodeSetBnd(pkd, pkdn,&bnd);
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
    double p1_r[3], p2_r[3];
    pkdNodeGetPos(pkd,p1,p1_r);
    pkdNodeGetPos(pkd,p2,p2_r);
    for (j=0;j<3;++j) kdn_r[j] = ifMass*(m1*p1_r[j] + m2*p2_r[j]);
    pkdNodeSetPos1(pkd,pkdn,kdn_r);
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

    /* Combine ball of balls */
    blitz::TinyVector<float,3> p1fBoBCenter = blitz::TinyVector<float,3>(p1->fBoBxCenter,p1->fBoByCenter,p1->fBoBzCenter);
    blitz::TinyVector<float,3> p2fBoBCenter = blitz::TinyVector<float,3>(p2->fBoBxCenter,p2->fBoByCenter,p2->fBoBzCenter);
    blitz::TinyVector<float,3> difference = p1fBoBCenter - p2fBoBCenter;
    blitz::TinyVector<float,3> direction = difference / sqrtf(blitz::dot(difference,difference));
    blitz::TinyVector<float,3> point1 = p1fBoBCenter + direction * sqrt(p1->fBoBr2);
    blitz::TinyVector<float,3> point2 = p2fBoBCenter - direction * sqrt(p2->fBoBr2);
    blitz::TinyVector<float,3> midpoint = (point1 - point2) / 2.0f + point2;
    pkdn->fBoBr2 = blitz::dot((point1 - point2) / 2.0f,(point1 - point2) / 2.0f);
    pkdn->fBoBxCenter = midpoint[0];
    pkdn->fBoByCenter = midpoint[1];
    pkdn->fBoBzCenter = midpoint[2];
    /* Combine marked flag */
    pkdn->bHasMarked = p1->bHasMarked || p2->bHasMarked;
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
	double kdn_r[3], p1_r[3], p2_r[3];
	pkdNodeGetPos(pkd,pkdn,kdn_r);
	pkdNodeGetPos(pkd,p1,p1_r);
	pkdNodeGetPos(pkd,p2,p2_r);
	*pkdNodeMom(pkd,pkdn) = *pkdNodeMom(pkd,p1);
	x = (float)(p1_r[0] - kdn_r[0]);
	y = (float)(p1_r[1] - kdn_r[1]);
	z = (float)(p1_r[2] - kdn_r[2]);
	momShiftFmomr(pkdNodeMom(pkd,pkdn),p1->bMax,x,y,z);

	momRescaleFmomr(pkdNodeMom(pkd,pkdn),pkdn->bMax,p1->bMax);

	mom = *pkdNodeMom(pkd,p2);
	x = (float)(p2_r[0] - kdn_r[0]);
	y = (float)(p2_r[1] - kdn_r[1]);
	z = (float)(p2_r[2] - kdn_r[2]);
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

void pkdTreeBuild(PKD pkd,int nBucket, int nGroup, uint32_t uRoot,uint32_t uTemp, double ddHonHLimit) {
#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Tree Build");
    __itt_string_handle* shMySubtask = __itt_string_handle_create("My SubTask");
#endif
#ifdef USE_ITT
     __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif

    /*
    ** The KDN at "uRoot" (e.g., ROOT) is already setup (pLower and pUpper are correct)
    ** For more information look a pkdDumpTrees and the Initialize*() routines above.
    */

    if (uTemp==0) BuildTemp(pkd,uRoot,nBucket,nGroup,HUGE_VAL);
    else  BuildFromTemplate(pkd,uRoot,nBucket,nGroup,uTemp);
    Create(pkd,uRoot,ddHonHLimit);

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

/*
** The array iGrpOffset[i] passed in here must have 2*pkd->nGroups entries!
*/
void pkdGroupOrder(PKD pkd,uint32_t *iGrpOffset) {
    uint32_t i,gid,iTree;
    uint32_t *iGrpEnd = &iGrpOffset[pkd->nGroups]; /* tricky because the 0th element is not used! */

    /* Count the number of particles in each group */
    for (i=0;i<=pkd->nGroups;++i) iGrpOffset[i] = 0;
    for (i=0;i<pkd->nLocal;++i) {
	auto p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	++iGrpOffset[gid+1];
	}
    iGrpOffset[1] = 0;
    /* Calculate starting offsets for particles in a group */
    for(i=2; i<=pkd->nGroups;++i) {
	iGrpOffset[i] += iGrpOffset[i-1];
	iGrpEnd[i-1] = iGrpOffset[i];
	}
    /* Reorder the particles into group order */
    for (iTree=1;iTree<pkd->nGroups;++iTree) {
	for (i=iGrpOffset[iTree];i<iGrpEnd[iTree];) {
	    auto p = pkdParticle(pkd,i);
	    gid = pkdGetGroup(pkd,p);
	    if (!gid) gid = pkd->nGroups;
	    if (gid == iTree) ++i;
	    else {
		PARTICLE *p2 = pkdParticle(pkd,iGrpOffset[gid]++);
		pkdSwapParticle(pkd,p,p2);
		}
	    }
	}
    }


void pkdTreeBuildByGroup(PKD pkd, int nBucket, int nGroup) {
    PARTICLE *p;
    KDN *pNode;
    double r[3], dMin[3], dMax[3];
    int i,j,k,n,gid,gid2,iRoot;
    int iTree;

    /*
    ** Should use the pkdGroupOrder function to just reorder the particles 
    ** without building the trees first. This is useful for some of the 
    ** groupstats functions.
    */
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
	int *iGrpOffset = CAST(int *,malloc(sizeof(int)*(pkd->nGroups+1)));
	int *iGrpEnd = CAST(int *,malloc(sizeof(int)*(pkd->nGroups+1)));

	/* Count the number of particles in each group */
	for (i=0; i<=pkd->nGroups;++i) iGrpOffset[i] = 0;
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    gid = pkdGetGroup(pkd,p);
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
		gid = pkdGetGroup(pkd,p);
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
	    pkdGetPos1(pkd,p,r);
	    gid = pkdGetGroup(pkd,p);
	    if (gid==0) break;
	    iRoot = pkd->hopGroups[gid].iTreeRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    
	    pNode->iLower = 0;
	    pNode->bGroup = 1;
	    assert(pNode->pLower == i);

	    for (j=0;j<3;++j) dMin[j] = dMax[j] = r[j];
	    for(p = pkdParticle(pkd,++i); i<pkd->nLocal && pkdGetGroup(pkd,p)==gid; ++i) {
		pkdGetPos1(pkd,p,r);
		pkdMinMax(r,dMin,dMax);
		}
	    pkdNodeSetBndMinMax(pkd,pNode,dMin,dMax);
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
	    pNode->bGroup = 1;
	    n = pNode->pUpper;
	    for(i=pNode->pLower; i<=n; ) {
		p = pkdParticle(pkd,i);
		pkdGetPos1(pkd,p,r);
		gid = pkdGetGroup(pkd,p);
		if (gid) {
		    assert(gid==gid2);
		    if (i==pNode->pLower) {
			for (j=0;j<3;++j) dMin[j] = dMax[j] = r[j];
			}
		    else {
			for (j=0;j<3;++j) r[j] = r[j];
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
	    Bound bnd(dMin,dMax);
	    pkdNodeSetBnd(pkd, pNode, &bnd);
	    }
	}


    for(gid=1; gid<pkd->nGroups; ++gid)
	if (pkd->hopGroups[gid].bNeedGrav)
	    BuildTemp(pkd,pkd->hopGroups[gid].iTreeRoot,nBucket,nGroup,HUGE_VAL);
    for(gid=1; gid<pkd->nGroups; ++gid)
	if (pkd->hopGroups[gid].bNeedGrav)
	    Create(pkd,pkd->hopGroups[gid].iTreeRoot,0.0);
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
    double xr = com[0];
    double yr = com[1];
    double zr = com[2];
    double x,y,z;
    double fMass;
    MOMC mc;
    KDN *kdn = pkdTreeNode(pkd,uRoot);
    int i = kdn->pLower;
    if (kdn->pLower > kdn->pUpper) momClearMomc(pmom);
    else {
	p = pkdParticle(pkd,i);
	pkdGetPos3(pkd,p,x,y,z);
	x -= xr;
	y -= yr;
	z -= zr;
	fMass = pkdMass(pkd,p);
	momMakeMomc(pmom,fMass,x,y,z);
	for (++i;i<=kdn->pUpper;++i) {
	    p = pkdParticle(pkd,i);
	    pkdGetPos3(pkd,p,x,y,z);
	    fMass = pkdMass(pkd,p);
	    x -= xr;
	    y -= yr;
	    z -= zr;
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

void pkdTreeUpdateMarkedFlags(PKD pkd,uint32_t uRoot) {
    KDN *c,*cLow,*cUp;
    int id,idLo,iCellLo,idUp,iCellUp;
    PARTICLE *p;
    int pj;

    c = pkdTreeNode(pkd,uRoot);

    if (c->bGroup) {
        for (pj=c->pLower;pj<=c->pUpper;++pj) {
		    p = pkdParticle(pkd,pj);
            if (p->bMarked) {
                c->bHasMarked = 1;
                break;
            }
        }
    } else {
        pkdGetChildCells(c,id,idLo,iCellLo,idUp,iCellUp);
        pkdTreeUpdateMarkedFlags(pkd,iCellLo);
        pkdTreeUpdateMarkedFlags(pkd,iCellUp);
        cLow = pkdTreeNode(pkd,iCellLo);
        cUp = pkdTreeNode(pkd,iCellUp);
        c->bHasMarked = cLow->bHasMarked || cUp->bHasMarked;
    }
    }