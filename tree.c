#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <assert.h>
#include "pkd.h"
#include "moments.h"
#include "floattype.h"

typedef struct pLite {
    FLOAT r[3];
    int i;
    unsigned int iActive;
    } PLITE;

/*
** M is the bucket size.
*/
int BuildTemp(PKD pkd,int M,int bSqueeze,int bExcludeVeryActive)
    {
    PLITE *p;
    PLITE t;
    FLOAT fSplit;
    FLOAT fMin[3],fMax[3];
    BND *pbnd;
    int *S;		/* this is the stack */
    int s,ns;
    int iNode,iLeft,iRight;
    int nNodes,nMaxNodes;
    int d,i,j;
    int nr,nl;
    int nBucket = 0;
    PARTICLE Temp;
    int iNew,iTemp;

    /*
    ** Allocate and initialize the temporary particles.
    */
    p = malloc(pkd->nLocal*sizeof(PLITE));
    assert(p != NULL);
    for (i=0;i<pkd->nLocal;++i) {
	for (j=0;j<3;++j) p[i].r[j] = pkd->pStore[i].r[j];
	p[i].i = i;
	p[i].iActive = pkd->pStore[i].iActive;
	}
    /*
    ** Allocate stack!
    */
    ns = floor(log(((double)(pkd->nLocal+1))/(M+1))/log(2.0));
    if (ns < 1) ns = 1;	/* want to allocate something! */	
    s = 0;
    S = malloc(ns*sizeof(int));
    assert(S != NULL);
    /*
    ** Allocate initial node storage.
    */
    nMaxNodes = 10000;
    pkd->kdTemp = malloc(nMaxNodes*sizeof(KDT));
    assert(pkd->kdTemp != NULL);
    /*
    ** Make sure we don't have buckets which are larger than the 
    ** pointer arrays for actives and inactives!
    */
    if (M > pkd->nMaxBucketActive) {
	pkd->nMaxBucketActive = M;
	pkd->piActive = realloc(pkd->piActive,pkd->nMaxBucketActive*sizeof(PARTICLE *));
	mdlassert(pkd->mdl,pkd->piActive != NULL);
	pkd->piInactive = realloc(pkd->piInactive,pkd->nMaxBucketActive*sizeof(PARTICLE *));
	mdlassert(pkd->mdl,pkd->piInactive != NULL);
	}
    /*
    ** If we want to split out very active particles from the tree
    ** we do it here, collecting them in Node index 0 as a bucket
    ** with possibly many particles.
    */
    pkd->nVeryActive = 0;
    if (bExcludeVeryActive) {
	/*
	** Now start the partitioning of the particles.
	*/
	i = 0;
	j = pkd->nLocal - 1;
	while (i <= j) {
	    if (p[i].iActive & TYPE_VERYACTIVE) ++i;
	    else break;
	    }
	while (i <= j) {
	    if (!(p[j].iActive & TYPE_VERYACTIVE)) --j;
	    else break;
	    }
	if (i < j) {
	    t = p[i];
	    p[i] = p[j];
	    p[j] = t;
	    while (1) {
		while (p[++i].iActive & TYPE_VERYACTIVE);
		while (!(p[--j].iActive & TYPE_VERYACTIVE));
		if (i < j) {
		    t = p[i];
		    p[i] = p[j];
		    p[j] = t;
		    }
		else break;
		}
	    }
	pkd->nVeryActive = i;
	pkd->kdTemp[0].pLower = 0;
	pkd->kdTemp[0].pUpper = i - 1;
	pkd->kdTemp[0].iLower = 0;
	pkd->kdTemp[0].iParent = 0;
	pkd->kdTemp[0].nGas = 0;	    
	if (pkd->nVeryActive > pkd->nMaxBucketActive) {
	    pkd->nMaxBucketActive = pkd->nVeryActive;
	    pkd->piActive = realloc(pkd->piActive,pkd->nMaxBucketActive*sizeof(PARTICLE *));
	    mdlassert(pkd->mdl,pkd->piActive != NULL);
	    pkd->piInactive = realloc(pkd->piInactive,pkd->nMaxBucketActive*sizeof(PARTICLE *));
	    mdlassert(pkd->mdl,pkd->piInactive != NULL);
	    }
#ifdef GASOLINE
	/*
	** No support for gasoline here quite yet.
	*/
	assert(0);
#endif
	}
    /*
    ** First set up the root node.
    ** Allocate it, set it's bounds and pointers.
    */
    nNodes = 1;
    iNode = nNodes++;
    pkd->kdTemp[iNode].iLower = 0;
    pkd->kdTemp[iNode].iParent = 0;
    pkd->kdTemp[iNode].pLower = pkd->nVeryActive;
    pkd->kdTemp[iNode].pUpper = pkd->nLocal-1;
    i = pkd->kdTemp[iNode].pLower;
    for (j=0;j<3;++j) {
	fMin[j] = p[i].r[j];
	fMax[j] = p[i].r[j];
	}
    for (++i;i<=pkd->kdTemp[iNode].pUpper;++i) {
	for (j=0;j<3;++j) {
	    if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
	    else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
	    }
	}
    for (j=0;j<3;++j) {
	pkd->kdTemp[iNode].bnd.fCenter[j] = 0.5*(fMax[j] + fMin[j]);
	pkd->kdTemp[iNode].bnd.fMax[j] = 0.5*(fMax[j] - fMin[j]);
	}
    if (pkd->kdTemp[iNode].pUpper - pkd->kdTemp[iNode].pLower + 1 <= M) 
	goto DonePart;
    while (1) {
	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/
	pbnd = &pkd->kdTemp[iNode].bnd;
	if (pbnd->fMax[0] < pbnd->fMax[1]) {
	    if (pbnd->fMax[1] < pbnd->fMax[2]) d = 2;
	    else d = 1;
	    }
	else if (pbnd->fMax[0] < pbnd->fMax[2]) d = 2;
	else d = 0;
	fSplit = pbnd->fCenter[d];
	/*
	** Now start the partitioning of the particles about
	** fSplit on dimension given by d.
	*/
	i = pkd->kdTemp[iNode].pLower;
	j = pkd->kdTemp[iNode].pUpper;
	while (i <= j) {
	    if (p[i].r[d] < fSplit) ++i;
	    else break;
	    }
	while (i <= j) {
	    if (fSplit < p[j].r[d]) --j;
	    else break;
	    }
	if (i < j) {
	    t = p[i];
	    p[i] = p[j];
	    p[j] = t;
	    while (1) {
		while (p[++i].r[d] < fSplit);
		while (fSplit < p[--j].r[d]);
		if (i < j) {
		    t = p[i];
		    p[i] = p[j];
		    p[j] = t;
		    }
		else break;
		}
	    }
	nl = i - pkd->kdTemp[iNode].pLower;
	nr = pkd->kdTemp[iNode].pUpper - i + 1;
	if (nl > 0 && nr > 0) {
	    /*
	    ** Allocate 2 new tree nodes making sure that we have
	    ** allocated enough storage.
	    */
	    if (nNodes+2 > nMaxNodes) {
		nMaxNodes += 10000;
		/*
		** Allocate two extra locations to cover the next
		** two nodes we will need.
		*/
		pkd->kdTemp = realloc(pkd->kdTemp,nMaxNodes*sizeof(KDT));
		assert(pkd->kdTemp != NULL);
		}
	    iLeft = nNodes++;
	    pkd->kdTemp[iLeft].iParent = iNode;
	    pkd->kdTemp[iLeft].pLower = pkd->kdTemp[iNode].pLower;
	    pkd->kdTemp[iLeft].pUpper = i-1;
	    iRight = nNodes++;
	    assert(iRight & 1);
	    pkd->kdTemp[iRight].iParent = iNode;
	    pkd->kdTemp[iRight].pLower = i;
	    pkd->kdTemp[iRight].pUpper = pkd->kdTemp[iNode].pUpper;
	    pkd->kdTemp[iNode].iLower = iLeft;
	    if (bSqueeze) {
		/*
		** Calculate bounds for these 2 new cells.
		*/
		i = pkd->kdTemp[iLeft].pLower;
		for (j=0;j<3;++j) {
		    fMin[j] = p[i].r[j];
		    fMax[j] = p[i].r[j];
		    }
		for (++i;i<=pkd->kdTemp[iLeft].pUpper;++i) {
		    for (j=0;j<3;++j) {
			if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
			else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
			}
		    }
		for (j=0;j<3;++j) {
		    pkd->kdTemp[iLeft].bnd.fCenter[j] = 0.5*(fMax[j] + fMin[j]);
		    pkd->kdTemp[iLeft].bnd.fMax[j] = 0.5*(fMax[j] - fMin[j]);
		    }
		i = pkd->kdTemp[iRight].pLower;
		for (j=0;j<3;++j) {
		    fMin[j] = p[i].r[j];
		    fMax[j] = p[i].r[j];
		    }
		for (++i;i<=pkd->kdTemp[iRight].pUpper;++i) {
		    for (j=0;j<3;++j) {
			if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
			else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
			}
		    }
		for (j=0;j<3;++j) {
		    pkd->kdTemp[iRight].bnd.fCenter[j] = 0.5*(fMax[j] + fMin[j]);
		    pkd->kdTemp[iRight].bnd.fMax[j] = 0.5*(fMax[j] - fMin[j]);
		    }
		}
	    else {
		pkd->kdTemp[iLeft].bnd = pkd->kdTemp[iNode].bnd;
		pkd->kdTemp[iLeft].bnd.fMax[d] *= 0.5;
		pkd->kdTemp[iRight].bnd = pkd->kdTemp[iLeft].bnd;
		pkd->kdTemp[iLeft].bnd.fCenter[d] -= pkd->kdTemp[iLeft].bnd.fMax[d];
		pkd->kdTemp[iRight].bnd.fCenter[d] += pkd->kdTemp[iRight].bnd.fMax[d];
		}
	    }
	else {
	    /*
	    ** No nodes allocated, this can't happen if we
	    ** are squeezing (or shouldn't).
	    ** Change the bounds!
	    */
	    assert(!bSqueeze);
	    pbnd->fMax[d] *= 0.5;
	    if (nl > 0) {
		pbnd->fCenter[d] -= pbnd->fMax[d];
		iLeft = iNode;
		}
	    else {
		pbnd->fCenter[d] += pbnd->fMax[d];
		iRight = iNode;
		}
	    }
	/*
	** Now figure out which subfile to process next.
	*/
	if (nl > M && nr > M) {
	    if (nr > nl) {
		S[s++] = iRight;	/* push tr */
		iNode = iLeft;		/* process lower subfile */
		}
	    else {
		S[s++] = iLeft;	/* push tl */
		iNode = iRight;		/* process upper subfile */
		}
	    }
	else {
	    if (nl > M) {
		iNode = iLeft;		/* process lower subfile */
		}
	    else if (nl > 0) {
		pkd->kdTemp[iLeft].iLower = 0;
		++nBucket;
#ifdef GASOLINE
		/*
		** Now start the partitioning of the particles into (Gas,Other) order.
		*/
		i = pkd->kdTemp[iLeft].pLower;
		j = pkd->kdTemp[iLeft].pUpper;
		while (i <= j) {
		    if (p[i].iActive & TYPE_GAS) ++i;
		    else break;
		    }
		while (i <= j) {
		    if (!(p[j].iActive & TYPE_GAS)) --j;
		    else break;
		    }
		if (i < j) {
		    t = p[i];
		    p[i] = p[j];
		    p[j] = t;
		    while (1) {
			while (p[++i].iActive & TYPE_GAS);
			while (!(p[--j].iActive & TYPE_GAS));
			if (i < j) {
			    t = p[i];
			    p[i] = p[j];
			    p[j] = t;
			    }
			else break;
			}
		    }
		pkd->kdTemp[iLeft].nGas = i - pkd->kdTemp[iLeft].pLower;
#endif
		}
	    if (nr > M) {
		iNode = iRight;		/* process upper subfile */
		}
	    else if (nr > 0) {
		pkd->kdTemp[iRight].iLower = 0;
		++nBucket;
#ifdef GASOLINE
		/*
		** Now start the partitioning of the particles into (Gas,Other) order.
		*/
		i = pkd->kdTemp[iRight].pLower;
		j = pkd->kdTemp[iRight].pUpper;
		while (i <= j) {
		    if (p[i].iActive & TYPE_GAS) ++i;
		    else break;
		    }
		while (i <= j) {
		    if (!(p[j].iActive & TYPE_GAS)) --j;
		    else break;
		    }
		if (i < j) {
		    t = p[i];
		    p[i] = p[j];
		    p[j] = t;
		    while (1) {
			while (p[++i].iActive & TYPE_GAS);
			while (!(p[--j].iActive & TYPE_GAS));
			if (i < j) {
			    t = p[i];
			    p[i] = p[j];
			    p[j] = t;
			    }
			else break;
			}
		    }
		pkd->kdTemp[iRight].nGas = i - pkd->kdTemp[iRight].pLower;
#endif
		}
	    }
	if (nl <= M && nr <= M) {
	    if (s) iNode = S[--s];		/* pop tn */
	    else break;
	    }
	}
    DonePart:
    /*
    ** Now we move the particles in one go using the temporary
    ** particles which have been shuffled.
    */
    iTemp = 0;
    while (1) {
	Temp = pkd->pStore[iTemp];
	i = iTemp;
	iNew = p[i].i;
	while (iNew != iTemp) {
	    pkd->pStore[i] = pkd->pStore[iNew];
	    p[i].i = 0;
	    i = iNew;
	    iNew = p[i].i;
	    }
	pkd->pStore[i] = Temp;
	p[i].i = 0;
	while (!p[iTemp].i) {
	    if (++iTemp == pkd->nLocal) goto Done;
	    }
	}
    Done:
    free(p);
    free(S);
    /* printf("nBucket:%d AvgM:%g\n",nBucket,pkd->nLocal/(double)nBucket); */
    return nNodes;
    }


void Create(PKD pkd,FLOAT diCrit2,int bTempBound)
    {
    PARTICLE *p = pkd->pStore;
    KDT *t = pkd->kdTemp;
    KDN *c = pkd->kdNodes;
    KDN *pkdn,*pkdl,*pkdu;
    MOMR mom;
    FLOAT m,fMass,x,y,z,vx,vy,vz,ft,d2,d2Max,dih2;
    int iNode,pj,d,nDepth;

    nDepth = 1;
    pkd->nMaxDepth = 1;
    if (pkd->nVeryActive > 0) {
	/*
	** Copy bucket-0 which possibly contains a number of very active particles.
	** We don't need any moment or bounds information for the very active particles
	** (this will change on a much shorter timescale anyway) so we just keep the 
	** particle pointers in this special bucket-0.
	*/
	iNode = 0;
	pkdn = &c[iNode];
	pkdn->iLower = t[iNode].iLower;
	pkdn->iParent = t[iNode].iParent;
	pkdn->pLower = t[iNode].pLower;
	pkdn->pUpper = t[iNode].pUpper;
	}
    /*
    ** Create the root node of the real tree!
    */
    iNode = ROOT;
    while (1) {
	while (t[iNode].iLower) {
	    iNode = t[iNode].iLower;
	    ++nDepth;
	    /*
	    ** Is this the deepest in the tree so far?
	    */
	    if (nDepth > pkd->nMaxDepth) pkd->nMaxDepth = nDepth;
	    }
	/*
	** Now calculate all bucket quantities!
	** This includes M,CoM,Moments and special
	** bounds and iMaxRung.
	*/
	pkdn = &c[iNode];
	pkdn->iLower = t[iNode].iLower;
	pkdn->iParent = t[iNode].iParent;
	pkdn->pLower = t[iNode].pLower;
	pkdn->pUpper = t[iNode].pUpper;
	if (bTempBound) {
	    pkdn->bnd = t[iNode].bnd;
	    }
	else {
	    pj = pkdn->pLower;
	    for (d=0;d<3;++d) {
		ft = p[pj].r[d];
		pkdn->bnd.fCenter[d] = ft;
		pkdn->bnd.fMax[d] = ft;
		}
	    for (++pj;pj<=pkdn->pUpper;++pj) {
		for (d=0;d<3;++d) {
		    ft = p[pj].r[d];
		    if (ft < pkdn->bnd.fCenter[d])
			pkdn->bnd.fCenter[d] = ft;
		    else if (ft > pkdn->bnd.fMax[d])
			pkdn->bnd.fMax[d] = ft;
		    }
		}
	    for (d=0;d<3;++d) {
		ft = pkdn->bnd.fCenter[d];
		pkdn->bnd.fCenter[d] = 0.5*(pkdn->bnd.fMax[d] + ft);
		pkdn->bnd.fMax[d] = 0.5*(pkdn->bnd.fMax[d] - ft);
		}			
	    }
	pj = pkdn->pLower;
	p[pj].iBucket = iNode;
	m = p[pj].fMass;
	fMass = m;
	dih2 = m/(p[pj].fSoft*p[pj].fSoft);
	x = m*p[pj].r[0];
	y = m*p[pj].r[1];
	z = m*p[pj].r[2];
	vx = m*p[pj].v[0];
	vy = m*p[pj].v[1];
	vz = m*p[pj].v[2];
	pkdn->iActive = p[pj].iActive;
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p[pj].iBucket = iNode;
	    m = p[pj].fMass;
	    fMass += m;
	    dih2 += m/(p[pj].fSoft*p[pj].fSoft);
	    x += m*p[pj].r[0];
	    y += m*p[pj].r[1];
	    z += m*p[pj].r[2];		
	    vx += m*p[pj].v[0];
	    vy += m*p[pj].v[1];
	    vz += m*p[pj].v[2];		
	    pkdn->iActive |= p[pj].iActive;
	    }
	m = 1/fMass;
	pkdn->r[0] = m*x;
	pkdn->r[1] = m*y;
	pkdn->r[2] = m*z;
	pkdn->v[0] = m*vx;
	pkdn->v[1] = m*vy;
	pkdn->v[2] = m*vz;
	dih2 *= m;
	pkdn->fSoft2 = 1/dih2;
	/*
	** Now calculate the reduced multipole moment.
	*/
	pj = pkdn->pLower;
	x = p[pj].r[0] - pkdn->r[0];
	y = p[pj].r[1] - pkdn->r[1];
	z = p[pj].r[2] - pkdn->r[2];
	d2Max = momMakeMomr(&pkdn->mom,p[pj].fMass,x,y,z);
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    x = p[pj].r[0] - pkdn->r[0];
	    y = p[pj].r[1] - pkdn->r[1];
	    z = p[pj].r[2] - pkdn->r[2];
	    d2 = momMakeMomr(&mom,p[pj].fMass,x,y,z);
	    momAddMomr(&pkdn->mom,&mom);
	    /*
	    ** Update bounding ball and softened bounding ball.
	    */
	    d2Max = (d2 > d2Max)?d2:d2Max;
	    }
	/*
	** Now determine the opening radius for gravity.
	*/
	d2Max *= FOPEN_FACTOR*diCrit2;
	pkdn->fOpen2 = d2Max;
#ifdef GASOLINE
	pkdn->nGas = t[iNode].nGas;
#endif
	/*
	** Finished with the bucket, move onto the next one,
	** or to the parent.
	*/
	while (iNode & 1) {
	    iNode = c[iNode].iParent;
	    --nDepth;
	    if (!iNode) {
		assert(nDepth == 0);
		return;	/* exit point!!! */
		}
	    /*
	    ** Now combine quantities from each of the children (2) of
	    ** this cell to form the quantities for this cell.
	    ** First find the CoM, just like for the bucket.
	    */
	    c[iNode].iLower = t[iNode].iLower;
	    c[iNode].iParent = t[iNode].iParent;
	    c[iNode].pLower = t[iNode].pLower;
	    c[iNode].pUpper = t[iNode].pUpper;
	    pkdn = &c[iNode];
	    pkdl = &c[pkdn->iLower];
	    pkdu = &c[pkdn->iLower + 1];
	    if (bTempBound) {
		pkdn->bnd = t[iNode].bnd;
		pkdCombineCells(pkdn,pkdl,pkdu,0);
		}
	    else {
		pkdCombineCells(pkdn,pkdl,pkdu,1);
		}
	    pj = pkdn->pLower;
	    if (pkdn->pUpper - pj < NMAX_OPENCALC) {
		x = p[pj].r[0] - pkdn->r[0];
		y = p[pj].r[1] - pkdn->r[1];
		z = p[pj].r[2] - pkdn->r[2];
		d2Max = x*x + y*y + z*z;
		for (++pj;pj<=pkdn->pUpper;++pj) {
		    x = p[pj].r[0] - pkdn->r[0];
		    y = p[pj].r[1] - pkdn->r[1];
		    z = p[pj].r[2] - pkdn->r[2];
		    d2 = x*x + y*y + z*z;
		    d2Max = (d2 > d2Max)?d2:d2Max;
		    }
		/*
		** Now determine the opening radius for gravity.
		*/
		d2Max *= FOPEN_FACTOR*diCrit2;
		pkdn->fOpen2 = d2Max;
		}
	    else {
		CALCOPEN(pkdn,diCrit2);
		}
	    }
	++iNode;
	}
    }



void pkdCombineCells(KDN *pkdn,KDN *p1,KDN *p2,int bCombineBound)
    {
    MOMR mom;
    FLOAT m1,m2,x,y,z,ifMass;

    m1 = p1->mom.m;
    m2 = p2->mom.m;
    ifMass = 1/(m1 + m2);
    pkdn->r[0] = ifMass*(m1*p1->r[0] + m2*p2->r[0]);
    pkdn->r[1] = ifMass*(m1*p1->r[1] + m2*p2->r[1]);
    pkdn->r[2] = ifMass*(m1*p1->r[2] + m2*p2->r[2]);
    pkdn->v[0] = ifMass*(m1*p1->v[0] + m2*p2->v[0]);
    pkdn->v[1] = ifMass*(m1*p1->v[1] + m2*p2->v[1]);
    pkdn->v[2] = ifMass*(m1*p1->v[2] + m2*p2->v[2]);
    pkdn->fSoft2 = 1.0/(ifMass*(m1/p1->fSoft2 + m2/p2->fSoft2));
    pkdn->iActive = (p1->iActive | p2->iActive);
    /*
    ** Now calculate the reduced multipole moment.
    ** Shift the multipoles of each of the children
    ** to the CoM of this cell and add them up.
    */
    pkdn->mom = p1->mom;
    x = p1->r[0] - pkdn->r[0];
    y = p1->r[1] - pkdn->r[1];
    z = p1->r[2] - pkdn->r[2];
    momShiftMomr(&pkdn->mom,x,y,z);
    mom = p2->mom;
    x = p2->r[0] - pkdn->r[0];
    y = p2->r[1] - pkdn->r[1];
    z = p2->r[2] - pkdn->r[2];
    momShiftMomr(&mom,x,y,z);
    momAddMomr(&pkdn->mom,&mom);
    /*
    ** Combine the bounds!
    */
    if (bCombineBound) {
	BND_COMBINE(pkdn->bnd,p1->bnd,p2->bnd);
	}
#ifdef GASOLINE
    BND_COMBINE(pkdn->bndBall,p1->bndBall,p2->bndBall);
    pkdn->nGas = p1->nGas + p2->nGas;
#endif
    }


void pkdTreeBuild(PKD pkd,int nBucket,FLOAT diCrit2,KDN *pkdn, int bSqueeze, int bExcludeVeryActive)
    {
    int nCells;

    if (pkd->kdNodes) {
	/*
	** Close cell caching space and free up nodes.
	*/
	mdlFinishCache(pkd->mdl,CID_CELL);
	mdlFree(pkd->mdl,pkd->kdNodes);
	}

    pkdClearTimer(pkd,0);
    pkdStartTimer(pkd,0);
    nCells = BuildTemp(pkd,nBucket,bSqueeze,bExcludeVeryActive);
    pkdStopTimer(pkd,0);
    /* 
       printf("Temp Tree Build wallclock: %g secs\n",
       pkdGetWallClockTimer(pkd,0));
       printf("Number of Cells: %d\n",nCells);
    */
    /*
    ** Now allocate the cell storage using mdlMalloc!
    ** Allocate one extra because cell-id 0 is a dummy
    ** cell.
    */
    pkd->nNodes = nCells + 1;
    pkd->kdNodes = mdlMalloc(pkd->mdl,pkd->nNodes*sizeof(KDN));
    assert(pkd->kdNodes != NULL);
    /*
    ** Now create the real tree from the temporary tree.
    ** Last argument is fBallChange, need to set this from somewhere.
    */
    pkdClearTimer(pkd,0);
    pkdStartTimer(pkd,0);
    Create(pkd,diCrit2,bSqueeze);
    pkdStopTimer(pkd,0);
    /*
      printf("Create Tree wallclock: %g secs\n",
      pkdGetWallClockTimer(pkd,0));
      printf("nMaxDepth:%d\n",pkd->nMaxDepth);
    */
    /*
    ** Free up the temporary tree, why not!
    */
    free(pkd->kdTemp);
    pkd->kdTemp = NULL;
    /*
    ** Finally activate a read only cache for remote access.
    */
    mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkd->nNodes);
    /*
    ** Copy the root node for the top-tree construction.
    */
    *pkdn = pkd->kdNodes[ROOT];
    }


void pkdDistribCells(PKD pkd,int nCell,KDN *pkdn)
    {
    KDN *c;
    int i;

    if (pkd->kdTop != NULL) free(pkd->kdTop);
    pkd->kdTop = malloc(nCell*sizeof(KDN));
    assert(pkd->kdTop != NULL);
    c = pkd->kdTop;
    for (i=1;i<nCell;++i) {
	if (pkdn[i].pUpper) {
	    c[i] = pkdn[i];
	    if (pkdn[i].pLower == pkd->idSelf) pkd->iTopRoot = i;
	    }
	}
    }


/*
** Hopefully we can bypass this step once we figure out how to do the
** Multipole Ewald with reduced multipoles.
*/
void pkdCalcRoot(PKD pkd,MOMC *pmom)
    {
    PARTICLE *p = pkd->pStore;
    FLOAT xr = pkd->kdTop[ROOT].r[0];
    FLOAT yr = pkd->kdTop[ROOT].r[1];
    FLOAT zr = pkd->kdTop[ROOT].r[2];
    FLOAT x,y,z;
    MOMC mc;
    int i = 0;

    x = p[i].r[0] - xr;
    y = p[i].r[1] - yr;
    z = p[i].r[2] - zr;
    momMakeMomc(pmom,p[i].fMass,x,y,z);
    for (++i;i<pkd->nLocal;++i) {
	x = p[i].r[0] - xr;
	y = p[i].r[1] - yr;
	z = p[i].r[2] - zr;
	momMakeMomc(&mc,p[i].fMass,x,y,z);
	momAddMomc(pmom,&mc);
	}
    }


void pkdDistribRoot(PKD pkd,MOMC *pmom)
    {
    pkd->momRoot = *pmom;
    }


#ifdef GASOLINE
/*
** Currently we only allow the calculation of bounds-of-balls only in 
** Gasoline, although it could also be needed for collision detection.
** This is mainly in an attempt to reduce the gravity only part of the 
** code to the bare minimum.
*/
void pkdCalcBoundBall(PKD pkd,double fBallFactor,BND *pbnd)
    {
    PARTICLE *p = pkd->pStore;
    KDN *c = pkd->kdNodes;
    KDN *pkdn,*pkdl,*pkdu;
    FLOAT ft,fBall;
    int iNode,pj,d;

    iNode = ROOT;
    while (1) {
	while (c[iNode].iLower) iNode = c[iNode].iLower;
	/*
	** Now calculate all bucket quantities!
	*/
	pkdn = &c[iNode];
	pj = pkdn->pLower;
	fBall = fBallFactor*p[pj].fBall;
	for (d=0;d<3;++d) {
	    ft = p[pj].r[d];
	    pkdn->bndBall.fCenter[d] = ft - fBall;
	    pkdn->bndBall.fMax[d] = ft + fBall;
	    }
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    fBall = fBallFactor*p[pj].fBall;
	    for (d=0;d<3;++d) {
		ft = p[pj].r[d];
		if (ft - fBall < pkdn->bndBall.fCenter[d])
		    pkdn->bndBall.fCenter[d] = ft - fBall;
		if (ft + fBall > pkdn->bndBall.fMax[d])
		    pkdn->bndBall.fMax[d] = ft + fBall;
		}
	    }
	for (d=0;d<3;++d) {
	    ft = pkdn->bndBall.fCenter[d];
	    pkdn->bndBall.fCenter[d] = 0.5*(pkdn->bndBall.fMax[d] + ft);
	    pkdn->bndBall.fMax[d] = 0.5*(pkdn->bndBall.fMax[d] - ft);
	    }
	/*
	** Finished with the bucket, move onto the next one,
	** or to the parent.
	*/
	while (iNode & 1) {
	    iNode = c[iNode].iParent;
	    if (!iNode) {
		*pbnd = pkd->kdNodes[ROOT].bndBall;
		return;	/* exit point!!! */
		}
	    /*
	    ** Now combine quantities from each of the children (2) of
	    ** this cell to form the quantities for this cell.
	    ** First find the CoM, just like for the bucket.
	    */
	    pkdn = &c[iNode];
	    pkdl = &c[pkdn->iLower];
	    pkdu = &c[pkdn->iLower + 1];
	    BND_COMBINE(pkdn->bndBall,pkdl->bndBall,pkdu->bndBall);
	    }
	++iNode;
	}
    }


void pkdDistribBoundBall(PKD pkd,int nCell,BND *pbnd)
    {
    KDN *c;
    int i;

    assert(pkd->kdTop != NULL);
    c = pkd->kdTop;
    for (i=1;i<nCell;++i) {
	if (pbnd[i].fMax[0] > 0) {
	    c[i].bndBall = pbnd[i];
	    }
	}
    }


#endif /* of GASOLINE */
