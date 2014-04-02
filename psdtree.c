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
#include "psdtree.h"
#include "moments.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "smooth.h"
#include "pst.h"
#include "qsort.h"

#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#define X 0
#define V 1

void psdInitializeParticles(PKD pkd, BND *prbnd, BND *pvbnd) {
    KDN *pNode;
    int j;
    BND *rbnd, *vbnd;

    /*
    **It is only forseen that there are 4 reserved nodes at present 0-NULL, 1-ROOT, 2-UNUSED, 3-VAROOT.
    */
    pkd->nNodes = NRESERVED_NODES;

    pkd->nVeryActive = 0;

    /*
    ** Set up the root node.
    */
    pNode = pkdTreeNode(pkd,ROOT);
    pNode->iLower = 0;
    pNode->iParent = 0;
    pNode->pLower = 0;
    pNode->pUpper = pkd->nLocal - 1;

    rbnd = pkdNodeBnd(pkd,pNode);
    for (j=0;j<3;++j) {
	rbnd->fCenter[j] = prbnd->fCenter[j];
	rbnd->fMax[j]    = prbnd->fMax[j];
    }
    vbnd = pkdNodeVBnd(pkd,pNode);
    for (j=0;j<3;++j) {
	vbnd->fCenter[j] = pvbnd->fCenter[j];
	vbnd->fMax[j]    = pvbnd->fMax[j];
	}
    }

int psdTreeDepth(PKD pkd) {
    int iNode = ROOT;
    int nDepth = 1;
    int maxDepth = 1;
    while (1) {
	while (pkdTreeNode(pkd,iNode)->iLower) {
	    iNode = pkdTreeNode(pkd,iNode)->iLower;
	    ++nDepth;
	    if (nDepth > maxDepth) maxDepth = nDepth;
	    }
	while (iNode & 1) {
	    iNode = pkdTreeNode(pkd,iNode)->iParent;
	    --nDepth;
	    if (!iNode) {
		assert(nDepth == 0);
		return maxDepth;
		}
	    }
	iNode++;
	}

    return -1;
    }

struct KEY {
    int64_t keyx;
    int64_t keyv;
    float   fMass;
};

int keyx_compar(const void *a0, const void *b0) {
    struct KEY *a = (struct KEY *)a0;
    struct KEY *b = (struct KEY *)b0;
    int64_t ax = a->keyx;
    int64_t bx = b->keyx;
    if (ax > bx) return +1;
    if (ax < bx) return -1;
    return 0;
    }

int keyv_compar(const void *a0, const void *b0) {
    struct KEY *a = (struct KEY *)a0;
    struct KEY *b = (struct KEY *)b0;
    int64_t av = a->keyv;
    int64_t bv = b->keyv;
    if (av > bv) return +1;
    if (av < bv) return -1;
    return 0;
    }

/*
** Compute entropy. E = -\sum (p ln(p))
*/
static float Ex(struct KEY *keys, int N, float M) {
    int i;
    float s = 0;
    int64_t key = keys[0].keyx;
    float   p   = keys[0].fMass;
    for (i=1; i < N; i++) {
	if (keys[i].keyx != key) {
	    p /= M;
	    s -= p * logf(p);
	    p = 0;
	    key = keys[i].keyx;
	    }

	p += keys[i].fMass;
	}

    if (p != 0)
    {
	p /= M;
	s -= p * logf(p);
    }

    return s;
}

static float Ev(struct KEY *keys, int N, float M) {
    int i;
    float s = 0;
    int64_t key = keys[0].keyv;
    float   p   = keys[0].fMass;
    assert(M > 0);
    for (i=1; i < N; i++) {
	if (keys[i].keyv != key) {
	    p /= M;
	    s -= p * logf(p);
	    p = 0;
	    key = keys[i].keyv;
	    }

	p += keys[i].fMass;
	}

    if (p != 0)
    {
	p /= M;
	s -= p * logf(p);
    }

    return s;
}

void inline SAVE_BOUNDS(PKD pkd, PSX smx, KDN *pNode, BND **bnd)
{
    int j;
    PSMETRIC *b = smx->psm + pNode->pLower;

#if 0
    for (j=0;j<3;++j) assert(!isinf(bnd[0]->fMax[j]));
    for (j=0;j<3;++j) assert(!isinf(bnd[1]->fMax[j]));
    for (j=0;j<3;++j) assert(!isinf(bnd[0]->fCenter[j]));
    for (j=0;j<3;++j) assert(!isinf(bnd[1]->fCenter[j]));

    for (j=0;j<3;++j) assert(fabs(bnd[0]->fMax[j]) < 1e10);
    for (j=0;j<3;++j) assert(fabs(bnd[1]->fMax[j]) < 1e10);
    for (j=0;j<3;++j) assert(fabs(bnd[0]->fCenter[j]) < 1e10);
    for (j=0;j<3;++j) assert(fabs(bnd[1]->fCenter[j]) < 1e10);
#endif

    for (j=0;j<3;++j)
    {
	b->rscale[j] = 0;
	if (bnd[X]->fMax[j] != 0) 
	{
#ifdef NO_PSMETRIC
	    b->rscale[j] = 1.;
#else
	    b->rscale[j] = 1. / bnd[X]->fMax[j]; 
#endif
	}
    }

    for (j=0;j<3;++j)
    {
	b->vscale[j] = 0;
	if (bnd[V]->fMax[j] != 0) 
	{
#ifdef NO_PSMETRIC
	    b->vscale[j] = 0.;
#else
	    b->vscale[j] = 1. / bnd[V]->fMax[j];
#endif
	}
    }
}

void entropy(PKD pkd, float *e, KDN *pNode, BND **bnd, struct KEY *keys) {
    int s, i,j;
    int64_t bx, by, bz;
    float Mtotal = 0;
    int Npart = pNode->pUpper - pNode->pLower + 1;
    int64_t Nb;
    Nb = 1+cbrt(Npart/10.);
    int maxNb = 1 << ((sizeof(Nb)*8) / 3);
    double *v;

    if (Nb > maxNb) Nb = maxNb;
    assert(Nb >= 1);

    const FLOAT dXx = bnd[X]->fMax[0] > 0 ? (FLOAT)Nb / (2*bnd[X]->fMax[0]) : 0;
    const FLOAT dXy = bnd[X]->fMax[1] > 0 ? (FLOAT)Nb / (2*bnd[X]->fMax[1]) : 0;
    const FLOAT dXz = bnd[X]->fMax[2] > 0 ? (FLOAT)Nb / (2*bnd[X]->fMax[2]) : 0;
    const FLOAT dVx = bnd[V]->fMax[0] > 0 ? (FLOAT)Nb / (2*bnd[V]->fMax[0]) : 0;
    const FLOAT dVy = bnd[V]->fMax[1] > 0 ? (FLOAT)Nb / (2*bnd[V]->fMax[1]) : 0;
    const FLOAT dVz = bnd[V]->fMax[2] > 0 ? (FLOAT)Nb / (2*bnd[V]->fMax[2]) : 0;

    const FLOAT edgeXx = bnd[X]->fCenter[0] - bnd[X]->fMax[0];
    const FLOAT edgeXy = bnd[X]->fCenter[1] - bnd[X]->fMax[1];
    const FLOAT edgeXz = bnd[X]->fCenter[2] - bnd[X]->fMax[2];
    const FLOAT edgeVx = bnd[V]->fCenter[0] - bnd[V]->fMax[0];
    const FLOAT edgeVy = bnd[V]->fCenter[1] - bnd[V]->fMax[1];
    const FLOAT edgeVz = bnd[V]->fCenter[2] - bnd[V]->fMax[2];

    Mtotal = 0;
    for (j=pNode->pLower; j <= pNode->pUpper; j++) 
    {
	PARTICLE *p = pkdParticle(pkd,j);
	v = pkdVel(pkd, p);
	float mass = pkdMass(pkd,p);
	Mtotal += mass;

	bx  = (int64_t)((p->r[0]-edgeXx) * dXx); bx -= (bx == Nb);
	by  = (int64_t)((p->r[1]-edgeXy) * dXy); by -= (by == Nb);
	bz  = (int64_t)((p->r[2]-edgeXz) * dXz); bz -= (bz == Nb);
	keys[j].keyx = bx + Nb * (by + Nb*bz);

	bx  = (int64_t)((v[0]-edgeVx) * dVx); bx -= (bx == Nb);
	by  = (int64_t)((v[1]-edgeVy) * dVy); by -= (by == Nb);
	bz  = (int64_t)((v[2]-edgeVz) * dVz); bz -= (bz == Nb);

	keys[j].keyv = bx + Nb * (by + Nb*bz);
	keys[j].fMass = mass;
    }

    struct KEY * const pkey = keys + pNode->pLower;
#define cmpx(a,b) ((a)->keyx < (b)->keyx)
    QSORT(struct KEY, pkey, Npart, cmpx)
    e[X] = Ex(pkey, Npart, Mtotal);

#define cmpv(a,b) ((a)->keyv < (b)->keyv)
    QSORT(struct KEY, pkey, Npart, cmpv);
    e[V] = Ev(pkey, Npart, Mtotal);
}

void shrink_wrap(PKD pkd, BND **bnd, int iLower, int iUpper)
{
    int pj;
    int d;
    double ft;
    PARTICLE *p;
    double *v;
    int subspace;

    for (subspace = X; subspace <= V; subspace++)
    {
	pj = iLower;
	p = pkdParticle(pkd,pj);
	v = pkdVel(pkd, p);
	for (d=0;d<3;++d) {
	    ft = (subspace == X) ? p->r[d] : v[d];
	    bnd[subspace]->fCenter[d] = ft;
	    bnd[subspace]->fMax[d] = ft;
	    }
	for (++pj;pj<=iUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    v = pkdVel(pkd, p);
	    for (d=0;d<3;++d) {
		ft = (subspace == X) ? p->r[d] : v[d];
		if (ft < bnd[subspace]->fCenter[d])
		    bnd[subspace]->fCenter[d] = ft;
		else if (ft > bnd[subspace]->fMax[d])
		    bnd[subspace]->fMax[d] = ft;
		}
	    }
	for (d=0;d<3;++d) {
	    ft = bnd[subspace]->fCenter[d];
	    bnd[subspace]->fCenter[d] = 0.5*(bnd[subspace]->fMax[d] + ft);
	    bnd[subspace]->fMax[d] = 0.5*(bnd[subspace]->fMax[d] - ft);
	    }
    }
}

/*
** Build a tree that will be used to estimate the phase-space metric.
**
** M is the bucket size.
** This function assumes that the root node is correctly set up (particularly the bounds).
*/

#define TEMP_S_INCREASE 100
void BuildPsdTemp(PKD pkd, PSX smx, int iNode,int M, int maxNb) {
    KDN *pLeft, *pRight;
    KDN *pNode;
    BND *bnd[2];
    FLOAT fSplit;
    int *S;		/* this is the stack */
    char *D;	     /* stack for the last cut dimension */
    uint8_t lastd;
    uint8_t last_subspace;

    int iLeft,iRight;
    int s,d,i,j;
    PARTICLE *pi, *pj;
    int subspace;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;
    int backup = 0;
    int Ntotal;

    double *v;

    pNode = pkdTreeNode(pkd,iNode);

    Ntotal = pNode->pUpper - pNode->pLower + 1;

    struct KEY *keys = malloc((Ntotal+1) * sizeof(*keys)); assert(keys != NULL);

    assert(maxNb > 0);

    /*************************************************************************/
    NEW_STACK(S, TEMP_S_INCREASE);
    NEW_STACK(D, TEMP_S_INCREASE);
    /*************************************************************************/


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

    if (pNode->pUpper - pNode->pLower + 1 <= M)
	goto DonePart;


    /*************************************************************************/

    bnd[0] = pkdNodeBnd(pkd, pNode);
    bnd[1] = pkdNodeVBnd(pkd, pNode);

    assert( bnd[0]->fMax[0] > 0.0
	 || bnd[0]->fMax[1] > 0.0
	 || bnd[0]->fMax[2] > 0.0
	 || bnd[1]->fMax[0] > 0.0
	 || bnd[1]->fMax[1] > 0.0
	 || bnd[1]->fMax[2] > 0.0);

    subspace = X;

    assert(STACK_EMPTY(S));
    PUSH(S, iNode); PUSH(D, (subspace<<2) | 0);

    while (!STACK_EMPTY(S)) {
	assert(!STACK_EMPTY(D));
	iNode = POP(S);
	int q = POP(D);
	lastd = q & 0x03;
	last_subspace = (q & ~0x03) >> 2;

	if (iNode == 0) {
	    assert(q == -1);
	    assert(backup != 0);
	    pkd->nNodes = backup;
	    backup = 0;
	    pRight = pkdTreeNode(pkd,pkd->nNodes-1);
	    pLeft  = pkdTreeNode(pkd,pkd->nNodes-2);

	    /* If one of these isn't really a bucket it will 
	    ** be popped of the stack next and iLower will be 
	    ** set correctly.
	    */
	    pLeft->iLower  = 0; 
	    pRight->iLower = 0; 

	    continue;
	    }


	pNode = pkdTreeNode(pkd,iNode);
	bnd[0] = pkdNodeBnd(pkd, pNode);
	bnd[1] = pkdNodeVBnd(pkd, pNode);

	int Npart = pNode->pUpper - pNode->pLower + 1;

#if 1
	if (!( bnd[0]->fMax[0] > 0.0 &&
	       bnd[0]->fMax[1] > 0.0 &&
	       bnd[0]->fMax[2] > 0.0 
#ifndef NO_PSMETRIC
	       &&
	       bnd[1]->fMax[0] > 0.0 &&
	       bnd[1]->fMax[1] > 0.0 &&
	       bnd[1]->fMax[2] > 0.0 
#endif
	       ))
	    fprintf(stderr, "%i %i, %g %g %g %g %g %g\n",
	        iNode, Npart,
		bnd[0]->fMax[0],
		bnd[0]->fMax[1],
		bnd[0]->fMax[2],
		bnd[1]->fMax[0],
		bnd[1]->fMax[1],
		bnd[1]->fMax[2]);
#endif

	assert( bnd[0]->fMax[0] > 0.0 
	     || bnd[0]->fMax[1] > 0.0 
	     || bnd[0]->fMax[2] > 0.0
	     || bnd[1]->fMax[0] > 0.0 
	     || bnd[1]->fMax[1] > 0.0 
	     || bnd[1]->fMax[2] > 0.0 
		);

	assert( !isinf(bnd[0]->fMax[0])
	     || !isinf(bnd[0]->fMax[1])
	     || !isinf(bnd[0]->fMax[2])
	     || !isinf(bnd[1]->fMax[0])
	     || !isinf(bnd[1]->fMax[1])
	     || !isinf(bnd[1]->fMax[2])
		);


    /*************************************************************************/

	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/

#ifdef NO_PSMETRIC
	d = max_bnd_range(bnd[subspace]->fMax, 0, 3);
	assert(d != -1);
#else
	float e[2];


	double fb = 0.5 * pow(Ntotal, 1./6);

	if (Npart > 7)
	{
	    s = X;
	    for (i=0; i < 3; i++)
	    {
		double l_min = bnd[s]->fCenter[i] - bnd[s]->fMax[i];
		double l_max = bnd[s]->fCenter[i] + bnd[s]->fMax[i];

		double x_max = -HUGE_VAL;
		double x_min = HUGE_VAL;
		for (j=pNode->pLower; j <= pNode->pUpper; j++) 
		{
		    PARTICLE *p = pkdParticle(pkd,j);
		    if (p->r[i] > x_max) x_max = p->r[i];
		    if (p->r[i] < x_min) x_min = p->r[i];
		}

		double mean_sep = (x_max - x_min) / (Npart - 1);

		if ((l_max - x_max) > fb * mean_sep
		&&  (x_min - l_min) > fb * mean_sep)
		{
		    bnd[s]->fCenter[i] = (x_max + x_min) / 2;
		    bnd[s]->fMax[i]    = (x_max - x_min) / 2;
		}
	    }

	    s = V;
	    for (i=0; i < 3; i++)
	    {
		double l_min = bnd[s]->fCenter[i] - bnd[s]->fMax[i];
		double l_max = bnd[s]->fCenter[i] + bnd[s]->fMax[i];

		double x_max = -HUGE_VAL;
		double x_min = HUGE_VAL;
		for (j=pNode->pLower; j <= pNode->pUpper; j++) 
		{
		    PARTICLE *p = pkdParticle(pkd,j);
		    v = pkdVel(pkd, p);
		    if (v[i] > x_max) x_max = v[i];
		    if (v[i] < x_min) x_min = v[i];
		}

		double mean_sep = (x_max - x_min) / (Npart - 1);

		if ((l_max - x_max) > fb * mean_sep
		&&  (x_min - l_min) > fb * mean_sep)
		{
		    bnd[s]->fCenter[i] = (x_max + x_min) / 2;
		    bnd[s]->fMax[i] = (x_max - x_min) / 2;
		}
	    }
	}

	e[V] = e[X] = 0;
	if (Npart > 1)
	    entropy(pkd, e, pNode, bnd, keys);
	if (e[V] == e[X])
	{
	    subspace = !subspace;
	}
	else
	{
	    if (e[V] != 0 && e[V] < e[X])
	    {
		subspace = V;
	    }
	    else
	    {
		subspace = X;
	    }
	}

	d = max_bnd_range(bnd[subspace]->fMax, 0, 3);

	if (d == -1)
	{
	    subspace = !subspace;
	    d = max_bnd_range(bnd[subspace]->fMax, 0, 3);
	    assert(d != -1);
	}
#endif

	fSplit = bnd[subspace]->fCenter[d];

	assert(subspace == X || subspace == V);
	assert(0 <= d && d < 3);

    /*************************************************************************/


	/*
	** Now start the partitioning of the particles about
	** fSplit on dimension given by d.
	*/
	i = pNode->pLower;
	j = pNode->pUpper;
	pi = pkdParticle(pkd,i);
	pj = pkdParticle(pkd,j);
	if (subspace == X)
	{
	    PARTITION(i<j,i<=j,
		   pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		   pkdSwapParticle(pkd, pi,pj),
		   pi->r[d] < fSplit, pj->r[d] >= fSplit);
	}
	else
	{
	    PARTITION(i<j,i<=j,
		   pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		   pkdSwapParticle(pkd, pi,pj),
		   pkdVel(pkd, pi)[d] < fSplit, pkdVel(pkd,pj)[d] >= fSplit);
	}

	nl = i - pNode->pLower;
	nr = pNode->pUpper - i + 1;

    /*************************************************************************/

	if (nl > 0 && nr > 0) {
	    BND * rbnd[2], *lbnd[2];
	    /*
	    ** Allocate 2 new tree nodes making sure that we have
	    ** allocated enough storage.
	    */
	    if ( pkd->nNodes+2 > pkd->nMaxNodes ) {
		pkdExtendTree(pkd);
		}
	    iLeft = pkd->nNodes++;
	    pLeft = pkdTreeNode(pkd,iLeft);
	    pLeft->iParent = iNode;
	    pLeft->pLower = pNode->pLower;
	    pLeft->pUpper = i-1;
	    iRight = pkd->nNodes++;
	    pRight = pkdTreeNode(pkd,iRight);
	    assert(iRight & 1);
	    pRight->iParent = iNode;
	    pRight->pLower = i;
	    pRight->pUpper = pNode->pUpper;
	    pNode->iLower = iLeft;


	    rbnd[0] = pkdNodeBnd(pkd, pRight); rbnd[1] = pkdNodeVBnd(pkd, pRight);
	    lbnd[0] = pkdNodeBnd(pkd, pLeft);  lbnd[1] = pkdNodeVBnd(pkd, pLeft);


	    /*
	    ** Now deal with the bounds.
	    **
	    ** ADD SHRINK WRAPPING -- jpc 27.3.2010
	    **
	    */
	    for (s=0; s < 2; s++)
	    {
		for (j=0;j<3;++j) {
		    lbnd[s]->fCenter[j] = rbnd[s]->fCenter[j] = bnd[s]->fCenter[j]; 
		    lbnd[s]->fMax[j]    = rbnd[s]->fMax[j]    = bnd[s]->fMax[j];
		    }
	    }
	    s = subspace;
	    j = d;
	    rbnd[s]->fMax[j]    = lbnd[s]->fMax[j]   = 0.5*bnd[s]->fMax[j];
	    lbnd[s]->fCenter[j] = bnd[s]->fCenter[j] - lbnd[s]->fMax[j];
	    rbnd[s]->fCenter[j] = bnd[s]->fCenter[j] + rbnd[s]->fMax[j];

	    /*
	    ** Now figure out which subfile to process next.
	    */
	    lc = ((nl > M)); /* this condition means the left child is not a bucket */
	    rc = ((nr > M));
	    if (rc && lc) {
		assert(backup == 0);
		EXTEND_STACK(S); /* Allocate more stack if required */
		EXTEND_STACK(D); /* Allocate more stack if required */

		if (nr > nl) {
		    PUSH(S, iRight); PUSH(D, (subspace << 2) | d);
		    PUSH(S, iLeft);  PUSH(D, (subspace << 2) | d);
		    }
		else {
		    PUSH(S, iLeft);  PUSH(D, (subspace << 2) | d);
		    PUSH(S, iRight); PUSH(D, (subspace << 2) | d);
		    }
		}
	    else if (lc) {
		/*
		** Right must be a bucket in this case!
		*/
		assert(backup == 0);
		PUSH(S, iLeft); PUSH(D, (subspace << 2) | d);
		++nBucket;

		if (nr > 1) {
		    backup = pkd->nNodes;
		    PUSH(S, 0); PUSH(D, -1);
		    PUSH(S, iRight); PUSH(D, (subspace << 2) | d);
		    }
		else {
		    pRight->iLower = 0;
		    SAVE_BOUNDS(pkd,smx,pRight,rbnd);
		    }
		}
	    else if (rc) {
		/*
		** Left must be a bucket in this case!
		*/
		assert(backup == 0);
		PUSH(S, iRight); PUSH(D, (subspace << 2) | d);
		++nBucket;

		if (nl > 1) {
		    backup = pkd->nNodes;
		    PUSH(S, 0); PUSH(D, -1);
		    PUSH(S, iLeft); PUSH(D, (subspace << 2) | d);
		    }
		else {
		    pLeft->iLower = 0;
		    SAVE_BOUNDS(pkd,smx,pLeft,lbnd);
		    }
		}
	    else {
		/*
		** Both are buckets.
		*/

		/* 
		** This should only happen once when we want to begin decending into
		** the two buckets. backup==0 ensures this
		*/
		if (!(nr==1 && nl==1) && backup == 0) { 
		    backup = pkd->nNodes;
		    PUSH(S, 0); PUSH(D, -1);
		    ++nBucket;
		    ++nBucket;
		    }

		if (nr > 1) { PUSH(S, iRight); PUSH(D, (subspace << 2) | d); }
		if (nl > 1) { PUSH(S, iLeft);  PUSH(D, (subspace << 2) | d); }

		if (nr == 1) {pRight->iLower=0; SAVE_BOUNDS(pkd,smx,pRight,rbnd); }
		if (nl == 1) {pLeft->iLower=0;  SAVE_BOUNDS(pkd,smx,pLeft,lbnd); }
		}
	    }
	else {
	    /*
	    ** No nodes allocated, Change the bounds if needed!
	    */

	    assert((nr == 0) ^ (nl == 0));
	    int n = nr + nl; /* one of these will be zero */
	    lc = rc = 0;

	    bnd[subspace]->fMax[d] *= 0.5;
	    if (bnd[subspace]->fMax[d] <= 0)
	    {
		fprintf(stderr, "subspace %i  d %i\n", subspace, d);
		fprintf(stderr, "nr %i  nl %i  fMax[%i] %f\n", nr, nl, d, bnd[subspace]->fMax[d]);
		fprintf(stderr, "Node %i\n", pkd->idSelf);
		for (i=pNode->pLower; i <=pNode->pUpper; i++)
		{
		    PARTICLE *p = pkdParticle(pkd, i);
		    double *v = pkdVel(pkd, p);
		    fprintf(stderr, "%10ld] %e %e %e  %e %e %e\n", 
			(uint64_t)p->iOrder,
			p->r[0], p->r[1], p->r[2],
			v[0], v[1], v[2]);

		}
		if (n == 1) {
		    pNode->iLower = 0;
		    if (backup == 0) ++nBucket;
		}
		assert(0);
	    }
	    else
	    {
		assert(bnd[subspace]->fMax[d] > 0);
		if (nl > 0) bnd[subspace]->fCenter[d] -= bnd[subspace]->fMax[d];
		else	bnd[subspace]->fCenter[d] += bnd[subspace]->fMax[d];
		if (n > 1) { PUSH(S, iNode); PUSH(D, (subspace << 2) | d); }
		if (n == 1) {
		    SAVE_BOUNDS(pkd,smx,pNode,bnd);
		    pNode->iLower = 0;
		    if (backup == 0) ++nBucket;
		    }
		}
	    }

	}
DonePart:
    FREE_STACK(S);
    FREE_STACK(D);
    free(keys);
    }

void psdBuildTree(PKD pkd, PSX psx, struct inPSD *in, KDN *pkdn) {
    int iStart;

    assert(pkd->oNodeVBnd);
    assert(pkd->oVelocity);
    assert(pkd->oMass);

    if (pkd->nNodes > 0) {
	/*
	** Close cell caching space and free up nodes.
	*/
	mdlFinishCache(pkd->mdl,CID_CELL);
	}

#ifdef USE_BSC
    MPItrace_event(10000, 1 );
#endif

    pkdClearTimer(pkd,0);
    pkdStartTimer(pkd,0);

    psdInitializeParticles(pkd,&pkd->bnd, &pkd->vbnd);

    BuildPsdTemp(pkd,psx, ROOT,in->nBucket, 102400);

    pkd->nNodesFull = pkd->nNodes;
    iStart = 0;

    pkdStopTimer(pkd,0);
#ifdef USE_BSC
    MPItrace_event(10000, 0 );
#endif
    /*
    ** Finally activate a read only cache for remote access.
    */
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,
	pkd->iTreeNodeSize,pkd->nNodes);
    /*
    ** Copy the root node for the top-tree construction.
    */
    pkdCopyNode(pkd,pkdn,pkdTreeNode(pkd,ROOT));
    mdlFinishCache(pkd->mdl,CID_CELL);
    }

