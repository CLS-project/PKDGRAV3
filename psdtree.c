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

#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif
#include <sys/time.h>


#define X 0
#define V 1

void psdInitializeParticles(PKD pkd, BND *prbnd, BND *pvbnd) {
    KDN *pNode;
    int j;
    pBND rbnd, vbnd;

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

    pkdNodeBnd(pkd,pNode, &rbnd);
    pkdNodeVBnd(pkd,pNode, &vbnd);
    for (j=0;j<3;++j) {
	rbnd.fCenter[j] = prbnd->fCenter[j];
	rbnd.fMax[j]    = prbnd->fMax[j];
	vbnd.fCenter[j] = pvbnd->fCenter[j];
	vbnd.fMax[j]    = pvbnd->fMax[j];
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

typedef struct {
    int64_t key;
    float fMass;
} KEY;

static int key_compar(const void *a0, const void *b0) {
    KEY *a = (KEY *)a0;
    KEY *b = (KEY *)b0;
    return a->key - b->key;
    }

/*
** Compute entropy. E = -\sum (p ln(p))
*/
static inline float E(KEY *keys, int N, float M) {
    int i;
    float s = 0;
    int64_t key = keys[0].key;
    float   p   = keys[0].fMass;
    assert(M > 0);
    for (i=1; i < N; i++) {
	if (keys[i].key != key) {
	    p /= M;
	    s -= p * log(p);
	    p = 0;
	    key = keys[i].key;
	    }

	p += keys[i].fMass;
	}

    if (p != 0)
    {
	p /= M;
	s -= p * log(p);
    }

    return s;
}

void SAVE_BOUNDS(PKD pkd, PSX smx, KDN *pNode, pBND bnd[2])
{
    int j;
    PARTICLE *p = pkdParticle(pkd,pNode->pLower);
    PSMETRIC *b = smx->psm + pNode->pLower;

    for (j=0;j<3;++j) assert(!isinf(bnd[0].fMax[j]));
    for (j=0;j<3;++j) assert(!isinf(bnd[1].fMax[j]));
    for (j=0;j<3;++j) assert(!isinf(bnd[0].fCenter[j]));
    for (j=0;j<3;++j) assert(!isinf(bnd[1].fCenter[j]));

    for (j=0;j<3;++j) assert(fabs(bnd[0].fMax[j]) < 1e10);
    for (j=0;j<3;++j) assert(fabs(bnd[1].fMax[j]) < 1e10);
    for (j=0;j<3;++j) assert(fabs(bnd[0].fCenter[j]) < 1e10);
    for (j=0;j<3;++j) assert(fabs(bnd[1].fCenter[j]) < 1e10);

    p->fDensity = pkdMass(pkd, p);
    for (j=0;j<3;++j)
    {
	b->rscale[j] = 0;
	b->vscale[j] = 0;
	if (bnd[X].fMax[j] != 0) 
	{
	    p->fDensity /= 2*bnd[X].fMax[j];
	    b->rscale[j] = 1. / bnd[X].fMax[j]; 
	}
	if (bnd[V].fMax[j] != 0) 
	{
	    p->fDensity /= 2*bnd[V].fMax[j];
	    b->vscale[j] = 1. / bnd[V].fMax[j];
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
    pBND bnd[2];
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
    int Nb;
    int backup = 0;
    int Ntotal;

    double *v;

    pNode = pkdTreeNode(pkd,iNode);

    Ntotal = pNode->pUpper - pNode->pLower + 1;

    KEY *xkeys = malloc((Ntotal+1) * sizeof(*xkeys)); assert(xkeys != NULL);
    KEY *vkeys = malloc((Ntotal+1) * sizeof(*vkeys)); assert(vkeys != NULL);

    assert(maxNb > 0);
    //assert(maxNb <= (1<<10));

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

    pkdNodeBnd(pkd, pNode, &bnd[0]);
    pkdNodeVBnd(pkd, pNode, &bnd[1]);

    assert( bnd[0].fMax[0] > 0.0
	 || bnd[0].fMax[1] > 0.0
	 || bnd[0].fMax[2] > 0.0
	 || bnd[1].fMax[0] > 0.0
	 || bnd[1].fMax[1] > 0.0
	 || bnd[1].fMax[2] > 0.0);

    assert(STACK_EMPTY(S));
    PUSH(S, iNode); PUSH(D, 0);

    while (!STACK_EMPTY(S)) {

	//fprintf(stderr, "Iter %i\n", iter++);

	iNode = POP(S);

	if (iNode == 0) {
	    //assert(0);
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

	    if (STACK_EMPTY(S)) break; /* Can happen at the end of the tree */
	    iNode = POP(S);
	    assert(iNode != 0);
	    }

	assert(!STACK_EMPTY(D));
	int q = POP(D);
	lastd = q & 0x03;
	last_subspace = (q & 0x04) >> 2;

	pNode = pkdTreeNode(pkd,iNode);
	pkdNodeBnd(pkd, pNode, &bnd[0]);
	pkdNodeVBnd(pkd, pNode, &bnd[1]);

	int Npart = pNode->pUpper - pNode->pLower + 1;

	//fprintf(stderr, "(ROOT %i) %i %i\n", ROOT, iNode,Npart);

#if 0
	if (!( bnd[0].fMax[0] > 0.0 &&
	       bnd[0].fMax[1] > 0.0 &&
	       bnd[0].fMax[2] > 0.0 &&
	       bnd[1].fMax[0] > 0.0 &&
	       bnd[1].fMax[1] > 0.0 &&
	       bnd[1].fMax[2] > 0.0 
	       ))
	    fprintf(stderr, "%g %g %g %g %g %g\n",
		bnd[0].fMax[0],
		bnd[0].fMax[1],
		bnd[0].fMax[2],
		bnd[1].fMax[0],
		bnd[1].fMax[1],
		bnd[1].fMax[2]);
#endif

	assert( bnd[0].fMax[0] > 0.0 
	     || bnd[0].fMax[1] > 0.0 
	     || bnd[0].fMax[2] > 0.0
	     || bnd[1].fMax[0] > 0.0 
	     || bnd[1].fMax[1] > 0.0 
	     || bnd[1].fMax[2] > 0.0 
		);

	assert( !isinf(bnd[0].fMax[0])
	     || !isinf(bnd[0].fMax[1])
	     || !isinf(bnd[0].fMax[2])
	     || !isinf(bnd[1].fMax[0])
	     || !isinf(bnd[1].fMax[1])
	     || !isinf(bnd[1].fMax[2])
		);


    /*************************************************************************/

	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/
	FLOAT dx_inv[2][3];
	FLOAT edge[2][3];
	float e[2];
	float Mtotal = 0;
	int b[2][3];

	Nb = 1+cbrt(Npart);
	if (Nb > maxNb) Nb = maxNb;
	assert(Nb > 1);

	float fb = 0.5 * pow(Ntotal, 1./6);
	float mean_sep;

	int n_node = pNode->pUpper - pNode->pLower + 1;

	if (n_node > 7)
	{
	    s = 0;
	    for (i=0; i < 3; i++)
	    {
		float l_min = bnd[s].fCenter[i] - bnd[s].fMax[i];
		float l_max = bnd[s].fCenter[i] + bnd[s].fMax[i];

		float x_max = -HUGE_VAL;
		float x_min = HUGE_VAL;
		for (j=pNode->pLower; j <= pNode->pUpper; j++) 
		{
		    PARTICLE *p = pkdParticle(pkd,j);
		    if (p->r[i] > x_max) x_max = p->r[i];
		    if (p->r[i] < x_min) x_min = p->r[i];
		}

		mean_sep = (x_max - x_min) / (n_node - 1);

		if ((l_max - x_max) > fb * mean_sep
		&&  (x_min - l_min) > fb * mean_sep)
		{
		    bnd[s].fCenter[i] = (x_max + x_min) / 2;
		    bnd[s].fMax[i]    = (x_max - x_min) / 2;
		}
	    }

	    s = 1;
	    for (i=0; i < 3; i++)
	    {
		float l_min = bnd[s].fCenter[i] - bnd[s].fMax[i];
		float l_max = bnd[s].fCenter[i] + bnd[s].fMax[i];

		float x_max = -HUGE_VAL;
		float x_min = HUGE_VAL;
		for (j=pNode->pLower; j <= pNode->pUpper; j++) 
		{
		    PARTICLE *p = pkdParticle(pkd,j);
		    v = pkdVel(pkd, p);
		    if (v[i] > x_max) x_max = v[i];
		    if (v[i] < x_min) x_min = v[i];
		}

		float mean_sep = (x_max - x_min) / (n_node - 1);

		if ((l_max - x_max) > fb * mean_sep
		&&  (x_min - l_min) > fb * mean_sep)
		{
		    bnd[s].fCenter[i] = (x_max + x_min) / 2;
		    bnd[s].fMax[i] = (x_max - x_min) / 2;
		}
	    }
	}


	for (s=0; s < 2; s++)
	{
	    for (i=0; i < 3; i++)
	    {
		if (bnd[s].fMax[i] > 0)
		    dx_inv[s][i] = (FLOAT)Nb / (2*bnd[s].fMax[i]);
		else
		    dx_inv[s][i] = 0;

		edge[s][i]   = bnd[s].fCenter[i]-bnd[s].fMax[i]; 
	    }
	}

	for (j=pNode->pLower; j <= pNode->pUpper; j++) 
	{
	    PARTICLE *p = pkdParticle(pkd,j);
	    v = pkdVel(pkd, p);
	    float mass = pkdMass(pkd,p);
	    Mtotal += mass;

	    s = 0;
	    for (i=0; i < 3; i++)
	    {
		b[s][i] = (int)((p->r[i]-edge[s][i]) * dx_inv[s][i]); 
		b[s][i] -= (b[s][i] == Nb);
	    }

	    s = 1;
	    for (i=0; i < 3; i++)
	    {
		b[s][i] = (int)((v[i]-edge[s][i]) * dx_inv[s][i]); 
		b[s][i] -= (b[s][i] == Nb);
	    }

	    xkeys[j].key = b[0][0] + Nb * (b[0][1] + Nb*b[0][2]);  xkeys[j].fMass = mass;
	    vkeys[j].key = b[1][0] + Nb * (b[1][1] + Nb*b[1][2]);  vkeys[j].fMass = mass;
	}

	qsort(xkeys + pNode->pLower, Npart, sizeof(*xkeys), key_compar);
	qsort(vkeys + pNode->pLower, Npart, sizeof(*vkeys), key_compar);

	e[X] = E(xkeys + pNode->pLower, Npart, Mtotal);
	e[V] = E(vkeys + pNode->pLower, Npart, Mtotal);


	//fprintf(stderr, "Entropy %e %e\n", e[0], e[1]);

	if (e[V] != 0 
	&& (e[V] < e[X]  ||  (e[V] == e[X] && last_subspace == X)))
	{
	    subspace = V;
	}
	else
	{
	    subspace = X;
	}

	d = max_bnd_range(bnd[subspace].fMax, 0, 3);

	if (d == -1)
	{
	    subspace = !subspace;
	    d = max_bnd_range(bnd[subspace].fMax, 0, 3);
	    assert(d != -1);
	}

#if 0
	// simulate normal tree build
	//subspace = V;
	subspace = !last_subspace;
	d = max_bnd_range(bnd[subspace].fMax, 0, 3);
	assert(d != -1);
#endif
	fSplit = bnd[subspace].fCenter[d];

	//fprintf(stderr, "subspace %i  d %i\n", subspace, d);

	assert(subspace == X || subspace == V);
	assert(0 <= d && d < 3);

    /*************************************************************************/


	//pNode->rbnd.lastd = d;

	//printf("%g %g %g %g\n", fSplit, pNode->rbnd.fCenter[d]-pNode->rbnd.fMax[d], pNode->rbnd.fCenter[d], pNode->rbnd.fCenter[d]+pNode->rbnd.fMax[d]);
	//printf("      %g %g\n", fSplit, pNode->rbnd.fMax[d]);
	//assert((pNode->rbnd.fCenter[d]-pNode->rbnd.fMax[d]) <= fSplit && fSplit <= (pNode->rbnd.fCenter[d]+pNode->rbnd.fMax[d]));
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
	    PARTITION(pi<pj,pi<=pj,
		   pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		   pkdSwapParticle(pkd, pi,pj),
		   pi->r[d] < fSplit, pj->r[d] >= fSplit);
	}
	else
	{
	    PARTITION(pi<pj,pi<=pj,
		   pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		   pkdSwapParticle(pkd, pi,pj),
		   pkdVel(pkd, pi)[d] < fSplit, pkdVel(pkd,pj)[d] >= fSplit);
	}

	nl = i - pNode->pLower;
	nr = pNode->pUpper - i + 1;

	//printf("%i %.15f  [%g %g] (%i %i) %i %i %s\n", d, fSplit, rbnd.fCenter[d]-rbnd.fMax[d], rbnd.fCenter[d]+rbnd.fMax[d], nl, nr, iNode, pNode->iParent, backup == 0 ? "" : "*");

    /*************************************************************************/

	if (nl > 0 && nr > 0) {
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

	    pBND rbnd[2], lbnd[2];

	    pkdNodeBnd(pkd, pRight, &rbnd[0]); pkdNodeVBnd(pkd, pRight, &rbnd[1]);
	    pkdNodeBnd(pkd, pLeft,  &lbnd[0]); pkdNodeVBnd(pkd, pLeft,  &lbnd[1]);

	    //pRight->bnd.lastd = d;
	    //pLeft->rbnd.lastd = d;

	    //fprintf(stderr, "%i %i %i\n", iNode, iLeft, iRight);

	    /*
	    ** Now deal with the bounds.
	    **
	    ** ADD SHRINK WRAPPING -- jpc 27.3.2010
	    **
	    */
	    for (s=0; s < 2; s++)
	    {
		for (j=0;j<3;++j) {
		    if (s == subspace && j == d) {
			rbnd[s].fMax[j] = lbnd[s].fMax[j] = 0.5*bnd[s].fMax[j];
			lbnd[s].fCenter[j] = bnd[s].fCenter[j] - lbnd[s].fMax[j];
			rbnd[s].fCenter[j] = bnd[s].fCenter[j] + rbnd[s].fMax[j];
			}
		    else {
			lbnd[s].fCenter[j] = bnd[s].fCenter[j]; lbnd[s].fMax[j] = bnd[s].fMax[j];
			rbnd[s].fCenter[j] = bnd[s].fCenter[j]; rbnd[s].fMax[j] = bnd[s].fMax[j];
			}
		    //assert(lbnd[0].fMax[j] > 0.0);
		    //assert(rbnd[0].fMax[j] > 0.0);
		    }
	    }

	    //ls = max_side(lbnd[0].fMax);     // MAXSIDE(pLeft.bnd.fMax,ls);
	    //rs = max_side(rbnd[0].fMax);     // MAXSIDE(pRight.rbnd.fMax,rs);
	    /*
	    ** Now figure out which subfile to process next.
	    */
	    lc = ((nl > M)); // ||((nl > 1)&&(ls>PKD_MAX_CELL_SIZE))); /* this condition means the left child is not a bucket */
	    rc = ((nr > M)); // ||((nr > 1)&&(rs>PKD_MAX_CELL_SIZE)));
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

		if (nr > 1) {
		    backup = pkd->nNodes; //fprintf(stderr, "backup [%i]\n", backup); 
		    PUSH(S, 0);
		    PUSH(S, iRight); PUSH(D, (subspace << 2) | d);
		    ++nBucket;
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

		if (nl > 1) {
		    backup = pkd->nNodes; //fprintf(stderr, "backup [%i]\n", backup); 
		    PUSH(S, 0);
		    PUSH(S, iLeft); PUSH(D, (subspace << 2) | d);
		    ++nBucket;
		    }
		else {
		    pLeft->iLower = 0;
		    SAVE_BOUNDS(pkd,smx,pLeft,lbnd);
		    }
		}
	    else {
		/*
		** Both are buckets (we need to pop from the stack to get the next subfile.)
		*/

		if (!(nr==1 && nl==1) && backup == 0) { 
		    backup = pkd->nNodes; //fprintf(stderr, "backup [%i]\n", backup); 
		    PUSH(S, 0);
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

	    bnd[subspace].fMax[d] *= 0.5;
	    if (bnd[subspace].fMax[d] <= 0)
		fprintf(stderr, "subspace %i  d %i\n", subspace, d);
	    assert(bnd[subspace].fMax[d] > 0);
	    if (nl > 0) bnd[subspace].fCenter[d] -= bnd[subspace].fMax[d];
	    else	bnd[subspace].fCenter[d] += bnd[subspace].fMax[d];

	    //ls = max_side(pNode->rbnd.fMax); // MAXSIDE(pNode->rbnd.fMax,ls);
	    //lc = ((n > M)||((n > 1)&&(ls>PKD_MAX_CELL_SIZE))); /* this condition means the node is not a bucket */

	    if (n > 1) { PUSH(S, iNode); PUSH(D, (subspace << 2) | d); }

	    if (n == 1) {
		SAVE_BOUNDS(pkd,smx,pNode,bnd);
		pNode->iLower = 0;
		if (backup == 0) ++nBucket;
		}
	    }
	}
DonePart:
    FREE_STACK(S);
    FREE_STACK(D);
    free(xkeys);
    free(vkeys);
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
    }

