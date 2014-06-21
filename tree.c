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

#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif





static void InitializeParticles(PKD pkd,int bExcludeVeryActive,BND *pbnd) {
    PLITE *pLite = pkd->pLite;
    PLITE t;
    PARTICLE *p;
    KDN *pNode;
    BND *bnd;
    int i,j;

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) pLite[i].r[j] = p->r[j];
	pLite[i].i = i;
	pLite[i].uRung = (bExcludeVeryActive==2?p->bMarked:p->uRung);
	}
    /*
    **It is only forseen that there are 4 reserved nodes at present 0-NULL, 1-ROOT, 2-UNUSED, 3-VAROOT.
    */
    pkd->nNodes = NRESERVED_NODES;
    /*
    ** If we want to split out very active particles from the tree
    ** we do it here, collecting them in Node index 0 as a bucket
    ** with possibly many particles.
    */
    pkd->nVeryActive = 0;
    /* We want only marked particles */
    if (bExcludeVeryActive == 2) {
	i = 0;
	j = pkd->nLocal - 1;
	while (i <= j) {
	    if ( pLite[i].uRung ) ++i;
	    else break;
	    }
	while (i <= j) {
	    if ( !pLite[j].uRung ) --j;
	    else break;
	    }
	if (i < j) {
	    t = pLite[i];
	    pLite[i] = pLite[j];
	    pLite[j] = t;
	    while (1) {
		while (pLite[++i].uRung);
		while (!pLite[--j].uRung);
		if (i < j) {
		    t = pLite[i];
		    pLite[i] = pLite[j];
		    pLite[j] = t;
		    }
		else break;
		}
	    }
	pkd->nVeryActive = pkd->nLocal - i;
	}
    else if (bExcludeVeryActive) {
	/*
	** Now start the partitioning of the particles.
	*/
	i = 0;
	j = pkd->nLocal - 1;
	while (i <= j) {
	    if ( pLite[i].uRung <= pkdRungVeryActive(pkd) ) ++i;
	    else break;
	    }
	while (i <= j) {
	    if ( pLite[j].uRung > pkdRungVeryActive(pkd) ) --j;
	    else break;
	    }
	if (i < j) {
	    t = pLite[i];
	    pLite[i] = pLite[j];
	    pLite[j] = t;
	    while (1) {
		while ((pLite[++i].uRung <= pkdRungVeryActive(pkd)));
		while (pLite[--j].uRung > pkdRungVeryActive(pkd));
		if (i < j) {
		    t = pLite[i];
		    pLite[i] = pLite[j];
		    pLite[j] = t;
		    }
		else break;
		}
	    }
	pkd->nVeryActive = pkd->nLocal - i;

	pNode = pkdTreeNode(pkd,VAROOT);
        bnd = pkdNodeBnd(pkd, pNode);
	if (pkd->nVeryActive > 0)
	    /*
	    ** Set up the very active root node.
	    */
	    pNode->iLower = 0;
	pNode->iParent = 0;
	pNode->pLower = pkd->nLocal - pkd->nVeryActive;
	pNode->pUpper = pkd->nLocal - 1;
	for (j=0;j<3;++j) {
	    bnd->fCenter[j] = pbnd->fCenter[j];
	    bnd->fMax[j] = pbnd->fMax[j];
	    }
	}
    /*
    ** Set up the root node.
    */
    pNode = pkdTreeNode(pkd,ROOT);
    pNode->iLower = 0;
    pNode->iParent = 0;
    pNode->pLower = 0;
    pNode->pUpper = pkd->nLocal - pkd->nVeryActive - 1;
    bnd = pkdNodeBnd(pkd, pNode);
    for (j=0;j<3;++j) {
	bnd->fCenter[j] = pbnd->fCenter[j];
	bnd->fMax[j] = pbnd->fMax[j];
	}
    }

#define MIN_SRATIO    0.05

/*
** M is the bucket size.
** This function assumes that the root node is correctly set up (particularly the bounds).
*/
#define TEMP_S_INCREASE 100
void BuildTemp(PKD pkd,int iNode,int M) {
    PLITE *p = pkd->pLite;
    KDN *pNode = pkdTreeNode(pkd,iNode);
    BND *bnd,*lbnd,*rbnd;
    KDN *pLeft, *pRight;
    PLITE t;
    FLOAT fSplit;
    int *S;		/* this is the stack */
    int s,ns;
    int iLeft,iRight;
    int d,i,j;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;

    bnd = pkdNodeBnd(pkd,pNode);

    /*
    ** Allocate stack!
    */
    ns = TEMP_S_INCREASE;
    s = 0;
    S = CAST(int*,malloc(ns*sizeof(int)));
    assert(S != NULL);
    /*
    ** Make sure we don't have buckets which are larger than the
    ** pointer arrays for actives and inactives!
    */
    if (M > pkd->nMaxBucketActive) {
	pkd->nMaxBucketActive = M;
	pkd->piActive = CAST(PARTICLE **,realloc(pkd->piActive,pkd->nMaxBucketActive*sizeof(PARTICLE *)));
	mdlassert(pkd->mdl,pkd->piActive != NULL);
	pkd->piInactive = CAST(PARTICLE **,realloc(pkd->piInactive,pkd->nMaxBucketActive*sizeof(PARTICLE *)));
	mdlassert(pkd->mdl,pkd->piInactive != NULL);
	}

    if (pNode->pUpper - pNode->pLower + 1 <= M)
	goto DonePart;

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
	fSplit = bnd->fCenter[d];
	/*
	** Now start the partitioning of the particles about
	** fSplit on dimension given by d.
	*/
	i = pNode->pLower;
	j = pNode->pUpper;
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

	nl = i - pNode->pLower;
	nr = pNode->pUpper - i + 1;
	if (nl > 0 && nr > 0) {
	    /*
	    ** Allocate 2 new tree nodes making sure that we have
	    ** allocated enough storage.
	    */
	    pkdTreeAllocNodePair(pkd,&iLeft,&iRight);
	    pLeft = pkdTreeNode(pkd,iLeft);
	    pLeft->iParent = iNode;
	    pLeft->pLower = pNode->pLower;
	    pLeft->pUpper = i-1;
	    pRight = pkdTreeNode(pkd,iRight);
	    assert(iRight & 1);
	    pRight->iParent = iNode;
	    pRight->pLower = i;
	    pRight->pUpper = pNode->pUpper;
	    pNode->iLower = iLeft;

            lbnd = pkdNodeBnd(pkd, pLeft);
            rbnd = pkdNodeBnd(pkd, pRight);

	    /*
	    ** Now deal with the bounds.
	    */
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
		}
	    /*
	    ** Now figure out which subfile to process next.
	    */
	    lc = (nl > M); /* this condition means the left child is not a bucket */
	    rc = (nr > M);
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
	    /*
	    ** No nodes allocated, Change the bounds if needed!
	    */
	    if (d >= 0 && d < 3) bnd->fMax[d] *= 0.5;
	    if (nl > 0) {
		if (d >= 0 && d < 3) bnd->fCenter[d] -= bnd->fMax[d];
		lc = (nl > M); /* this condition means the node is not a bucket */
		if (!lc) {
		    pNode->iLower = 0;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    else {
		if (d >= 0 && d < 3) bnd->fCenter[d] += bnd->fMax[d];
		rc = (nr > M);
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
DonePart:
    free(S);
    }

/*
** If this is called with iStart being the index of the first very active particle
** then it reshuffles only the very actives. This is again a bit ugly, but will
** do for now.
*/
void ShuffleParticles(PKD pkd,int iStart) {
    PARTICLE *p, *pNew, *pNewer;
    int i,iNew,iNewer,iTemp;

    /*
    ** Now we move the particles in one go using the temporary
    ** particles which have been shuffled.
    */
    for(iTemp=iStart; iTemp<pkd->nLocal; ++iTemp) {
	if (pkd->pLite[iTemp].i==iTemp) continue;
	i = iTemp;
	p = pkdParticle(pkd,i);
	pkdSaveParticle(pkd,p);
	iNew = pkd->pLite[i].i;
	pNew = pkdParticle(pkd,iNew);
	while(iNew!=iTemp) {
	    pkdCopyParticle(pkd,p,pNew);
	    pkd->pLite[i].i = i;
	    i = iNew;
	    p = pNew;
	    iNew = pkd->pLite[i].i;
	    pNew = pkdParticle(pkd,iNew);
	    }
	pkd->pLite[i].i = i;
	pkdLoadParticle(pkd,p);
	}
    }

static double zeroV[3] = {0.0,0.0,0.0};
static float  zeroF[3] = {0.0,0.0,0.0};

void Create(PKD pkd,int iNode) {
    PARTICLE *p;
    KDN *pkdn,*pkdl,*pkdu;
    FMOMR mom;
    SPHBNDS *bn;
    BND *bnd;
    FLOAT m,fMass,fSoft,x,y,z,vx,vy,vz,ax,ay,az,ft,d2,d2Max,dih2,bmin,b;
    float *a;
    double *v;
    int pj,d,nDepth,ism;
    const int nMaxStackIncrease = 1;

    nDepth = 1;
    while (1) {
	while ((pkdn=pkdTreeNode(pkd,iNode))->iLower) {
	    ++nDepth;
	    iNode = pkdn->iLower;
	    }
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
		}
	    pkd->nMaxStack += nMaxStackIncrease;
	    }
	/*
	** Now calculate all bucket quantities!
	** This includes M,CoM,Moments and special
	** bounds and iMaxRung.
	*/
	pkdn = pkdTreeNode(pkd,iNode);
        bnd = pkdNodeBnd(pkd, pkdn);
	pkdn->nActive = 0;
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
	    ft = p->r[d];
	    bnd->fCenter[d] = ft;
	    bnd->fMax[d] = ft;
	    }
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    for (d=0;d<3;++d) {
		ft = p->r[d];
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
	x = m*p->r[0];
	y = m*p->r[1];
	z = m*p->r[2];
	vx = m*v[0];
	vy = m*v[1];
	vz = m*v[2];
	ax = m*a[0];
	ay = m*a[1];
	az = m*a[2];
	pkdn->uMinRung = pkdn->uMaxRung = p->uRung;
	pkdn->bDstActive = p->bDstActive;
	if (pkdIsActive(pkd,p)) ++pkdn->nActive;
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    a = pkd->oAcceleration ? pkdAccel(pkd,p) : zeroF;
	    m = pkdMass(pkd,p);
	    fSoft = pkdSoft(pkd,p);
	    v = pkd->oVelocity ? pkdVel(pkd,p) : zeroV;
	    fMass += m;
	    if (fSoft>dih2) dih2=fSoft;
	    x += m*p->r[0];
	    y += m*p->r[1];
	    z += m*p->r[2];
	    vx += m*v[0];
	    vy += m*v[1];
	    vz += m*v[2];
	    ax += m*a[0];
	    ay += m*a[1];
	    az += m*a[2];
	    if ( p->uRung > pkdn->uMaxRung ) pkdn->uMaxRung = p->uRung;
	    if ( p->uRung < pkdn->uMinRung ) pkdn->uMinRung = p->uRung;
	    if ( p->bDstActive ) pkdn->bDstActive = 1;
	    if (pkdIsActive(pkd,p)) ++pkdn->nActive;
	    }
	m = 1/fMass;
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
	    double *pVel = pkdNodeVel(pkd,pkdn);
	    pVel[0] = m*vx;
	    pVel[1] = m*vy;
	    pVel[2] = m*vz;
	    }
	if (pkd->oNodeAcceleration) {
	    double *pAcc = pkdNodeAccel(pkd,pkdn);
	    pAcc[0] = m*ax;
	    pAcc[1] = m*ay;
	    pAcc[2] = m*az;
	    }
	pkdn->fSoft2 = dih2*dih2;
	d2Max = 0.0;
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    x = p->r[0] - pkdn->r[0];
	    y = p->r[1] - pkdn->r[1];
	    z = p->r[2] - pkdn->r[2];
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
	/*
	** Now calculate the reduced multipole moment.
	** Note that we use the cell's openening radius as the scaling factor!
	*/
	if (pkd->oNodeMom) {
	    momClearFmomr(pkdNodeMom(pkd,pkdn));
	    for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		p = pkdParticle(pkd,pj);
		x = p->r[0] - pkdn->r[0];
		y = p->r[1] - pkdn->r[1];
		z = p->r[2] - pkdn->r[2];
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
		if (pkdIsGas(pkd,p)) {
		    /*
		    ** This first ball bound over all gas particles is only used for remote searching.
		    */
		    for (d=0;d<3;++d) bn->B.min[d] = fmin(bn->B.min[d],p->r[d] - (1+pkd->param.ddHonHLimit)*p->fBall);
		    for (d=0;d<3;++d) bn->B.max[d] = fmax(bn->B.max[d],p->r[d] + (1+pkd->param.ddHonHLimit)*p->fBall);
		    if (pkdIsActive(pkd,p)) {
			for (d=0;d<3;++d) bn->A.min[d] = fmin(bn->A.min[d],p->r[d]);
			for (d=0;d<3;++d) bn->A.max[d] = fmax(bn->A.max[d],p->r[d]);
		    }
		    else {
			for (d=0;d<3;++d) bn->BI.min[d] = fmin(bn->BI.min[d],p->r[d] - (1+pkd->param.ddHonHLimit)*p->fBall);
			for (d=0;d<3;++d) bn->BI.max[d] = fmax(bn->BI.max[d],p->r[d] + (1+pkd->param.ddHonHLimit)*p->fBall);
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
	while (iNode & 1) {
	    iNode = pkdTreeNode(pkd,iNode)->iParent;
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
		p = pkdParticle(pkd,pj);
		x = p->r[0] - pkdn->r[0];
		y = p->r[1] - pkdn->r[1];
		z = p->r[2] - pkdn->r[2];
		d2Max = x*x + y*y + z*z;
		for (++pj;pj<=pkdn->pUpper;++pj) {
		    p = pkdParticle(pkd,pj);
		    x = p->r[0] - pkdn->r[0];
		    y = p->r[1] - pkdn->r[1];
		    z = p->r[2] - pkdn->r[2];
		    d2 = x*x + y*y + z*z;
		    d2Max = (d2 > d2Max)?d2:d2Max;
		    }
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
    if (0xffffffffu - p1->nActive < p2->nActive) pkdn->nActive = 0xffffffffu; 
    else pkdn->nActive = p1->nActive + p2->nActive;

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
	x = p1->r[0] - pkdn->r[0];
	y = p1->r[1] - pkdn->r[1];
	z = p1->r[2] - pkdn->r[2];
	momShiftFmomr(pkdNodeMom(pkd,pkdn),p1->bMax,x,y,z);

	momRescaleFmomr(pkdNodeMom(pkd,pkdn),pkdn->bMax,p1->bMax);

	mom = *pkdNodeMom(pkd,p2);
	x = p2->r[0] - pkdn->r[0];
	y = p2->r[1] - pkdn->r[1];
	z = p2->r[2] - pkdn->r[2];
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


void pkdVATreeBuild(PKD pkd,int nBucket) {
    PARTICLE *p;
    int i,j,iStart;

    iStart = pkd->nLocal - pkd->nVeryActive;
    /*
    ** First initialize the very active temporary particles.
    */
    for (i=iStart;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) pkd->pLite[i].r[j] = p->r[j];
	pkd->pLite[i].i = i;
	pkd->pLite[i].uRung = p->uRung;
	}
    /*
    ** Then clear the VA tree by setting the node index back to one node past the end
    ** of the non VA tree.
    */
    pkd->nNodes = pkd->nNonVANodes;
    BuildTemp(pkd,VAROOT,nBucket);

    ShuffleParticles(pkd,iStart);

    Create(pkd,VAROOT);
    }


void pkdTreeBuild(PKD pkd,int nBucket,KDN *pkdn,int bExcludeVeryActive) {
    int iStart;

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

    InitializeParticles(pkd,bExcludeVeryActive,&pkd->bnd);

    BuildTemp(pkd,ROOT,nBucket);
    if (bExcludeVeryActive) {
	pkd->nNonVANodes = pkd->nNodes;
	}
    else {
	pkd->nNodesFull = pkd->nNodes;
	}
    iStart = 0;
    ShuffleParticles(pkd,iStart);
    Create(pkd,ROOT);

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

void pkdTreeBuildByGroup(PKD pkd, int nBucket) {
    PLITE *pLite = pkd->pLite;
    PARTICLE *p;
    KDN *pNode;
    BND *bnd;
    double dMin[3], dMax[3];
    int i,j,k,n,gid,gid2,iRoot;

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
    pkd->nVeryActive = 0;

    if (pkd->hopSavedRoots == 0) {
	/* Sort particle by group, but with group 0 at the end */
	int *iGrpOffset = malloc(sizeof(int)*(pkd->nGroups+1));
	for (i=0; i<=pkd->nGroups;++i) iGrpOffset[i] = 0;
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    gid = *pkdGroup(pkd,p);
	    ++iGrpOffset[gid+1];
	    }
	iGrpOffset[0] = iGrpOffset[1];
	iGrpOffset[1] = 0;
	for(i=2; i<=pkd->nGroups;++i) {
	    iGrpOffset[i] += iGrpOffset[i-1];
	    }
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    gid = *pkdGroup(pkd,p);
	    if (gid==0) k = iGrpOffset[pkd->nGroups]++;
	    else k = iGrpOffset[gid]++;
	    for (j=0;j<3;++j) pLite[k].r[j] = p->r[j];
	    pLite[k].i = i;
	    pLite[k].uRung = p->uRung;
	    pLite[k].uGroup = gid;
	    }
	i = 0;
	for(gid=1; gid<pkd->nGroups;++gid) {
	    pkdTreeAllocRootNode(pkd,&iRoot);
	    pkd->hopGroups[gid].iTreeRoot = iRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    pNode->pLower = i;
	    i = iGrpOffset[gid];
	    pNode->pUpper = i - 1;
	    }
	free(iGrpOffset);
	for (i=0;i<pkd->nLocal && pLite[i].uGroup;) {
	    gid = pLite[i].uGroup;
	    iRoot = pkd->hopGroups[gid].iTreeRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    bnd = pkdNodeBnd(pkd, pNode);
	    
	    pNode->iLower = 0;
	    pNode->iParent = 0;
	    assert(pNode->pLower == i);

	    for (j=0;j<3;++j) dMin[j] = dMax[j] = pLite[i].r[j];
	    for(++i; i<pkd->nLocal && pLite[i].uGroup==gid; ++i) {
		pkdMinMax(pLite[i].r,dMin,dMax);
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
	for (i=0;i<pkd->nLocal;++i) pLite[i].i = i;
	for(gid2=1; gid2<pkd->nGroups;++gid2) {
	    if (!pkd->hopGroups[gid2].bNeedGrav) continue;
	    iRoot = pkd->hopGroups[gid2].iTreeRoot;
	    pNode = pkdTreeNode(pkd,iRoot);
	    pNode->iLower = 0;
	    pNode->iParent = 0;
	    bnd = pkdNodeBnd(pkd, pNode);
	    n = pNode->pUpper;
	    for(i=k=pNode->pLower; i<=pNode->pUpper; ++i) {
		p = pkdParticle(pkd,i);
		gid = *pkdGroup(pkd,p);
		if (gid) {
		    assert(gid==gid2);
		    for (j=0;j<3;++j) pLite[k].r[j] = p->r[j];
		    pLite[k].i = i;
		    pLite[k].uRung = p->uRung;
		    pLite[k].uGroup = gid;
		    if (k==pNode->pLower) {
			for (j=0;j<3;++j) dMin[j] = dMax[j] = pLite[k].r[j];
			}
		    else {
			pkdMinMax(pLite[k].r,dMin,dMax);
			}
		    ++k;
		    }
		else {
		    for (j=0;j<3;++j) pLite[n].r[j] = p->r[j];
		    pLite[n].i = i;
		    pLite[n].uRung = p->uRung;
		    pLite[n].uGroup = gid;
		    --n;
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
	    BuildTemp(pkd,pkd->hopGroups[gid].iTreeRoot,nBucket);
    ShuffleParticles(pkd,0);
    for(gid=1; gid<pkd->nGroups; ++gid)
	if (pkd->hopGroups[gid].bNeedGrav)
	    Create(pkd,pkd->hopGroups[gid].iTreeRoot);
    /*
    ** Finally activate a read only cache for remote access.
    */
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,
	pkd->iTreeNodeSize,pkd->nNodes);


    }

void pkdDistribCells(PKD pkd,int nCell,KDN *pkdn) {
    KDN *pSrc, *pDst;
    int i;

    pkdAllocateTopTree(pkd,nCell);
    for (i=1;i<nCell;++i) {
	pSrc = pkdNode(pkd,pkdn,i);
	if (pSrc->pUpper) {
	    pDst = pkdTopNode(pkd,i);
	    pkdCopyNode(pkd,pDst,pSrc);
	    if (pDst->pLower == pkd->idSelf) pkd->iTopRoot = i;
	    }
	}
    }


/*
** Hopefully we can bypass this step once we figure out how to do the
** Multipole Ewald with reduced multipoles.
*/
void pkdCalcRoot(PKD pkd,MOMC *pmom) {
    PARTICLE *p;
    FLOAT xr = pkdTopNode(pkd,ROOT)->r[0];
    FLOAT yr = pkdTopNode(pkd,ROOT)->r[1];
    FLOAT zr = pkdTopNode(pkd,ROOT)->r[2];
    FLOAT x,y,z;
    FLOAT fMass;
    MOMC mc;
    int i = 0;

    p = pkdParticle(pkd,i);
    x = p->r[0] - xr;
    y = p->r[1] - yr;
    z = p->r[2] - zr;
    fMass = pkdMass(pkd,p);
    momMakeMomc(pmom,fMass,x,y,z);
    for (++i;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	fMass = pkdMass(pkd,p);
	x = p->r[0] - xr;
	y = p->r[1] - yr;
	z = p->r[2] - zr;
	momMakeMomc(&mc,fMass,x,y,z);
	momAddMomc(pmom,&mc);
	}
    }


void pkdDistribRoot(PKD pkd,MOMC *pmom) {
    pkd->ew.mom = *pmom;
    }


void pkdTreeNumSrcActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi) {
    KDN *kdn;
    PARTICLE *p;
    int iNode,pj;

    kdn = pkdTreeNode(pkd,iNode = ROOT);
    while (1) {
	while (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iNode = kdn->iLower);
	    }
	/*
	** We have to test each particle of the bucket for activity.
	*/
	kdn->nActive = 0;
	for (pj=kdn->pLower;pj<=kdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    if (pkdIsSrcActive(p,uRungLo,uRungHi)) ++kdn->nActive;
	}
	while (iNode & 1) {
	    kdn = pkdTreeNode(pkd,iNode = kdn->iParent);
	    if (!iNode) return;	/* exit point!!! */
	    kdn->nActive = pkdTreeNode(pkd,kdn->iLower)->nActive
		+ pkdTreeNode(pkd,kdn->iLower + 1)->nActive;
	}
	kdn = pkdTreeNode(pkd,++iNode);
    }
}


void pkdBoundWalk(PKD pkd,BND *pbnd,uint8_t uRungLo,uint8_t uRungHi,uint32_t *pnActive,uint32_t *pnContained) {
    KDN *kdn;
    PARTICLE *p;
    double d;
    int iNode,pj;    
    BND *bnd;

    *pnActive = 0;
    *pnContained = 0;
    kdn = pkdTreeNode(pkd,iNode = ROOT);
    bnd = pkdNodeBnd(pkd, kdn);
    while (1) {
	d = fabs(pbnd->fCenter[0] - bnd->fCenter[0]) - pbnd->fMax[0];
	if (d - bnd->fMax[0] > 0) goto NoIntersect;
	else if (d + bnd->fMax[0] <= 0) {
	    d = fabs(pbnd->fCenter[1] - bnd->fCenter[1]) - pbnd->fMax[1];
	    if (d - bnd->fMax[1] > 0) goto NoIntersect;
	    else if (d + bnd->fMax[1] <= 0) {
		d = fabs(pbnd->fCenter[2] - bnd->fCenter[2]) - pbnd->fMax[2];
		if (d - bnd->fMax[2] > 0) goto NoIntersect;
		else if (d + bnd->fMax[2] <= 0) goto Contained;
		}
	    else {
		d = fabs(pbnd->fCenter[2] - bnd->fCenter[2]) - pbnd->fMax[2];
		if (d - bnd->fMax[2] > 0) goto NoIntersect;
		}
	    }
	else {
	    d = fabs(pbnd->fCenter[1] - bnd->fCenter[1]) - pbnd->fMax[1];
	    if (d - bnd->fMax[1] > 0) goto NoIntersect;
	    d = fabs(pbnd->fCenter[2] - bnd->fCenter[2]) - pbnd->fMax[2];
	    if (d - bnd->fMax[2] > 0) goto NoIntersect;
	    }	
	/*
	** We have an intersection to test!
	*/
	if (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iNode = kdn->iLower);
	    continue;
	    }
	else {
	    /*
	    ** We have to test each active particle of the bucket for containment.
	    */
	    for (pj=kdn->pLower;pj<=kdn->pUpper;++pj) {
		p = pkdParticle(pkd,pj);
		if (fabs(pbnd->fCenter[0] - p->r[0]) - pbnd->fMax[0] > 0) continue;
		if (fabs(pbnd->fCenter[1] - p->r[1]) - pbnd->fMax[1] > 0) continue;
		if (fabs(pbnd->fCenter[2] - p->r[2]) - pbnd->fMax[2] > 0) continue;
		/*
		** This particle is contained.
		*/
		*pnContained += 1;
		if (pkdIsSrcActive(p,uRungLo,uRungHi)) *pnActive += 1;
		}
	    }

    Contained:
	/*
	** Cell is contained within the bounds.
	*/
	*pnContained += (kdn->pUpper - kdn->pLower + 1);
	*pnActive += kdn->nActive;  /* this must be set with SrcActive for all cells first */
    NoIntersect:
	while (iNode & 1) {
	    kdn = pkdTreeNode(pkd,iNode = kdn->iParent);
	    if (!iNode) return;    /* exit point */
	}
	kdn = pkdTreeNode(pkd,++iNode);
    }
}
