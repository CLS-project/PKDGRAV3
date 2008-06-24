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
#include <sys/time.h>



void InitializeParticles(PKD pkd,int bExcludeVeryActive,BND *pbnd) {
    PLITE *pLite = pkd->pLite;
    PLITE t;
    PARTICLE *p;
    int i,j;

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) pLite[i].r[j] = p->r[j];
	pLite[i].i = i;
	pLite[i].uRung = p->uRung;
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
    if (bExcludeVeryActive) {
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

	if (pkd->nVeryActive > 0)
	    /*   printf("%d:nVeryActive = %d\n",mdlSelf(pkd->mdl),pkd->nVeryActive);*/
	    /*
	    ** Set up the very active root node.
	    */
	    pkd->kdNodes[VAROOT].iLower = 0;
	pkd->kdNodes[VAROOT].iParent = 0;
	pkd->kdNodes[VAROOT].pLower = pkd->nLocal - pkd->nVeryActive;
	pkd->kdNodes[VAROOT].pUpper = pkd->nLocal - 1;
	for (j=0;j<3;++j) {
	    pkd->kdNodes[VAROOT].bnd.fCenter[j] = pbnd->fCenter[j];
	    pkd->kdNodes[VAROOT].bnd.fMax[j] = pbnd->fMax[j];
	    }
	}
    /*
    ** Set up the root node.
    */
    pkd->kdNodes[ROOT].iLower = 0;
    pkd->kdNodes[ROOT].iParent = 0;
    pkd->kdNodes[ROOT].pLower = 0;
    pkd->kdNodes[ROOT].pUpper = pkd->nLocal - pkd->nVeryActive - 1;
    for (j=0;j<3;++j) {
	pkd->kdNodes[ROOT].bnd.fCenter[j] = pbnd->fCenter[j];
	pkd->kdNodes[ROOT].bnd.fMax[j] = pbnd->fMax[j];
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
    PLITE t;
    FLOAT fSplit;
    FLOAT ls,rs;
    int *S;		/* this is the stack */
    int s,ns;
    int iLeft,iRight;
    int d,i,j;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;

    /*
    ** Allocate stack!
    */
    ns = TEMP_S_INCREASE;
    s = 0;
    S = malloc(ns*sizeof(int));
    assert(S != NULL);
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

    if (pkd->kdNodes[iNode].pUpper - pkd->kdNodes[iNode].pLower + 1 <= M)
	goto DonePart;

    assert( pkd->kdNodes[iNode].bnd.fMax[0] > 0.0 ||
	    pkd->kdNodes[iNode].bnd.fMax[1] > 0.0 ||
	    pkd->kdNodes[iNode].bnd.fMax[2] > 0.0 );
    while (1) {
	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/
	if (pkd->kdNodes[iNode].bnd.fMax[0] < pkd->kdNodes[iNode].bnd.fMax[1]) {
	    if (pkd->kdNodes[iNode].bnd.fMax[1] < pkd->kdNodes[iNode].bnd.fMax[2]) d = 2;
	    else d = 1;
	    }
	else if (pkd->kdNodes[iNode].bnd.fMax[0] < pkd->kdNodes[iNode].bnd.fMax[2]) d = 2;
	else d = 0;
	fSplit = pkd->kdNodes[iNode].bnd.fCenter[d];
	/*
	** Now start the partitioning of the particles about
	** fSplit on dimension given by d.
	*/
	i = pkd->kdNodes[iNode].pLower;
	j = pkd->kdNodes[iNode].pUpper;
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

	nl = i - pkd->kdNodes[iNode].pLower;
	nr = pkd->kdNodes[iNode].pUpper - i + 1;
	if (nl > 0 && nr > 0) {
	    /*
	    ** Allocate 2 new tree nodes making sure that we have
	    ** allocated enough storage.
	    */
	    assert(pkd->nNodes+2 <= pkd->nMaxNodes);
	    iLeft = pkd->nNodes++;
	    pkd->kdNodes[iLeft].iParent = iNode;
	    pkd->kdNodes[iLeft].pLower = pkd->kdNodes[iNode].pLower;
	    pkd->kdNodes[iLeft].pUpper = i-1;
	    iRight = pkd->nNodes++;
	    assert(iRight & 1);
	    pkd->kdNodes[iRight].iParent = iNode;
	    pkd->kdNodes[iRight].pLower = i;
	    pkd->kdNodes[iRight].pUpper = pkd->kdNodes[iNode].pUpper;
	    pkd->kdNodes[iNode].iLower = iLeft;
	    /*
	    ** Now deal with the bounds.
	    */
	    for (j=0;j<3;++j) {
		if (j == d) {
		    pkd->kdNodes[iRight].bnd.fMax[j] = pkd->kdNodes[iLeft].bnd.fMax[j] = 0.5*pkd->kdNodes[iNode].bnd.fMax[j];
		    pkd->kdNodes[iLeft].bnd.fCenter[j] = pkd->kdNodes[iNode].bnd.fCenter[j] - pkd->kdNodes[iLeft].bnd.fMax[j];
		    pkd->kdNodes[iRight].bnd.fCenter[j] = pkd->kdNodes[iNode].bnd.fCenter[j] + pkd->kdNodes[iRight].bnd.fMax[j];
		    }
		else {
		    pkd->kdNodes[iLeft].bnd.fCenter[j] = pkd->kdNodes[iNode].bnd.fCenter[j];
		    pkd->kdNodes[iLeft].bnd.fMax[j] = pkd->kdNodes[iNode].bnd.fMax[j];
		    pkd->kdNodes[iRight].bnd.fCenter[j] = pkd->kdNodes[iNode].bnd.fCenter[j];
		    pkd->kdNodes[iRight].bnd.fMax[j] = pkd->kdNodes[iNode].bnd.fMax[j];
		    }
		}
	    MAXSIDE(pkd->kdNodes[iLeft].bnd.fMax,ls);
	    MAXSIDE(pkd->kdNodes[iRight].bnd.fMax,rs);
	    /*
	    ** Now figure out which subfile to process next.
	    */
	    lc = ((nl > M)||((nl > 1)&&(ls>PKD_MAX_CELL_SIZE))); /* this condition means the left child is not a bucket */
	    rc = ((nr > M)||((nr > 1)&&(rs>PKD_MAX_CELL_SIZE)));
	    if (rc && lc) {
		/* Allocate more stack if required */
		if ( s+1 >= ns ) {
		    assert( s+1 == ns );
		    ns += TEMP_S_INCREASE;
		    S = realloc(S,ns*sizeof(int));
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
		pkd->kdNodes[iRight].iLower = 0;
		++nBucket;
		}
	    else if (rc) {
		/*
		** Left must be a bucket in this case!
		*/
		iNode = iRight;   /* process upper subfile */
		pkd->kdNodes[iLeft].iLower = 0;
		++nBucket;
		}
	    else {
		/*
		** Both are buckets (we need to pop from the stack to get the next subfile.
		*/
		pkd->kdNodes[iLeft].iLower = 0;
		++nBucket;
		pkd->kdNodes[iRight].iLower = 0;
		++nBucket;
		if (s) iNode = S[--s];		/* pop tn */
		else break;
		}
	    }
	else {
	    /*
	    ** No nodes allocated, Change the bounds if needed!
	    */
	    if (d >= 0 && d < 3) pkd->kdNodes[iNode].bnd.fMax[d] *= 0.5;
	    if (nl > 0) {
		if (d >= 0 && d < 3) pkd->kdNodes[iNode].bnd.fCenter[d] -= pkd->kdNodes[iNode].bnd.fMax[d];
		MAXSIDE(pkd->kdNodes[iNode].bnd.fMax,ls);
		lc = ((nl > M)||((nl > 1)&&(ls>PKD_MAX_CELL_SIZE))); /* this condition means the node is not a bucket */
		if (!lc) {
		    pkd->kdNodes[iNode].iLower = 0;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    else {
		if (d >= 0 && d < 3) pkd->kdNodes[iNode].bnd.fCenter[d] += pkd->kdNodes[iNode].bnd.fMax[d];
		MAXSIDE(pkd->kdNodes[iRight].bnd.fMax,rs);
		rc = ((nr > M)||((nr > 1)&&(rs>PKD_MAX_CELL_SIZE)));
		if (!rc) {
		    pkd->kdNodes[iNode].iLower = 0;
		    ++nBucket;
		    if (s) iNode = S[--s];		/* pop tn */
		    else break;
		    }
		}
	    }
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
    iTemp = iStart;
    while (1) {
	p = pkdParticle(pkd,iTemp);
	pkdSaveParticle(pkd,p);
	i = iTemp;
	iNew = pkd->pLite[i].i;
	while (iNew != iTemp) {
	    pNew = pkdParticle(pkd,iNew);
	    iNewer = pkd->pLite[iNew].i;
	    pNewer = pkdParticle(pkd,iNewer);
	    /* Particles are being shuffled here in a non-linear order.
	    ** Being smart humans, we can tell the CPU where the next chunk
	    ** of data can be found.  The limit is 8 outstanding prefetches
	    ** (according to the Opteron Guide).
	    */
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
	    __builtin_prefetch((char *)(pkd->pLite+iNewer)
			       + offsetof(struct pLite,i), 1, 3 );
	    __builtin_prefetch((char *)(pNewer)+0,1,0);
#ifndef __ALTIVEC__
	    __builtin_prefetch((char *)(pNewer)+64,1,0);
#endif
#endif
	    pkdCopyParticle(pkd,p,pNew);
	    pkd->pLite[i].i = 0;
	    i = iNew;
	    p = pkdParticle(pkd,i);
	    iNew = pkd->pLite[i].i;
	    }
	pkdLoadParticle(pkd,p);
	pkd->pLite[i].i = 0;
	while (!pkd->pLite[iTemp].i) {
	    if (++iTemp == pkd->nLocal) return;
	    }
	}
    }

static double zeroV[3] = {0.0,0.0,0.0};

void Create(PKD pkd,int iNode,FLOAT diCrit2,double dTimeStamp) {
    PARTICLE *p;
    KDN *c = pkd->kdNodes;
    KDN *pkdn,*pkdl,*pkdu;
    MOMR mom;
    FLOAT m,fMass,fSoft,x,y,z,vx,vy,vz,ax,ay,az,ft,d2,d2Max,dih2,b;
    float *a;
    double *v;
    int pj,d,nDepth,ism;
    const int nMaxStackIncrease = 1;

    assert(pkd->oVelocity); /* Validate memory model */
    assert(pkd->oAcceleration); /* Validate memory model */

    nDepth = 1;
    while (1) {
	while (c[iNode].iLower) {
	    iNode = c[iNode].iLower;
	    ++nDepth;
	    /*
	    ** Is this the deepest in the tree so far? We might need to have more stack
	    ** elements for the tree walk!
	    */
	    if (nDepth > pkd->nMaxStack/*+1*/) {
		pkd->S = realloc(pkd->S,(pkd->nMaxStack+nMaxStackIncrease)*sizeof(CSTACK));
		assert(pkd->S != NULL);
		for (ism=pkd->nMaxStack;ism<(pkd->nMaxStack+nMaxStackIncrease);++ism) {
		    pkd->S[ism].Check = malloc(pkd->nMaxCheck*sizeof(CELT));
		    assert(pkd->S[ism].Check != NULL);
		    }
		pkd->nMaxStack += nMaxStackIncrease;
		}
	    }
	/*
	** Now calculate all bucket quantities!
	** This includes M,CoM,Moments and special
	** bounds and iMaxRung.
	*/
	pkdn = &c[iNode];
	pj = pkdn->pLower;
	p = pkdParticle(pkd,pj);
	for (d=0;d<3;++d) {
	    ft = p->r[d];
	    pkdn->bnd.fCenter[d] = ft;
	    pkdn->bnd.fMax[d] = ft;
	    }
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    for (d=0;d<3;++d) {
		ft = p->r[d];
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
	pj = pkdn->pLower;
	p = pkdParticle(pkd,pj);
	a = pkdAccel(pkd,p);
	m = pkdMass(pkd,p);
	fSoft = pkdSoft(pkd,p);
	if ( pkd->oVelocity ) v = pkdVel(pkd,p);
	else v = zeroV;
	fMass = m;
	dih2 = m/(fSoft*fSoft);
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
	for (++pj;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    a = pkdAccel(pkd,p);
	    m = pkdMass(pkd,p);
	    fSoft = pkdSoft(pkd,p);
	    if ( pkd->oVelocity ) v = pkdVel(pkd,p);
	    else v = zeroV;
	    fMass += m;
	    dih2 += m/(fSoft*fSoft);
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
	    }
	m = 1/fMass;
	pkdn->r[0] = m*x;
	pkdn->r[1] = m*y;
	pkdn->r[2] = m*z;
	pkdn->v[0] = m*vx;
	pkdn->v[1] = m*vy;
	pkdn->v[2] = m*vz;
	pkdn->a[0] = m*ax;
	pkdn->a[1] = m*ay;
	pkdn->a[2] = m*az;
	dih2 *= m;
	pkdn->fSoft2 = 1/dih2;
	/*
	** Now calculate the reduced multipole moment.
	*/
	momClearMomr(&pkdn->mom);
	d2Max = 0;
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    x = p->r[0] - pkdn->r[0];
	    y = p->r[1] - pkdn->r[1];
	    z = p->r[2] - pkdn->r[2];
	    if (pkdIsSrcActive(p,0,MAX_RUNG)) {
		m = pkdMass(pkd,p);
		d2 = momMakeMomr(&mom,m,x,y,z);
		momAddMomr(&pkdn->mom,&mom);
		}
	    else {
		d2 = x*x + y*y + z*z;
		}
	    /*
	    ** Update bounding ball and softened bounding ball.
	    */
	    d2Max = (d2 > d2Max)?d2:d2Max;
	    }
	/*
	** Now determine the opening radius for gravity.
	*/
	MAXSIDE(pkdn->bnd.fMax,b);
#ifdef CLASSICAL_FOPEN
#ifdef LOCAL_EXPANSION
	pkdn->fOpen = sqrt(FOPEN_FACTOR*d2Max*diCrit2);
#else
	pkdn->fOpen2 = FOPEN_FACTOR*d2Max*diCrit2;
#endif
#else
#ifdef LOCAL_EXPANSION
	pkdn->fOpen = b*sqrt(diCrit2);
	if (pkdn->fOpen < sqrt(d2Max)) pkdn->fOpen = sqrt(d2Max);
#else
	pkdn->fOpen2 = b*b*diCrit2;
	if (pkdn->fOpen2 < d2Max) pkdn->fOpen2 = d2Max;
#endif
#endif
	/*
	** Set the timestamp for the node.
	*/
	pkdn->dTimeStamp = dTimeStamp;
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
	    pkdn = &c[iNode];
	    pkdl = &c[pkdn->iLower];
	    pkdu = &c[pkdn->iLower + 1];
	    pkdCombineCells(pkdn,pkdl,pkdu);
	    pj = pkdn->pLower;
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
		MAXSIDE(pkdn->bnd.fMax,b);
#ifdef CLASSICAL_FOPEN
#ifdef LOCAL_EXPANSION
		pkdn->fOpen = sqrt(FOPEN_FACTOR*d2Max*diCrit2);
#else
		pkdn->fOpen2 = FOPEN_FACTOR*d2Max*diCrit2;
#endif
#else
#ifdef LOCAL_EXPANSION
		pkdn->fOpen = b*sqrt(diCrit2);
		if (pkdn->fOpen < sqrt(d2Max)) pkdn->fOpen = sqrt(d2Max);
#else
		pkdn->fOpen2 = b*b*diCrit2;
		if (pkdn->fOpen2 < d2Max) pkdn->fOpen2 = d2Max;
#endif
#endif
		}
	    else {
		CALCOPEN(pkdn,diCrit2);
		}
	    }
	++iNode;
	}
    }



void pkdCombineCells(KDN *pkdn,KDN *p1,KDN *p2) {
    MOMR mom;
    FLOAT m1,m2,x,y,z,ifMass;
    FLOAT r1[3],r2[3];
    int j;

    for (j=0;j<3;++j) {
	r1[j] = p1->r[j];
	r2[j] = p2->r[j];
	}
    if (p1->dTimeStamp > p2->dTimeStamp) {
	/*
	** Shift r2 to be time synchronous to p1->dTimeStamp.
	*/
	pkdn->dTimeStamp = p1->dTimeStamp;
	assert(0);  /* need to code the shift y*/
	}
    else if (p1->dTimeStamp < p2->dTimeStamp) {
	/*
	** Shift r1 to be time synchronous to p2->dTimeStamp.
	*/
	pkdn->dTimeStamp = p2->dTimeStamp;
	assert(0); /* need to code the shift */
	}
    else {
	/*
	** Both child cells are time synchronous so we don't need to
	** shift either of them and we can also use the timestamp of either.
	*/
	pkdn->dTimeStamp = p1->dTimeStamp;
	}
    m1 = p1->mom.m;
    m2 = p2->mom.m;
    ifMass = 1/(m1 + m2);
    /*
    ** In the case where a cell has all its particles source inactive mom.m == 0, which is ok, but we
    ** still need a reasonable center in order to define opening balls in the tree code.
    */
    if ( m1==0.0 || m2 == 0.0 ) {
	ifMass = 1.0;
	m1 = m2 = 0.5;
	}
    for (j=0;j<3;++j) {
	pkdn->r[j] = ifMass*(m1*r1[j] + m2*r2[j]);
	pkdn->v[j] = ifMass*(m1*p1->v[j] + m2*p2->v[j]);
	pkdn->a[j] = ifMass*(m1*p1->a[j] + m2*p2->a[j]);
	}
    pkdn->fSoft2 = 1.0/(ifMass*(m1/p1->fSoft2 + m2/p2->fSoft2));
    pkdn->uMinRung = p1->uMinRung < p2->uMinRung ? p1->uMinRung : p2->uMinRung;
    pkdn->uMaxRung = p1->uMaxRung > p2->uMaxRung ? p1->uMaxRung : p2->uMaxRung;
    pkdn->bDstActive = p1->bDstActive || p2->bDstActive;

    /*
    ** Now calculate the reduced multipole moment.
    ** Shift the multipoles of each of the children
    ** to the CoM of this cell and add them up.
    */
    pkdn->mom = p1->mom;
    x = r1[0] - pkdn->r[0];
    y = r1[1] - pkdn->r[1];
    z = r1[2] - pkdn->r[2];
    momShiftMomr(&pkdn->mom,x,y,z);
    mom = p2->mom;
    x = r2[0] - pkdn->r[0];
    y = r2[1] - pkdn->r[1];
    z = r2[2] - pkdn->r[2];
    momShiftMomr(&mom,x,y,z);
    momAddMomr(&pkdn->mom,&mom);
    BND_COMBINE(pkdn->bnd,p1->bnd,p2->bnd);
    }


void pkdVATreeBuild(PKD pkd,int nBucket,FLOAT diCrit2,double dTimeStamp) {
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

    Create(pkd,VAROOT,diCrit2,dTimeStamp);
    }


void pkdTreeBuild(PKD pkd,int nBucket,FLOAT diCrit2,KDN *pkdn,int bExcludeVeryActive,double dTimeStamp) {
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
    Create(pkd,ROOT,diCrit2,dTimeStamp);

    pkdStopTimer(pkd,0);
#ifdef USE_BSC
    MPItrace_event(10000, 0 );
#endif
    /*
    ** Finally activate a read only cache for remote access.
    */
    mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkd->nNodes);
    /*
    ** Copy the root node for the top-tree construction.
    */
    *pkdn = pkd->kdNodes[ROOT];
    }


void pkdDistribCells(PKD pkd,int nCell,KDN *pkdn) {
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
void pkdCalcRoot(PKD pkd,MOMC *pmom) {
    PARTICLE *p;
    FLOAT xr = pkd->kdTop[ROOT].r[0];
    FLOAT yr = pkd->kdTop[ROOT].r[1];
    FLOAT zr = pkd->kdTop[ROOT].r[2];
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
    pkd->momRoot = *pmom;
    }


void pkdTreeNumSrcActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi) {
    KDN *c = pkd->kdNodes;
    PARTICLE *p;
    int iNode,pj;

    iNode = ROOT;
    while (1) {
	while (c[iNode].iLower) {
	    iNode = c[iNode].iLower;
	    }
	/*
	** We have to test each particle of the bucket for activity.
	*/
	c[iNode].nActive = 0;
	for (pj=c[iNode].pLower;pj<=c[iNode].pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    if (pkdIsSrcActive(p,uRungLo,uRungHi)) ++c[iNode].nActive;
	}
	while (iNode & 1) {
	    iNode = c[iNode].iParent;
	    if (!iNode) return;	/* exit point!!! */
	    c[iNode].nActive = c[c[iNode].iLower].nActive + c[c[iNode].iLower + 1].nActive;
	}
	++iNode;
    }
}


void pkdBoundWalk(PKD pkd,BND *pbnd,uint8_t uRungLo,uint8_t uRungHi,uint32_t *pnActive,uint32_t *pnContained) {
    KDN *c = pkd->kdNodes;
    PARTICLE *p;
    double d;
    int iNode,pj;    

    *pnActive = 0;
    *pnContained = 0;
    iNode = ROOT;
    while (1) {
	d = fabs(pbnd->fCenter[0] - c[iNode].bnd.fCenter[0]) - pbnd->fMax[0];
	if (d - c[iNode].bnd.fMax[0] > 0) goto NoIntersect;
	else if (d + c[iNode].bnd.fMax[0] <= 0) {
	    d = fabs(pbnd->fCenter[1] - c[iNode].bnd.fCenter[1]) - pbnd->fMax[1];
	    if (d - c[iNode].bnd.fMax[1] > 0) goto NoIntersect;
	    else if (d + c[iNode].bnd.fMax[1] <= 0) {
		d = fabs(pbnd->fCenter[2] - c[iNode].bnd.fCenter[2]) - pbnd->fMax[2];
		if (d - c[iNode].bnd.fMax[2] > 0) goto NoIntersect;
		else if (d + c[iNode].bnd.fMax[2] <= 0) goto Contained;
		}
	    else {
		d = fabs(pbnd->fCenter[2] - c[iNode].bnd.fCenter[2]) - pbnd->fMax[2];
		if (d - c[iNode].bnd.fMax[2] > 0) goto NoIntersect;
		}
	    }
	else {
	    d = fabs(pbnd->fCenter[1] - c[iNode].bnd.fCenter[1]) - pbnd->fMax[1];
	    if (d - c[iNode].bnd.fMax[1] > 0) goto NoIntersect;
	    d = fabs(pbnd->fCenter[2] - c[iNode].bnd.fCenter[2]) - pbnd->fMax[2];
	    if (d - c[iNode].bnd.fMax[2] > 0) goto NoIntersect;
	    }	
	/*
	** We have an intersection to test!
	*/
	if (c[iNode].iLower) {
	    iNode = c[iNode].iLower;
	    continue;
	    }
	else {
	    /*
	    ** We have to test each active particle of the bucket for containment.
	    */
	    for (pj=c[iNode].pLower;pj<=c[iNode].pUpper;++pj) {
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
	*pnContained += (c[iNode].pUpper - c[iNode].pLower + 1);
	*pnActive += c[iNode].nActive;  /* this must be set with SrcActive for all cells first */
    NoIntersect:
	while (iNode & 1) {
	    iNode = c[iNode].iParent;
	    if (!iNode) return;    /* exit point */
	}
	++iNode;
    }
}
