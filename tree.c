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



void InitializeParticles(PKD pkd,int bExcludeVeryActive) {
    PLITE *p = pkd->pLite;
    PLITE t;
    int i,j;

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	for (j=0;j<3;++j) p[i].r[j] = pkd->pStore[i].r[j];
	p[i].i = i;
	p[i].uRung = pkd->pStore[i].iRung;
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
	    if ( p[i].uRung <= pkdRungVeryActive(pkd) ) ++i;
	    else break;
	    }
	while (i <= j) {
	    if ( p[j].uRung > pkdRungVeryActive(pkd) ) --j;
	    else break;
	    }
	if (i < j) {
	    t = p[i];
	    p[i] = p[j];
	    p[j] = t;
	    while (1) {
		while ((p[++i].uRung <= pkdRungVeryActive(pkd)));
		while (p[--j].uRung > pkdRungVeryActive(pkd));
		if (i < j) {
		    t = p[i];
		    p[i] = p[j];
		    p[j] = t;
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
	}
    /*
    ** Set up the root node.
    */
    pkd->kdNodes[ROOT].iLower = 0;
    pkd->kdNodes[ROOT].iParent = 0;
    pkd->kdNodes[ROOT].pLower = 0;
    pkd->kdNodes[ROOT].pUpper = pkd->nLocal - pkd->nVeryActive - 1;
    }

#define MIN_SRATIO    0.05

/*
** M is the bucket size.
*/
void BuildTemp(PKD pkd,int iNode,int M,int bSqueeze) {
    PLITE *p = pkd->pLite;
    PLITE t;
    FLOAT fSplit,sRatio;
    FLOAT fMin[3],fMax[3];
    FLOAT ls,rs;
    BND *pbnd;
    int *S;		/* this is the stack */
    int s,ns;
    int iLeft,iRight;
    int d,i,j;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;
    int nActive = 0;

    //assert(bSqueeze > 0);
    /*
    ** Allocate stack!
    */
    ns = floor(log(((double)(pkd->nLocal+1))/(M+1))/log(2.0));
    if (ns < 1) ns = 1;	/* want to allocate something! */	
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
    
#ifdef GASOLINE
    /*
    ** No support for gasoline here quite yet.
    */
    assert(0);
#endif
    /*
    ** First set up the root node.
    ** Allocate it, set it's bounds and pointers.
    ** Also find out how many active particles there are.
    */
    i = pkd->kdNodes[iNode].pLower;
    for (j=0;j<3;++j) {
	fMin[j] = p[i].r[j];
	fMax[j] = p[i].r[j];
	if (pkdIsRungActive(pkd,p[i].uRung)) ++nActive;
	}
    for (++i;i<=pkd->kdNodes[iNode].pUpper;++i) {
	for (j=0;j<3;++j) {
	    if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
	    else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
	    }
	if (pkdIsRungActive(pkd,p[i].uRung)) ++nActive;
	}
    for (j=0;j<3;++j) {
	pkd->kdNodes[iNode].bnd.fCenter[j] = 0.5*(fMax[j] + fMin[j]);
	pkd->kdNodes[iNode].bnd.fMax[j] = 0.5*(fMax[j] - fMin[j]);
	}
    if (pkd->kdNodes[iNode].pUpper - pkd->kdNodes[iNode].pLower + 1 <= M) 
	goto DonePart;

    while (1) {
	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/
	pbnd = &pkd->kdNodes[iNode].bnd;
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
    JumpInFromActiveInactive:
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
	    if (bSqueeze) {
		/*
		** Calculate bounds for these 2 new cells.
		*/
		i = pkd->kdNodes[iLeft].pLower;
		for (j=0;j<3;++j) {
		    fMin[j] = p[i].r[j];
		    fMax[j] = p[i].r[j];
		    }
		for (++i;i<=pkd->kdNodes[iLeft].pUpper;++i) {
		    for (j=0;j<3;++j) {
			if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
			else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
			}
		    }
		for (j=0;j<3;++j) {
		    pkd->kdNodes[iLeft].bnd.fCenter[j] = 0.5*(fMax[j] + fMin[j]);
		    pkd->kdNodes[iLeft].bnd.fMax[j] = 0.5*(fMax[j] - fMin[j]);
		    }
		i = pkd->kdNodes[iRight].pLower;
		for (j=0;j<3;++j) {
		    fMin[j] = p[i].r[j];
		    fMax[j] = p[i].r[j];
		    }
		for (++i;i<=pkd->kdNodes[iRight].pUpper;++i) {
		    for (j=0;j<3;++j) {
			if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
			else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
			}
		    }
		for (j=0;j<3;++j) {
		    pkd->kdNodes[iRight].bnd.fCenter[j] = 0.5*(fMax[j] + fMin[j]);
		    pkd->kdNodes[iRight].bnd.fMax[j] = 0.5*(fMax[j] - fMin[j]);
		    }
		}
	    else {
		pkd->kdNodes[iLeft].bnd = pkd->kdNodes[iNode].bnd;
		pkd->kdNodes[iRight].bnd = pkd->kdNodes[iLeft].bnd;
		if (d >= 0 && d < 3) {    /* otherwise we have done some other kind of splitting and we don't cut the bounds */
		  pkd->kdNodes[iLeft].bnd.fMax[d] *= 0.5;
		  pkd->kdNodes[iLeft].bnd.fCenter[d] -= pkd->kdNodes[iLeft].bnd.fMax[d];
		  pkd->kdNodes[iRight].bnd.fCenter[d] += pkd->kdNodes[iRight].bnd.fMax[d];
		}
		}
	    MAXSIDE(pkd->kdNodes[iLeft].bnd.fMax,ls);
	    MAXSIDE(pkd->kdNodes[iRight].bnd.fMax,rs);
	    }
	else {
	    /*
	    ** No nodes allocated, this can't happen if we
	    ** are squeezing (or shouldn't).
	    ** Change the bounds!
	    */
	    assert(!bSqueeze);
	    if (d >= 0 && d < 3) pbnd->fMax[d] *= 0.5;
	    if (nl > 0) {
	      if (d >= 0 && d < 3) pbnd->fCenter[d] -= pbnd->fMax[d];
	      iLeft = iNode;
	      MAXSIDE(pkd->kdNodes[iLeft].bnd.fMax,ls);
	      rs = 0.0;
	    }
	    else {
	      if (d >= 0 && d < 3) pbnd->fCenter[d] += pbnd->fMax[d];
	      iRight = iNode;
	      MAXSIDE(pkd->kdNodes[iRight].bnd.fMax,rs);
	      ls = 0.0;
	    }
	    }
	/*
	** Now figure out which subfile to process next.
	*/
	lc = ((nl > M)||(ls>1e-2));
	rc = ((nr > M)||(rs>1e-2));
	MAXSIDE(pkd->kdNodes[iLeft].bnd.fMax,ls);
	MAXSIDE(pkd->kdNodes[iRight].bnd.fMax,rs);

	if (rc && lc) {
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
	    if (lc) {
		iNode = iLeft;		/* process lower subfile */
		}
	    else if (nl > 0) {
		pkd->kdNodes[iLeft].iLower = 0;
		++nBucket;
		}
	    if (rc) {
		iNode = iRight;		/* process upper subfile */
		}
	    else if (nr > 0) {
		pkd->kdNodes[iRight].iLower = 0;
		++nBucket;
		}
	    }
	if ( !rc && !lc ) {
	    if (s) iNode = S[--s];		/* pop tn */
	    else break;
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
    PARTICLE Temp;
    int i,iNew,iNewer,iTemp;

    /*
    ** Now we move the particles in one go using the temporary
    ** particles which have been shuffled.
    */
    iTemp = iStart;
    while (1) {
	Temp = pkd->pStore[iTemp];
	i = iTemp;
	iNew = pkd->pLite[i].i;
	iNewer = pkd->pLite[iNew].i;
	while (iNew != iTemp) {
	    /* Particles are being shuffled here in a non-linear order.
	    ** Being smart humans, we can tell the CPU where the next chunk
	    ** of data can be found.  The limit is 8 outstanding prefetches
	    ** (according to the Opteron Guide).
	    */
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
	    __builtin_prefetch((char *)(pkd->pLite+iNewer)
		+ offsetof(struct pLite,i), 1, 3 );
	    __builtin_prefetch((char *)(pkd->pStore+iNewer)+0,1,0);
#ifndef __ALTIVEC__
	    __builtin_prefetch((char *)(pkd->pStore+iNewer)+64,1,0);
#endif
	    __builtin_prefetch((char *)(pkd->pStore+iNewer)+128,1,0);
#ifndef __ALTIVEC__
	    __builtin_prefetch((char *)(pkd->pStore+iNewer)+192,1,0);

#endif
#endif
#ifdef xx__SSE__
	    _mm_prefetch((char *)(pkd->pLite+iNewer)+offsetof(struct pLite,i),
			 _MM_HINT_T0 );
	    _mm_prefetch((char *)(pkd->pStore+iNewer)+0,_MM_HINT_NTA);
	    _mm_prefetch((char *)(pkd->pStore+iNewer)+64,_MM_HINT_NTA);
	    _mm_prefetch((char *)(pkd->pStore+iNewer)+128,_MM_HINT_NTA);
	    _mm_prefetch((char *)(pkd->pStore+iNewer)+192,_MM_HINT_NTA);
#endif
#ifdef xx__ALTIVEC__
	__asm__ __volatile__ ("dcbt 0, %0"::"r"((char *)(pkd->pLite+iNewer)+offsetof(struct pLite,i)));
	__asm__ __volatile__ ("dcbt 0, %0"::"r"((char *)(pkd->pStore+iNewer)+0));
	__asm__ __volatile__ ("dcbt 0, %0"::"r"((char *)(pkd->pStore+iNewer)+128));
#endif
	    pkd->pStore[i] = pkd->pStore[iNew];
	    pkd->pLite[i].i = 0;
	    i = iNew;
	    iNew = pkd->pLite[i].i;
	    iNewer = pkd->pLite[iNew].i;
	    }
	pkd->pStore[i] = Temp;
	pkd->pLite[i].i = 0;
	while (!pkd->pLite[iTemp].i) {
	    if (++iTemp == pkd->nLocal) return;
	    }
	}
    }


void Create(PKD pkd,int iNode,FLOAT diCrit2,double dTimeStamp,int bTempBound) {
    PARTICLE *p = pkd->pStore;
    KDN *c = pkd->kdNodes;
    KDN *pkdn,*pkdl,*pkdu;
    MOMR mom;
    FLOAT m,fMass,x,y,z,vx,vy,vz,ax,ay,az,ft,d2,d2Max,dih2,b;
    int pj,d,nDepth,ism;
    const int nMaxStackIncrease = 1;

    nDepth = 1;
    while (1) {
	while (c[iNode].iLower) {
	    iNode = c[iNode].iLower;
	    ++nDepth;
	    /*
	    ** Is this the deepest in the tree so far? We might need to have more stack
	    ** elements for the tree walk!
	    */
	    if (nDepth > pkd->nMaxStack+1) {
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
	if (!bTempBound) {
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
	ax = m*p[pj].a[0];
	ay = m*p[pj].a[1];
	az = m*p[pj].a[2];
	pkdn->iActive = p[pj].iActive;
	pkdn->uMinRung = pkdn->uMaxRung = p[pj].iRung;
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
	    ax += m*p[pj].a[0];
	    ay += m*p[pj].a[1];
	    az += m*p[pj].a[2];
	    pkdn->iActive |= p[pj].iActive;
	    if ( p[pj].iRung > pkdn->uMaxRung ) pkdn->uMaxRung = p[pj].iRung;
	    if ( p[pj].iRung < pkdn->uMinRung ) pkdn->uMinRung = p[pj].iRung;
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
	    pkdn = &c[iNode];
	    pkdl = &c[pkdn->iLower];
	    pkdu = &c[pkdn->iLower + 1];
	    pkdCombineCells(pkdn,pkdl,pkdu,!bTempBound);
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



void pkdCombineCells(KDN *pkdn,KDN *p1,KDN *p2,int bCombineBound) {
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
      ** Both child cells are tiume synchronous so we don't need to 
      ** shift either of them and we can also use the timestamp of either.
      */
      pkdn->dTimeStamp = p1->dTimeStamp;
    }
    m1 = p1->mom.m;
    m2 = p2->mom.m;
    ifMass = 1/(m1 + m2);
    for (j=0;j<3;++j) {
      pkdn->r[j] = ifMass*(m1*r1[j] + m2*r2[j]);
      pkdn->v[j] = ifMass*(m1*p1->v[j] + m2*p2->v[j]);
      pkdn->a[j] = ifMass*(m1*p1->a[j] + m2*p2->a[j]);
    }
    pkdn->fSoft2 = 1.0/(ifMass*(m1/p1->fSoft2 + m2/p2->fSoft2));
    pkdn->iActive = (p1->iActive | p2->iActive);
    pkdn->uMinRung = p1->uMinRung < p2->uMinRung ? p1->uMinRung : p2->uMinRung;
    pkdn->uMaxRung = p1->uMaxRung > p2->uMaxRung ? p1->uMaxRung : p2->uMaxRung;

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
    /*
    ** Combine the bounds!
    */
    if (bCombineBound) {
	BND_COMBINE(pkdn->bnd,p1->bnd,p2->bnd);
        }
    }


void pkdVATreeBuild(PKD pkd,int nBucket,FLOAT diCrit2,int bSqueeze,double dTimeStamp) {
    int i,j,iStart;

    iStart = pkd->nLocal - pkd->nVeryActive;
    /*
    ** First initialize the very active temporary particles.
    */
    for (i=iStart;i<pkd->nLocal;++i) {
	for (j=0;j<3;++j) pkd->pLite[i].r[j] = pkd->pStore[i].r[j];
	pkd->pLite[i].i = i;
	pkd->pLite[i].uRung = pkd->pStore[i].iRung;
	}
    /*
    ** Then clear the VA tree by setting the node index back to one node past the end
    ** of the non VA tree.
    */
    pkd->nNodes = pkd->nNonVANodes;
    BuildTemp(pkd,VAROOT,nBucket,bSqueeze);

    ShuffleParticles(pkd,iStart);
    
    Create(pkd,VAROOT,diCrit2,dTimeStamp,bSqueeze);
    }


void pkdTreeBuild(PKD pkd,int nBucket,FLOAT diCrit2,KDN *pkdn,int bSqueeze,int bExcludeVeryActive,double dTimeStamp) {
    int iStart,nNodesEst;

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

    InitializeParticles(pkd,bExcludeVeryActive);

    BuildTemp(pkd,ROOT,nBucket,bSqueeze);
    if (bExcludeVeryActive) {
	pkd->nNonVANodes = pkd->nNodes;
	}
    else {
	pkd->nNodesFull = pkd->nNodes;
	}
    iStart = 0;
    ShuffleParticles(pkd,iStart);
    Create(pkd,ROOT,diCrit2,dTimeStamp,bSqueeze);

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


void pkdDistribBoundBall(PKD pkd,int nCell,BND *pbnd) {
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
