#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *psdtree_module_id = "$Id";

#ifdef USE_PSD

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

#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif
#include <sys/time.h>



void psdInitializeParticles(PKD pkd, BND *pbnd) {
    PLITE *pLite = pkd->pLite;
    PLITE t;
    PARTICLE *p;
    KDN *pNode;
    int i,j;

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) pLite[i].r[j] = p->r[j];
	for (j=0;j<3;++j) pLite[i].v[j] = pkdVel(pkd,p)[j];
        pLite[i].fMass = pkdMass(pkd,p);
	pLite[i].i = i;
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
    /*
    ** Set up the root node.
    */
    pNode = pkdTreeNode(pkd,ROOT);
    pNode->iLower = 0;
    pNode->iParent = 0;
    pNode->pLower = 0;
    pNode->pUpper = pkd->nLocal - pkd->nVeryActive - 1;
    pBND bnd = pkdNodeBnd(pkd,pNode);
    for (j=0;j<6;++j) {
	bnd.fCenter[j] = pbnd->fCenter[j];
	bnd.fMax[j]    = pbnd->fMax[j];
	}
    }

#if 1

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
    int32_t key;
    float fMass;
} KEY;

int key_compar(const void *a0, const void *b0) {
    KEY *a = (KEY *)a0;
    KEY *b = (KEY *)b0;
    return a->key - b->key;
    }

static inline float E(KEY *keys, int N, float M) {
    int i;
    float s = 0;
    int32_t key = keys[0].key;
    float   p   = keys[0].fMass;
    assert(M > 0);
    for (i=1; i < N; i++) {
        if (keys[i].key != key) {
            p /= M;
            s -= p * log10(p);
            p = 0;
            key = keys[i].key;
            }

        p += keys[i].fMass;
        }

    s -= p * log10(p);

    return s;
}

/*
** M is the bucket size.
** This function assumes that the root node is correctly set up (particularly the bounds).
*/

#if 0
static int max_bnd(const double *fMax)
{
    int i;
    int d = 0;
    double max = fMax[0];
    for (i=1; i < 6; i++)
        if (fMax[i] > max) { max = fMax[i]; d = i; }
    return d;
}


static int max_side(const double *fMax)
{
    return 2 * max_bnd(fMax);
}
#endif

#define SAVE_BOUNDS(node, bnd) do {\
    PSMETRIC *b = pkdPsMetric(pkd, pkdParticle(pkd,p[(node)->pLower].i)); \
    for (j=0;j<6;++j) { \
        b->scale[j] = bnd.fMax[j]; \
        } \
} while (0)

#define CUBICCELLS 0
#define CUBICCELLS2 1
#define LONGDIMSPLIT 0
#define RHOSPLIT 1
#define MEANSPLIT 2
#define ROTATESPLIT 3
#define SPLIT RHOSPLIT
//#define SPLIT LONGDIMSPLIT
#define TEMP_S_INCREASE 100
void BuildPsdTemp(PKD pkd,int iNode,int M, int maxNb) {
    PLITE6 *p = UNION_CAST(pkd->pLite, PLITE*, PLITE6*);
    PLITE6 t;
    KDN *pLeft, *pRight;
    KDN *pNode = pkdTreeNode(pkd,iNode);
    pBND bnd = pkdNodeBnd(pkd, pNode);
    FLOAT fSplit;
    FLOAT ls,rs;
    int *S;		/* this is the stack */
    int *D;             /* stack for the last cut dimension */
    int iLeft,iRight;
    int d,i,j;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;
    int Nb;
    int backup = 0;
    int lastd;



#if SPLIT==RHOSPLIT
#if CUBICCELLS
    float *rho[2];
    int rho_size = maxNb*maxNb*maxNb;
    rho[0] = malloc(rho_size * sizeof(float)); assert(rho[0] != NULL);
    rho[1] = malloc(rho_size * sizeof(float)); assert(rho[1] != NULL);
    memset(rho[0], 0, rho_size * sizeof(float));
    memset(rho[1], 0, rho_size * sizeof(float));
#elif CUBICCELLS2
    int Ntotal = pNode->pUpper - pNode->pLower + 1;

    KEY *xkeys = malloc((Ntotal+1) * sizeof(*xkeys));
    KEY *vkeys = malloc((Ntotal+1) * sizeof(*vkeys));
    
#else
    float *rho[6];
    size_t rho_size = maxNb;
    for (i=0; i<6; i++)
        rho[i] = malloc(rho_size * sizeof(float));
#endif
#endif

    assert(maxNb > 0);
    if (CUBICCELLS2) assert(maxNb <= (1<<10));


    printf("BuildPsdTemp M=%i\n", M);

    /*
    ** Allocate stack!
    */
    NEW_STACK(S, TEMP_S_INCREASE);
    NEW_STACK(D, TEMP_S_INCREASE);

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

    assert( bnd.fMax[0] > 0.0 
         || bnd.fMax[1] > 0.0 
         || bnd.fMax[2] > 0.0 
	     || bnd.fMax[3] > 0.0 
         || bnd.fMax[4] > 0.0 
         || bnd.fMax[5] > 0.0 
            );

    assert(STACK_EMPTY(S));
    PUSH(S, iNode); PUSH(D, 0);

    while (!STACK_EMPTY(S)) {

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

        lastd = POP(D);


	pNode = pkdTreeNode(pkd,iNode);
        bnd   = pkdNodeBnd(pkd, pNode);

        int Npart = pNode->pUpper - pNode->pLower + 1;

        //fprintf(stderr, "%i %i\n", iNode,Npart);

        if (!( bnd.fMax[0] > 0.0 &&
               bnd.fMax[1] > 0.0 &&
               bnd.fMax[2] > 0.0 &&
               bnd.fMax[3] > 0.0 &&
               bnd.fMax[4] > 0.0 &&
               bnd.fMax[5] > 0.0 
               ))
            fprintf(stderr, "%g %g %g %g %g %g\n",
                bnd.fMax[0],
                bnd.fMax[1],
                bnd.fMax[2],
                bnd.fMax[3],
                bnd.fMax[4],
                bnd.fMax[5]);

        assert( bnd.fMax[0] > 0.0 
             || bnd.fMax[1] > 0.0 
             || bnd.fMax[2] > 0.0
             || bnd.fMax[3] > 0.0 
             || bnd.fMax[4] > 0.0 
             || bnd.fMax[5] > 0.0 
                );

	/*
	** Begin new stage!
	** Calculate the appropriate fSplit.
	** Pick longest dimension and split it in half.
	*/
#if SPLIT==LONGDIMSPLIT
        d = max_bnd(pNode->bnd.fMax);
	fSplit = pNode->bnd.fCenter[d];
#elif SPLIT==ROTATESPLIT
        d = pNode->bnd.lastd = (pNode->bnd.lastd+1) % 6;
	fSplit = pNode->bnd.fCenter[d];
#elif SPLIT==MEANSPLIT
        fSplit = 0;
        if (d < 3) {
            for (i=pNode->pLower; i <= pNode->pUpper; i++) {
                fSplit += p[i].r[d];
            }
        }
        else {
            for (i=pNode->pLower; i <= pNode->pUpper; i++) {
                fSplit += p[i].v[d-3];
            }
        }

        fSplit /= pNode->pUpper - pNode->pLower + 1;
#elif SPLIT==RHOSPLIT

#if CUBICCELLS
        FLOAT dx_inv[6];
        float e[2];
        float Mtotal = 0;
        int k;
        FLOAT edge[6];
        int b[6];

        Nb = Npart;
        if (Nb > maxNb) Nb = maxNb;
        assert(Nb > 1);

        int Nb3 = Nb*Nb*Nb;

        for (i=0; i < 6; i++)
        {
            dx_inv[i] = (FLOAT)Nb / (2*bnd.fMax[i]);
            edge[i] = bnd.fCenter[i]-bnd.fMax[i]; 
        }

        for (j=pNode->pLower; j <= pNode->pUpper; j++) 
        {
            //float mass = pkdMass(pkd, pkdParticle(pkd,p[j].i));
            float mass = p[j].fMass;
            Mtotal += mass;


            b[0] = (int)((p[j].r[0]-edge[0]) * dx_inv[0]); b[0] -= (b[0] == Nb); //assert(0 <= i0 && i0 < Nb);
            b[1] = (int)((p[j].r[1]-edge[1]) * dx_inv[1]); b[1] -= (b[1] == Nb); //assert(0 <= i1 && i1 < Nb);
            b[2] = (int)((p[j].r[2]-edge[2]) * dx_inv[2]); b[2] -= (b[2] == Nb); //assert(0 <= i2 && i2 < Nb);

            b[3] = (int)((p[j].r[3]-edge[3]) * dx_inv[3]); b[3] -= (b[3] == Nb); //assert(0 <= i3 && i3 < Nb);
            b[4] = (int)((p[j].r[4]-edge[4]) * dx_inv[4]); b[4] -= (b[4] == Nb); //assert(0 <= i4 && i4 < Nb);
            b[5] = (int)((p[j].r[5]-edge[5]) * dx_inv[5]); b[5] -= (b[5] == Nb); //assert(0 <= i5 && i5 < Nb);

            rho[0][b[0] + Nb * (b[1] + Nb*b[2])] += mass;
            rho[1][b[3] + Nb * (b[4] + Nb*b[5])] += mass;
        }

        e[0] = 0;
        e[1] = 0;
        for (i=0; i < Nb3; i++)
            if (rho[0][i]) { rho[0][i] /= Mtotal; e[0] -= rho[0][i] * log10(rho[0][i]); rho[0][i] = 0;}

        for (i=0; i < Nb3; i++)
            if (rho[1][i]) { rho[1][i] /= Mtotal; e[1] -= rho[1][i] * log10(rho[1][i]); rho[1][i] = 0;}

        int dmin = 0;
        if ((e[0] >  e[1])
        ||  (e[0] == e[1] && lastd < 3))
        //||  (e[0] == e[1] && pNode->bnd.lastd < 3))
        {
            dmin = 3;
        }

        d = max_bnd_range(bnd.fMax, dmin, dmin+3);
        fSplit = bnd.fCenter[d];
#elif CUBICCELLS2
        FLOAT dx_inv[6];
        FLOAT edge[6];
        float e[2];
        float Mtotal = 0;
        int b[6];

        Nb = Npart;
        if (Nb > maxNb) Nb = maxNb;
        assert(Nb > 1);

        for (i=0; i < 6; i++)
        {
            dx_inv[i] = (FLOAT)Nb / (2*bnd.fMax[i]);
            edge[i]   = bnd.fCenter[i]-bnd.fMax[i]; 
        }

        for (j=pNode->pLower; j <= pNode->pUpper; j++) 
        {
            float mass = p[j].fMass;
            Mtotal += mass;

            for (i=0; i < 6; i++)
                b[i] = (int)((p[j].r[i]-edge[i]) * dx_inv[i]); b[i] -= (b[i] == Nb); //assert(0 <= b[i] && b[i] < Nb)

            xkeys[j].key = b[0] + Nb * (b[1] + Nb*b[2]);  xkeys[j].fMass = mass;
            vkeys[j].key = b[3] + Nb * (b[4] + Nb*b[5]);  vkeys[j].fMass = mass;
        }

        qsort(xkeys + pNode->pLower, Npart, sizeof(*xkeys), key_compar);
        qsort(vkeys + pNode->pLower, Npart, sizeof(*vkeys), key_compar);

        e[0] = E(xkeys + pNode->pLower, Npart, Mtotal);
        e[1] = E(vkeys + pNode->pLower, Npart, Mtotal);

        int dmin = 0;
        if ((e[0] >  e[1])
        ||  (e[0] == e[1] && lastd < 3))
        //||  (e[0] == e[1] && pNode->bnd.lastd < 3))
        {
            dmin = 3;
        }

        d = max_bnd_range(bnd.fMax, dmin, dmin+3);
        if (d == -1)
        {
            if (dmin==3)
                d = max_bnd_range(bnd.fMax, 0, 3);
            else
                d = max_bnd_range(bnd.fMax, 3, 6);
            assert(d != -1);
        }
        fSplit = bnd.fCenter[d];
#else
        FLOAT dx[6];
        float e[6];
        float Mtotal = 0;

        Nb = Npart;
        assert(Nb > 1);
        if (Nb > maxNb) Nb = maxNb;

        for (i=0; i < 6; i++)
        {
            memset(rho[i], 0, Nb * sizeof(float));
            dx[i] = 2*pNode->bnd.fMax[i] / (FLOAT)Nb;
            assert(dx[i] > 0);
        }

        for (j=pNode->pLower; j <= pNode->pUpper; j++) 
        {
            float mass = pkdMass(pkd, pkdParticle(pkd,p[j].i));
            Mtotal += mass;
            for (i=0; i < 3; i++) 
            {
                FLOAT edge = pNode->bnd.fCenter[i]-pNode->bnd.fMax[i];
                assert(p[j].r[i] >= edge);
                int rx = (int)((p[j].r[i]-edge) / dx[i]); if (rx == Nb) rx = Nb-1;
                if (rx < 0 || rx >= Nb)
                        fprintf(stderr, "%.15g %.15g %g %i %i %.15g\n", p[j].r[i], edge, dx[i], Nb, rx, pNode->bnd.fCenter[i]+pNode->bnd.fMax[i]);
                assert(rx >= 0);
                assert(rx < Nb);
                rho[i][rx] += mass;
            }
            for (i=3; i < 6; i++) 
            {
                FLOAT edge = pNode->bnd.fCenter[i]-pNode->bnd.fMax[i];
                int rx = (int)((p[j].v[i-3]-edge) / dx[i]); if (rx == Nb) rx = Nb-1;
                assert(rx >= 0);
                assert(rx < Nb);
                rho[i][rx] += mass;
            }
        }
        for (i=0; i < 6; i++)
        {
            e[i] = 0;
            for (j=0; j < Nb; j++)
            {
                assert(!isnan(rho[i][j]));
                rho[i][j] /= Mtotal;
                //fprintf(stderr, "rho[%i][%i]=%g\n", i, j, rho[i][j]);
                assert(rho[i][j] >= 0);
                if (rho[i][j])
                    e[i] -= rho[i][j] * log10(rho[i][j]);
                assert(!isnan(e[i]));
            }
        }

        int allsame = 1;
        float emin = e[0];
        d = 0;
        for (i=1; i < 6; i++)
        {
            if (0 < e[i] && e[i] < emin) { emin=e[i]; d = i; }
            else if (e[i] != emin) allsame = 0;
        }

        if (allsame || emin == 0)
            d = (pNode->bnd.lastd+1) % 6;

        fSplit = pNode->bnd.fCenter[d];
#endif

#else
#error "Unknown split defined."
#endif

        //pNode->bnd.lastd = d;

        //printf("%g %g %g %g\n", fSplit, pNode->bnd.fCenter[d]-pNode->bnd.fMax[d], pNode->bnd.fCenter[d], pNode->bnd.fCenter[d]+pNode->bnd.fMax[d]);
        //printf("      %g %g\n", fSplit, pNode->bnd.fMax[d]);
        //assert((pNode->bnd.fCenter[d]-pNode->bnd.fMax[d]) <= fSplit && fSplit <= (pNode->bnd.fCenter[d]+pNode->bnd.fMax[d]));
	/*
	** Now start the partitioning of the particles about
	** fSplit on dimension given by d.
	*/
	i = pNode->pLower;
	j = pNode->pUpper;
        PARTITION(i<j,i<=j,
               ++i,--j,
               (t = p[i], p[i]=p[j], p[j] = t),
               p[i].r[d] < fSplit, p[j].r[d] >= fSplit);

	nl = i - pNode->pLower;
	nr = pNode->pUpper - i + 1;

        //printf("%i %.15f  [%g %g] (%i %i) %i %i %s\n", d, fSplit, bnd.fCenter[d]-bnd.fMax[d], bnd.fCenter[d]+bnd.fMax[d], nl, nr, iNode, pNode->iParent, backup == 0 ? "" : "*");

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

            pBND rbnd = pkdNodeBnd(pkd, pRight);
            pBND lbnd = pkdNodeBnd(pkd, pLeft);

            //pRight->bnd.lastd = d;
            //pLeft->bnd.lastd = d;

            //fprintf(stderr, "%i %i %i\n", iNode, iLeft, iRight);

	    /*
	    ** Now deal with the bounds.
            **
            ** ADD SHRINK WRAPPING -- jpc 27.3.2010
            **
	    */
	    for (j=0;j<6;++j) {
		if (j == d) {
		    rbnd.fMax[j]    = lbnd.fMax[j] = 0.5*bnd.fMax[j];
		    lbnd.fCenter[j] = bnd.fCenter[j] - lbnd.fMax[j];
		    rbnd.fCenter[j] = bnd.fCenter[j] + rbnd.fMax[j];
		    }
		else {
		    lbnd.fCenter[j] = bnd.fCenter[j];
		    lbnd.fMax[j]    = bnd.fMax[j];
		    rbnd.fCenter[j] = bnd.fCenter[j];
		    rbnd.fMax[j]    = bnd.fMax[j];
		    }
                //assert(lbnd.fMax[j] > 0.0);
                //assert(rbnd.fMax[j] > 0.0);
		}

            ls = max_side(lbnd.fMax);     // MAXSIDE(pLeft.bnd.fMax,ls);
            rs = max_side(rbnd.fMax);     // MAXSIDE(pRight.bnd.fMax,rs);
	    /*
	    ** Now figure out which subfile to process next.
	    */
	    lc = ((nl > M)||((nl > 1)&&(ls>PKD_MAX_CELL_SIZE))); /* this condition means the left child is not a bucket */
	    rc = ((nr > M)||((nr > 1)&&(rs>PKD_MAX_CELL_SIZE)));
	    if (rc && lc) {
                assert(backup == 0);
                EXTEND_STACK(S); /* Allocate more stack if required */
                EXTEND_STACK(D); /* Allocate more stack if required */

#if 1
                if (nr > nl) {
                    PUSH(S, iRight); PUSH(D, d);
                    PUSH(S, iLeft);  PUSH(D, d);
                    }
                else {
                    PUSH(S, iLeft);  PUSH(D, d);
                    PUSH(S, iRight); PUSH(D, d);
                    }
#else
		if (nr > nl) {
		    S[s++] = iRight;	/* push tr */
		    iNode = iLeft;		/* process lower subfile */
		    }
		else {
		    S[s++] = iLeft;	/* push tl */
		    iNode = iRight;		/* process upper subfile */
		    }
#endif
		}
	    else if (lc) {
		/*
		** Right must be a bucket in this case!
		*/
                assert(backup == 0);
                PUSH(S, iLeft); PUSH(D, d);

                if (nr > 1) {
                    backup = pkd->nNodes; //fprintf(stderr, "backup [%i]\n", backup); 
                    PUSH(S, 0);
                    PUSH(S, iRight); PUSH(D, d);
                    ++nBucket;
                    }
                else {
                    pRight->iLower = 0;
                    for (j=0;j<6;++j) assert(!isinf(rbnd.fMax[j]));
                    SAVE_BOUNDS(pRight,rbnd);
                    }
		}
	    else if (rc) {
		/*
		** Left must be a bucket in this case!
		*/
                assert(backup == 0);
                PUSH(S, iRight); PUSH(D, d);

                if (nl > 1) {
                    backup = pkd->nNodes; //fprintf(stderr, "backup [%i]\n", backup); 
                    PUSH(S, 0);
                    PUSH(S, iLeft); PUSH(D, d);
                    ++nBucket;
                    }
                else {
                    pLeft->iLower = 0;
                    for (j=0;j<6;++j) assert(!isinf(lbnd.fMax[j]));
                    SAVE_BOUNDS(pLeft,lbnd);
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

                if (nr > 1) { PUSH(S, iRight); PUSH(D, d); }
                if (nl > 1) { PUSH(S, iLeft);  PUSH(D, d); }

                if (nr == 1) {pRight->iLower=0; SAVE_BOUNDS(pRight,rbnd); for (j=0;j<6;++j) assert(!isinf(lbnd.fMax[j])); }
                if (nl == 1) {pLeft->iLower=0;  SAVE_BOUNDS(pLeft,lbnd); for (j=0;j<6;++j) assert(!isinf(lbnd.fMax[j]));}
		}
	    }
	else {
	    /*
	    ** No nodes allocated, Change the bounds if needed!
	    */

            assert(nr == 0 || nl == 0);
            int n = nr + nl; /* one of these will be zero */
            lc = rc = 0;
            if (0 <= d && d < 6) {
                float xx = bnd.fMax[d];
                bnd.fMax[d] *= 0.5;
                if (!(bnd.fMax[d] > 0))
                    fprintf(stderr, "%f\n", d);
                assert(bnd.fMax[d] > 0);
	        if (nl > 0) bnd.fCenter[d] -= bnd.fMax[d];
	        else        bnd.fCenter[d] += bnd.fMax[d];
		}

            //ls = max_side(pNode->bnd.fMax); // MAXSIDE(pNode->bnd.fMax,ls);
	    //lc = ((n > M)||((n > 1)&&(ls>PKD_MAX_CELL_SIZE))); /* this condition means the node is not a bucket */

            if (n > 1) { PUSH(S, iNode); PUSH(D, d); }

            if (n == 1) {
                for (j=0;j<6;++j) assert(!isinf(bnd.fMax[j]));
                SAVE_BOUNDS(pNode,bnd);
                pNode->iLower = 0;
                if (backup == 0) ++nBucket;
                }
	    }
	}
DonePart:
    FREE_STACK(S);
    FREE_STACK(D);
#if SPLIT==RHOSPLIT
#if CUBICCELLS
    free(rho[0]);
    free(rho[1]);
#elif CUBICCELLS2
    free(xkeys);
    free(vkeys);
#else
    for (i=0; i<6; i++)
        free(rho[i]);
#endif
#endif
    }
#endif

#if 0
static double zeroV[3] = {0.0,0.0,0.0};
static float  zeroF[3] = {0.0,0.0,0.0};
#endif

#if 0
void pkdCombineCells(PKD pkd,KDN *pkdn,KDN *p1,KDN *p2) {
    MOMR mom;
    SPHBNDS *b1,*b2,*bn;
    FLOAT m1,m2,x,y,z,ifMass;
    FLOAT r1[3],r2[3];
    int j;

    for (j=0;j<3;++j) {
	r1[j] = p1->r[j];
	r2[j] = p2->r[j];
	}
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
    for (j=0;j<3;++j) {
	pkdn->r[j] = ifMass*(m1*r1[j] + m2*r2[j]);
	if (pkd->oNodeVelocity)
	    pkdNodeVel(pkd,pkdn)[j]
		= ifMass*(m1*pkdNodeVel(pkd,p1)[j] + m2*pkdNodeVel(pkd,p2)[j]);
	if (pkd->oNodeAcceleration)
	    pkdNodeAccel(pkd,pkdn)[j]
		= ifMass*(m1*pkdNodeAccel(pkd,p1)[j] + m2*pkdNodeAccel(pkd,p2)[j]);
	}
    if(p1->fSoft2 == 0.0 || p2->fSoft2 == 0.0)
	pkdn->fSoft2 = 0.0;
    else
	pkdn->fSoft2 = 1.0/(ifMass*(m1/p1->fSoft2 + m2/p2->fSoft2));
    pkdn->uMinRung = p1->uMinRung < p2->uMinRung ? p1->uMinRung : p2->uMinRung;
    pkdn->uMaxRung = p1->uMaxRung > p2->uMaxRung ? p1->uMaxRung : p2->uMaxRung;
    pkdn->bDstActive = p1->bDstActive || p2->bDstActive;
    if (0xffffffffu - p1->nActive < p2->nActive) pkdn->nActive = 0xffffffffu; 
    else pkdn->nActive = p1->nActive + p2->nActive;
    /*
    ** Now calculate the reduced multipole moment.
    ** Shift the multipoles of each of the children
    ** to the CoM of this cell and add them up.
    */
    if (pkd->oNodeMom) {
	*pkdNodeMom(pkd,pkdn) = *pkdNodeMom(pkd,p1);
	x = r1[0] - pkdn->r[0];
	y = r1[1] - pkdn->r[1];
	z = r1[2] - pkdn->r[2];
	momShiftMomr(pkdNodeMom(pkd,pkdn),x,y,z);
	mom = *pkdNodeMom(pkd,p2);
	x = r2[0] - pkdn->r[0];
	y = r2[1] - pkdn->r[1];
	z = r2[2] - pkdn->r[2];
	momShiftMomr(&mom,x,y,z);
	momAddMomr(pkdNodeMom(pkd,pkdn),&mom);
	}
    PSDBND_COMBINE(pkdn->bnd,p1->bnd,p2->bnd);
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
#endif


#if 0
void pkdVATreeBuild(PKD pkd,int nBucket,FLOAT diCrit2) {
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

    Create(pkd,VAROOT,diCrit2);
    }
#endif

void psdBuildTree(PKD pkd,int nBucket,FLOAT diCrit2,KDN *pkdn) {
    int iStart;

    assert(pkd->oNodeBnd6);
    assert(pkd->oPsMetric);
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

    psdInitializeParticles(pkd,&pkd->bnd);

#if CUBICCELLS
    BuildPsdTemp(pkd,ROOT,nBucket, 128);
#elif CUBICCELLS2
    BuildPsdTemp(pkd,ROOT,nBucket, 1024);
#else
    BuildPsdTemp(pkd,ROOT,nBucket, 1000000);
#endif

    pkd->nNodesFull = pkd->nNodes;
    iStart = 0;
    ShuffleParticles(pkd,iStart);

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

#if 0
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
#endif

#if 0
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
#endif

#if 0
void pkdDistribRoot(PKD pkd,MOMC *pmom) {
    pkd->momRoot = *pmom;
    }
#endif


#if 0
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
#endif


#if 0
void pkdBoundWalk(PKD pkd,PSDBND *pbnd,uint8_t uRungLo,uint8_t uRungHi,uint32_t *pnActive,uint32_t *pnContained) {
    KDN *kdn;
    PARTICLE *p;
    double d;
    int iNode,pj;    

    *pnActive = 0;
    *pnContained = 0;
    kdn = pkdTreeNode(pkd,iNode = ROOT);
    while (1) {
	d = fabs(pbnd->fCenter[0] - kdn->bnd.fCenter[0]) - pbnd->fMax[0];
	if (d - kdn->bnd.fMax[0] > 0) goto NoIntersect;
	else if (d + kdn->bnd.fMax[0] <= 0) {
	    d = fabs(pbnd->fCenter[1] - kdn->bnd.fCenter[1]) - pbnd->fMax[1];
	    if (d - kdn->bnd.fMax[1] > 0) goto NoIntersect;
	    else if (d + kdn->bnd.fMax[1] <= 0) {
		d = fabs(pbnd->fCenter[2] - kdn->bnd.fCenter[2]) - pbnd->fMax[2];
		if (d - kdn->bnd.fMax[2] > 0) goto NoIntersect;
		else if (d + kdn->bnd.fMax[2] <= 0) goto Contained;
		}
	    else {
		d = fabs(pbnd->fCenter[2] - kdn->bnd.fCenter[2]) - pbnd->fMax[2];
		if (d - kdn->bnd.fMax[2] > 0) goto NoIntersect;
		}
	    }
	else {
	    d = fabs(pbnd->fCenter[1] - kdn->bnd.fCenter[1]) - pbnd->fMax[1];
	    if (d - kdn->bnd.fMax[1] > 0) goto NoIntersect;
	    d = fabs(pbnd->fCenter[2] - kdn->bnd.fCenter[2]) - pbnd->fMax[2];
	    if (d - kdn->bnd.fMax[2] > 0) goto NoIntersect;
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
#endif

#endif /* USE_PSD */

