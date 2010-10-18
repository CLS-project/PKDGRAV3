#ifndef PSDTREE_HINCLUDED
#define PSDTREE_HINCLUDED

#include <stdint.h>
#include <string.h>

#include "mdl.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

/*
** These macros implement the index manipulation tricks we use for moving
** around in the top-tree. Note: these do NOT apply to the local-tree!
*/
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define SIBLING(i) 	(i^1)
#define PARENT(i)	(i>>1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

/*
** This is useful for debugging the very-active force calculation.
*/
#define A_VERY_ACTIVE  1


#define MAX_TIMERS		10

/*
** This constant is used to limit the size of a cell.
** Was #define PKD_MAX_CELL_SIZE (1e-2), but in this version of the code it is
** no longer needed!
*/
#define PKD_MAX_CELL_SIZE 1e20

#define PKD_MAX_CLASSES 256
#define MAX_RUNG     63

/*
** Here we define some special reserved nodes. Node-0 is a sentinel or null node, node-1
** is here defined as the ROOT of the local tree (or top tree), node-2 is unused and
** node-3 is the root node for the very active tree.
*/
#define VAROOT          3
#define ROOT		1
#define NRESERVED_NODES MAX_RUNG+1

#define PSDBND_COMBINE(b,b1,b2)\
{\
	int BND_COMBINE_j;\
	for (BND_COMBINE_j=0;BND_COMBINE_j<6;++BND_COMBINE_j) {\
		FLOAT BND_COMBINE_t1,BND_COMBINE_t2,BND_COMBINE_max,BND_COMBINE_min;\
		BND_COMBINE_t1 = (b1).fCenter[BND_COMBINE_j] + (b1).fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2).fCenter[BND_COMBINE_j] + (b2).fMax[BND_COMBINE_j];\
		BND_COMBINE_max = (BND_COMBINE_t1 > BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		BND_COMBINE_t1 = (b1).fCenter[BND_COMBINE_j] - (b1).fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2).fCenter[BND_COMBINE_j] - (b2).fMax[BND_COMBINE_j];\
		BND_COMBINE_min = (BND_COMBINE_t1 < BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		(b).fCenter[BND_COMBINE_j] = 0.5*(BND_COMBINE_max + BND_COMBINE_min);\
		(b).fMax[BND_COMBINE_j] = 0.5*(BND_COMBINE_max - BND_COMBINE_min);\
		}\
	}

typedef struct {
    double fCenter[6];
    double fMax[6];
    double size;
    } PSDBND;

#define PSDMINDIST(bnd,pos,min2) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<6;++BND_j) {\
	BND_dMin = fabs((bnd).fCenter[BND_j] - (pos)[BND_j]) - (bnd).fMax[BND_j]; \
	if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin;			\
	}\
    }

static inline int IN_PSDBND(const FLOAT *R,const PSDBND *b) {
    int i;
    for( i=0; i<6; i++ )
	if ( R[i]<b->fCenter[i]-b->fMax[i] || R[i]>=b->fCenter[i]+b->fMax[i] )
	    return 0;
    return 1;
    }


/*
** General partition macro
** LT,LE: Compare less-than/less-than or equal
** ii,dj: Increment i and decrement j
** SWAP: Swap the i'th and j'th element
** LOWER,UPPER: comparison predicates
** e.g.,
** PARTICLE *pi = pkdParticle(pkd,i);
** PARTICLE *pj = pkdParticle(pkd,j);
**    PARTITION(pi<pj,pi<=pj,
**              pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
**              pkdSwapParticle(pkd,pi,pj),
**	        pi->r[d] >= fSplit,pj->r[d] < fSplit);
** When finished, the 'I' variable points to the first element of
** the upper partition (or one past the end).
** NOTE: Because this function supports tiled data structures,
**       the LT, followed by "if (LE)" needs to remain this way.
*/
#define PARTITION(LT,LE,INCI,DECJ,SWAP,LOWER,UPPER)	\
    {							\
    while ((LT) && (LOWER)) { INCI; }			\
    if ((LE) && (LOWER)) { INCI; }			\
    else {						\
	while ((LT) && (UPPER)) { DECJ; }		\
	while (LT) {					\
		    { SWAP; }				\
		    do { DECJ; } while (UPPER);	\
		    do { INCI; } while (LOWER);		\
	    }						\
	}						\
    }

typedef struct psdNode {
    double r[6];
    PSDBND bnd;
    int iLower;
    int iParent;
    int pLower;		/* also serves as thread id for the LTT */
    int pUpper;		/* pUpper < 0 indicates no particles in tree! */
#ifdef LOCAL_EXPANSION
    FLOAT fOpen;
#else
    FLOAT fOpen2;
#endif
    FLOAT fSoft2;
    uint32_t nActive; /* local active count used for walk2 */
    uint8_t uMinRung;
    uint8_t uMaxRung;
    uint8_t bSrcActive;
    uint8_t bDstActive;
    } PSDNODE;

#ifdef NEW_TREE
#define MAX_NBUCKET 5

typedef struct psdNew {
    float s;  /* scale of the cell, if -ve then it indicates a bucket! */
    union celltype {
	struct cell {
	    double dSplit;
	    uint32_t iLower;
	    uint32_t iUpper;
	    uint16_t idLower;
	    uint16_t idUpper;
	    } c;
	struct bucket {
	    uint32_t iPart[MAX_NBUCKET];
	    } b;
	};
    double r[6];
    } PSDNEW;
#endif

static inline int max_bnd(const double *fMax)
{
    int i;
    int d = 0;
    double max = fMax[0];
    for (i=1; i < 6; i++)
        if (fMax[i] > max) { max = fMax[i]; d = i; }
    return d;
}

static inline int max_bnd_range(const double *fMax, const int m, const int M)
{
    int i;
    int d = m;
    double max = fMax[m];
    for (i=m+1; i < M; i++)
        if (fMax[i] > max) { max = fMax[i]; d = i; }
    if (max == 0) 
        return -1;
    return d;
}

static inline int min_bnd(const double *fMax)
{
    int i;
    int d = 0;
    double min = fMax[0];
    for (i=1; i < 6; i++)
        if (fMax[i] < min) { min = fMax[i]; d = i; }
    return d;
}

static inline int max_side(const double *fMax)
{
    return 2 * max_bnd(fMax);
}

static inline int min_side(const double *fMax)
{
    return 2 * min_bnd(fMax);
}

void psdBuildTree(PKD pkd,int nBucket,FLOAT diCrit2,KDN *pkdn);
int psdTreeDepth(PKD pkd);

#endif
