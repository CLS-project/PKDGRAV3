#ifndef KNN6D_INCLUDED
#define KNN6D_INCLUDED

#include "pkd.h"
#include "psd.h"

typedef struct pq6Node {
    struct pq6Node *pqLoser;
    struct pq6Node *pqFromInt;
    struct pq6Node *pqFromExt;
    struct pq6Node *pqWinner;    /* Only used when building initial tree */
    PARTICLE *pPart;
    FLOAT dr[6];
    FLOAT fDist2;
    int iIndex;
    int iPid;
    } PQ6;

/*
** Used to keep track of the state of things from a previous call to knn()
*/
struct knn_data
{
    FLOAT rscale[3], vscale[3];
    FLOAT rLast[3], vLast[3];
};


typedef struct knn6dContext
{
    int nnListSize;
    int nnListMax;
    int pqSize;
    PQ6 *pq;
    PQ6 *pqTEMP;
    PSMETRIC *psm;
    struct knn_data data;
    FLOAT fBall2;
    int bPeriodic;

    PARTICLE pSentinel;
    /*
    ** Flags to mark local particles which are inactive either because they
    ** are source inactive or because they are already present in the prioq.
    ** In this extra array is also space for a queue of particles, needed 
    ** for the fast gas routines or for friends-of-friends.
    ** This will point to the pLite array, so it will be destroyed after a tree
    ** build or domain decomposition.
    */
    struct smExtraArray *ea;
    /*
    ** Hash table to indicate whether a remote particle is already present in the 
    ** priority queue.
    */
    int nHash;  /* should be a prime number > nSmooth */
    struct hashElement *pHash;
    struct hashElement *pFreeHash;
    /*
     ** Two stacks for the search algorithm.
     */
    int *S;
    FLOAT *Smin;
    /*
     ** Also need the two stacks for the search
     ** within the top tree.
     */
    int *ST;
    FLOAT *SminT;
} * KNN6D;

#define PQ6_INIT(pq,n)\
{\
    int j;\
    if ((n) == 1) {\
	(pq)[0].pqFromInt = NULL;\
	(pq)[0].pqFromExt = NULL;\
	}\
    else {\
	for (j=0;j<(n);++j) {\
	    if (j < 2) (pq)[j].pqFromInt = NULL;\
	    else (pq)[j].pqFromInt = &(pq)[j>>1];\
	    (pq)[j].pqFromExt = &(pq)[(j+(n))>>1];\
	    }\
	}\
    }


#define PQ6_BUILD(pq,n,q)\
{\
    int i,j;\
    PQ6 *t,*lt;\
    for (j=(n)-1;j>0;--j) {\
	i = (j<<1);\
	if (i < (n)) t = (pq)[i].pqWinner;\
	else t = &(pq)[i-(n)];\
	++i;\
	if (i < (n)) lt = (pq)[i].pqWinner;\
	else lt = &(pq)[i-(n)];\
	if (t->fDist2 < lt->fDist2) {\
	    (pq)[j].pqLoser = t;\
	    (pq)[j].pqWinner = lt;\
	    }\
	else {\
	    (pq)[j].pqLoser = lt;\
	    (pq)[j].pqWinner = t;\
	    }\
	}\
    if ((n) == 1) (q) = (pq);\
    else (q) = (pq)[1].pqWinner;\
    }


#define PQ6_REPLACE(q)\
{\
    PQ6 *t,*lt;\
    t = (q)->pqFromExt;\
    while (t) {\
	if (t->pqLoser->fDist2 > (q)->fDist2) {\
	    lt = t->pqLoser;\
	    t->pqLoser = (q);\
	    (q) = lt;\
	    }\
	t = t->pqFromInt;\
	}\
    }

#if 0

#define PSMINDIST(bnd,r,v,rscale,vscale, min_out) MINDIST(bnd[0],r,min_out)

#else

#define _PSMINDIST(bnd,pos,min2, scale) do {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<3;++BND_j) {\
	BND_dMin = (fabs((bnd)->fCenter[BND_j] - (pos)[BND_j]) - (bnd)->fMax[BND_j]) * scale[BND_j]; \
	BND_dMin *= (BND_dMin > 0); \
	(min2) += BND_dMin*BND_dMin; \
	}\
    } while (0)

#define PSMINDIST(bnd,r,v,rscale,vscale, min_out) do { \
    FLOAT rmin, vmin=0; \
    _PSMINDIST(bnd[0], r, rmin, rscale); \
    _PSMINDIST(bnd[1], v, vmin, vscale); \
    min_out = rmin + vmin; \
} while(0)
#endif

void knn6d(PKD pkd, KNN6D knn, int pid, float *fBall, int first_time);
void knn6dGather(PKD pkd, KNN6D knn, float fBall, int pid, int first_time);

int knn6dInitialize(PKD pkd, KNN6D knn, int nSmooth,int bPeriodic);
void knn6dFree(PKD pkd, KNN6D  smx);
void knn6dFinish(PKD pkd, KNN6D  psx);

#endif
