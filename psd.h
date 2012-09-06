#ifndef PSD_INCLUDED
#define PSD_INCLUDED

#include "pkd.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#include "listcomp.h"

#define NNLIST_INCREMENT    200	/* number of extra neighbor elements added to nnList */

typedef struct pq6Node {
    struct pq6Node *pqLoser;
    struct pq6Node *pqFromInt;
    struct pq6Node *pqFromExt;
    struct pq6Node *pqWinner;    /* Only used when building initial tree */
    PARTICLE *pPart;
    FLOAT fDist2;
    FLOAT dr[6];
    int iIndex;
    int iPid;
    } PQ6;


typedef PQ6 NN6;

struct pkdContext;

typedef struct psfParameters {
    int bComove;
    double dTime;
    double H;
    double a;
    double dDeltaT;
    double dTau2;
    double dVTau2;
    int bTauAbs;
    int nMinMembers;
    int nBins;
    int iCenterType;
    int bLogBins;
    FLOAT binFactor;
    FLOAT fMinRadius;
    struct pkdContext * pkd; /* useful for diagnostics, etc. */
    } PSF;

#if 0
struct hashElement {
    void *p;
    struct hashElement *coll;
};

struct smExtraArray {
    uint32_t iIndex;
    char bInactive;
    char bDone;
};
#endif

struct bridge
{
    int32_t local_gid;
    int32_t remote_gid;
    int iPid;
    int64_t iIndex;
    int done;
};

typedef struct psContext {
    struct pkdContext * pkd;
    PARTICLE pSentinel;
    int nSmooth;
    int bPeriodic;
    PQ6 *pq;
    PSMETRIC *psm;
    struct bridge *bridges;
    int nBridges;
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
    int nnListSize;
    int nnListMax;
    NN6 *nnList;
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

    } * PSX;

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

/******************************************************************************/

/* Standard M_4 Kernel */
#define BALL2(a) ((a)->fBall*(a)->fBall)
#define KERNEL(ak,ar2) { \
	ak = 2.0 - sqrt(ar2); \
	if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
	else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
	else ak = 0.0;\
	}


int psdInitialize(PSX smx, struct pkdContext * pkd, PSF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType, int initCache);
void psdFinish(PSX smx, PSF *smf);
void psdSmooth(PSX smx, PSF *smf);
void psdSmoothLink(PSX smx, PSF *smf);
int psdJoinBridges(PSX psx, PSF *smf);
int psdCountLocalGroups(PSX psx);
void psdMergeNoisyGroups(PSX psx);
void psdSetGlobalId(PSX psx);

#endif

