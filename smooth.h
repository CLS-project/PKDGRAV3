
#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED
#define SMOOTH_H_MODULE_ID "$Id$"

#include "pkd.h"
#include "smoothfcn.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#define NNLIST_INCREMENT	200		/* number of extra neighbor elements added to nnList */


struct hashElement {
    PARTICLE *p;
    struct hashElement *coll;
};


typedef struct smContext {
    PKD pkd;
    PARTICLE pSentinel;
    void (*fcnSmooth)(PARTICLE *,int,NN *,SMF *);
    void (*fcnPost)(void *,PARTICLE *,SMF *);
    int nSmooth;
    int nQueue;
    int bPeriodic;
    PQ *pq;
    /*
    ** Flags to mark local particles which are inactive either because they
    ** are source inactive or because they are already present in the prioq.
    */
    char *bInactive;
    /*
    ** Hash table to indicate whether a remote particle is already present in the 
    ** priority queue.
    */
    int nHash;  /* should be a prime number > nSmooth */
    struct hashElement *pHash;
    struct hashElement *pFreeHash;
    int nnListSize;
    int nnListMax;
    NN *nnList;
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
    } * SMX;


/*
** Assumes that p does not already occur in the hash table!!!
*/
inline void smHashAdd(SMX smx,PARTICLE *p) {
    struct hashElement *t;
    uint32_t i = p%smx->nHash;
    if (!smx->pHash[i].p) {
	smx->pHash[i] = p;
    }
    else {
	t = smx->pFreeHash;
	assert(t != NULL);
	smx->pFreeHash = t->coll;
	t->coll = smx->pHash[i].coll;
	smx->pHash[i].coll = t;
	t->p = p;
    }
}

/*
** Assumes that p is definitely in the hash table!!!
*/
inline void smHashDel(SMX smx,PARTICLE *p) {
    struct hashElement *t,*tt;
    uint32_t i = p%smx->nHash;
    if (!smx->pHash[i].coll) {
	/*
	** It has to be the first element.
	*/
	smx->pHash[i].p = NULL;
    }
    else if (smx->pHash[i].p == p) {
	/*
	** It is the first element, but there are others!
	*/
	t = smx->pHash[i].coll;
	smx->pHash[i].coll = t->coll;
	smx->pHash[i].p = t->p;
	t->coll = smx->pFreeHash;
	smx->pFreeHash = t;
    }
    else {
	tt = &smx->pHash[i];
	while (tt->next->p != p) tt = tt->next;
	t = tt->coll;
	tt->coll = t->coll; /* unlink */
	t->coll = smx->pFreeHash;
	smx->pFreeHash = t;	
    }
}


inline int smHashPresent(SMX smx,PARTICLE *p) {
    struct hashElement *t;
    uint32_t i = p%smx->nHash;

    if (smx->pHash[i].p == p) return 1;
    t = smx->pHash[i].coll;
    while (t) {
	if (t->p == p) return 1;
	else t = t->coll;
    }
    return 0;
}



int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,
		 int bPeriodic,int bSymmetric,int iSmoothType);

void smFinish(SMX,SMF *);
void smSmooth(SMX,SMF *);
void smReSmooth(SMX,SMF *);

void smFof(SMX smx, SMF *smf);
int smGroupMerge(SMF *smf, int bPeriodic);
int smGroupProfiles(SMX smx, SMF *smf,int nTotalGroups);

#endif
