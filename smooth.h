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
    void *p;
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


int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,
		 int bPeriodic,int bSymmetric,int iSmoothType);

void smFinish(SMX,SMF *);
void smSmooth(SMX,SMF *);
void smReSmooth(SMX,SMF *);

void smFof(SMX smx, SMF *smf);
int smGroupMerge(SMF *smf, int bPeriodic);
int smGroupProfiles(SMX smx, SMF *smf,int nTotalGroups);

#endif
