#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "listcomp.h"
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

struct smExtraArray {
    uint32_t iIndex;
    char bInactive;
    char bDone;
};

typedef struct smContext {
    PKD pkd;
    PARTICLE pSentinel;
    void (*fcnSmooth)(PARTICLE *,int,NN *,SMF *);
    void (*fcnPost)(void *,PARTICLE *,SMF *);
    int nSmooth;
    int nQueue;
    int bPeriodic;
    int bOwnCache;
    FLOAT rLast[3]; /* For the snake */
    PQ *pq;
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
    ** Similarly, this will also use a portion of the pLite.
    */
    struct smGroupArray *ga;
    /*
    ** Flags to mark local particles which are finished in the processing
    ** of the smFastGas routine. They have updated their densities.
    */
    char *bDone;
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
    /*
    ** Context for nearest neighbor lists.
    */
    LCODE lcmp;
    } * SMX;


int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,
		 int bPeriodic,int bSymmetric,int iSmoothType);
int smInitializeRO(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,
		   int bPeriodic,int iSmoothType);
void smFinish(SMX,SMF *);
void smSmoothInitialize(SMX smx);
void smSmoothFinish(SMX smx);
void smSmoothSingle(SMX smx,SMF *smf,PARTICLE *p);
void smSmooth(SMX,SMF *);
void smReSmooth(SMX,SMF *);

void smFastGasPhase1(SMX smx,SMF *smf);
void smFastGasPhase2(SMX smx,SMF *smf);
void pkdFastGasCleanup(PKD pkd);  /* frees up the neighbor lists */

void smFof(SMX smx, SMF *smf);
int smGroupMerge(SMF *smf, int bPeriodic);
int smGroupProfiles(SMX smx, SMF *smf,int nTotalGroups);

void smHopLink(SMX smx,SMF *smf);
int smHopJoin(PKD pkd);
#endif
