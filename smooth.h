#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "pkd.h"
#include "smoothfcn.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#define NNLIST_INCREMENT	100		/* number of extra neighbor elements added to nnList */

#define PQ_LOAD_FACTOR 0.50


typedef struct pqNode {
	struct pqNode *pqLoser;
	struct pqNode *pqFromInt;
	struct pqNode *pqFromExt;
	struct pqNode *pqWinner;	/* Only used when building initial tree */
	PARTICLE *pPart;
	FLOAT fDist2;
	FLOAT dx;
	FLOAT dy;
	FLOAT dz;
	int bRemote;
	} PQ;


typedef struct smContext {
	PKD pkd;
	int bLowhFix;
	double dfBall2OverSoft2;
	void (*fcnSmooth)(PARTICLE *,int,NN *,SMF *);
	void (*fcnPost)(PARTICLE *,SMF *);
	int nSmooth;
	int nQueue;
	int bGasOnly;
	int bPeriodic;
	PQ *pq;
	int nnListSize;
	int nnListMax;
	NN *nnList;
	int *nnbRemote;
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


#define PQ_INIT(pq,n)\
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


#define PQ_BUILD(pq,n,q)\
{\
    int i,j;\
	PQ *t,*lt;\
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


#define PQ_REPLACE(q)\
{\
	PQ *t,*lt;\
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


int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bGasOnly,
		 int bPeriodic,int bSymmetric,int iSmoothType,
		 double dfBall2OverSoft2);

void smFinish(SMX,SMF *);
void smSmooth(SMX,SMF *);
void smReSmooth(SMX,SMF *);

void smFof(SMX smx, int nFOFsDone, SMF *smf);
int smGroupMerge(SMF *smf,int bPeriodic);
int smGroupProfiles(SMX smx, SMF *smf,int bPeriodic,int nTotalGroups,int bLogBins,int nFOFsDone);
FLOAT corrPos(FLOAT com, FLOAT r,FLOAT l);
int CmpParticleGroupIds(const void *v1,const void *v2);
int CmpRMs(const void *v1,const void *v2);
int CmpProtoGroups(const void *v1,const void *v2);
FLOAT phase_dist(PARTICLE p1,PARTICLE p2,double h);

#endif
