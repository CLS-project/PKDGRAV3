#ifndef PSD_INCLUDED
#define PSD_INCLUDED

#include "pkd.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#include "listcomp.h"
#include "knn6d.h"

#define NNLIST_INCREMENT    200	/* number of extra neighbor elements added to nnList */

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
    int nSmooth;
    int bPeriodic;
    struct bridge *bridges;
    int nBridges;
    struct knn6dContext *knn;
    PSMETRIC *psm;
    } * PSX;

/******************************************************************************/


/* Standard M_4 Kernel */
#define BALL2(a) ((a)->fBall*(a)->fBall)
#define KERNEL(ak,ar2) { \
	ak = 2.0 - sqrt(ar2); \
	if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
	else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
	else ak = 0.0;\
	}


int psdInitialize(PKD pkd, PSX psx, int nSmooth,int bPeriodic);
void psdFinish(PKD pkd, PSX smx);
void psdSmooth(PKD pkd, PSX smx);
void psdSmoothLink(PKD pkd, PSX smx);
int psdJoinBridges(PKD pkd, PSX psx);
int psdCountLocalGroups(PKD pkd);
void psdMergeNoisyGroups(PKD pkd, PSX psx);
void psdAssignGlobalIds(PKD pkd, int offs, int count);
void psdUpdateGroupProperties(PKD pkd);
void psdSetGlobalId(PKD pkd);
int psdJoinGroupBridges(PKD pkd, PSX psx);
#endif

