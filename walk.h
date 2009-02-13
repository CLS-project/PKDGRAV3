#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED
#define WALK_H_MODULE_ID "$Id$"

#include "pkd.h"

static inline int pkdIsCellActive(KDN *c,uint8_t uRungLo,uint8_t uRungHi) {
    return (uRungLo <= c->uMaxRung && uRungHi >= c->uMinRung) && c->bDstActive;
    }

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,
		int bVeryActive,double *pdFlop,double *pdPartSum,double *pdCellSum);
#endif
