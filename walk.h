#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"

/*
** This is really means that the cell MIGHT have an active particle, but it 
** is not for certain, since the contained DstActive particles may not lie in 
** the rung range. To be certain that there are actually active particles
** one has to look at all particles of this cell (recursively walking the 
** subcells).
*/
static inline int pkdIsCellActive(KDN *c,uint8_t uRungLo,uint8_t uRungHi) {
    return (uRungLo <= c->uMaxRung && uRungHi >= c->uMinRung) && c->bDstActive;
    }

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,int bKickClose,int bKickOpen,
    vel_t *dtClose,vel_t *dtOpen,double *dtLCDrift,double *dtLCKick,double dLookbackFac,double dLookbackFacLCP,
    double dAccFac,double dTime,int nReps,int bEwald,int nGroup, int iRoot1, int iRoot2,
    int iVARoot, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum);

int pkdGravWalkGroups(PKD pkd,double dTime,int nGroup, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum);

#endif
