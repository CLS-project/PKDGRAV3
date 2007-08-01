#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"

int pkdGravWalk(PKD pkd,double dTime,int nReps,int bEwald,int bEwaldKick,int bVeryActive,double fEwCut,
				double *pdFlop,double *pdPartSum,double *pdCellSum);

#endif
