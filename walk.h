#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bEwald,
		int bVeryActive,double *pdFlop,double *pdPartSum,double *pdCellSum);
#endif
