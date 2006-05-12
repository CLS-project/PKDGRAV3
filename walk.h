#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"

int pkdGravWalk(PKD pkd,int nReps,int bEwald,double fEwCut,
				double *pdFlop,double *pdPartSum,double *pdCellSum);
void pkdParticleWalk(PKD pkd,int iParticle,int nReps, ILP **pilp, int *pnPart, int *pnMaxPart, ILC **pilc, int *pnCell, int *pnMaxCell, ILPB **pilpb, int *pnPartBucket, int *pnMaxPartBucket);

#endif
