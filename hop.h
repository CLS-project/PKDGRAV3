#ifndef HOP_H
#define HOP_H
#include "smooth.h"

int smHopLink(SMX smx,SMF *smf);
int smHopJoin(SMX smx,SMF *smf,double dHopTau,int *nLocal);
int pkdHopFinishUp(PKD pkd, int nMinGroupSize, int bPeriodic, double *dPeriod);
void pkdHopAssignGID(PKD pkd);
void pkdHopTreeBuild(PKD pkd,int nBucket);
int pkdHopUnbind(PKD pkd,double dTime,int nMinGroupSize, int bPeriodic, double *dPeriod);
void pkdHopSendStats(PKD pkd);
int pkdGravWalkHop(PKD pkd,double dTime,int nGroup, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum);

int smNewFof(SMX smx,SMF *smf);

#endif
