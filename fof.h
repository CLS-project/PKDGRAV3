#ifndef FOF_INCLUDED
#define FOF_INCLUDED

#include "pkd.h"

void pkdNewFof(PKD pkd,double dTau2,int nMinMembers);
int pkdFofPhases(PKD pkd);
uint64_t pkdFofFinishUp(PKD pkd,int nMinGroupSize);

#endif

