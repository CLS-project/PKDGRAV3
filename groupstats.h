#ifndef GROUPSTATS_INCLUDED
#define GROUPSTATS_INCLUDED
#include "pkd.h"

void pkdHopSendStats(PKD pkd);
void pkdCalculateGroupStats(PKD pkd, int bPeriodic,double dPeriod[3],double rEnvironment[2]);

#endif
