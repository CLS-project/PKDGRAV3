#ifndef GROUP_H
#define GROUP_H

#include "pkd.h"

struct smGroupArray {
    remoteID id;       /* iPid, iIndex */
    int32_t iGid;      /* local group ID */
    union {
	int32_t iNewGid;   /* new local group ID */
	uint32_t iLink;    /* link to remote groups */
	};
    };

int pkdCombineDuplicateGroupIds(PKD pkd, int nGroups, struct smGroupArray *ga,int bIndexIsGID);
void pkdPurgeSmallGroups(PKD pkd,int nMinGroupSize, int bPeriodic, double *dPeriod);
void pkdHopSendStats(PKD pkd);
void pkdHopAssignGID(PKD pkd,uint64_t iStartGID);
int pkdHopCountGID(PKD pkd);

#endif
