#ifndef GROUP_H
#define GROUP_H

#include "pkd.h"

struct smGroupArray {
    remoteID id;       /* iPid, iIndex */
    int32_t iGid;      /* local group ID */
    union {
	int32_t iNewGid;   /* new local group ID */
	uint32_t iLink;    /* link to remote groups */
	uint32_t nTotal;   /* count of total number of particles in the group */
	};
    };

int pkdGroupCombineDuplicateIds(PKD pkd,int nGroups, struct smGroupArray *ga,int bIndexIsGID);
int pkdGroupRelocate(PKD pkd,int nGroups,struct smGroupArray *ga);
int pkdPurgeSmallGroups(PKD pkd,int nGroups, struct smGroupArray *ga,int nMinGroupSize);
void pkdGroupAssignGID(PKD pkd,uint64_t iStartGID);
int pkdGroupCountGID(PKD pkd);

#endif
