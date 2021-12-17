/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GROUP_H
#define GROUP_H

#include "pkd.h"

struct smGroupArray {
    remoteID id;       /* iPid, iIndex */
    int32_t iGid;      /* local group ID */
    union {
        int32_t iNewGid;   /* new local group ID */
        uint32_t iLink;    /* link to remote groups (used during fof) */
        uint32_t nTotal;   /* count of total number of particles in the group (set by pkdGroupCounts) */
    };
    float minPot;
    uint32_t iMinPart;
};

#ifdef __cplusplus
extern "C" {
#endif
int pkdGroupCombineDuplicateIds(PKD pkd,int nGroups, struct smGroupArray *ga,int bIndexIsGID);
int pkdGroupRelocate(PKD pkd,int nGroups,struct smGroupArray *ga);
int pkdPurgeSmallGroups(PKD pkd,int nGroups, struct smGroupArray *ga,int nMinGroupSize);
int pkdGroupCounts(PKD pkd,int nGroups, struct smGroupArray *ga);
#ifdef __cplusplus
}
#endif
#endif
