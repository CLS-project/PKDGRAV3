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
    return (uRungLo <= c->uMaxRung && uRungHi >= c->uMinRung);
}
#ifdef __cplusplus
extern "C" {
#endif
/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
    double dTime,int nReps,int bEwald,int nGroup, int iRoot1, int iRoot2,
    int iVARoot, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum,SPHOptions *SPHoptions);

int pkdGravWalkGroups(PKD pkd,double dTime,int nGroup, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum);
#ifdef __cplusplus
}
#endif

#endif
