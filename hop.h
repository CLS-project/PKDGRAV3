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

#ifndef HOP_H
#define HOP_H
#include "smooth/smooth.h"
#include "group.h"
#ifdef __cplusplus
extern "C" {
#endif
int smHopLink(SMX smx,SMF *smf);
int smHopJoin(SMX smx,SMF *smf,double dHopTau,int *nLocal);
int pkdHopFinishUp(PKD pkd, int nMinGroupSize, int bPeriodic, double *dPeriod);
void pkdHopTreeBuild(PKD pkd,int nBucket,int nGroup);
int pkdHopUnbind(PKD pkd,double dTime,int nMinGroupSize, int bPeriodic, double *dPeriod);
int pkdGravWalkHop(PKD pkd,double dTime,int nGroup, double dThetaMin,double *pdFlop
,double *pdPartSum,double *pdCellSum);
#ifdef __cplusplus
}
#endif

#endif
