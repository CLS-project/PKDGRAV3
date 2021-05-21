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

#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED

#include "pkd.h"
#include "moments.h"
#include "smooth.h"

static inline double softmassweight(double m1,double h12,double m2,double h22) {
    double tmp = h12*h22;
    if (m1 == 0.0) return(h22);
    if (m2 == 0.0) return(h12);
    if (tmp > 0.0) return((m1+m2)*tmp/(h22*m1+h12*m2));
    else return(0.0);
    }

void pkdGravStartEwald(PKD pkd);
void pkdGravFinishEwald(PKD pkd);

int pkdGravInteract(PKD pkd,
    struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
    KDN *pBucket,LOCR *pLoc,ILP ilp,ILC ilc,
    float dirLsum,float normLsum,int bEwald,double *pdFlop,
    SMX smx,SMF *smf,int iRoot1,int iRoot2,uint64_t SPHoptions);

void pkdParticleWorkDone(workParticle *work);

#ifdef TIMESTEP_CRITICAL
double pkdRho1(double rhopmaxlocal, double summ, double dir, double x, double y, double z, double vx, double vy, double vz, double EccFacMax);
#endif

#endif



