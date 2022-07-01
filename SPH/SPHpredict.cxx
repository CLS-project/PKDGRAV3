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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif

#include "SPHpredict.h"

float getDtPredDrift(struct pkdKickParameters *kick, int bMarked, int uRungLo, int uRung) {
    if (uRung < uRungLo) {
        return kick->dtPredDrift[uRung];
    }
    else {
        if (bMarked) {
            return - kick->dtOpen[uRung];
        }
        else {
            return kick->dtClose[uRung];
        }
    }
}

void SPHpredictOnTheFly(PKD pkd, PARTICLE *p, struct pkdKickParameters *kick, int uRungLo, float *vpred, float *P, float *cs, SPHOptions *SPHoptions) {
    NEWSPHFIELDS *pNewSph = pkdNewSph(pkd,p);
    float dtPredDrift = getDtPredDrift(kick,p->bMarked,uRungLo,p->uRung);
    const float *ap = pkdAccel(pkd,p);
    const vel_t *v = pkdVel(pkd,p);
    vpred[0] = v[0] + dtPredDrift * ap[0];
    vpred[1] = v[1] + dtPredDrift * ap[1];
    vpred[2] = v[2] + dtPredDrift * ap[2];
    if ((SPHoptions->VelocityDamper > 0.0) & p->bMarked) {
        vpred[0] /= 1.0 - kick->dtClose[p->uRung] * SPHoptions->VelocityDamper;
        vpred[1] /= 1.0 - kick->dtClose[p->uRung] * SPHoptions->VelocityDamper;
        vpred[2] /= 1.0 - kick->dtClose[p->uRung] * SPHoptions->VelocityDamper;
    }
    if (SPHoptions->doSPHForces) {
        if (SPHoptions->doOnTheFlyPrediction) {
            *P = SPHEOSPCofRhoU(pkd,pkdDensity(pkd,p),pNewSph->u + dtPredDrift * pNewSph->uDot,cs,pkdiMat(pkd,p),SPHoptions);
        }
        else {
            *P = pNewSph->P;
            *cs = pNewSph->cs;
        }
    }
}

void SPHpredictInDensity(PKD pkd, PARTICLE *p, struct pkdKickParameters *kick, int uRungLo, float *P, float *cs, SPHOptions *SPHoptions) {
    NEWSPHFIELDS *pNewSph = pkdNewSph(pkd,p);
    float dtPredDrift = getDtPredDrift(kick,0,uRungLo,p->uRung);
    *P = SPHEOSPCofRhoU(pkd,pkdDensity(pkd,p),pNewSph->u + dtPredDrift * pNewSph->uDot,cs,pkdiMat(pkd,p),SPHoptions);
}
