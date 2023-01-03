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

void SPHpredictOnTheFly(PKD pkd, PARTICLE *p, struct pkdKickParameters *kick, int uRungLo, float *vpred, float *P, float *cs, float *T, SPHOptions *SPHoptions) {
    NEWSPHFIELDS *pNewSph = pkdNewSph(pkd,p);
    float dtPredDrift = getDtPredDrift(kick,p->bMarked,uRungLo,p->uRung);
    const float *ap = pkdAccel(pkd,p);
    const vel_t *v = pkdVel(pkd,p);
    if (SPHoptions->doConsistentPrediction) {
        vpred[0] = pNewSph->vpredx;
        vpred[1] = pNewSph->vpredy;
        vpred[2] = pNewSph->vpredz;
    }
    else {
        vpred[0] = v[0] + dtPredDrift * ap[0];
        vpred[1] = v[1] + dtPredDrift * ap[1];
        vpred[2] = v[2] + dtPredDrift * ap[2];
        if ((SPHoptions->VelocityDamper > 0.0) && p->bMarked) {
            vpred[0] *= exp(kick->dtOpen[p->uRung] * SPHoptions->VelocityDamper);
            vpred[1] *= exp(kick->dtOpen[p->uRung] * SPHoptions->VelocityDamper);
            vpred[2] *= exp(kick->dtOpen[p->uRung] * SPHoptions->VelocityDamper);
        }
    }
    if (SPHoptions->doSPHForces || SPHoptions->doDensityCorrection) {
        if (SPHoptions->doOnTheFlyPrediction) {
            float uPred = 0.0f;
            if (SPHoptions->useIsentropic && !(pkdiMat(pkd,p) == 0 && SPHoptions->useBuiltinIdeal)) {
                if (p->bMarked) {
                    // undo kick
                    uPred = pNewSph->u + dtPredDrift * pNewSph->uDot;
                }
                else {
                    // undo kick
                    uPred = pNewSph->u + kick->dtPredISPHUndoOpen[p->uRung] * pNewSph->uDot;
                    // new opening kick
                    uPred += kick->dtPredISPHOpen[p->uRung] * pNewSph->uDot;
                    // isentropic evolution
                    uPred = SPHEOSIsentropic(pkd,pNewSph->oldRho,uPred,pkdDensity(pkd,p),pkdiMat(pkd,p),SPHoptions);
                    // new closing kick
                    uPred += kick->dtPredISPHClose[p->uRung] * pNewSph->uDot;
                }
            }
            else {
                uPred = pNewSph->u + dtPredDrift * pNewSph->uDot;
            }
            *P = SPHEOSPCTofRhoU(pkd,pkdDensity(pkd,p),uPred,cs,T,pkdiMat(pkd,p),SPHoptions);
        }
        else {
            *P = pNewSph->P;
            *cs = pNewSph->cs;
            if (T) *T = pNewSph->T;
        }
    }
}

void SPHpredictInDensity(PKD pkd, PARTICLE *p, struct pkdKickParameters *kick, int uRungLo, float *P, float *cs, float *T, SPHOptions *SPHoptions) {
    // CAREFUL!! When this is called, p->bMarked does not mean "has been kicked", but it is a fastgas marker
    NEWSPHFIELDS *pNewSph = pkdNewSph(pkd,p);
    if (SPHoptions->doUConversion && SPHoptions->doInterfaceCorrection) {
        *T = pNewSph->u;
        *P = SPHEOSPofRhoT(pkd, pkdDensity(pkd,p), pNewSph->u, pkdiMat(pkd,p), SPHoptions);
    }
    else {
        float dtPredDrift = getDtPredDrift(kick,0,uRungLo,p->uRung);
        float uPred = 0.0f;
        if (SPHoptions->useIsentropic && !(pkdiMat(pkd,p) == 0 && SPHoptions->useBuiltinIdeal)) {
            // undo kick
            uPred = pNewSph->u + kick->dtPredISPHUndoOpen[p->uRung] * pNewSph->uDot;
            // new opening kick
            uPred += kick->dtPredISPHOpen[p->uRung] * pNewSph->uDot;
            // isentropic evolution
            uPred = SPHEOSIsentropic(pkd,pNewSph->oldRho,uPred,pkdDensity(pkd,p),pkdiMat(pkd,p),SPHoptions);
            // new closing kick
            uPred += kick->dtPredISPHClose[p->uRung] * pNewSph->uDot;
        }
        else {
            uPred = pNewSph->u + dtPredDrift * pNewSph->uDot;
        }
        *P = SPHEOSPCTofRhoU(pkd,pkdDensity(pkd,p),uPred,cs,T,pkdiMat(pkd,p),SPHoptions);
        if (SPHoptions->doConsistentPrediction) {
            const vel_t *v = pkdVel(pkd,p);
            const float *ap = pkdAccel(pkd,p);
            pNewSph->vpredx = v[0] + dtPredDrift * ap[0];
            pNewSph->vpredy = v[1] + dtPredDrift * ap[1];
            pNewSph->vpredz = v[2] + dtPredDrift * ap[2];
        }
    }
}
