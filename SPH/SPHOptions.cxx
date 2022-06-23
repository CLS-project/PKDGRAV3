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

#include "SPHOptions.h"
#include "parameters.h"
#include <stdio.h>
#include "basetype.h"

SPHOptions initializeSPHOptions(struct parameters param, CSM csm, double dTime) {
    SPHOptions SPHoptions;
    SPHoptions.fKernelTarget = param.fKernelTarget;
    SPHoptions.epsilon = 0.01f;
    SPHoptions.alpha = param.dConstAlpha;
    SPHoptions.beta = param.dConstBeta;
    SPHoptions.EtaCourant = param.dEtaCourant;
    SPHoptions.gamma = param.dConstGamma;
    SPHoptions.TuFac = param.units.dGasConst/(param.dConstGamma - 1)/param.dMeanMolWeight;
    SPHoptions.FastGasFraction = param.dFastGasFraction;
    SPHoptions.VelocityDamper = param.dVelocityDamper;
    SPHoptions.nSmooth = param.nSmooth;
    SPHoptions.ballSizeLimit = param.dBallSizeLimit;
    SPHoptions.fBallFactor = 1.1f;
    SPHoptions.dKpcUnit = param.units.dKpcUnit;
    SPHoptions.dMsolUnit = param.units.dMsolUnit;
    SPHoptions.dMeanMolWeight = param.dMeanMolWeight;
    SPHoptions.nPredictRung = 0;
    SPHoptions.nRungCorrection = 2;
    if (csm->val.bComove) {
        SPHoptions.a = csmTime2Exp(csm,dTime);
        SPHoptions.H = csmTime2Hub(csm,dTime);
    }
    else {
        SPHoptions.a = 1.0f;
        SPHoptions.H = 0.0f;
    }
    SPHoptions.doGravity = 0;
    SPHoptions.doDensity = 0;
    SPHoptions.doSPHForces = 0;
    SPHoptions.doUConversion = 0;
    SPHoptions.doSetDensityFlags = 0;
    SPHoptions.dofBallFactor = 0;
    SPHoptions.useNumDen = 0;
    SPHoptions.useIsentropic = param.bGasIsentropic;
    SPHoptions.useBuiltinIdeal = param.bGasBuiltinIdeal;
    SPHoptions.useDensityFlags = 0;
    SPHoptions.doOnTheFlyPrediction = 0;
    SPHoptions.kernelType = param.iKernelType;
    return SPHoptions;
}

void copySPHOptions(SPHOptions *source, SPHOptions *target) {
    target->fKernelTarget = source->fKernelTarget;
    target->epsilon = source->epsilon;
    target->alpha = source->alpha;
    target->beta = source->beta;
    target->EtaCourant = source->EtaCourant;
    target->gamma = source->gamma;
    target->TuFac = source->TuFac;
    target->FastGasFraction = source->FastGasFraction;
    target->VelocityDamper = source->VelocityDamper;
    target->nSmooth = source->nSmooth;
    target->ballSizeLimit = source->ballSizeLimit;
    target->fBallFactor = source->fBallFactor;
    target->dKpcUnit = source->dKpcUnit;
    target->dMsolUnit = source->dMsolUnit;
    target->dMeanMolWeight = source->dMeanMolWeight;
    target->nPredictRung = source->nPredictRung;
    target->nRungCorrection = source->nRungCorrection;
    target->a = source->a;
    target->H = source->H;
    target->doGravity = 0;
    target->doDensity = 0;
    target->doSPHForces = 0;
    target->doUConversion = 0;
    target->doSetDensityFlags = 0;
    target->dofBallFactor = 0;
    target->useNumDen = source->useNumDen;
    target->useIsentropic = source->useIsentropic;
    target->useBuiltinIdeal = source->useBuiltinIdeal;
    target->useDensityFlags = 0;
    target->doOnTheFlyPrediction = source->doOnTheFlyPrediction;
    target->kernelType = source->kernelType;
}

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
