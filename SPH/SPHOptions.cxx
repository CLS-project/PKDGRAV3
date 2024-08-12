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
#include "pkd_parameters.h"
#include <stdio.h>
#include "basetype.h"

SPHOptions initializeSPHOptions(pkd_parameters &parameters, CSM csm, double dTime) {
    UNITS units(parameters,csm->val.h);
    SPHOptions SPHoptions;
    SPHoptions.fKernelTarget = parameters.get_nSmooth();
    SPHoptions.epsilon = 0.01f;
    SPHoptions.alpha = parameters.get_dConstAlpha();
    SPHoptions.beta = parameters.get_dConstBeta();
    SPHoptions.EtaCourant = parameters.get_dEtaCourant();
    SPHoptions.EtauDot = parameters.get_dEtauDot();
    SPHoptions.EtaSdot = parameters.get_dEtaSdot();
    SPHoptions.timeStepSmin = parameters.get_dTimeStepSmin();
    SPHoptions.gamma = parameters.get_dConstGamma();
    SPHoptions.TuFac = units.dGasConst/(SPHoptions.gamma - 1)/parameters.get_dMeanMolWeight();
    SPHoptions.FastGasFraction = parameters.get_dFastGasFraction();
    SPHoptions.VelocityDamper = 0.0f;
    auto dDelta = parameters.get_dDelta();
    auto dVelocityDamper = parameters.get_dVelocityDamper();
    auto dVelocityDamperEnd = parameters.get_dVelocityDamperEnd();
    auto dVelocityDamperEndTime = parameters.get_dVelocityDamperEndTime();
    if (dDelta > 0.0 && dVelocityDamper > 0.0) {
        if ((dVelocityDamperEnd > 0.0) && (dVelocityDamperEndTime > 0.0)) {
            if (dTime <= dVelocityDamperEndTime) {
                SPHoptions.VelocityDamper = 2.0 / dDelta * pow(10.0, log10(dVelocityDamperEnd/dVelocityDamper)/dVelocityDamperEndTime*dTime + log10(dVelocityDamper));
            }
        }
        else {
            SPHoptions.VelocityDamper = 2.0 / dDelta * dVelocityDamper;
        }
        printf("Velocity Damper active, VelocityDamper = %g\n", SPHoptions.VelocityDamper * dDelta * 0.5);
    }
    SPHoptions.nSmooth = parameters.get_nSmooth();
    SPHoptions.ballSizeLimit = parameters.get_dBallSizeLimit();
    SPHoptions.fBallFactor = 1.1f;
    SPHoptions.dKpcUnit = units.dKpcUnit;
    SPHoptions.dMsolUnit = units.dMsolUnit;
    SPHoptions.dMeanMolWeight = parameters.get_dMeanMolWeight();
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
    SPHoptions.doDensityCorrection = 0;
    SPHoptions.doSPHForces = 0;
    SPHoptions.doUConversion = 0;
    SPHoptions.doSetDensityFlags = 0;
    SPHoptions.dofBallFactor = 1;
    SPHoptions.useIsentropic = parameters.get_bGasIsentropic();
    SPHoptions.useBuiltinIdeal = parameters.get_bGasBuiltinIdeal();
    SPHoptions.useDensityFlags = 0;
    SPHoptions.doOnTheFlyPrediction = parameters.get_bGasOnTheFlyPrediction();
    SPHoptions.doInterfaceCorrection = parameters.get_bGasInterfaceCorrection();
    SPHoptions.doSetNNflags = 0;
    SPHoptions.useNNflags = 0;
    SPHoptions.doConsistentPrediction = parameters.get_bGasConsistentPrediction();
    SPHoptions.kernelType = parameters.get_iKernelType();
    SPHoptions.doCentrifugal = parameters.get_bCentrifugal();
    SPHoptions.CentrifugalT0 = parameters.get_dCentrifT0();
    SPHoptions.CentrifugalT1 = parameters.get_dCentrifT1();
    SPHoptions.CentrifugalOmega0 = parameters.get_dCentrifOmega0();
    SPHoptions.doExtensiveILPTest = parameters.get_bGasDoExtensiveILPTest();
    SPHoptions.doShearStrengthModel = parameters.get_bShearStrengthModel();
    SPHoptions.evolveDensity = parameters.get_bGasEvolveDensity();
    return SPHoptions;
}

void copySPHOptions(SPHOptions *source, SPHOptions *target) {
    target->fKernelTarget = source->fKernelTarget;
    target->epsilon = source->epsilon;
    target->alpha = source->alpha;
    target->beta = source->beta;
    target->EtaCourant = source->EtaCourant;
    target->EtauDot = source->EtauDot;
    target->EtaSdot = source->EtaSdot;
    target->timeStepSmin = source->timeStepSmin;
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
    target->doDensityCorrection = 0;
    target->doSPHForces = 0;
    target->doUConversion = 0;
    target->doSetDensityFlags = 0;
    target->dofBallFactor = 1;
    target->useIsentropic = source->useIsentropic;
    target->useBuiltinIdeal = source->useBuiltinIdeal;
    target->useDensityFlags = 0;
    target->doOnTheFlyPrediction = source->doOnTheFlyPrediction;
    target->doInterfaceCorrection = source->doInterfaceCorrection;
    target->doSetNNflags = source->doSetNNflags;
    target->useNNflags = source->useNNflags;
    target->doConsistentPrediction = source->doConsistentPrediction;
    target->kernelType = source->kernelType;
    target->doCentrifugal = source->doCentrifugal;
    target->CentrifugalT0 = source->CentrifugalT0;
    target->CentrifugalT1 = source->CentrifugalT1;
    target->CentrifugalOmega0 = source->CentrifugalOmega0;
    target->doExtensiveILPTest = source->doExtensiveILPTest;
    target->doShearStrengthModel = source->doShearStrengthModel;
    target->evolveDensity = source->evolveDensity;
}

void copySPHOptionsGPU(SPHOptions *source, SPHOptionsGPU *target) {
    target->kernelType = source->kernelType;
    target->doInterfaceCorrection = source->doInterfaceCorrection;
    target->useIsentropic = source->useIsentropic;
    target->epsilon = source->epsilon;
    target->alpha = source->alpha;
    target->beta = source->beta;
    target->EtaCourant = source->EtaCourant;
    target->a = source->a;
    target->H = source->H;
}

float calculateInterfaceCorrectionPrefactor(float nSmooth,int kernelType) {
    float alpha = 0.0f;
    switch (kernelType) {
    case 0: {
        if (nSmooth < 10.0f) nSmooth = 10.0f; // Below 10, the function decreases
        alpha = -7.44731686e+02f / (nSmooth * nSmooth) + 1.42956727e+02f / nSmooth + 4.46685213e+00f;
        break;
    }
    case 1: {
        if (nSmooth < 3.0f) nSmooth = 3.0f; // Below 3, the function decreases
        alpha = -2.27180786e+02f / (nSmooth * nSmooth) + 1.55942641e+02f / nSmooth + 4.53619157e+00f;
        break;
    }
    case 2: {
        alpha = 1.81941313e+03f / (nSmooth * nSmooth) + 1.73546428e+02f / nSmooth + 4.71450606e+00f;
        break;
    }
    case 3: {
        alpha = 6.00324255e+03f / (nSmooth * nSmooth) + 1.71766754e+02f / nSmooth + 4.92489934e+00f;
        break;
    }
    default: assert(0);
    }
    return alpha;
}
