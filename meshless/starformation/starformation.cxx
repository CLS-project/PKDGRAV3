#include "starformation.h"
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#include "eEOS/eEOS.h"
#include "master.h"

/* MSR layer
 */
void MSR::SetStarFormationParam() {
    calc.dSFThresholdu = parameters.get_dSFThresholdT() * dTuFacPrimNeutral;

    if (csm->val.bComove) {
        // As usual, in code units the critical density is 1.0
        calc.dSFThresholdOD = parameters.get_dSFMinOverDensity() * csm->val.dOmegab;
    }
    else {
        calc.dSFThresholdOD = 0.0;
    }

    const double dnHToRho = MHYDR / parameters.get_dInitialH() / units.dGmPerCcUnit;
    calc.dSFThresholdDen = parameters.get_dSFThresholdDen() * dnHToRho; // Code physical density now
    const double Msolpcm2 = 1. / units.dMsolUnit * pow(units.dKpcUnit * 1e3, 2);
    calc.dSFnormalizationKS = parameters.get_dSFnormalizationKS() / units.dMsolUnit * units.dSecUnit / SECONDSPERYEAR * pow(units.dKpcUnit, 2) *
                              pow(Msolpcm2, -parameters.get_dSFindexKS());
}

int MSR::ValidateStarFormationParam() {
    if (!parameters.get_bDoEOS()) {
        print_warning("WARNING: Star formation is active but no eEOS is selected!\n");
    }
    return 1;
}

void MSR::StarForm(double dTime, double dDelta, int iRung) {
    struct inStarForm in;
    struct outStarForm out;

    TimerStart(TIMER_STARFORM);

    const double a = csmTime2Exp(csm, dTime);

    in.dHubble = Comove() ? csmTime2Hub(csm, dTime) : 0.0;
    in.dScaleFactor = a;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.dSFindexKS = parameters.get_dSFindexKS();
    in.dSFnormalizationKS = calc.dSFnormalizationKS;
    in.dConstGamma = parameters.get_dConstGamma();
    in.dSFGasFraction = parameters.get_dSFGasFraction();
    in.dSFThresholdDen = calc.dSFThresholdDen * a * a * a; // Send in comoving units
    in.dSFThresholdOD = calc.dSFThresholdOD;
    in.dSFThresholdu = calc.dSFThresholdu;
    in.dSFEfficiency = parameters.get_dSFEfficiency();
#ifdef HAVE_METALLICITY
    in.bSFThresholdDenSchaye2004 = parameters.get_bSFThresholdDenSchaye2004();
#endif
    in.bCCSNFeedback = parameters.get_bCCSNFeedback();
    in.bSNIaFeedback = parameters.get_bSNIaFeedback();
    in.dSNFBEfficiency = parameters.get_dSNFBEfficiency();
    in.dSNFBMaxEff = parameters.get_dSNFBMaxEff();
    in.dSNFBEffIndex = parameters.get_dSNFBEffIndex();
    in.dSNFBEffnH0 = calc.dSNFBEffnH0;
    in.eEOS = eEOSparam(parameters, calc);
#ifdef STELLAR_EVOLUTION
    in.dSNIaMaxMass = parameters.get_dSNIaMaxMass();
    in.dCCSNMinMass = parameters.get_dCCSNMinMass();
    in.dCCSNMaxMass = parameters.get_dCCSNMaxMass();
#endif

    if (parameters.get_bVDetails()) printf("Star Form (rung: %d) ... ", iRung);

    ActiveRung(iRung, 1);
    mdl->RunService<ServiceStarForm>(sizeof(in), &in, &out, sizeof(out));
    if (parameters.get_bVDetails()) printf("%d Stars formed with mass %g, %d gas deleted\n", out.nFormed, out.dMassFormed, out.nDeleted);
    massFormed += out.dMassFormed;
    starFormed += out.nFormed;

    nGas -= out.nFormed;
    nStar += out.nFormed;

    TimerStop(TIMER_STARFORM);
    auto dsec = TimerGet(TIMER_STARFORM);
    printf("Star Formation Calculated, Wallclock: %f secs\n\n", dsec);
}
