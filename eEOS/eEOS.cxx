#include "master.h"

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
void MSR::SetEOSParam() {
#ifdef EEOS_POLYTROPE
    param.dEOSPolyFlooru = param.dEOSPolyFloorTemp * dTuFac * (param.dMeanMolWeight / 1.2285);

    if (csm->val.bComove) {
        assert(csm->val.dOmegab > 0.);
        // If in PKDGRAV3 units, this should always be unity
        const double rhoCrit0 = 3. * csm->val.dHubble0 * csm->val.dHubble0 /
                                (8. * M_PI);
        param.dEOSPolyFloorOD = rhoCrit0 * csm->val.dOmegab *
                                param.dEOSMinOverDensity;
    }
    else {
        param.dEOSPolyFloorOD = 0.0;
    }
#endif

    if (!param.bRestart) {
#ifdef EEOS_POLYTROPE
        param.dEOSPolyFloorIndex -= 1.;
        const double dnHToRho = MHYDR / param.dInitialH / param.units.dGmPerCcUnit;
        param.dEOSPolyFloorDen *= dnHToRho; // Code density
#endif
#ifdef EEOS_JEANS
        param.dEOSNJeans = pow(param.dEOSNJeans, 0.666666);
#endif
    }
}
#endif
