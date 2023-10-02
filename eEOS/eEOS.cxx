#include "master.h"

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
void MSR::SetEOSParam() {
    const double dHydFrac = param.dInitialH;
    const double dnHToRho = MHYDR / dHydFrac / units.dGmPerCcUnit;
    param.dEOSFloorDen = param.dEOSFloornH*dnHToRho;
    param.dEOSFlooru = param.dEOSFloorTemp*dTuFacPrimNeutral;
    if (csm->val.bComove)
        param.dEOSFloorMinBaryonOD = param.dEOSFloorMinOD*csm->val.dOmegab;
    else
        param.dEOSFloorMinBaryonOD = 0.;
#ifdef EEOS_POLYTROPE
    param.dEOSPolyFloorExponent = param.dEOSPolyFloorIndex-1.;
    param.dEOSPolyFloorDen =  param.dEOSPolyFloornH*dnHToRho;
    param.dEOSPolyFlooru = param.dEOSPolyFloorTemp*dTuFacPrimNeutral;
    if (csm->val.bComove)
        param.dEOSPolyFloorMinBaryonOD = param.dEOSPolyFloorMinOD*csm->val.dOmegab;
    else
        param.dEOSPolyFloorMinBaryonOD = 0.;
#endif
}


int MSR::ValidateEOSParam() {
    if (!prmSpecified(prm, "dOmegab") && prmSpecified(prm, "dEOSFloorMinOD")) {
        fprintf(stderr,"ERROR: dEOSFloorMinOD is specified but dOmegab is not set\n");
        return 0;
    }
#ifdef EEOS_POLYTROPE
    if (!prmSpecified(prm, "dOmegab") && prmSpecified(prm, "dEOSPolyFloorMinOD")) {
        fprintf(stderr,"ERROR: dEOSPolyFloorMinOD is specified but dOmegab is not set\n");
        return 0;
    }
#endif
    return 1;
}

#endif
