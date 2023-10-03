#include "master.h"

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)

eEOSparam::eEOSparam(class pkd_parameters &parameters,struct CALC &calc) {
#ifdef EEOS_POLYTROPE
    dPolyFloorMinOD = calc.dEOSPolyFloorMinBaryonOD;
    dPolyFloorExponent = calc.dEOSPolyFloorExponent;
    dPolyFloorDen = calc.dEOSPolyFloorDen;
    dPolyFlooru = calc.dEOSPolyFlooru;
#endif
#ifdef EEOS_JEANS
    dNJeans = parameters.get_dEOSNJeans();
#endif
    dFlooru = calc.dEOSFlooru;
    dFloorDen = calc.dEOSFloorDen;
    dFloorMinOD = calc.dEOSFloorMinBaryonOD;
}

void MSR::SetEOSParam() {
    const double dHydFrac = parameters.get_dInitialH();
    const double dnHToRho = MHYDR / dHydFrac / units.dGmPerCcUnit;
    calc.dEOSFloorDen = parameters.get_dEOSFloornH()*dnHToRho;
    calc.dEOSFlooru = parameters.get_dEOSFloorTemp()*dTuFacPrimNeutral;
    if (csm->val.bComove)
        calc.dEOSFloorMinBaryonOD = parameters.get_dEOSFloorMinOD()*csm->val.dOmegab;
    else
        calc.dEOSFloorMinBaryonOD = 0.;
#ifdef EEOS_POLYTROPE
    calc.dEOSPolyFloorExponent = parameters.get_dEOSPolyFloorIndex()-1.;
    calc.dEOSPolyFloorDen =  parameters.get_dEOSPolyFloornH()*dnHToRho;
    calc.dEOSPolyFlooru = parameters.get_dEOSPolyFloorTemp()*dTuFacPrimNeutral;
    if (csm->val.bComove)
        calc.dEOSPolyFloorMinBaryonOD = parameters.get_dEOSPolyFloorMinOD()*csm->val.dOmegab;
    else
        calc.dEOSPolyFloorMinBaryonOD = 0.;
#endif
}

int MSR::ValidateEOSParam() {
    if (!parameters.has_dOmegab() && parameters.has_dEOSFloorMinOD()) {
        fprintf(stderr,"ERROR: dEOSFloorMinOD is specified but dOmegab is not set\n");
        return 0;
    }
#ifdef EEOS_POLYTROPE
    if (!parameters.has_dOmegab() && parameters.has_dEOSPolyFloorMinOD()) {
        fprintf(stderr,"ERROR: dEOSPolyFloorMinOD is specified but dOmegab is not set\n");
        return 0;
    }
#endif
    return 1;
}

#endif
