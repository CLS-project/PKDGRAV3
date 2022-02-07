#include "master.h"

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
void MSR::SetEOSParam() {
    const double dHydFrac = param.dInitialH;
    const double dnHToRho = MHYDR / dHydFrac / param.units.dGmPerCcUnit;
#ifdef EEOS_POLYTROPE
    param.dEOSPolyFloorIndex -= 1.;
    param.dEOSPolyFloorDen *=  dnHToRho; // Code density
    param.dEOSPolyFlooru *= dTuFac; // Code internal energy
    // per unit mass
#endif
#ifdef EEOS_JEANS
    param.dEOSNJeans = pow(param.dEOSNJeans, 0.666666);
#endif
}
#endif
