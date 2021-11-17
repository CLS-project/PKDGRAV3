#include "eEOS/eEOS.h"

void msrSetEOSParam(MSR msr){
    const double dHydFrac = msr->param.dInitialH;
    const double dnHToRho = MHYDR / dHydFrac / msr->param.dGmPerCcUnit;
#ifdef EEOS_POLYTROPE
    msr->param.dEOSPolyFloorIndex -= 1.;
    msr->param.dEOSPolyFloorDen *=  dnHToRho; // Code density
    msr->param.dEOSPolyFlooru *= msr->param.dTuFac; // Code internal energy 
                                                    // per unit mass
#endif
#ifdef EEOS_JEANS
    msr->param.dEOSNJeans = pow(msr->param.dEOSNJeans, 0.666666);
#endif
}
