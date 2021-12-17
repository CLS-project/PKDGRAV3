
#include "smooth/smooth.h"
#include "pkd.h"

#ifdef __cplusplus
extern "C" {
#endif
int pkdPlaceBHSeed(PKD pkd, double dTime, double dScaleFactor,
                   uint8_t uRungMax, double dDenMin, double dBHMhaloMin,
                   double dTau, double dInitialH, double dBHSeedMass);
#ifdef __cplusplus
}
#endif
