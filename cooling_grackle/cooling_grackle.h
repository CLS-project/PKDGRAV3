#ifdef __cplusplus
extern "C" {
#endif
#include "pkd.h"
#include "pst.h"
#include <grackle.h>

void pkdGrackleInit(PKD pkd, int bComove, double dScaleFactor, char *achCoolingTable, UNITS units);
void pkdGrackleUpdate(PKD pkd, double dScaleFactor, char *achCoolingTable, UNITS units);
void pkdGrackleCooling(PKD pkd, PARTICLE* p, double pDelta, double dTuFac);

#ifdef __cplusplus
}
#endif
