#include "master.h"
#include "pkd.h"
#include "pst.h"
#include <grackle.h>

void msrGrackleInit(MSR msr, int bComove, double dScaleFactor);
void pkdGrackleInit(PKD pkd, int bComove, double dScaleFactor);

void pkdGrackleUpdate(PKD pkd, double dScaleFactor);

void pkdGrackleCooling(PKD pkd, PARTICLE* p, double pDelta);

