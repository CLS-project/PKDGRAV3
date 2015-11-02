#ifndef EWALD_HINCLUDED
#define EWALD_HINCLUDED

#include "pkd.h"

double pkdParticleEwald(PKD pkd, double *r, float *pa, float *pPot);
void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut);

#endif
