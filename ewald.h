#ifndef EWALD_HINCLUDED
#define EWALD_HINCLUDED

#include "pkd.h"

#ifdef __cplusplus
extern "C" {
#endif
double pkdParticleEwald(PKD pkd, double *r, float *pa, float *pPot,double *pdFlopSingle, double *pdFlopDouble);
void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut);
#ifdef __cplusplus
    }
#endif

#endif
