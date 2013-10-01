#ifndef EWALD_HINCLUDED
#define EWALD_HINCLUDED

#include "pkd.h"

int pkdParticleEwald(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    PARTICLE *p, float *pa, float *pPot);
int pkdParticleEwaldSIMD(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    PARTICLE *p, float *pa, float *pPot);

void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut);

#endif
