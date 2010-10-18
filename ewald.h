#ifndef EWALD_HINCLUDED
#define EWALD_HINCLUDED

#include "pkd.h"

int pkdParticleEwald(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,PARTICLE *p);
void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut);

#endif
