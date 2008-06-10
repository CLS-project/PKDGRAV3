#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

typedef struct smfParameters {
    double H;
    double a;
    double dDeltaT;
#ifdef SYMBA
    double dSunMass;
#endif
    double dTau2;
    double dVTau2;
    int bTauAbs;
    int nMinMembers;
    int nMinProfile;
    int nBins;
    int bUsePotmin;
    FLOAT fContrast;
    FLOAT Delta;
    FLOAT binFactor;
    PKD pkd; /* useful for diagnostics, etc. */
    } SMF;


typedef struct nNeighbor {
    PARTICLE *pPart;
    FLOAT fDist2;
    FLOAT dx;
    FLOAT dy;
    FLOAT dz;
    int iPid;
    int iIndex;
    } NN;

#define SMX_NULL                            0
void NullSmooth(PARTICLE *,int,NN *,SMF *);

#define SMX_DENSITY				1
void initDensity(void *,void *);
void combDensity(void *,void *,void *);
void Density(PARTICLE *,int,NN *,SMF *);
void DensitySym(PARTICLE *,int,NN *,SMF *);

#define SMX_MEANVEL				2
void initMeanVel(void *,void *);
void combMeanVel(void *,void *,void *);
void MeanVel(PARTICLE *,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,int,NN *,SMF *);

#define SMX_FOF			25
void initGroupIds(void *,void *p);
void initGroupMerge(void *,void *g);
void combGroupMerge(void *,void *g1, void *g2);
void initGroupBins(void *,void *b);
void combGroupBins(void *,void *b1, void *b2);

#define SMX_RELAXATION		26
void AddRelaxation(PARTICLE *,int,NN *,SMF *);

#ifdef SYMBA
#define SMX_SYMBA               27
void DrmininDrift(PARTICLE *,int,NN *,SMF *);
#endif

#endif










