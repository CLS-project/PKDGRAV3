#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED
#define SMOOTHFCN_H_MODULE_ID "$Id$"

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
    int nBins;
    int iCenterType;
    int bLogBins;
    FLOAT binFactor;
    FLOAT fMinRadius;
    /* Gas */
    double alpha;
    double beta;
    double gamma;
    int bViscosityLimiter;
    /* diffusion */
    double dMetalDiffusionCoeff;
    double dThermalDiffusionCoeff;
    int bConstantDiffusion;
    /* star form */
    /*
    double dMinMassFrac;
    int bShortCoolShutoff;
    int bSNTurnOffCooling;
    int bSmallSNSmooth;
    double dSecUnit;
    double dGmUnit;
    struct snContext sn;
    */
    /* end starform */   
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

#define SMX_DIVV				3
void initDivv(void *,void *);
void combDivv(void *,void *,void *);
void Divv(PARTICLE *,int,NN *,SMF *);
void DivvSym(PARTICLE *,int,NN *,SMF *);

#define SMX_VELDISP2				4
void initVelDisp2(void *,void *);
void combVelDisp2(void *,void *,void *);
void VelDisp2(PARTICLE *,int,NN *,SMF *);
void VelDisp2Sym(PARTICLE *,int,NN *,SMF *);

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










