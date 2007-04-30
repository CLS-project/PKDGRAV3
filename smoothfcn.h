#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

typedef struct smfParameters {
	double H;
	double a;
#ifdef RELAXATION
	double dDeltaT;
#endif
	double dTau2;
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
void initDensity(void *);
void combDensity(void *,void *);
void Density(PARTICLE *,int,NN *,SMF *);
void DensitySym(PARTICLE *,int,NN *,SMF *);

#define SMX_MARKDENSITY			6
void initParticleMarkDensity(void *);
void initMarkDensity(void *);
void combMarkDensity(void *,void *);
void MarkDensity(PARTICLE *,int,NN *,SMF *);
void MarkDensitySym(PARTICLE *,int,NN *,SMF *);

#define SMX_MARKIIDENSITY		16
void initParticleMarkIIDensity(void *);
void initMarkIIDensity(void *);
void combMarkIIDensity(void *,void *);
void MarkIIDensity(PARTICLE *,int,NN *,SMF *);
void MarkIIDensitySym(PARTICLE *,int,NN *,SMF *);

#define SMX_MARK				17
void initMark(void *);
void combMark(void *,void *);

#define SMX_MEANVEL				2
void initMeanVel(void *);
void combMeanVel(void *,void *);
void MeanVel(PARTICLE *,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,int,NN *,SMF *);

#define SMX_FOF			25
void initGroupIds(void *p);
void initGroupMerge(void *g);
void combGroupMerge(void *g1, void *g2);
void initGroupBins(void *b);
void combGroupBins(void *b1, void *b2);

#ifdef RELAXATION	
#define SMX_RELAXATION		26
void AddRelaxation(PARTICLE *,int,NN *,SMF *);
#endif

#endif










