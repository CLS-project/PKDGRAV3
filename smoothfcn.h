/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"

typedef struct smfParameters {
    int bComove;
    double dTime;
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
    double binFactor;
    double fMinRadius;
    /* Gas */
    double alpha;
    double beta;
    double gamma;
    double dDelta;
    double dEtaCourant;
    int iViscosityLimiter;
    /* diffusion */
    int iDiffusion;
    double dMetalDiffusionCoeff;
    double dThermalDiffusionCoeff;
    /* star form */
    double SFdESNPerStarMass;
    double SFdtCoolingShutoff;
    double SFdtFeedbackDelay;
    double SFdMassLossPerStarMass;
    double SFdZMassPerStarMass;
    double SFdFBFac;
    /* end starform */   
    PKD pkd; /* useful for diagnostics, etc. */
    remoteID hopParticleLink;
    int bDone;
    float *pfDensity;
    } SMF;


typedef struct pqNode {
    struct pqNode *pqLoser;
    struct pqNode *pqFromInt;
    struct pqNode *pqFromExt;
    struct pqNode *pqWinner;	/* Only used when building initial tree */
    PARTICLE *pPart;
    double fDist2;
    double dx;
    double dy;
    double dz;
    int iIndex;
    int iPid;
    } PQ;


typedef PQ NN;

#define PQ_INIT(pq,n)\
{\
	int j;\
	if ((n) == 1) {\
		(pq)[0].pqFromInt = NULL;\
		(pq)[0].pqFromExt = NULL;\
		}\
	else {\
	    for (j=0;j<(n);++j) {\
		    if (j < 2) (pq)[j].pqFromInt = NULL;\
		    else (pq)[j].pqFromInt = &(pq)[j>>1];\
		    (pq)[j].pqFromExt = &(pq)[(j+(n))>>1];\
		    }\
	    }\
    }


#define PQ_BUILD(pq,n,q)\
{\
    int i,j;\
	PQ *t,*lt;\
	for (j=(n)-1;j>0;--j) {\
		i = (j<<1);\
		if (i < (n)) t = (pq)[i].pqWinner;\
		else t = &(pq)[i-(n)];\
		++i;\
		if (i < (n)) lt = (pq)[i].pqWinner;\
		else lt = &(pq)[i-(n)];\
		if (t->fDist2 < lt->fDist2) {\
			(pq)[j].pqLoser = t;\
			(pq)[j].pqWinner = lt;\
			}\
		else {\
			(pq)[j].pqLoser = lt;\
			(pq)[j].pqWinner = t;\
			}\
		}\
    if ((n) == 1) (q) = (pq);\
	else (q) = (pq)[1].pqWinner;\
	}


#define PQ_REPLACE(q)\
{\
	PQ *t,*lt;\
	t = (q)->pqFromExt;\
	while (t) {\
		if (t->pqLoser->fDist2 > (q)->fDist2) {\
			lt = t->pqLoser;\
			t->pqLoser = (q);\
			(q) = lt;\
			}\
		t = t->pqFromInt;\
		}\
	}


#define SMX_NULL                            0
void NullSmooth(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DENSITY				1
void initDensity(void *,void *);
void combDensity(void *,void *,void *);
void Density(PARTICLE *,float fBall,int,NN *,SMF *);
void DensitySym(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_MEANVEL				2
void initMeanVel(void *,void *);
void combMeanVel(void *,void *,void *);
void MeanVel(PARTICLE *,float fBall,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DIVV				3
void initDivv(void *,void *);
void combDivv(void *,void *,void *);
void Divv(PARTICLE *,float fBall,int,NN *,SMF *);
void DivvSym(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_VELDISP2				4
void initVelDisp2(void *,void *);
void combVelDisp2(void *,void *,void *);
void VelDisp2(PARTICLE *,float fBall,int,NN *,SMF *);
void VelDisp2Sym(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DENDVDX				5
void initDenDVDX(void *,void *);
void combDenDVDX(void *,void *,void *);
void DenDVDX(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_SPHFORCES				6
void initSphForcesParticle(void *,void *);
void initSphForces(void *,void *);
void combSphForces(void *,void *,void *);
void SphForces(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DIST_DELETED_GAS                    7
void initDistDeletedGas(void *,void *p1);
void combDistDeletedGas(void *,void *p1,void *p2);
void DistDeletedGas(PARTICLE *,float fBall, int, NN *, SMF *);

#define SMX_DIST_SN_ENERGY                      8
void initDistSNEnergy(void *,void *p1);
void combDistSNEnergy(void *,void *p1,void *p2);
void DistSNEnergy(PARTICLE *p, float, int, NN *, SMF *);

#define SMX_PRINTNN                            9
void PrintNN(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_FOF			25
void initGroupIds(void *,void *p);
void initGroupMerge(void *,void *g);
void combGroupMerge(void *,void *g1, void *g2);
void initGroupBins(void *,void *b);
void combGroupBins(void *,void *b1, void *b2);

#define SMX_RELAXATION		26
void AddRelaxation(PARTICLE *,float fBall,int,NN *,SMF *);

#ifdef SYMBA
#define SMX_SYMBA               27
void DrmininDrift(PARTICLE *,float fBall,int,NN *,SMF *);
#endif

#define SMX_DENSITY_F1          28
void DensityF1(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DENSITY_M3          29
void DensityM3(PARTICLE *,float fBall,int,NN *,SMF *);
#define SMX_GRADIENT_M3         30
void LinkGradientM3(PARTICLE *,float fBall,int,NN *,SMF *);
#define SMX_HOP_LINK            31
void LinkHopChains(PARTICLE *,float fBall,int,NN *,SMF *);
#endif










