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
    int nSmooth;
    int iMaxRung;
    int bComove;
    int bDoGravity;
    UNITS units;
    double dTime;
    double H;
    double a;
    double dDeltaT;
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
    double gamma;
    double dDelta;
    double dEtaCourant;
    PKD pkd; /* useful for diagnostics, etc. */
    remoteID hopParticleLink;
    int bDone;
    float *pfDensity;
    /* IA: Meshless hydro */
    int bMeshlessHydro;
    int bIterativeSmoothingLength;
    int bUpdateBall;
    int nBucket;
    double dCFLacc;
    double dConstGamma;
    double dhMinOverSoft;
    double dNeighborsStd;
    struct eEOSparam eEOS;
#ifdef FEEDBACK
    double dSNFBDu;
    double dCCSNFBDelay;
    double dCCSNFBSpecEnergy;
    double dSNIaFBDelay;
    double dSNIaFBSpecEnergy;
#endif
#ifdef BLACKHOLES
    double dBHFBEff;
    double dBHFBEcrit;
    double dBHAccretionEddFac;
    double dBHAccretionAlpha;
    double dBHAccretionCvisc;
    int bBHFeedback;
    int bBHAccretion;
#endif
#ifdef STELLAR_EVOLUTION
    double dWindSpecificEkin;
    double dSNIaNorm;
    double dSNIaScale;
    double dCCSNMinMass;
#endif
} SMF;

typedef struct pqNode {
    struct pqNode *pqLoser;
    struct pqNode *pqFromInt;
    struct pqNode *pqFromExt;
    struct pqNode *pqWinner;    /* Only used when building initial tree */
    PARTICLE *pPart;
    double fDist2;
    blitz::TinyVector<double,3> dr;
    float fBall;
    int bMarked;
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

#define SMX_DENSITY             1
void initDensity(void *,void *);
void combDensity(void *,void *,const void *);
void Density(PARTICLE *,float fBall,int,NN *,SMF *);
void DensitySym(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_HYDRO_DENSITY     40
#define SMX_HYDRO_DENSITY_FINAL 41
#define SMX_HYDRO_GRADIENT    42
#define SMX_HYDRO_FLUX        43
#define SMX_HYDRO_STEP        44
#define SMX_HYDRO_FLUX_VEC    45

#define SMX_SN_FEEDBACK       50

#define SMX_BH_MERGER         55
#define SMX_BH_DRIFT          56
#define SMX_BH_STEP           57

#ifdef STELLAR_EVOLUTION
    #define SMX_CHEM_ENRICHMENT   60
#endif

#define SMX_PRINTNN                            9
void PrintNN(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DENSITY_F1          28
void DensityF1(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_DENSITY_M3          29
void DensityM3(PARTICLE *,float fBall,int,NN *,SMF *);
#define SMX_GRADIENT_M3         30
void LinkGradientM3(PARTICLE *,float fBall,int,NN *,SMF *);
#define SMX_HOP_LINK            31
void LinkHopChains(PARTICLE *,float fBall,int,NN *,SMF *);

#define SMX_BALL                33
void initBall(void *,void *);
void BallSmooth(PARTICLE *,float fBall,int,NN *,SMF *);
#endif
