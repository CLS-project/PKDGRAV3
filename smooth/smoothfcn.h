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
    int bStricterSlopeLimiter;
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

enum SMOOTH_TYPE {
    SMX_NULL,
    SMX_DENSITY,
    SMX_DENSITY_F1,
    SMX_DENSITY_M3,
    SMX_GRADIENT_M3,
    SMX_HOP_LINK,
    SMX_BALL,
    SMX_PRINTNN,
    SMX_HYDRO_DENSITY,
    SMX_HYDRO_DENSITY_FINAL,
    SMX_HYDRO_GRADIENT,
    SMX_HYDRO_FLUX,
    SMX_HYDRO_STEP,
    SMX_HYDRO_FLUX_VEC,
    SMX_SN_FEEDBACK,
    SMX_BH_MERGER,
    SMX_BH_EVOLVE,
    SMX_BH_STEP,
    SMX_BH_GASPIN,
    SMX_CHEM_ENRICHMENT,
};

extern void init_shared_gas_decorator(PKD pkd, void *dst);
extern void init_shared_bh_decorator(PKD pkd, void *dst);

// Smooth operates on particles, so we use a particleCACHEhelper.
// Additionally we can call a function to initialize the local data
// before we open the cache. Most often this is not needed, but
// occasionally we loop over the local particles calling init.
class SmoothHelper : public particleCACHEhelper {
public:
    SmoothHelper(pkdContext *pkd,bool bModify=false) : particleCACHEhelper(pkd,bModify) {}
    virtual void initialize_local() {}
};

template<typename P,typename F=P>
class SmoothHelperPacked : public SmoothHelper {
protected:
    virtual void    pack(void *dst, const void *src) = 0;
    virtual void  unpack(void *dst, const void *src, const void *key) = 0;
public:
    using pack_type = P;
    using flush_type = F;
    explicit SmoothHelperPacked(pkdContext *pkd, bool bModify=false) : SmoothHelper(pkd,bModify) {
        nPackSize = sizeof(P);
        nFlushSize = sizeof(F);
    }
};

// This is a simple cache helper that does not allow symmetric smoothing.
// If bSymmetric is true, it will assert. For symmetric smoothing, the
// helper requires at least a combine function.
class SmoothHelperAsymmetric : public SmoothHelper {
public:
    SmoothHelperAsymmetric(pkdContext *pkd,bool bSymmetric=false)  : SmoothHelper(pkd,false) {
        assert(!bSymmetric);
    }
};

//#define SMX_NULL                            0
void NullSmooth(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
class SmoothHelperNull : public SmoothHelper {
public:
    SmoothHelperNull(pkdContext *pkd) : SmoothHelper(pkd) {}
};

//#define SMX_DENSITY             1
void Density(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
void DensitySym(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
class SmoothHelperDensity : public SmoothHelper {
    virtual void init(void *pIn) override final;
    virtual void combine(void *dst, const void *src, const void *key) override final;
    virtual void initialize_local() override final;
public:
    SmoothHelperDensity(pkdContext *pkd,bool bSymmetric);
};

//#define SMX_PRINTNN                            9
void PrintNN(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);

//#define SMX_DENSITY_F1          28
void DensityF1(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
using SmoothHelperDensityF1 = SmoothHelperAsymmetric;

//#define SMX_DENSITY_M3          29
void DensityM3(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
class SmoothHelperDensityM3 : public SmoothHelper {
public:
    SmoothHelperDensityM3(pkdContext *pkd,bool bSymmetric=false) : SmoothHelper(pkd,false) {
        assert(!bSymmetric);
        assert(pkd->particles.present(PKD_FIELD::oGroup));
    }
};
//#define SMX_GRADIENT_M3         30
void LinkGradientM3(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
using SmoothHelperGradientM3 = SmoothHelperAsymmetric;
//#define SMX_HOP_LINK            31
void LinkHopChains(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
using SmoothHelperHopChains = SmoothHelperAsymmetric;

//#define SMX_BALL                33
void BallSmooth(particleStore::ParticleReference &p,float fBall,int,NN *,SMF *);
class SmoothHelperBallSmooth : public SmoothHelper {
    virtual void init(void *pIn) override final;
    virtual void initialize_local() override final;
public:
    SmoothHelperBallSmooth(pkdContext *pkd);
};
#endif
