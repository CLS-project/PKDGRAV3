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

#ifndef BASETYPE_H
#define BASETYPE_H
#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <stdint.h>
#include "gravity/ilp.h"
#include "gravity/ilc.h"
#include "SPH/SPHOptions.h"
#include "units.h"

#define IORDERBITS 43
#define IORDERMAX ((((uint64_t) 1)<<IORDERBITS)-1)

#define IRUNGBITS 6
#define IRUNGMAX ((1<<IRUNGBITS)-1)

typedef float vel_t;
typedef double pos_t;

/*
** Costs: ADD/SUB/MUL/AND/CMP 1 cycle
** RSQRT: 7 cycles (latency)
** DIV: 35 cycles (latency)
*/
#define COST_FLOP_PP 53 // +-*: 38 AND/CMP:8 rsqrt:1
#define COST_FLOP_PC 215 // +-*: 206 AND/CMP:2 rsqrt:1
#define COST_FLOP_EWALD 386 // +-*: 306 AND/CMP:3 div:2 rsqrt:1
#define COST_FLOP_HLOOP 62 // +-*:46  AND/CMP:16 div: rsqrt:0
#define COST_FLOP_SOFT 15
#define COST_FLOP_OPEN 97 // +-*:55  AND/CMP:42

/*
** This is an important base type. Alter with care, or even better, leave it alone.
*/
typedef struct {
    int32_t  iPid;      /* A processor */
    int32_t  iIndex;    /* Index of item on the processor */
} remoteID;

/* Regular particle with order and all the goodies */
typedef struct particle {
uint64_t  uRung      :  IRUNGBITS;
    uint64_t  bMarked    :  1;
uint64_t  uNewRung   :  IRUNGBITS;  /* Optional with bNewKDK + bMemUnordered */
    uint64_t  iClass     :  8;          /* Optional with bMemUnordered */
uint64_t  iOrder     :  IORDERBITS; /* Optional with bMemUnordered */
} PARTICLE;
static_assert(sizeof(PARTICLE)==sizeof(uint64_t));

/* Abbreviated particle header with group id */
#define IGROUPBITS (32-IRUNGBITS-1)
#define IGROUPMAX ((1<<IGROUPBITS)-1)

typedef struct uparticle {
uint32_t  uRung      :  IRUNGBITS;
    uint32_t  bMarked    :  1;
uint32_t  iGroup     :  IGROUPBITS;
} UPARTICLE;
static_assert(sizeof(UPARTICLE)==sizeof(uint32_t));

#define PP_CUDA_MEMORY_LIMIT (2*1024*1024)

typedef struct {
    float r[3];
    float a[3];
    float fSmooth2;
    float fDensity;
    /* SPH fields follow */
    float fBall;
    float Omega;
    float v[3];
    float rho; /* fDensity above is used for different stuff, calculated with different kernel etc */
    float P;
    float c;
    int32_t species;
    int isTooLarge;
    /*    float v[3];*/
    /*    float fMass;*/
    /*    float fSoft;*/
} PINFOIN;

typedef struct {
    float a[3];
    float fPot;
    float dirsum, normsum;
    float rhopmax;
    /* SPH fields follow */
    float rho, drhodfball, nden, dndendfball, fBall, nSmooth;
    float uDot, divv, dtEst, maxRung;
} PINFOOUT;


#if defined(USE_SIMD)
    typedef float ewaldFloatType;
#else
    typedef double ewaldFloatType;
#endif
typedef union {
    ewaldFloatType *f;
#if defined(USE_SIMD) && !defined(__CUDACC__)
    v_sf *p;
#endif
} ewaldFloat;

typedef struct {
    ewaldFloat hx,hy,hz;
    ewaldFloat hCfac,hSfac;
} EwaldTable;
struct EwaldVariables {
    momFloat r[3]; /* Center of mass of the box */
    MOMC mom; /* moment of the box */
    momFloat fEwCut2,fInner2,alpha,ialpha,alpha2,k1,ka,Lbox;
    momFloat Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
    int nMaxEwhLoop;
    int nEwLoopInner, nEwhLoop;
    int nReps,nEwReps;
};

struct pkdTimestepParameters {
    double dDelta, dEta;
    double dAccFac, dRhoFac, dPreFacRhoLoc;
    uint8_t uRungLo,uRungHi,uMaxRung;
    uint8_t bGravStep;
    int iTimeStepCrit;
    int nPartRhoLoc;
    double nPartColl;
    double dEccFacMax;
};

struct pkdKickParameters {
    int bKickClose, bKickOpen;
    vel_t dtClose[IRUNGMAX+1];
    vel_t dtOpen[IRUNGMAX+1];
    vel_t dtPredDrift[IRUNGMAX+1];
};

struct pkdLightconeParameters {
    double dtLCDrift[IRUNGMAX+1];
    double dtLCKick[IRUNGMAX+1];
    double dLookbackFac;
    double dLookbackFacLCP;
    double dBoxSize;
    int bLightConeParticles;
};

struct inChemCompInit {
    double dInitialH;
#ifdef HAVE_HELIUM
    double dInitialHe;
#endif
#ifdef HAVE_CARBON
    double dInitialC;
#endif
#ifdef HAVE_NITROGEN
    double dInitialN;
#endif
#ifdef HAVE_OXYGEN
    double dInitialO;
#endif
#ifdef HAVE_NEON
    double dInitialNe;
#endif
#ifdef HAVE_MAGNESIUM
    double dInitialMg;
#endif
#ifdef HAVE_SILICON
    double dInitialSi;
#endif
#ifdef HAVE_IRON
    double dInitialFe;
#endif
#ifdef HAVE_METALLICITY
    double dInitialMetallicity;
#endif
};

struct inEndTimestep {
    UNITS units;
#ifdef GRACKLE
    char achCoolingTable[256];
#endif
    double dTime;
    double dDelta;
    double dConstGamma;
    double dTuFac;
#ifdef STAR_FORMATION
    double dSFMinOverDensity;
#endif
#ifdef COOLING
    double dCoolingFloorDen;
    double dCoolingFlooru;
#endif
#ifdef EEOS_POLYTROPE
    double dEOSPolyFloorIndex;
    double dEOSPolyFloorDen;
    double dEOSPolyFlooru;
#endif
#ifdef EEOS_JEANS
    double dEOSNJeans;
#endif
#ifdef BLACKHOLES
    double dBHRadiativeEff;
#endif
#ifdef STELLAR_EVOLUTION
    int bChemEnrich;
#endif
};

#ifdef STAR_FORMATION
struct inStarForm {
    double dDenMin;
    double dDelta;
    double dTime;
    double dScaleFactor;

    double dSFindexKS;
    double dSFnormalizationKS;
    double dConstGamma;
    double dSFGasFraction;
    double dSFThresholdu;
    double dSFEfficiency;

    double dEOSPolyFloorIndex;
    double dEOSPolyFloorDen;
    double dEOSPolyFlooru;

#ifdef FEEDBACK
    double dSNFBEfficiency;
    double dSNFBMaxEff;
    double dSNFBEffnH0;
    double dSNFBEffIndex;
#endif

#ifdef STELLAR_EVOLUTION
    double dSNIaMaxMass;
    double dCCSNMinMass;
    double dCCSNMaxMass;
#endif
};
#endif


/*
** Accumulates the work for a set of particles
*/
typedef struct {
    PARTICLE **pPart;
    uint32_t *iPart;
    PINFOIN *pInfoIn;
    PINFOOUT *pInfoOut;
    double dFlop;
    double c[3];
    int nP;
    int nRefs;
    void *ctx;

    struct pkdKickParameters *kick;
    struct pkdLightconeParameters *lc;
    struct pkdTimestepParameters *ts;
    double dirLsum, normLsum;

    double dTime;

    double dFlopSingleCPU;
    double dFlopSingleGPU;
    double dFlopDoubleCPU;
    double dFlopDoubleGPU;

    SPHOptions *SPHoptions;
    ilpList *ilp;
    int bGravStep;
} workParticle;

#define EWALD_ALIGN 64
#define EWALD_MASK (EWALD_ALIGN-1)
typedef struct {
    momFloat X[EWALD_ALIGN];
    momFloat Y[EWALD_ALIGN];
    momFloat Z[EWALD_ALIGN];
} gpuEwaldInput;

typedef struct {
    momFloat X[EWALD_ALIGN];
    momFloat Y[EWALD_ALIGN];
    momFloat Z[EWALD_ALIGN];
    momFloat Pot[EWALD_ALIGN];
    momFloat FlopDouble[EWALD_ALIGN];
} gpuEwaldOutput;
#endif
