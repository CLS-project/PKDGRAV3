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

#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <stdint.h>
#include <string.h>

#include "mdl.h"
#ifdef USE_CUDA
    #include "cuda/cudautil.h"
#endif
#include "gravity/ilp.h"
#include "gravity/ilc.h"
#include "gravity/cl.h"
#include "gravity/moments.h"
#include "cosmo.h"
#include "units.h"
#include "io/fio.h"
#ifdef USE_GRAFIC
    #include "grafic.h"
#endif
#include "basetype.h"
#include "io/iomodule.h"
#include "SPHOptions.h"
#include "core/bound.h"
#ifdef GRACKLE
    #include <grackle.h>
#endif
#include "chemistry.h"

#define CAST(T,V) reinterpret_cast<T>(V)

typedef uint_fast32_t local_t; /* Count of particles locally (per processor) */
typedef uint_fast64_t total_t; /* Count of particles globally (total number) */

static inline int d2i(double d)  {
    return (int)d;
}

static inline int64_t d2u64(double d) {
    return (uint64_t)d;
}

#define INTEGER_FACTOR 0x80000000u
#define pkdDblToIntPos(pkd,d) (int32_t)((d)*INTEGER_FACTOR)
#define pkdIntPosToDbl(pkd,pos) ((pos)*(1.0/INTEGER_FACTOR))

/*
** Handy type punning macro.
*/
#define UNION_CAST(x, sourceType, destType) \
    (((union {sourceType a; destType b;})x).b)
/*
** The following sort of definition should really be in a global
** configuration header file -- someday...
*/

#define CID_PARTICLE    0
#define CID_CELL    1
#define CID_PARTICLE2   8
#define CID_CELL2   9
#define CID_HEALPIX     7
#define CID_GROUP   2
#define CID_SHRINK      3
#define CID_RM      3
#define CID_BIN     4
#define CID_SHAPES  5
#define CID_PK          2
#define CID_GridLinFx   2
#define CID_GridLinFy   10
#define CID_GridLinFz   11
#define CID_PNG         2
#define CID_SADDLE_BUF  3
#define CID_TREE_ROOT   3

/*
** This is useful for debugging the very-active force calculation.
*/
#define A_VERY_ACTIVE  1


#define MAX_TIMERS      10

/*
** Memory models.  Each is a bit mask that indicates that additional fields should be
** added to the particle structure.
*/
#define PKD_MODEL_VELOCITY     (1<<0)  /* Velocity Required */
#define PKD_MODEL_ACCELERATION (1<<1)  /* Acceleration Required */
#define PKD_MODEL_POTENTIAL    (1<<2)  /* Potential Required */
#define PKD_MODEL_GROUPS       (1<<3)  /* Group profiling */
#define PKD_MODEL_RELAXATION   (1<<5)  /* Trace relaxation */
#define PKD_MODEL_MASS         (1<<6)  /* Mass for each particle */
#define PKD_MODEL_DENSITY      (1<<7)  /* Density for each particle */
#define PKD_MODEL_BALL         (1<<8)  /* Ball for each particle */
#define PKD_MODEL_SOFTENING    (1<<9)  /* Softening for each particle */
#define PKD_MODEL_VELSMOOTH    (1<<10) /* Velocity Smoothing */
#define PKD_MODEL_SPH          (1<<11) /* Sph Fields */
#define PKD_MODEL_STAR         (1<<12) /* Star Fields */
#define PKD_MODEL_PARTICLE_ID  (1<<13) /* Particles have a unique ID */
#define PKD_MODEL_UNORDERED    (1<<14) /* Particles do not have an order */
#define PKD_MODEL_INTEGER_POS  (1<<15) /* Particles do not have an order */
#define PKD_MODEL_BH           (1<<16) /* BH fields */

#define PKD_MODEL_NODE_MOMENT  (1<<24) /* Include moment in the tree */
#define PKD_MODEL_NODE_ACCEL   (1<<25) /* mean accel on cell (for grav step) */
#define PKD_MODEL_NODE_VEL     (1<<26) /* center of mass velocity for cells */
#define PKD_MODEL_NODE_SPHBNDS (1<<27) /* Include 3 extra bounds in tree */

#define PKD_MODEL_NODE_BND     (1<<28) /* Include normal bounds in tree */
#define PKD_MODEL_NODE_VBND    (1<<29) /* Include velocity bounds in tree for phase-space density*/

#define PKD_MODEL_NEW_SPH      (1ULL<<32) /* New Sph Fields */

typedef struct {
    double rscale[3];
    double vscale[3];
} PSMETRIC;

#define PKD_MAX_CLASSES 256
#define MAX_RUNG     63

/*
** Here we define some special reserved nodes. Node-0 is a sentinel or null node, node-1
** is here defined as the ROOT of the local tree (or top tree), node-2 is unused and
** node-3 is the root node "fixed" tree.
*/
#define FIXROOT         3
#define VAROOT          3
#define ROOT        1
#define NRESERVED_NODES MAX_RUNG+1

struct PARTCLASS {
    float       fMass;    /* Particle mass */
    float       fSoft;    /* Current softening */
    FIO_SPECIES eSpecies; /* Species: dark, star, etc. */

    bool operator <(const PARTCLASS &b) const {
        if ( fMass < b.fMass ) return true;
        else if ( fMass > b.fMass ) return false;
        else if ( fSoft < b.fSoft ) return true;
        else if ( fSoft > b.fSoft ) return false;
        else return eSpecies < b.eSpecies;
    }
    bool operator==(const PARTCLASS &b) const {
        return fMass==b.fMass && fSoft==b.fSoft && eSpecies==b.eSpecies;
    }
    PARTCLASS() = default;
    PARTCLASS(float fMass,float fSoft,FIO_SPECIES eSpecies)
        : fMass(fMass), fSoft(fSoft), eSpecies(eSpecies) {}
};
static_assert(std::is_trivial<PARTCLASS>());

typedef struct velsmooth {
    float vmean[3];
    float divv;
    float veldisp2;
} VELSMOOTH;


typedef double myreal;




typedef struct sphfields {
    char *pNeighborList; /* pointer to nearest neighbor list - compressed */
    double vPred[3];

    float c;        /* sound speed */
#ifndef OPTIM_REMOVE_UNUSED
    float u;            /* thermal energy */
    float uPred;    /* predicted thermal energy */
    float uDot;
    float divv;
    float BalsaraSwitch;    /* Balsara viscosity reduction */
    float fMetals;      /* mass fraction in metals, a.k.a, Z - tipsy output variable */

    /* diffusion */
    float diff;
    float fMetalsPred;
    float fMetalsDot;
#endif //OPTIM_REMOVE_UNUSED

    /* IA: B matrix to 'easily' reconstruct faces'. Reminder: it is symmetric */
    double B[6];

    /* IA: Condition number for pathological configurations */
    myreal Ncond;

    /* IA: Gradients */
    myreal gradRho[3];
    myreal gradVx[3];
    myreal gradVy[3];
    myreal gradVz[3];
    myreal gradP[3];

    /* IA: last time this particle's primitve variables were updated */
    myreal lastUpdateTime;
    myreal lastAcc[3];
    myreal lastMom[3];
    myreal lastE;
    myreal lastUint;
    myreal lastHubble; // TODO: Maybe there is a more intelligent way to avoid saving this...
#ifndef USE_MFM
    myreal lastDrDotFrho[3];
#endif
    float lastMass;

    /* IA: normalization factor (Eq 7 Hopkins 2015) at the particle position */
    double omega;

    /* IA: Fluxes */
    myreal Frho;
    myreal Fmom[3];
    myreal Fene;

#ifndef USE_MFM
    double drDotFrho[3];
#endif
    /* IA: Conserved variables */
    double mom[3];
    double E;
    /* IA: Internal energy, which is evolved in parallel and used for updating the pressure if we are in a cold flow */
    double Uint;

#ifdef ENTROPY_SWITCH
    double S;
    double lastS;
    double maxEkin;
#endif

    /* IA: Primitive variables */
    double P;

    /* IA: fBall from the last iteration. Used for the bisection algorithm */
    //float fLastBall;
    /* IA: Number of neighbors correspoding to that fBall */
    //int nLastNeighs;

    /* IA: TODO temporarly */
    //uint8_t uNewRung;

#ifdef STAR_FORMATION
    myreal SFR;
#endif

    float afElemMass[ELEMENT_COUNT];

#ifdef COOLING
    myreal lastCooling;
    float cooling_dudt;
#endif

#ifdef HAVE_METALLICITY
    float fMetalMass;
#endif

#ifdef STELLAR_EVOLUTION
    float afReceivedMom[3];
    float fReceivedMass;
    float fReceivedE;
#endif


#ifdef FEEDBACK
    float fAccFBEnergy;
#endif


    uint8_t uWake;

} SPHFIELDS;

typedef struct newsphfields {
    float Omega;        /* Correction factor */
    float divv;         /* Divergence of v */
    float u;            /* Thermodynamical variable, can be T, A(s) or u */
    float uDot;         /* Derivative of the thermodynamical variable */
} NEWSPHFIELDS;

typedef struct starfields {
    double omega;
#ifdef STELLAR_EVOLUTION
    float afElemAbun[ELEMENT_COUNT]; /* Formation abundances */
    float fMetalAbun;            /* Formation metallicity */
    float fInitialMass;
    float fLastEnrichTime;
    float fLastEnrichMass;
    int iLastEnrichMass;
    float fNextEnrichTime;
    struct {
        int oZ;
        float fDeltaZ;
    } CCSN, AGB, Lifetime;
    float fSNIaOnsetTime;
#endif

    float fTimer;  /* Time of formation */
    float fSNEfficiency;
    int hasExploded; /* Has exploded as a supernova? */
} STARFIELDS;

typedef struct blackholefields {
    PARTICLE *pLowPot;
    double omega;
    double dInternalMass;
    double newPos[3];
    double lastUpdateTime;
    double dAccretionRate;
    double dEddingtonRatio;
    double dFeedbackRate;
    double dAccEnergy;
    float fTimer;    /* Time of formation */
} BHFIELDS;

#ifdef OPTIM_UNION_EXTRAFIELDS
typedef union extrafields {
    SPHFIELDS sph;
    STARFIELDS star;
    BHFIELDS bh;
} EXTRAFIELDS;
#endif

typedef struct partLightCone {
    float pos[3];
    float vel[3];
#ifdef POTENTIAL_IN_LIGHTCONE
    float pot;
#endif
} LIGHTCONEP;

/*
** General partition macro
** LT,LE: Compare less-than/less-than or equal
** ii,dj: Increment i and decrement j
** SWAP: Swap the i'th and j'th element
** LOWER,UPPER: comparison predicates
** e.g.,
** PARTICLE *pi = pkdParticle(pkd,i);
** PARTICLE *pj = pkdParticle(pkd,j);
**    PARTITION(pi<pj,pi<=pj,
**              pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
**              pkdSwapParticle(pkd,pi,pj),
**          pi->r[d] >= fSplit,pj->r[d] < fSplit);
** When finished, the 'I' variable points to the first element of
** the upper partition (or one past the end).
** NOTE: Because this function supports tiled data structures,
**       the LT, followed by "if (LE)" needs to remain this way.
*/
#define PARTITION(LT,LE,INCI,DECJ,SWAP,LOWER,UPPER) \
    {                           \
    while ((LT) && (LOWER)) { INCI; }           \
    if ((LE) && (LOWER)) { INCI; }          \
    else {                      \
    while ((LT) && (UPPER)) { DECJ; }       \
    while (LT) {                    \
            { SWAP; }               \
            do { DECJ; } while (UPPER); \
            do { INCI; } while (LOWER);     \
        }                       \
    }                       \
    }

#define NEW_STACK(S,inc) \
    int S ## __s=0, S ## __ns=inc, S ## __inc=inc; \
    do { S = malloc((S ## __ns)*sizeof(*S)); assert(S != NULL); } while (0)
#define FREE_STACK(S) do { free(S); } while (0)

#define PUSH(S,v) do { S[(S ## __s)++] = (v); } while (0)
#define POP(S) (S[--(S ## __s)])
#define STACK_EMPTY(S) (S ## __s == 0)
#define EXTEND_STACK(S) do { \
    if ( (S ## __s)+1 >= (S ## __ns) ) { \
        assert( (S ## __s)+1 == (S ## __ns) ); \
        (S ## __ns) += (S ## __inc); \
        S = realloc(S,(S ## __ns)*sizeof(*S)); \
        } \
} while (0)
#define CLEAR_STACK(S) do { S ## __s=0; } while (0)
#define STACK_SIZE(S) (S ## __s)

/* This macro give the id (processor) and index of the two child cells */
#define pkdGetChildCells(c,id,idLower,idxLower,idUpper,idxUpper)    \
    do {                                \
        idxLower = c->iLower;   /* This is always true */       \
    idLower = idUpper = id; /* Default */               \
    if (c->bTopTree) {                      \
        /* We may point off node, but then nodes are adjacent */    \
        if (c->bRemote) {                       \
        idLower = c->pLower;                    \
        idUpper = idLower;                  \
        idxUpper = idxLower+1;                  \
        }                           \
        else {                          \
        idxUpper = c->pUpper;                   \
        }                           \
        }                               \
    else { idxUpper = idxLower+1; }                 \
    } while(0)                          \

typedef struct kdNode {
    int pLower;          /* also serves as thread id for the LTT */
    int pUpper;          /* pUpper < 0 indicates no particles in tree! */
    uint32_t iLower;         /* Local lower node (or remote processor w/bRemote=1) */
    uint32_t iDepth     : 15;
    uint32_t bGroup     : 1;
    uint32_t iSplitDim  : 2;
    uint32_t uMinRung   : 6;
    uint32_t uMaxRung   : 6;
    uint32_t bTopTree   : 1; /* This is a top tree node: pLower,pUpper are node indexes */
    uint32_t bRemote    : 1; /* children are remote */
    float bMax;
    float fSoft2;
#if SPHBALLOFBALLS
    float fBoBr2;       /* Ball of Balls radius squared */
    float fBoBxCenter;
    float fBoByCenter;
    float fBoBzCenter;
#endif
#if SPHBOXOFBALLS
    float fBoBxMin;
    float fBoBxMax;
    float fBoByMin;
    float fBoByMax;
    float fBoBzMin;
    float fBoBzMax;
#endif
    uint64_t bHasMarked : 1;         /* flag if node has a marked particle, there are still 31 bit left*/
} KDN;

typedef struct sphBounds {
    struct minmaxBound {
        double min[3];
        double max[3];
    } A,B,BI;
} SPHBNDS;


#define NMAX_OPENCALC   1000

#define MAXSIDE(fMax,b) {\
    if ((fMax)[0] > (fMax)[1]) {\
    if ((fMax)[0] > (fMax)[2]) b = 2.0*(fMax)[0];\
    else b = 2.0*(fMax)[2];\
    }\
    else {\
    if ((fMax)[1] > (fMax)[2]) b = 2.0*(fMax)[1];\
    else b = 2.0*(fMax)[2];\
    }\
    }

#define MINSIDE(fMax,b) {\
    if ((fMax)[0] < (fMax)[1]) {\
    if ((fMax)[0] < (fMax)[2]) b = 2.0*(fMax)[0];\
    else b = 2.0*(fMax)[2];\
    }\
    else {\
    if ((fMax)[1] < (fMax)[2]) b = 2.0*(fMax)[1];\
    else b = 2.0*(fMax)[2];\
    }\
    }

#define MINDIST(bnd,pos,min2) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;                 \
    for (BND_j=0;BND_j<3;++BND_j) {                 \
    BND_dMin = fabs((bnd)->fCenter[BND_j] - (pos)[BND_j]) - (bnd)->fMax[BND_j]; \
    if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin;          \
    }\
    }
#define MAXDIST(bnd,pos,max2) {                 \
    double BND_dMax;                            \
    int BND_j;                              \
    (max2) = 0;                                 \
    for (BND_j=0;BND_j<3;++BND_j) {                     \
    BND_dMax = fabs((bnd)->fCenter[BND_j] - (pos)[BND_j]) + (bnd)->fMax[BND_j];     \
    (max2) += BND_dMax*BND_dMax;                    \
    }                               \
    }

#define CALCOPEN(pkdn,minside) {                    \
        double CALCOPEN_d2 = 0;                     \
    double CALCOPEN_b;                      \
    const BND CALCOPEN_bnd = pkdNodeGetBnd(pkd, pkdn);      \
    double CALCOPEN_r[3];                       \
    pkdNodeGetPos(pkd, (pkdn), CALCOPEN_r);             \
    MAXDIST(&CALCOPEN_bnd,CALCOPEN_r,CALCOPEN_d2)           \
    MAXSIDE(CALCOPEN_bnd.fMax,CALCOPEN_b);              \
    if (CALCOPEN_b < minside) CALCOPEN_b = minside;         \
    if (CALCOPEN_b*CALCOPEN_b < CALCOPEN_d2) CALCOPEN_b = sqrt(CALCOPEN_d2); \
    (pkdn)->bMax = CALCOPEN_b;                  \
    }

#if (0)
#define CALCOPEN(pkdn) {                        \
        double CALCOPEN_d2 = 0;                     \
    const BND CALCOPEN_bnd = pkdNodeGetBnd(pkd, pkdn);      \
    MAXDIST(&CALCOPEN_bnd,(pkdn)->r,CALCOPEN_d2)            \
        CALCOPEN_d2 = sqrt(CALCOPEN_d2);      \
        if (CALCOPEN_d2 < (pkdn)->bMax) (pkdn)->bMax = CALCOPEN_d2;   \
    }
#endif

/*
** Components required for tree walking.
*/
typedef struct CheckStack {
    clList *cl;
    LOCR L;
    float dirLsum;
    float normLsum;
    int iNodeIndex;
} CSTACK;

/*
** components required for time-step calculation (only grav.c)
*/

typedef struct RhoLocalArray {
    double d2;
    double m;
} RHOLOCAL;

#if defined(USE_SIMD) && !defined(__CUDACC__)
typedef struct {
    struct PMOMC {
        v_df m;
        v_df xx,yy,xy,xz,yz;
        v_df xxx,xyy,xxy,yyy,xxz,yyz,xyz;
        v_df xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
        v_df zz;
        v_df xzz,yzz,zzz;
        v_df xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
    } ewm;
    struct PEWALDVARS {
        v_df fEwCut2,fInner2,alpha,alpha2,ialpha,k1,ka;
        v_df Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
    } ewp;
} ewaldSIMD;
#endif

/*
** This is the temporary group table used when Grasshopping.
** We eventually contruct a proper table.
*/
typedef remoteID GHtmpGroupTable;

typedef struct {
    uint64_t nTotal;      /* Total particles in this group */
    remoteID id;          /* Owner (or myself) */
    remoteID rmt;
    uint32_t iGlobalId;   /* Global unique group id */
    uint32_t nLocal;      /* Local to this processor */
    uint32_t iTreeRoot;   /* Our local tree root */
    uint32_t iAllRoots;
    uint16_t nRemote;     /* Number of remote partners */
    uint8_t  bNeedGrav : 1;
    uint8_t  bComplete : 1;
    float fMass;
    float fRMSRadius;
    double dEnergy;
    double rref[3];
    double ravg[3];
    double rmin[3];
    double rmax[3];
    double rcom[3];
    double vcom[3];
} HopGroupTable;

typedef struct {
    remoteID key;
    remoteID name;
    uint32_t iLink;
} FOFRemote;

typedef struct {
    float rPot[3];
    float minPot;
    float rcen[3];
    float rcom[3];
    float vcom[3];
    float angular[3];
    float inertia[6];
    float sigma;
    float rMax;
    float fMass;
    float fEnvironDensity0;
    float fEnvironDensity1;
    float rHalf;
    int nBH;
    int nStar;
    int nGas;
    int nDM;
} TinyGroupTable;

/* IA: For the BH seeding, we need a small part of the FoF information
 *  to survive for a while
 */
typedef struct {
    double rPot[3];
    float fMass;
    int nBH;
} VeryTinyGroupTable;

//typedef struct {
//    } SmallGroupTable;

struct remote_root_id {
    int iPid;
    int iLocalGroupId;
    int iLocalRootId;
};

struct saddle_point_group {
    int iGlobalId;
    int iLocalId;
    int iPid;
    double fDensity;
};

struct saddle_point {
    /* Information about the particle that is the saddle point */
    double fDensity;
    int iLocalId;
    int iPid;

    /* The group that owns the saddle point */
    struct saddle_point_group owner;
    /* The other group joined by the saddle point */
    struct saddle_point_group nbr;
    struct saddle_point_group parent;
};

struct saddle_point_buffer {
    /* Destination Group Id. I.e., who should get the following saddle points. */
    int iLocalId;
    int iGlobalId;

    int nSaddlePoints;
    struct saddle_point sp[32];
};

struct saddle_point_list {
    /* Number of saddle points in the list */
    int n;
    /* Size of allocated array */
    int size;

    struct saddle_point *sp;
    struct saddle_point_buffer *buf;
};

struct tree_root_cache_msg {
    struct remote_root_id tr;
    int iPid;
    int iLocalId;
};

struct psGroup {
    int iPid;
    int iLocalId;
    int iGlobalId;
    /*-----Unique-Local-Data-----*/
    int bridge;
    int dup;
    int nLocal;
    int nTreeRoots;
    struct remote_root_id *treeRoots;
    int nSaddlePoints;
    int *sp;
    /*-----Shared-Data-----------*/
    uint64_t nTotal;
    double fDensity;
    double fMass;
    double fRMSRadius;
    double r[3], rcom[3];
    double v[3], vcom[3];
    double fMass_com;
};

struct psGroupTable {
    int nGroups;
    struct psGroup *pGroup;
};

typedef struct groupBin {
    double fRadius;
    int nMembers;
    double fMassInBin;
} FOFBIN;

typedef struct profileBin {
    double dRadius;
    double dMassInBin;
    double dVolume;
    double L[3];
    double vel_radial;
    double vel_radial_sigma;
    double vel_tang_sigma;
    uint64_t nParticles;
} PROFILEBIN;

typedef struct shapesBin {
    double com[3];
    double dMassEnclosed;
    double dInertia[3][3];
    double ell_matrix[3][3];
    double ell_center[3];
} SHAPESBIN;

/*#define PKD_GROUP_SIZE 256*/

typedef struct {
    uint64_t nActiveBelow;
    uint64_t nActiveAbove;
    uint64_t nTotalBelow;
    uint64_t nTotalAbove;
} ORBCOUNT;

typedef struct {
    uint32_t nGrouped;
    uint32_t nUngrouped;
    float fPotential;
} healpixData;

enum PKD_FIELD {
    oPosition,
    oVelocity, /* Three vel_t */
    oAcceleration, /* Three float */
    oPotential, /* One float */
    oGroup, /* One int32 */
    oMass, /* One float */
    oSoft, /* One float */
    oDensity, /* One float */
    oBall, /* One float */
    oSph, /* Sph structure */
    oNewSph, /* NewSph structure */
    oStar, /* Star structure */
    oBH, /* BH structure */
    oRelaxation,
    oVelSmooth,
    oRungDest, /* Destination processor for each rung */
    oParticleID,

    MAX_PKD_FIELD
};

class pkdContext {
public:
    explicit pkdContext(
        mdl::mdlClass *mdl,int nStore,uint64_t nMinTotalStore,uint64_t nMinEphemeral,uint32_t nEphemeralBytes,
        int nTreeBitsLo, int nTreeBitsHi,
        int iCacheSize,int iWorkQueueSize,int iCUDAQueueSize,double *fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,uint64_t nBH,
        uint64_t mMemoryModel, int bLightCone, int bLightConeParticles);
    virtual ~pkdContext();

protected:  // Support for memory models
    int NodeAddStruct    (int n);   // Add a NODE structure: assume double alignment
    int NodeAddDouble    (int n=1); // Add n doubles to the node structure
    int NodeAddFloat     (int n=1); // Add n floats to the node structure
    int NodeAddInt64     (int n=1); // Add n 64-bit integers to the node structure
    int NodeAddInt32     (int n=1); // Add n 32-bit integers to the node structure
    int ParticleAddStruct(int n);   // Add a structure: assume double alignment
    int ParticleAddDouble(int n=1); // Add n doubles to the particle structure
    int ParticleAddFloat (int n=1); // Add n floats to the particle structure
    int ParticleAddInt64 (int n=1); // Add n 64-bit integers to the particle structure
    int ParticleAddInt32 (int n=1); // Add n 32-bit integers to the particle structure
protected:
    PARTICLE *pStorePRIVATE = nullptr;
    PARTICLE *pTempPRIVATE = nullptr;
    int nStore = 0; // Maximum local particles
    int nLocal = 0; // Current number of particles
    size_t nParticleAlign = 0, iParticle32 = 0;
    size_t iParticleSize = 0; // Size (in bytes) of a single PARTICLE
    size_t iTreeNodeSize = 0; // Size (in bytes) of a tree node
    uint32_t nEphemeralBytes = 0; /* per-particle */
public:
    int FreeStore() { return nStore; }
    int Local() { return nLocal; }
    int SetLocal(int n) {return (nLocal=n);}
    int AddLocal(int n) {return (nLocal+=n);}
    auto EphemeralBytes() {return nEphemeralBytes; }
    static constexpr auto MaxNodeSize() { return sizeof(KDN) + 2*sizeof(BND) + sizeof(FMOMR) + 6*sizeof(double) + sizeof(SPHBNDS); }
    auto NodeSize() { return iTreeNodeSize; }
    auto ParticleSize() {return iParticleSize; }
    auto ParticleMemory() { return (ParticleSize() + EphemeralBytes()) * (FreeStore()+1); }
    auto ParticleGet(void *pBase, int i) {
        auto v = static_cast<char *>(pBase);
        return reinterpret_cast<PARTICLE *>(v + ((uint64_t)i)*ParticleSize());
    }
    PARTICLE *ParticleBase() { return pStorePRIVATE; }
    PARTICLE *Particle(int i) { return ParticleGet(ParticleBase(),i); }

    void SaveParticle(PARTICLE *a) { memcpy(pTempPRIVATE,a,ParticleSize()); }
    void LoadParticle(PARTICLE *a) { memcpy(a,pTempPRIVATE,ParticleSize()); }
    void CopyParticle(PARTICLE *a, PARTICLE *b) { memcpy(a,b,ParticleSize()); }

protected:
    KDN **kdNodeListPRIVATE; /* BEWARE: KDN is actually variable length! */
    int nTreeBitsLo;
    int nTreeBitsHi;
    int iTreeMask;
    int nTreeTiles;
    int nTreeTilesReserved;
    int nNodes, nMaxNodes;
    auto TreeBase(int iTile=0) { return kdNodeListPRIVATE[iTile]; }
public:
    auto Nodes() const { return nNodes; }
    /*[[deprecated]]*/ void SetNodeCount(int n) { nNodes = n; }
    auto Node(KDN *pBase,int iNode) {
        return reinterpret_cast<KDN *>(reinterpret_cast<char *>(pBase)+NodeSize()*iNode);
    }
    auto TreeNode(int iNode) {
        return Node(TreeBase(iNode>>nTreeBitsLo),iNode&iTreeMask);
    }
    void ExtendTree();
    size_t TreeMemory() { return nTreeTiles * (1<<nTreeBitsLo) * NodeSize(); }
    auto TreeAlignNode() {
        if (nNodes&1) ++nNodes;
        return nNodes;
    }
    auto TreeAllocNode(int n=1) {
        int iNode = nNodes;
        nNodes += n;
        while (nNodes > nMaxNodes) ExtendTree();
        return iNode;
    }
    void TreeAllocNodePair(int *iLeft, int *iRight) {
        *iLeft = TreeAllocNode();
        *iRight = TreeAllocNode();
    }
    auto TreeAllocRootNode() {
        if ((nNodes&1)==0) ++nNodes;
        return TreeAllocNode();
    }

public:
    mdl::mdlClass *mdl;
    auto Self()    const { return mdl->Self(); }
    auto Threads() const { return mdl->Threads(); }

public:
    int nRejects;
    int nActive;
    uint64_t nRung[IRUNGMAX+1];
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nBH;
    double fPeriod[3];
    int iTopTree[NRESERVED_NODES];
    int nNodesFull;     /* number of nodes in the full tree (including very active particles) */
    BND bnd;
    BND vbnd;
    /*
    ** Light cone variables.
    */
    double lcOffset0[184];
    double lcOffset1[184];
    double lcOffset2[184];
    asyncFileInfo afiLightCone;
    LIGHTCONEP *pLightCone;
    int nLightCone, nLightConeMax;
    int64_t nHealpixPerDomain;
    int64_t nSideHealpix;
    healpixData *pHealpixData;

    std::vector<PARTCLASS> ParticleClasses;
    float fSoftFix;
    float fSoftFac;
    float fSoftMax;
    void *pLite;
    /*
    ** Advanced memory models
    */
    int bNoParticleOrder;
    int bIntegerPosition;
    int oFieldOffset[MAX_PKD_FIELD];

    /*
    ** Advanced memory models - Tree Nodes
    */
    int oNodePosition; /* Three double or int32_t (if bIntegerPosition) */
    int oNodeVelocity; /* Three vel_t */
    int oNodeAcceleration; /* Three doubles */
    int oNodeSoft;
    int oNodeMom; /* an FMOMR */
    int oNodeBnd;
    int oNodeSphBounds; /* Three Bounds */
    int oNodeVBnd; /* Velocity bounds */
#ifdef OPTIM_REORDER_IN_NODES
    int oNodeNgas;
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
    int oNodeNstar;
#endif
    int oNodeNbh;
#endif

    /*
    ** Tree walk variables.
    */
    int nMaxStack;
    CSTACK *S;
    ilpList ilp;
    ilcList ilc;
    ilcList ill;
    clList::free_list clFreeList;
    clList *cl;
    clList *clNew;
    double dFlop;
    double dFlopSingleCPU, dFlopDoubleCPU;
    double dFlopSingleGPU, dFlopDoubleGPU;
    /*
    ** Opening angle table for mass weighting.
    */
    float fiCritTheta;
    /* Potential Energy for when potential is not in the particle */
    double dEnergyU;
    /* also put kinetic energy here to calculate it on the fly */
    double dEnergyT;
    double dEnergyW;
    double dEnergyF[3];
    double dEnergyL[3];

    /*
    ** New activation methods
    */
    uint8_t uMinRungActive;
    uint8_t uMaxRungActive;

    /*
    ** Ewald summation setup.
    */
    struct EwaldVariables ew;
    EwaldTable ewt;
#ifdef USE_SIMD_EWALD
    ewaldSIMD es;
#endif

    struct psGroupTable psGroupTable;

    int nGroups, nLocalGroups;
    struct smGroupArray *ga;
    BND bndInterior;  /* this gets calculated at the start of Fof for now but should be done elsewhere */
    /*
    ** Some variables needed for pkdNewFof().
    */
    uint32_t iRemoteGroup,nMaxRemoteGroups;
    FOFRemote *tmpFofRemote;
    TinyGroupTable *tinyGroupTable;
    VeryTinyGroupTable *veryTinyGroupTable;

    GHtmpGroupTable *tmpHopGroups;
    HopGroupTable *hopGroups;
    int *hopRootIndex;
    int hopSavedRoots;
    remoteID *hopRoots;

    struct saddle_point_list saddle_points;
    int nRm;
    int nMaxRm;
    int nBins;

    FOFBIN *groupBin;

    PROFILEBIN *profileBins;

    CSM csm;

#ifdef COOLING
    // IA: we add here the needed cooling information available to all procs
    struct cooling_function_data *cooling;
    struct cooling_tables *cooling_table;
#endif
#ifdef GRACKLE
    chemistry_data *grackle_data;
    chemistry_data_storage *grackle_rates;
    grackle_field_data *grackle_field;
    code_units *grackle_units;
#endif
#ifdef STELLAR_EVOLUTION
    struct StellarEvolutionData *StelEvolData;
#endif

#ifdef USE_CUDA
    CudaClient *cudaClient;
#endif
#ifdef MDL_FFTW
    MDLFFT fft;
#endif

    SPHOptions SPHoptions;

};
typedef pkdContext *PKD;

#if defined(USE_SIMD) && defined(__SSE2__) && 0
#define pkdMinMax(dVal,dMin,dMax) do {                  \
    (dMin)[0] = _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd((dMin)[0]),_mm_set_sd((dVal)[0]))); \
    (dMin)[1] = _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd((dMin)[1]),_mm_set_sd((dVal)[1]))); \
    (dMin)[2] = _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd((dMin)[2]),_mm_set_sd((dVal)[2]))); \
    (dMax)[0] = _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd((dMax)[0]),_mm_set_sd((dVal)[0]))); \
    (dMax)[1] = _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd((dMax)[1]),_mm_set_sd((dVal)[1]))); \
    (dMax)[2] = _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd((dMax)[2]),_mm_set_sd((dVal)[2]))); \
    } while(0)
#else
#define pkdMinMax(dVal,dMin,dMax) {\
    (dMin)[0] = (dVal)[0] < (dMin)[0] ? (dVal)[0] : (dMin)[0];  \
    (dMin)[1] = (dVal)[1] < (dMin)[1] ? (dVal)[1] : (dMin)[1];      \
    (dMin)[2] = (dVal)[2] < (dMin)[2] ? (dVal)[2] : (dMin)[2];      \
    (dMax)[0] = (dVal)[0] > (dMax)[0] ? (dVal)[0] : (dMax)[0];      \
    (dMax)[1] = (dVal)[1] > (dMax)[1] ? (dVal)[1] : (dMax)[1];      \
    (dMax)[2] = (dVal)[2] > (dMax)[2] ? (dVal)[2] : (dMax)[2];      \
    }
#endif
#define pkdMinMax6(dVal0,dVal1,dMin,dMax) {\
    (dMin)[0] = (dVal0)[0] < (dMin)[0] ? (dVal0)[0] : (dMin)[0];    \
    (dMin)[1] = (dVal0)[1] < (dMin)[1] ? (dVal0)[1] : (dMin)[1];    \
    (dMin)[2] = (dVal0)[2] < (dMin)[2] ? (dVal0)[2] : (dMin)[2];    \
    (dMin)[3] = (dVal1)[0] < (dMin)[3] ? (dVal1)[0] : (dMin)[3];    \
    (dMin)[4] = (dVal1)[1] < (dMin)[4] ? (dVal1)[1] : (dMin)[4];    \
    (dMin)[5] = (dVal1)[2] < (dMin)[5] ? (dVal1)[2] : (dMin)[5];    \
    (dMax)[0] = (dVal0)[0] > (dMax)[0] ? (dVal0)[0] : (dMax)[0];    \
    (dMax)[1] = (dVal0)[1] > (dMax)[1] ? (dVal0)[1] : (dMax)[1];    \
    (dMax)[2] = (dVal0)[2] > (dMax)[2] ? (dVal0)[2] : (dMax)[2];    \
    (dMax)[3] = (dVal1)[0] > (dMax)[3] ? (dVal1)[0] : (dMax)[3];    \
    (dMax)[4] = (dVal1)[1] > (dMax)[4] ? (dVal1)[1] : (dMax)[4];    \
    (dMax)[5] = (dVal1)[2] > (dMax)[5] ? (dVal1)[2] : (dMax)[5];    \
    }

static inline int pkdIsRungRange(PARTICLE *p,uint8_t uRungLo,uint8_t uRungHi) {
    return ((p->uRung >= uRungLo)&&(p->uRung <= uRungHi));
}

static inline int pkdIsRungActive(PKD pkd, uint8_t uRung ) {
    return uRung >= pkd->uMinRungActive && uRung <= pkd->uMaxRungActive;
}
static inline int pkdIsActive(PKD pkd, PARTICLE *p ) {
    return pkdIsRungActive(pkd,p->uRung);
}

/*
** A tree node is of variable size.  The following routines are used to
** access individual fields.
*/
static inline void pkdCopyNode(PKD pkd, KDN *a, KDN *b) {
    memcpy(a,b,pkd->NodeSize());
}
static inline void *pkdNodeField( KDN *n, int iOffset ) {
    char *v = (char *)n;
    /*assert(iOffset);*/ /* Remove this for better performance */
    return (void *)(v + iOffset);
}
static inline FMOMR *pkdNodeMom(PKD pkd,KDN *n) {
    return CAST(FMOMR *,pkdNodeField(n,pkd->oNodeMom));
}
static inline vel_t *pkdNodeVel( PKD pkd, KDN *n ) {
    return CAST(vel_t *,pkdNodeField(n,pkd->oNodeVelocity));
}
static inline float *pkdNodeAccel( PKD pkd, KDN *n ) {
    return CAST(float *,pkdNodeField(n,pkd->oNodeAcceleration));
}
static inline SPHBNDS *pkdNodeSphBounds( PKD pkd, KDN *n ) {
    return CAST(SPHBNDS *,pkdNodeField(n,pkd->oNodeSphBounds));
}

static inline void pkdNodeGetPos(PKD pkd,KDN *n,double *r) {
    if (pkd->bIntegerPosition) {
        int32_t *pr = CAST(int32_t *,pkdNodeField(n,pkd->oNodePosition));
        r[0] = pkdIntPosToDbl(pkd,pr[0]);
        r[1] = pkdIntPosToDbl(pkd,pr[1]);
        r[2] = pkdIntPosToDbl(pkd,pr[2]);
    }
    else {
        double *pr = CAST(double *,pkdNodeField(n,pkd->oNodePosition));
        r[0] = pr[0];
        r[1] = pr[1];
        r[2] = pr[2];
    }
}
static inline void pkdNodeSetPos3(PKD pkd,KDN *n,double x, double y, double z) {
    if (pkd->bIntegerPosition) {
        int32_t *pr = CAST(int32_t *,pkdNodeField(n,pkd->oNodePosition));
        pr[0] = pkdDblToIntPos(pkd,x);
        pr[1] = pkdDblToIntPos(pkd,y);
        pr[2] = pkdDblToIntPos(pkd,z);
    }
    else {
        double *pr = CAST(double *,pkdNodeField(n,pkd->oNodePosition));
        pr[0] = x;
        pr[1] = y;
        pr[2] = z;
    }
}
static inline void pkdNodeSetPos1(PKD pkd,KDN *n,const double *r) {
    pkdNodeSetPos3(pkd,n,r[0],r[1],r[2]);
}

static inline BND *pkdNodeBndPRIVATE( PKD pkd, KDN *n ) {
    return CAST(BND *,pkdNodeField(n,pkd->oNodeBnd));
}
static inline BND pkdNodeGetBnd( PKD pkd, KDN *n ) {
    if (pkd->bIntegerPosition) {
        IBND *ibnd = (IBND *)pkdNodeField(n,pkd->oNodeBnd);
        BND bnd;
        bnd.fCenter[0] = pkdIntPosToDbl(pkd,ibnd->fCenter[0]);
        bnd.fCenter[1] = pkdIntPosToDbl(pkd,ibnd->fCenter[1]);
        bnd.fCenter[2] = pkdIntPosToDbl(pkd,ibnd->fCenter[2]);
        bnd.fMax[0] = pkdIntPosToDbl(pkd,ibnd->fMax[0]);
        bnd.fMax[1] = pkdIntPosToDbl(pkd,ibnd->fMax[1]);
        bnd.fMax[2] = pkdIntPosToDbl(pkd,ibnd->fMax[2]);
        return bnd;
    }
    else return *pkdNodeBndPRIVATE(pkd,n);
}
static inline void pkdNodeSetBnd( PKD pkd, KDN *n, const BND *bnd ) {
    if (pkd->bIntegerPosition) {
        IBND *ibnd = (IBND *)pkdNodeField(n,pkd->oNodeBnd);
        ibnd->fCenter[0] = pkdDblToIntPos(pkd,bnd->fCenter[0]);
        ibnd->fCenter[1] = pkdDblToIntPos(pkd,bnd->fCenter[1]);
        ibnd->fCenter[2] = pkdDblToIntPos(pkd,bnd->fCenter[2]);
        ibnd->fMax[0] = pkdDblToIntPos(pkd,bnd->fMax[0]);
        ibnd->fMax[1] = pkdDblToIntPos(pkd,bnd->fMax[1]);
        ibnd->fMax[2] = pkdDblToIntPos(pkd,bnd->fMax[2]);
    }
    else {
        *pkdNodeBndPRIVATE(pkd,n) = *bnd;
    }
}

static inline void pkdNodeSetBndMinMax( PKD pkd, KDN *n, double *dMin, double *dMax ) {
    BND bnd;
    int j;
    for (j=0; j<3; ++j) {
        bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
        bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
    }
    pkdNodeSetBnd(pkd,n,&bnd);
}

#ifdef OPTIM_REORDER_IN_NODES
static inline int pkdNodeNgas( PKD pkd, KDN *n) {
    return *CAST(int *, pkdNodeField(n, pkd->oNodeNgas));
}
static inline void pkdNodeSetNgas(  PKD pkd, KDN *n, int ngas) {
    *CAST(int *, pkdNodeField(n, pkd->oNodeNgas)) = ngas;
}
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
static inline int pkdNodeNstar( PKD pkd, KDN *n) {
    return *CAST(int *, pkdNodeField(n, pkd->oNodeNstar));
}
static inline void pkdNodeSetNstar(  PKD pkd, KDN *n, int nstar) {
    *CAST(int *, pkdNodeField(n, pkd->oNodeNstar)) = nstar;
}

#endif
static inline int pkdNodeNbh( PKD pkd, KDN *n) {
    return *CAST(int *, pkdNodeField(n, pkd->oNodeNbh));
}
static inline void pkdNodeSetNBH(  PKD pkd, KDN *n, int nbh) {
    *CAST(int *, pkdNodeField(n, pkd->oNodeNbh)) = nbh;
}
#endif

static inline BND *pkdNodeVBnd( PKD pkd, KDN *n ) {
    return CAST(BND *,pkdNodeField(n,pkd->oNodeVBnd));
}

void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize);
//static inline KDN *pkdTopNode(PKD pkd,int iNode) {
//    assert(0); // no top tree now
//    }
/*
** The size of a particle is variable based on the memory model.
** The following three routines must be used instead of accessing pStore
** directly.  pkdParticle will return a pointer to the i'th particle.
** The Size and Base functions are intended for cache routines; no other
** code should care about sizes of the particle structure.
*/
static inline void pkdSwapParticle(PKD pkd, PARTICLE *a, PARTICLE *b) {
    pkd->SaveParticle(a);
    pkd->CopyParticle(a,b);
    pkd->LoadParticle(b);
}

static inline const void *pkdFieldRO( const PARTICLE *p, int iOffset ) {
    const char *v = (const char *)p;
    /*assert(iOffset);*/ /* Remove this for better performance */
    return (const void *)(v + iOffset);
}

static inline void *pkdField( PARTICLE *p, int iOffset ) {
    char *v = (char *)p;
    /*assert(iOffset);*/ /* Remove this for better performance */
    return (void *)(v + iOffset);
}

static inline int32_t *pkdInt32( PARTICLE *p, int iOffset ) {
    char *v = (char *)p;
    return (int32_t *)(v + iOffset);
}

static inline int pkdIsFieldPresent(PKD pkd, enum PKD_FIELD field) {
    assert(field>=0 && field<MAX_PKD_FIELD);
    return pkd->oFieldOffset[field] != 0;
}

static inline int32_t pkdGetGroup( PKD pkd, const PARTICLE *p ) {
    if (pkd->bNoParticleOrder) return ((const UPARTICLE *)p)->iGroup;
    assert(pkd->oFieldOffset[oGroup]);
    return CAST(const int32_t *, pkdFieldRO(p,pkd->oFieldOffset[oGroup]))[0];
}

static inline void pkdSetGroup( PKD pkd, PARTICLE *p, uint32_t gid ) {
    if (pkd->bNoParticleOrder) ((UPARTICLE *)p)->iGroup = gid;
    else if (pkd->oFieldOffset[oGroup]) CAST(int32_t *, pkdField(p,pkd->oFieldOffset[oGroup]))[0] = gid;
}

static inline float pkdDensity( PKD pkd, const PARTICLE *p ) {
    assert(pkd->oFieldOffset[oDensity]);
    return * CAST(const float *, pkdFieldRO(p,pkd->oFieldOffset[oDensity]));
}
static inline void pkdSetDensity( PKD pkd, PARTICLE *p, float fDensity ) {
    if (pkd->oFieldOffset[oDensity]) *CAST(float *, pkdField(p,pkd->oFieldOffset[oDensity])) = fDensity;
}

static inline float pkdBall( PKD pkd, PARTICLE *p ) {
    assert(pkd->oFieldOffset[oBall]);
    return *CAST(float *, pkdField(p,pkd->oFieldOffset[oBall]));
}
static inline void pkdSetBall(PKD pkd, PARTICLE *p, float fBall) {
    if (pkd->oFieldOffset[oBall]) *CAST(float *, pkdField(p,pkd->oFieldOffset[oBall])) = fBall;
}


/* Here is the new way of getting mass and softening */
static inline float pkdMass( PKD pkd, PARTICLE *p ) {
    if ( pkd->oFieldOffset[oMass] ) {
        float *pMass = CAST(float *,pkdField(p,pkd->oFieldOffset[oMass]));
        return *pMass;
    }
    else if (pkd->bNoParticleOrder) return pkd->ParticleClasses[0].fMass;
    else return pkd->ParticleClasses[p->iClass].fMass;
}
static inline float pkdSoft0( PKD pkd, PARTICLE *p ) {
    if ( pkd->oFieldOffset[oSoft] ) {
        float *pSoft = CAST(float *,pkdField(p,pkd->oFieldOffset[oSoft]));
        return *pSoft;
    }
    else if (pkd->bNoParticleOrder) return pkd->ParticleClasses[0].fSoft;
    else return pkd->ParticleClasses[p->iClass].fSoft;
}
static inline float pkdSoft( PKD pkd, PARTICLE *p ) {
    float fSoft;
    if ( pkd->fSoftFix >= 0.0 ) fSoft = pkd->fSoftFix;
    else fSoft = pkdSoft0(pkd,p);
    fSoft *= pkd->fSoftFac;
    if ( fSoft > pkd->fSoftMax ) fSoft = pkd->fSoftMax;
    return fSoft;
}
static inline FIO_SPECIES pkdSpecies( PKD pkd, PARTICLE *p ) {
    if (pkd->bNoParticleOrder) return pkd->ParticleClasses[0].eSpecies;
    else return pkd->ParticleClasses[p->iClass].eSpecies;
}

/*
** Integerized coordinates: signed integer -0x7fffffff to +0x7fffffff
** We assume a periodic box of width 1 so a simple multiple will convert.
** The situation is more complicated with non-periodic boxes, or for boxes
** with a different period so this is not currently supported.
*/
#define pkdPosRaw(pkd,p,d) (CAST(int32_t *,pkdField(p,pkd->oFieldOffset[oPosition]))[d])
#define pkdSetPosRaw(pkd,p,d,v) (CAST(int32_t *,pkdField(p,pkd->oFieldOffset[oPosition]))[d]) = (v)

static inline double pkdPos(PKD pkd,PARTICLE *p,int d) {
    if (pkd->bIntegerPosition) return pkdIntPosToDbl(pkd,CAST(int32_t *,pkdField(p,pkd->oFieldOffset[oPosition]))[d]);
    else return CAST(double *,pkdField(p,pkd->oFieldOffset[oPosition]))[d];
}
static inline void pkdSetPos(PKD pkd,PARTICLE *p,int d,double v) {
    if (pkd->bIntegerPosition) CAST(int32_t *,pkdField(p,pkd->oFieldOffset[oPosition]))[d] = pkdDblToIntPos(pkd,v);
    else CAST(double *,pkdField(p,pkd->oFieldOffset[oPosition]))[d] = v;
}
#define pkdGetPos3(pkd,p,d1,d2,d3) ((d1)=pkdPos(pkd,p,0),(d2)=pkdPos(pkd,p,1),(d3)=pkdPos(pkd,p,2))
#define pkdGetPos1(pkd,p,d) pkdGetPos3(pkd,p,(d)[0],(d)[1],(d)[2])

#if defined(__AVX__) && defined(USE_SIMD)
static inline __m128i pkdGetPosRaw(PKD pkd,PARTICLE *p) {
    return _mm_loadu_si128((__m128i *)pkdField(p,pkd->oFieldOffset[oPosition]));
}
static inline __m256d pkdGetPos(PKD pkd,PARTICLE *p) {
    return _mm256_mul_pd(_mm256_cvtepi32_pd(pkdGetPosRaw(pkd,p)),_mm256_set1_pd(1.0/INTEGER_FACTOR));
}
#endif

static inline vel_t *pkdVel( PKD pkd, PARTICLE *p ) {
    return CAST(vel_t *,pkdField(p,pkd->oFieldOffset[oVelocity]));
}
static inline float *pkdAccel( PKD pkd, PARTICLE *p ) {
    return CAST(float *,pkdField(p,pkd->oFieldOffset[oAcceleration]));
}
static inline const float *pkdAccelRO( PKD pkd, const PARTICLE *p ) {
    return CAST(const float *,pkdFieldRO(p,pkd->oFieldOffset[oAcceleration]));
}
static inline float *pkdPot( PKD pkd, PARTICLE *p ) {
    return pkd->oFieldOffset[oPotential] ? CAST(float *,pkdField(p,pkd->oFieldOffset[oPotential])) : NULL;
}
static inline uint16_t *pkdRungDest( PKD pkd, PARTICLE *p ) {
    return CAST(uint16_t *,pkdField(p,pkd->oFieldOffset[oRungDest]));
}
static inline uint64_t *pkdParticleID( PKD pkd, PARTICLE *p ) {
    return CAST(uint64_t *,pkdField(p,pkd->oFieldOffset[oParticleID]));
}
/* Sph variables */
static inline SPHFIELDS *pkdSph( PKD pkd, PARTICLE *p ) {
#if defined(OPTIM_UNION_EXTRAFIELDS) && defined(DEBUG_UNION_EXTRAFIELDS)
    assert( pkdSpecies(pkd,p)==FIO_SPECIES_SPH);
#endif
    return ((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]));
}
/* NewSph variables */
static inline NEWSPHFIELDS *pkdNewSph( PKD pkd, PARTICLE *p ) {
    return ((NEWSPHFIELDS *) pkdField(p,pkd->oFieldOffset[oNewSph]));
}
static inline const SPHFIELDS *pkdSphRO( PKD pkd, const PARTICLE *p ) {
#if defined(OPTIM_UNION_EXTRAFIELDS) && defined(DEBUG_UNION_EXTRAFIELDS)
    assert( pkdSpecies(pkd,p)==FIO_SPECIES_SPH);
#endif
    return ((const SPHFIELDS *) pkdFieldRO(p,pkd->oFieldOffset[oSph]));
}
static inline const NEWSPHFIELDS *pkdNewSphRO( PKD pkd, const PARTICLE *p ) {
    return ((const NEWSPHFIELDS *) pkdFieldRO(p,pkd->oFieldOffset[oNewSph]));
}
static inline STARFIELDS *pkdStar( PKD pkd, PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
    assert( pkdSpecies(pkd,p)==FIO_SPECIES_STAR);
#endif //DEBUG
    return ((STARFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]));
#else
    return ((STARFIELDS *) pkdField(p,pkd->oFieldOffset[oStar]));
#endif
}
static inline const STARFIELDS *pkdStarRO( PKD pkd, PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
    assert( pkdSpecies(pkd,p)==FIO_SPECIES_STAR);
#endif //DEBUG
    return ((const STARFIELDS *) pkdFieldRO(p,pkd->oFieldOffset[oSph]));
#else
    return ((const STARFIELDS *) pkdFieldRO(p,pkd->oFieldOffset[oStar]));
#endif
}
static inline BHFIELDS *pkdBH( PKD pkd, PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
    assert( pkdSpecies(pkd,p)==FIO_SPECIES_BH);
#endif //DEBUG
    return ((BHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]));
#else
    return ((BHFIELDS *) pkdField(p,pkd->oFieldOffset[oBH]));
#endif
}
static inline const BHFIELDS *pkdBHRO( PKD pkd, PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
    assert( pkdSpecies(pkd,p)==FIO_SPECIES_BH);
#endif //DEBUG
    return ((const BHFIELDS *) pkdFieldRO(p,pkd->oFieldOffset[oSph]));
#else
    return ((const BHFIELDS *) pkdFieldRO(p,pkd->oFieldOffset[oBH]));
#endif
}
static inline double *pkd_vPred( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->vPred[0]);
}
#ifndef OPTIM_REMOVE_UNUSED
static inline float *pkd_u( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->u);
}
static inline float *pkd_uPred( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->uPred);
}
static inline float *pkd_uDot( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->uDot);
}
static inline float *pkd_c( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->c);
}
static inline float *pkd_divv( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->divv);
}
static inline float *pkd_fMetals( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->fMetals);
}
static inline float *pkd_diff( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->diff);
}
static inline float *pkd_fMetalsDot( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->fMetalsDot);
}
static inline float *pkd_fMetalsPred( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->fMetalsPred);
}
#endif //OPTIM_REMOVE_UNUSED
static inline char **pkd_pNeighborList( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oFieldOffset[oSph]))->pNeighborList);
}

static inline float *pkd_Timer( PKD pkd, PARTICLE *p ) {
    return &(((STARFIELDS *) pkdField(p,pkd->oFieldOffset[oStar]))->fTimer);
}

static inline int pkdIsDeleted(PKD pkd,PARTICLE *p) {
    return (pkdSpecies(pkd,p) == FIO_SPECIES_LAST);
}

static inline int pkdIsNew(PKD pkd,PARTICLE *p) {
    return (p->iOrder == IORDERMAX);
}


/*
** From tree.c:
*/
void pkdVATreeBuild(PKD pkd,int nBucket);
void pkdTreeBuild(PKD pkd,int nBucket,int nGroup,uint32_t uRoot,uint32_t uTemp,double ddHonHLimit);
uint32_t pkdDistribTopTree(PKD pkd, uint32_t uRoot, uint32_t nTop, KDN *pTop, int allocateMemory);
void pkdOpenCloseCaches(PKD pkd,int bOpen,int bFixed);
void pkdTreeInitMarked(PKD pkd);
void pkdDumpTrees(PKD pkd,int bOnlyVA,uint8_t uRungDD);
void pkdCombineCells1(PKD,KDN *pkdn,KDN *p1,KDN *p2);
void pkdCombineCells2(PKD,KDN *pkdn,KDN *p1,KDN *p2);
void pkdCalcRoot(PKD,uint32_t,double *,MOMC *);
void pkdDistribRoot(PKD,double *,MOMC *);
void pkdGroupOrder(PKD pkd,uint32_t *iGrpOffset);
void pkdTreeBuildByGroup(PKD pkd, int nBucket, int nGroup);

/*
** From pkd.c:
*/
size_t pkdClCount(PKD pkd);
size_t pkdClMemory(PKD pkd);
size_t pkdIllMemory(PKD pkd);
size_t pkdIlcMemory(PKD pkd);
size_t pkdIlpMemory(PKD pkd);
void pkdReadFIO(PKD pkd,FIO fio,uint64_t iFirst,int nLocal,double dvFac, double dTuFac);
void pkdSetSoft(PKD pkd,double dSoft);
void pkdSetSmooth(PKD pkd,double dSmooth);
void pkdSetCrit(PKD pkd,double dCrit);
void pkdCalcBound(PKD,BND *);
void pkdEnforcePeriodic(PKD,BND *);
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul);
int pkdWeight(PKD,int,double,int,int,int,int *,int *,double *,double *);
void pkdCountVA(PKD,int,double,int *,int *);
double pkdTotalMass(PKD pkd);
uint8_t pkdGetMinDt(PKD pkd);
void   pkdSetGlobalDt(PKD pkd, uint8_t minDt);
int pkdLowerPart(PKD,int,double,int,int);
int pkdUpperPart(PKD,int,double,int,int);
int pkdWeightWrap(PKD,int,double,double,int,int,int,int *,int *);
int pkdLowerPartWrap(PKD,int,double,double,int,int);
int pkdUpperPartWrap(PKD,int,double,double,int,int);
int pkdLowerOrdPart(PKD,uint64_t,int,int);
int pkdUpperOrdPart(PKD,uint64_t,int,int);
int pkdActiveOrder(PKD);

void pkdOrbBegin(PKD pkd, int nRungs);
int pkdOrbSelectRung(PKD pkd, int iRung);
void pkdOrbUpdateRung(PKD pkd);
void pkdOrbFinish(PKD pkd);
void pkdOrbSplit(PKD pkd,int iDomain);
int pkdOrbRootFind(
    PKD pkd,double dFraction,uint64_t nLowerMax, uint64_t nUpperMax,
    double dReserveFraction, BND *bnd, double *dSplitOut, int *iDim);
void pkdOrbUpdateRung(PKD pkd);
/*#define PEANO_HILBERT_KEY_MAX 0x3ffffffffffull*/ /* 2d */
#define PEANO_HILBERT_KEY_MAX 0x7fffffffffffffffull /* 3d */
void pkdPeanoHilbertDecomp(PKD pkd, int nRungs, int iMethod);
void pkdRungOrder(PKD pkd, int iRung, total_t *nMoved);
int pkdColRejects(PKD,int);
int pkdColRejects_Old(PKD,int,double,double,int);

int pkdSwapRejects(PKD,int);
int pkdSwapSpace(PKD);
int pkdActive(PKD);
int pkdInactive(PKD);

int pkdColOrdRejects(PKD,uint64_t,int);
void pkdLocalOrder(PKD,uint64_t iMinOrder,uint64_t iMaxOrder);
void pkdCheckpoint(PKD pkd,const char *fname);
void pkdRestore(PKD pkd,const char *fname);
void pkdWriteHeaderFIO(PKD pkd, FIO fio, double dScaleFactor, double dTime,
                       uint64_t nDark, uint64_t nGas, uint64_t nStar, uint64_t nBH,
                       double dBoxSize, double h, int nProcessors, UNITS units);
uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac,double dTuFac,BND *bnd);
void pkdWriteFromNode(PKD pkd,int iNode, FIO fio,double dvFac,double dTuFac,BND *bnd);
void pkdWriteViaNode(PKD pkd, int iNode);
char *pkdPackArray(PKD pkd,int iSize,void *vBuff,int *piIndex,int n,int field,int iUnitSize,double dvFac,int bMarked);
void pkdSendArray(PKD pkd, int iNode, int field, int iUnitSize,double dvFac,int bMarked);
void *pkdRecvArray(PKD pkd,int iNode, void *pDest, int iUnitSize);
void pkdGravAll(PKD pkd,
                struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                double dTime,int nReps,int bPeriodic,
                int bEwald,int nGroup,int iRoot1, int iRoot2,
                double fEwCut,double fEwhCut,double dThetaMin,SPHOptions *SPHoptions,
                uint64_t *pnActive,
                double *pdPart,double *pdPartNumAccess,double *pdPartMissRatio,
                double *pdCell,double *pdCellNumAccess,double *pdCellMissRatio,
                double *pdFlop,uint64_t *pnRung);
void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double *L,double *F,double *W);
void pkdProcessLightCone(PKD pkd,PARTICLE *p,float fPot,double dLookbackFac,double dLookbackFacLCP,
                         double dDriftDelta,double dKickDelta,double dBoxSize,int bLightConeParticles);
void pkdGravEvalPP(const PINFOIN &Part, ilpTile &tile, PINFOOUT &Out );
void pkdDensityEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions);
void pkdSPHForcesEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions);
void pkdGravEvalPC(const PINFOIN &pPart, ilcTile &tile, PINFOOUT &pOut );
void pkdDrift(PKD pkd,int iRoot,double dTime,double dDelta,double,double,int bDoGas);
void pkdEndTimestepIntegration(PKD pkd, struct inEndTimestep in);
#ifdef OPTIM_REORDER_IN_NODES
    void pkdReorderWithinNodes(PKD pkd);
#endif
void pkdKickKDKOpen(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void pkdKickKDKClose(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void pkdKick(PKD pkd,double dTime,double dDelta,int bDoGas,double,double,double,uint8_t uRungLo,uint8_t uRungHi);
void pkdKickTree(PKD pkd,double dTime,double dDelta,double,double,double,int iRoot);
void pkdSwapAll(PKD pkd, int idSwap);
void pkdInitCosmology(PKD pkd,struct csmVariables *cosmo);
void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdActiveRung(PKD pkd, int iRung, int bGreater);
void pkdCountRungs(PKD pkd,uint64_t *nRungs);
void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
                  double dDelta, int iMaxRung,
                  double dEta,double dVelFac,double dAccFac,
                  int bDoGravity,int bEpsAcc,double dhMinOverSoft);
void pkdCooling(PKD pkd,double,double,int,int,int,int);
void pkdChemCompInit(PKD pkd, struct inChemCompInit in);
void pkdSphStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
                double dDelta, int iMaxRung,double dEta, double dAccFac, double dEtaUDot);
#define CORRECTENERGY_IN 1
#define CORRECTENERGY_OUT 2
#define CORRECTENERGY_SPECIAL 3
#define CORRECTENERGY_IC_MESHLESS 4
void pkdCorrectEnergy(PKD pkd, double dTuFac, double z, double dTime, int iType );
void pkdDensityStep(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, int iMaxRung, double dDelta, double dEta, double dRhoFac);
int pkdUpdateRung(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
                  uint8_t uRung,int iMaxRung,uint64_t *nRungCount);
uint8_t pkdDtToRung(double dT, double dDelta, uint8_t uMaxRung);
static inline uint8_t pkdDtToRungInverse(float fT, float fiDelta, uint8_t uMaxRung) {
    union {
        float f;
        struct {
            uint32_t mantisa : 23;
            uint32_t exponent : 8;
            uint32_t sign : 1;
        } ieee;
    } T;
    int iRung;
    T.f = fabsf(fiDelta)*fT;
    if (T.f>=1.0) return 0;
    iRung = 126 - T.ieee.exponent; /* -log2(d) */
    if (iRung > uMaxRung) return uMaxRung;
    else return iRung;
}
int pkdOrdWeight(PKD pkd,uint64_t iOrdSplit,int iSplitSide,int iFrom,int iTo,
                 int *pnLow,int *pnHigh);
void pkdDeleteParticle(PKD pkd, PARTICLE *p);
void pkdNewParticle(PKD pkd, PARTICLE *p);
int pkdResetTouchRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdIsGas(PKD,PARTICLE *);
int pkdIsDark(PKD,PARTICLE *);
int pkdIsStar(PKD,PARTICLE *);
int pkdIsBH(PKD,PARTICLE *);
void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
                  int *nDeltaStar);
void pkdNewOrder(PKD pkd, int nStart);

struct outGetNParts {
    total_t n;
    total_t nGas;
    total_t nDark;
    total_t nStar;
    total_t nBH;
    total_t nMaxOrder;
};

void pkdMoveDeletedParticles(PKD pkd, total_t *n, total_t *nGas, total_t *nDark, total_t *nStar, total_t *nBH);
void pkdGetNParts(PKD pkd, struct outGetNParts *out );
void pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar, int nBH);
void pkdInitRelaxation(PKD pkd);

#ifdef USE_GRAFIC
void pkdGenerateIC(PKD pkd, GRAFICCTX gctx, int iDim,
                   double fSoft, double fMass, int bCannonical);
#endif
int pkdGetClasses( PKD pkd, int nMax, PARTCLASS *pClass );
void pkdSetClasses( PKD pkd, int n, PARTCLASS *pClass, int bUpdate );
void pkdSetClass( PKD pkd, float fMass, float fSoft, FIO_SPECIES eSpecies, PARTICLE *p );

int pkdCountSelected(PKD pkd);
int pkdSelSpecies(PKD pkd,uint64_t mSpecies, int setIfTrue, int clearIfFalse);
int pkdSelGroup(PKD pkd, int iGroup, int setIfTrue, int clearIfFalse);
int pkdSelActive(PKD pkd, int setIfTrue, int clearIfFalse);
int pkdSelBlackholes(PKD pkd, int setIfTrue, int clearIfFalse);
int pkdSelMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse );
int pkdSelById(PKD pkd,uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse );
int pkdSelPhaseDensity(PKD pkd,double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse );
int pkdSelBox(PKD pkd,double *dCenter, double *dSize, int setIfTrue, int clearIfFalse );
int pkdSelSphere(PKD pkd,double *r, double dRadius, int setIfTrue, int clearIfFalse );
int pkdSelCylinder(PKD pkd,double *dP1, double *dP2, double dRadius, int setIfTrue, int clearIfFalse );

int pkdDeepestPot(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
                  double *r, float *fPot);
void pkdProfile(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
                const double *dCenter, const double *dRadii, int nBins,
                const double *com, const double *vcm, const double *L);
void pkdCalcDistance(PKD pkd, double *dCenter, int bPeriodic);
uint_fast32_t pkdCountDistance(PKD pkd, double r2i, double r2o );
void pkdCalcCOM(PKD pkd, double *dCenter, double dRadius, int bPeriodic,
                double *com, double *vcm, double *L,
                double *M, uint64_t *N);
void pkdCalcMtot(PKD pkd, double *M, uint64_t *N);
void pkdTreeUpdateFlagBounds(PKD pkd,uint32_t uRoot,SPHOptions *SPHoptions);
#ifdef MDL_FFTW
void pkdAssignMass(PKD pkd, uint32_t iLocalRoot, int iAssignment, int iGrid, float dDelta);
void pkdInterlace(PKD pkd, int iGridTarget, int iGridSource);
float getLinAcc(PKD pkd, MDLFFT fft,int cid, double r[3]);
void pkdSetLinGrid(PKD pkd,double a0, double a, double a1, double dBSize, int nGrid, int iSeed,
                   int bFixed, float fPhase);
void pkdMeasureLinPk(PKD pkd, int nGrid, double dA, double dBoxSize,
                     int nBins,  int iSeed, int bFixed, float fPhase,
                     double *fK, double *fPower, uint64_t *nPower);
#endif
void pkdOutPsGroup(PKD pkd,char *pszFileName,int iType);

void pkdLightConeOpen(PKD pkd, const char *fname,int nSideHealpix);
void pkdLightConeClose(PKD pkd, const char *healpixname);
void pkdLightCone(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
                  double dLookbackFac,double dLookbackFacLCP,
                  double *dtLCDrift,double *dtLCKick);
void pkdLightConeVel(PKD pkd,double dBoxSize);

struct outGetParticles { /* Array of these */
    uint64_t id;
    float    mass, phi;
    float    r[3], v[3];
};
int pkdGetParticles(PKD pkd, int nIn, uint64_t *ID, struct outGetParticles *out);

#ifdef USE_CUDA
    extern int CUDAinitWorkPP( void *vpp, void *vwork );
    extern int CUDAcheckWorkPP( void *vpp, void *vwork );
    extern int CUDAinitWorkPC( void *vpp, void *vwork );
    extern int CUDAcheckWorkPC( void *vpp, void *vwork );
    extern int CUDAinitWorkEwald( void *vpp, void *vwork );
    extern int CUDAcheckWorkEwald( void *vpp, void *vwork );
    void *pkdCudaClientInitialize(PKD pkd);
#endif
#ifdef USE_CL
    extern int CLinitWorkEwald( void *vcl, void *ve, void *vwork );
    int CLcheckWorkEwald( void *ve, void *vwork );
    extern void clEwaldInit(void *cudaCtx, struct EwaldVariables *ewIn, EwaldTable *ewt );
#endif

#define vec_sub(r,a,b) do {\
    int i;\
    for (i=0; i<3; i++) (r)[i] = (a)[i] - (b)[i];   \
} while(0)

#define vec_add_const_mult(r,a,c,b) do {\
    int i;\
    for (i=0; i<3; i++) (r)[i] = (a)[i] + (c) * (b)[i]; \
} while(0)

#define matrix_vector_mult(b,mat,a) do {\
    int i;\
    for (i=0; i<3; i++) {\
        int j;\
    (b)[i] = 0.0;                   \
        for (j=0; j<3; j++) (b)[i] += (mat)[i][j] * (a)[j]; \
    }\
} while(0)

static inline double dot_product(const double *a,const double *b) {
    double r = 0.0;
    int i;
    for (i=0; i<3; i++) r += a[i]*b[i];
    return r;
}

#define cross_product(r,a,b) do {\
    (r)[0] = (a)[1] * (b)[2] - (a)[2] * (b)[1] ;    \
    (r)[1] = (a)[2] * (b)[0] - (a)[0] * (b)[2] ;    \
    (r)[2] = (a)[0] * (b)[1] - (a)[1] * (b)[0] ;    \
} while(0)

#define mat_transpose(mat,trans_mat) do {\
    int i;               \
    for (i=0; i<3; i++) {           \
    int j;                  \
        for (j=0; j<3; j++) {           \
            (trans_mat)[i][j] = (mat)[j][i];    \
        }                   \
    }                   \
} while(0)


#endif
