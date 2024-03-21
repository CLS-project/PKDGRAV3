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
#include <vector>

#include "mdl.h"
#ifdef USE_CUDA
    #include "cuda/cudautil.h"
#endif
#ifdef USE_METAL
    #include "metal/metal.h"
#endif
#include "gravity/ilp.h"
#include "gravity/ilc.h"
#include "gravity/cl.h"
#include "gravity/moments.h"
#include "cosmo.h"
#include "units.h"
#include "io/fio.h"
#include "basetype.h"
#include "core/treenode.h"
#include "io/iomodule.h"
#include "SPH/SPHOptions.h"
#ifdef HAVE_EOSLIB_H
    #include <EOSlib.h>
#endif
#include "core/bound.h"
#include "eEOS/eEOS_struct.h"
#include "chemistry.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef GRACKLE
#include <grackle.h>
#endif

#ifdef __cplusplus
}
#endif

typedef uint_fast32_t local_t; /* Count of particles locally (per processor) */
typedef uint_fast64_t total_t; /* Count of particles globally (total number) */

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

#define MAX_TIMERS      10

/*
** Memory models.  Each is a bit mask that indicates that additional fields should be
** added to the particle structure.
*/
#define PKD_MODEL_VELOCITY     (1<<0)  /* Velocity Required */
#define PKD_MODEL_ACCELERATION (1<<1)  /* Acceleration Required */
#define PKD_MODEL_POTENTIAL    (1<<2)  /* Potential Required */
#define PKD_MODEL_GROUPS       (1<<3)  /* Group profiling */
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
#define PKD_MODEL_GLOBALGID    (1<<17) /* Global group identifier per particle */

#define PKD_MODEL_NODE_MOMENT  (1<<24) /* Include moment in the tree */
#define PKD_MODEL_NODE_ACCEL   (1<<25) /* mean accel on cell (for grav step) */
#define PKD_MODEL_NODE_VEL     (1<<26) /* center of mass velocity for cells */
#define PKD_MODEL_NODE_SPHBNDS (1<<27) /* Include 3 extra bounds in tree */

#define PKD_MODEL_NODE_BND     (1<<28) /* Include normal bounds in tree */
#define PKD_MODEL_NODE_VBND    (1<<29) /* Include velocity bounds in tree for phase-space density*/
#define PKD_MODEL_NODE_BOB     (1<<30) /* Ball of Balls, or Box of Balls */

#define PKD_MODEL_NEW_SPH      (1ULL<<32) /* New Sph Fields */

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

typedef struct partLightCone {
    uint64_t id;
    float pos[3];
    float vel[3];
    float pot;
} LIGHTCONEP;

#define NMAX_OPENCALC   1000

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

#if defined(USE_SIMD) && !defined(__CUDACC__)
typedef struct {
    struct PMOMC {
        dvec m;
        dvec xx,yy,xy,xz,yz;
        dvec xxx,xyy,xxy,yyy,xxz,yyz,xyz;
        dvec xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
        dvec zz;
        dvec xzz,yzz,zzz;
        dvec xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
    } ewm;
    struct PEWALDVARS {
        dvec fEwCut2,fInner2,alpha,alpha2,ialpha,k1,ka;
        dvec Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
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
    blitz::TinyVector<double,3> rref;
    blitz::TinyVector<double,3> ravg;
    blitz::TinyVector<double,3> rmin;
    blitz::TinyVector<double,3> rmax;
    blitz::TinyVector<double,3> rcom;
    blitz::TinyVector<double,3> vcom;
} HopGroupTable;

typedef struct {
    remoteID key;
    remoteID name;
    uint32_t iLink;
} FOFRemote;

typedef struct {
    blitz::TinyVector<float,3> rPot;
    float minPot;
    blitz::TinyVector<float,3> rcen;
    blitz::TinyVector<float,3> rcom;
    blitz::TinyVector<float,3> vcom;
    blitz::TinyVector<float,3> angular;
    blitz::TinyVector<float,6> inertia;
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
    uint64_t iGlobalGid;
} TinyGroupTable;

/* IA: For the BH seeding, we need a small part of the FoF information
 *  to survive for a while
 */
typedef struct {
    blitz::TinyVector<double,3>  rPot;
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
    blitz::TinyVector<double,3> L;
    double vel_radial;
    double vel_radial_sigma;
    double vel_tang_sigma;
    uint64_t nParticles;
} PROFILEBIN;

typedef struct shapesBin {
    blitz::TinyVector<double,3> com;
    double dMassEnclosed;
    blitz::TinyVector<blitz::TinyVector<double,3>,3> dInertia;
    blitz::TinyVector<blitz::TinyVector<double,3>,3> ell_matrix;
    blitz::TinyVector<double,3> ell_center;
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

class pkdContext {
public:
    explicit pkdContext(
        mdl::mdlClass *mdl,int nStore,uint64_t nMinTotalStore,uint64_t nMinEphemeral,uint32_t nEphemeralBytes,
        int nTreeBitsLo, int nTreeBitsHi,
        int iCacheSize,int iCacheMaxInflight,int iWorkQueueSize,const blitz::TinyVector<double,3> &fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,uint64_t nBH,
        uint64_t mMemoryModel);
    virtual ~pkdContext();

protected:  // Support for memory models
    int NodeAddStruct    (int n);   // Add a NODE structure: assume double alignment
    int NodeAddDouble    (int n=1); // Add n doubles to the node structure
    int NodeAddFloat     (int n=1); // Add n floats to the node structure
    int NodeAddInt64     (int n=1); // Add n 64-bit integers to the node structure
    int NodeAddInt32     (int n=1); // Add n 32-bit integers to the node structure
protected:
    PARTICLE *pTempPRIVATE = nullptr;
    uint32_t nEphemeralBytes = 0; /* per-particle */
    void *storageBase;
    uint64_t storageSize;
public:
    particleStore particles;
    treeStore tree;
    int FreeStore() { return particles.FreeStore(); }
    int Local() { return particles.Local(); }
    int SetLocal(int n);
    int AddLocal(int n);
    auto EphemeralBytes() {return nEphemeralBytes; }
    static constexpr auto MaxNodeSize() { return sizeof(KDN) + 2*sizeof(Bound) + sizeof(FMOMR) + 6*sizeof(double) + sizeof(SPHBNDS); }
    auto NodeSize() { return tree.ElementSize(); }
    auto ParticleMemory() { return (particles.ParticleSize() + EphemeralBytes()) * (FreeStore()+1); }
    PARTICLE *Particle(int i) { return particles.Element(i); }

    void SaveParticle(PARTICLE *a) {
        memcpy(pTempPRIVATE,a,particles.ParticleSize());
    }
    void LoadParticle(PARTICLE *a) {
        memcpy(a,pTempPRIVATE,particles.ParticleSize());
    }
    void CopyParticle(PARTICLE *a, PARTICLE *b) {
        memcpy(a,b,particles.ParticleSize());
    }

public:
    auto Nodes() const {
        return tree.Nodes();
    }
    auto Node(KDN *pBase,int iNode) {
        return tree.Element(pBase,iNode);
    }
    auto TreeNode(int iNode) {
        return tree[iNode];
    }
    auto TreeAlignNode() {
        return tree.AlignNode();
    }
    auto TreeAllocNode(int n=1) {
        return tree.AllocNode(n);
    }
    void TreeAllocNodePair(int *iLeft, int *iRight) {
        *iLeft = TreeAllocNode();
        *iRight = TreeAllocNode();
    }
    auto TreeAllocRootNode() {
        if ((Nodes()&1)==0) TreeAllocNode();
        return TreeAllocNode();
    }
// I/O
public:
    void Restore(uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd);

// Rockstar Analysis
protected:
    uint64_t nRsElements;
public:
    void RsHaloIdStart(uint64_t nElements,bool bAppend=false);
    void RsHaloIdFinish(uint64_t nElements);
    void RsHaloIdRead(uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd);
    void RsIdStart(uint64_t nElements,bool bAppend=false);
    void RsIdFinish(uint64_t nElements);
    void RsIdRead(uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd);
    void RsIdSave(int iGroup,const std::string &filename,int iSegment,int nSegment);
    void RsReorder(uint64_t *pOrds);
    void RsExtract(const std::string &filename, int iSegment, uint64_t *pOrds);

public:
    mdl::mdlClass *mdl;
    auto Self()    const {
        return mdl->Self();
    }
    auto Threads() const {
        return mdl->Threads();
    }

public:
    int nRejects;
    int nActive;
    blitz::TinyVector<uint64_t,MAX_RUNG+1> nRung;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nBH;
    blitz::TinyVector<double,3> fPeriod;
    int iTopTree[NRESERVED_NODES];
    int nNodesFull;     /* number of nodes in the full tree (including very active particles) */
    Bound bnd;
    Bound vbnd;
    /*
    ** Light cone variables.
    */
    int nLayerMax;
    int *nBoxLC;
    double *lcOffset0;
    double *lcOffset1;
    double *lcOffset2;
    asyncFileInfo afiLightCone;
    LIGHTCONEP *pLightCone;
    int nLightCone, nLightConeMax;
    int64_t nHealpixPerDomain;
    int64_t nSideHealpix;
    healpixData *pHealpixData;
    gsl_spline *interp_scale; // interpolation table for 1/a given r in the lightcone
    gsl_interp_accel *interp_accel = gsl_interp_accel_alloc();

    void *pLite;
    /*
    ** Advanced memory models
    */
    int bNoParticleOrder;
    int bIntegerPosition;

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
    int nWpPending;
    uint64_t nTilesTotal, nTilesCPU;
    /*
    ** Opening angle table for mass weighting.
    */
    float fiCritTheta;
    /* Potential Energy for when potential is not in the particle */
    double dEnergyU;
    /* also put kinetic energy here to calculate it on the fly */
    double dEnergyT;
    double dEnergyW;
    blitz::TinyVector<double,3> dEnergyF, dEnergyL;

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
    Bound bndInterior;  /* this gets calculated at the start of Fof for now but should be done elsewhere */
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
#ifdef USE_METAL
    MetalClient *metalClient;
#endif
#ifdef MDL_FFTW
    MDLFFT fft;
#endif

    SPHOptions SPHoptions;
#ifdef HAVE_EOSLIB_H
    EOSmaterial *materials[EOS_N_MATERIAL_MAX] = {NULL};
#endif
public:
    auto &NodeBOB(KDN *n) {
        return tree.BOB(n);
    }
    auto &NodeBOB(const KDN *n) {
        return tree.BOB(n);
    }
public:
    void ActiveRung(int iRung, bool bGreater) {
        uMinRungActive = iRung;
        uMaxRungActive = bGreater ? 255 : iRung;
        particles.ActiveRung(uMinRungActive,uMaxRungActive);
    }
    template<typename T>
    auto Ephemeral() {
        return static_cast<T *>(pLite);
    }
};
using PKD = pkdContext *;

static inline int pkdIsRungActive(PKD pkd, uint8_t uRung ) {
    return uRung >= pkd->uMinRungActive && uRung <= pkd->uMaxRungActive;
}
static inline int pkdIsActive(PKD pkd, PARTICLE *p ) {
    return pkdIsRungActive(pkd,p->uRung);
}

void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize);

#if defined(FEEDBACK) || defined(BLACKHOLES)
static inline void pkdAddFBEnergy(PKD pkd, particleStore::ParticleReference &p, meshless::FIELDS *psph, double dConstGamma) {
#ifndef OLD_FB_SCHEME
    psph->Uint += psph->fAccFBEnergy;
    psph->E += psph->fAccFBEnergy;
#ifdef ENTROPY_SWITCH
    psph->S += psph->fAccFBEnergy*(dConstGamma-1.) *
               pow(p.density(), -dConstGamma+1);
#endif
    psph->fAccFBEnergy = 0.0;
#endif //OLD_FB_SCHEME
}
#endif

/*
** The size of a particle is variable based on the memory model.
** The following three routines must be used instead of accessing pStore
** directly.  pkdParticle will return a pointer to the i'th particle.
** The Size and Base functions are intended for cache routines; no other
** code should care about sizes of the particle structure.
*/

static inline int32_t pkdGetGroup( PKD pkd, const PARTICLE *p ) {
    return pkd->particles.group(p);
}
static inline int64_t pkdGetGlobalGid( PKD pkd, const PARTICLE *p ) {
    return pkd->particles.global_gid(p);
}
static inline void pkdSetGroup( PKD pkd, PARTICLE *p, uint32_t gid ) {
    pkd->particles.set_group(p,gid);
}
static inline void pkdSetGlobalGid( PKD pkd, PARTICLE *p, uint64_t gid ) {
    pkd->particles.global_gid(p) = gid;
}

static inline float pkdDensity( PKD pkd, const PARTICLE *p ) {
    return pkd->particles.density(p);
}
static inline void pkdSetDensity( PKD pkd, PARTICLE *p, float fDensity ) {
    pkd->particles.set_density(p,fDensity);
}

/* Here is the new way of getting mass and softening */
static inline float pkdMass( PKD pkd, PARTICLE *p ) {
    return pkd->particles.mass(p);
}
static inline float pkdSoft0( PKD pkd, PARTICLE *p ) {
    return pkd->particles.soft0(p);
}
static inline float pkdSoft( PKD pkd, PARTICLE *p ) {
    return pkd->particles.soft(p);
}
static inline FIO_SPECIES pkdSpecies( PKD pkd, PARTICLE *p ) {
    return pkd->particles.species(p);
}
static inline int pkdiMat( PKD pkd, PARTICLE *p ) {
    return pkd->particles.iMat(p);
}

static inline double pkdPos(PKD pkd,PARTICLE *p,int d) {
    if (pkd->bIntegerPosition) return pkdIntPosToDbl(pkd,pkd->particles.get<int32_t[3]>(p,PKD_FIELD::oPosition)[d]);
    else return pkd->particles.get<double[3]>(p,PKD_FIELD::oPosition)[d];
}
static inline void pkdSetPos(PKD pkd,PARTICLE *p,int d,double v) {
    if (pkd->bIntegerPosition) pkd->particles.get<int32_t[3]>(p,PKD_FIELD::oPosition)[d] = pkdDblToIntPos(pkd,v);
    else pkd->particles.get<double[3]>(p,PKD_FIELD::oPosition)[d] = v;
}
#define pkdGetPos3(pkd,p,d1,d2,d3) ((d1)=pkdPos(pkd,p,0),(d2)=pkdPos(pkd,p,1),(d3)=pkdPos(pkd,p,2))
#define pkdGetPos1(pkd,p,d) pkdGetPos3(pkd,p,(d)[0],(d)[1],(d)[2])

#if defined(__AVX__) && defined(USE_SIMD)
static inline __m128i pkdGetPosRaw(PKD pkd,PARTICLE *p) {
    return _mm_loadu_si128(&pkd->particles.get<__m128i>(p,PKD_FIELD::oPosition));
}
static inline __m256d pkdGetPos(PKD pkd,PARTICLE *p) {
    return _mm256_mul_pd(_mm256_cvtepi32_pd(pkdGetPosRaw(pkd,p)),_mm256_set1_pd(1.0/INTEGER_FACTOR));
}
#endif

static inline int pkdIsDeleted(PKD pkd,PARTICLE *p) {
    return (pkdSpecies(pkd,p) == FIO_SPECIES_UNKNOWN);
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
void pkdCombineCells1(PKD,treeStore::NodePointer pkdn,treeStore::NodePointer p1,treeStore::NodePointer p2);
void pkdCombineCells2(PKD,treeStore::NodePointer pkdn,treeStore::NodePointer p1,treeStore::NodePointer p2);
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
void pkdEnforcePeriodic(PKD,Bound);
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul);
int pkdWeight(PKD,int,double,int,int,int,int *,int *,double *,double *);
void pkdCountVA(PKD,int,double,int *,int *);
double pkdTotalMass(PKD pkd);
uint8_t pkdGetMinDt(PKD pkd);
void pkdSetGlobalDt(PKD pkd, uint8_t minDt);
int pkdLowerPart(PKD,int,double,int,int);
int pkdUpperPart(PKD,int,double,int,int);
int pkdWeightWrap(PKD,int,double,double,int,int,int,int *,int *);
int pkdLowerPartWrap(PKD,int,double,double,int,int);
int pkdUpperPartWrap(PKD,int,double,double,int,int);
int pkdLowerOrdPart(PKD,uint64_t,int,int);
int pkdUpperOrdPart(PKD,uint64_t,int,int);
int pkdActiveOrder(PKD);

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
void pkdWriteHeaderFIO(PKD pkd, FIO fio, double dScaleFactor, double dTime,
                       uint64_t nDark, uint64_t nGas, uint64_t nStar, uint64_t nBH,
                       double dBoxSize, double h, int nProcessors, UNITS units);
uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac,double dTuFac,Bound bnd);
void pkdWriteFromNode(PKD pkd,int iNode, FIO fio,double dvFac,double dTuFac,Bound bnd);
void pkdWriteViaNode(PKD pkd, int iNode);
char *pkdPackArray(PKD pkd,int iSize,void *vBuff,int *piIndex,int n,PKD_FIELD field,int iUnitSize,double dvFac,int bMarked);
void pkdSendArray(PKD pkd, int iNode, PKD_FIELD field, int iUnitSize,double dvFac,int bMarked);
void *pkdRecvArray(PKD pkd,int iNode, void *pDest, int iUnitSize);
void pkdGravAll(PKD pkd,
                struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                double dTime,int nReps,int bPeriodic,int bGPU,
                int bEwald,int iRoot1, int iRoot2,
                double fEwCut,double fEwhCut,double dThetaMin,SPHOptions *SPHoptions,
                uint64_t *pnActive,
                double *pdPart,double *pdPartNumAccess,double *pdPartMissRatio,
                double *pdCell,double *pdCellNumAccess,double *pdCellMissRatio,
                double *pdFlop,uint64_t *pnRung);
void pkdCalcEandL(PKD pkd,double &T,double &U,double &Eth,blitz::TinyVector<double,3> &L,blitz::TinyVector<double,3> &F,double &W);
void pkdProcessLightCone(PKD pkd,PARTICLE *p,float fPot,double dLookbackFac,double dLookbackFacLCP,
                         double dDriftDelta,double dKickDelta,double dBoxSize,int bLightConeParticles,
                         blitz::TinyVector<double,3> hlcp,double tanalpha2);
void pkdGravEvalPP(const PINFOIN &Part, ilpTile &tile, PINFOOUT &Out );
void pkdDensityEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions);
void pkdDensityCorrectionEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions);
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
void pkdInitCosmology(PKD pkd, struct csmVariables *cosmo);
void pkdInitLightcone(PKD pkd,int bBowtie,int bLightConeParticles,double dBoxSize,double dRedshiftLCP,double alphaLCP,blitz::TinyVector<double,3> hLCP);
void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
                  double dDelta, int iMaxRung,
                  double dEta,double dVelFac,double dAccFac,
                  int bDoGravity,int bEpsAcc,double dhMinOverSoft);
void pkdCooling(PKD pkd,double,double,int,int,int,int);
void pkdChemCompInit(PKD pkd, struct inChemCompInit in);
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
void pkdDeleteParticle(PKD pkd, particleStore::ParticleReference &p);
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

int pkdDeepestPot(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
                  double *r, float *fPot);
void pkdProfile(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
                const blitz::TinyVector<double,3> &dCenter, const double *dRadii, int nBins,
                const blitz::TinyVector<double,3> &com, const blitz::TinyVector<double,3> &vcm,
                const blitz::TinyVector<double,3> &L);
void pkdCalcDistance(PKD pkd, blitz::TinyVector<double,3> dCenter, int bPeriodic);
uint_fast32_t pkdCountDistance(PKD pkd, double r2i, double r2o );
void pkdCalcCOM(PKD pkd, const blitz::TinyVector<double,3> &dCenter, double dRadius, int bPeriodic,
                blitz::TinyVector<double,3> &com, blitz::TinyVector<double,3> &vcm,
                blitz::TinyVector<double,3> &L,
                double &M, uint64_t &N);
void pkdCalcMtot(PKD pkd, double *M, uint64_t *N);
void pkdResetCOM(PKD pkd, blitz::TinyVector<double,3> r_com, blitz::TinyVector<double,3> v_com);
void pkdInitializeEOS(PKD pkd);
void pkdCopySPHOptionsToDevice(PKD pkd, SPHOptions *SPHoptions, int bGPU);
void pkdUpdateGasValues(PKD pkd, struct pkdKickParameters *kick, SPHOptions *SPHoptions);
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
void pkdSetupInterpScale(PKD pkd,double dBoxSize,double mrMax);

struct outGetParticles { /* Array of these */
    uint64_t id;
    float    mass, phi;
    blitz::TinyVector<float,3> r, v;
};
int pkdGetParticles(PKD pkd, int nIn, uint64_t *ID, struct outGetParticles *out);

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
