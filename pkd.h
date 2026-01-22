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

#include <cstdint>
#include <stdint.h>
#include <string.h>
#include <optional>
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
#include "core/ephemeral.h"
#include "core/integerize.h"
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

#define CID_PARTICLE 0
#define CID_CELL 1
#define CID_PARTICLE2 8
#define CID_CELL2 9
#define CID_HEALPIX 7
#define CID_GROUP 2
#define CID_SHRINK 3
#define CID_RM 3
#define CID_BIN 4
#define CID_SHAPES 5
#define CID_PK 2
#define CID_GridLinFx 2
#define CID_GridLinFy 10
#define CID_GridLinFz 11
#define CID_PNG 2
#define CID_SADDLE_BUF 3
#define CID_TREE_ROOT 3

#define MAX_TIMERS 10

/*
** Memory models.  Each is a bit mask that indicates that additional fields should be
** added to the particle structure.
*/
#define PKD_MODEL_VELOCITY (1 << 0)     /* Velocity Required */
#define PKD_MODEL_ACCELERATION (1 << 1) /* Acceleration Required */
#define PKD_MODEL_POTENTIAL (1 << 2)    /* Potential Required */
#define PKD_MODEL_GROUPS (1 << 3)       /* Group profiling */
#define PKD_MODEL_MASS (1 << 6)         /* Mass for each particle */
#define PKD_MODEL_DENSITY (1 << 7)      /* Density for each particle */
#define PKD_MODEL_BALL (1 << 8)         /* Ball for each particle */
#define PKD_MODEL_SOFTENING (1 << 9)    /* Softening for each particle */
#define PKD_MODEL_VELSMOOTH (1 << 10)   /* Velocity Smoothing */
#define PKD_MODEL_MFM (1 << 11)         /* Sph Fields */
#define PKD_MODEL_MFV (1 << 12)         /* Sph Fields */
#define PKD_MODEL_STAR (1 << 13)        /* Star Fields */
#define PKD_MODEL_PARTICLE_ID (1 << 14) /* Particles have a unique ID */
#define PKD_MODEL_UNORDERED (1 << 15)   /* Particles do not have an order */
#define PKD_MODEL_INTEGER_POS (1 << 16) /* Particles do not have an order */
#define PKD_MODEL_BH (1 << 17)          /* BH fields */
#define PKD_MODEL_GLOBALGID (1 << 18)   /* Global group identifier per particle */
#define PKD_MODEL_NEWRUNG (1 << 19)     /* New rung field */

#define PKD_MODEL_NODE_MOMENT (1 << 24)  /* Include moment in the tree */
#define PKD_MODEL_NODE_ACCEL (1 << 25)   /* mean accel on cell (for grav step) */
#define PKD_MODEL_NODE_VEL (1 << 26)     /* center of mass velocity for cells */
#define PKD_MODEL_NODE_SPHBNDS (1 << 27) /* Include 3 extra bounds in tree */

#define PKD_MODEL_NODE_BND (1 << 28)  /* Include normal bounds in tree */
#define PKD_MODEL_NODE_VBND (1 << 29) /* Include velocity bounds in tree for phase-space density*/
#define PKD_MODEL_NODE_BOB (1 << 30)  /* Ball of Balls, or Box of Balls */

#define PKD_MODEL_NEW_SPH (1ULL << 32) /* New Sph Fields */
/* P and c if not doing onthefly prediction */
/* oldRho if doing ISPH */
/* expImb2 if doing interface correction */
/* T if doing interface correction and not onthefly prediction*/
/* vpred if doing consistent prediction */
#define PKD_MODEL_STRENGTH (1ULL << 38) /* New Sph Fields */
/* Spred if doing consistent prediction */
#define PKD_MODEL_PEAKPRESSURE (1ULL << 40)

#define MAX_RUNG 63

/*
** Here we define some special reserved nodes. Node-0 is a sentinel or null node, node-1
** is here defined as the ROOT of the local tree (or top tree), node-2 is unused and
** node-3 is the root node "fixed" tree.
*/
#define FIXROOT 3
#define VAROOT 3
#define ROOT 1
#define NRESERVED_NODES MAX_RUNG + 1

typedef struct partLightCone {
    uint64_t id;
    float pos[3];
    float vel[3];
    float pot;
} LIGHTCONEP;

#define NMAX_OPENCALC 1000

/*
** Components required for tree walking.
*/
struct CSTACK {
    using list_type = typename clList::list_type; // the boost::intrusive::list alias
    // CSTACK() = default;
    explicit CSTACK(list_type& freeList) : cl(freeList) {}
    // clList is non-copyable (intrusive); but movable
    CSTACK(const CSTACK&) = delete;
    CSTACK& operator=(const CSTACK&) = delete;
    CSTACK(CSTACK&&) noexcept = default;
    CSTACK& operator=(CSTACK&&) noexcept = default;

    clList cl;
    LOCR L;
    float dirLsum;
    float normLsum;
    int iNodeIndex;
};

#if defined(USE_SIMD) && !defined(__CUDACC__)
typedef struct {
    struct PMOMC {
        dvec m;
        dvec xx, yy, xy, xz, yz;
        dvec xxx, xyy, xxy, yyy, xxz, yyz, xyz;
        dvec xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
        dvec zz;
        dvec xzz, yzz, zzz;
        dvec xxzz, xyzz, xzzz, yyzz, yzzz, zzzz;
    } ewm;
    struct PEWALDVARS {
        dvec fEwCut2, fInner2, alpha, alpha2, ialpha, k1, ka;
        dvec Q4xx, Q4xy, Q4xz, Q4yy, Q4yz, Q4zz, Q4, Q3x, Q3y, Q3z, Q2;
    } ewp;
} ewaldSIMD;
#endif

/*
** This is the temporary group table used when Grasshopping.
** We eventually contruct a proper table.
*/
typedef remoteID GHtmpGroupTable;

typedef struct {
    uint64_t nTotal; /* Total particles in this group */
    remoteID id;     /* Owner (or myself) */
    remoteID rmt;
    uint32_t iGlobalId; /* Global unique group id */
    uint32_t nLocal;    /* Local to this processor */
    uint32_t iTreeRoot; /* Our local tree root */
    uint32_t iAllRoots;
    uint16_t nRemote; /* Number of remote partners */
    uint8_t bNeedGrav : 1;
    uint8_t bComplete : 1;
    float fMass;
    float fRMSRadius;
    double dEnergy;
    blitz::TinyVector<double, 3> rref;
    blitz::TinyVector<double, 3> ravg;
    blitz::TinyVector<double, 3> rmin;
    blitz::TinyVector<double, 3> rmax;
    blitz::TinyVector<double, 3> rcom;
    blitz::TinyVector<double, 3> vcom;
} HopGroupTable;

typedef struct {
    remoteID key;
    remoteID name;
    uint32_t iLink;
} FOFRemote;

typedef struct {
    blitz::TinyVector<float, 3> rPot;
    float minPot;
    blitz::TinyVector<float, 3> rcen;
    blitz::TinyVector<float, 3> rcom;
    blitz::TinyVector<float, 3> vcom;
    blitz::TinyVector<float, 3> angular;
    blitz::TinyVector<float, 6> inertia;
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
    blitz::TinyVector<double, 3> rPot;
    float fMass;
    int nBH;
} VeryTinyGroupTable;

// typedef struct {
//     } SmallGroupTable;

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

    struct saddle_point* sp;
    struct saddle_point_buffer* buf;
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
    struct remote_root_id* treeRoots;
    int nSaddlePoints;
    int* sp;
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
    struct psGroup* pGroup;
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
    blitz::TinyVector<double, 3> L;
    double vel_radial;
    double vel_radial_sigma;
    double vel_tang_sigma;
    uint64_t nParticles;
} PROFILEBIN;

typedef struct shapesBin {
    blitz::TinyVector<double, 3> com;
    double dMassEnclosed;
    blitz::TinyVector<blitz::TinyVector<double, 3>, 3> dInertia;
    blitz::TinyVector<blitz::TinyVector<double, 3>, 3> ell_matrix;
    blitz::TinyVector<double, 3> ell_center;
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

namespace base_type {
using force_t = float;
using force_vector_t = blitz::TinyVector<force_t, 3>;
using force_array_t = blitz::Array<force_t, 3>;
using position_t = double;
using position_vector_t = blitz::TinyVector<position_t, 3>;
using shape_vector_t = blitz::TinyVector<int, 3>;
using mass_t = float;
using mass_array_t = blitz::Array<mass_t, 3>;
} // namespace base_type

class pkdContext : public Integerize {
  public:
    explicit pkdContext(mdl::mdlClass* mdl, int nStore, int nDecoratorStore, int nBucket, uint64_t nMinTotalStore, uint64_t nMinEphemeral,
                        uint32_t nEphemeralBytes, int nTreeBitsLo, int nTreeBitsHi, int iCacheSize, int iCacheMaxInflight, int iWorkQueueSize,
                        mdl::MPI_RECV_MODE mpi_recv_mode, int mpi_recv_count, bool mpi_RDMA, const blitz::TinyVector<double, 3>& fPeriod, uint64_t nDark,
                        uint64_t nGas, uint64_t nStar, uint64_t nBH, uint64_t mMemoryModel, uint32_t nIntegerFactor);
    virtual ~pkdContext();
    void set_factor(std::uint32_t factor);

  public:
    mdl::mdlClass* mdl;
    auto Self() const { return mdl->Self(); }
    auto Threads() const { return mdl->Threads(); }

  protected:                      // Support for memory models
    int NodeAddStruct(int n);     // Add a NODE structure: assume double alignment
    int NodeAddDouble(int n = 1); // Add n doubles to the node structure
    int NodeAddFloat(int n = 1);  // Add n floats to the node structure
    int NodeAddInt64(int n = 1);  // Add n 64-bit integers to the node structure
    int NodeAddInt32(int n = 1);  // Add n 32-bit integers to the node structure
  protected:
    // Storage Segments: there are three different segments for each thread:
    // 1. The main particle store, which contains the particles and their properties.
    // 2. The ephemeral store, which contains temporary data for each particle.
    // 3. The overflow store, which is used if the grids for the would not fit into the ephermeral store.
    //    The storage in this segment is reused for (some of) the tree nodes after the IC generation is done.
    static constexpr auto nSegments = 4;
    blitz::TinyVector<void*, nSegments> pSegments;
    blitz::TinyVector<size_t, nSegments> nSegmentSizeBytes;

  public:
    auto SegmentSizeBytes(int i) const { return nSegmentSizeBytes[i]; }
    template <typename T = void>
    T* Segment(int i) const {
        return reinterpret_cast<T*>(pSegments[i]);
    }
    template <typename T>
    auto Ephemeral() {
        return Segment<T>(1);
    }
    std::optional<ephemeral_buffer_resource> ephemeral;

  protected:
    PARTICLE* pTempPRIVATE = nullptr;
    size_t nEphemeralBytes = 0; /* per-particle */
    size_t storageSize;

  public:
    particleStore particles;
    decoratorStore decorators;
    treeStore tree;
    int FreeStore() { return particles.capacity(); }
    int FreeDecoratorStore() { return decorators.capacity(); }
    // I do not like at all the name "Local" for indicating the number of particles...
    int Local() { return particles.size(); }
    int DecoratorLocal() { return decorators.size(); }
    int SetLocal(int n);
    int AddLocal(int n);
    auto EphemeralElementBytes() { return nEphemeralBytes; }
    auto EphemeralSize() { return 1ul * EphemeralElementBytes() * (FreeStore() + 1); }
    static constexpr auto MaxNodeSize() { return sizeof(KDN) + 2 * sizeof(Bound) + sizeof(FMOMR) + 6 * sizeof(double) + sizeof(SPHBNDS); }
    auto NodeSize() { return tree.ElementSize(); }
    auto DecoratorSize() { return decorators.ElementSize(); }
    auto DecoratorMemory() { return decorators.ElementSize() * (decorators.capacity() + 1); }
    auto ParticleMemory() { return (particles.ParticleSize() + EphemeralElementBytes()) * (FreeStore() + 1); }
    PARTICLE* Particle(int i) { return particles.Element(i); }
    decorators::FIELDS* Decorator(int i) { return decorators.Element(i); }

    void SaveParticle(PARTICLE* a) { memcpy(pTempPRIVATE, a, particles.ParticleSize()); }
    void LoadParticle(PARTICLE* a) { memcpy(a, pTempPRIVATE, particles.ParticleSize()); }
    void CopyParticle(PARTICLE* a, PARTICLE* b) { memcpy(a, b, particles.ParticleSize()); }
    void CopyDecorator(decorators::FIELDS* a, decorators::FIELDS* b) { memcpy(a, b, decorators.DecoratorSize()); }

  public:
    auto Nodes() const { return tree.size(); }
    auto Node(KDN* pBase, int iNode) { return tree.Element(pBase, iNode); }
    auto TreeNode(int iNode) { return tree[iNode]; }
    auto TreeAlignNode() { return tree.align(); }
    auto TreeAllocNode(int n = 1) { return tree.append(n); }
    void TreeAllocNodePair(int* iLeft, int* iRight) {
        *iLeft = TreeAllocNode();
        *iRight = TreeAllocNode();
    }
    auto TreeAllocRootNode() {
        if ((Nodes() & 1) == 0) TreeAllocNode();
        return TreeAllocNode();
    }
    // I/O
  public:
    void Restore(uint64_t iElement, const std::string& filename, uint64_t iBeg, uint64_t iEnd);
    void RestoreDecorator(uint64_t iElement, const std::string& filename, uint64_t iBeg, uint64_t iEnd);
    void GroupStatsSave(int iGroup, const std::string& filename, int iSegment, int nSegment);

    // Rockstar Analysis
  protected:
    uint64_t nRsElements;

  public:
    void RsHaloIdStart(uint64_t nElements, bool bAppend = false);
    void RsHaloIdFinish(uint64_t nElements);
    void RsHaloIdRead(uint64_t iElement, const std::string& filename, uint64_t iBeg, uint64_t iEnd);
    void RsIdStart(uint64_t nElements, bool bAppend = false);
    void RsIdFinish(uint64_t nElements);
    void RsIdRead(uint64_t iElement, const std::string& filename, uint64_t iBeg, uint64_t iEnd);
    void RsIdSave(int iGroup, const std::string& filename, int iSegment, int nSegment);
    void RsReorder(const uint64_t* pOrds);
    void RsExtract(const std::string& filename, int iSegment, uint64_t* pOrds);

    // Lightcone
  protected:
    void flushLightCone(int octant);
    void addToLightCone(double dvFac, double* r, float fPot, PARTICLE* p, int bParticleOutput);

  public:
    void LightConeOpen(const char* name_template, int nSideHealpix);
    void LightConeClose(const char* healpixname_template);
    void ProcessLightCone(PARTICLE* p, float fPot, double dLookbackFac, double dLookbackFacLCP, double dDriftDelta, double dKickDelta, double dBoxSize,
                          int bLightConeParticles, blitz::TinyVector<double, 3> hlcp, double tanalpha2);
    void LightConeVel(double dBoxSize);

  protected:
    /// @brief Partition particles array between pLower and pUpper (inclusive).
    /// @tparam pivot_t Allows for floating point or integer positions.
    /// @param pkd pkdContext object.
    /// @param pLower Starting index of partition.
    /// @param pUpper End index of partition (inclusive).
    /// @param d Dimension to split in.
    /// @param pivot Values lower then pivot are sorted to the left of the array, values larger then pivot are sorted to the right of the array.
    /// @return The index of the split.
    template <class pivot_t>
    auto PartPart(int pLower, int pUpper, int d, pivot_t pivot) -> int {
        auto pi = particles.begin() + pLower;
        auto pj = particles.begin() + pUpper;
        return std::partition(pi, pj + 1, [pivot, d](auto& p) { return p.template raw_position<pivot_t>(d) < pivot; }) - particles.begin();
    }

    void BuildTemp(int iNode, int bucketSize, int nGroup, double maxBucketWidth);
    void Create(int iRoot, double ddHonHLimit);
    void BuildFromTemplate(int iNode, int M, int nGroup, int iTemplate);
    void TreeUpdateFlagBoundsRecurse(uint32_t uRoot, const SPHOptions* SPHoptions);
    auto InitializeRootCommon(uint32_t uRoot) {
        auto pRoot = tree[uRoot];
        pRoot->set_local();
        pRoot->set_bound(bnd);
        return pRoot;
    }

    void fetch_forces(int cid, int nGrid, base_type::force_array_t& forces, const base_type::shape_vector_t& lower);
    void flush_masses(int nGrid, const base_type::mass_array_t& masses, const base_type::shape_vector_t& lower, int fold);
    void FofRemoteSearch(double dTau2, int bPeriodic, int nReplicas, int nBucket);
    auto FofGatherLocal(std::vector<int>& S, double fBall2, blitz::TinyVector<double, 3> r, uint32_t iGroup, uint32_t iTail, uint32_t* Fifo,
                        int* pbCurrFofContained, const Bound::coord_type fMinFofContained, const Bound::coord_type fMaxFofContained) -> uint32_t;
    auto GetDistance2(particleStore::Particle& p, blitz::TinyVector<double, 3> dCenter, int bPeriodic) -> double;
    auto GatherMass(remoteID* S, double fBall, blitz::TinyVector<double, 3> r, int bPeriodic, blitz::TinyVector<double, 3> dPeriod) -> double;
    auto gatherMass(remoteID* S, double fBall2, double ri2, double ro2, blitz::TinyVector<double, 3> r) -> double;
    auto getCell(int iCell, int id) { return (id == Self()) ? tree[iCell] : tree[mdl->Fetch<KDN>(CID_CELL, iCell, id)]; }

  public:
    void TreeBuild(int nBucket, int nGroup, uint32_t uRoot, uint32_t uTemp, double ddHonHLimit);
    void TreeBuildByGroup(int nBucket, int nGroup);
    void TreeUpdateFlagBounds(uint32_t uRoot, const SPHOptions* SPHoptions);
    void TreeInitMarked();
    auto DistribTopTree(uint32_t uRoot, uint32_t nTop, const KDN* pTop, int allocateMemory) -> uint32_t;
    void DistribRoot(const double*, const MOMC*);
    void DumpTrees(int bOnlyVA, uint8_t uRungDD);
    void CombineCells1(treeStore::NodePointer pkdn, treeStore::NodePointer p1, treeStore::NodePointer p2);
    void CombineCells2(treeStore::NodePointer pkdn, treeStore::NodePointer p1, treeStore::NodePointer p2);
    void GroupOrder(uint32_t* iGrpOffset);
    auto ActiveOrder() -> uint64_t;

    void UpdateGasValues(struct pkdKickParameters* kick, SPHOptions* SPHoptions);
    void ResetCOM(blitz::TinyVector<double, 3> r_com, blitz::TinyVector<double, 3> v_com);

    auto TotalMass() -> double;
    [[maybe_unused]] auto GridCreateFFT(int nGrid) -> std::tuple<size_t, int, int, int, int>;
    void GridDeleteFFT();
    void SetLinearGrid(double a0, double a, double a1, double dBSize, int nGrid, int iSeed, int bFixed, float fPhase);
    void GenerateLinearGrid(double a, double a_next, double Lbox, int iSeed, int bFixed, float fPhase, int bRho);
    void LinearKick(vel_t dtOpen, vel_t dtClose, int iAssignment = 3);
    void MeasureLinearPk(int nGrid, double dA, double dBoxSize, int nBins, int iSeed, int bFixed, float fPhase, double* fK, double* fPower, uint64_t* nPower);
    void AddLinearSignal(int iGrid, int iSeed, bool bFixed, float fPhase, double Lbox, double a);
    void AssignMass(uint32_t iLocalRoot, int iAssignment, int iGrid, float fDelta, int fold);
    void DensityContrast(double dTotalMass, int iGrid, bool k = false);
    void Interlace(int iGridTarget, int iGridSource);
    void WindowCorrection(int iAssignment, int iGrid);
    void GridBinK(int nBins, int iGrid, std::span<double> fK, std::span<double> fPower, std::span<uint64_t> nPower);
    void BispectrumSelect(int iGridTarget, int iGridSource, double kmin, double kmax);
    auto BispectrumCalculate(int iGrid1, int iGrid2, int iGrid3) -> double;
    void BispectrumNormalize(int iGridTarget, double kmin, double kmax);
    void NewFof(double dTau2, int nMinMembers, int bPeriodic, int nReplicas, int nBucket);
    auto FofPhases() -> bool;
    auto FofFinishUp(int nMinGroupSize) -> uint64_t;
    auto HopFinishUp(int nMinGroupSize, int bPeriodic, blitz::TinyVector<double, 3> dPeriod) -> int;
    void HopTreeBuild(int nBucket, int nGroup);
    auto HopGravity(double dTime, double dThetaMin, double* pdFlop, double* pdPartSum, double* pdCellSum) -> int;
    auto HopUnbind(double dTime, int nMinGroupSize, int bPeriodic, blitz::TinyVector<double, 3> dPeriod) -> int;
    void CalculateGroupStats(int bPeriodic, blitz::TinyVector<double, 3> dPeriod, std::array<double, 2> rEnvironment, uint64_t iGlobalStart);
    void Profile(uint8_t uRungLo, uint8_t uRungHi, const blitz::TinyVector<double, 3>& dCenter, const double* dRadii, int nBins,
                 const blitz::TinyVector<double, 3>& com, const blitz::TinyVector<double, 3>& vcm, const blitz::TinyVector<double, 3>& L);
    void CalcDistance(blitz::TinyVector<double, 3> dCenter, int bPeriodic);
    void CalcCOM(const blitz::TinyVector<double, 3>& dCenter, double dRadius, int bPeriodic, blitz::TinyVector<double, 3>& com,
                 blitz::TinyVector<double, 3>& vcm, blitz::TinyVector<double, 3>& L, double& M, uint64_t& N);
    uint8_t MaxRung();
    void SetGlobalRung(uint8_t minDt);
    void DensityStep(uint8_t uRungLo, uint8_t uRungHi, int iMaxRung, double dDelta, double dEta, double dRhoFac);
    void AccelStep(uint8_t uRungLo, uint8_t uRungHi, double dDelta, int iMaxRung, double dEta, double dVelFac, double dAccFac, int bDoGravity, int bEpsAcc,
                   double dhMinOverSoft);
    void Kick(double dTime, double dDelta, int bDoGas, double, double, double, uint8_t uRungLo, uint8_t uRungHi);
    void Drift(int iRoot, double dTime, double dDelta, double, double, int bDoGas, int bGasEvolveDensity);

    void WriteViaNode(int iNode);
    void BHInit(uint8_t uRungMax);
    auto PlaceBHSeed(double dTime, double dScaleFactor, uint8_t uRungMax, double dDenMin, double dBHMhaloMin, double dTau, double dBHSeedMass) -> int;
    void BHReposition();
    void BHAccretion(double dScaleFactor);
    void ChemCompInit(double H, double He, double C, double N, double O, double Ne, double Mg, double Si, double Fe, double Metallicity);

#ifdef GRACKLE
    void GrackleInit(int bComove, double dScaleFactor, char* achCoolingTable, UNITS units);
    void GrackleUpdate(double dScaleFactor, char* achCoolingTable, UNITS units);
    void GrackleCooling(particleStore::ParticleReference& p, double pDelta, double dTuFac);
#endif

    void swap_classes(int id);
#ifdef OPTIM_REORDER_IN_NODES
    void ReorderInNodes();
#endif

  public:
    int nRejects;
    int nDecoratorRejects;
    int nActive;
    blitz::TinyVector<uint64_t, MAX_RUNG + 1> nRung;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nBH;
    blitz::TinyVector<double, 3> fPeriod;
    int iTopTree[NRESERVED_NODES];
    int nNodesFull; /* number of nodes in the full tree (including very active particles) */
    Bound bnd;
    Bound vbnd;
    /*
    ** Light cone variables.
    */
    int nLayerMax;
    std::vector<int> nBoxLC;
    std::vector<double> lcOffset0, lcOffset1, lcOffset2;
    asyncFileInfo afiLightCone[8];
    int nLightConeMax, mLightConeOctants;
    blitz::TinyVector<LIGHTCONEP*, 8> pLightCone;
    blitz::TinyVector<int, 8> nLightCone;
    int64_t nHealpixPerDomain;
    int64_t nSideHealpix;
    healpixData* pHealpixData;
    gsl_spline* interp_scale = nullptr; // interpolation table for 1/a given r in the lightcone
    gsl_interp_accel* interp_accel = nullptr;

    /*
    ** Advanced memory models
    */
    int bNoParticleOrder;
    int bIntegerPosition;
    bool bStars;      // Stellar particles are allowed
    bool bBlackholes; // Blackhole particles are allowed

    /*
    ** Tree walk variables.
    */
    std::vector<CSTACK> S;
    ilpList ilp;
    ilcList ilc;
    ilcList ill;
    clList::free_list clFreeList;
    clList cl;
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
    blitz::TinyVector<double, 3> dEnergyF, dEnergyL;

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
    struct smGroupArray* ga;
    Bound bndInterior; /* this gets calculated at the start of Fof for now but should be done elsewhere */
    /*
    ** Some variables needed for NewFof().
    */
    uint32_t iRemoteGroup, nMaxRemoteGroups;
    FOFRemote* tmpFofRemote;
    TinyGroupTable* tinyGroupTable;
    VeryTinyGroupTable* veryTinyGroupTable;

    GHtmpGroupTable* tmpHopGroups;
    HopGroupTable* hopGroups;
    int* hopRootIndex;
    int hopSavedRoots;
    remoteID* hopRoots;

    struct saddle_point_list saddle_points;
    int nRm;
    int nMaxRm;
    int nBins;

    FOFBIN* groupBin;

    PROFILEBIN* profileBins;

    CSM csm;

    // IA: we add here the needed cooling information available to all procs
    struct cooling_function_data* cooling = nullptr;
    struct cooling_tables* cooling_table = nullptr;
#ifdef GRACKLE
    chemistry_data* grackle_data = nullptr;
    chemistry_data_storage* grackle_rates = nullptr;
    grackle_field_data* grackle_field = nullptr;
    code_units* grackle_units = nullptr;
#endif
    struct StellarEvolutionData* StelEvolData = nullptr;

#ifdef USE_CUDA
    std::unique_ptr<CudaClient> cudaClient;
#endif
#ifdef USE_METAL
    std::unique_ptr<MetalClient> metalClient;
#endif
#ifdef MDL_FFTW
    MDLFFT fft;
#endif

    SPHOptions SPHoptions;
#ifdef HAVE_EOSLIB_H
    EOSmaterial* materials[EOS_N_MATERIAL_MAX] = {NULL};
#endif
  public:
    auto& NodeBOB(KDN* n) { return tree.BOB(n); }
    auto& NodeBOB(const KDN* n) { return tree.BOB(n); }

  public:
    void ActiveRung(int iRung, bool bGreater) {
        uMinRungActive = iRung;
        uMaxRungActive = bGreater ? 255 : iRung;
        particles.ActiveRung(uMinRungActive, uMaxRungActive);
    }
};
using PKD = pkdContext*;

#ifdef FIXED_SIZE_TREES
class treeCACHEhelper : public mdl::arrayCACHEhelper {
  protected:
    pkdContext* const pkd;

  public:
    explicit treeCACHEhelper(pkdContext* pkd) : arrayCACHEhelper(pkd->tree.ElementSize(), pkd->tree, pkd->Local(), false), pkd(pkd) {}
};
#else
class treeCACHEhelper : public mdl::CACHEhelper {
  protected:
    pkdContext* const pkd;

  public:
    explicit treeCACHEhelper(pkdContext* pkd) : CACHEhelper(pkd->NodeSize(), pkd->tree.size(), false), pkd(pkd) {}
    virtual void* get_element(uint32_t i) override final { return pkd->tree[i]; }
};
#endif

class particleCACHEhelper : public mdl::arrayCACHEhelper {
  protected:
    pkdContext* const pkd;

  public:
    explicit particleCACHEhelper(pkdContext* pkd, bool bModify = false)
        : arrayCACHEhelper(pkd->particles.ParticleSize(), pkd->particles, pkd->Local(), bModify), pkd(pkd) {}
};

static inline int pkdIsRungActive(PKD pkd, uint8_t uRung) {
    return uRung >= pkd->uMinRungActive && uRung <= pkd->uMaxRungActive;
}
static inline int pkdIsActive(PKD pkd, PARTICLE* p) {
    return pkdIsRungActive(pkd, p->uRung);
}

/*
** The size of a particle is variable based on the memory model.
** The following three routines must be used instead of accessing pStore
** directly.  pkdParticle will return a pointer to the i'th particle.
** The Size and Base functions are intended for cache routines; no other
** code should care about sizes of the particle structure.
*/

static inline FIO_SPECIES pkdSpecies(PKD pkd, PARTICLE* p) {
    return pkd->particles.species(p);
}
static inline int pkdiMat(PKD pkd, PARTICLE* p) {
    return pkd->particles.iMat(p);
}

static inline int pkdIsDeleted(PKD pkd, PARTICLE* p) {
    return (pkdSpecies(pkd, p) == FIO_SPECIES_UNKNOWN);
}

static inline int pkdIsNew(PKD pkd, PARTICLE* p) {
    return (p->iOrder == IORDERMAX);
}

/*
** From tree.c:
*/
void pkdVATreeBuild(PKD pkd, int nBucket);
void pkdOpenCloseCaches(PKD pkd, int bOpen, int bFixed);

/*
** From pkd.c:
*/
size_t pkdClCount(PKD pkd);
size_t pkdClMemory(PKD pkd);
size_t pkdIllMemory(PKD pkd);
size_t pkdIlcMemory(PKD pkd);
size_t pkdIlpMemory(PKD pkd);
void pkdReadFIO(PKD pkd, FIO fio, uint64_t iFirst, int nLocal, double dvFac, double dTuFac, FIO_SPECIES s);
void pkdSetSoft(PKD pkd, double dSoft);
void pkdSetSmooth(PKD pkd, double dSmooth);
void pkdSetCrit(PKD pkd, double dCrit);
void pkdEnforcePeriodic(PKD, Bound);
void pkdPhysicalSoft(PKD pkd, double dSoftMax, double dFac, int bSoftMaxMul);
int pkdWeight(PKD, int, double, int, int, int, int*, int*, double*, double*);
void pkdCountVA(PKD, int, double, int*, int*);
int pkdLowerPart(PKD, int, double, int, int);
int pkdUpperPart(PKD, int, double, int, int);
int pkdWeightWrap(PKD, int, double, double, int, int, int, int*, int*);
int pkdLowerPartWrap(PKD, int, double, double, int, int);
int pkdUpperPartWrap(PKD, int, double, double, int, int);
int pkdLowerOrdPart(PKD, uint64_t, int, int);
int pkdUpperOrdPart(PKD, uint64_t, int, int);

/*#define PEANO_HILBERT_KEY_MAX 0x3ffffffffffull*/  /* 2d */
#define PEANO_HILBERT_KEY_MAX 0x7fffffffffffffffull /* 3d */
void pkdPeanoHilbertDecomp(PKD pkd, int nRungs, int iMethod);
void pkdRungOrder(PKD pkd, int iRung, total_t* nMoved);
std::pair<int, int> pkdColRejects(PKD, int);
std::pair<int, int> pkdColRejects_Old(PKD, int, double, double, int);

std::pair<int, int> pkdSwapRejects(PKD, int);
int pkdSwapSpace(PKD);
int pkdSwapDecoratorSpace(PKD);
int pkdActive(PKD);
int pkdInactive(PKD);

std::pair<int, int> pkdColOrdRejects(PKD, uint64_t, int);
void pkdLocalOrder(PKD, uint64_t iMinOrder, uint64_t iMaxOrder);
void pkdCheckpoint(PKD pkd, const char* fname, const char* fname_decorator);
void pkdWriteHeaderFIO(PKD pkd, FIO fio, double dScaleFactor, double dTime, fioSpeciesList nTotal, fioSpeciesList nThisFile, double dBoxSize, double h,
                       int nProcessors, UNITS units);
uint32_t pkdWriteFIO(PKD pkd, FIO fio, double dvFac, double dTuFac, Bound bnd, fioSpeciesList nWritten);
void pkdWriteFromNode(PKD pkd, int iNode, FIO fio, double dvFac, double dTuFac, Bound bnd, fioSpeciesList nWritten);
auto pkdPackArray(PKD pkd, int iSize, void* vBuff, int* piIndex, int n, PKD_FIELD field, int iUnitSize, double dvFac, int bMarked) -> std::byte*;
void pkdSendArray(PKD pkd, int iNode, PKD_FIELD field, int iUnitSize, double dvFac, int bMarked);
void* pkdRecvArray(PKD pkd, int iNode, void* pDest, int iUnitSize);
void pkdGravAll(PKD pkd, struct pkdKickParameters* kick, struct pkdLightconeParameters* lc, struct pkdTimestepParameters* ts, double dTime, int nReps,
                int bPeriodic, int bGPU, int bEwald, int iRoot1, int iRoot2, double fEwCut, double fEwhCut, double dThetaMin, SPHOptions* SPHoptions,
                uint64_t* pnActive, double* pdPart, double* pdPartNumAccess, double* pdPartMissRatio, double* pdCell, double* pdCellNumAccess,
                double* pdCellMissRatio, double* pdFlop, uint64_t* pnRung);
void pkdCalcEandL(PKD pkd, double& T, double& U, double& Eth, blitz::TinyVector<double, 3>& L, blitz::TinyVector<double, 3>& F, double& W);
void pkdGravEvalPP(const PINFOIN& Part, ilpTile& tile, PINFOOUT& Out);
void pkdDensityEval(const PINFOIN& Part, ilpTile& tile, PINFOOUT& Out, SPHOptions* SPHoptions);
void pkdDensityCorrectionEval(const PINFOIN& Part, ilpTile& tile, PINFOOUT& Out, SPHOptions* SPHoptions);
void pkdSPHForcesEval(const PINFOIN& Part, ilpTile& tile, PINFOOUT& Out, SPHOptions* SPHoptions);
void pkdGravEvalPC(const PINFOIN& pPart, ilcTile& tile, PINFOOUT& pOut);
void pkdEndTimestepIntegration(PKD pkd, struct inEndTimestep in);
void pkdReorderDecorator(PKD pkd);
void pkdSwapAll(PKD pkd, int idSwap);
void pkdInitCosmology(PKD pkd, const struct csmVariables* cosmo);
void pkdInitLightcone(PKD pkd, int bBowtie, int bLightConeParticles, double dBoxSize, double dRedshiftLCP, double alphaLCP, blitz::TinyVector<double, 3> hLCP,
                      int mOctants);
void pkdZeroNewRung(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdCooling(PKD pkd, double, double, int, int, int, int);
#define CORRECTENERGY_IN 1
#define CORRECTENERGY_OUT 2
#define CORRECTENERGY_SPECIAL 3
#define CORRECTENERGY_IC_MESHLESS 4
int pkdUpdateRung(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung, int iMaxRung, uint64_t* nRungCount);
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
    T.f = fabsf(fiDelta) * fT;
    if (T.f >= 1.0) return 0;
    iRung = 126 - T.ieee.exponent; /* -log2(d) */
    if (iRung > uMaxRung) return uMaxRung;
    else return iRung;
}
int pkdOrdWeight(PKD pkd, uint64_t iOrdSplit, int iSplitSide, int iFrom, int iTo, int* pnLow, int* pnHigh);
void pkdDeleteParticle(PKD pkd, particleStore::ParticleReference& p);
int pkdIsGas(PKD, PARTICLE*);
int pkdIsDark(PKD, PARTICLE*);
int pkdIsStar(PKD, PARTICLE*);
int pkdIsBH(PKD, PARTICLE*);

int pkdDeepestPot(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, double* r, float* fPot);
void pkdInitializeEOS(PKD pkd);
void pkdCopySPHOptionsToDevice(PKD pkd, SPHOptions* SPHoptions, int bGPU);
#ifdef MDL_FFTW
float getLinAcc(PKD pkd, MDLFFT fft, int cid, double r[3]);
#endif
void pkdOutPsGroup(PKD pkd, char* pszFileName, int iType);

void pkdLightCone(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, double dLookbackFac, double dLookbackFacLCP, double* dtLCDrift, double* dtLCKick);
void pkdSetupInterpScale(PKD pkd, double dBoxSize, double mrMax);

#define matrix_vector_mult(b, mat, a)                                                                                                                          \
    do {                                                                                                                                                       \
        int i;                                                                                                                                                 \
        for (i = 0; i < 3; i++) {                                                                                                                              \
            int j;                                                                                                                                             \
            (b)[i] = 0.0;                                                                                                                                      \
            for (j = 0; j < 3; j++) (b)[i] += (mat)[i][j] * (a)[j];                                                                                            \
        }                                                                                                                                                      \
    } while (0)

#endif
