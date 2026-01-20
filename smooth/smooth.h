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

#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include <memory>
#include <sys/types.h>
#include <vector>
#include <boost/unordered/unordered_flat_set.hpp>
#include "pkd.h"
#include "smooth/smoothfcn.h"
#include "group/group.h"
#include "blitz/array.h"

#define NNLIST_INCREMENT 200 /* number of extra neighbor elements added to nnList */

// To use the new Smooth interace you must provide a "helper" class that defines:
// - How to smooth (derived from SmoothFunction) by providing a Smooth() method
// - (optionally) How to access the particle cache (derived from particleCACHEhelper).
// Normally you derive a new class from SmoothHelper or more specialised classes like
// SmoothHelperDensity, SmoothHelperAsymmetric, etc. These classes provide default
// implementations for smoothing density, asymmetric smoothing, etc.
//
// When it comes time to call Smooth(), ReSmooth(), or NodeSmooth(), you also
// need to provide a "policy" class that defines which particles are smoothed
// and which particles are accepted as neighbors. See the SmoothEverything class
// below for an example. You can derive new classes from that one to provide
// filtering based on particle type (gas, BH, etc) and/or active status.

// Policy classes for filtering particles during node smoothing
struct SmoothEverything {
    // Returns the open range [start, end) of particles in the node to process
    // Certain optimisations (see below) allow us to avoid checking particle types
    static inline auto range(treeStore::NodePointer node) { return std::pair{node->lower(), node->upper() + 1}; }
    // Filter function to select nodes of interest (here, all particles)
    // If the filter is done in range(), this can be a no-op (always true)
    static constexpr auto may_have_active(treeStore::NodePointer&) -> bool { return true; }
    // This deterines if the particle is active (including type checks)
    // Normally is_active() is used, but certain optimisations allow us to avoid that check
    static constexpr auto is_active(particleStore::ParticleReference& p) -> bool { return true; }
    // Returns the open range [start, end) of particles in the node to process for neighbours
    // Certain optimisations (see below) allow us to avoid checking particle types
    static inline auto range_neighbour(treeStore::NodePointer node) { return std::pair{node->lower(), node->upper() + 1}; }
    // Filter function to select neighbors, e.g., particles types of interest (gas, BH, etc)
    // If the filter is done in range(), this can be a no-op (always true)
    static constexpr auto accept_neighbour(particleStore::ParticleReference& p) -> bool { return true; }
    // The particles to be smoothed are added to the list of sinks, but we may need to
    // check again (for example during BH mergers we may have deleted particles)
    static constexpr auto accept_sink(particleStore::ParticleReference& p) -> bool { return true; }
    // The fBall factor to apply during node smoothing
    static constexpr double fBall_factor = 1.0;
    // True if the ball should be updated
    static constexpr bool update_ball = false;
    // When evaluating sinks in ReSmoothNode we need a cutoff based on ball^2 (normally)
    static constexpr auto ball2_sink(particleStore::ParticleReference& p) -> float { return p.ball() * p.ball(); }
    // In resmoothSinks we may need to filter neighbors based on distance and other criteria
    static auto inline accept_ball(PKD pkd, PARTICLE* P, double fDist2, blitz::TinyVector<double, 3> dr_node, const NN& nn) -> bool { return true; }
};

struct SmoothOnlyActive : SmoothEverything {
    static constexpr auto is_active(particleStore::ParticleReference& p) -> bool { return p.is_active(); }
};

// This will select only gas particles as sinks and neighbors but does not check care if they are active.
// As an opmisiation, only the range of particles in the node are checked if this option is enabled.
struct SmoothOnlyGas : SmoothEverything {
#ifdef OPTIM_REORDER_IN_NODES
    static inline auto range(treeStore::NodePointer node) {
        auto start = node->lower();
        auto end = start + node->Ngas();
        return std::pair{start, end};
    }
    static inline auto range_neighbour(treeStore::NodePointer node) { return range(node); }
    static inline auto may_have_active(treeStore::NodePointer& node) -> bool { return node->Ngas() > 0; }
    // The raange() method filters particles, so is_active() will always check only gas particles
    static inline auto is_gas(particleStore::ParticleReference& p) -> bool { return true; }
#else
    static inline auto is_gas(particleStore::ParticleReference& p) -> bool { return p.is_gas(); }
#endif
    static inline auto is_active(particleStore::ParticleReference& p) -> bool { return is_gas(p); }
    static inline auto accept_neighbour(particleStore::ParticleReference& p) -> bool { return p.is_gas(); }
};

// This does the same as SmoothOnlyGas but also checks if the particles are active
struct SmoothOnlyGasActive : SmoothOnlyGas {
    static inline auto is_active(particleStore::ParticleReference& p) -> bool { return is_gas(p) && p.is_active(); }
};

// For iterative density smoothing, the active particles are first marked and we use that flag here
struct SmoothOnlyGasOptimised : SmoothOnlyGasActive {
#ifdef OPTIM_AVOID_IS_ACTIVE
    static inline auto is_active(particleStore::ParticleReference& p) -> bool { return is_gas(p) && p.marked(); }
#endif
};

struct smContext {
    PKD pkd;
    PARTICLE* pSentinel = nullptr;
    std::shared_ptr<SmoothFunction> function;
    std::shared_ptr<particleCACHEhelper> cacheHelper;
    int nSmooth = 0;
    bool bPeriodic = false;
    blitz::TinyVector<double, 3> rLast{0.0, 0.0, 0.0}; /* For the snake */
    std::vector<PQ> vpq;
    /*
    ** Hash table to indicate whether a remote particle is already present in the
    ** priority queue.
    */
    boost::unordered_flat_set<void*> hash;
    std::vector<NN> nnList;
    /*
     ** Also need the stacks for the tree search
     */
    std::vector<uint32_t> S;
    struct stStack {
        int id;
        int iCell;
        double min;
    };
    std::vector<stStack> ST;

    int iSmoothType = 0;
    static auto makeHelper(PKD pkd, SMF* smf, int bSymmetric, int iSmoothType) -> std::shared_ptr<SmoothHelper>;

    void create();
    void finish();

    auto pqSearch(PQ* pq, blitz::TinyVector<double, 3> r, int iRoot) -> PQ*;
    void Gather(double fBall2, blitz::TinyVector<double, 3> r);
    auto getCell(int iCell, int id) { return (id == pkd->Self()) ? pkd->tree[iCell] : pkd->tree[pkd->mdl->Fetch<KDN>(CID_CELL, iCell, id)]; }
    void Release(std::span<NN> list);

    // for node-based smoothing
#ifdef OPTIM_SMOOTH_NODE
    template <class RangePolicy>
    auto getNodeSinks(treeStore::NodePointer node, std::vector<PARTICLE*>& sinks, RangePolicy info = {}) -> Bound;

    template <class RangePolicy>
    void resmoothSinks(std::span<PARTICLE*> sinks, Bound bnd_node, RangePolicy info = {});
    template <class RangePolicy>
    void buildInteractionList(Bound bnd_node, blitz::TinyVector<double, 3> r, RangePolicy info = {});
    template <class RangePolicy>
    void buildInteractionList(Bound bnd_node, RangePolicy info = {});
#endif
  public:
    smContext(PKD pkdIn, int iSmoothType, SMF* smf, bool bSymmetric, int nSmooth, bool bPeriodic);
    smContext(PKD pkdIn, std::shared_ptr<SmoothFunction> helper, int nSmooth, bool bPeriodic);

    template <class T>
        requires std::is_base_of_v<SmoothFunction, T> && std::is_base_of_v<particleCACHEhelper, T>
    smContext(PKD pkdIn, std::shared_ptr<T> helper, int nSmooth, bool bPeriodic)
        : smContext(pkdIn, std::static_pointer_cast<SmoothFunction>(helper), nSmooth, bPeriodic) {
        cacheHelper = std::static_pointer_cast<particleCACHEhelper>(helper);
        helper->initialize_local();
        // Start particle caching space (cell cache is already active).
        pkd->mdl->CacheInitialize(CID_PARTICLE, helper);
    }

    virtual ~smContext();
    auto SmoothSingle(particleStore::ParticleReference& p, int iRoot1 = ROOT, int iRoot2 = 0) -> float;
    auto ReSmoothSingle(particleStore::ParticleReference& p, double fBall) -> bool;

    template <class RangePolicy = ::SmoothEverything>
    void Smooth();
    template <class RangePolicy = ::SmoothEverything>
    auto ReSmooth() -> uint64_t;
#ifdef OPTIM_SMOOTH_NODE
    template <class RangePolicy = ::SmoothEverything>
    auto ReSmoothNode(RangePolicy info = {}) -> uint64_t;
#endif
};

template <class RangePolicy>
void smContext::Smooth() {
    create();
    for (auto& p : pkd->particles) {
        if (RangePolicy::is_active(p)) {
            auto fBall = SmoothSingle(p);
            if constexpr (RangePolicy::update_ball) p.set_ball(fBall);
            pkd->mdl->CacheCheck();
        }
    }
    finish();
}

template <class RangePolicy>
auto smContext::ReSmooth() -> uint64_t {
    uint64_t nSmoothed = 0;
    for (auto& p : pkd->particles) {
        if (RangePolicy::is_active(p) && RangePolicy::accept_sink(p)) {
            ReSmoothSingle(p, p.ball());
            ++nSmoothed;
            pkd->mdl->CacheCheck();
        }
    }
    return nSmoothed;
}

#ifdef OPTIM_SMOOTH_NODE
/* IA: In this version, we loop over the buckets, rather than over the particles.
 *
 * For each bucket, we look for all the surroiding buckets that may interact
 * with any particle in said bucket.
 *
 * Then, we put all those particles (including the own bucket)
 * in an interaction list.
 */
template <class RangePolicy>
auto smContext::ReSmoothNode(RangePolicy info) -> uint64_t {
    uint64_t nSmoothed = 0;
    // Here we store the pointers to the particles whose interaction need to be computed
    std::vector<PARTICLE*> sinks;
    sinks.reserve(256); // At most, the bucket size
    for (auto i = NRESERVED_NODES; i < pkd->Nodes() - 1; ++i) {
        auto node = pkd->tree[i];
        if (node->is_bucket() && RangePolicy::may_have_active(node)) {
            auto bnd_node = getNodeSinks<RangePolicy>(node, sinks, info); // ball factor is 1.0
            if (sinks.empty()) continue;
            buildInteractionList<RangePolicy>(bnd_node, info);
            resmoothSinks<RangePolicy>(sinks, bnd_node, info);
            nSmoothed += sinks.size();
        }
    }
    return nSmoothed;
}

template <class RangePolicy>
auto smContext::getNodeSinks(treeStore::NodePointer node, std::vector<PARTICLE*>& sinks, RangePolicy info) -> Bound {
    sinks.clear();
    sinks.reserve(256); // At most, the bucket size
    auto bnd_node = node->bound();

    blitz::TinyVector<double, 3> fMax_shrink{0.};
    auto [pStart, pEnd] = RangePolicy::range(node);
    for (auto pj = pStart; pj < pEnd; ++pj) {
        auto p = pkd->particles[pj];
        if (RangePolicy::is_active(p)) {
            fMax_shrink = blitz::max(fMax_shrink, abs(p.position() - bnd_node.center()) + p.ball() * RangePolicy::fBall_factor);
            sinks.push_back(&p);
        }
    }
    if (!sinks.empty()) bnd_node.shrink(fMax_shrink);
    return bnd_node;
}

template <class RangePolicy>
void smContext::buildInteractionList(Bound bnd_node, RangePolicy info) {
    auto r = bnd_node.center();
    nnList.clear();
    if (bPeriodic) {
        blitz::TinyVector<int, 3> iStart{floor((r - bnd_node.apothem()) / pkd->fPeriod + 0.5)};
        blitz::TinyVector<int, 3> iEnd{floor((r + bnd_node.apothem()) / pkd->fPeriod + 0.5)};
        for (int ix = iStart[0]; ix <= iEnd[0]; ++ix) {
            r[0] = bnd_node.center(0) - ix * pkd->fPeriod[0];
            for (int iy = iStart[1]; iy <= iEnd[1]; ++iy) {
                r[1] = bnd_node.center(1) - iy * pkd->fPeriod[1];
                for (int iz = iStart[2]; iz <= iEnd[2]; ++iz) {
                    r[2] = bnd_node.center(2) - iz * pkd->fPeriod[2];
                    buildInteractionList<RangePolicy>(bnd_node, r, info);
                }
            }
        }
    }
    else {
        buildInteractionList<RangePolicy>(bnd_node, r, info);
    }
}

template <class RangePolicy>
void smContext::buildInteractionList(Bound bnd_node, blitz::TinyVector<double, 3> r, RangePolicy) {
    auto mdl = pkd->mdl;
    int id, iCell;

    // We look for the biggest node that encloses the needed domain
    // We can only take advantage of this if we are are in the original cell
    auto kdn = getCell(iCell = pkd->iTopTree[ROOT], id = pkd->Self());

    //  Now we start the walk as usual
    ST.clear();
    while (1) {
        auto bnd = kdn->bound();
        if (any(abs(bnd.center() - r) - bnd.apothem() - bnd_node.apothem() > 0.)) goto NoIntersect;

        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            int idUpper, iUpper;
            std::tie(iCell, id, iUpper, idUpper) = kdn->get_child_cells(id);
            kdn = getCell(iCell, id);
            ST.push_back({.id = idUpper, .iCell = iUpper});
            continue;
        }
        else {
            if (id == pkd->Self()) {
                auto [pStart, pEnd] = RangePolicy::range_neighbour(kdn);
                for (auto pj = pStart; pj < pEnd; ++pj) {
                    auto p = pkd->particles[pj];
                    if (RangePolicy::accept_neighbour(p)) {
                        const blitz::TinyVector<double, 3> dr{r - p.position()};
                        if (all(abs(dr) <= bnd_node.apothem())) {
                            nnList.push_back(
                                {.pPart = &p, .fDist2 = dot(dr, dr), .dr = dr, .fBall = p.ball(), .bMarked = p.marked(), .iIndex = pj, .iPid = id});
                        }
                    }
                }
            }
            else {
                auto [pStart, pEnd] = RangePolicy::range_neighbour(kdn);
                for (auto pj = pStart; pj < pEnd; ++pj) {
                    auto p = pkd->particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
                    if (RangePolicy::accept_neighbour(p)) {
                        const blitz::TinyVector<double, 3> dr{r - p.position()};
                        if (all(abs(dr) <= bnd_node.apothem())) {
                            // The actual particle information fill be acquired later
                            nnList.push_back(
                                {.pPart = nullptr, .fDist2 = dot(dr, dr), .dr = dr, .fBall = p.ball(), .bMarked = p.marked(), .iIndex = pj, .iPid = id});
                        }
                    }
                }
            }
        }
    NoIntersect:
        if (!ST.empty()) {
            auto& top = ST.back();
            id = top.id;
            iCell = top.iCell;
            kdn = getCell(iCell, id);
            ST.pop_back();
        }
        else break;
    }
}

// IA: Now we should have inside nnList all the particles in the
//  bucket (sinks) and those of which can interact with them
//  from other buckets (smx->nnList)
//
// We just have to proceed to compute the correct dx, dy, dz and
// pass that nnList to the smoothfcn routine
//
// However, we have different options to do so:
// 1) Naive: we pass the whole nnList
//
// 2) Sorting: we could follow Gonnet, 2007 (10.1002/jcc.20563) to
// reduce the number of distance computations, but this is troublesome
// in our case because:
//      a) we need to compute the distance anyway for sorting
//      b) we could sort relative to the cell, but this is suboptimal
//      c) we are not computing cell-cell interactions, so there is
//            no well-defined axis that could be used for projection
//
template <class RangePolicy>
void smContext::resmoothSinks(std::span<PARTICLE*> sinks, Bound bnd_node, RangePolicy) {
    std::vector<NN> nnList_p;
    nnList_p.reserve(nnList.size());

    for (auto& P : sinks) {
        auto partj = pkd->particles[P];
        if (!RangePolicy::accept_sink(partj)) continue;
        const float pBall2 = RangePolicy::ball2_sink(partj);
        const blitz::TinyVector<double, 3> dr_node{bnd_node.center() - partj.position()};

        nnList_p.clear();
        for (auto& nn : nnList) {
            if (P == nn.pPart) continue;
            const blitz::TinyVector<double, 3> dr{nn.dr - dr_node};
            const double fDist2 = dot(dr, dr);
            if (fDist2 < pBall2 && RangePolicy::accept_ball(pkd, P, fDist2, dr_node, nn)) {
                nnList_p.push_back({.pPart = (nn.iPid == pkd->Self()) ? nn.pPart : pkd->mdl->Acquire<PARTICLE>(CID_PARTICLE, nn.iIndex, nn.iPid),
                                    .fDist2 = fDist2,
                                    .dr = dr,
                                    .iIndex = nn.iIndex,
                                    .iPid = nn.iPid});
            }
        }
        function->Smooth(partj, partj.ball(), nnList_p);
        Release(nnList_p);
    }
}

#endif // OPTIM_SMOOTH_NODE
using SMX = smContext*;
void smSmooth(SMX, SMF*);
auto smReSmooth(SMX, int) -> uint64_t;

#endif
