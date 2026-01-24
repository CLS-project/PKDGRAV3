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

#include "pkd_config.h"
#include "smooth/smoothfcn.h"
#include <array>
#include <cmath>

#include <new>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "hydro/hydro.h"
#include "starformation/feedback.h"
#include "blackhole/merger.h"
#include "blackhole/evolve.h"
#include "blackhole/step.h"
#include "blackhole/drift.h"
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#include <sys/stat.h>

using blitz::abs;
using blitz::all;
using blitz::any;
using blitz::dot;
using blitz::floor;
using blitz::TinyVector;

// This delegating constructor will go away once we convert the functions
smContext::smContext(PKD pkdIn, int iSmoothType, SMF* smf, bool bSymmetric, int nSmooth, bool bPeriodic, bool bMakeCache)
    : smContext(pkdIn, makeHelper(pkdIn, smf, bSymmetric, iSmoothType), nSmooth, bPeriodic, bMakeCache) {
    this->iSmoothType = iSmoothType;
}

smContext::smContext(PKD pkdIn, std::shared_ptr<SmoothHelper> helper, int nSmooth, bool bPeriodic, bool bMakeCache)
    : pkd(pkdIn), function(helper), nSmooth(nSmooth), bPeriodic(bPeriodic), iSmoothType(-1) {
    pSentinel = static_cast<PARTICLE*>(malloc(pkd->particles.ParticleSize()));

    helper->initialize_local();
    /*
    ** Start particle caching space (cell cache is already active).
    */
    if (bMakeCache) {
        bOwnCache = 1;
        pkd->mdl->CacheInitialize(CID_PARTICLE, helper);
    }
    else bOwnCache = 0;
    /*
    ** Allocate Nearest-Neighbor List.
    */
    nnList.reserve(NNLIST_INCREMENT);
    /*
    ** Allocate priority queue.
    */
    vpq.resize(nSmooth);
    PQ_INIT(vpq.data(), nSmooth);
    // Setup hash table.
    hash.max_load_factor(0.7f);
    hash.reserve(nSmooth); // expected size

    // Allocate special stacks for searching within the tree.
    // 1024 is more than enough.
    S.reserve(1024);
    ST.reserve(1024);
    /*
    ** Set up the sentinel particle with some very far away distance.
    ** This is used to initially load the priority queue and all references
    ** to this particle should end up being replaced in the priority queue
    ** as long as there are nSmooth particles.
    */
    auto sentinel = pkd->particles[pSentinel];
    if (pkd->bIntegerPosition) sentinel.raw_position<int32_t>() = TinyVector<int32_t, 3>(INT32_MAX, INT32_MAX, INT32_MAX);
    else sentinel.raw_position<double>() = TinyVector<double, 3>(HUGE_VAL, HUGE_VAL, HUGE_VAL);
}

smContext::~smContext() {

    /*
    ** Stop particle caching space.
    */
    if (bOwnCache) pkd->mdl->FinishCache(CID_PARTICLE);
    free(pSentinel);
}

// This function should go away once we convert the functions
auto smContext::makeHelper(PKD pkd, SMF* smf, int bSymmetric, int iSmoothType) -> std::shared_ptr<SmoothHelper> {
    switch (iSmoothType) {
    case SMX_NULL:
        return std::make_shared<SmoothHelperNull>(pkd);
    case SMX_DENSITY:
        if (bSymmetric) return std::make_shared<SmoothHelperDensitySymmetric>(pkd);
        else return std::make_shared<SmoothHelperDensity>(pkd);
    case SMX_DENSITY_F1:
        assert(!bSymmetric);
        return std::make_shared<SmoothHelperDensityF1>(pkd);
    case SMX_DENSITY_M3:
        assert(!bSymmetric);
        assert(pkd->particles.present(PKD_FIELD::oGroup));
        return std::make_shared<SmoothHelperDensityM3>(pkd);
    case SMX_GRADIENT_M3:
        assert(!bSymmetric);
        return std::make_shared<SmoothHelperGradientM3>(pkd);
    case SMX_HOP_LINK:
        assert(!bSymmetric);
        return std::make_shared<SmoothHelperHopChains>(pkd);
    case SMX_PRINTNN:
        return std::make_shared<SmoothHelperPrintNN>(pkd);
    case SMX_HYDRO_DENSITY:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroDensity>(pkd, smf->dNeighborsStd, smf->bIterativeSmoothingLength);
    case SMX_HYDRO_DENSITY_FINAL:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroDensityFinal>(pkd);
    case SMX_HYDRO_GRADIENT:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroGradients>(pkd, smf->bConditionedBarthJespersen, smf->bStricterSlopeLimiter);
#ifdef OPTIM_FLUX_VEC
    case SMX_HYDRO_FLUX_VEC:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroFluxes>(pkd, smf->dTime, smf->dDelta, smf->a, smf->H, smf->dConstGamma, smf->eEOS);
#endif
    case SMX_HYDRO_STEP:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroStep>(pkd, smf->dEtaCourant, smf->a, smf->dCFLacc, smf->dDelta, smf->bDoGravity);
    case SMX_SN_FEEDBACK:
        return std::make_shared<SmoothHelperSNFeedback>(pkd, smf->dTime, smf->dCCSNFBDelay, smf->dSNIaFBDelay, smf->dSNFBDu, smf->dCCSNFBSpecEnergy,
                                                        smf->dSNIaFBSpecEnergy, smf->dConstGamma);
    case SMX_BH_MERGER:
        return std::make_shared<SmoothHelperBHmerger>(pkd);
    case SMX_BH_STEP:
        return std::make_shared<SmoothHelperBHstep>(pkd);
    case SMX_BH_EVOLVE:
        return std::make_shared<SmoothHelperBHEvolve>(pkd, smf->a, smf->dBHAccretionAlpha, smf->dBHAccretionCvisc, smf->dBHAccretionEddFac, smf->dConstGamma,
                                                      smf->dBHFBEff, smf->dBHFBEcrit, smf->bBHAccretion, smf->bBHFeedback);
    case SMX_BH_GASPIN:
        return std::make_shared<SmoothHelperBHGasPin>(pkd, smf->a);
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        return std::make_shared<SmoothHelperChemEnrich>(pkd, smf->dTime, smf->dWindSpecificEkin, smf->dCCSNMinMass, smf->dSNIaNorm, smf->dSNIaScale);
#endif
    case SMX_BALL:
        assert(pkd->particles.present(PKD_FIELD::oBall));
        return std::make_shared<SmoothHelperBallSmooth>(pkd);
    }
    return nullptr;
}

static auto getCell(PKD pkd, int iCell, int id) {
    return (id == pkd->Self()) ? pkd->tree[iCell] : pkd->tree[pkd->mdl->Fetch<KDN>(CID_CELL, iCell, id)];
}

auto smContext::pqSearch(PQ* pq, TinyVector<double, 3> r, int iRoot) -> PQ* {
    auto mdl = pkd->mdl;
    int idSelf = pkd->Self();
    int iCell, id;
    double fDist2;
    TinyVector<double, 3> dr;

    /* Start at the root node of the tree */
    auto kdn = getCell(pkd, pkd->iTopTree[iRoot], id = idSelf);
    ST.clear();
    while (1) {
        while (kdn->is_cell()) {
            auto [iLower, idLower, iUpper, idUpper] = kdn->get_child_cells(id);
            kdn = getCell(pkd, iLower, idLower);
            auto min1 = kdn->bound().mindist(r);
            kdn = getCell(pkd, iUpper, idUpper);
            auto min2 = kdn->bound().mindist(r);
            if (min1 < min2) {
                if (min1 > pq->fDist2) goto NoIntersect;
                ST.push_back({.id = idUpper, .iCell = iUpper, .min = min2});
                id = idLower;
                iCell = iLower;
                kdn = getCell(pkd, iCell, id);
            }
            else {
                if (min2 > pq->fDist2) goto NoIntersect;
                ST.push_back({.id = idLower, .iCell = iLower, .min = min1});
                id = idUpper;
                iCell = iUpper;
            }
        }
        /* Now at a bucket */
        if (id == idSelf) {
#ifdef OPTIM_REORDER_IN_NODES
            auto pEnd = (function->gasOnly()) ? kdn->lower() + kdn->Ngas() : kdn->upper() + 1;
#else
            auto pEnd = kdn->upper() + 1;
#endif
            for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                auto p = pkd->particles[pj];
                if (!p.marked()) continue;
#ifndef OPTIM_REORDER_IN_NODES
                if (function->gasOnly() && !p.is_gas()) continue;
#endif
                dr = r - p.position();
                fDist2 = dot(dr, dr);
                if (fDist2 <= pq->fDist2) {
                    if (pq->iPid == idSelf) {
                        pkd->particles[pq->iIndex].set_marked(true);
                    }
                    else {
                        hash.erase(pq->pPart);
                        mdl->Release(CID_PARTICLE, pq->pPart);
                        pq->iPid = idSelf;
                    }
                    pq->pPart = &p;
                    pq->fDist2 = fDist2;
                    pq->dr = dr;
                    pq->iIndex = pj;
                    p.set_marked(false); /* de-activate a particle that enters the queue */
                    PQ_REPLACE(pq);
                }
            }
        }
        else {
#ifdef OPTIM_REORDER_IN_NODES
            auto pEnd = (function->gasOnly()) ? kdn->lower() + kdn->Ngas() : kdn->upper() + 1;
#else
            auto pEnd = kdn->upper() + 1;
#endif
            for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                auto p = pkd->particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
#ifndef OPTIM_REORDER_IN_NODES
                if (function->gasOnly() && !p.is_gas()) continue;
#endif
                if (hash.contains(&p)) continue;
                dr = r - p.position();
                fDist2 = dot(dr, dr);
                if (fDist2 <= pq->fDist2) {
                    if (pq->iPid == idSelf) {
                        pkd->particles[pq->iIndex].set_marked(true);
                    }
                    else {
                        hash.erase(pq->pPart);
                        mdl->Release(CID_PARTICLE, pq->pPart);
                    }
                    auto p = pkd->particles[mdl->Acquire<PARTICLE>(CID_PARTICLE, pj, id)];
                    pq->pPart = &p;
                    pq->fDist2 = fDist2;
                    pq->dr = dr;
                    pq->iIndex = pj;
                    pq->iPid = id;
                    hash.insert(&p);
                    PQ_REPLACE(pq);
                }
            }
        }

    NoIntersect:
        while (!ST.empty() && ST.back().min > pq->fDist2) ST.pop_back();
        if (!ST.empty()) {
            auto& top = ST.back();
            id = top.id;
            iCell = top.iCell;
            kdn = getCell(pkd, iCell, id);
            ST.pop_back();
        }
        else return pq;
    }
}

void smContext::initialize() {
    // Initialize the priority queue first.
    auto sentinel = pkd->particles[pSentinel];
    for (auto& q : vpq) {
        q.pPart = pSentinel;
        q.iIndex = pkd->Local();
        q.iPid = pkd->Self();
        q.dr = sentinel.position();
        q.fDist2 = dot(q.dr, q.dr);
    }
    rLast = 0.0;
}

void smContext::finish() {
    // Release acquired pointers and source-reactivate particles in prioq.
    for (auto& q : vpq) {
        if (q.iPid == pkd->Self()) {
            pkd->particles[q.iIndex].set_marked(true);
        }
        else {
            hash.erase(q.pPart);
            pkd->mdl->Release(CID_PARTICLE, q.pPart);
        }
    }
}

auto smContext::SmoothSingle(particleStore::ParticleReference& p, int iRoot1, int iRoot2) -> float {
    if (!function->shouldSmooth(p)) return 0.0f; // No smoothing needed for this particle
    int ix, iy, iz;
    TinyVector<double, 3> r;
    double fBall;
    int j;
    PQ* pq;

    auto p_r = p.position();

    /*
    ** Correct distances and rebuild priority queue.
    */
    if (bPeriodic) {
        for (j = 0; j < 3; ++j) {
            if (p_r[j] > rLast[j] + 0.5 * pkd->fPeriod[j]) rLast[j] += pkd->fPeriod[j];
            else if (p_r[j] < rLast[j] - 0.5 * pkd->fPeriod[j]) rLast[j] -= pkd->fPeriod[j];
        }
    }

    for (auto& q : vpq) {
        q.dr += p_r - rLast;
        q.fDist2 = dot(q.dr, q.dr);
    }
    r = p_r;
    rLast = p_r;
    PQ_BUILD(vpq.data(), vpq.size(), pq);
    pq = pqSearch(pq, r, iRoot1);
    if (iRoot2) pq = pqSearch(pq, r, iRoot2);
    /*
    ** Search in replica boxes if it is required.
    */
    //    printf("fPeriod %f %f %f \n",pkd->fPeriod[0],pkd->fPeriod[1],pkd->fPeriod[2] );
    //    printf("x %f y %f z %f \n", p_r[0], p_r[1], p_r[2]);
    if (bPeriodic) {
        fBall = sqrt(pq->fDist2);
        TinyVector<int, 3> iStart = floor((p_r - fBall) / pkd->fPeriod + 0.5);
        TinyVector<int, 3> iEnd = floor((p_r + fBall) / pkd->fPeriod + 0.5);
        for (ix = iStart[0]; ix <= iEnd[0]; ++ix) {
            r[0] = p_r[0] - ix * pkd->fPeriod[0];
            for (iy = iStart[1]; iy <= iEnd[1]; ++iy) {
                r[1] = p_r[1] - iy * pkd->fPeriod[1];
                for (iz = iStart[2]; iz <= iEnd[2]; ++iz) {
                    r[2] = p_r[2] - iz * pkd->fPeriod[2];
                    if (ix || iy || iz) {
                        pq = pqSearch(pq, r, iRoot1);
                        if (iRoot2) pq = pqSearch(pq, r, iRoot2);
                    }
                }
            }
        }
    }
    fBall = sqrt(pq->fDist2);

    /* IA: I do not fully understand this fBall, so I will compute my own kernel length such that it encloses
     * all the particles in the neighbor list. This means that h > 0.5*max(dist). I have taken 0.501 as a safe
     * value, because 0.5 would exclude the furthest particle(s) */
    //    int i;
    //    fBall = 0.0;
    //    for (i=0; i<nSmooth; ++i){
    //       if (fBall < vpq[i].fDist2) fBall = vpq[i].fDist2;
    //    }
    //    fBall = 0.50*sqrt(fBall);

    /*
    ** Apply smooth funtion to the neighbor list.
    */
    function->Smooth(p, fBall, vpq);
    return fBall;
}

void smContext::Smooth() {
    initialize();
    SmoothParticles();
    finish();
}

void smContext::SmoothParticles() {
    for (auto& p : pkd->particles) {
        SmoothSingle(p);
        pkd->mdl->CacheCheck();
    }
}

void smSmooth(SMX smx, SMF* smf) {
    PKD pkd = smx->pkd;
    float fBall;

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (auto& p : pkd->particles) p.set_marked(true);
    smx->initialize();
    switch (smx->iSmoothType) {
    case SMX_BH_EVOLVE:
    case SMX_BH_STEP:
    case SMX_BH_GASPIN:
        for (auto& p : pkd->particles) {
            if (p.is_bh()) {
                smx->SmoothSingle(p);
            }
        }
        break;
    case SMX_SN_FEEDBACK:
        for (auto& p : pkd->particles) {
            if (p.is_star()) {
                auto& star = p.star();

                if (!star.bCCSNFBDone && ((smf->dTime - star.fTimer) > smf->dCCSNFBDelay)) {
                    smx->SmoothSingle(p);
                }

                if (!star.bSNIaFBDone && ((smf->dTime - star.fTimer) > smf->dSNIaFBDelay)) {
                    smx->SmoothSingle(p);
                }
            }
        }
        break;
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        for (auto& p : pkd->particles) {
            if (p.is_star()) {
                auto& star = p.star();
                if (smf->dTime > star.fNextEnrichTime) {
                    smx->SmoothSingle(p);
                }
            }
        }
        break;
#endif
    default:
        for (auto& p : pkd->particles) {
            if (!smf->bMeshlessHydro) {
                smx->SmoothSingle(p);
                // p.set_ball(smx->SmoothSingle(p));
            }
            else {
                if (p.is_active()) {
                    fBall = smx->SmoothSingle(p);
                    if (smf->bUpdateBall) {
                        p.set_ball(fBall);
                    }
                }
            }
            /*
            ** Call mdlCacheCheck to make sure we are making progress!
            */
            pkd->mdl->CacheCheck();
        }
    }
    smx->finish();
}

void smContext::Gather(double fBall2, TinyVector<double, 3> r) {
    auto mdl = pkd->mdl;
    int idSelf = pkd->Self();
    int iCell, id;
    int pj, pEnd;

    auto kdn = getCell(pkd, iCell = pkd->iTopTree[ROOT], id = idSelf);

    ST.clear();
    while (1) {
        auto min2 = kdn->bound().mindist(r);
        if (min2 > fBall2) {
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            int idUpper, iUpper;
            std::tie(iCell, id, iUpper, idUpper) = kdn->get_child_cells(id);
            kdn = getCell(pkd, iCell, id);
            ST.push_back({.id = idUpper, .iCell = iUpper, .min = 0.0});
            continue;
        }
        else {
            if (id == pkd->Self()) {
                pEnd = kdn->upper();
                for (pj = kdn->lower(); pj <= pEnd; ++pj) {
                    auto p = pkd->particles[pj];
                    if (!p.is_gas()) continue;
                    TinyVector<double, 3> dr{r - p.position()};
                    double fDist2 = dot(dr, dr);
                    if (fDist2 <= fBall2) {
                        nnList.push_back({.pPart = &p, .fDist2 = fDist2, .dr = dr, .iIndex = pj, .iPid = idSelf});
                    }
                }
            }
            else {
                pEnd = kdn->upper();
                for (pj = kdn->lower(); pj <= pEnd; ++pj) {
                    auto p = pkd->particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
                    if (!p.is_gas()) continue;
                    TinyVector<double, 3> dr{r - p.position()};
                    double fDist2 = dot(dr, dr);
                    if (fDist2 <= fBall2) {
                        nnList.push_back({.pPart = mdl->Acquire<PARTICLE>(CID_PARTICLE, pj, id), .fDist2 = fDist2, .dr = dr, .iIndex = pj, .iPid = id});
                    }
                }
            }
        }
    NoIntersect:
        if (!ST.empty()) {
            auto& top = ST.back();
            id = top.id;
            iCell = top.iCell;
            kdn = getCell(pkd, iCell, id);
            ST.pop_back();
        }
        else break;
    }
}

void smDoGatherLocal(SMX smx, double fBall2, TinyVector<double, 3> r, void (*Do)(SMX, PARTICLE*, double)) {
    PKD pkd = smx->pkd;
    int iCell;

    std::array<int, 1024> S;
    size_t sp = 0;
    auto kdn = pkd->tree[iCell = ROOT];
    while (1) {
        auto min2 = kdn->bound().mindist(r);
        if (min2 > fBall2) {
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            S[sp++] = kdn->rchild();
            kdn = pkd->tree[iCell = kdn->lchild()];
            continue;
        }
        else {
            for (auto& p : *kdn) {
                TinyVector<double, 3> dr{r - p.position()};
                double fDist2 = dot(dr, dr);
                if (fDist2 <= fBall2) {
                    Do(smx, &p, fDist2);
                }
            }
        }
    NoIntersect:
        if (sp) {
            kdn = pkd->tree[iCell = S[--sp]];
        }
        else break;
    }
}

auto smContext::ReSmoothSingle(particleStore::ParticleReference& p, double fBall) -> bool {
    if (!function->shouldSmooth(p)) return false; // No smoothing needed for this particle
    TinyVector<double, 3> r;

    auto R = p.position();

    nnList.clear();
    /*
    ** Note for implementing SLIDING PATCH, the offsets for particles are
    ** negative here, reflecting the relative +ve offset of the simulation
    ** volume.
    */
    if (bPeriodic) {
        TinyVector<int, 3> iStart = floor((R - fBall) / pkd->fPeriod + 0.5);
        TinyVector<int, 3> iEnd = floor((R + fBall) / pkd->fPeriod + 0.5);
        for (int ix = iStart[0]; ix <= iEnd[0]; ++ix) {
            r[0] = R[0] - ix * pkd->fPeriod[0];
            for (int iy = iStart[1]; iy <= iEnd[1]; ++iy) {
                r[1] = R[1] - iy * pkd->fPeriod[1];
                for (int iz = iStart[2]; iz <= iEnd[2]; ++iz) {
                    r[2] = R[2] - iz * pkd->fPeriod[2];
                    Gather(fBall * fBall, r);
                }
            }
        }
    }
    else {
        Gather(fBall * fBall, R);
    }
    /*
    ** Apply smooth funtion to the neighbor list.
    */
    function->Smooth(p, fBall, nnList);
    /*
    ** Release acquired pointers.
    */
    for (auto& nn : nnList) {
        if (nn.iPid != pkd->Self()) {
            pkd->mdl->Release(CID_PARTICLE, nn.pPart);
        }
    }
    return true;
}

auto smContext::ReSmoothParticles() -> uint64_t {
    uint64_t nSmoothed = 0;
    for (auto& p : pkd->particles) {
        ReSmoothSingle(p, p.ball());
        ++nSmoothed;
        pkd->mdl->CacheCheck();
    }
    return nSmoothed;
}

uint64_t smReSmooth(SMX smx, int iSmoothType) {
    PKD pkd = smx->pkd;
    uint64_t nSmoothed = 0;

    switch (iSmoothType) {
    case SMX_HYDRO_DENSITY:
        for (auto& p : pkd->particles) {
            if (p.is_active() && p.marked() && p.is_gas()) {
                smx->ReSmoothSingle(p, p.ball());
                nSmoothed++;
            }
        }
        break;

    case SMX_HYDRO_STEP:
        for (auto& p : pkd->particles) {
            if (p.is_gas()) {
                if (p.is_active()) {
                    smx->ReSmoothSingle(p, p.ball());
                    nSmoothed++;
                }
            }
        }
        break;

    case SMX_SN_FEEDBACK:
    case SMX_CHEM_ENRICHMENT:
    case SMX_BH_EVOLVE:
    case SMX_BH_STEP:
    case SMX_BH_GASPIN:
        assert(0);
        break;

    default:
        for (auto& p : pkd->particles) {
            if (p.is_active() && p.is_gas()) {
                smx->ReSmoothSingle(p, p.ball());
                nSmoothed++;
            }
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
uint64_t smReSmoothNode(SMX smx, SMF* smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    auto mdl = pkd->mdl;
    uint64_t nSmoothed = 0;

    smx->nnList.clear();

    std::vector<NN> nnList_p;
    nnList_p.reserve(NNLIST_INCREMENT);

    // Here we store the pointers to the particles whose interaction
    // need to be computed
    std::vector<PARTICLE*> sinks;
    sinks.reserve(smf->nBucket); // At most, the bucket size

    /* For allowing vectorization, it is better to use an structure of
     *  arrays rather than an array of structures.
     *
     *  In our case, the structure is just an array of pointers to the locations
     *  in the buffer where a given array of variables starts.
     *
     *  Having everything in the same buffer can be advantageous as it should
     *  be in cache, but it was painful to code...
     */

    double fBall_factor = 1.;
    if (iSmoothType == SMX_HYDRO_DENSITY) fBall_factor *= 1.2; // An small margin is kept in case fBall needs to increase

    for (auto i = NRESERVED_NODES; i < pkd->Nodes() - 1; ++i) {
        auto node = pkd->tree[i];
        if (node->is_bucket()) {

            // Prepare the interaction list

            auto bnd_node = node->bound();

            TinyVector<double, 3> r{bnd_node.center()};

            // First, we add all the particles whose interactions need to be computed
            sinks.clear();
    #ifdef OPTIM_REORDER_IN_NODES
            auto pEnd = node->lower() + node->Ngas();
    #else // OPTIM_REORDER_IN_NODES
            auto pEnd = node->upper() + 1;
    #endif

            TinyVector<double, 3> fMax_shrink{0.};
            for (auto pj = node->lower(); pj < pEnd; ++pj) {
                auto p = pkd->particles[pj];

    #ifdef OPTIM_AVOID_IS_ACTIVE
                int pIsActive = p.marked();
    #else
                int pIsActive = p.is_active();
    #endif

                if (pIsActive) {
                    if ((iSmoothType == SMX_HYDRO_DENSITY) || (iSmoothType == SMX_HYDRO_DENSITY_FINAL)) {

    #ifndef OPTIM_AVOID_IS_ACTIVE
                        if (!p.marked()) continue;
    #endif

                        // If there is only feedback, once the particle explodes
                        // there is no need to updated its fBall.
                        // However, if we have stellar evolution, fBall needs to be
                        // always updated.
                        if (!pkd->StelEvolData)
                            if (p.is_star() && p.star().bCCSNFBDone && p.star().bSNIaFBDone) continue;

    #ifndef OPTIM_REORDER_IN_NODES
                        // Explicit check of the types
                        if (!p.is_gas() && !p.is_star()) continue;
    #endif
                    }
                    else {
    #ifndef OPTIM_REORDER_IN_NODES
                        // Explicit check of the types
                        if (!p.is_gas()) continue;
    #endif
                    } // SMX_HYDRO_DENSITY

                    fMax_shrink = blitz::max(fMax_shrink, abs(p.position() - bnd_node.center()) + p.ball() * fBall_factor);

                    sinks.push_back(&p);
                }
            }
            // There are no elligible particles in this bucket, go to the next
            if (sinks.empty()) continue;

            smx->nnList.clear();

            bnd_node.shrink(fMax_shrink);

            if (smx->bPeriodic) {
                TinyVector<int, 3> iStart{floor((r - bnd_node.apothem()) / pkd->fPeriod + 0.5)};
                TinyVector<int, 3> iEnd{floor((r + bnd_node.apothem()) / pkd->fPeriod + 0.5)};
                for (int ix = iStart[0]; ix <= iEnd[0]; ++ix) {
                    r[0] = bnd_node.center(0) - ix * pkd->fPeriod[0];
                    for (int iy = iStart[1]; iy <= iEnd[1]; ++iy) {
                        r[1] = bnd_node.center(1) - iy * pkd->fPeriod[1];
                        for (int iz = iStart[2]; iz <= iEnd[2]; ++iz) {
                            r[2] = bnd_node.center(2) - iz * pkd->fPeriod[2];
                            buildInteractionList(smx, node, bnd_node, r, ix, iy, iz);
                        }
                    }
                }
            }
            else {
                buildInteractionList(smx, node, bnd_node, r, 0, 0, 0);
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

            // For the smoothing length determination we can bypass the typical
            //  flow of calling fcnsmooth, as probably we have gathered more
            //  neighbours than needed and thus the iterative procedure should be
            //  faster
            if (iSmoothType == SMX_HYDRO_DENSITY) {
                hydroDensity_node(pkd, smf, bnd_node, sinks, smx->nnList.data(), smx->nnList.size());
            }
            else {
                for (auto& P : sinks) {
                    assert(P);
                    auto partj = pkd->particles[P];
                    const float pBall2 = partj.ball() * partj.ball();
                    const TinyVector<double, 3> dr_node{bnd_node.center() - partj.position()};

                    nnList_p.clear();
                    for (auto& nn : smx->nnList) {
                        if (P == nn.pPart) continue;
                        const TinyVector<double, 3> dr{nn.dr - dr_node};
                        const double fDist2 = dot(dr, dr);
                        if (fDist2 < pBall2) {

                            // Reasons not to compute this interaction
                            if (iSmoothType == SMX_HYDRO_FLUX || iSmoothType == SMX_HYDRO_FLUX_VEC) {

                                const auto& qBall = nn.fBall;
                                if (qBall * qBall < fDist2) continue;

    #ifdef OPTIM_AVOID_IS_ACTIVE
                                int qIsActive = nn.bMarked;
    #else
                                abort(); // not supported
    #endif

    #ifdef OPTIM_NO_REDUNDANT_FLUXES
                                if (qIsActive) {
                                    if (pkd->Self() == nn.iPid) {
                                        if (P < nn.pPart) continue;
                                    }
                                    else if (pkd->Self() < nn.iPid) continue;
                                }
    #endif
                            }

                            nnList_p.push_back({.pPart = (nn.iPid == pkd->Self()) ? nn.pPart : mdl->Acquire<PARTICLE>(CID_PARTICLE, nn.iIndex, nn.iPid),
                                                .fDist2 = fDist2,
                                                .dr = dr,
                                                .iIndex = nn.iIndex,
                                                .iPid = nn.iPid});
                        }
                    }
                    smx->function->Smooth(partj, partj.ball(), nnList_p);
                    for (auto& nn : nnList_p) {
                        if (nn.iPid != pkd->Self()) {
                            mdl->Release(CID_PARTICLE, nn.pPart);
                        }
                    }
                }
            }

            nSmoothed += sinks.size();
        }
    }
    return nSmoothed;
}

void buildInteractionList(SMX smx, KDN* node, Bound bnd_node, TinyVector<double, 3> r, int ix, int iy, int iz) {
    PKD pkd = smx->pkd;
    auto mdl = pkd->mdl;
    int id, iCell;

    // We look for the biggest node that encloses the needed domain
    // We can only take advantage of this if we are are in the original cell
    auto kdn = getCell(pkd, iCell = pkd->iTopTree[ROOT], id = pkd->Self());

    //  Now we start the walk as usual
    smx->ST.clear();
    while (1) {
        auto bnd = kdn->bound();
        if (any(abs(bnd.center() - r) - bnd.apothem() - bnd_node.apothem() > 0.)) goto NoIntersect;

        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            int idUpper, iUpper;
            std::tie(iCell, id, iUpper, idUpper) = kdn->get_child_cells(id);
            kdn = getCell(pkd, iCell, id);
            smx->ST.push_back({.id = idUpper, .iCell = iUpper});
            continue;
        }
        else {
            if (id == pkd->Self()) {
    #ifdef OPTIM_REORDER_IN_NODES
                auto pEnd = kdn->lower() + kdn->Ngas();
    #else
                auto pEnd = kdn->upper() + 1;
    #endif
                // printf("pEnd %d \n", pEnd);
                for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                    auto p = pkd->particles[pj];
    #ifndef OPTIM_REORDER_IN_NODES
                    if (!p.is_gas()) continue;
    #endif
                    const TinyVector<double, 3> dr{r - p.position()};
                    if (all(abs(dr) <= bnd_node.apothem())) {
                        smx->nnList.push_back(
                            {.pPart = &p, .fDist2 = dot(dr, dr), .dr = dr, .fBall = p.ball(), .bMarked = p.marked(), .iIndex = pj, .iPid = id});
                    }
                }
            }
            else {
    #ifdef OPTIM_REORDER_IN_NODES
                auto pEnd = kdn->lower() + kdn->Ngas();
    #else
                auto pEnd = kdn->upper() + 1;
    #endif
                for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                    auto p = pkd->particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
    #ifndef OPTIM_REORDER_IN_NODES
                    if (!p.is_gas()) continue;
    #endif
                    const TinyVector<double, 3> dr{r - p.position()};
                    if (all(abs(dr) <= bnd_node.apothem())) {
                        // The actual particle information fill be acquired later
                        smx->nnList.push_back(
                            {.pPart = nullptr, .fDist2 = dot(dr, dr), .dr = dr, .fBall = p.ball(), .bMarked = p.marked(), .iIndex = pj, .iPid = id});
                    }
                }
            }
        }
    NoIntersect:
        if (!smx->ST.empty()) {
            auto& top = smx->ST.back();
            id = top.id;
            iCell = top.iCell;
            kdn = getCell(pkd, iCell, id);
            smx->ST.pop_back();
        }
        else break;
    }
}
#endif // OPTIM_SMOOTH_NODE
