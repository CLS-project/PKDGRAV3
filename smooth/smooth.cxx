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

#include "core/bound.h"
#include "pkd_config.h"
#include "smooth/smoothfcn.h"
#include <array>
#include <cmath>

#include <new>
#include <utility>
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
smContext::smContext(PKD pkdIn, int iSmoothType, SMF* smf, bool bSymmetric, int nSmooth, bool bPeriodic)
    : smContext(pkdIn, makeHelper(pkdIn, smf, bSymmetric, iSmoothType), nSmooth, bPeriodic) {
    this->iSmoothType = iSmoothType;
}

smContext::smContext(PKD pkdIn, std::shared_ptr<SmoothFunction> helper, int nSmooth, bool bPeriodic)
    : pkd(pkdIn), function(std::move(helper)), nSmooth(nSmooth), bPeriodic(bPeriodic), iSmoothType(-1) {
    pSentinel = static_cast<PARTICLE*>(malloc(pkd->particles.ParticleSize()));
    /*
    ** Allocate Nearest-Neighbor List.
    */
    nnList.reserve(nSmooth);
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
    // Stop particle caching space, but only if it was created by this class.
    if (cacheHelper) pkd->mdl->FinishCache(CID_PARTICLE);
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
    case SMX_HYDRO_GRADIENT:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroGradients>(pkd, smf->bConditionedBarthJespersen, smf->bStricterSlopeLimiter);
    case SMX_HYDRO_STEP:
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
        return std::make_shared<SmoothHelperHydroStep>(pkd, smf->dEtaCourant, smf->a, smf->dCFLacc, smf->dDelta, smf->bDoGravity);
    case SMX_SN_FEEDBACK:
        return std::make_shared<SmoothHelperSNFeedback>(pkd, smf->dTime, smf->dCCSNFBDelay, smf->dSNIaFBDelay, smf->dSNFBDu, smf->dCCSNFBSpecEnergy,
                                                        smf->dSNIaFBSpecEnergy, smf->dConstGamma);
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

void smContext::Release(std::span<NN> list) {
    for (auto& nn : list) {
        if (nn.iPid != pkd->Self()) {
            pkd->mdl->Release(CID_PARTICLE, nn.pPart);
        }
    }
}

auto smContext::pqSearch(PQ* pq, TinyVector<double, 3> r, int iRoot) -> PQ* {
    auto mdl = pkd->mdl;
    int idSelf = pkd->Self();
    int iCell, id;
    double fDist2;
    TinyVector<double, 3> dr;

    /* Start at the root node of the tree */
    auto kdn = getCell(pkd->iTopTree[iRoot], id = idSelf);
    ST.clear();
    while (1) {
        while (kdn->is_cell()) {
            auto [iLower, idLower, iUpper, idUpper] = kdn->get_child_cells(id);
            kdn = getCell(iLower, idLower);
            auto min1 = kdn->bound().mindist(r);
            kdn = getCell(iUpper, idUpper);
            auto min2 = kdn->bound().mindist(r);
            if (min1 < min2) {
                if (min1 > pq->fDist2) goto NoIntersect;
                ST.push_back({.id = idUpper, .iCell = iUpper, .min = min2});
                id = idLower;
                iCell = iLower;
                kdn = getCell(iCell, id);
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
            auto [pStart, pEnd] = function->range_neighbour(kdn);
            for (auto pj = pStart; pj < pEnd; ++pj) {
                auto p = pkd->particles[pj];
                if (!p.marked() || !function->accept_neighbour(p)) continue;
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
            auto [pStart, pEnd] = function->range_neighbour(kdn);
            for (auto pj = pStart; pj < pEnd; ++pj) {
                auto p = pkd->particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
                if (!function->accept_neighbour(p)) continue;
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
            kdn = getCell(iCell, id);
            ST.pop_back();
        }
        else return pq;
    }
}

void smContext::create() {
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

void smSmooth(SMX smx, SMF* smf) {
    PKD pkd = smx->pkd;
    float fBall;

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (auto& p : pkd->particles) p.set_marked(true);
    smx->create();
    switch (smx->iSmoothType) {
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

    auto kdn = getCell(iCell = pkd->iTopTree[ROOT], id = idSelf);

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
            kdn = getCell(iCell, id);
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
            kdn = getCell(iCell, id);
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
    Release(nnList);
    return true;
}

uint64_t smReSmooth(SMX smx, int iSmoothType) {
    PKD pkd = smx->pkd;
    uint64_t nSmoothed = 0;

    switch (iSmoothType) {
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
