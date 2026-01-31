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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include "hop.h"
#include "groupstats.h"
#include "smooth/smooth.h"
#include <algorithm>

class HelperJoinLoops : public mdl::arrayCACHEhelper {
    bool bDone = true;
    void combine(void* dst, const void* src, const void* key) final {
        auto g1 = static_cast<GHtmpGroupTable*>(dst);
        const auto g2 = static_cast<const GHtmpGroupTable*>(src);
        if (g1->iPid > g2->iPid || (g1->iPid == g2->iPid && g1->iIndex > g2->iIndex)) {
            g1->iPid = g2->iPid;
            g1->iIndex = g2->iIndex;
            bDone = false;
        }
    }

  public:
    [[nodiscard]] auto done() const -> bool { return bDone; }
    explicit HelperJoinLoops(PKD pkd) : arrayCACHEhelper(pkd->tmpHopGroups, pkd->nGroups, true) {}
};

class HelperSetArc : public particleCACHEhelper {
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
        const auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];
        if (p2.is_marked()) p1.set_marked(true);
        assert(p1.group() == p2.group());
    }

  public:
    explicit HelperSetArc(PKD pkd) : particleCACHEhelper(pkd, true) {}
};

/*
** Follow the link for the given particle
*/
static int traverseLink(PKD pkd, struct smGroupArray* ga, int pi, int* iPid2, int* iIndex2, int* iMinPartPid, int* iMinPartIndex) {
    auto mdl = pkd->mdl;
    int gid2;
    GHtmpGroupTable* g2;
    int bDone = 0;
    int bRemote = *iPid2 != pkd->Self();

    auto p = (!bRemote) ? pkd->particles[*iIndex2] : pkd->particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, *iIndex2, *iPid2)];
    gid2 = p.group();
    g2 = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, gid2, *iPid2);

    /* Update the minimum group id in case we need it below */
    if (iMinPartPid && (*iPid2 < *iMinPartPid || (*iPid2 == *iMinPartPid && *iIndex2 < *iMinPartIndex))) {
        *iMinPartPid = *iPid2;
        *iMinPartIndex = *iIndex2;
    }

    /*
    ** This group points to the same node, so it is a terminal group.
    ** The particle ids may not match, but the group will point to
    ** a particle that points back to this group.
    */
    if (g2->iPid == *iPid2) {
        ga[pi].id.iPid = g2->iPid;
        ga[pi].id.iIndex = g2->iIndex;
        bDone = 1;
    }
    /*
    ** If we have linked back to our group on our node, then we are done.
    ** We choice the lowest gid we have seen so far. We also have to mark
    ** the remote particle that started this, because it starts an arc!
    */
    else if (*iPid2 == mdl->Self() && gid2 == pi) { /*LOOP*/
        auto p2 = pkd->particles[pkd->mdl->Acquire<PARTICLE>(CID_PARTICLE, ga[pi].id.iIndex, ga[pi].id.iPid)];
        p2.set_marked(true);
        pkd->mdl->Release(CID_PARTICLE, &p2);
        assert(ga[pi].id.iPid != pkd->Self());
        ga[pi].id.iPid = *iMinPartPid;
        ga[pi].id.iIndex = *iMinPartIndex;
        bDone = 1;
    }

    /* Follow the link - ignored if bDone is set. */
    *iPid2 = g2->iPid;
    *iIndex2 = g2->iIndex;

    pkd->mdl->Release(CID_GROUP, g2);

    return bDone;
}

/*
** Link particles based on density gradients
** After we finish, all "chains" will be complete globally across all domains.
** Chains still need to be joined if their terminating loops are close together.
** Prerequisites: Density and fBall from smooth
*/
static int smHopLink(SMX smx, std::shared_ptr<SmoothHelperGradientM3> helper) {
    PKD pkd = smx->pkd;
    auto mdl = pkd->mdl;
    int pi, j, gid, nLoop, nSpur;
    int nGroups, nLocal, nRemote;
    int iIndex1, iIndex2, iMinPartIndex, iPid1, iPid2, iMinPartPid, iParticle;
    struct smGroupArray* ga;
    /* Next particle in the chain or -1 if at the end in which case */
    /* the next particle will be remote and can be found in the group array */
    uint32_t* pl;

    pkd->ga = ga = pkd->Ephemeral<smGroupArray>();
    pl = (uint32_t*)((pkd->Ephemeral<char>() + 1ul * pkd->Local() * pkd->EphemeralElementBytes()) - (pkd->Local() + 1) * sizeof(uint32_t));
    assert((uint32_t*)ga < pl);

    ga[0].iGid = 0;
    ga[0].id.iPid = mdl->Self();
    ga[0].id.iIndex = -1; /* Sentinel */
    for (pi = 0; pi < pkd->Local(); ++pi) {
        auto p = pkd->particles[pi];
        p.set_group(0);     // Ungrouped
        p.set_marked(true); // Used by smooth to determine active particles
    }
    smx->create();

    pkd->nGroups = 0;
    pkd->nLocalGroups = 0;
    pkd->tmpHopGroups = NULL;

    /*
    ** Run through all of the particles, and find a particle to link to
    ** by calculating the gradient and looking in that direction. We make
    ** a link to this particle, then follow the link. Once we have reached
    ** the end of the chain, we consider the next particle (which may have
    ** been already done in which case we skip it).
    */
    nGroups = 1;
    nSpur = nLoop = nRemote = 0;
    for (pi = 0; pi < pkd->Local(); ++pi) {
        auto p = pkd->particles[iParticle = pi];
        if (p.group() > 0) continue; /* Already done (below) */
        assert((uint32_t*)&ga[nGroups + 1] < pl);
        ga[nGroups].iGid = nGroups;
        for (;;) {
            p.set_group(nGroups);
            /* Need: fDensity and fBall but pl[iParticle] is not yet used */
            /* fDensity MUST persist but fBall could be stashed into pl[] */
            smx->ReSmoothSingle(p, p.ball());
            /* Done: fBall */
            auto link = helper->link();
            ga[nGroups].id = link;
            /*
            ** Call mdlCacheCheck to make sure we are making progress!
            */
            mdl->CacheCheck();
            /*
            ** If we link to a remote node, then this group chain is done for now.
            ** Later we will merge remote groups.
            **/
            if (link.iPid != pkd->Self()) {
                ++nGroups;
                ++nRemote;
                break;
            }
            else {
                iParticle = pl[iParticle] = link.iIndex;
                p = pkd->particles[iParticle];
                gid = p.group();
                /* We link to ourselves: this forms a new "loop" */
                if (gid == nGroups) {
                    ++nGroups;
                    ++nLoop;
                    break;
                }
                /* Ok, some other group: merge this "spur" on to it */
                else if (gid > 0) {
                    for (j = pi; j != iParticle; j = pl[j]) {
                        p = pkd->particles[j];
                        p.set_group(gid);
                    }
                    ++nSpur;
                    break;
                }
            }
            /* None of the above? Just keep following the links. */
        }
    }
    smx->finish();
    /*
    ** All chains now terminate at a specific *particle*, either local or remote.
    ** Spurs were automatically joined, and point to their corresponding "loop".
    ** This "loop" particle is guaranteed to be part of the terminating loop.
    ** Remote particles are still pending, so we combine duplicates so there are
    ** fewer remote searches to perform. This results in substantially fewer chains.
    */
    mdl->FinishCache(CID_PARTICLE);
    nGroups = pkdGroupCombineDuplicateIds(pkd, nGroups, ga, 0);
    mdl->CacheInitialize(CID_PARTICLE, std::make_shared<particleCACHEhelper>(pkd));
    /*
    ** Now lets go looking for remote loops and see if we are part of them.
    ** The group table still contains PARTICLE links, so we use it to
    ** accelerate the search. We will find three cases:
    **  1. We find a terminating loop. Update the link to this loop.
    **  2. We are part of a terminating loop. Update the link to the lowest
    **     loop id and we are done. All members will do this.
    **  3. We encounter a loop that we are not a part of. Save any member
    **     of the loop. Later, we will fetch the correct id from it.
    */

    /* We will update bMarked (an arc) remotely, so we need a combiner cache */
    mdl->FinishCache(CID_PARTICLE);
    for (auto& p : pkd->particles) p.set_marked(false); // Not an arc (yet)
    mdl->CacheInitialize(CID_PARTICLE, std::make_shared<HelperSetArc>(pkd));
    pkd->tmpHopGroups = mdl->Malloc<GHtmpGroupTable>(nGroups * sizeof(GHtmpGroupTable));
    assert(pkd->tmpHopGroups != NULL);
    for (pi = 1; pi < nGroups; ++pi) {
        pkd->tmpHopGroups[pi].iPid = ga[pi].id.iPid;
        pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
    }
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(pkd->tmpHopGroups, nGroups));
    for (pi = 1; pi < nGroups; ++pi) {
        /* Completely local loop: any particle can start the "arc" */
        if (ga[pi].id.iPid == pkd->Self()) {
            pkd->particles[ga[pi].id.iIndex].set_marked(true);
        }
        /* Remote: may form a loop, or we may simply be a spur */
        else {
            iMinPartPid = iPid1 = iPid2 = ga[pi].id.iPid;
            iMinPartIndex = iIndex1 = iIndex2 = ga[pi].id.iIndex;
            for (;;) {
                /* We use Floyd's cycle-finding algorithm here */
                if (traverseLink(pkd, ga, pi, &iPid2, &iIndex2, &iMinPartPid, &iMinPartIndex)) break;
                if (traverseLink(pkd, ga, pi, &iPid2, &iIndex2, &iMinPartPid, &iMinPartIndex)) break;
                traverseLink(pkd, ga, pi, &iPid1, &iIndex1, NULL, NULL);
                /*
                ** If we see a loop here, then we are not a part of it (we are a spur).
                ** We point ourselves at one of the loop members so we can later
                ** retrieve the proper terminating particle. We can be pointing
                ** to our node of course, but to a different group.
                */
                if (iPid1 == iPid2 && iIndex1 == iIndex2) {
                    ga[pi].id.iPid = iPid1;
                    ga[pi].id.iIndex = iIndex1;
                    break;
                }
            }
        }
    }
    mdl->FinishCache(CID_GROUP);
    /*
    ** By now, we have marked the particle that starts each arc on the node,
    ** and exactly one particle if the group is a local loop.
    ** We now want to mark the rest of the particles in the arcs and loops.
    */
    mdl->FinishCache(CID_PARTICLE);
    for (pi = 0; pi < pkd->Local(); ++pi) {
        auto p = pkd->particles[pi];
        if (!p.marked()) continue;
        iIndex1 = pl[pi];
        while (iIndex1 != -1 && !(p = pkd->particles[iIndex1]).marked()) {
            p.set_marked(true);
            iIndex1 = pl[iIndex1];
        }
    }
    /* Done: pl[] is no longer needed */
    pl = NULL;

    nGroups = pkdGroupCombineDuplicateIds(pkd, nGroups, ga, 0);
    /*
    ** Now handle deferred groups (ones pointing to, but not part of, loops).
    */
    for (pi = 1; pi < nGroups; ++pi) {
        pkd->tmpHopGroups[pi].iPid = ga[pi].id.iPid;
        pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
    }
    mdl->CacheInitialize(CID_PARTICLE, std::make_shared<particleCACHEhelper>(pkd));
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(pkd->tmpHopGroups, nGroups));
    nSpur = 0;
    for (pi = 1; pi < nGroups; ++pi) {
        auto p1 = mdl->Acquire<PARTICLE>(CID_PARTICLE, ga[pi].id.iIndex, ga[pi].id.iPid);
        auto P1 = pkd->particles[p1];
        int gid1 = P1.group();
        auto g1 = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, gid1, ga[pi].id.iPid);
        auto p2 = mdl->Acquire<PARTICLE>(CID_PARTICLE, g1->iIndex, g1->iPid);
        auto P2 = pkd->particles[p2];
        int gid2 = P2.group();
        auto g2 = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, gid2, g1->iPid);
        assert(g2->iPid == g1->iPid);
        ga[pi].id.iPid = g2->iPid;
        ga[pi].id.iIndex = g2->iIndex;
        mdl->Release(CID_GROUP, g2);
        mdl->Release(CID_PARTICLE, p2);
        mdl->Release(CID_GROUP, g1);
        mdl->Release(CID_PARTICLE, p1);
    }
    mdl->FinishCache(CID_GROUP);

    /* If we end up point on-node, then we might have to combine the groups */
    for (pi = 1; pi < nGroups; ++pi) {
        if (ga[pi].id.iPid == pkd->Self()) {
            gid = pkd->particles[ga[pi].id.iIndex].group();
            if (gid != pi) {
                assert(ga[gid].id.iPid == pkd->Self());
                ga[pi].iGid = pi;
                ga[pi].id.iPid = ga[gid].id.iPid;
                ga[pi].id.iIndex = ga[gid].id.iIndex;
            }
        }
    }
    mdl->FinishCache(CID_PARTICLE);
    nGroups = pkdGroupCombineDuplicateIds(pkd, nGroups, ga, 0);
    mdl->CacheInitialize(CID_PARTICLE, std::make_shared<particleCACHEhelper>(pkd));

    /*
    ** Whew. The situation for groups is now one of the following.
    **  1. Groups that point locally form a local terminal loop or arc,
    **     so it must point to itself.
    **  2. Groups that point to a remote group must point to a remote
    **     group that points to itself.
    ** All particles that form a loop (either local or remote) have
    ** been marked (bMarked=1).
    */

    /* Turn particle IDs into group IDs */
    for (pi = 1; pi < nGroups; ++pi) {
        pkd->tmpHopGroups[pi].iPid = ga[pi].id.iPid;
        pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
    }
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(pkd->tmpHopGroups, nGroups));

    for (pi = 1; pi < nGroups; ++pi) {
        if (ga[pi].id.iPid == pkd->Self()) {
            gid = pkd->particles[ga[pi].id.iIndex].group();
            assert(gid == pi); /* I point to myself! */
        }
        else {
            auto p = pkd->particles[mdl->Acquire<PARTICLE>(CID_PARTICLE, ga[pi].id.iIndex, ga[pi].id.iPid)];
            gid = p.group();
            mdl->Release(CID_PARTICLE, &p);
            auto g = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, gid, ga[pi].id.iPid);
            /* The remote particle must point to itself */
            assert(g->iPid == ga[pi].id.iPid && g->iIndex == ga[pi].id.iIndex);
            mdl->Release(CID_GROUP, g);
        }
        ga[pi].id.iIndex = gid;
        assert(ga[pi].id.iPid != pkd->Self() || ga[pi].id.iIndex == pi);
    }
    mdl->FinishCache(CID_GROUP);

    pkd->nGroups = nGroups;
    nLocal = nRemote = 0;
    for (pi = 1; pi < pkd->nGroups; ++pi) {
        if (ga[pi].id.iPid == mdl->Self()) nLocal++;
        else nRemote++;
        pkd->tmpHopGroups[pi].iPid = ga[pi].id.iPid;
        pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
    }
    pkd->nLocalGroups = nLocal;

    /* tmpHopGroups is still allocated: freed in pkdHopFinishUp */
    return nLocal;
}

auto ServiceHopLink::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast</*const*/ input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto helper = std::make_shared<SmoothHelperGradientM3>(pkd);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
    *out = smHopLink(&smx, helper);
    return sizeof(output);
}

/*
** Join together chains if their terminating loops are withing a linking length.
**
** OPTIMIZATION: we should mark individual groups as done/not done so we can skip
**               subsequent iterations for those groups.
*/
static int smHopJoin(SMX smx, double dHopTau, int* nLocal) {
    PKD pkd = smx->pkd;
    auto mdl = pkd->mdl;
    auto pRoot = pkd->tree[ROOT];
    struct smGroupArray* ga = pkd->ga;
    int pi;

    auto join_helper = std::make_shared<HelperJoinLoops>(pkd);
    mdl->CacheInitialize(CID_GROUP, join_helper);
    /*
    ** Check all particles that are part of a loop or an arc (they are marked).
    ** We have contructed a tree with only marked particles, so just check those.
    */
    for (auto& p : *pRoot) {
        assert(p.marked());
        float fBall = (dHopTau < 0.0) ? -dHopTau * p.soft() : dHopTau;
        fBall = std::max(fBall, p.ball() * 0.5f);
        smx->ReSmoothSingle(p, fBall);
    }
    mdl->FinishCache(CID_GROUP);

    /* Follow the chains to the end */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(pkd->tmpHopGroups, pkd->nGroups));
    for (pi = 1; pi < pkd->nGroups; ++pi) {
        int iPid = pkd->Self();
        int iIndex = pi;
        int iNextPid = pkd->tmpHopGroups[iIndex].iPid;
        int iNextIndex = pkd->tmpHopGroups[iIndex].iIndex;
        while (iPid != iNextPid || iIndex != iNextIndex) {
            auto g = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, iNextIndex, iNextPid);
            iPid = iNextPid;
            iIndex = iNextIndex;
            iNextPid = g->iPid;
            iNextIndex = g->iIndex;
            mdl->Release(CID_GROUP, g);
        }
        pkd->tmpHopGroups[pi].iPid = iPid;
        pkd->tmpHopGroups[pi].iIndex = iIndex;
    }
    mdl->FinishCache(CID_GROUP);

    /* Merge duplicates for the next round */
    for (pi = 1; pi < pkd->nGroups; ++pi) {
        ga[pi].iGid = pi;
        ga[pi].id.iPid = pkd->tmpHopGroups[pi].iPid;
        ga[pi].id.iIndex = pkd->tmpHopGroups[pi].iIndex;
    }
    mdl->FinishCache(CID_PARTICLE);
    pkd->nGroups = pkdGroupCombineDuplicateIds(pkd, pkd->nGroups, ga, 1);
    mdl->CacheInitialize(CID_PARTICLE, std::make_shared<particleCACHEhelper>(pkd));
    int nRemote = 0;
    *nLocal = 0;
    for (pi = 1; pi < pkd->nGroups; ++pi) {
        if (ga[pi].id.iPid == mdl->Self()) (*nLocal)++;
        else nRemote++;
        pkd->tmpHopGroups[pi].iPid = ga[pi].id.iPid;
        pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
    }
    pkd->nLocalGroups = *nLocal;
    return join_helper->done();
}

auto ServiceHopJoin::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast</*const*/ input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto helper = std::make_shared<SmoothHelperHopChains>(pkd);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
    int nLocal;
    out->bDone = smHopJoin(&smx, in->dHopTau, &nLocal) && helper->done();
    out->nGroups = nLocal;
    return sizeof(output);
}

auto ServiceHopFinishUp::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vout.size() == sizeof(output));
    *out = pst->pkd()->HopFinishUp(in->nMinGroupSize, in->bPeriodic, in->fPeriod);
    pst->data_emplace<ServiceGroupStats::InfoPST>(ServiceGroupStats::InfoPST{.nLowerGroups = 0});
    return sizeof(output);
}

auto ServiceHopFinishUp::Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST pst) -> int {
    auto out1 = reinterpret_cast<output*>(vout1.data());
    assert(vout1.size() == sizeof(output));
    auto out2 = reinterpret_cast<const output*>(vout2.data());
    assert(vout2.size() == sizeof(output));
    pst->data_emplace<ServiceGroupStats::InfoPST>(ServiceGroupStats::InfoPST{.nLowerGroups = *out1});
    *out1 += *out2;
    return sizeof(output);
}

auto pkdContext::HopFinishUp(int nMinGroupSize, int bPeriodic, blitz::TinyVector<double, 3> dPeriod) -> int {
    struct smGroupArray* ga = Ephemeral<smGroupArray>();
    int i;

    for (i = 1; i < nGroups; ++i) {
        ga[i].iGid = i;
        ga[i].id.iPid = tmpHopGroups[i].iPid;
        ga[i].id.iIndex = tmpHopGroups[i].iIndex;
    }
    mdl->Free(tmpHopGroups);
    tmpHopGroups = NULL;

    nGroups = pkdPurgeSmallGroups(this, nGroups, ga, nMinGroupSize);
    hopSavedRoots = 0;
    return nLocalGroups;
}

class HelperHopGetRoots : public mdl::arrayCACHEhelper {
    class pkdContext* pkd;
    virtual void init(void* v) override final {
        auto g = static_cast<HopGroupTable*>(v);
        g->rmt.iPid = -1;
        g->rmt.iIndex = 0;
        g->iGlobalId = 0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<HopGroupTable*>(dst);
        const auto g2 = static_cast<const HopGroupTable*>(src);
        if (g2->iGlobalId) {
            int gid = g1 - pkd->hopGroups;
            int iRoot = pkd->hopRootIndex[gid]++;
            assert(gid < pkd->nLocalGroups); /* We should only get data for local groups */
            assert(pkd->hopRoots[iRoot].iPid == -1 && pkd->hopRoots[iRoot].iIndex == -1);
            pkd->hopRoots[iRoot].iPid = g2->rmt.iPid;
            pkd->hopRoots[iRoot].iIndex = g2->rmt.iIndex;
        }
    }

  public:
    explicit HelperHopGetRoots(pkdContext* pkd) : arrayCACHEhelper(pkd->hopRoots, pkd->nGroups, true), pkd(pkd) {}
};

auto ServiceHopTreeBuild::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->pkd();
    pkd->HopTreeBuild(in->nBucket, in->nGroup);
    return 0;
}

void pkdContext::HopTreeBuild(int nBucket, int nGroup) {
    int i, gid, iRoot;
    int nRootsTotal;

    /*
    ** The remote groups should be sorted by processor and remote group id,
    ** but we will rely on this behaviour later so we verify it now.
    */
    for (i = 2 + nLocalGroups; i < nGroups; ++i) {
        assert(hopGroups[i - 1].id.iPid < hopGroups[i].id.iPid ||
               (hopGroups[i - 1].id.iPid == hopGroups[i].id.iPid && hopGroups[i - 1].id.iIndex < hopGroups[i].id.iIndex));
    }

    /* Setup the buffer for tree roots */
    hopRootIndex = new int[nGroups + 1];
    /* Calculate the index into the tree roots table */
    nRootsTotal = 0;
    for (i = 1; i < nGroups; ++i) {
        hopRootIndex[i] = nRootsTotal;
        hopGroups[i].iAllRoots = nRootsTotal;
        nRootsTotal += 1 + hopGroups[i].nRemote;
    }

    hopRoots = new remoteID[nRootsTotal + 1];

    /* Invalidate all entries (for debug checks) */
    for (iRoot = 0; iRoot <= nRootsTotal; ++iRoot) hopRoots[iRoot].iPid = hopRoots[iRoot].iIndex = -1;

    /* We build a tree of each group of particles - hopGroups[].iTreeRoot is set */
    TreeBuildByGroup(nBucket, nGroup);

    /* Add ourselves for local groups */
    for (gid = 1; gid < nLocalGroups; ++gid) {
        iRoot = hopRootIndex[gid]++;
        assert(iRoot < nRootsTotal);
        assert(hopRoots[iRoot].iPid == -1 && hopRoots[iRoot].iIndex == -1);
        hopRoots[iRoot].iPid = mdl->Self();
        hopRoots[iRoot].iIndex = i;
    }

    /*
    ** Let the group master know about the remote tree roots. We can use the
    ** combiner cache because there is at most 1 root per group per node.
    */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<HelperHopGetRoots>(this));
    for (i = 1 + nLocalGroups; i < nGroups; ++i) {
        assert(hopGroups[i].id.iPid != mdl->Self());
        //  if (hopNumRoots[i] == 0) continue;
        auto g = mdl->Acquire<HopGroupTable>(CID_GROUP, hopGroups[i].id.iIndex, hopGroups[i].id.iPid);
        g->rmt.iPid = mdl->Self();
        g->rmt.iIndex = hopGroups[i].iTreeRoot;
        g->iGlobalId = i;
        hopGroups[i].iAllRoots = g->iAllRoots; /* We will fix this later! */
        mdl->Release(CID_GROUP, g);
    }
    mdl->FinishCache(CID_GROUP);

    /*
    ** Now we just need to fetch all of the roots for remote groups
    */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(hopRoots, nRootsTotal));

    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        int iIndex = hopGroups[gid].iAllRoots;
        for (i = 0; i < hopGroups[gid].nRemote + 1; ++i) {
            auto pRoot = mdl->Fetch<remoteID>(CID_GROUP, iIndex + i, hopGroups[gid].id.iPid);
            int iRoot = hopRootIndex[gid]++;
            assert(hopRoots[iRoot].iPid == -1 && hopRoots[iRoot].iIndex == -1);
            hopRoots[iRoot].iPid = pRoot->iPid;
            hopRoots[iRoot].iIndex = pRoot->iIndex;
        }
    }
    mdl->FinishCache(CID_GROUP);

    /* Fix the indexes */
    nRootsTotal = 0;
    for (gid = 1; gid < nGroups; ++gid) {
        hopGroups[gid].iAllRoots = nRootsTotal;
        nRootsTotal += 1 + hopGroups[gid].nRemote;
    }

    for (iRoot = 1; iRoot < nRootsTotal; ++iRoot) assert(hopRoots[iRoot].iPid >= 0 && hopRoots[iRoot].iIndex > 0);

    for (gid = 1; gid < nGroups; ++gid) {
        int iRootUs = hopRootIndex[gid];
        for (iRoot = hopRootIndex[gid]; iRoot < hopRootIndex[gid + 1]; ++iRoot) {
            assert(iRoot < nRootsTotal);
            assert(hopRoots[iRoot].iPid >= 0 && hopRoots[iRoot].iIndex > 0);
            /* Move our local root to the start */
            if (hopRoots[iRoot].iPid == mdl->Self() && iRoot > iRootUs) {
                remoteID t;
                assert(hopRoots[iRootUs].iPid != mdl->Self());
                t = hopRoots[iRootUs];
                hopRoots[iRootUs] = hopRoots[iRoot];
                hopRoots[iRoot] = t;
            }
        }
    }
    delete[] hopRootIndex;
    hopRootIndex = NULL;
    /* free: We have allocated pkd->hopRoots -- needed for gravity */
}

auto ServiceHopGravity::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->pkd();
    double dFlop = 0, dPartSum = 0, dCellSum = 0;
    pkd->HopGravity(in->dTime, in->dTheta, &dFlop, &dPartSum, &dCellSum);
    return 0;
}

auto ServiceHopUnbind::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->pkd();
    out->nEvaporated = pkd->HopUnbind(in->dTime, in->nMinGroupSize, in->bPeriodic, in->fPeriod);
    out->nGroups = pkd->nLocalGroups;
    return sizeof(output);
}

struct EE { /*EnergyElement*/
    double dPot;
    double dSqrtKin;
    double dTot;
    int i;
    bool operator<(const EE& other) const { return dTot < other.dTot; }
};

class HelperMaxHopEnergy : public mdl::arrayCACHEhelper {
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<HopGroupTable*>(dst);
        const auto g2 = static_cast<const HopGroupTable*>(src);
        if (g2->dEnergy > g1->dEnergy) g1->dEnergy = g2->dEnergy;
    }

  public:
    explicit HelperMaxHopEnergy(pkdContext* pkd) : arrayCACHEhelper(pkd->hopGroups, pkd->nGroups, true) {}
};

auto pkdContext::HopUnbind(double dTime, int nMinGroupSize, int bPeriodic, blitz::TinyVector<double, 3> dPeriod) -> int {
    auto ee = Ephemeral<EE>();
    int gid, iRoot, i, j;
    int nEvaporated;
    double a = csmTime2Exp(csm, dTime);
    double ia = 1.0 / a;
    double a2 = a * a;
    double dv, dv2;
    double dEnergy;

    mdl->CacheInitialize(CID_GROUP, std::make_shared<HelperMaxHopEnergy>(this));
    for (gid = 1; gid < nGroups; ++gid) {
        if (hopGroups[gid].bComplete) continue;
        iRoot = hopGroups[gid].iTreeRoot;
        auto pNode = tree[iRoot];
        // n = pNode->pUpper - pNode->pLower + 1;
        /* Calculate kinetic energy & total energy */
        for (auto& P : *pNode) {
            assert(P.group() == gid);
            auto v = P.velocity();
            dv2 = 0.0;
            for (j = 0; j < 3; ++j) {
                dv = v[j] - hopGroups[gid].vcom[j];
                dv2 += dv * dv;
            }
            ee[i].i = i;
            ee[i].dPot = P.potential() * ia;
            if (P.is_gas() && particles.present(PKD_FIELD::oSph)) { /* TODO: is this correct? */
                const auto& sph = P.sph();
                ee[i].dPot += sph.Uint / P.mass();
            }
            ee[i].dSqrtKin = a * sqrt(0.5 * dv2);
            ee[i].dTot = 0.5 * a2 * dv2 + ee[i].dPot;
        }
        std::sort(ee + pNode->lower(), ee + pNode->upper() + 1);
        i = pNode->count() > 100 ? pNode->upper() - 100 : pNode->lower();
        while (i <= pNode->upper() && ee[i].dTot <= 0.0) ++i;
        if (i > pNode->upper()) dEnergy = 0.0;
        else dEnergy = ee[i].dTot;
        hopGroups[gid].dEnergy = dEnergy;

        if (hopGroups[gid].id.iPid != Self()) {
            auto g = mdl->Acquire<HopGroupTable>(CID_GROUP, hopGroups[gid].id.iIndex, hopGroups[gid].id.iPid);
            g->dEnergy = dEnergy;
            mdl->Release(CID_GROUP, g);
        }
    }
    mdl->FinishCache(CID_GROUP);

    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(hopGroups, nGroups));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        if (hopGroups[gid].bComplete) continue;
        assert(hopGroups[gid].id.iPid != Self());
        auto g = mdl->Acquire<HopGroupTable>(CID_GROUP, hopGroups[gid].id.iIndex, hopGroups[gid].id.iPid);
        assert(hopGroups[gid].id.iIndex == g->id.iIndex && hopGroups[gid].id.iPid == g->id.iPid);
        hopGroups[gid].dEnergy = g->dEnergy;
        mdl->Release(CID_GROUP, g);
    }
    mdl->FinishCache(CID_GROUP);
    nEvaporated = 0;
    for (gid = 1; gid < nGroups; ++gid) {
        if (hopGroups[gid].bComplete) continue;
        iRoot = hopGroups[gid].iTreeRoot;
        dEnergy = hopGroups[gid].dEnergy;
        auto pNode = tree[iRoot];
        // n = pNode->count();
        for (i = pNode->upper(); i >= pNode->lower() && ee[i].dTot >= dEnergy; --i) {
            auto P = particles[ee[i].i];
            P.set_group(0);
            ++nEvaporated;
        }
        if (i == pNode->upper()) hopGroups[gid].bComplete = 1;
        /* Move evaporated particles to the end */
        else
            for (i = pNode->lower(); i <= pNode->upper();) {
                auto P = particles[i];
                if (P.group()) ++i;
                else {
                    auto Q = particles[pNode->upper()--];
                    swap(P, Q);
                }
            }
    }

    assert(0); /* the purge small groups below must work on ga... */
    nGroups = pkdPurgeSmallGroups(this, nGroups, ga, nMinGroupSize);

    return nEvaporated;
}
