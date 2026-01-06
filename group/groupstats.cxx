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
#include <algorithm>
#include <math.h>
#include "groupstats.h"
#include "group.h"

using blitz::dot;
using blitz::TinyVector;

auto ServiceGroupStats::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->lcl()->pkd;
    pkd->CalculateGroupStats(in->bPeriodic, in->dPeriod, in->rEnvironment, in->iGlobalStart);
    return 0;
}
auto ServiceGroupStats::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto mdl = pst->mdl();
    input inUpper = *in;
    inUpper.iGlobalStart += pst->nLowerGroups; /* for the upper subset the global id starting point is more by the number of lower groups */
    auto rID = ReqService(pst, inUpper);
    Traverse(pst->lower(), vin, vout);
    return mdl->GetReply(rID, vout);
}

#if 0
double illinois(double (*func)(double),double r,double s,double xacc,int *pnIter) {
    const int maxIter = 100;
    double t,fr,fs,ft,phis,phir,gamma;
    int i;

    fr = func(r);
    fs = func(s);
    assert(fr*fs < 0);
    t = (s*fr - r*fs)/(fr - fs);
    for (i=0; i < maxIter && fabs(t-s) > xacc; ++i) {

        ft = func(t);

        if (ft*fs < 0) {
            /*
            ** Unmodified step.
            */
            r = s;
            s = t;
            fr = fs;
            fs = ft;
        }
        else {
            /*
            ** Modified step to make sure we do not retain the
            ** endpoint r indefinitely.
            */
    #if 1
            phis = ft/fs;
            phir = ft/fr;
            gamma = 1 - (phis/(1-phir));  /* method 3 */
            if (gamma < 0) gamma = 0.5;
    #else
            gamma = 0.5;    /* illinois */
    #endif
            fr *= gamma;
            s = t;
            fs = ft;
        }
        t = (s*fr - r*fs)/(fr - fs);
        /*
            printf("%d %d r:%.14f fr:%.14f s:%.14f fs:%.14f t:%.14f\n",
                   i,r,fr,s,fs,t);
        */
    }
    *pnIter = i;
    return (t);
}
#endif

class MinPotHelper : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto g = static_cast<TinyGroupTable*>(dst);
        g->minPot = HUGE_VALF; /* make sure it is always set */
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<TinyGroupTable*>(dst);
        auto g2 = static_cast<const TinyGroupTable*>(src);
        if (g2->minPot < g1->minPot) {
            g1->minPot = g2->minPot;
            g1->rPot = g2->rPot;
        }
    }

  public:
    MinPotHelper(TinyGroupTable* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

class HelperTinyGroup : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto g = static_cast<TinyGroupTable*>(dst);
        g->fMass = 0.0;
        g->nBH = 0;
        g->nDM = 0;
        g->nGas = 0;
        g->nStar = 0;
        g->rcom = 0.0;
        g->vcom = 0.0;
        g->sigma = 0.0;
        g->rMax = 0.0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<TinyGroupTable*>(dst);
        auto g2 = static_cast<const TinyGroupTable*>(src);
        float x;

        g1->fMass += g2->fMass;
        x = g2->fMass / g1->fMass;
        g1->rcom = (1 - x) * g1->rcom + x * g2->rcom;
        g1->vcom = (1 - x) * g1->vcom + x * g2->vcom;
        g1->sigma = (1 - x) * g1->sigma + x * g2->sigma;
        if (g2->rMax > g1->rMax) g1->rMax = g2->rMax;
        g1->nBH += g2->nBH;
        g1->nDM += g2->nDM;
        g1->nGas += g2->nGas;
        g1->nStar += g2->nStar;
    }

  public:
    HelperTinyGroup(TinyGroupTable* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

class HelperTinyRmax : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto g = static_cast<TinyGroupTable*>(dst);
        g->rMax = 0.0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<TinyGroupTable*>(dst);
        auto g2 = static_cast<const TinyGroupTable*>(src);
        if (g2->rMax > g1->rMax) g1->rMax = g2->rMax;
    }

  public:
    HelperTinyRmax(TinyGroupTable* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

struct ShrinkStruct {
    double fMass;
    TinyVector<double, 3> rcom;
    int nEnclosed;
};

class HelperShrink : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto g = static_cast<ShrinkStruct*>(dst);

        g->fMass = 0.0;
        g->rcom = 0.0;
        g->nEnclosed = 0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<ShrinkStruct*>(dst);
        auto g2 = static_cast<const ShrinkStruct*>(src);
        float x;

        g1->fMass += g2->fMass;
        x = g2->fMass / g1->fMass;
        g1->rcom = (1 - x) * g1->rcom + x * g2->rcom;
        g1->nEnclosed += g2->nEnclosed;
    }

  public:
    HelperShrink(ShrinkStruct* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

class HelperAngular : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto g = static_cast<TinyGroupTable*>(dst);
        g->angular = 0.0;
        g->inertia = 0.0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<TinyGroupTable*>(dst);
        auto g2 = static_cast<const TinyGroupTable*>(src);
        g1->angular += g2->angular;
        g1->inertia += g2->inertia;
    }

  public:
    HelperAngular(TinyGroupTable* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

struct MassRadius {
    float fMass;
    float dr;
    auto operator==(const MassRadius& other) const { return dr == other.dr; }
    auto operator<=>(const MassRadius& other) const { return dr <=> other.dr; }
};

class HelperRoot : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto g = static_cast<MassRadius*>(dst);
        g->fMass = 0.0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto g1 = static_cast<MassRadius*>(dst);
        auto g2 = static_cast<const MassRadius*>(src);
        g1->fMass += g2->fMass;
    }

  public:
    HelperRoot(MassRadius* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

auto pkdContext::gatherMass(remoteID* S, double fBall2, double ri2, double ro2, TinyVector<double, 3> r) -> double {
    double min2, max2;
    int iCell, id;
    int sp = 0;
    double dx, dy, dz, fDist2;
    int pj, pEnd;
    double fMass = 0;

    auto kdn = getCell(iCell = iTopTree[ROOT], id = Self());
    while (1) {
        Bound bnd = kdn->bound();
        min2 = bnd.mindist(r);
        if (min2 > ri2) goto NoIntersect;
        max2 = bnd.maxdist(r);
        if (max2 <= ro2) {
            fMass += kdn->mass();
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            int idUpper, iUpper;
            std::tie(iCell, id, iUpper, idUpper) = kdn->get_child_cells(id);
            kdn = getCell(iCell, id);
            S[sp].iPid = idUpper;
            S[sp].iIndex = iUpper;
            ++sp;
            continue;
        }
        else {
            if (id == Self()) {
                pEnd = kdn->upper();
                for (pj = kdn->lower(); pj <= pEnd; ++pj) {
                    auto p = particles[pj];
                    auto p_r = p.position();
                    dx = r[0] - p_r[0];
                    dy = r[1] - p_r[1];
                    dz = r[2] - p_r[2];
                    fDist2 = dx * dx + dy * dy + dz * dz;
                    if (fDist2 <= fBall2) {
                        fMass += p.mass();
                    }
                }
            }
            else {
                pEnd = kdn->upper();
                for (pj = kdn->lower(); pj <= pEnd; ++pj) {
                    auto p = particles[mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
                    auto p_r = p.position();
                    dx = r[0] - p_r[0];
                    dy = r[1] - p_r[1];
                    dz = r[2] - p_r[2];
                    fDist2 = dx * dx + dy * dy + dz * dz;
                    if (fDist2 <= fBall2) {
                        fMass += p.mass();
                    }
                }
            }
        }
    NoIntersect:
        if (sp) {
            --sp;
            id = S[sp].iPid;
            iCell = S[sp].iIndex;
            kdn = getCell(iCell, id);
        }
        else break;
    }
    return (fMass);
}

#if 0
static double gatherLocalMass(PKD pkd,remoteID *S,double fBall2,double ri2,double ro2,TinyVector<double,3> r) {
    double min2,max2;
    int iCell;
    int sp = 0;
    double dx,dy,dz,fDist2;
    double fMass = 0;

    auto kdn = pkd->tree[iCell=ROOT];
    while (1) {
        Bound bnd = kdn.bound();
        min2 = bnd.mindist(r);
        if (min2 > ri2) goto NoIntersect;
        max2 = bnd.maxdist(r);
        if (max2 <= ro2) {
            fMass += kdn.mass();
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->iLower) {
            kdn = pkd->tree[iCell = kdn->iLower];
            S[sp++].iIndex = iCell+1;
            continue;
        }
        else {
            for (auto &p : kdn) {
                auto p_r = p.position();
                dx = r[0] - p_r[0];
                dy = r[1] - p_r[1];
                dz = r[2] - p_r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                    fMass += p.mass();
                }
            }
        }
NoIntersect:
        if (sp) kdn = pkd->tree[iCell = S[--sp].iIndex];
        else break;
    }
    return (fMass);
}
#endif

auto pkdContext::GatherMass(remoteID* S, double fBall, TinyVector<double, 3> r, int bPeriodic, TinyVector<double, 3> dPeriod) -> double {
    TinyVector<double, 3> rp;
    int ix, iy, iz, nReps;
    double fMass = 0;
    double fBall2 = fBall * fBall;
    double ri, ri2, ro2;

    assert(tree.present(KDN_FIELD::oNodeMom));
    ri = 0.9 * fBall; /* we use an approximate radius but make sure it is unbiased by volume */
    ri2 = ri * ri;
    ro2 = pow(2 * fBall2 * fBall - ri2 * ri, 2.0 / 3.0);
    if (bPeriodic) nReps = 1;
    else nReps = 0;
    for (ix = -nReps; ix <= nReps; ++ix) {
        rp[0] = r[0] + ix * dPeriod[0];
        for (iy = -nReps; iy <= nReps; ++iy) {
            rp[1] = r[1] + iy * dPeriod[1];
            for (iz = -nReps; iz <= nReps; ++iz) {
                rp[2] = r[2] + iz * dPeriod[2];
                fMass += gatherMass(S, fBall2, ri2, ro2, rp);
            }
        }
    }
    return (fMass);
}

struct RootFindingTable {
    float a;
    float b;
    int bConverged;
    int iter; /* this is only used for a diagnostic */
    int gid;
};

void pkdContext::CalculateGroupStats(int bPeriodic, TinyVector<double, 3> dPeriod, std::array<double, 2> rEnvironment, uint64_t iGlobalStart) {
    const int bDoShrinkingSphere = 0;
    double fMass;
    TinyVector<double, 3> r;
    double v2;
    double diVol0, diVol1;
    float r2, rMax;
    int i, j, gid, n;
    int nLocalGroups;
    double* dAccumulate;
    MassRadius *mr, *mrFree, *rootFunction;
    remoteID* S;
    uint32_t *iGrpOffset, *iGrpEnd;
    int nRootFind, bIncomplete, nMaxIter, iter;
    int *mrIndex, iRoot, *bRemoteDone;
    RootFindingTable* rootFindingTable;
    ShrinkStruct* shrink;
    int bShrink;
    double f2;
    int bLocal; /* a flip-flop which toggles when the first remote group is encountered, used as a consistency check only */

    assert(nGroups * (sizeof(*ga) + sizeof(*tinyGroupTable) + sizeof(ShrinkStruct)) < 1ul * EphemeralElementBytes() * FreeStore());
    tinyGroupTable = (TinyGroupTable*)(&ga[nGroups]);
    dAccumulate = (double*)(&tinyGroupTable[nGroups]);
    /*
    ** Initialize the table.
    */
    nLocalGroups = 0;
    /*
    ** Global id of ungrouped particles is always 0!
    */
    tinyGroupTable[0].iGlobalGid = 0;
    bLocal = 1;
    for (gid = 0; gid < nGroups; ++gid) {
        if (ga[gid].id.iPid == Self() && gid) {
            assert(bLocal == 1); /* make sure all the local groups come first in the ga array!!! */
            tinyGroupTable[gid].iGlobalGid = iGlobalStart + nLocalGroups++;
        }
        else if (gid) {
            bLocal = 0; /* at this point should follow all remote groups */
        }
        /*
        ** Here we assume that the ga table was setup prior to calling
        ** gravity so that we have filled in the minimum potential particle
        ** for each group.
        */
        tinyGroupTable[gid].minPot = ga[gid].minPot;
        if (gid) {
            assert(ga[gid].minPot < FLOAT_MAXVAL);
            assert(ga[gid].iMinPart < Local());
            auto p = particles[ga[gid].iMinPart];
            tinyGroupTable[gid].rPot = p.position();
            tinyGroupTable[gid].rcen = 0.0; /* relative to rPot */
        }
        tinyGroupTable[gid].rMax = 0;
        tinyGroupTable[gid].nBH = 0;
        tinyGroupTable[gid].nDM = 0;
        tinyGroupTable[gid].nGas = 0;
        tinyGroupTable[gid].nStar = 0;
        dAccumulate[4 * gid] = 0;
        dAccumulate[4 * gid + 1] = 0;
        dAccumulate[4 * gid + 2] = 0;
        dAccumulate[4 * gid + 3] = 0;
    }
    /*
    ** Make sure we got a consistent number of local groups from after FofFinishUp().
    */
    assert(this->nLocalGroups == nLocalGroups);
    /*
    ** First determine the minimum potential particle for each group.
    ** This will be the reference position for the group as well.
    ** In this first pass we also sum up the v2 and vcom locally.
    */
    for (i = 0; i < Local(); ++i) {
        auto p = particles[i];
        gid = p.group();
        fMass = p.mass();
        if (!gid) continue;
        auto v = p.velocity();
        v2 = 0;
        for (j = 0; j < 3; ++j) {
            dAccumulate[4 * gid + j] += fMass * v[j];
            v2 += v[j] * v[j];
        }
        dAccumulate[4 * gid + 3] += fMass * v2;
    }
    for (gid = 1; gid < nGroups; ++gid) {
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].vcom[j] = dAccumulate[4 * gid + j];
            dAccumulate[4 * gid + j] = 0;
        }
        tinyGroupTable[gid].sigma = dAccumulate[4 * gid + 3];
        dAccumulate[4 * gid + 3] = 0;
    }
    /*
    ** Now look at remote group particle to see if we have a lower potential.
    ** Note that we only expose the local groups to the cache! This allows us
    ** to make any desired update to the remote group entries of the table.
    */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<MinPotHelper>(tinyGroupTable, nLocalGroups + 1));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        assert(ga[gid].id.iPid != Self());
        auto g = mdl->VirtualFetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
        /*
        ** We update the remote center, the combiner makes sure the minimum minPot is set.
        ** We don't have to check the virtual fetch value in this case since we fetch each
        ** remote group only once (there cannot be more than one reference) otherwise we
        ** would have to check the fetched value.
        */
        g->minPot = tinyGroupTable[gid].minPot;
        for (j = 0; j < 3; ++j) {
            g->rPot[j] = tinyGroupTable[gid].rPot[j];
        }
    }
    mdl->FinishCache(CID_GROUP);
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(tinyGroupTable, nLocalGroups + 1));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        auto g = mdl->Fetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
        /*
        ** We update the local version of the remote center, here we can just copy it.
        */
        tinyGroupTable[gid].minPot = g->minPot;
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].rPot[j] = g->rPot[j];
        }
        /*
        ** Read the remote group's global id.
        */
        tinyGroupTable[gid].iGlobalGid = g->iGlobalGid;
    }
    mdl->FinishCache(CID_GROUP);
    /*
    ** Now we can set the GlobalGid of every particle (if this field is present)!
    */
    if (particles.present(PKD_FIELD::oGlobalGid)) {
        for (i = 0; i < Local(); ++i) {
            auto p = particles[i];
            gid = p.group();
            p.global_gid() = tinyGroupTable[gid].iGlobalGid;
        }
    }
    /*
    ** Now find rcom and rMax and total mass locally.
    */
    TinyVector<double, 3> dHalf;
    if (bPeriodic) dHalf = 0.5 * dPeriod;
    else dHalf = std::numeric_limits<double>::max();
    for (i = 0; i < Local(); ++i) {
        auto p = particles[i];
        gid = p.group();
        fMass = p.mass();
        int pSpecies = p.species();
        switch (pSpecies) {
        case FIO_SPECIES_BH:
            tinyGroupTable[gid].nBH++;
            break;
        case FIO_SPECIES_DARK:
            tinyGroupTable[gid].nDM++;
            break;
        case FIO_SPECIES_SPH:
            tinyGroupTable[gid].nGas++;
            break;
        case FIO_SPECIES_STAR:
            tinyGroupTable[gid].nStar++;
            break;
        default: // Maybe a deleted particle, skip it
            continue;
        }
        dAccumulate[4 * gid + 3] += fMass;
        if (gid > 0) {
            r = p.position() - tinyGroupTable[gid].rPot;
            r2 = 0.0;
            for (j = 0; j < 3; ++j) {
                if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
                dAccumulate[4 * gid + j] += fMass * r[j];
                r2 += r[j] * r[j];
            }
            rMax = sqrtf(r2);
            if (rMax > tinyGroupTable[gid].rMax) tinyGroupTable[gid].rMax = rMax;
        }
    }
    for (gid = 1; gid < nGroups; ++gid) {
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].rcom[j] = dAccumulate[4 * gid + j] / dAccumulate[4 * gid + 3];
            tinyGroupTable[gid].vcom[j] /= dAccumulate[4 * gid + 3];
            dAccumulate[4 * gid + j] = 0;
        }
        tinyGroupTable[gid].sigma /= dAccumulate[4 * gid + 3]; /* "sigma" is actually still v^2 */
        tinyGroupTable[gid].fMass = dAccumulate[4 * gid + 3];
        dAccumulate[4 * gid + 3] = 0; /* if we still need masses later we should not clear this field */
    }
    /*
    ** Now accumulate totals globally
    */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<HelperTinyGroup>(tinyGroupTable, nLocalGroups + 1));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        auto g = mdl->VirtualFetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
        g->fMass += tinyGroupTable[gid].fMass;
        g->nBH += tinyGroupTable[gid].nBH;
        g->nDM += tinyGroupTable[gid].nDM;
        g->nGas += tinyGroupTable[gid].nGas;
        g->nStar += tinyGroupTable[gid].nStar;
        for (j = 0; j < 3; j++) {
            g->rcom[j] += tinyGroupTable[gid].rcom[j];
            g->vcom[j] += tinyGroupTable[gid].vcom[j];
        }
        g->sigma += tinyGroupTable[gid].sigma;
        if (tinyGroupTable[gid].rMax > g->rMax) g->rMax = tinyGroupTable[gid].rMax;
    }
    mdl->FinishCache(CID_GROUP);
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(tinyGroupTable, nLocalGroups + 1));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        auto g = mdl->Fetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
        tinyGroupTable[gid].fMass = g->fMass;
        tinyGroupTable[gid].nBH = g->nBH;
        tinyGroupTable[gid].nDM = g->nDM;
        tinyGroupTable[gid].nGas = g->nGas;
        tinyGroupTable[gid].nStar = g->nStar;
        for (j = 0; j < 3; j++) {
            tinyGroupTable[gid].rcom[j] = g->rcom[j];
            tinyGroupTable[gid].vcom[j] = g->vcom[j];
            tinyGroupTable[gid].sigma = g->sigma;
        }
        tinyGroupTable[gid].rMax = g->rMax;
    }
    mdl->FinishCache(CID_GROUP);

    for (gid = 1; gid < nGroups; ++gid) {
        v2 = 0;
        for (j = 0; j < 3; j++) {
            /*
            ** If we want absolute rcom instead of relative to the minimum potential then we
            ** need to uncomment the line below.
            ** tinyGroupTable[gid].rcom[j] += tinyGroupTable[gid].rPot[j];
            */
            v2 += tinyGroupTable[gid].vcom[j] * tinyGroupTable[gid].vcom[j];
        }
        /*
        ** Now convert from mass weighted average v2 to sigma.
        */
        tinyGroupTable[gid].sigma -= v2;
        tinyGroupTable[gid].sigma = sqrtf(tinyGroupTable[gid].sigma);
    }
    if (bDoShrinkingSphere) {
        /*
        ** Set the center of each group to be the center of mass, at least initially...
        */
        for (gid = 1; gid < nGroups; ++gid) {
            for (j = 0; j < 3; ++j) {
                tinyGroupTable[gid].rcen[j] = tinyGroupTable[gid].rcom[j];
            }
            tinyGroupTable[gid].rMax = 0; /* reset this here since we will recalculate it about the new center */
        }
        /*
        ** Now find rMax about the center (which is currently the center of mass).
        */
        for (i = 0; i < Local(); ++i) {
            auto p = particles[i];
            gid = p.group();
            if (gid > 0) {
                r = p.position() - tinyGroupTable[gid].rPot;
                for (j = 0; j < 3; ++j) {
                    if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                    else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
                    /*
                    ** Make the position relative to the current best center.
                    ** This starts as being the center of mass of the entire group.
                    */
                }
                r -= tinyGroupTable[gid].rcen;
                r2 = dot(r, r);
                rMax = sqrtf(r2);
                if (rMax > tinyGroupTable[gid].rMax) tinyGroupTable[gid].rMax = rMax;
            }
        }
        /*
        ** Now find the largest rmax for the local groups globally.
        */
        mdl->CacheInitialize(CID_GROUP, std::make_shared<HelperTinyRmax>(tinyGroupTable, nLocalGroups + 1));
        for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
            auto g = mdl->VirtualFetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
            g->rMax = tinyGroupTable[gid].rMax;
        }
        mdl->FinishCache(CID_GROUP);
        mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(tinyGroupTable, nLocalGroups + 1));
        for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
            auto g = mdl->Fetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
            tinyGroupTable[gid].rMax = g->rMax;
        }
        mdl->FinishCache(CID_GROUP);

        shrink = (ShrinkStruct*)dAccumulate;                                        /* overlay the shrink structure on the accumulators */
        for (gid = 0; gid < nGroups; ++gid) shrink[gid].nEnclosed = ga[gid].nTotal; /* all particles in the group */

        /*
        ** We also may want the "shrinking sphere" center for many applications.
        ** At this point we have a sphere of maximum radius enclosing all the particles
        ** as well as an initial guess of the center. We now shrink the sphere by
        ** 2.5% (a la Power, et. al. 2003) and find a new center of mass. We continue this
        ** process until we reach a minumum number of particles within the sphere (~1000).
        */
        f2 = 1.0;
        while (1) {
            /*
            ** Start of the main shrinking iteration.
            */
            bShrink = 0; /* by default we will not shrink anything */
            /*
            ** Clear the sphere totals
            */
            for (gid = 0; gid < nGroups; ++gid) {
                for (j = 0; j < 3; ++j) shrink[gid].rcom[j] = 0.0;
                shrink[gid].fMass = 0.0;
                if (shrink[gid].nEnclosed > 1000) {
                    shrink[gid].nEnclosed = 0;
                    bShrink = 1;
                }
                else {
                    shrink[gid].nEnclosed = -1; /* marks not to shrink this group */
                }
            }
            shrink[0].nEnclosed = -1; /* make sure ungrouped particles (gid==0) are not involved */
            if (!bShrink) break;      /* exits the main shrinking iteration...we are done */
            f2 *= 0.975 * 0.975;      /* make the shrink factor 2.5 % smaller */
            for (i = 0; i < Local(); ++i) {
                auto p = particles[i];
                gid = p.group();
                if (shrink[gid].nEnclosed < 0) continue; /* skip this particle */
                r = p.position() - tinyGroupTable[gid].rPot;
                for (j = 0; j < 3; ++j) {
                    if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                    else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
                }
                r -= tinyGroupTable[gid].rcen;
                r2 = dot(r, r);
                /*
                ** If the particle lies within the shrinking sphere, then add use it to calculate the center of mass.
                */
                if (r2 < f2 * tinyGroupTable[gid].rMax * tinyGroupTable[gid].rMax) {
                    fMass = p.mass();
                    for (j = 0; j < 3; ++j) shrink[gid].rcom[j] += fMass * (r[j] + tinyGroupTable[gid].rcen[j]);
                    shrink[gid].fMass += fMass;
                    shrink[gid].nEnclosed += 1; /* count number of particles inside rmax */
                }
            }
            for (gid = 1; gid < nGroups; ++gid) {
                for (j = 0; j < 3; ++j) shrink[gid].rcom[j] /= shrink[gid].fMass;
            }
            /*
            ** Now accumulate totals globally, we just sum the 5 fields for each group.
            */
            mdl->CacheInitialize(CID_GROUP, std::make_shared<HelperShrink>(shrink, nLocalGroups + 1));
            for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
                if (shrink[gid].nEnclosed < 0) continue;
                auto g = mdl->VirtualFetch<ShrinkStruct>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
                for (j = 0; j < 3; j++) {
                    g->rcom[j] = shrink[gid].rcom[j];
                }
                g->fMass = shrink[gid].fMass;
                g->nEnclosed = shrink[gid].nEnclosed;
            }
            mdl->FinishCache(CID_GROUP);
            mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(shrink, nLocalGroups + 1));
            for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
                if (shrink[gid].nEnclosed < 0) continue;
                auto g = mdl->Fetch<ShrinkStruct>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
                for (j = 0; j < 3; j++) {
                    shrink[gid].rcom[j] = g->rcom[j];
                }
                shrink[gid].fMass = g->fMass;
                shrink[gid].nEnclosed = g->nEnclosed;
            }
            mdl->FinishCache(CID_GROUP);
            for (gid = 1; gid < nGroups; ++gid) {
                if (shrink[gid].nEnclosed > 0) {
                    assert(shrink[gid].fMass > 0.0);
                    for (j = 0; j < 3; ++j) {
                        tinyGroupTable[gid].rcen[j] = shrink[gid].rcom[j];
                    }
                }
            }
        }
        shrink = NULL; /* make sure we can't use it anymore by accident */
    } /* end of if (doShrinkingSphere) */
    /*
    ** Calculate angular momentum about the center of the group.
    */
    for (i = 0; i < Local(); ++i) {
        auto p = particles[i];
        gid = p.group();
        if (gid > 0) {
            fMass = p.mass();
            auto v = p.velocity();
            r = p.position() - tinyGroupTable[gid].rPot;
            for (j = 0; j < 3; ++j) {
                if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
            }
            r -= tinyGroupTable[gid].rcen;
            dAccumulate[4 * gid + 0] += fMass * (r[1] * v[2] - r[2] * v[1]);
            dAccumulate[4 * gid + 1] += fMass * (r[2] * v[0] - r[0] * v[2]);
            dAccumulate[4 * gid + 2] += fMass * (r[0] * v[1] - r[1] * v[0]);
        }
    }
    for (gid = 1; gid < nGroups; ++gid) {
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].angular[j] = dAccumulate[4 * gid + j];
            dAccumulate[4 * gid + j] = 0;
        }
    }
    /*
    ** Calculate the moment of inertia (we do this in 2 stages to save accumulator memory.
    */
    for (i = 0; i < Local(); ++i) {
        auto p = particles[i];
        gid = p.group();
        if (gid > 0) {
            fMass = p.mass();
            r = p.position() - tinyGroupTable[gid].rPot;
            for (j = 0; j < 3; ++j) {
                if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
            }
            r -= tinyGroupTable[gid].rcen;
            dAccumulate[4 * gid + 0] += fMass * r[0] * r[0];
            dAccumulate[4 * gid + 1] += fMass * r[0] * r[1];
            dAccumulate[4 * gid + 2] += fMass * r[0] * r[2];
        }
    }
    for (gid = 1; gid < nGroups; ++gid) {
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].inertia[j] = dAccumulate[4 * gid + j];
            dAccumulate[4 * gid + j] = 0;
        }
    }
    for (i = 0; i < Local(); ++i) {
        auto p = particles[i];
        gid = p.group();
        if (gid > 0) {
            fMass = p.mass();
            for (j = 1; j < 3; ++j) { /* careful, this loop skips r[0]! */
                r[j] = p.position(j) - tinyGroupTable[gid].rPot[j];
                if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
                r[j] -= tinyGroupTable[gid].rcen[j];
            }
            dAccumulate[4 * gid + 0] += fMass * r[1] * r[1];
            dAccumulate[4 * gid + 1] += fMass * r[1] * r[2];
            dAccumulate[4 * gid + 2] += fMass * r[2] * r[2];
        }
    }
    for (gid = 1; gid < nGroups; ++gid) {
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].inertia[j + 3] = dAccumulate[4 * gid + j];
            dAccumulate[4 * gid + j] = 0;
        }
    }
    /*
    ** Now accumulate angular momentum globally
    */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<HelperAngular>(tinyGroupTable, nLocalGroups + 1));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        auto g = mdl->VirtualFetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
        for (j = 0; j < 3; j++) {
            g->angular[j] += tinyGroupTable[gid].angular[j];
        }
        for (j = 0; j < 6; ++j) {
            g->inertia[j] += tinyGroupTable[gid].inertia[j];
        }
    }
    mdl->FinishCache(CID_GROUP);
    /*
    ** Scoop environment mass (note the full tree must be present).
    ** At the very least we need to have masses in the cells, either as part of the FMOMR structure
    ** or a seperate mass field!
    */
    S = (remoteID*)(&tinyGroupTable[nGroups]);
    mdl->CacheInitialize(CID_PARTICLE, std::make_shared<particleCACHEhelper>(this));
    diVol0 = 3.0 / 4.0 * M_1_PI * pow(rEnvironment[0], -3);
    diVol1 = 3.0 / 4.0 * M_1_PI * pow(rEnvironment[1], -3);
    for (gid = 1; gid <= nLocalGroups; ++gid) {
        for (j = 0; j < 3; ++j) r[j] = tinyGroupTable[gid].rPot[j];
        if (rEnvironment[0] > 0.0) {
            tinyGroupTable[gid].fEnvironDensity0 = GatherMass(S, rEnvironment[0], r, bPeriodic, dPeriod) * diVol0;
        }
        else tinyGroupTable[gid].fEnvironDensity0 = 0.0;
        if (rEnvironment[1] > 0.0) {
            tinyGroupTable[gid].fEnvironDensity1 = GatherMass(S, rEnvironment[1], r, bPeriodic, dPeriod) * diVol1;
        }
        else tinyGroupTable[gid].fEnvironDensity1 = 0.0;
    }
    mdl->FinishCache(CID_PARTICLE);
    /*
    ** Fetch remote group properties (we probably only need rMax to be fetched here).
    */
    mdl->CacheInitialize(CID_GROUP, std::make_shared<mdl::arrayCACHEhelper>(tinyGroupTable, nLocalGroups + 1));
    for (gid = 1 + nLocalGroups; gid < nGroups; ++gid) {
        auto g = mdl->Fetch<TinyGroupTable>(CID_GROUP, ga[gid].id.iIndex, ga[gid].id.iPid);
        tinyGroupTable[gid].fMass = g->fMass;
        for (j = 0; j < 3; ++j) {
            tinyGroupTable[gid].rcom[j] = g->rcom[j];
            tinyGroupTable[gid].vcom[j] = g->vcom[j];
            tinyGroupTable[gid].angular[j] = g->angular[j];
        }
        for (j = 0; j < 6; ++j) {
            tinyGroupTable[gid].inertia[j] = g->inertia[j];
        }
        tinyGroupTable[gid].sigma = g->sigma;
        tinyGroupTable[gid].rMax = g->rMax;
        tinyGroupTable[gid].fEnvironDensity0 = g->fEnvironDensity0;
        tinyGroupTable[gid].fEnvironDensity1 = g->fEnvironDensity1;
    }
    mdl->FinishCache(CID_GROUP);
    /*
    ** Now find the half mass radii of the groups.
    ** Start by reordering particles into group order.
    ** Note that we must not access tree cells until the next tree build!
    */
    mdl->FinishCache(CID_CELL);
    iGrpOffset = (uint32_t*)(&tinyGroupTable[nGroups]);
    mrIndex = (int*)(&iGrpOffset[2 * nGroups]); /* this stores the ending index of each group in the mr below */
    mr = (MassRadius*)(&mrIndex[nGroups]);      /* we must be careful how much of this we use! */
    GroupOrder(iGrpOffset);
    iGrpEnd = &iGrpOffset[nGroups];
    iGrpEnd[0] = 0;
    /*
    ** First do all purely local groups. We can do these one at a time without
    ** worrying about remote particles. We should use the bounds determined by
    ** fof previously.
    */
    mrFree = mr;
    nRootFind = 0;
    mrIndex[0] = 0;
    for (gid = 1; gid <= nLocalGroups; ++gid) {
        rMax = std::min(tinyGroupTable[gid].rMax,
                        static_cast<float>(blitz::min(bndInterior.apothem() - blitz::abs(bndInterior.center() - tinyGroupTable[gid].rPot))));
        n = 0;
        for (i = iGrpEnd[gid - 1]; i < iGrpEnd[gid]; ++i, ++n) {
            auto p = particles[i];
            assert(gid == p.group()); /* make sure it is really in group order */
            r = p.position() - tinyGroupTable[gid].rPot;
            for (j = 0; j < 3; ++j) {
                if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
            }
            r -= tinyGroupTable[gid].rcen;
            r2 = dot(r, r);
            mrFree[n].fMass = p.mass();
            mrFree[n].dr = sqrtf(r2);
        }
        std::sort(mrFree, mrFree + n);
        fMass = mrFree[0].fMass;
        for (i = 1; i < n; ++i) {
            fMass += mrFree[i].fMass;
            mrFree[i].fMass = fMass;
        }
        /*
        ** If the radius in this sum ever exceeds rMax, then we have to resolve this
        ** group's rHalf by using particles from the other remote processors.
        */
        bIncomplete = 0;
        for (i = 0; i < n; ++i) {
            if (mrFree[i].dr > rMax) {
                bIncomplete = 1;
                break;
            }
            else if (mrFree[i].fMass > 0.5 * tinyGroupTable[gid].fMass) break;
        }
        tinyGroupTable[gid].rHalf = (i > 0) ? mrFree[--i].dr : 0;
        if (bIncomplete) {
            /*
            ** Compact the mrFree array.
            */
            for (j = 0; i < n; ++i, ++j) mrFree[j] = mrFree[i];
            /*
            ** We need to do some root finding on this group in order to determine the
            ** correct rHalf. The stored rHalf is only a lower bound at present in this
            ** case. Only the rootTable is exposed to the remote processors via the cache.
            */
            mrIndex[gid] = mrFree - mr + j;
            /*
            ** Set the new mrFree to unallocated space.
            */
            mrFree = mr + mrIndex[gid];
            ++nRootFind;
        }
        else {
            mrIndex[gid] = mrFree - mr;
        }
    }
    for (gid = nLocalGroups + 1; gid < nGroups; ++gid) {
        n = 0;
        for (i = iGrpEnd[gid - 1]; i < iGrpEnd[gid]; ++i, ++n) {
            auto p = particles[i];
            assert(gid == p.group()); /* make sure it is really in group order */
            r = p.position() - tinyGroupTable[gid].rPot;
            for (j = 0; j < 3; ++j) {
                if (r[j] < -dHalf[j]) r[j] += dPeriod[j];
                else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
            }
            r -= tinyGroupTable[gid].rcen;
            r2 = dot(r, r);
            mrFree[n].fMass = p.mass();
            mrFree[n].dr = sqrtf(r2);
        }
        std::sort(mrFree, mrFree + n);
        fMass = mrFree[i].fMass;
        for (i = 1; i < n; ++i) {
            fMass += mrFree[i].fMass;
            mrFree[i].fMass = fMass;
        }
        mrIndex[gid] = mrFree - mr + n;
        /*
        ** Set the new mrFree to unallocated space.
        */
        mrFree = mr + mrIndex[gid];
    }
    /*
    ** Set up a table which defines the rootFunction of all enclosed
    ** group mass at a given dr. Also set up extra table needed to keep
    ** track of lower and upper bounds for root finding after all the
    ** saved mass-radius values.
    ** We will only need nRootFind entries in this table.
    */
    iRoot = 0;
    rootFunction = mrFree;
    bRemoteDone = (int*)(&rootFunction[nLocalGroups + 1]);
    for (gid = nLocalGroups + 1; gid < nGroups; ++gid) bRemoteDone[gid - (nLocalGroups + 1)] = 0;
    rootFindingTable = (RootFindingTable*)(&bRemoteDone[nGroups - (nLocalGroups + 1)]);
    for (gid = 0; gid <= nLocalGroups; ++gid) {
        if (mrIndex[gid] > mrIndex[gid - 1]) {
            rootFindingTable[iRoot].bConverged = 0;
            rootFindingTable[iRoot].gid = gid;
            rootFindingTable[iRoot].iter = 65536;
            rootFindingTable[iRoot].a = tinyGroupTable[gid].rHalf;
            rootFindingTable[iRoot].b = tinyGroupTable[gid].rMax;
            ++iRoot;
        }
        else {
            rootFunction[gid].dr = -1.0; /* marks no root finding needed */
        }
        rootFunction[gid].fMass = 0;
    }
    assert(iRoot == nRootFind);
    /*
    ** Now everything is set up for the root finding phases.
    */
    nMaxIter = 16;
    for (iter = 0; iter < nMaxIter; ++iter) {
        /*
        ** First calculate the local mass contained in a trial dr.
        */
        for (iRoot = 0; iRoot < nRootFind; ++iRoot) {
            gid = rootFindingTable[iRoot].gid;
            if (!rootFindingTable[iRoot].bConverged) {
                rootFunction[gid].dr = 0.5 * (rootFindingTable[iRoot].a + rootFindingTable[iRoot].b);
                /*
                ** Binary search for rootFunction[gid].dr
                */
                int ia = mrIndex[gid - 1];
                int ib = mrIndex[gid] - 1;
                while (ia < ib) {
                    int ic = (ia + ib + 1) / 2;
                    if (mr[ic].dr > rootFunction[gid].dr) ib = ic - 1;
                    else ia = ic;
                }
                assert(ia == ib);
                rootFunction[gid].fMass = (rootFunction[gid].dr > mr[ia].dr) ? mr[ia].fMass : 0;
            }
            else {
                rootFunction[gid].dr = -1.0; /* communicates to the remote that nothing needs to be done for this gid */
            }
        }
        /*
        ** Now scatter the enclosed local mass to remote groups at their trial dr.
        */
        mdl->CacheInitialize(CID_RM, std::make_shared<HelperRoot>(rootFunction, nLocalGroups + 1));
        for (gid = nLocalGroups + 1; gid < nGroups; ++gid) {
            if (bRemoteDone[gid - (nLocalGroups + 1)]) continue;
            auto pmr = mdl->Acquire<MassRadius>(CID_RM, ga[gid].id.iIndex, ga[gid].id.iPid);
            if (pmr->dr > 0) {
                /*
                ** Binary search for pmr->dr!
                */
                int ia = mrIndex[gid - 1];
                int ib = mrIndex[gid] - 1;
                while (ia < ib) {
                    int ic = (ia + ib + 1) / 2;
                    if (mr[ic].dr > pmr->dr) ib = ic - 1;
                    else ia = ic;
                }
                assert(ia == ib);
                pmr->fMass += (pmr->dr > mr[ia].dr) ? mr[ia].fMass : 0;
            }
            else {
                bRemoteDone[gid - (nLocalGroups + 1)] = 1; /* This will even bypass the mdlAcquire in the next iteration. */
            }
            mdl->Release(CID_RM, pmr);
        }
        mdl->FinishCache(CID_RM);
        /*
        ** Now the root functions are complete.
        */
        for (iRoot = 0; iRoot < nRootFind; ++iRoot) {
            if (!rootFindingTable[iRoot].bConverged) {
                gid = rootFindingTable[iRoot].gid;
                if (rootFunction[gid].fMass > (0.5 + 0.6 / ga[gid].nTotal) * tinyGroupTable[gid].fMass) rootFindingTable[iRoot].b = rootFunction[gid].dr;
                else if (rootFunction[gid].fMass < (0.5 - 0.6 / ga[gid].nTotal) * tinyGroupTable[gid].fMass) rootFindingTable[iRoot].a = rootFunction[gid].dr;
                else {
                    rootFindingTable[iRoot].bConverged = 1;
                    tinyGroupTable[gid].rHalf = rootFunction[gid].dr;
                    rootFindingTable[iRoot].iter = iter;
                }
            }
        }
    }
    /*
    ** Important to really make sure pkd->nLocalGroups is set correctly before output!
    */
    this->nLocalGroups = nLocalGroups;

    /*
    ** IA: We copy part of the stats to a even smaller structure that will survive later in memory
    **  such that it can be used for planting the BH seeds
    */
    if (veryTinyGroupTable != NULL) free(veryTinyGroupTable);
    veryTinyGroupTable = (VeryTinyGroupTable*)malloc((1 + nLocalGroups) * sizeof(VeryTinyGroupTable));
    for (gid = 1; gid <= nLocalGroups; gid++) {
        veryTinyGroupTable[gid].rPot[0] = tinyGroupTable[gid].rPot[0];
        veryTinyGroupTable[gid].rPot[1] = tinyGroupTable[gid].rPot[1];
        veryTinyGroupTable[gid].rPot[2] = tinyGroupTable[gid].rPot[2];
        veryTinyGroupTable[gid].fMass = tinyGroupTable[gid].fMass;
        veryTinyGroupTable[gid].nBH = tinyGroupTable[gid].nBH;
    }

    /*
    ** Clear ga so that later gravity calculations don't set the minimum potential again.
    ** This means we cannot directly output ga and could use this space for something else.
    */
    ga = NULL;
}
