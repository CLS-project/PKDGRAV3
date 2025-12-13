#ifndef WALK_HYDRO_H
#define WALK_HYDRO_H
#include "../TraversePST.h"
#include "basewalk.h"
#include "blitz/tinyvec2.h"
#include "core/simd.h"
#include "gravity/moments.h"
#include "gravity/vmoments.h"
#include "interaction/gravity.h"
#include "pkd.h"
#include "smooth/smooth.h"
#include "walk/stats.h"

class HydroWalkTraits {
public:
    struct name {
        static constexpr std::string_view value = "HydroWalk";
    };

    // For easier porting, we just duplicate the options,
    // but they should be logically separated probably
    struct Options {
        double dTheta; // needed for the opening criteria (could be removed)
        int iRoot1;
        int iRoot2;
        int nReps;
        float ballSizeLimit;
        float fBallFactor;
        int dofBallFactor;

        uint8_t uRungLo;
        uint8_t uRungHi;
    };

    struct output {
        STAT sLocal;
        STAT sActive;
        STAT sPart;
        STAT sPartNumAccess;
        STAT sPartMissRatio;
        STAT sCell;
        STAT sCellNumAccess;
        STAT sCellMissRatio;
        STAT sFlop;
#ifdef INSTRUMENT
        STAT sComputing;
        STAT sWaiting;
        STAT sSynchronizing;
#endif
#ifdef __linux__
        STAT sRSS;
        STAT sFreeMemory;
#endif
        double dFlopSingleCPU;
        double dFlopDoubleCPU;
        double dFlopSingleGPU;
        double dFlopDoubleGPU;
        uint64_t nActive;
        uint64_t nTilesTotal;
        uint64_t nTilesCPU;
    };
};

template <typename Interaction, typename Derived>
class HydroWalkImpl : public BaseTreeWalk<HydroWalkTraits, Interaction,
    HydroWalkImpl<Interaction, Derived>> {

    template <class T> friend class BaseTreeWalkInterface;
    using base = BaseTreeWalk<HydroWalkTraits, Interaction,
          HydroWalkImpl<Interaction, Derived>>;

public:
    explicit HydroWalkImpl(PST pst)
        : BaseTreeWalk<HydroWalkTraits, Interaction, HydroWalkImpl>(pst) {}

private:
    double dPartSum;

protected:
    void initCLAccumulatedVariablesImpl() {};

    void beforeCLProcessingImpl(treeStore::NodePointer k) {};

    enum outcome {
        KEEP_ON_CL,
        ADD_PARTICLES_PP_LIST,
        OPEN_BUCKET,
        OPEN_CELL,
        ACCEPT_MULTIPOLE,
        ACCEPT_LOCAL_EXPANSION,
        REMOVE_CELL,
    };

    // We could really try to optimize this for our case, but for now it is ok
    void openOutcomeBlockImpl(treeStore::NodePointer k, int n, clBlock &blk) {
        const float walk_min_multipole = 3;
        fmask T0, T1, T2, T3, T4, T6, T7;
        i32v P1, P2, P3, P4;
        fvec xc, yc, zc, dx, dy, dz, d2, diCrit, cOpen, cOpen2, d2Open, mink2,
             minbnd2, fourh2;
        fvec distk2, distc2;
        fvec box1xMin, box1xMax, box1yMin, box1yMax, box1zMin, box1zMax,
             box2xMin, box2xMax, box2yMin, box2yMax, box2zMin, box2zMax;
        int i;
        i32v iOpen, iOpenA, iOpenB;
        fvec k_xCenter, k_yCenter, k_zCenter, k_xMax, k_yMax, k_zMax;
        fvec k_xMinBnd, k_yMinBnd, k_zMinBnd, k_xMaxBnd, k_yMaxBnd, k_zMaxBnd;
        fvec k_x, k_y, k_z, k_bMax, k_Open;
        fmask k_notgrp;
        auto k_r = k->position();

        assert(k->mass() > 0.0f);

        diCrit = 1.0f / this->options->dTheta;

        auto kbnd = k->bound();
        k_xMinBnd = kbnd.lower(0);
        k_yMinBnd = kbnd.lower(1);
        k_zMinBnd = kbnd.lower(2);
        k_xMaxBnd = kbnd.upper(0);
        k_yMaxBnd = kbnd.upper(1);
        k_zMaxBnd = kbnd.upper(2);
        k_xCenter = kbnd.center(0);
        k_yCenter = kbnd.center(1);
        k_zCenter = kbnd.center(2);
        k_xMax = kbnd.apothem(0);
        k_yMax = kbnd.apothem(1);
        k_zMax = kbnd.apothem(2);
        k_x = k_r[0];
        k_y = k_r[1];
        k_z = k_r[2];
        k_bMax = k->bMax();
        k_notgrp = cvt_fvec(i32v(k->is_group())) == 0.0;
        k_Open = 1.5f * k_bMax * diCrit;
        const auto &SPHbob = k->have_BOB() ? k->BOB() : SPHBOB();
        fvec k_fBoBr2 = SPHbob.fBoBr * SPHbob.fBoBr;
        fvec blk_fBoBr2 = 0.0f;
        distk2 = HUGE_VALF;
        distc2 = HUGE_VALF;

        fmask bIntersectSource, bIntersectSink;
        bIntersectSource = 0;
        bIntersectSink = 0;

        n = (n + fvec::mask()) / fvec::width(); // Now number of blocks
        for (i = 0; i < n; ++i) {
            fourh2 = blk.fourh2.v[i];

            distk2 = 0.0f;
            dx = SPHbob.fBoBCenter[0] - blk.xCenter.v[i] - blk.xOffset.v[i] -
                 blk.xMax.v[i];
            distk2 += maskz_mov(dx > 0, dx * dx);
            dx = blk.xCenter.v[i] + blk.xOffset.v[i] - blk.xMax.v[i] -
                 SPHbob.fBoBCenter[0];
            distk2 += maskz_mov(dx > 0, dx * dx);

            dx = SPHbob.fBoBCenter[1] - blk.yCenter.v[i] - blk.yOffset.v[i] -
                 blk.yMax.v[i];
            distk2 += maskz_mov(dx > 0, dx * dx);
            dx = blk.yCenter.v[i] + blk.yOffset.v[i] - blk.yMax.v[i] -
                 SPHbob.fBoBCenter[1];
            distk2 += maskz_mov(dx > 0, dx * dx);

            dx = SPHbob.fBoBCenter[2] - blk.zCenter.v[i] - blk.zOffset.v[i] -
                 blk.zMax.v[i];
            distk2 += maskz_mov(dx > 0, dx * dx);
            dx = blk.zCenter.v[i] + blk.zOffset.v[i] - blk.zMax.v[i] -
                 SPHbob.fBoBCenter[2];
            distk2 += maskz_mov(dx > 0, dx * dx);
            fmask bIntersectSourceBall = distk2 < k_fBoBr2;

            box1xMin = SPHbob.fBoBMin[0];
            box1xMax = SPHbob.fBoBMax[0];
            box1yMin = SPHbob.fBoBMin[1];
            box1yMax = SPHbob.fBoBMax[1];
            box1zMin = SPHbob.fBoBMin[2];
            box1zMax = SPHbob.fBoBMax[2];
            box2xMin = blk.xCenter.v[i] + blk.xOffset.v[i] - blk.xMax.v[i];
            box2xMax = blk.xCenter.v[i] + blk.xOffset.v[i] + blk.xMax.v[i];
            box2yMin = blk.yCenter.v[i] + blk.yOffset.v[i] - blk.yMax.v[i];
            box2yMax = blk.yCenter.v[i] + blk.yOffset.v[i] + blk.yMax.v[i];
            box2zMin = blk.zCenter.v[i] + blk.zOffset.v[i] - blk.zMax.v[i];
            box2zMax = blk.zCenter.v[i] + blk.zOffset.v[i] + blk.zMax.v[i];
            fmask bIntersectSourceBox =
                (box1xMin < box2xMax) & (box2xMin < box1xMax) &
                (box1yMin < box2yMax) & (box2yMin < box1yMax) &
                (box1zMin < box2zMax) & (box2zMin < box1zMax);

            bIntersectSource = bIntersectSourceBall & bIntersectSourceBox;
            blk_fBoBr2 = blk.fBoBr.v[i] * blk.fBoBr.v[i];
            distc2 = 0.0f;
            dx = k_xMinBnd - blk.fBoBxCenter.v[i] - blk.xOffset.v[i];
            distc2 += maskz_mov(dx > 0, dx * dx);
            dx = blk.fBoBxCenter.v[i] + blk.xOffset.v[i] - k_xMaxBnd;
            distc2 += maskz_mov(dx > 0, dx * dx);

            dx = k_yMinBnd - blk.fBoByCenter.v[i] - blk.yOffset.v[i];
            distc2 += maskz_mov(dx > 0, dx * dx);
            dx = blk.fBoByCenter.v[i] + blk.yOffset.v[i] - k_yMaxBnd;
            distc2 += maskz_mov(dx > 0, dx * dx);

            dx = k_zMinBnd - blk.fBoBzCenter.v[i] - blk.zOffset.v[i];
            distc2 += maskz_mov(dx > 0, dx * dx);
            dx = blk.fBoBzCenter.v[i] + blk.zOffset.v[i] - k_zMaxBnd;
            distc2 += maskz_mov(dx > 0, dx * dx);
            fmask bIntersectSinkBall = distc2 < blk_fBoBr2;

            box1xMin = k_xMinBnd;
            box1xMax = k_xMaxBnd;
            box1yMin = k_yMinBnd;
            box1yMax = k_yMaxBnd;
            box1zMin = k_zMinBnd;
            box1zMax = k_zMaxBnd;
            box2xMin = blk.fBoBxMin.v[i] + blk.xOffset.v[i];
            box2xMax = blk.fBoBxMax.v[i] + blk.xOffset.v[i];
            box2yMin = blk.fBoByMin.v[i] + blk.yOffset.v[i];
            box2yMax = blk.fBoByMax.v[i] + blk.yOffset.v[i];
            box2zMin = blk.fBoBzMin.v[i] + blk.zOffset.v[i];
            box2zMax = blk.fBoBzMax.v[i] + blk.zOffset.v[i];
            fmask bIntersectSinkBox =
                (box1xMin < box2xMax) & (box2xMin < box1xMax) &
                (box1yMin < box2yMax) & (box2yMin < box1yMax) &
                (box1zMin < box2zMax) & (box2zMin < box1zMax);

            bIntersectSink = bIntersectSinkBall & bIntersectSinkBox;

            xc = blk.x.v[i] + blk.xOffset.v[i];
            yc = blk.y.v[i] + blk.yOffset.v[i];
            zc = blk.z.v[i] + blk.zOffset.v[i];
            dx = k_x - xc;
            dy = k_y - yc;
            dz = k_z - zc;
            d2 = dx * dx + dy * dy + dz * dz;
            cOpen = blk.cOpen.v[i];
            cOpen2 = cOpen * cOpen;
            d2Open = cOpen + k_Open;
            d2Open = d2Open * d2Open;

            dx = abs(xc - k_xCenter) - k_xMax;
            dy = abs(yc - k_yCenter) - k_yMax;
            dz = abs(zc - k_zCenter) - k_zMax;

            dx = maskz_mov(dx > 0, dx);
            dy = maskz_mov(dy > 0, dy);
            dz = maskz_mov(dz > 0, dz);
            mink2 = dx * dx + dy * dy + dz * dz;
            minbnd2 = 0.0f;

            dx =
                k_xMinBnd - blk.xCenter.v[i] - blk.xOffset.v[i] - blk.xMax.v[i];
            minbnd2 += maskz_mov(dx > 0, dx * dx);
            dx =
                blk.xCenter.v[i] + blk.xOffset.v[i] - blk.xMax.v[i] - k_xMaxBnd;
            minbnd2 += maskz_mov(dx > 0, dx * dx);

            dx =
                k_yMinBnd - blk.yCenter.v[i] - blk.yOffset.v[i] - blk.yMax.v[i];
            minbnd2 += maskz_mov(dx > 0, dx * dx);
            dx =
                blk.yCenter.v[i] + blk.yOffset.v[i] - blk.yMax.v[i] - k_yMaxBnd;
            minbnd2 += maskz_mov(dx > 0, dx * dx);

            dx =
                k_zMinBnd - blk.zCenter.v[i] - blk.zOffset.v[i] - blk.zMax.v[i];
            minbnd2 += maskz_mov(dx > 0, dx * dx);
            dx =
                blk.zCenter.v[i] + blk.zOffset.v[i] - blk.zMax.v[i] - k_zMaxBnd;
            minbnd2 += maskz_mov(dx > 0, dx * dx);

            T0 = blk.m.v[i] > 0.0f;
            T1 = (d2 > d2Open) & (minbnd2 > fourh2) & ~bIntersectSource &
                 ~bIntersectSink;
            T2 = cvt_fvec(blk.iLower.v[i]) == 0.0;
            T3 = (walk_min_multipole > cvt_fvec(blk.nc.v[i])) |
                 (mink2 <= cOpen2);
            T4 = (minbnd2 > fourh2) & ~bIntersectSource & ~bIntersectSink;
            T6 = cOpen > k_Open;
            T7 = k_notgrp;
            iOpenA = mask_mov(i32v(static_cast<int>(OPEN_CELL)), T2, i32v(static_cast<int>(ADD_PARTICLES_PP_LIST)));
            iOpenB = mask_mov(mask_mov(iOpenA, T4, i32v(static_cast<int>(ACCEPT_MULTIPOLE))), T3,
                              iOpenA);
            P1 = mask_mov(i32v(static_cast<int>(OPEN_CELL)), T2, i32v(static_cast<int>(OPEN_BUCKET)));
            P2 = mask_mov(iOpenB, T7, i32v(static_cast<int>(KEEP_ON_CL)));
            P3 = mask_mov(P2, T6, P1);
            P4 = mask_mov(P3, T1, i32v(static_cast<int>(ACCEPT_LOCAL_EXPANSION)));
            iOpen = mask_mov(i32v(static_cast<int>(REMOVE_CELL)), T0, P4);
            blk.iOpen.v[i] = iOpen;
        }

        double dFlopOpening = COST_FLOP_OPEN * n;
        this->walkStats.dFlop += dFlopOpening;
        this->walkStats.dFlopSingleCPU += dFlopOpening;
    };

    void selectOutcomeImpl(treeStore::NodePointer k) {
        blitz::TinyVector<float, 3> fOffset;
        blitz::TinyVector<double, 3> r, dx, c_r, dOffset;
        treeStore::NodePointer c(this->pkd->tree);
        float fMass, fSoft;
        int id, iCell, iSib, iCheckCell, iCheckLower, iCidPart;
        int j, pj;
        int jTile;
        auto k_r = k->position();
        for (auto &tile : this->pkd->cl) {
            auto nBlocks = tile.count() / tile.width;
            for (auto iBlock = 0; iBlock <= nBlocks; ++iBlock) {
                int n = iBlock < nBlocks ? tile.width
                        : tile.count() - nBlocks * tile.width;
                auto &blk = tile[iBlock];
                for (jTile = 0; jTile < n; ++jTile) {
                    switch (blk.iOpen[jTile]) {
                    case KEEP_ON_CL:
                        /*
                        ** This checkcell stays on the checklist.
                        */
                        this->clSib->append(blk, jTile);
                        break;
                    case ADD_PARTICLES_PP_LIST:
                        /*
                        ** This checkcell's particles are added
                        to the P-P list.
                        */
                        iCheckCell = blk.iCell[jTile];
                        id = blk.idCell[jTile];
                        if (iCheckCell < 0) {
                            pj = -1 - iCheckCell;
                            iCidPart = blk.iCache[jTile] == CID_CELL
                                       ? CID_PARTICLE
                                       : CID_PARTICLE2;
                            auto p = (id == this->pkd->Self())
                                     ? this->pkd->particles[pj]
                                     : this->pkd->particles
                                     [this->pkd->mdl->template Fetch<
                                                               PARTICLE>(iCidPart, pj, id)];
                            r = p.position();
                            if (!this->bReferenceFound) {
                                this->bReferenceFound = 1;
                                this->interaction.setReferenceIL(r);
                            }
                            auto v = p.velocity();
                            this->interaction.addParticleIL(
                                r[0] + blk.xOffset[jTile],
                                r[1] + blk.yOffset[jTile],
                                r[2] + blk.zOffset[jTile], blk.m[jTile],
                                blk.fourh2[jTile], &p);
                        }
                        else {
                            assert(id >= 0);
                            c = (id == this->pkd->Self())
                                ? this->pkd->tree[iCheckCell]
                                : this->pkd->tree[this->pkd->mdl
                                                  ->template Fetch<KDN>(
                                                      blk.iCache[jTile],
                                                      iCheckCell, id)];
                            iCidPart = blk.iCache[jTile] == CID_CELL
                                       ? CID_PARTICLE
                                       : CID_PARTICLE2;
                            if (!this->bReferenceFound) {
                                this->bReferenceFound = 1;
                                auto p =
                                    (id == this->pkd->Self())
                                    ? this->pkd->particles[c->lower()]
                                    : this->pkd->particles
                                    [this->pkd->mdl->template Fetch<
                                                              PARTICLE>(iCidPart,
                                                                        c->lower(), id)];
                                r = p.position();
                                this->interaction.setReferenceIL(r);
                            }
                            for (pj = c->lower(); pj <= c->upper(); ++pj) {
                                auto p =
                                    (id == this->pkd->Self())
                                    ? this->pkd->particles[pj]
                                    : this->pkd->particles
                                    [this->pkd->mdl->template Fetch<
                                                              PARTICLE>(iCidPart, pj, id)];
                                fMass = p.mass();
                                fSoft = p.soft();
                                r = p.position();
                                this->interaction.addParticleIL(
                                    r[0] + blk.xOffset[jTile],
                                    r[1] + blk.yOffset[jTile],
                                    r[2] + blk.zOffset[jTile], fMass,
                                    4 * fSoft * fSoft, &p);
                            }
                        }
                        break;
                    case OPEN_BUCKET:
                        /*
                        ** Now I am trying to open a bucket,
                        which means I add each particle back on
                        the
                        ** checklist with a cell size of zero.
                        */
                        iCheckCell = blk.iCell[jTile];
                        assert(iCheckCell >= 0);
                        id = blk.idCell[jTile];
                        fOffset[0] = blk.xOffset[jTile];
                        fOffset[1] = blk.yOffset[jTile];
                        fOffset[2] = blk.zOffset[jTile];
                        assert(id >= 0);
                        c = (id == this->pkd->Self())
                            ? this->pkd->tree[iCheckCell]
                            : this->pkd->tree
                            [this->pkd->mdl->template Fetch<KDN>(
                                                          blk.iCache[jTile], iCheckCell, id)];
                        iCidPart = blk.iCache[jTile] == CID_CELL
                                   ? CID_PARTICLE
                                   : CID_PARTICLE2;
                        for (pj = c->lower(); pj <= c->upper(); ++pj) {
                            auto p = (id == this->pkd->Self())
                                     ? this->pkd->particles[pj]
                                     : this->pkd->particles
                                     [this->pkd->mdl->template Fetch<
                                                               PARTICLE>(iCidPart, pj, id)];
                            r = p.position();
                            fMass = p.mass();
                            fSoft = p.soft();
                            float fBallFactor = (this->options->dofBallFactor) ? this->options->fBallFactor : 1.0f;
                            float limitedBallSize = std::min(this->options->ballSizeLimit,fBallFactor * p.ball());
                            this->clNew->append(
                                blk.iCache[jTile], id, -1 - pj, 0, 0, 0, 0, 1,
                                0.0, fMass, 4.0f * fSoft * fSoft, r, fOffset,
                                Bound(r, r), SPHBOB(r, limitedBallSize));
                        }
                        break;
                    case OPEN_CELL:
                        /*
                        ** Open the cell.
                        ** We could do a prefetch here for
                        non-local
                        ** cells.
                        */
                        iCheckCell = blk.iCell[jTile];
                        assert(iCheckCell >= 0);
                        iCheckLower = blk.iLower[jTile];
                        assert(iCheckLower > 0);

                        fOffset[0] = blk.xOffset[jTile];
                        fOffset[1] = blk.yOffset[jTile];
                        fOffset[2] = blk.zOffset[jTile];

                        this->addChild(this->pkd, blk.iCache[jTile],
                                       *this->clNew, blk.iLower[jTile],
                                       blk.idLower[jTile], fOffset);
                        this->addChild(this->pkd, blk.iCache[jTile],
                                       *this->clNew, blk.iUpper[jTile],
                                       blk.idUpper[jTile], fOffset);

                        break;
                    case ACCEPT_MULTIPOLE:
                        break;
                    case ACCEPT_LOCAL_EXPANSION:
                        break;
                    case REMOVE_CELL:
                        /*
                        ** This checkcell is removed from the
                        checklist since it has zero/negative mass.
                        */
                        break;
                    default:
                        assert(0);
                    }
                }
            }
        }
    };

    bool testAddNodeCLImpl(treeStore::NodePointer k) {
        return k->is_marked();
    }

    int interact_wrapperImpl(treeStore::NodePointer k, int iRoot, int iRoot2) {
        auto nActiveInt = this->interaction.interact(k, iRoot, iRoot2);

        if (nActiveInt) {
            // This is NOT a good approximation in this case, as a particle
            // does not interact with all the IL (most of them may be outside
            // of its smoothing kernel)
            dPartSum += nActiveInt * this->interaction.countIL();
        }
        return nActiveInt;
    }

    void initStatsImpl() {
        dPartSum = 0.0;
        base::initStatsImpl();
    }

    // We reuse the gravity statistics, but maybe this can be
    // part of the default implementation?
    void computeStatsImpl() {
        // Flop and tile count
        this->out->sFlop.dSum = 1e-9 * this->walkStats.dFlop;
        this->out->dFlopSingleCPU = 1e-9 * this->walkStats.dFlopSingleCPU;
        this->out->dFlopDoubleCPU = 1e-9 * this->walkStats.dFlopDoubleCPU;
        this->out->dFlopSingleGPU = 1e-9 * this->walkStats.dFlopSingleGPU;
        this->out->dFlopDoubleGPU = 1e-9 * this->walkStats.dFlopDoubleGPU;
        this->out->nTilesTotal = this->walkStats.nTilesTotal;
        this->out->nTilesCPU = this->walkStats.nTilesCPU;
        this->out->nActive = this->walkStats.nActive;

        double dActive = static_cast<double>(this->walkStats.nActive);
        // Set dSum for the rest of statistics
        if (this->walkStats.nActive) {
            this->out->sCell.dSum = 0.0;
            this->out->sPart.dSum = dPartSum / dActive;
            this->out->sPartNumAccess.dSum =
                this->pkd->mdl->NumAccess(CID_PARTICLE) / dActive;
            this->out->sCellNumAccess.dSum =
                this->pkd->mdl->NumAccess(CID_CELL) / dActive;
        }
        else {
            this->out->sCell.dSum = 0;
            this->out->sPart.dSum = 0;
            this->out->sPartNumAccess.dSum = 0;
            this->out->sCellNumAccess.dSum = 0;
        }
        this->out->sCellMissRatio.dSum =
            100.0 * this->pkd->mdl->MissRatio(CID_CELL);
        this->out->sPartMissRatio.dSum =
            100.0 * this->pkd->mdl->MissRatio(CID_PARTICLE);

        this->out->sComputing.dSum = this->pkd->mdl->TimeComputing();
        this->out->sWaiting.dSum = this->pkd->mdl->TimeWaiting();
        this->out->sSynchronizing.dSum = this->pkd->mdl->TimeSynchronizing();

        this->out->sLocal.dSum = this->pkd->Local();
        this->out->sActive.dSum = dActive;

#ifdef __linux__
        FILE *fp;
        char buffer[512], *save, *f, *v;
        fp = fopen("/proc/self/stat", "r");
        if (fp != NULL) {
            if (fgets(buffer, sizeof(buffer), fp) == NULL)
                buffer[0] = '\0';
            fclose(fp);
            f = strtok_r(buffer, " ", &save);
            for (auto i = 0; i <= 36 && f; i++) {
                switch (i) {
                case 23:
#ifdef HAVE_GETPAGESIZE
                    this->out->sRSS.dSum =
                        atol(f) * 1.0 * getpagesize() / 1024 / 1024 / 1024;
#else
                    this->out->sRSS.dSum = atol(f) * 1.0 / 1024 / 1024;
#endif
                    break;
                default:
                    break;
                }
                f = strtok_r(NULL, " ", &save);
            }
        }
        else
            this->out->sRSS.dSum = 0;

        fp = fopen("/proc/meminfo", "r");
        this->out->sFreeMemory.dSum = 0;
        if (fp != NULL) {
            while (fgets(buffer, sizeof(buffer), fp)) {
                f = strtok_r(buffer, ":", &save);
                v = strtok_r(NULL, " ", &save);
                if (strcmp(f, "MemFree") == 0 || strcmp(f, "Buffers") == 0 ||
                        strcmp(f, "Cached") == 0) {
                    this->out->sFreeMemory.dSum -= atol(v) * 1.0 / 1024 / 1024;
                }
            }
            fclose(fp);
        }
#endif
        /*
        ** Initialize statistics tracking variables assuming dSum is set.
        */
#ifdef __linux__
        this->initStat(&this->out->sRSS, this->pkd->Self());
        this->initStat(&this->out->sFreeMemory, this->pkd->Self());
#endif
        this->initStat(&this->out->sLocal, this->pkd->Self());
        this->initStat(&this->out->sActive, this->pkd->Self());
        this->initStat(&this->out->sPart, this->pkd->Self());
        this->initStat(&this->out->sPartNumAccess, this->pkd->Self());
        this->initStat(&this->out->sPartMissRatio, this->pkd->Self());
        this->initStat(&this->out->sCell, this->pkd->Self());
        this->initStat(&this->out->sCellNumAccess, this->pkd->Self());
        this->initStat(&this->out->sCellMissRatio, this->pkd->Self());
        this->initStat(&this->out->sFlop, this->pkd->Self());
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
        this->initStat(&this->out->sComputing, this->pkd->Self());
        this->initStat(&this->out->sWaiting, this->pkd->Self());
        this->initStat(&this->out->sSynchronizing, this->pkd->Self());
#endif
        /*
        ** If there were no actives then we need to zero statistics counts
        ** for those quantities which are based on per active particles.
        */
        if (!this->out->nActive) {
            this->out->sPart.n = 0;
            this->out->sPartNumAccess.n = 0;
            this->out->sPartMissRatio.n = 0;
            this->out->sCell.n = 0;
            this->out->sCellNumAccess.n = 0;
            this->out->sCellMissRatio.n = 0;
        }
    }

public:
    virtual int CombineImpl(void *vout, void *vout2) {
        auto tmp = static_cast<typename HydroWalkTraits::output *>(vout2);
        auto out = static_cast<typename HydroWalkTraits::output *>(vout);
#ifdef __linux__
        this->combStat(&this->out->sRSS, &tmp->sRSS);
        this->combStat(&this->out->sFreeMemory, &tmp->sFreeMemory);
#endif
        this->combStat(&this->out->sLocal, &tmp->sLocal);
        this->combStat(&this->out->sActive, &tmp->sActive);
        this->combStat(&this->out->sPart, &tmp->sPart);
        this->combStat(&this->out->sPartNumAccess, &tmp->sPartNumAccess);
        this->combStat(&this->out->sPartMissRatio, &tmp->sPartMissRatio);
        this->combStat(&this->out->sCell, &tmp->sCell);
        this->combStat(&this->out->sCellNumAccess, &tmp->sCellNumAccess);
        this->combStat(&this->out->sCellMissRatio, &tmp->sCellMissRatio);
        this->combStat(&this->out->sFlop, &tmp->sFlop);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
        this->combStat(&this->out->sComputing, &tmp->sComputing);
        this->combStat(&this->out->sWaiting, &tmp->sWaiting);
        this->combStat(&this->out->sSynchronizing, &tmp->sSynchronizing);
#endif
        /*
        ** and the number of actives processed in this gravity call.
        */
        this->out->nActive += tmp->nActive;
        this->out->dFlopSingleCPU += tmp->dFlopSingleCPU;
        this->out->dFlopDoubleCPU += tmp->dFlopDoubleCPU;
        this->out->dFlopSingleGPU += tmp->dFlopSingleGPU;
        this->out->dFlopDoubleGPU += tmp->dFlopDoubleGPU;
        this->out->nTilesTotal += tmp->nTilesTotal;
        this->out->nTilesCPU += tmp->nTilesCPU;

        return sizeof(out);
    };
};
template <typename Interaction>
class HydroWalk : public HydroWalkImpl<Interaction, HydroWalk<Interaction>> {
public:
    explicit HydroWalk(PST pst)
        : HydroWalkImpl<Interaction, HydroWalk<Interaction>>(pst) {}
};
#endif
