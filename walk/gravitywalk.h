#ifndef WALK_GRAVITY_H
#define WALK_GRAVITY_H
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

class GravityWalkTraits {
public:
    struct name {
        static constexpr std::string_view value = "GravityWalk";
    };

    // For easier porting, we just duplicate the options,
    // but they should be logically separated probably
    using Options = GravityInteractionTraits::Options;

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
        uint64_t nRung[IRUNGMAX + 1];
        uint64_t nTilesTotal;
        uint64_t nTilesCPU;
    };
};

template <typename Interaction, typename Derived>
class GravityWalkImpl
    : public BaseTreeWalk<GravityWalkTraits, Interaction,
      GravityWalkImpl<Interaction, Derived>> {

    // Note here to but the template outside to be "friend" with all the
    // template specializations of the base class
    template <class T> friend class BaseTreeWalkInterface;
    using base = BaseTreeWalk<GravityWalkTraits, Interaction,
          GravityWalkImpl<Interaction, Derived>>;

public:
    explicit GravityWalkImpl(PST pst)
        : BaseTreeWalk<GravityWalkTraits, Interaction, GravityWalkImpl>(pst) {}

    int ServiceImpl(PST pst, void *vin, int nIn, void *vout, int nOut) {
        this->pkd = pst->plcl->pkd;
        this->options = static_cast<GravityWalkTraits::Options *>(vin);
        this->out = static_cast<GravityWalkTraits::output *>(vout);
        assert(sizeof(*this->out) == nOut);
        assert(sizeof(*this->options) == nIn);
        // We reuse the same options for the interaction
        this->interaction.setup(this->pkd, this->options, &this->walkStats);

        // options->copyToDevice(); // or something similar
        this->setupCache();

        this->initStats();

        this->walk();

        this->computeStats();

        this->pkd->mdl->FinishCache(CID_PARTICLE);
        return sizeof(*this->out);
    }

private:
    blitz::TinyVector<float, 3> a;
    float dirLsum, normLsum;
    FMOMR monoPole;
    LOCR L;
    SMX smx;
    SMF smf;
    double dFlop;
    int nWpPending;

protected:
    void initCLAccumulatedVariablesImpl() {
        momClearFmomr(&monoPole);
        momClearLocr(&L);
        dirLsum = 0;
        normLsum = 0;
        a = 0;
    };

    void beforeCLProcessingImpl(treeStore::NodePointer k) {
#ifdef USE_SIMD_FMM
        a = (this->options->ts.bGravStep) ? k->acceleration() : 0.0f;
#else
        float imaga;
        if (this->options->ts.bGravStep) {
            a = k->acceleration();
            imaga = 1.0 / sqrtf(blitz::dot(a, a));
        }
#endif
#ifdef USE_SIMD_FMM
        this->interaction.clearILL();
        // this->pkd->ill.clear();
#endif
    };
    // This would be a good moment to eliminate the magic numbers
    // and use an enum instead..
    enum outcome {
        KEEP_ON_CL,
        ADD_PARTICLES_PP_LIST,
        OPEN_BUCKET,
        OPEN_CELL,
        ACCEPT_MULTIPOLE,
        ACCEPT_LOCAL_EXPANSION,
        REMOVE_CELL,
    };

    // Mock opening criteria, only for compiling tests
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
            iOpenA = mask_mov(i32v(OPEN_CELL), T2, i32v(ADD_PARTICLES_PP_LIST));
            iOpenB = mask_mov(mask_mov(iOpenA, T4, i32v(ACCEPT_MULTIPOLE)), T3,
                              iOpenA);
            P1 = mask_mov(i32v(OPEN_CELL), T2, i32v(OPEN_BUCKET));
            P2 = mask_mov(iOpenB, T7, i32v(0));
            P3 = mask_mov(P2, T6, P1);
            P4 = mask_mov(P3, T1, i32v(ACCEPT_LOCAL_EXPANSION));
            iOpen = mask_mov(i32v(REMOVE_CELL), T0, P4);
            blk.iOpen.v[i] = iOpen;
        }

        double dFlop = COST_FLOP_OPEN * n;
        dFlop += dFlop;
        this->walkStats.dFlopSingleCPU += dFlop;
    };

    void selectOutcomeImpl(clBlock blk, int outcome, int jTile) {
        blitz::TinyVector<float, 3> fOffset;
        blitz::TinyVector<double, 3> r, dx, k_r, c_r, dOffset;
        treeStore::NodePointer c(this->pkd->tree);
        float fMass, fSoft;
        int id, iCell, iSib, iCheckCell, iCheckLower, iCidPart;
        int j, pj;
        switch (outcome) {
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
                iCidPart = blk.iCache[jTile] == CID_CELL ? CID_PARTICLE
                           : CID_PARTICLE2;
                auto p =
                    (id == this->pkd->Self())
                    ? this->pkd->particles[pj]
                    : this->pkd->particles[this->pkd->mdl->template Fetch<
                                                                    PARTICLE>(iCidPart, pj, id)];
                r = p.position();
                if (!this->bReferenceFound) {
                    this->bReferenceFound = 1;
                    this->interaction.setReferenceIL(r);
                }
                auto v = p.velocity();
                this->pkd->ilp.append(
                    r[0] + blk.xOffset[jTile], r[1] + blk.yOffset[jTile],
                    r[2] + blk.zOffset[jTile], blk.m[jTile], blk.fourh2[jTile],
                    v[0], v[1], v[2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, p.rung(),
                    0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
            }
            else {
                assert(id >= 0);
                auto c =
                    (id == this->pkd->Self())
                    ? this->pkd->tree[iCheckCell]
                    : this->pkd->tree[this->pkd->mdl->template Fetch<KDN>(
                                                                   blk.iCache[jTile], iCheckCell, id)];
                iCidPart = blk.iCache[jTile] == CID_CELL ? CID_PARTICLE
                           : CID_PARTICLE2;
                if (!this->bReferenceFound) {
                    this->bReferenceFound = 1;
                    auto p =
                        (id == this->pkd->Self())
                        ? this->pkd->particles[c->lower()]
                        : this->pkd
                        ->particles[this->pkd->mdl->template Fetch<
                                                             PARTICLE>(iCidPart, c->lower(), id)];
                    r = p.position();
                    this->interaction.setReferenceIL(r);
                }
                for (pj = c->lower(); pj <= c->upper(); ++pj) {
                    auto p =
                        (id == this->pkd->Self())
                        ? this->pkd->particles[pj]
                        : this->pkd
                        ->particles[this->pkd->mdl->template Fetch<
                                                             PARTICLE>(iCidPart, pj, id)];
                    fMass = p.mass();
                    fSoft = p.soft();
                    r = p.position();
                    auto v = p.velocity();
                    this->pkd->ilp.append(
                        r[0] + blk.xOffset[jTile], r[1] + blk.yOffset[jTile],
                        r[2] + blk.zOffset[jTile], fMass, 4 * fSoft * fSoft,
                        v[0], v[1], v[2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0,
                        p.rung(), 0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                        0.0f);
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
                : this->pkd->tree[this->pkd->mdl->template Fetch<KDN>(
                                                               blk.iCache[jTile], iCheckCell, id)];
            iCidPart =
                blk.iCache[jTile] == CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
            for (pj = c->lower(); pj <= c->upper(); ++pj) {
                auto p =
                    (id == this->pkd->Self())
                    ? this->pkd->particles[pj]
                    : this->pkd->particles[this->pkd->mdl->template Fetch<
                                                                    PARTICLE>(iCidPart, pj, id)];
                r = p.position();
                fMass = p.mass();
                fSoft = p.soft();
                // assert(0);
                this->clNew->append(blk.iCache[jTile], id, -1 - pj, 0, 0, 0, 0,
                                    1, 0.0, fMass, 4.0f * fSoft * fSoft, r,
                                    fOffset, Bound(r, r), SPHBOB(r, 0.0));
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

            this->addChild(this->pkd, blk.iCache[jTile], *this->clNew,
                           blk.iLower[jTile], blk.idLower[jTile], fOffset);
            this->addChild(this->pkd, blk.iCache[jTile], *this->clNew,
                           blk.iUpper[jTile], blk.idUpper[jTile], fOffset);

            break;
        case ACCEPT_MULTIPOLE:
            /*
            ** Accept multipole!
            ** Interact += Moment(c);
            */
            iCheckCell = blk.iCell[jTile];
            assert(iCheckCell >= 0);
            id = blk.idCell[jTile];
            c = (id == this->pkd->Self())
                ? this->pkd->tree[iCheckCell]
                : this->pkd->tree[this->pkd->mdl->template Fetch<KDN>(
                                                               blk.iCache[jTile], iCheckCell, id)];
            r[0] = blk.x[jTile] + blk.xOffset[jTile];
            r[1] = blk.y[jTile] + blk.yOffset[jTile];
            r[2] = blk.z[jTile] + blk.zOffset[jTile];
            if (!this->bReferenceFound) {
                this->bReferenceFound = 1;
                this->interaction.setReferenceIL(r);
            }
            this->interaction.addCellIL(r, &c->moment(), c->bMax());
            // pkd->ilc.append(r[0], r[1], r[2], &c->moment(), c->bMax());

            break;
        case ACCEPT_LOCAL_EXPANSION:
            /*
            ** Local expansion accepted!
            */
            iCheckCell = blk.iCell[jTile];
            if (iCheckCell < 0) {
                fOffset[0] = blk.xOffset[jTile];
                fOffset[1] = blk.yOffset[jTile];
                fOffset[2] = blk.zOffset[jTile];
                dx[0] = k_r[0] - (blk.x[jTile] + blk.xOffset[jTile]);
                dx[1] = k_r[1] - (blk.y[jTile] + blk.yOffset[jTile]);
                dx[2] = k_r[2] - (blk.z[jTile] + blk.zOffset[jTile]);
#ifdef USE_SIMD_FMM
                monoPole.m = blk.m[jTile];

                this->interaction.addLocalIL(dx, &monoPole, 1.0);
                // pkd->ill.append((float)dx[0], (float)dx[1], (float)dx[2],
                // &monoPole,  1.0);
#else
                double d2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                double dir = 1.0 / sqrt(d2);
                /* monoPole.m = blk.m[jTile];*/
                /* *pdFlop +=
                momLocrAddFmomr5cm(&L,&monoPole,0.0,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);*/
                double tax, tay, taz;
                dFlop = momLocrAddMono5(&L, blk.m[jTile], dir, dx[0], dx[1],
                                        dx[2], &tax, &tay, &taz);
                this->dFlop += dFlop;
                this->walkStats.dFlopDoubleCPU += dFlop;
                if (this->options->ts.bGravStep) {
                    float adotai = a[0] * tax + a[1] * tay + a[2] * taz;
                    if (adotai > 0) {
                        adotai *= imaga;
                        dirLsum += dir * adotai * adotai;
                        normLsum += adotai * adotai;
                    }
                }
#endif
            }
            else {
                id = blk.idCell[jTile];
                dOffset[0] = blk.xOffset[jTile];
                dOffset[1] = blk.yOffset[jTile];
                dOffset[2] = blk.zOffset[jTile];
                c = (id == this->pkd->Self())
                    ? this->pkd->tree[iCheckCell]
                    : this->pkd->tree[this->pkd->mdl->template Fetch<KDN>(
                                                                   blk.iCache[jTile], iCheckCell, id)];
                c_r = c->position();
#ifdef USE_SIMD_FMM
                dx = k_r - (c_r + dOffset);
                this->interaction.addLocalIL(dx, &c->moment(), c->bMax());
                // pkd->ill.append((float)dx[0], (float)dx[1], (float)dx[2],
                //                 &c->moment(), c->bMax());
#else
                dx = k_r - (c_r + dOffset);
                double d2 = blitz::dot(dx, dx);
                double dir = 1.0 / sqrt(d2);
                double tax, tay, taz;
                dFlop =
                    momLocrAddFmomr5cm(&L, &c->moment(), c->bMax(), dir, dx[0],
                                       dx[1], dx[2], &tax, &tay, &taz);
                this->dFlop += dFlop;
                this->walkStats.dFlopDoubleCPU += dFlop;
                if (options->ts.bGravStep) {
                    float adotai = a[0] * tax + a[1] * tay + a[2] * taz;
                    if (adotai > 0) {
                        adotai *= imaga;
                        dirLsum += dir * adotai * adotai;
                        normLsum += adotai * adotai;
                    }
                }
#endif
            }
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

    void calculateLocalExpansionImpl(treeStore::NodePointer k) {
        FLOCR Lf;
#ifdef USE_SIMD_FMM
        // Need to get the scaling factor correct here
        if (this->pkd->ill.count()) {
            float v = k->bMax();
            dFlop = momFlocrSetVFmomr5cm(&Lf, v, this->pkd->ill, a.data(),
                                         &dirLsum, &normLsum);
            this->dFlop += dFlop;
            this->walkStats.dFlopSingleCPU += dFlop;
            momLocrAddFlocr(&this->L, &Lf, v);
        }
#endif
    };

    void updateStackMomentTensorImpl(treeStore::NodePointer c, double xParent,
                                     double yParent, double zParent) {
        this->clStack.back().L = L;
        this->clStack.back().dirLsum = dirLsum;
        this->clStack.back().normLsum = normLsum;
        auto c_r = c->position();
        momShiftLocr(&this->clStack.back().L, c_r[0] - xParent,
                     c_r[1] - yParent, c_r[2] - zParent);
    }

    void updateMomentTensorFromStackImpl() {
        L = this->clStack.back().L;
        dirLsum = this->clStack.back().dirLsum;
        normLsum = this->clStack.back().normLsum;
    }

    void shiftToNextNodeImpl(blitz::TinyVector<double, 3> k_r, double xParent,
                             double yParent, double zParent) {
        auto dFlopMom = momShiftLocr(&this->L, k_r[0] - xParent,
                                     k_r[1] - yParent, k_r[2] - zParent);
        dFlop += dFlopMom;
        this->walkStats.dFlopDoubleCPU += dFlopMom;
    }

    bool testAddNodeCLImpl(treeStore::NodePointer k) {
        return (k->min_rung() <= this->options->ts.uRungHi &&
                k->max_rung() >= this->options->ts.uRungLo);
    }

    int interact_wrapperImpl(treeStore::NodePointer k, int iRoot, int iRoot2) {
        // Prepare the internal state of the interaction
        this->interaction.dirLsum = dirLsum;
        this->interaction.normLsum = normLsum;
        this->interaction.pLoc = &L;
        this->interaction.smx = smx;
        this->interaction.smf = smf;

        // The gravity can make use of the smx/smf structs
        return this->interaction.interact(k, iRoot, iRoot2);
    }

    void initWalkImpl() {
        this->pkd->dEnergyU = 0.0;
        this->pkd->dEnergyT = 0.0;
        this->pkd->dEnergyW = 0.0;
        this->pkd->dEnergyF[0] = this->pkd->dEnergyF[1] =
                                     this->pkd->dEnergyF[2] = 0.0;
        this->pkd->dEnergyL[0] = this->pkd->dEnergyL[1] =
                                     this->pkd->dEnergyL[2] = 0.0;

        this->pkd->fiCritTheta = 1.0f / this->options->dTheta;

        assert(this->pkd->tree.present(KDN_FIELD::oNodeMass));
        if (this->options->ts.bGravStep) {
            assert(this->pkd->tree.present(KDN_FIELD::oNodeAcceleration));
            if (this->options->ts.iTimeStepCrit == 1) {
                assert(this->pkd->tree.present(KDN_FIELD::oNodeVelocity));
                assert(this->pkd->particles.present(PKD_FIELD::oVelocity));
            }
        }

        /*
        ** Setup smooth for calculating local densities when a particle has
        *too few P-P interactions.
        */
        if (this->options->ts.bGravStep) {
            smInitializeRO(&smx, this->pkd, &smf, this->options->ts.nPartRhoLoc,
                           this->options->bPeriodic, SMX_DENSITY_F1);
            smSmoothInitialize(smx);
        }
        else
            smx = NULL;
    }

    void doneWalkImpl() {
        if (smx) {
            smSmoothFinish(smx);
            smFinish(smx, &smf);
        }
        for (auto i = 0; i <= IRUNGMAX; ++i)
            this->out->nRung[i] = this->pkd->nRung[i];
    }

    void initStatsImpl() {
        dFlop = 0.0;
        base::initStatsImpl();
    }

    void computeStatsImpl() {
        // Flop and tile count
        this->out->dFlopSingleCPU = 1e-9 * this->walkStats.dFlopSingleCPU;
        this->out->dFlopDoubleCPU = 1e-9 * this->walkStats.dFlopDoubleCPU;
        this->out->dFlopSingleGPU = 1e-9 * this->walkStats.dFlopSingleGPU;
        this->out->dFlopDoubleGPU = 1e-9 * this->walkStats.dFlopDoubleGPU;
        this->out->nTilesTotal = this->walkStats.nTilesTotal;
        this->out->nTilesCPU = this->walkStats.nTilesCPU;
        this->out->nActive = this->nActive;

        double dActive = static_cast<double>(this->nActive);
        // Set dSum for the rest of statistics
        this->out->sFlop.dSum = this->dFlop;
        if (this->nActive) {
            this->out->sCell.dSum = this->dCellSum / dActive;
            this->out->sPart.dSum = this->dPartSum / dActive;
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
        this->out->sActive.dSum = (double)this->out->nActive;

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
        auto tmp = static_cast<typename GravityWalkTraits::output *>(vout2);
        auto out = static_cast<typename GravityWalkTraits::output *>(vout);
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
        ** Combine the rung counts for the next timestep...
        */
        for (auto i = this->options->uRungLo; i <= IRUNGMAX; ++i)
            this->out->nRung[i] += tmp->nRung[i];
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
class GravityWalk
    : public GravityWalkImpl<Interaction, GravityWalk<Interaction>> {
public:
    explicit GravityWalk(PST pst)
        : GravityWalkImpl<Interaction, GravityWalk<Interaction>>(pst) {}
};
#endif
