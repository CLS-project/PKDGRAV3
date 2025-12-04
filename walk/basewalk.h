#ifndef BASE_WALK_H
#define BASE_WALK_H
#include "../TraversePST.h"
#include "gravity/cl.h"
#include "pkd.h"
#include "pkd_config.h"
#include "smooth/smooth.h"

// A little help from ChatGPT for the string concatenation at
// compile time... but probably with newer C++ standards can
// be easier
template <typename A, typename B> struct Concat {
    static constexpr auto value = [] {
        constexpr auto a = A::value;
        constexpr auto b = B::value;
        std::array<char, a.size() + b.size() + 1> buf{}; // +1 for '\0'

        for (size_t i = 0; i < a.size(); ++i)
            buf[i] = a[i];
        for (size_t i = 0; i < b.size(); ++i)
            buf[a.size() + i] = b[i];

        buf[a.size() + b.size()] = '\0';
        return buf;
    }();
};

// As with the interaction, we define traits for the
// derived tree walks to that we can access their types
// from the base class.
template <typename Traits, typename Interaction>
class BaseTreeWalk : public TraverseCombinePST {
    // We choose to derive TraverseCombinePST to allow gathering
    // of information of the underlying walk/interactions, such
    // as the statistics shown for the gravity.
public:
    using Options = typename Traits::Options;
    using output = typename Traits::output;
    using name = typename Traits::name;

    // For passing both the tree walk options and the
    // interaction options at the PST level, we pack them
    // in a single struct for convenience
    struct OptionPack {
        Options walk;
        typename Interaction::Options inter;
    };

    // The service name is building by concatenating
    // at compile time the names of the walk+interaction
    // defined in their respective traits
    static constexpr auto service_name =
        Concat<typename Traits::name, typename Interaction::name>::value;

    explicit BaseTreeWalk(PST pst)
        : TraverseCombinePST(pst, service_name.data(),
                             sizeof(OptionPack), sizeof(output)) {
        // cl = new CL();
    };

protected:
    PKD pkd;
    SMX smx;
    Interaction interaction;
    Options options;

private:
    // IA: We still use the clList defined in PKD, I do not know why
    // is done that way and if we can simply define and allocate
    // one for us here.
    // clList *cl;
    // The same goes with the stack... pkd->S

    // Any variable that must be accumulated during the walk
    double dPartSum;
    double dCellSum;
    int nActive;
    double dFlopTotal;

    void addChild(PKD pkd, int iCache, clList &cl, int iChild, int id,
                  blitz::TinyVector<float, 3> fOffset) {
        auto c = (id == pkd->Self())
                 ? pkd->tree[iChild]
                 : pkd->tree[pkd->mdl->Fetch<KDN>(iCache, iChild, id)];
        auto nc =
            (c->is_remote() || c->is_top_tree()) ? 1000000000 : c->count();
        auto cOpen = c->bMax() * pkd->fiCritTheta;
        auto c_r = c->position();
        auto cbnd = c->bound();
        auto [iLower, idLower, iUpper, idUpper] = c->get_child_cells(id);
        const auto &SPHbob = c->have_BOB() ? c->BOB() : SPHBOB();
        cl.append(iCache, id, iChild, idLower, iLower, idUpper, iUpper, nc,
                  cOpen, c->mass(), 4.0f * c->fSoft2(), c_r, fOffset, cbnd,
                  SPHbob);
    }

    void processCheckList(PKD pkd, SMX smx, SMF smf, int iRoot, int iRoot2,
                          Options options) {
        treeStore::NodePointer c(pkd->tree);
        int id, iCell, iSib, iCheckCell, iCheckLower;
        FLOCR Lf;
        float vpred[3];
        blitz::TinyVector<float, 3> a;
        blitz::TinyVector<double, 3> r, k_r, c_r, dOffset, dx;
        double xParent, yParent, zParent;
        blitz::TinyVector<float, 3> fOffset;
        float dirLsum, normLsum;
        float fMass, fSoft;
        int j, jTile, pj, nCLActive;
        static const blitz::TinyVector<float, 3> fZero3(0, 0, 0);
        int iCidPart;
        int bReferenceFound;
        double dFlop;

        pkd->dFlop = 0.0; /* Flops are accumulated here! */
        pkd->S.clear();
        pkd->S.reserve(40); // Start with a sensible stack size

        /*
        ** Clear monopole sentinel and local expansion and the timestepping
        *sums.
        */
        nCLActive = 0;
        initCLAccumulatedVariables();
        // momClearFmomr(&monoPole);
        // momClearLocr(&L);
        // dirLsum = 0;
        // normLsum = 0;
        // a = 0;

        /*
        ** We are now going to work on the local tree.
        ** Make iCell point to the root of the tree again.
        */
        auto k = pkd->tree[iCell = iRoot];
        k_r = k->position();
        while (1) {
#ifdef ILP_ILC_CAN_BE_NON_EMPTY
            // IA: are these two used?
            ILPTILE tile;
            ILCTILE ctile;
            blitz::TinyVector<double, 3> reference;
            double d2c;
            int iLower, iUpper;
            /*
            ** Find the next active particle that will be encountered in the
            *walk algorithm below
            ** in order to set a good initial center for the P-P interaction
            *list.
            */
            auto kFind = k;
            while (kFind->iLower) {
                int iCellDescend;
                kFind = pkd->tree[iCellDescend = kFind->iLower];
                if (kFind->uMinRung > uRungHi || kFind->uMaxRung < uRungLo) {
                    /*
                    ** Move onto processing the sibling.
                    */
                    kFind = pkd->tree[++iCellDescend];
                }
            }
            for (auto &p : *kfind) {
                if (!p.is_active())
                    continue;
                reference = p.position();
                goto found_it;
            }
            printf("%d: TREE ERROR\n", pkd->Self());
            assert(0); /* We didn't find an active particle */
found_it:
            reference -= interaction.getReferenceIL();
            d2c = blitz::dot(reference, reference);
            if (d2c > 1e-5) {
                // Correct all remaining PP/PC interactions to this new center.
                interaction.recenterIL(reference);
                // pkd->ilp.recenter(reference);
                // pkd->ilhydro.recenter(reference);
                // pkd->ilc.recenter(reference);
            }
            bReferenceFound = 1;
#else
            bReferenceFound = 0;
#endif
            while (1) {
                /*
                ** Process the Checklist.
                */
                beforeCLProcessing(k);
                // #ifdef USE_SIMD_FMM
                //             a = (ts->bGravStep) ? k->acceleration() : 0.0f;
                // #else
                //             float imaga;
                //             if (ts->bGravStep) {
                //                 a = k->acceleration();
                //                 imaga = 1.0 / sqrtf(blitz::dot(a,a));
                //             }
                // #endif
                // #ifdef USE_SIMD_FMM
                //             pkd->ill.clear();
                // #endif
                /*
                ** For cells which will remain on the checklist for a further
                *deeper level of
                ** level of the tree we will be using the stack cl (even if no
                *push happens later).
                */
                // IA: any reason why the stack/cl are on the PKD level? Can
                // they be placed as private variables of TreeWalk?
                clList clSib(pkd->clFreeList);
                do {
                    for (auto &tile : pkd->cl) {
                        openOutcomeSIMD(k, tile);
                    }
                    clList clNew(pkd->clFreeList);
                    for (auto &tile : pkd->cl) {
                        auto nBlocks = tile.count() / tile.width;
                        for (auto iBlock = 0; iBlock <= nBlocks; ++iBlock) {
                            int n = iBlock < nBlocks
                                    ? tile.width
                                    : tile.count() - nBlocks * tile.width;
                            auto &blk = tile[iBlock];
                            for (jTile = 0; jTile < n; ++jTile) {
                                // IA: performance consideration for having this
                                // here rather than outside the loop? It should
                                // be inlined anyway right?
                                selectOutcome(blk, blk.iOpen[jTile]);
                                // printf("%d\n", blk.iOpen[jTile]);
                                // switch (blk.iOpen[jTile]) {
                                // case 0:
                                //     /*
                                //     ** This checkcell stays on the checklist.
                                //     */
                                //     pkd->S[iStack+1].cl->append(blk,jTile);
                                //     break;
                                // case 1:
                                //     /*
                                //     ** This checkcell's particles are added
                                //     to the P-P list.
                                //     */
                                //     iCheckCell = blk.iCell[jTile];
                                //     id = blk.idCell[jTile];
                                //     if (iCheckCell < 0) {
                                //         pj = -1 - iCheckCell;
                                //         iCidPart =
                                //         blk.iCache[jTile]==CID_CELL ?
                                //         CID_PARTICLE : CID_PARTICLE2; if
                                //         (SPHoptions->doSetDensityFlags ||
                                //         SPHoptions->doSetNNflags) {
                                //             auto p = (id == pkd->Self()) ?
                                //             pkd->particles[pj]
                                //             :
                                //             pkd->particles[pkd->mdl->Acquire<PARTICLE>(iCidPart,pj,id)];
                                //             if
                                //             (SPHoptions->doSetDensityFlags)
                                //             p.set_marked(true); if
                                //             (SPHoptions->doSetNNflags)
                                //             p.set_NN_flag(true); if (id !=
                                //             pkd->Self())
                                //             pkd->mdl->Release(iCidPart,&p);
                                //         }
                                //         else {
                                //             auto p = (id == pkd->Self()) ?
                                //             pkd->particles[pj]
                                //                      :
                                //                      pkd->particles[pkd->mdl->Fetch<PARTICLE>(iCidPart,pj,id)];
                                //             r = p.position();
                                //             if (!bReferenceFound) {
                                //                 bReferenceFound=1;
                                //                 pkd->ilp.setReference(r);
                                //                 pkd->ilhydro.setReference(r);
                                //                 pkd->ilc.setReference(r);
                                //             }
                                //             if
                                //             (pkd->particles.present(PKD_FIELD::oNewSph))
                                //             {
                                //                 const auto &NewSph =
                                //                 p.newsph(); float Omega =
                                //                 pkdIsGas(pkd, &p) ?
                                //                 NewSph.Omega : 1.0f; float P
                                //                 = 0.0f; float cs = 0.0f;
                                //                 float T = 0.0f; float expImb2
                                //                 = pkdIsGas(pkd, &p) ?
                                //                 NewSph.expImb2 : 1.0f; float
                                //                 Sxx = 0.0f; float Syy = 0.0f;
                                //                 float Sxy = 0.0f;
                                //                 float Sxz = 0.0f;
                                //                 float Syz = 0.0f;
                                //                 SPHpredictOnTheFly(pkd, p,
                                //                 kick,
                                //                 SPHoptions->nPredictRung,
                                //                 vpred, &P, &cs, &T, &Sxx,
                                //                 &Syy, &Sxy, &Sxz, &Syz,
                                //                 SPHoptions); pkd->ilp.append(
                                //                     r[0] +
                                //                     blk.xOffset[jTile], r[1]
                                //                     + blk.yOffset[jTile],
                                //                     r[2] +
                                //                     blk.zOffset[jTile],
                                //                     blk.m[jTile],
                                //                     blk.fourh2[jTile],
                                //                     vpred[0], vpred[1],
                                //                     vpred[2], pkdIsGas(pkd,
                                //                     &p) ? p.ball() : 1.0f,
                                //                     Omega, pkdIsGas(pkd, &p)
                                //                     ? p.density() : 1.0f, P,
                                //                     cs, p.species(),
                                //                     p.rung(), p.imaterial(),
                                //                     T, expImb2, pkdIsGas(pkd,
                                //                     &p) ? 1.0f : 0.0f, Sxx,
                                //                     Syy, Sxy, Sxz, Syz);
                                //             } else if (p.have_mfm()) {
                                //                 auto &psph = p.sph();

                                //                float qomega = psph.omega;
                                //                float qspheom= psph.spheom;
                                //                float qH     = p.ball();
                                //                float qlast  =
                                //                psph.lastUpdateTime; float qDt
                                //                = ts->dDelta/(1<<p.rung());
                                //                float qBXX   = psph.B[XX];
                                //                float qBXY   = psph.B[XY];
                                //                float qBXZ   = psph.B[XZ];
                                //                float qBYY   = psph.B[YY];
                                //                float qBYZ   = psph.B[YZ];
                                //                float qBZZ   = psph.B[ZZ];
                                //                float qDens  = p.density();
                                //                float qvx= p.velocity()[0];
                                //                float qvy= p.velocity()[1];
                                //                float qvz= p.velocity()[2];
                                //                float qP     = psph.P;
                                //                float qS;
                                //#ifdef ENTROPY_SWITCH
                                //                qS    = psph.S;
                                //#endif
                                //                float qgradRhox =
                                //                psph.gradRho[0]; float
                                //                qgradRhoy = psph.gradRho[1];
                                //                float qgradRhoz =
                                //                psph.gradRho[2]; float qgradPx
                                //                = psph.gradP[0]; float qgradPy
                                //                = psph.gradP[1]; float qgradPz
                                //                = psph.gradP[2]; float
                                //                qgradVxx = psph.gradVx[0];
                                //                float qgradVxy =
                                //                psph.gradVx[1]; float qgradVxz
                                //                = psph.gradVx[2]; float
                                //                qgradVyx = psph.gradVy[0];
                                //                float qgradVyy =
                                //                psph.gradVy[1]; float qgradVyz
                                //                = psph.gradVy[2]; float
                                //                qgradVzx = psph.gradVz[0];
                                //                float qgradVzy =
                                //                psph.gradVz[1]; float qgradVzz
                                //                = psph.gradVz[2]; float
                                //                qlastAccx = psph.lastAcc[0];
                                //                float qlastAccy =
                                //                psph.lastAcc[1]; float
                                //                qlastAccz = psph.lastAcc[2];
                                //                //printf("Adding %ld\n",
                                //                p.ParticleID());
                                //                pkd->ilhydro.append(
                                //                    r[0] + blk.xOffset[jTile],
                                //                    r[1] + blk.yOffset[jTile],
                                //                    r[2] + blk.zOffset[jTile],
                                //                    qH, qlast, qDt, qomega,
                                //                    qspheom,
                                // qBXX, qBXY, qBXZ, qBYY, qBYZ, qBZZ, qDens,
                                // qvx, qvy, qvz, qP, qS,
                                //// qgradRhox, qgradRhoy, qgradRhoz, qgradPx,
                                /// qgradPy, qgradPz,
                                /// qgradVxx, / qgradVxy, qgradVxz, qgradVyx,
                                /// qgradVyy, qgradVyz, qgradVzx, qgradVzy, /
                                /// qgradVzz, qlastAccx, qlastAccy, qlastAccz,
                                // pkdIsGas(pkd, &p) ? 1.0f : 0.0f,
                                // pkdIsActive(pkd, &p) ? 1.0f : 0.0f,
                                // p.marked() ? 1.0f : 0.0f
                                //                     );
                                //             }
                                //             else {
                                //                 auto v = p.velocity();
                                //                 pkd->ilp.append(
                                //                     r[0] +
                                //                     blk.xOffset[jTile], r[1]
                                //                     + blk.yOffset[jTile],
                                //                     r[2] +
                                //                     blk.zOffset[jTile],
                                //                     blk.m[jTile],
                                //                     blk.fourh2[jTile], v[0],
                                //                     v[1], v[2], 0.0f, 0.0f,
                                //                     0.0f, 0.0f, 0.0f, 0,
                                //                     p.rung(), 0, 0.0f, 0.0f,
                                //                     0.0f, 0.0f, 0.0f, 0.0f,
                                //                     0.0f, 0.0f);
                                //             }
                                //         }
                                //     }
                                //     else {
                                //         assert(id >= 0);
                                //         c = (id == pkd->Self()) ?
                                //         pkd->tree[iCheckCell]
                                //             :
                                //             pkd->tree[pkd->mdl->Fetch<KDN>(blk.iCache[jTile],iCheckCell,id)];
                                //         iCidPart =
                                //         blk.iCache[jTile]==CID_CELL ?
                                //         CID_PARTICLE : CID_PARTICLE2; if
                                //         (!bReferenceFound) {
                                //             bReferenceFound=1;
                                //             auto p = (id == pkd->Self()) ?
                                //             pkd->particles[c->lower()] :
                                //             ((SPHoptions->doSetDensityFlags
                                //             || SPHoptions->doSetNNflags) ?
                                //             pkd->particles[pkd->mdl->Acquire<PARTICLE>(iCidPart,c->lower(),id)]
                                //             :
                                //             pkd->particles[pkd->mdl->Fetch<PARTICLE>(iCidPart,c->lower(),id)]);
                                //             r = p.position();
                                //             pkd->ilp.setReference(r);
                                //             pkd->ilhydro.setReference(r);
                                //             pkd->ilc.setReference(r);
                                //             if ((id != pkd->Self()) &&
                                //             (SPHoptions->doSetDensityFlags ||
                                //             SPHoptions->doSetNNflags))
                                //             pkd->mdl->Release(iCidPart,&p);
                                //         }
                                //         for (pj=c->lower(); pj<=c->upper();
                                //         ++pj) {
                                //             if (SPHoptions->doSetDensityFlags
                                //             || SPHoptions->doSetNNflags) {
                                //                 auto p = (id == pkd->Self())
                                //                 ? pkd->particles[pj] :
                                //                 pkd->particles[pkd->mdl->Acquire<PARTICLE>(iCidPart,pj,id)];
                                //                 if
                                //                 (SPHoptions->doSetDensityFlags)
                                //                 p.set_marked(true); if
                                //                 (SPHoptions->doSetNNflags)
                                //                 p.set_NN_flag(true); if (id
                                //                 != pkd->Self())
                                //                 pkd->mdl->Release(iCidPart,&p);
                                //             }
                                //             else {
                                //                 auto p = (id == pkd->Self())
                                //                 ? pkd->particles[pj]
                                //                          :
                                //                          pkd->particles[pkd->mdl->Fetch<PARTICLE>(iCidPart,pj,id)];
                                //                 fMass = p.mass();
                                //                 fSoft = p.soft();
                                //                 r = p.position();
                                //                 if (p.have_newsph()) {
                                //                     const auto &NewSph =
                                //                     p.newsph(); float Omega =
                                //                     pkdIsGas(pkd, &p) ?
                                //                     NewSph.Omega : 1.0f;
                                //                     float P = 0.0f; float cs
                                //                     = 0.0f; float T = 0.0f;
                                //                     float expImb2 =
                                //                     pkdIsGas(pkd, &p) ?
                                //                     NewSph.expImb2 : 1.0f;
                                //                     float Sxx = 0.0f; float
                                //                     Syy = 0.0f; float Sxy =
                                //                     0.0f; float Sxz = 0.0f;
                                //                     float Syz = 0.0f;
                                //                     SPHpredictOnTheFly(pkd,
                                //                     p, kick,
                                //                     SPHoptions->nPredictRung,
                                //                     vpred, &P, &cs, &T, &Sxx,
                                //                     &Syy, &Sxy, &Sxz, &Syz,
                                //                     SPHoptions);
                                //                     pkd->ilp.append(
                                //                         r[0] +
                                //                         blk.xOffset[jTile],
                                //                         r[1] +
                                //                         blk.yOffset[jTile],
                                //                         r[2] +
                                //                         blk.zOffset[jTile],
                                //                         fMass, 4*fSoft*fSoft,
                                //                         vpred[0], vpred[1],
                                //                         vpred[2],
                                //                         pkdIsGas(pkd, &p) ?
                                //                         p.ball() : 1.0f,
                                //                         Omega, pkdIsGas(pkd,
                                //                         &p) ? p.density()
                                //                         : 1.0f, P, cs,
                                //                         p.species(),
                                //                         p.rung(),
                                //                         p.imaterial(), T,
                                //                         expImb2,
                                //                         pkdIsGas(pkd, &p)
                                //                         ? 1.0f : 0.0f, Sxx,
                                //                         Syy, Sxy, Sxz, Syz);
                                //                 } else if (p.have_mfm()) {
                                //                     auto &psph = p.sph();

                                //                    float qomega = psph.omega;
                                //                    float qspheom=
                                //                    psph.spheom; float qH =
                                //                    p.ball(); float qlast  =
                                //                    psph.lastUpdateTime; float
                                //                    qDt    =
                                //                    ts->dDelta/(1<<p.rung());
                                //                    float qBXX   = psph.B[XX];
                                //                    float qBXY   = psph.B[XY];
                                //                    float qBXZ   = psph.B[XZ];
                                //                    float qBYY   = psph.B[YY];
                                //                    float qBYZ   = psph.B[YZ];
                                //                    float qBZZ   = psph.B[ZZ];
                                //                    float qDens  =
                                //                    p.density(); float qvx=
                                //                    p.velocity()[0]; float
                                //                    qvy= p.velocity()[1];
                                //                    float qvz=
                                //                    p.velocity()[2]; float qP
                                //                    = psph.P; float qS;
                                //#ifdef ENTROPY_SWITCH
                                //                    qS    = psph.S;
                                //#endif
                                //                    float qgradRhox =
                                //                    psph.gradRho[0]; float
                                //                    qgradRhoy =
                                //                    psph.gradRho[1]; float
                                //                    qgradRhoz =
                                //                    psph.gradRho[2]; float
                                //                    qgradPx = psph.gradP[0];
                                //                    float qgradPy =
                                //                    psph.gradP[1]; float
                                //                    qgradPz = psph.gradP[2];
                                //                    float qgradVxx =
                                //                    psph.gradVx[0]; float
                                //                    qgradVxy = psph.gradVx[1];
                                //                    float qgradVxz =
                                //                    psph.gradVx[2]; float
                                //                    qgradVyx = psph.gradVy[0];
                                //                    float qgradVyy =
                                //                    psph.gradVy[1]; float
                                //                    qgradVyz = psph.gradVy[2];
                                //                    float qgradVzx =
                                //                    psph.gradVz[0]; float
                                //                    qgradVzy = psph.gradVz[1];
                                //                    float qgradVzz =
                                //                    psph.gradVz[2]; float
                                //                    qlastAccx =
                                //                    psph.lastAcc[0]; float
                                //                    qlastAccy =
                                //                    psph.lastAcc[1]; float
                                //                    qlastAccz =
                                //                    psph.lastAcc[2];
                                //                    //printf("Adding %ld\n",
                                //                    p.ParticleID());
                                //                    pkd->ilhydro.append(
                                //                        r[0] +
                                //                        blk.xOffset[jTile],
                                //                        r[1] +
                                //                        blk.yOffset[jTile],
                                //                        r[2] +
                                //                        blk.zOffset[jTile],
                                //                        qH, qlast, qDt,
                                //                        qomega, qspheom,
                                // qBXX, qBXY, qBXZ, qBYY, qBYZ, qBZZ, qDens,
                                // qvx, qvy, qvz, qP, qS,
                                //// qgradRhox, qgradRhoy, qgradRhoz, qgradPx,
                                /// qgradPy, qgradPz,
                                /// qgradVxx, / qgradVxy, qgradVxz, qgradVyx,
                                /// qgradVyy, qgradVyz, qgradVzx, qgradVzy, /
                                /// qgradVzz, qlastAccx, qlastAccy, qlastAccz,
                                // pkdIsGas(pkd, &p) ? 1.0f : 0.0f,
                                // pkdIsActive(pkd, &p) ? 1.0f : 0.0f,
                                // p.marked() ? 1.0f : 0.0f
                                //                         );
                                //                 } else {
                                //                     assert(0);
                                //                     auto v = p.velocity();
                                //                     pkd->ilp.append(
                                //                         r[0] +
                                //                         blk.xOffset[jTile],
                                //                         r[1] +
                                //                         blk.yOffset[jTile],
                                //                         r[2] +
                                //                         blk.zOffset[jTile],
                                //                         fMass, 4*fSoft*fSoft,
                                //                         v[0], v[1], v[2],
                                //                         0.0f, 0.0f, 0.0f,
                                //                         0.0f, 0.0f, 0,
                                //                         p.rung(), 0, 0.0f,
                                //                         0.0f, 0.0f, 0.0f,
                                //                         0.0f, 0.0f, 0.0f,
                                //                         0.0f);
                                //                 }
                                //             }
                                //         }
                                //     }
                                //     break;
                                // case 2:
                                //     /*
                                //     ** Now I am trying to open a bucket,
                                //     which means I add each particle back on
                                //     the
                                //     ** checklist with a cell size of zero.
                                //     */
                                //     iCheckCell = blk.iCell[jTile];
                                //     assert(iCheckCell>=0);
                                //     id = blk.idCell[jTile];
                                //     fOffset[0] = blk.xOffset[jTile];
                                //     fOffset[1] = blk.yOffset[jTile];
                                //     fOffset[2] = blk.zOffset[jTile];
                                //     assert(id >= 0);
                                //     c = (id == pkd->Self()) ?
                                //     pkd->tree[iCheckCell]
                                //         :
                                //         pkd->tree[pkd->mdl->Fetch<KDN>(blk.iCache[jTile],iCheckCell,id)];
                                //     iCidPart = blk.iCache[jTile]==CID_CELL ?
                                //     CID_PARTICLE : CID_PARTICLE2;
                                //     //assert(0);
                                //     for (pj=c->lower(); pj<=c->upper(); ++pj)
                                //     {
                                //         auto p = (id == pkd->Self()) ?
                                //         pkd->particles[pj] :
                                //         ((SPHoptions->doSetDensityFlags ||
                                //         SPHoptions->doSetNNflags) ?
                                //         pkd->particles[pkd->mdl->Acquire<PARTICLE>(iCidPart,pj,id)]
                                //         :
                                //         pkd->particles[pkd->mdl->Fetch<PARTICLE>(iCidPart,pj,id)]);
                                //         r = p.position();
                                //         fMass = p.mass();
                                //         fSoft = p.soft();
                                //         if (p.have_newsph() || p.have_mfm())
                                //         {
                                //             //assert(0);
                                //             float fBallFactor =
                                //             (SPHoptions->dofBallFactor) ?
                                //             SPHoptions->fBallFactor : 1.0f;
                                //             float limitedBallSize =
                                //             std::min(SPHoptions->ballSizeLimit,fBallFactor
                                //             * p.ball());
                                //             pkd->clNew->append(blk.iCache[jTile],id,-1
                                //             -
                                //             pj,0,0,0,0,1,0.0,fMass,4.0f*fSoft*fSoft,
                                //                                r, // center
                                //                                of mass
                                //                                fOffset, //
                                //                                fOffset
                                //                                Bound(r,r), //
                                //                                zero size box
                                //                                at r
                                //                                SPHBOB(r,limitedBallSize));
                                //         }
                                //         else {
                                //             //assert(0);
                                //             pkd->clNew->append(blk.iCache[jTile],id,-1
                                //             -
                                //             pj,0,0,0,0,1,0.0,fMass,4.0f*fSoft*fSoft,
                                //                                r, // center
                                //                                of mass
                                //                                fOffset, //
                                //                                fOffset
                                //                                Bound(r,r), //
                                //                                zero size box
                                //                                at r
                                //                                SPHBOB(r,0.0));
                                //         }
                                //         if ((id != pkd->Self()) &&
                                //         (SPHoptions->doSetDensityFlags ||
                                //         SPHoptions->doSetNNflags))
                                //         pkd->mdl->Release(iCidPart,&p);
                                //     }
                                //     break;
                                // case 3:
                                //     /*
                                //     ** Open the cell.
                                //     ** We could do a prefetch here for
                                //     non-local
                                //     ** cells.
                                //     */
                                //     iCheckCell = blk.iCell[jTile];
                                //     assert(iCheckCell >= 0); iCheckLower =
                                //     blk.iLower[jTile]; assert(iCheckLower >
                                //     0);

                                //    fOffset[0] = blk.xOffset[jTile];
                                //    fOffset[1] = blk.yOffset[jTile];
                                //    fOffset[2] = blk.zOffset[jTile];

                                //    addChild(pkd,blk.iCache[jTile],pkd->clNew,blk.iLower[jTile],blk.idLower[jTile],fOffset);
                                //    addChild(pkd,blk.iCache[jTile],pkd->clNew,blk.iUpper[jTile],blk.idUpper[jTile],fOffset);

                                //    break;
                                // case 4:
                                //    /*
                                //    ** Accept multipole!
                                //    ** Interact += Moment(c);
                                //    */
                                //    if (SPHoptions->doGravity) {
                                //        iCheckCell = blk.iCell[jTile];
                                //        assert(iCheckCell>=0);
                                //        id = blk.idCell[jTile];
                                //        c = (id == pkd->Self()) ?
                                //        pkd->tree[iCheckCell]
                                //            :
                                //            pkd->tree[pkd->mdl->Fetch<KDN>(blk.iCache[jTile],iCheckCell,id)];
                                //        r[0] = blk.x[jTile] +
                                //        blk.xOffset[jTile]; r[1] =
                                //        blk.y[jTile] + blk.yOffset[jTile];
                                //        r[2] = blk.z[jTile] +
                                //        blk.zOffset[jTile]; if
                                //        (!bReferenceFound) {
                                //            bReferenceFound=1;
                                //            pkd->ilp.setReference(r);
                                //            pkd->ilhydro.setReference(r);
                                //            pkd->ilc.setReference(r);
                                //        }
                                //        pkd->ilc.append(r[0],r[1],r[2],&c->moment(),c->bMax());
                                //    }
                                //    break;
                                // case 8:
                                //    /*
                                //    ** Local expansion accepted!
                                //    */
                                //    if (SPHoptions->doGravity) {
                                //        iCheckCell = blk.iCell[jTile];
                                //        if (iCheckCell<0) {
                                //            fOffset[0] = blk.xOffset[jTile];
                                //            fOffset[1] = blk.yOffset[jTile];
                                //            fOffset[2] = blk.zOffset[jTile];
                                //            dx[0] = k_r[0] - (blk.x[jTile] +
                                //            blk.xOffset[jTile]); dx[1] =
                                //            k_r[1] - (blk.y[jTile] +
                                //            blk.yOffset[jTile]); dx[2] =
                                //            k_r[2] - (blk.z[jTile] +
                                //            blk.zOffset[jTile]);
                                //#ifdef USE_SIMD_FMM
                                //            monoPole.m = blk.m[jTile];
                                //            pkd->ill.append((float)dx[0],(float)dx[1],(float)dx[2],&monoPole,1.0);
                                //#else
                                //            double d2 = dx[0]*dx[0] +
                                //            dx[1]*dx[1] + dx[2]*dx[2]; double
                                //            dir = 1.0/sqrt(d2);
                                //            /* monoPole.m = blk.m[jTile];*/
                                //            /* *pdFlop +=
                                //            momLocrAddFmomr5cm(&L,&monoPole,0.0,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);*/
                                //            double tax,tay,taz;
                                //            dFlop =
                                //            momLocrAddMono5(&L,blk.m[jTile],dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
                                //            *pdFlop += dFlop;
                                //            pkd->dFlopDoubleCPU += dFlop;
                                //            if (ts->bGravStep) {
                                //                float adotai = a[0]*tax +
                                //                a[1]*tay + a[2]*taz; if
                                //                (adotai > 0) {
                                //                    adotai *= imaga;
                                //                    dirLsum +=
                                //                    dir*adotai*adotai;
                                //                    normLsum += adotai*adotai;
                                //                }
                                //            }
                                //#endif
                                //        }
                                //        else {
                                //            id = blk.idCell[jTile];
                                //            dOffset[0] = blk.xOffset[jTile];
                                //            dOffset[1] = blk.yOffset[jTile];
                                //            dOffset[2] = blk.zOffset[jTile];
                                //            c = (id == pkd->Self()) ?
                                //            pkd->tree[iCheckCell]
                                //                :
                                //                pkd->tree[pkd->mdl->Fetch<KDN>(blk.iCache[jTile],iCheckCell,id)];
                                //            c_r = c->position();
                                //#ifdef USE_SIMD_FMM
                                //            dx = k_r - (c_r + dOffset);
                                //            pkd->ill.append((float)dx[0],(float)dx[1],(float)dx[2],&c->moment(),c->bMax());
                                //#else
                                //            dx = k_r - (c_r + dOffset);
                                //            double d2 = blitz::dot(dx,dx);
                                //            double dir = 1.0/sqrt(d2);
                                //            double tax,tay,taz;
                                //            dFlop =
                                //            momLocrAddFmomr5cm(&L,&c->moment(),c->bMax(),dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
                                //            *pdFlop += dFlop;
                                //            pkd->dFlopDoubleCPU += dFlop;
                                //            if (ts->bGravStep) {
                                //                float adotai = a[0]*tax +
                                //                a[1]*tay + a[2]*taz; if
                                //                (adotai > 0) {
                                //                    adotai *= imaga;
                                //                    dirLsum +=
                                //                    dir*adotai*adotai;
                                //                    normLsum += adotai*adotai;
                                //                }
                                //            }
                                //#endif
                                //        }
                                //    }
                                //    break;
                                // case 10:
                                //    /*
                                //    ** This checkcell is removed from the
                                //    checklist since it has zero/negative mass.
                                //    */
                                //    break;
                                // default:
                                //    assert(0);
                                //}
                            }
                        }
                    } /* end of CL_LOOP */
                    pkd->cl = std::move(clNew);
                } while (pkd->cl.count());
                /*
                ** Now calculate the local expansion.
                */
                afterCLProcessing();
                // #ifdef USE_SIMD_FMM
                //             // Need to do something here with dirLsum and
                //             normLsum for GravStep
                //             // Need to get the scaling factor correct here
                //             if (pkd->ill.count()) {
                //                 float v = k->bMax();
                //                 dFlop =
                //                 momFlocrSetVFmomr5cm(&Lf,v,pkd->ill,a.data(),&dirLsum,&normLsum);
                //                 *pdFlop += dFlop;
                //                 pkd->dFlopSingleCPU += dFlop;
                //                 momLocrAddFlocr(&L,&Lf,v);
                //             }
                // #endif
                /*
                ** Done processing of the Checklist.
                ** Now prepare to proceed to the next deeper
                ** level of the tree.
                */
                if (k->is_group())
                    break; /* A bucket is ALWAYS a group */
                xParent = k_r[0];
                yParent = k_r[1];
                zParent = k_r[2];
                for (j = 0; j < 3; ++j)
                    fOffset[j] = 0.0f;
                iCell = k->lchild();
                k = pkd->tree[iCell];
                k_r = k->position();
                for (auto &p : *k) {
                    if (p.is_marked())
                        assert(k->is_marked());
                }
                /*
                ** Check iCell is active. We eventually want to just to a
                ** rung check here when we start using tree repair, but
                ** for now this is just as good.
                */
                if (testAddNodeCL(k)) {
                    // if ((k->min_rung()<=ts->uRungHi && k->max_rung() >=
                    // ts->uRungLo) || (SPHoptions->useDensityFlags &&
                    // k->is_marked()) || (SPHoptions->useNNflags && k->is_NN())
                    // || (SPHoptions->doMFM && k->is_marked()) ) {
                    /*
                    ** iCell is active, continue processing it.
                    ** Put the sibling onto the checklist.
                    */
                    iSib = iCell + 1;
                    c = pkd->tree[iSib];
                    /*
                    ** Test whether the sibling is active as well.
                    ** If not we don't push it onto the stack, but we
                    ** have to be careful to not pop the stack when we
                    ** hit the sibling.
                    */
                    if (testAddNodeCL(c)) {
                        // if ((c->min_rung()<=ts->uRungHi && c->max_rung() >=
                        // ts->uRungLo)
                        // || (SPHoptions->useDensityFlags && c->is_marked()) ||
                        // (SPHoptions->useNNflags && c->is_NN()) ||
                        // (SPHoptions->doMFM && c->is_marked())) {
                        /*
                        ** Sibling is active so we need to clone the checklist!
                        */
                        pkd->cl = clSib;
                        /*
                        ** Sibling is active as well.
                        ** Push Checklist for the sibling onto the stack
                        ** before proceeding deeper in the tree.
                        */
                        pkd->S.emplace_back(pkd->clFreeList);
                        pkd->S.back().cl = std::move(clSib);
                        pkd->S.back().iNodeIndex = iSib;

                        /*
                        ** At this point, the ILP is normally empty if you never
                        *do P-P except when reaching a bucket.
                        ** Softened multi-poles are also an exception.
                        */
                        // assert(pkd->ilp.count() == 0);
                        // assert(pkd->ilhydro.count() == 0);
                        // assert(pkd->ilc.count() == 0);

                        /*
                        ** Note here we already have the correct elements in
                        *S[iStack] (iStack+1 was used previously), just need to
                        *add one.
                        */
                        // IA:Are these variables even used?
                        // If so, add a hook
                        // pkd->S.back().L = L;
                        // pkd->S.back().dirLsum = dirLsum;
                        // pkd->S.back().normLsum = normLsum;
                        // c_r = c->position();
                        // momShiftLocr(&pkd->S.back().L,
                        //              c_r[0] - xParent,
                        //              c_r[1] - yParent,
                        //              c_r[2] - zParent);
                    }
                    else {
                        /*
                        ** Otherwise we can simple grab it.
                        */
                        pkd->cl = std::move(clSib);
                    }
                }
                else {
                    /*
                    ** Skip iCell, but add it to the Checklist.
                    ** No need to push anything onto the stack here.
                    ** We can simply grab it the checklist from the Stack.
                    */
                    pkd->cl = std::move(clSib);
                    /*
                    ** Move onto processing the sibling.
                    */
                    k = pkd->tree[++iCell];
                    k_r = k->position();
                }
                // IA: Again, not sure about what is going on here
                // dFlopMom = momShiftLocr(&L,k_r[0] - xParent,
                //                      k_r[1] - yParent,
                //                      k_r[2] - zParent);
                // dFlop += dFlopMom;
                // pkd->dFlopDoubleCPU += dFlopMom;
            }
            /*
            ** Now the interaction list should be complete and the
            ** Checklist should be empty! Calculate gravity on this
            ** Bucket!
            */
            auto [nActiveInt, flops] = interact_wrapper(k, iRoot, iRoot2);
            // nActiveInt = pkdGravInteract(pkd,kick,lc,ts,
            //                           k,&L,pkd->ilp,pkd->ilc,dirLsum,normLsum,bEwald,pdFlop,
            //                           smx, &smf, iRoot, iRoot2, SPHoptions,
            //                           bGPU);

            /*
            ** Update the limit for a shift of the center here based on the
            *opening radius of this
            ** cell (the one we just evaluated).
            */
            if (nActiveInt) {
                /*
                ** Here we used to set the weights of particles based on the
                *work done, but now we just assume that
                ** all particles cost the same in domain decomposition, so we
                *really don't need to set anything here.
                */
                dPartSum += nActiveInt * pkd->ilp.count();
                dCellSum += nActiveInt * pkd->ilc.count();
                nCLActive += nActiveInt;
            }
            /* Get the next cell to process from the stack */
            if (pkd->S.empty())
                goto doneCheckList;
            k = pkd->tree[iCell = pkd->S.back().iNodeIndex];
            k_r = k->position();
            /*
            ** Pop the Checklist from the top of the stack,
            ** also getting the state of the interaction list.
            */
            interaction.clearIL();
            // pkd->ilp.clear();
            // pkd->ilhydro.clear();
            // pkd->ilc.clear();
            /*
            ** Grab the checklist from the stack.
            */
            assert(pkd->cl.count() == 0);
            pkd->cl = std::move(pkd->S.back().cl);
            // IA: again these variables
            // L = pkd->S.back().L;
            // dirLsum = pkd->S.back().dirLsum;
            // normLsum = pkd->S.back().normLsum;
            pkd->S.pop_back();
        }
doneCheckList:
        doneChecklist();
        // #ifdef USE_CUDA
        //     pkd->cudaClient->flushCUDA();
        //     if (SPHoptions->doDensity) {
        //         /*
        //         ** In the density pass, a wp can requeue itself to perform
        //         the Newton iterations.
        //         ** To make sure that after the thread has reached the end of
        //         the tree walk
        //         ** the work that is still being resent by the last few wps
        //         will still be treated,
        //         ** we not only have to check (as done in mdlCompleteAllWork),
        //         ** but also flush the queues repeatedly.
        //         ** Once this while loop exits, ALL density work has been
        //         done.
        //         */
        //         // while (pkd->mdl->gpu.flushCompleted())
        //         pkd->cudaClient->flushCUDA(); while (pkd->nWpPending) {
        //             pkd->mdl->gpu.flushCompleted();
        //             while (!pkd->cudaClient->wps.empty()) {
        //                 pkdParticleWorkDone(pkd->cudaClient->wps.front());
        //                 pkd->cudaClient->wps.pop();
        //             }
        //             pkd->cudaClient->flushCUDA();
        //         }
        //     }
        // #endif
        // #ifdef USE_METAL
        //     pkd->metalClient->flushMETAL();
        // #endif
        //     pkd->mdl->CompleteAllWork();
        //     *pdFlop += pkd->dFlop; /* Accumulate work flops (notably Ewald)
        //     */
        // This is now a private variable
        // return nTotActive;
        nActive += nCLActive;
    }

    void doneChecklist() {
#ifdef USE_CUDA
        pkd->cudaClient->flushCUDA();
#endif
#ifdef USE_METAL
        pkd->metalClient->flushMETAL();
#endif
        pkd->mdl->CompleteAllWork();
        dFlopTotal += pkd->dFlop; /* Accumulate work flops (notably Ewald) */
    }

private:
    // Convenient hooks when processing the CL to add gravity-related
    // operations.
    void initCLAccumulatedVariables() {};
    void beforeCLProcessing(treeStore::NodePointer k) {};
    void afterCLProcessing() {};

    // Opening criteria
    void openOutcomeBlock(treeStore::NodePointer k, int n, clBlock &blk) {};
    void openOutcomeSIMD(treeStore::NodePointer k, clTile &tile) {
        // Do the full blocks, followed by the last (probably not full) block
        auto nBlocks = tile.count() / tile.width;
        for (auto iBlock = 0; iBlock < nBlocks; ++iBlock) {
            openOutcomeBlock(k, tile.width, tile[iBlock]);
        }
        openOutcomeBlock(k, tile.count() - nBlocks * tile.width, tile[nBlocks]);
    }

    // Based on openOutcomeBlock do something.
    // We split it so allow for child classes
    // to reuse the same criteria easily
    void selectOutcome(clBlock blk, int outcome) {};

    // Whether a node is eligible for adding
    // it to the interaction list
    bool testAddNodeCL(treeStore::NodePointer k) { return true; };

    // pkdGravWalk
    void walk() {
        int id;
        float fOffset[3];
        int ix, iy, iz;
        // Not sure how to handle this.. this seems to be used
        // only for the sph density?
        SMX smx;
        SMF smf;
        int iTop1, iTop2;

        initWalk();

        iTop1 = pkd->iTopTree[options.iLocalRoot1];
        iTop2 = pkd->iTopTree[options.iLocalRoot2];
        id = pkd->Self();

        /*
        ** Walk tree 1 against trees 1 (and optionally 2) if there are active
        *particles
        */
        auto k = pkd->tree[iTop1];
        if (k->count() && testAddNodeCL(k)) {
            // if (k->count() && (k->is_active(ts->uRungLo,ts->uRungHi) ||
            // (SPHoptions->useDensityFlags && k->is_marked())  ||
            // (SPHoptions->useNNflags && k->is_NN()) || ( SPHoptions->doMFM &&
            // k->is_marked()) )) {
            /*
            ** Initially we set our cell pointer to
            ** point to the top tree.
            */
            interaction.clearIL();
            // pkd->ilp.clear();
            // pkd->ilhydro.clear();
            // pkd->ilc.clear();
            pkd->cl.clear();

            /*
            ** Add all replicas of the entire box to the Checklist.
            ** We add at least one box (0,0,0). The root cell is alway on
            *processor 0.
            */
            for (ix = -options.nReps; ix <= options.nReps; ++ix) {
                fOffset[0] = ix * pkd->fPeriod[0];
                for (iy = -options.nReps; iy <= options.nReps; ++iy) {
                    fOffset[1] = iy * pkd->fPeriod[1];
                    for (iz = -options.nReps; iz <= options.nReps; ++iz) {
                        fOffset[2] = iz * pkd->fPeriod[2];
                        addChild(pkd, CID_CELL, pkd->cl, iTop1, id, fOffset);
#ifndef SINGLE_CACHES
                        if (iTop2 > 0)
                            addChild(pkd, CID_CELL2, pkd->cl, iTop2, id,
                                     fOffset);
#else
                        if (iTop2 > 0)
                            addChild(pkd, CID_CELL, pkd->cl, iTop2, id,
                                     fOffset);
#endif
                    }
                }
            }
            processCheckList(pkd, smx, smf, iTop1, iTop2, options);
        }
        doneWalk();
    }

    // initWalkGrav
    void initWalk() {};

    // doneWalkGrav
    void doneWalk() {};

    // Allow derived classes to start a cache with a different helper to allow
    // for packing, for example
    void setupCache() {
        pkd->mdl->CacheInitialize(CID_PARTICLE,
                                  std::make_shared<particleCACHEhelper>(pkd));
    }

    // Initialize flop/counts statistics for latter analisis
    void initStats() {};

    // Process the gatheres statistics during the walk (which should be
    // available as private variables)
    void computeStats() {};

    // We abstract this a bit to allow passing different interaction lists
    // depending
    std::tuple<int, double> interact_wrapper(treeStore::NodePointer k,
            int iRoot, int iRoot2) {
        return interaction.interact(k, iRoot, iRoot2);
    }

    static_assert(std::is_void<Options>() ||
                  std::is_standard_layout<Options>());
    static_assert(std::is_void<output>() || std::is_standard_layout<output>());

public:
    // pkdGravAll
    int Service(PST pst, void *vin, int nIn, void *vout, int nOut) {
        this->pkd = pst->plcl->pkd;
        auto *optpack = static_cast<OptionPack *>(vin);
        typename Interaction::Options opt = optpack->inter;
        this->options = optpack->walk;
        output *out = static_cast<output *>(vout);
        assert(sizeof(this->options) == nIn);
        interaction.setup(this->pkd, opt);

        // options.copyToDevice(); // or something similar
        setupCache();

        initStats();

        walk();

        computeStats();

        pkd->mdl->FinishCache(CID_PARTICLE);
        // This kind of thing should go inside doneWalk
        // for (auto i = 0; i <= IRUNGMAX; ++i)
        //   nRung[i] = pkd->nRung[i];
        return sizeof(output);
    }

    // pstGravity combination of statistics
    virtual int Combine(void *vout, void *vout2) { return sizeof(output); };
};
#endif
