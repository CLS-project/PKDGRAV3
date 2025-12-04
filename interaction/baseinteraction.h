#ifndef BASE_INTERACTION_H
#define BASE_INTERACTION_H
#include "pkd.h"
#include "pkd_config.h"

// Question, how to make this scheme to work on smooths
// like the gradients that does:
// 1. loop over ngb to gather some number
// 2. do something with that number
// 3. loop again using that value to gather other values
// 4. do the final computation for the sink particle based
//      on those gathererd values
//
// I guess that you can do multiple queues

// Another important question:
// The p-p force is computed twice in the current
// implementation (afaik), any idea about how to
// improve this?

// And related to this:
// It seems that there is always a gather approach,
// rather than scatter. This is convenient because
// there is no need for combiner caches and locking.
// But requires to mark those particles that, although
// not active, will interact with an active and compute
// their interactions. Still better to use the gather?
//
// In my prototype implementation of MFM on GPUs I the
// marking using a normal ReSmooth, as doing it within
// the gravity loop (as with NewSPH) was way slower in
// my tests.

// I'm using traits classes for the interactions.
// These contains all the types information used for
// the class, and is passed here as a template parameter
// so that the base class can use the derived class'
// types. All this at compile-time, so no run-time
// overhead due to virtualization or anything.
template <typename Traits> class BaseInteraction {

public:
    // Access the derived classes types via trait class
    using IL = typename Traits::IL;
    using Options = typename Traits::Options;
    using PINFOIN = typename Traits::PINFOIN;
    using PINFOOUT = typename Traits::PINFOOUT;
    using name = typename Traits::name;
    using ResultInteraction = typename Traits::template ResultInteraction<fvec>;

    // Default constructor without argument is needed when
    // building the service. We do a delayed initialization
    // later in setup(), which is called at Service runtime,
    // not construction.
    BaseInteraction() { il = new IL(); };

    void setup(PKD pkd, Options options) {
        this->pkd = pkd;
        this->options = options;
    }

    // Generic interfaces to the interaction list(s).
    // If the interaction has more than one, these
    // functions should be overloaded accordingly.
    void clearIL() { il->clear(); };
    void recenterIL(blitz::TinyVector<double, 3> ref) { il->recenter(ref); };
    void setReferenceIL(blitz::TinyVector<double, 3> ref) {
        il->setReference(ref);
    };
    blitz::TinyVector<double, 3> getReferenceIL() {
        return il->getReference();
    };

    // Test that the IL(s) has been populated correctly, only
    // for debugging runs
    void extensiveILTest() {};

protected:
    PKD pkd;
    Options options;

    // By default we only have one interaction list,
    // but derived classes can manually add more
    // if required
    IL *il;

    struct workParticle {
        // Use a constructor for convenience
        workParticle(PKD pkd, treeStore::NodePointer k) {
            nRefs = 1; /* I am using it currently */
            ctx = pkd;
            dFlop = 0.0;
            dFlopSingleCPU = dFlopSingleGPU = 0.0;
            dFlopDoubleCPU = dFlopDoubleGPU = 0.0;
            nP = 0;

            auto knP = k->count();
            pPart = new PARTICLE*[knP];
            iPart = new uint32_t[knP];
            pInfoIn = new PINFOIN[knP];
            pInfoOut = new PINFOOUT[knP];
            c[0] = il->getReference(0);
            // assert(wp->c[0] == ilc.getReference(0));
            c[1] = il->getReference(1);
            // assert(wp->c[1] == ilc.getReference(1));
            c[2] = il->getReference(2);
            // assert(wp->c[2] == ilc.getReference(2));
        };
        IL *il;

        PARTICLE **pPart;
        uint32_t *iPart;
        PINFOIN *pInfoIn;
        PINFOOUT *pInfoOut;
        double dFlop;
        blitz::TinyVector<double, 3> c;
        int nP;
        int nRefs;
        void *ctx;

        double dFlopSingleCPU;
        double dFlopSingleGPU;
        double dFlopDoubleCPU;
        double dFlopDoubleGPU;
    };

    // This hook allows or explicit cloning of the IL
    // if needed (it seems to be used when bDensity & bGPU)
    void cloneIL(workParticle *wp) { wp->il = il; }

    // This would wrap all the queues or any processing
    // required by a derived class.
    // By default just queue it.
    void processWorkParticle(workParticle *wp) { queueInteraction(wp, il); };

    // Acts on own workParticle, overloaded by the
    // different choices in pkdParticleWorkDone
    void particleWorkDone(workParticle *wp) {};

    // queuePP, queueMFM,...
    void queueInteraction(workParticle *wp, IL *il) {};

    // Condition for not processing a particle
    bool particleSkipCondition(particleStore::ParticlePointer p) {
        return !p->is_rung_range(options.uRungLo, options.uRungHi);
    };

    void aggregateResult(PINFOOUT &out, ResultInteraction &result) = delete;
    void evalInteraction(const PINFOIN &Part, typename IL::tile_type &tile,
                         PINFOOUT &out) = delete;

public:
    // pkdGravInteract
    // Returns the number of interactions and the flops
    std::tuple<int, double> interact(treeStore::NodePointer pkdn, int iRoot1,
                                     int iRoot2) {
        float fBall;
        int i, nSoft, nActive;
        int nP;

        auto kdn_r = pkdn->position();

        /*
        ** Now process the two interaction lists for each active particle.
        */
        nActive = 0;
        nSoft = 0;

        /* Collect the bucket particle information */
        auto wp = new workParticle(pkd, pkdn);
        assert(wp != NULL);

        ++pkd->nWpPending;
        // wp->nRefs = 1; /* I am using it currently */
        // wp->ctx = pkd;
        // wp->bGPU = bGPU;
        // wp->dFlop = 0.0;
        // wp->dFlopSingleCPU = wp->dFlopSingleGPU = 0.0;
        // wp->dFlopDoubleCPU = wp->dFlopDoubleGPU = 0.0;
        // wp->nP = 0;

        // nP = pkdn->count();
        // wp->pPart = new PARTICLE *[nP];
        // wp->iPart = new uint32_t[nP];
        // wp->pInfoIn = new PINFOIN[nP];
        // wp->pInfoOut = new PINFOOUT[nP];
        // wp->c[0] = il.getReference(0);// assert(wp->c[0] ==
        // ilc.getReference(0)); wp->c[1] = il.getReference(1);//
        // assert(wp->c[1] == ilc.getReference(1)); wp->c[2] =
        // il.getReference(2);// assert(wp->c[2] == ilc.getReference(2));

        // wp->options = options;

        for (auto &p : *pkdn) {
            if (particleSkipCondition(&p))
                continue;
            // if (SPHoptions->doMFM) {
            //     if (!p.is_marked()) continue;
            // }else{
            //     if (!p.is_rung_range(ts->uRungLo,ts->uRungHi)
            //         && !(SPHoptions->useDensityFlags && p.marked())
            //         && !(SPHoptions->useNNflags && p.NN_flag()))
            //         continue;
            // }

            auto r = p.position();

            nP = wp->nP++;
            wp->pPart[nP] = &p;
            wp->iPart[nP] = &p - pkd->particles.begin();
            wp->pInfoIn[nP].init(&p);
            wp->pInfoOut[nP].init(&p);

            //                 wp->pInfoIn[nP].r = pkd->ilc.get_dr(r);

            //                 if (p.have_newsph()) {
            //                     const auto &NewSph = p.newsph();
            //                     wp->pInfoIn[nP].fBall = pkdIsGas(pkd, &p) ?
            //                     p.ball() : 1.0f; wp->pInfoIn[nP].isTooLarge =
            //                     0; wp->pInfoIn[nP].Omega = pkdIsGas(pkd, &p)
            //                     ? NewSph.Omega : 1.0f; wp->pInfoIn[nP].iMat =
            //                     p.imaterial(); SPHpredictOnTheFly(pkd, p,
            //                     kick, wp->SPHoptions->nPredictRung,
            //                     wp->pInfoIn[nP].v, &wp->pInfoIn[nP].P,
            //                     &wp->pInfoIn[nP].cs, NULL,
            //                     &wp->pInfoIn[nP].Sxx, &wp->pInfoIn[nP].Syy,
            //                     &wp->pInfoIn[nP].Sxy, &wp->pInfoIn[nP].Sxz,
            //                     &wp->pInfoIn[nP].Syz, SPHoptions);
            //                     wp->pInfoIn[nP].rho = pkdIsGas(pkd, &p) ?
            //                     p.density() : 1.0f; wp->pInfoIn[nP].isGas =
            //                     pkdIsGas(pkd, &p) ? 1.0f : 0.0f;

            //                     wp->pInfoOut[nP].rho = 0.0f;
            //                     wp->pInfoOut[nP].drhodfball = 0.0f;
            //                     wp->pInfoOut[nP].nden = 0.0f;
            //                     wp->pInfoOut[nP].dndendfball = 0.0f;
            //                     wp->pInfoOut[nP].fBall = 0.0f;
            //                     wp->pInfoOut[nP].nSmooth = 0.0f;
            //                     wp->pInfoOut[nP].imbalanceX = 0.0f;
            //                     wp->pInfoOut[nP].imbalanceY = 0.0f;
            //                     wp->pInfoOut[nP].imbalanceZ = 0.0f;
            //                     wp->pInfoOut[nP].uDot = 0.0f;
            //                     wp->pInfoOut[nP].divv = 0.0f;
            //                     wp->pInfoOut[nP].dtEst = 1e14f;
            //                     wp->pInfoOut[nP].maxRung = 0.0f;
            //                     wp->pInfoOut[nP].corrT = 0.0f;
            //                     wp->pInfoOut[nP].corrP = 0.0f;
            //                     wp->pInfoOut[nP].corr = 0.0f;
            //                     if (SPHoptions->doShearStrengthModel) {
            //                         wp->pInfoOut[nP].dvxdx = 0.0f;
            //                         wp->pInfoOut[nP].dvxdy = 0.0f;
            //                         wp->pInfoOut[nP].dvxdz = 0.0f;
            //                         wp->pInfoOut[nP].dvydx = 0.0f;
            //                         wp->pInfoOut[nP].dvydy = 0.0f;
            //                         wp->pInfoOut[nP].dvydz = 0.0f;
            //                         wp->pInfoOut[nP].dvzdx = 0.0f;
            //                         wp->pInfoOut[nP].dvzdy = 0.0f;
            //                         wp->pInfoOut[nP].dvzdz = 0.0f;
            //                         wp->pInfoOut[nP].Cinvxx = 0.0f;
            //                         wp->pInfoOut[nP].Cinvxy = 0.0f;
            //                         wp->pInfoOut[nP].Cinvxz = 0.0f;
            //                         wp->pInfoOut[nP].Cinvyx = 0.0f;
            //                         wp->pInfoOut[nP].Cinvyy = 0.0f;
            //                         wp->pInfoOut[nP].Cinvyz = 0.0f;
            //                         wp->pInfoOut[nP].Cinvzx = 0.0f;
            //                         wp->pInfoOut[nP].Cinvzy = 0.0f;
            //                         wp->pInfoOut[nP].Cinvzz = 0.0f;
            //                     }
            //                 }
            //                 if (p.have_mfm()) {
            //                     if (pkdIsGas(pkd, &p)) {
            //                         const auto &psph = p.sph();
            //                         wp->pInfoIn[nP].fBall = pkdIsGas(pkd, &p)
            //                         ? p.ball() : 1.0f; wp->pInfoIn[nP].isGas
            //                         = pkdIsGas(pkd, &p) ? 1.0f : 0.0f;
            //                         wp->pInfoIn[nP].isActive =
            //                         pkdIsActive(pkd, &p) ? 1.0f : 0.0f;
            //                         wp->pInfoIn[nP].isMarked = p.marked()
            //                         ? 1.0f : 0.0f; wp->pInfoIn[nP].Omega =
            //                         psph.omega; wp->pInfoIn[nP].rho =
            //                         pkdIsGas(pkd, &p) ? p.density() : 1.0f;
            //                         wp->pInfoIn[nP].spheom= psph.spheom;
            //                         wp->pInfoIn[nP].last  =
            //                         psph.lastUpdateTime; wp->pInfoIn[nP].Dt =
            //                         wp->ts->dDelta/(1<<p.rung());
            //                         wp->pInfoIn[nP].BXX   = psph.B[XX];
            //                         wp->pInfoIn[nP].BXY   = psph.B[XY];
            //                         wp->pInfoIn[nP].BXZ   = psph.B[XZ];
            //                         wp->pInfoIn[nP].BYY   = psph.B[YY];
            //                         wp->pInfoIn[nP].BYZ   = psph.B[YZ];
            //                         wp->pInfoIn[nP].BZZ   = psph.B[ZZ];
            //                         wp->pInfoIn[nP].v[0] = p.velocity()[0];
            //                         wp->pInfoIn[nP].v[1] = p.velocity()[1];
            //                         wp->pInfoIn[nP].v[2] = p.velocity()[2];
            //                         wp->pInfoIn[nP].P  = psph.P;
            //                         wp->pInfoIn[nP].S = 0.0f; //psph.S;
            //                         wp->pInfoIn[nP].gradRhox =
            //                         psph.gradRho[0]; wp->pInfoIn[nP].gradRhoy
            //                         = psph.gradRho[1];
            //                         wp->pInfoIn[nP].gradRhoz =
            //                         psph.gradRho[2]; wp->pInfoIn[nP].gradPx =
            //                         psph.gradP[0]; wp->pInfoIn[nP].gradPy =
            //                         psph.gradP[1]; wp->pInfoIn[nP].gradPz =
            //                         psph.gradP[2]; wp->pInfoIn[nP].gradVxx =
            //                         psph.gradVx[0]; wp->pInfoIn[nP].gradVxy =
            //                         psph.gradVx[1]; wp->pInfoIn[nP].gradVxz =
            //                         psph.gradVx[2]; wp->pInfoIn[nP].gradVyx =
            //                         psph.gradVy[0]; wp->pInfoIn[nP].gradVyy =
            //                         psph.gradVy[1]; wp->pInfoIn[nP].gradVyz =
            //                         psph.gradVy[2]; wp->pInfoIn[nP].gradVzx =
            //                         psph.gradVz[0]; wp->pInfoIn[nP].gradVzy =
            //                         psph.gradVz[1]; wp->pInfoIn[nP].gradVzz =
            //                         psph.gradVz[2]; wp->pInfoIn[nP].lastAccx
            //                         = psph.lastAcc[0];
            //                         wp->pInfoIn[nP].lastAccy =
            //                         psph.lastAcc[1]; wp->pInfoIn[nP].lastAccz
            //                         = psph.lastAcc[2]; wp->pInfoOut[nP].Frho
            //                         = 0.0f; wp->pInfoOut[nP].Fmom[0] = 0.0f;
            //                         wp->pInfoOut[nP].Fmom[1] = 0.0f;
            //                         wp->pInfoOut[nP].Fmom[2] = 0.0f;
            //                         wp->pInfoOut[nP].Fene = 0.0f;
            //                         wp->pInfoOut[nP].FS = 0.0f;
            //                     }
            //                 }

            // wp->pInfoOut[nP].a[0] = 0.0f;
            // wp->pInfoOut[nP].a[1] = 0.0f;
            // wp->pInfoOut[nP].a[2] = 0.0f;
            // wp->pInfoOut[nP].fPot = 0.0f;
            // wp->pInfoOut[nP].dirsum = dirLsum;
            // wp->pInfoOut[nP].normsum = normLsum;
            // wp->pInfoOut[nP].rhopmax = 0.0f;

            /*
            ** Calculate local density and kernel smoothing
            ** length for dynamical time-stepping
            */
            // IA: Maybe we need some type of hook here
            // for this case?
            // if (SPHoptions->doGravity) {
            //     if (ts->bGravStep) {
            //         /*
            //         ** Calculate local density using smooth; this is fast
            //         because the particles are
            //         ** likely to be cached already because they will be on
            //         the P-P list.
            //         */
            //         smf->pfDensity = &wp->pInfoIn[nP].fDensity;
            //         fBall = smSmoothSingle(smx,smf,p,iRoot1,iRoot2);
            //         wp->pInfoIn[nP].fSmooth2 = fBall * fBall;
            //     }
            //     else {
            //         /*
            //         ** We are not using GravStep!
            //         */
            //         wp->pInfoIn[nP].fSmooth2 = 0.0;
            //     }
            // }
        }

        /*
        ** If CUDA is used, and we are doing density on the GPU,
        ** we need to clone the ilp,
        ** otherwise we use the reference to the global ilp.
        */
        cloneIL(wp);
        // #ifdef USE_CUDA
        //             if (SPHoptions->doDensity && bGPU) {
        //                 wp->ilp = new ilpList;
        //                 wp->ilp->clone(pkd->ilp);
        //                 wp->ilhydro->clone(pkd->ilhydro);
        //                 wp->ilp->setReference(pkd->ilp.getReference());
        //                 wp->ilhydro->setReference(pkd->ilhydro.getReference());
        //             }
        //             else {
        // #endif
        //                 wp->ilp = &ilp;
        // #ifdef USE_CUDA
        //             }
        // #endif

        nActive += wp->nP;

        if (options.doExtensiveTest) {
            extensiveILTest();
        }
        // if (SPHoptions->doExtensiveILPTest && !(SPHoptions->doSetDensityFlags
        // || SPHoptions->doSetNNflags)) {
        //     extensiveILPTest(pkd, wp, pkd->ilp);
        // }
        // extensiveILHTest(pkd, wp, pkd->ilhydro);

        processWorkParticle(wp);
        /*
        ** Evaluate the local expansion.
        */
        // if (pLoc && SPHoptions->doGravity) {
        //     for ( i=0; i<wp->nP; i++ ) {
        //         momFloat ax,ay,az, dPot;
        //         blitz::TinyVector<double,3> r = wp->c + wp->pInfoIn[i].r;

        //         /*
        //         ** Evaluate local expansion.
        //         */
        //         blitz::TinyVector<double,3> dr = r - kdn_r;
        //         dPot = 0;
        //         ax = 0;
        //         ay = 0;
        //         az = 0;

        //         momEvalLocr(pLoc,dr[0],dr[1],dr[2],&dPot,&ax,&ay,&az);

        //         wp->pInfoOut[i].fPot += dPot;
        //         wp->pInfoOut[i].a[0] += ax;
        //         wp->pInfoOut[i].a[1] += ay;
        //         wp->pInfoOut[i].a[2] += az;
        //     }
        // }

        // if (SPHoptions->doGravity) {
        //     /*
        //     ** Evaluate the P-C interactions
        //     */
        //     queuePC( pkd,  wp, pkd->ilc, ts->bGravStep, bGPU);

        //     /*
        //     ** Evaluate the P-P interactions
        //     */
        //     queuePP( pkd, wp, pkd->ilp, ts->bGravStep, bGPU);
        // }

        // if (SPHoptions->doMFM) {
        //     queueMFM( pkd, wp, pkd->ilhydro, bGPU);
        // }

        // if (SPHoptions->doDensity) {
        //     /*
        //     ** Evaluate the Density on the P-P interactions
        //     */
        //     queueDensity( pkd, wp, pkd->ilp, bGPU);
        // }

        // if (SPHoptions->doDensityCorrection) {
        //     /*
        //     ** Evaluate the weighted averages of P and T
        //     */
        //     queueDensityCorrection( pkd, wp, pkd->ilp, bGPU);
        // }

        // if (SPHoptions->doSPHForces) {
        //     /*
        //     ** Evaluate the SPH forces on the P-P interactions
        //     */
        //     queueSPHForces( pkd, wp, pkd->ilp, bGPU);
        //     if (SPHoptions->doCentrifugal) {
        //         addCentrifugalAcceleration(pkd, wp);
        //     }
        //     if (wp->SPHoptions->doExternalAcceleration) {
        //         for (int i=0; i< wp->nP; i++) {
        //             wp->pInfoOut[i].a[0] +=
        //             wp->SPHoptions->externalAccelerationX;
        //             wp->pInfoOut[i].a[1]
        //             += wp->SPHoptions->externalAccelerationY;
        //             wp->pInfoOut[i].a[2] +=
        //             wp->SPHoptions->externalAccelerationZ;
        //         }
        //     }
        // }

        // if (SPHoptions->doGravity) {
        //     /*
        //     ** Calculate the Ewald correction for this particle, if it is
        //     required.
        //     */
        //     if (bEwald) {
        //         queueEwald( pkd,  wp, bGPU );
        //     }
        // }

        // #ifdef TIMESTEP_CRITICAL
        // if (SPHoptions->doGravity) {
        //     for ( i=0; i<wp->nP; i++ ) {
        //         double *c = wp->c;
        //         float *in = wp->pInfoIn[i].r;
        //         r[0] = c[0] + in[0];
        //         r[1] = c[1] + in[1];
        //         r[2] = c[2] + in[2];
        //         float fMass = p.mass();
        //         float fSoft = p.soft();
        //         const auto &v = p.velocity();
        //         /*
        //         ** Set value for time-step, note that we have the current
        //         ewald acceleration
        //         ** in this now as well!
        //         */
        //         if (ts->bGravStep && wp->ts->iTimeStepCrit == 1) {
        //             double vx,vy,vz;
        //             float d2,dir,dir2;
        //             float fx,fy,fz;
        //             float rhopmax,rhopmaxlocal;
        //             float summ;
        //             float fourh2;
        //             int n;
        //             ILPTILE tile;

        //             /*
        //             ** GravStep if iTimeStepCrit =
        //             ** 0: Mean field regime for dynamical time
        //             (normal/standard setting)
        //             ** 1: Gravitational scattering regime for dynamical time
        //             with eccentricity correction
        //             */
        //             rhopmax = 0.0;
        //             ILP_LOOP(ilp,tile) {
        //                 int blk,prt;
        //                 for ( blk=0; blk<=tile->lstTile.nBlocks; ++blk ) {
        //                     n = (blk==tile->lstTile.nBlocks ?
        //                     tile->lstTile.nInLast : ilp->lst.nPerBlock); for
        //                     (prt=0; prt<n; ++prt) {
        //                         if (p.order() < wp->ts->nPartColl ||
        //                         tile->xtr[blk].iOrder.i[prt] <
        //                         wp->ts->nPartColl)
        //                         {
        //                             fx = r[0] - ilp->cx +
        //                             tile->blk[blk].dx.f[prt]; fy = r[1] -
        //                             ilp->cy
        //                             + tile->blk[blk].dy.f[prt]; fz = r[2] -
        //                             ilp->cz + tile->blk[blk].dz.f[prt]; d2 =
        //                             fx*fx + fy*fy + fz*fz; if (p.order() ==
        //                             tile->xtr[blk].iOrder.i[prt]) continue;
        //                             fourh2 =
        //                             softmassweight(fMass,4*fSoft*fSoft,
        //                                                     tile->blk[blk].m.f[prt],tile->blk[blk].fourh2.f[prt]);
        //                             if (d2 > fourh2) {
        //                                 SQRT1(d2,dir);
        //                                 dir2 = dir*dir*dir;
        //                             }
        //                             else {
        //                                 /*
        //                                 ** This uses the Dehnen K1 kernel
        //                                 function now, it's fast!
        //                                 */
        //                                 SQRT1(fourh2,dir);
        //                                 dir2 = dir*dir;
        //                                 d2 *= dir2;
        //                                 dir2 *= dir;
        //                                 d2 = 1 - d2;
        //                                 dir *= 1.0f + d2*(0.5f +
        //                                 d2*(3.0f/8.0f + d2*(45.0f/32.0f)));
        //                                 dir2 *= 1.0f + d2*(1.5f +
        //                                 d2*(135.0f/16.0f));
        //                             }
        //                             summ = fMass+tile->blk[blk].m.f[prt];
        //                             rhopmaxlocal = summ*dir2;
        //                             vx = v[0] - tile->xtr[blk].vx.d[prt];
        //                             vy = v[1] - tile->xtr[blk].vy.d[prt];
        //                             vz = v[2] - tile->xtr[blk].vz.d[prt];
        //                             rhopmaxlocal =
        //                             pkdRho1(rhopmaxlocal,summ,dir,
        //                                                    fx,fy,fz,vx,vy,vz,wp->ts->dEccFacMax);
        //                             rhopmax = (rhopmaxlocal >
        //                             rhopmax)?rhopmaxlocal:rhopmax;
        //                         }
        //                     }
        //                 }
        //             }
        //             wp->pInfoOut[i].rhopmax = rhopmax;
        //         }
        //     } /* end of i-loop cells & particles */
        // }
        // #endif

        particleWorkDone(wp);

        //*pdFlop += nActive*(pkd->ilp.count()*COST_FLOP_PP +
        // pkd->ilc.count()*COST_FLOP_PC) + nSoft*15; return (nActive);
        return std::tuple<int, double>(nActive, wp->dFlop);
    };
};

#endif
