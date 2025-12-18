#ifndef BASE_INTERACTION_H
#define BASE_INTERACTION_H
#include "cuda/cudautil.h"
#include "gpu/pppcdata.h"
#include "pkd.h"
#include "pkd_config.h"
#include "basetype.h"
#include "walk/stats.h"

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

template <typename Traits, class Derived> class BaseInteractionInterface {
public:
    using PINFOIN = typename Traits::PINFOIN;
    using PINFOOUT = typename Traits::PINFOOUT;
    using IL = typename Traits::IL;
    using ILtile = typename Traits::IL::tile_type;
    using interWorkParticle = typename Traits::WORKPARTICLE;
    inline void setup(PKD pkd, typename Traits::Options *options,
                      WalkStats *walkStats) {
        derived().setupImpl(pkd, options, walkStats);
    };
    inline void clearIL() { derived().clearILImpl(); };
    inline int countIL() { return derived().countILImpl(); };
    inline void setReferenceIL(blitz::TinyVector<double, 3> ref) {
        derived().setReferenceILImpl(ref);
    };
    inline blitz::TinyVector<double, 3> getReferenceIL() {
        return derived().getReferenceILImpl();
    };
    void addParticleIL(double x, double y, double z, float fMass, float fourh2, particleStore::ParticlePointer p) {
        derived().addParticleILImpl(x,y,z,fMass,fourh2,p);
    }
    inline void extensiveILTest(interWorkParticle *wp) { derived().extensiveILTestImpl(wp); };
    inline void initializeInfoIn(interWorkParticle *wp, int nP,
                                 particleStore::ParticlePointer p) {
        derived().initializeInfoInImpl(wp, nP, p);
    };
    inline void initializeInfoOut(interWorkParticle *wp, int nP,
                                  particleStore::ParticlePointer p) {
        derived().initializeInfoOutImpl(wp, nP, p);
    };
    inline void cloneIL(interWorkParticle *wp) { derived().cloneILImpl(wp); };
    inline void setupCache() { setupCacheImpl(); };
    inline void queueInteraction(interWorkParticle *wp) {
        derived().queueInteractionImpl(wp);
    };
    inline int queueInteractionGPU(interWorkParticle *wp, ILtile tile) {
        return derived().queueInteractionGPUImpl(wp, tile);
    };
    inline void processWorkParticle(interWorkParticle *wp) {
        derived().processWorkParticleImpl(wp);
    };
    inline void particleWorkDone(interWorkParticle *wp) {
        derived().particleWorkDoneImpl(wp);
    };
    inline bool particleSkipCondition(particleStore::ParticlePointer p) {
        return derived().particleSkipConditionImpl(p);
    };
    inline void evalInteraction(const PINFOIN &Part,
                                ILtile &tile,
                                PINFOOUT &out) {
        derived().evalInteractionImpl(Part, tile, out);
    };
    inline int interact(treeStore::NodePointer pkdn,
                        int iRoot1, int iRoot2) {
        return derived().interactImpl(pkdn, iRoot1, iRoot2);
    };
    static void trampoline(void *ctx, workParticleVariants wp) {
        static_cast<Derived *>(ctx)->particleWorkDone(*std::get_if<interWorkParticle *>(&wp));
    }
    void registerCudaCallback(CudaClient &cuda) {
        cuda.ctx = this;
        cuda.particleWorkDone = &BaseInteractionInterface::trampoline;
    }

protected:
    inline void setupImpl(PKD pkd, typename Traits::Options *options,
                          WalkStats *walkStats) {
        derived().setupImpl(pkd, options, walkStats);
    };
    inline void clearILImpl() { derived().clearILImpl(); };
    inline int countILImpl() { return derived().countILImpl(); };
    inline void setReferenceILImpl(blitz::TinyVector<double, 3> ref) {
        derived().setReferenceILImpl(ref);
    };
    inline blitz::TinyVector<double, 3> getReferenceILImpl() {
        return derived().getReferenceILImpl();
    };
    void addParticleILImpl(double x, double y, double z, float fMass, float fourh2, particleStore::ParticlePointer p) {
        assert(0);
    }
    inline void extensiveILTestImpl(interWorkParticle *wp) { derived().extensiveILTestImpl(wp); };
    inline void initializeInfoInImpl(interWorkParticle *wp, int nP,
                                     particleStore::ParticlePointer p) {
        derived().initializeInfoInImpl(wp, nP, p);
    };
    inline void initializeInfoOutImpl(interWorkParticle *wp, int nP,
                                      particleStore::ParticlePointer p) {
        derived().initializeInfoOutImpl(wp, nP, p);
    };
    inline void cloneILImpl(interWorkParticle *wp) { derived().cloneILImpl(wp); };
    inline void processWorkParticleImpl(interWorkParticle *wp) {
        derived().processWorkParticleImpl(wp);
    };
    inline void setupCacheImpl() { derived().setupCacheImpl(); };
    inline int queueInteractionGPUImpl(interWorkParticle *wp, ILtile tile) {
        return derived().queueInteractionGPUImpl(wp, tile);
    };
    inline void particleWorkDoneImpl(interWorkParticle *wp) {
        derived().particleWorkDoneImpl(wp);
    };
    inline bool particleSkipConditionImpl(particleStore::ParticlePointer p) {
        return derived().particleSkipConditionImpl(p);
    };
    inline void evalInteractionImpl(const PINFOIN &Part,
                                    ILtile &tile,
                                    PINFOOUT &out) {
        derived().evalInteractionImpl(Part, tile, out);
    };
    inline std::tuple<int, double> interactImpl(treeStore::NodePointer pkdn,
            int iRoot1, int iRoot2) {
        return derived().interactImpl(pkdn, iRoot1, iRoot2);
    };
    inline Derived &derived() { return static_cast<Derived &>(*this); };
};

// I'm using traits classes for the interactions.
// These contains all the types information used for
// the class, and is passed here as a template parameter
// so that the base class can use the derived class'
// types. All this at compile-time, so no run-time
// overhead due to virtualization or anything.
template <typename Traits, class Derived>
class BaseInteraction : public BaseInteractionInterface<Traits, Derived> {
    using impl = BaseInteractionInterface<Traits, Derived>;

public:
    // Access the derived classes types via trait class
    using IL = typename Traits::IL;
    using ILtile = typename Traits::IL::tile_type;
    using Options = typename Traits::Options;
    using PINFOIN = typename Traits::PINFOIN;
    using PINFOOUT = typename Traits::PINFOOUT;
    using name = typename Traits::name;
    using ResultInteraction = typename Traits::template ResultInteraction<fvec>;
    using interWorkParticle = typename impl::interWorkParticle;

    // Default constructor without argument is needed when
    // building the service. We do a delayed initialization
    // later in setup(), which is called at Service runtime,
    // not construction.
    BaseInteraction() {};

    void setupImpl(PKD pkd, Options *options, WalkStats *stats) {
        this->pkd = pkd;
        this->registerCudaCallback(*pkd->cudaClient);
        this->options = options;
        this->walkStats = stats;
        this->nTilesTotal = 0;
        this->nTilesCPU = 0;
    }

    // Generic interfaces to the interaction list(s).
    // If the interaction has more than one, these
    // functions should be overloaded accordingly.
    void clearILImpl() { il->clear(); };
    int countILImpl() { return il->count(); };
    void setReferenceILImpl(blitz::TinyVector<double, 3> ref) {
        il->setReference(ref);
    };
    blitz::TinyVector<double, 3> getReferenceILImpl() {
        return il->getReference();
    };

    // Test that the IL(s) has been populated correctly, only
    // for debugging runs
    void extensiveILTestImpl(interWorkParticle *wp) {};

protected:
    PKD pkd;
    Options *options;
    WalkStats *walkStats;
    int nTilesTotal;
    int nTilesCPU;

    blitz::TinyVector<double, 3> kdn_r;

    // By default we only have one interaction list,
    // but derived classes can manually add more
    // if required
    IL *il;

    // Allow for some flexibility when initializing this at it
    // could contain information of the interaction (e.g., IL)
    void initializeInfoInImpl(interWorkParticle *wp, int nP,
                              particleStore::ParticlePointer p) {
        wp->pInfoIn[nP].init(p);
    };
    void initializeInfoOutImpl(interWorkParticle *wp, int nP,
                               particleStore::ParticlePointer p) {
        wp->pInfoOut[nP].init(p);
    };

    // This hook allows or explicit cloning of the IL
    // if needed (it seems to be used when bDensity & bGPU)
    void cloneILImpl(interWorkParticle *wp) { }

    // Allow derived classes to start a cache with a different
    // helper to allow for packing, for example
    void setupCacheImpl() {
        pkd->mdl->CacheInitialize(CID_PARTICLE,
                                  std::make_shared<particleCACHEhelper>(pkd));
    }


    // This would wrap all the queues or any processing
    // required by a derived class.
    // By default just queue it.
    void processWorkParticleImpl(interWorkParticle *wp) {
        this->queueInteraction(wp);
    };

    int queueInteractionGPUImpl(interWorkParticle *wp, ILtile tile) {
        return 0;
        // if (options->bGPU) {
        // #ifdef USE_CUDA
        //     if (pkd->cudaClient->queuePP(wp, tile))
        //         continue;
        // #endif
        // }
    }

    void queueInteractionImpl(interWorkParticle *wp) {
        for (auto &tile : *this->il) {
            ++this->nTilesTotal;
            if (this->queueInteractionGPU(wp, tile))
                continue;
            ++this->nTilesCPU;
            for (auto i = 0; i < wp->nP; ++i) {
                this->evalInteraction(wp->pInfoIn[i], tile, wp->pInfoOut[i]);
            }
        }
    }

    // Acts on own interWorkParticle, overloaded by the
    // different choices in pkdParticleWorkDone
    void particleWorkDoneImpl(interWorkParticle *wp) {
        this->walkStats->dFlopSingleCPU += wp->dFlopSingleCPU;
        this->walkStats->dFlopSingleGPU += wp->dFlopSingleGPU;
        this->walkStats->dFlopDoubleCPU += wp->dFlopDoubleCPU;
        this->walkStats->dFlopDoubleGPU += wp->dFlopDoubleGPU;
        this->walkStats->nTilesTotal += this->nTilesTotal;
        this->walkStats->nTilesCPU += this->nTilesCPU;
        this->walkStats->nActive += wp->nP;
        this->walkStats->dFlop += wp->dFlop;
    };

    // Condition for not processing a particle
    bool particleSkipConditionImpl(particleStore::ParticlePointer p) {
        return !p->is_rung_range(options->uRungLo, options->uRungHi);
    };

    void evalInteractionImpl(const PINFOIN &Part, ILtile &tile,
                             PINFOOUT &out) {
        assert(0); // evalInteraction should be overloaded
    };

public:
    // pkdGravInteract
    int interactImpl(treeStore::NodePointer pkdn,
                     int iRoot1, int iRoot2) {
        int nActive = 0;
        int nP;

        this->kdn_r = pkdn->position();

        /*
        ** Now process the two interaction lists for each active particle.
        */

        /* Collect the bucket particle information */
        auto wp = new interWorkParticle(pkd, pkdn->count(), il);
        assert(wp != NULL);

        ++pkd->nWpPending;

        for (auto &p : *pkdn) {
            if (this->particleSkipCondition(&p))
                continue;

            auto r = p.position();

            nP = wp->nP++;
            wp->pPart[nP] = &p;
            wp->iPart[nP] = &p - pkd->particles.begin();
            this->initializeInfoIn(wp, nP, &p);
            this->initializeInfoOut(wp, nP, &p);

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
        ** we need to clone the il,
        ** otherwise we use the reference to the global il.
        */
        this->cloneIL(wp);

        nActive += wp->nP;

        if (options->doExtensiveTest) {
            this->extensiveILTest(wp);
        }

        this->processWorkParticle(wp);


        this->particleWorkDone(wp);

        --pkd->nWpPending;
        delete wp;
        return nActive;
    };
};

#endif
