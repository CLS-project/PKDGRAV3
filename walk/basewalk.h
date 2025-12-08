#ifndef BASE_WALK_H
#define BASE_WALK_H
#include "../TraversePST.h"
#include "gravity/cl.h"
#include "pkd.h"
#include "pkd_config.h"
#include "smooth/smooth.h"

template <auto BasePtr, auto DerivedPtr> constexpr bool is_overridden() {
    using BaseType = decltype(BasePtr);
    using DerivedType = decltype(DerivedPtr);

    // force identical pointer-to-member type
    using Common = BaseType;

    return static_cast<Common>(DerivedPtr) != static_cast<Common>(BasePtr);
}

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

template <class Derived>
class BaseTreeWalkInterface : public TraverseCombinePST {
  public:
    explicit BaseTreeWalkInterface(PST pst, std::string_view service_name, int nInBytes, int nOutBytes)
        : TraverseCombinePST(pst, service_name, nInBytes, nOutBytes) {
    };

    inline void addChild(PKD pkd, int iCache, clList& cl, int iChild, int id,
                         blitz::TinyVector<float, 3> fOffset) {
        derived().addChildImpl(pkd, iCache, cl, iChild, id, fOffset);
    }
    inline void processCheckList(PKD pkd, int iRoot, int iRoot2) {
        derived().processCheckListImpl(pkd, iRoot, iRoot2);
    }
    inline void doneChecklist() { derived().doneChecklistImpl(); };
    inline void initCLAccumulatedVariables() {
        initCLAccumulatedVariablesImpl();
    };
    inline void beforeCLProcessing(treeStore::NodePointer k) {
        beforeCLProcessingImpl(k);
    };
    inline void calculateLocalExpansion(treeStore::NodePointer k) {
        calculateLocalExpansionImpl(k);
    }
    inline void updateStackMomentTensor(treeStore::NodePointer c,
                                        double xParent, double yParent,
                                        double zParent) {
        updateStackMomentTensorImpl(c, xParent, yParent, zParent);
    };
    inline void shiftToNextNode(blitz::TinyVector<double, 3> k_r,
                                double xParent, double yParent,
                                double zParent) {
        shiftToNextNodeImpl(k_r, xParent, yParent, zParent);
    };
    inline void updateMomentTensorFromStack() {
        updateMomentTensorFromStackImpl();
    };
    inline void openOutcomeBlock(treeStore::NodePointer k, int n,
                                 clBlock& blk) {
        openOutcomeBlockImpl(k, n, blk);
    };
    inline void selectOutcome(clBlock blk, int outcome, int jTile) {
        selectOutcomeImpl(blk, outcome, jTile);
    };
    inline bool testAddNodeCL(treeStore::NodePointer k) {
        return testAddNodeCLImpl(k);
    };
    inline void walk() { walkImpl(); };
    inline void initWalk() { initWalkImpl(); };
    inline void doneWalk() { doneWalkImpl(); };
    inline void setupCache() { setupCacheImpl(); };
    inline void initStats() { initStatsImpl(); };
    inline void computeStats() { computeStatsImpl(); };
    inline std::tuple<int, double> interact_wrapper(treeStore::NodePointer k,
                                                    int iRoot, int iRoot2) {
        return interact_wrapperImpl(k, iRoot, iRoot2);
    };
    int Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
        return ServiceImpl(pst, vin, nIn, vout, nOut);
    };
    virtual int Combine(void* vout, void* vout2) {
      return CombineImpl(vout, vout2);
    }

  protected:
    inline void addChildImpl(PKD pkd, int iCache, clList& cl, int iChild,
                             int id, blitz::TinyVector<float, 3> fOffset) {
        derived().addChildImpl(pkd, iCache, cl, iChild, id, fOffset);
    };
    inline void processCheckListImpl(PKD pkd, int iRoot, int iRoot2) {
        derived().processCheckListImpl(pkd, iRoot, iRoot2);
    }
    inline void doneChecklistImpl() { derived().doneChecklistImpl(); };
    inline void initCLAccumulatedVariablesImpl() {
        derived().initCLAccumulatedVariablesImpl();
    };
    inline void beforeCLProcessingImpl(treeStore::NodePointer k) {
        derived().beforeCLProcessingImpl(k);
    };
    inline void calculateLocalExpansionImpl(treeStore::NodePointer k) {
        derived().calculateLocalExpansionImpl(k);
    }
    inline void updateStackMomentTensorImpl(treeStore::NodePointer c,
                                            double xParent, double yParent,
                                            double zParent) {
        derived().updateStackMomentTensorImpl(c, xParent, yParent, zParent);
    };
    inline void shiftToNextNodeImpl(blitz::TinyVector<double, 3> k_r,
                                    double xParent, double yParent,
                                    double zParent) {
        derived().shiftToNextNodeImpl(k_r, xParent, yParent, zParent);
    };
    inline void updateMomentTensorFromStackImpl() {
        derived().updateMomentTensorFromStackImpl();
    };
    inline void openOutcomeBlockImpl(treeStore::NodePointer k, int n,
                                     clBlock& blk) {
        derived().openOutcomeBlockImpl(k, n, blk);
    };
    inline void selectOutcomeImpl(clBlock blk, int outcome, int jTile) {
        derived().selectOutcomeImpl(blk, outcome, jTile);
    };
    inline bool testAddNodeCLImpl(treeStore::NodePointer k) {
        return derived().testAddNodeCLImpl(k);
    };
    inline void walkImpl() { derived().walkImpl(); };
    inline void initWalkImpl() { derived().initWalkImpl(); };
    inline void doneWalkImpl() { derived().doneWalkImpl(); };
    inline void setupCacheImpl() { derived().setupCacheImpl(); };
    inline void initStatsImpl() { derived().initStatsImpl(); };
    inline void computeStatsImpl() { derived().computeStatsImpl(); };
    inline std::tuple<int, double>
    interact_wrapperImpl(treeStore::NodePointer k, int iRoot, int iRoot2) {
        return derived().interact_wrapperImpl(k, iRoot, iRoot2);
    };
    inline int ServiceImpl(PST pst, void* vin, int nIn, void* vout, int nOut) {
        return derived().ServiceImpl(pst, vin, nIn, vout, nOut);
    };
    inline int CombineImpl(void* vout, void* vout2) {
      return derived().CombineImpl(vout, vout2);
    }
    inline Derived& derived() { return static_cast<Derived&>(*this); }

};

// As with the interaction, we define traits for the
// derived tree walks to that we can access their types
// from the base class.
template <typename Traits, typename Interaction, class Derived>
class BaseTreeWalk : public BaseTreeWalkInterface<Derived> {
    // We choose to derive TraverseCombinePST to allow gathering
    // of information of the underlying walk/interactions, such
    // as the statistics shown for the gravity.
    using impl = BaseTreeWalkInterface<Derived>;
    friend impl;

  public:
    using Options = typename Traits::Options;
    using output = typename Traits::output;
    using name = typename Traits::name;
    using Base = BaseTreeWalk<Traits, Interaction, Derived>;

    // For passing both the tree walk options and the
    // interaction options at the PST level, we pack them
    // in a single struct for convenience
    struct OptionPack {
        Options walk;
        typename Interaction::Options inter;
    };

    // The service name is built by concatenating
    // at compile time the names of the walk+interaction
    // defined in their respective traits
    static constexpr auto service_name =
        Concat<typename Traits::name, typename Interaction::name>::value;

    explicit BaseTreeWalk(PST pst)
        : BaseTreeWalkInterface<Derived>(pst, service_name.data(), sizeof(OptionPack),
                             sizeof(output)) {
        // cl = new CL();
        clStack.reserve(40); // Start with a sensible stack size
    };

  protected:
    PKD pkd;
    SMX smx;
    Interaction interaction;
    Options* options;
    output* out;

    clList* clSib;
    clList* clNew;

    int bReferenceFound;

    void initStat(STAT* ps, int id) {
        ps->idMax = id;
        ps->dMax = ps->dSum;
        ps->dSum2 = ps->dSum * ps->dSum;
        ps->n = 1; /* sometimes we will need to zero this */
    }

    void combStat(STAT* ps, STAT* pa) {
        if (pa->dMax > ps->dMax) {
            ps->idMax = pa->idMax;
            ps->dMax = pa->dMax;
        }
        ps->dSum += pa->dSum;
        ps->dSum2 += pa->dSum2;
        ps->n += pa->n;
    }

  protected:
    // IA: We still use the clList defined in PKD, I do not know why
    // is done that way and if we can simply define and allocate
    // one for us here.
    // clList *cl;
    std::vector<CSTACK> clStack;

    // Any variable that must be accumulated during the walk
    double dPartSum;
    double dCellSum;
    int nActive;
    double dFlopTotal;

    void addChildImpl(PKD pkd, int iCache, clList& cl, int iChild, int id,
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
        const auto& SPHbob = c->have_BOB() ? c->BOB() : SPHBOB();
        cl.append(iCache, id, iChild, idLower, iLower, idUpper, iUpper, nc,
                  cOpen, c->mass(), 4.0f * c->fSoft2(), c_r, fOffset, cbnd,
                  SPHbob);
    }

    void processCheckListImpl(PKD pkd, int iRoot, int iRoot2) {
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
        double dFlop;

        pkd->dFlop = 0.0; /* Flops are accumulated here! */
        clStack.clear();

        /*
        ** Clear monopole sentinel and local expansion and the timestepping
        *sums.
        */
        nCLActive = 0;
        impl::initCLAccumulatedVariables();

        /*
        ** We are now going to work on the local tree.
        ** Make iCell point to the root of the tree again.
        */
        auto k = pkd->tree[iCell = iRoot];
        k_r = k->position();
        printf("stuck??\n");
        while (1) {
            bReferenceFound = 0;
            while (1) {
                /*
                ** Process the Checklist.
                */
                impl::beforeCLProcessing(k);
                /*
                ** For cells which will remain on the checklist for a further
                *deeper level of
                ** level of the tree we will be using the stack cl (even if no
                *push happens later).
                */
                this->clSib = new clList(pkd->clFreeList);
                do {
                    for (auto& tile : pkd->cl) {
                        openOutcomeSIMD(k, tile);
                    }
                    this->clNew = new clList(pkd->clFreeList);
                    for (auto& tile : pkd->cl) {
                        auto nBlocks = tile.count() / tile.width;
                        for (auto iBlock = 0; iBlock <= nBlocks; ++iBlock) {
                            int n = iBlock < nBlocks
                                        ? tile.width
                                        : tile.count() - nBlocks * tile.width;
                            auto& blk = tile[iBlock];
                            for (jTile = 0; jTile < n; ++jTile) {
                              impl::selectOutcome(blk, blk.iOpen[jTile], jTile);
                            }
                        }
                    } /* end of CL_LOOP */
                    pkd->cl = std::move(*clNew);
                } while (pkd->cl.count());
                /*
                ** Now calculate the local expansion.
                */
                impl::calculateLocalExpansion(k);
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
                for (auto& p : *k) {
                    if (p.is_marked())
                        assert(k->is_marked());
                }
                /*
                ** Check iCell is active. We eventually want to just to a
                ** rung check here when we start using tree repair, but
                ** for now this is just as good.
                */
                if (impl::testAddNodeCL(k)) {
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
                    if (impl::testAddNodeCL(c)) {
                        // if ((c->min_rung()<=ts->uRungHi && c->max_rung() >=
                        // ts->uRungLo)
                        // || (SPHoptions->useDensityFlags && c->is_marked()) ||
                        // (SPHoptions->useNNflags && c->is_NN()) ||
                        // (SPHoptions->doMFM && c->is_marked())) {
                        /*
                        ** Sibling is active so we need to clone the checklist!
                        */
                        pkd->cl = *clSib;
                        /*
                        ** Sibling is active as well.
                        ** Push Checklist for the sibling onto the stack
                        ** before proceeding deeper in the tree.
                        */
                        clStack.emplace_back(pkd->clFreeList);
                        clStack.back().cl = std::move(*clSib);
                        clStack.back().iNodeIndex = iSib;

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
                        impl::updateStackMomentTensor(c, xParent, yParent, zParent);
                    } else {
                        /*
                        ** Otherwise we can simple grab it.
                        */
                        pkd->cl = std::move(*clSib);
                    }
                } else {
                    /*
                    ** Skip iCell, but add it to the Checklist.
                    ** No need to push anything onto the stack here.
                    ** We can simply grab it the checklist from the Stack.
                    */
                    pkd->cl = std::move(*clSib);
                    /*
                    ** Move onto processing the sibling.
                    */
                    k = pkd->tree[++iCell];
                    k_r = k->position();
                }

                    impl::shiftToNextNode(k_r, xParent, yParent, zParent);
            }
            /*
            ** Now the interaction list should be complete and the
            ** Checklist should be empty! Calculate gravity on this
            ** Bucket!
            */
            auto [nActiveInt, flops] = impl::interact_wrapper(k, iRoot, iRoot2);

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
            if (clStack.empty())
                goto doneCheckList;
            k = pkd->tree[iCell = clStack.back().iNodeIndex];
            k_r = k->position();
            /*
            ** Pop the Checklist from the top of the stack,
            ** also getting the state of the interaction list.
            */
            interaction.clearIL();
            /*
            ** Grab the checklist from the stack.
            */
            assert(pkd->cl.count() == 0);
            pkd->cl = std::move(clStack.back().cl);
            impl::updateMomentTensorFromStack();
            clStack.pop_back();
        }
    doneCheckList:
        assert(0);
        impl::doneChecklist();
        nActive += nCLActive;
    }

    void doneChecklistImpl() {
#ifdef USE_CUDA
        pkd->cudaClient->flushCUDA();
#endif
#ifdef USE_METAL
        pkd->metalClient->flushMETAL();
#endif
        pkd->mdl->CompleteAllWork();
        dFlopTotal += pkd->dFlop; /* Accumulate work flops (notably Ewald) */
    }

  protected:
    // Convenient hooks when processing the CL to add gravity-related
    // operations.
    void initCLAccumulatedVariablesImpl() {}

    void beforeCLProcessingImpl(treeStore::NodePointer k) {}

    void calculateLocalExpansionImpl(treeStore::NodePointer k) {}

    void updateStackMomentTensorImpl(treeStore::NodePointer c, double xParent,
                                     double yParent, double zParent) {}

    void updateMomentTensorFromStackImpl() {}

    void shiftToNextNodeImpl(blitz::TinyVector<double, 3> k_r, double xParent,
                             double yParent, double zParent) {}

    // Opening criteria
    void openOutcomeBlockImpl(treeStore::NodePointer k, int n, clBlock& blk) {}

    void openOutcomeSIMD(treeStore::NodePointer k, clTile& tile) {
        // Do the full blocks, followed by the last (probably not full) block
        auto nBlocks = tile.count() / tile.width;
        for (auto iBlock = 0; iBlock < nBlocks; ++iBlock) {
            impl::openOutcomeBlock(k, tile.width, tile[iBlock]);
        }
        impl::openOutcomeBlock(k, tile.count() - nBlocks * tile.width,
                               tile[nBlocks]);
    }

    // Based on openOutcomeBlock do something.
    // We split it so allow for child classes
    // to reuse the same criteria easily
    void selectOutcomeImpl(clBlock blk, int outcome, int jTile) {}

    // Whether a node is eligible for adding
    // it to the interaction list
    bool testAddNodeCLImpl(treeStore::NodePointer k) { return true; }

    // pkdGravWalk
    void walkImpl() {
        int id;
        float fOffset[3];
        int ix, iy, iz;
        int iTop1, iTop2;

        impl::initWalk();

        iTop1 = pkd->iTopTree[options->iRoot1];
        iTop2 = pkd->iTopTree[options->iRoot2];
        id = pkd->Self();

        /*
        ** Walk tree 1 against trees 1 (and optionally 2) if there are active
        *particles
        */
        auto k = pkd->tree[iTop1];
        if (k->count() && impl::testAddNodeCL(k)) {
            /*
            ** Initially we set our cell pointer to
            ** point to the top tree.
            */
            interaction.clearIL();
            pkd->cl.clear();

            /*
            ** Add all replicas of the entire box to the Checklist.
            ** We add at least one box (0,0,0). The root cell is alway on
            *processor 0.
            */
            for (ix = -options->nReps; ix <= options->nReps; ++ix) {
                fOffset[0] = ix * pkd->fPeriod[0];
                for (iy = -options->nReps; iy <= options->nReps; ++iy) {
                    fOffset[1] = iy * pkd->fPeriod[1];
                    for (iz = -options->nReps; iz <= options->nReps; ++iz) {
                        fOffset[2] = iz * pkd->fPeriod[2];
                        impl::addChild(pkd, CID_CELL, pkd->cl, iTop1, id,
                                       fOffset);
#ifndef SINGLE_CACHES
                        if (iTop2 > 0)
                            impl::addChild(pkd, CID_CELL2, pkd->cl, iTop2, id,
                                           fOffset);
#else
                        if (iTop2 > 0)
                            impl::addChild(pkd, CID_CELL, pkd->cl, iTop2, id,
                                           fOffset);
#endif
                    }
                }
            }
            impl::processCheckList(pkd, iTop1, iTop2);
        }
        impl::doneWalk();
    }

    // initWalkGrav
    void initWalkImpl() {}

    // doneWalkGrav
    void doneWalkImpl() {}

    // Allow derived classes to start a cache with a different helper to allow
    // for packing, for example
    void setupCacheImpl() {
        pkd->mdl->CacheInitialize(CID_PARTICLE,
                                  std::make_shared<particleCACHEhelper>(pkd));
    }

    // Initialize flop/counts statistics for latter analisis
    void initStatsImpl() {}

    // Process the gatheres statistics during the walk (which should be
    // available as private variables)
    void computeStatsImpl() {}

    // We abstract this a bit to allow passing different interaction lists
    std::tuple<int, double> interact_wrapperImpl(treeStore::NodePointer k,
                                                 int iRoot, int iRoot2) {

        return interaction.interact(k, iRoot, iRoot2);
    }

    static_assert(std::is_void<Options>() ||
                  std::is_standard_layout<Options>());
    static_assert(std::is_void<output>() || std::is_standard_layout<output>());

  public:
    // pkdGravAll
    int ServiceImpl(PST pst, void* vin, int nIn, void* vout, int nOut) {
        this->pkd = pst->plcl->pkd;
        auto* optpack = static_cast<OptionPack*>(vin);
        typename Interaction::Options opt = optpack->inter;
        this->options = &optpack->walk;
        this->out = static_cast<output*>(vout);
        assert(sizeof(this->options) == nIn);
        interaction.setup(this->pkd, &opt);

        // options->copyToDevice(); // or something similar
        impl::setupCache();

        impl::initStats();

        impl::walk();

        impl::computeStats();

        pkd->mdl->FinishCache(CID_PARTICLE);
        // This kind of thing should go inside doneWalk
        return sizeof(output);
    }

    // pstGravity combination of statistics
    virtual int CombineImpl(void* vout, void* vout2) {
        return sizeof(output);
    }
};
#endif
