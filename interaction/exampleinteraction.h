#include "baseinteraction.h"

// A interaction definition is separated in three blocks inside
// this file. The first block is the interaction "Traits".
// These define the datatypes used for the interactions and its
// name.
//
// The second block is the evaluation of the interaction, which
// contains all the interesting physics.
// This is done outside of the class, and in most complex
// scenarios will probably be outside of this header file,
// but for the simple example here it can stay here
// The interface for this part has not changed drastically,
// thus already working interactions should be easy to migrate.
//
// Lastly, we have the interaction class itself.
// It overloads certain functions (without "virtual"), so there
// is no performance overhead at run-time.
// In this simple example there is nothing fancy going on there.

/*
 * INTERACTION TRAITS
 */
class ExampleInteractionTraits {
public:
    struct name {
        static constexpr std::string_view value = "ExampleInteraction";
    };

    struct Options {
        float time;
        float delta;

        bool testoption;

        uint8_t uRungLo;
        uint8_t uRungHi;

        bool doExtensiveTest;
        bool bGPU;
    };

    typedef struct PINFOIN {
        blitz::TinyVector<float, 3> r;

        void init(particleStore::ParticlePointer p) { r = p->position(); }
    } PINFOIN;

    typedef struct PINFOOUT {
        int n;

        void init(particleStore::ParticlePointer p) { n = 0; }
    } PINFOOUT;

    using IL = ilpList;

    template <class F> struct ResultInteraction {
        F n;

        void zero() { n = 0; }
        ResultInteraction operator+=(const ResultInteraction a) {
            n += a.n;
            return *this;
        }
    };
};

/*
 * INTERACTION EVALUATION
 */
template <class F, class M>
/*PP_CUDA_BOTH*/ ExampleInteractionTraits::ResultInteraction<F>
Evaluate(F fx, F fy, F fz, F Idx, F Idy, F Idz, bool testoption) {
    F dx = fx + Idx;
    F dy = fy + Idy;
    F dz = fz + Idz;
    F dr = sqrt(dx * dx + dy * dy + dz * dz);

    ExampleInteractionTraits::ResultInteraction<F> z;
    z.zero();

    M mask = dr > 0.0;
    F one = 1.0;
    mask_mov(z.n, mask, one);

    return z;
};

template <typename BLOCK>
struct ilist::EvalBlock<ExampleInteractionTraits::ResultInteraction<fvec>,
           BLOCK> {
    typedef ExampleInteractionTraits::ResultInteraction<fvec> result_type;
    const fvec fx, fy, fz;
    const ExampleInteractionTraits::Options *const options;
    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy, fvec fz,
              ExampleInteractionTraits::Options *options)
        : fx(fx), fy(fy), fz(fz), options(options) {}

    result_type operator()(int n, BLOCK &blk) {
        n /= fvec::width(); // Now number of blocks

        result_type result;
        result.zero();
        for (auto i = 0; i < n; ++i) {

            result +=
                Evaluate<fvec, fmask>(fx, fy, fz, blk.dx.v[i], blk.dy.v[i],
                                      blk.dz.v[i], options->testoption);
        }
        return result;
    }
};

/*
 * INTERACTION CLASS
 */
class ExampleInteraction : public BaseInteraction<ExampleInteractionTraits> {
private:
    // How the result of the interaction is aggregated to the
    // per-particle output information
    void aggregateResult(PINFOOUT &out, ResultInteraction &result) {
        out.n += int(hadd(result.n));
    };

    void queueInteraction(workParticle *wp, IL *il) {
        for (auto &tile : *il) {
            ++pkd->nTilesTotal;
            // TODO: how to do cuda queues now?
            // if (options->bGPU) {
            // #ifdef USE_CUDA
            //     if (pkd->cudaClient->queueMFM(wp,tile)) continue;
            // #endif
            // }
            ++pkd->nTilesCPU;
            for (auto i = 0; i < wp->nP; ++i) {
                evalInteraction(wp->pInfoIn[i], tile, wp->pInfoOut[i]);
            }
        }
    }

    void evalInteraction(const PINFOIN &Part, IL::tile_type &tile,
                         PINFOOUT &out) {
        ilist::EvalBlock<ExampleInteraction::ResultInteraction, IL::block_type>
        eval(Part.r[0], Part.r[1], Part.r[2], &options);
        auto result = EvalTile(tile, eval);
        aggregateResult(out, result);
    }
};
