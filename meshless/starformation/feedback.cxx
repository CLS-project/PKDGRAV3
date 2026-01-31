#include <numeric>
#include "hydro/hydro.h"
#include "feedback.h"
#include "imf.h"
struct SmoothSNFeedback : smooth::policy::OnlyStar {
    double dTime;
    double dCCSNFBDelay;
    double dSNIaFBDelay;

    inline auto active(particleStore::ParticleReference& p) -> bool {
        if (p.is_star()) {
            auto& star = p.star();
            if (!star.bCCSNFBDone && ((dTime - star.fTimer) > dCCSNFBDelay)) return true;
            if (!star.bSNIaFBDone && ((dTime - star.fTimer) > dSNIaFBDelay)) return true;
        }
        return false;
    }
};

auto ServiceSNFeedback::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->pkd();
    auto helper = std::make_shared<SmoothHelperSNFeedback>(pkd, in->dTime, in->dCCSNFBDelay, in->dSNIaFBDelay, in->dSNFBDu, in->dCCSNFBSpecEnergy,
                                                           in->dSNIaFBSpecEnergy, in->dConstGamma);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
    *out = smx.Smooth<SmoothSNFeedback>({
        .dTime = in->dTime,
        .dCCSNFBDelay = in->dCCSNFBDelay,
        .dSNIaFBDelay = in->dSNIaFBDelay,
    });
    return sizeof(output);
}

static inline void snFeedback(PKD pkd, PARTICLE* pIn, std::span<NN> nnList, const float fNgbTotMass, const double dDeltau, const double dStarMass,
                              const double dSpecEnergy, const float fEfficiency, const double dConstGamma);

/* Function that will be called with the information of all the neighbors.
 * The helper snFeedback computes the probability of explosion and adds
 * the energy to the gas particles in nnList
 */

void SmoothHelperSNFeedback::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    const auto fNgbTotMass =
        std::accumulate(nnList.begin(), nnList.end(), 0.0f, [this](auto a, const auto& b) { return a + pkd->particles[b.pPart].mass(); }) - p.mass();

    auto& star = p.star();

#ifdef STELLAR_EVOLUTION
    const double dStarMass = star.fInitialMass;
#else
    const double dStarMass = p.mass();
#endif

    if (!star.bCCSNFBDone && ((dTime - star.fTimer) > dCCSNFBDelay)) {
        snFeedback(pkd, &p, nnList, fNgbTotMass, dSNFBDu, dStarMass, dCCSNFBSpecEnergy, star.fSNEfficiency, dConstGamma);
        star.bCCSNFBDone = 1;
    }

    if (!star.bSNIaFBDone && ((dTime - star.fTimer) > dSNIaFBDelay)) {
        snFeedback(pkd, &p, nnList, fNgbTotMass, dSNFBDu, dStarMass, dSNIaFBSpecEnergy, star.fSNEfficiency, dConstGamma);
        star.bSNIaFBDone = 1;
    }
}

void SmoothHelperSNFeedback::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->fMass = p2.mass();
#ifdef ENTROPY_SWITCH
        p1->fDensity = p2.density();
#endif
    }
}

void SmoothHelperSNFeedback::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_adj_decorator();
        p1.set_position(p2->position);
        p1.set_mass(p2->fMass);
#ifdef ENTROPY_SWITCH
        p1.set_density(p2->fDensity);
#endif
    }
}

void SmoothHelperSNFeedback::init(void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];

    if (p.is_gas()) {
        auto& sph = p.sph();

        sph.fAccFBEnergy = 0.;
    }
}

void SmoothHelperSNFeedback::flush(void* dst, const void* src) {
    auto p1 = static_cast<flush_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    if (p2.is_gas()) {
        const auto& sph = p2.sph();

        p1->fAccFBEnergy = sph.fAccFBEnergy;
    }
}

void SmoothHelperSNFeedback::combine(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const flush_type*>(src);

    if (p1.is_gas()) {
        auto& sph = p1.sph();

        sph.fAccFBEnergy += p2->fAccFBEnergy;
    }
}

static inline void snFeedback(PKD pkd, PARTICLE* pIn, std::span<NN> nnList, const float fNgbTotMass, const double dDeltau, const double dStarMass,
                              const double dSpecEnergy, const float fEfficiency, const double dConstGamma) {

    const double dProb = fEfficiency * dSpecEnergy * dStarMass / (dDeltau * fNgbTotMass);
    const double dCorrFac = dProb > 1.0 ? dProb : 1.0;

    for (auto& nn : nnList) {
        if (nn.pPart == pIn) continue;

        if (rand() < RAND_MAX * dProb) {
            auto q = pkd->particles[nn.pPart];
            auto& qsph = q.sph();
            const double dEnergyInput = dCorrFac * dDeltau * q.mass();

            qsph.fAccFBEnergy += dEnergyInput;
        }
    }
}
