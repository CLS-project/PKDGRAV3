#include "pkd.h"
#include "smooth/smooth.h"

static inline void pkdAddFBEnergy(PKD pkd, particleStore::ParticleReference& p, decorators::MFM* psph, double dConstGamma) {
    psph->Uint += psph->fAccFBEnergy;
    psph->E += psph->fAccFBEnergy;
#ifdef ENTROPY_SWITCH
    psph->S += psph->fAccFBEnergy * (dConstGamma - 1.) * pow(p.density(), -dConstGamma + 1);
#endif
    psph->fAccFBEnergy = 0.0;
}

struct snFeedbackPack {
    blitz::TinyVector<double, 3> position;
    float fMass;
#ifdef ENTROPY_SWITCH
    float fDensity;
#endif
    uint8_t iClass;
};

struct snFeedbackFlush {
    float fAccFBEnergy;
};

class SmoothHelperSNFeedback : public SmoothHelperPacked<snFeedbackPack, snFeedbackFlush> {
    double dTime, dCCSNFBDelay, dSNIaFBDelay, dSNFBDu, dCCSNFBSpecEnergy, dSNIaFBSpecEnergy, dConstGamma;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    void pack(void* dst, const void* src) final;
    void unpack(void* dst, const void* src, const void* key) final;
    void init(void* dst) final;
    void flush(void* dst, const void* src) final;
    void combine(void* dst, const void* src, const void* key) final;
    void init_shared(void* dst) final { init_shared_gas_decorator(pkd, dst); }

  public:
    SmoothHelperSNFeedback(pkdContext* pkd, double dTime, double dCCSNFBDelay, double dSNIaFBDelay, double dSNFBDu, double dCCSNFBSpecEnergy,
                           double dSNIaFBSpecEnergy, double dConstGamma)
        : SmoothHelperPacked<snFeedbackPack, snFeedbackFlush>(pkd, true), dTime(dTime), dCCSNFBDelay(dCCSNFBDelay), dSNIaFBDelay(dSNIaFBDelay),
          dSNFBDu(dSNFBDu), dCCSNFBSpecEnergy(dCCSNFBSpecEnergy), dSNIaFBSpecEnergy(dSNIaFBSpecEnergy), dConstGamma(dConstGamma) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
    auto range_neighbour(treeStore::NodePointer node) -> std::pair<int, int> override { return SmoothOnlyGas::range_neighbour(node); }
    auto accept_neighbour(particleStore::ParticleReference& p) -> bool override { return SmoothOnlyGas::accept_neighbour(p); }
};

static inline float SNFeedbackEfficiency(float Z, float rho, double dSNFBEfficiency, double dSNFBMaxEff, double dSNFBEffIndex, double dSNFBEffnH0) {
    if (dSNFBMaxEff > 0.0) {
        const double den = 1.0 + pow(Z / 0.00127, dSNFBEffIndex) * pow(rho / dSNFBEffnH0, -dSNFBEffIndex);
        return dSNFBEfficiency + (dSNFBMaxEff - dSNFBEfficiency) / den;
    }
    else {
        return dSNFBEfficiency;
    }
}
