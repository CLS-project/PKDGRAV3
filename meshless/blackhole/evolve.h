#include "smooth/smooth.h"
#include "pkd.h"

#define NOT_ACCRETED -1

struct bhEvolvePack {
    blitz::TinyVector<double, 3> position;
    blitz::TinyVector<double, 3> velocity;
    float c;
    float fMass;
#ifdef ENTROPY_SWITCH
    float fDensity;
#endif
    uint8_t iClass;
};

struct bhEvolveFlush {
#ifdef OLD_FB_SCHEME
    double E;
    double Uint;
    #ifdef ENTROPY_SWITCH
    double S;
    #endif
#else // OLD_FB_SCHEME
    float fAccFBEnergy;
#endif
    int iAccPid;
    int iAccIndex;
};

class SmoothHelperBHEvolve : public SmoothHelperPacked<bhEvolvePack, bhEvolveFlush> {
    double a_, dBHAccretionAlpha, dBHAccretionCvisc, dBHAccretionEddFac, dConstGamma, dBHFBEff, dBHFBEcrit;
    bool bBHAccretion, bBHFeedback;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;
    virtual void init(void* dst) override final;
    virtual void flush(void* dst, const void* src) override final;
    virtual void combine(void* dst, const void* src, const void* key) override final;
    virtual void init_shared(void* dst) override final { init_shared_bh_decorator(pkd, dst); }

  public:
    SmoothHelperBHEvolve(pkdContext* pkd, double a, double dBHAccretionAlpha, double dBHAccretionCvisc, double dBHAccretionEddFac, bool bBHAccretion,
                         bool bBHFeedback, double dConstGamma, double dBHFBEff, double dBHFBEcrit)
        : SmoothHelperPacked<bhEvolvePack, bhEvolveFlush>(pkd, true), a_(a), dBHAccretionAlpha(dBHAccretionAlpha), dBHAccretionCvisc(dBHAccretionCvisc),
          dBHAccretionEddFac(dBHAccretionEddFac), dConstGamma(dConstGamma), dBHFBEff(dBHFBEff), dBHFBEcrit(dBHFBEcrit), bBHAccretion(bBHAccretion),
          bBHFeedback(bBHFeedback) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
    auto range_neighbour(treeStore::NodePointer node) -> std::pair<int, int> override { return SmoothOnlyGas::range_neighbour(node); }
    auto accept_neighbour(particleStore::ParticleReference& p) -> bool override { return SmoothOnlyGas::accept_neighbour(p); }
};

void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference& p, double dTime, double dDelta, double dBHRadiativeEff);
