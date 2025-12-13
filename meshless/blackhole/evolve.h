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

void smBHEvolve(particleStore::ParticleReference& p, float fBall, int nSmooth, NN* nnList, SMF* smf);
class SmoothHelperBHEvolve : public SmoothHelperPacked<bhEvolvePack, bhEvolveFlush> {
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;
    virtual void init(void* dst) override final;
    virtual void flush(void* dst, const void* src) override final;
    virtual void combine(void* dst, const void* src, const void* key) override final;
    virtual void init_shared(void* dst) override final { init_shared_bh_decorator(pkd, dst); }

  public:
    SmoothHelperBHEvolve(pkdContext* pkd) : SmoothHelperPacked<bhEvolvePack, bhEvolveFlush>(pkd, true) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};

void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference& p, double dTime, double dDelta, double dBHRadiativeEff);
