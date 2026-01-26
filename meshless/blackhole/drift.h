#include "smooth/smooth.h"
#include "pkd.h"

struct bhGasPinPack {
    blitz::TinyVector<double, 3> position;
    blitz::TinyVector<double, 3> velocity;
    float fMass;
    float fPotential;
    uint8_t iClass;
};

class SmoothHelperBHGasPin : public SmoothHelperPacked<bhGasPinPack> {
    double a_;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;

  public:
    SmoothHelperBHGasPin(pkdContext* pkd, double a) : SmoothHelperPacked<bhGasPinPack>(pkd, false), a_(a) {}
    auto gasOnly() -> bool override { return true; }
};

void pkdBHReposition(PKD pkd);
