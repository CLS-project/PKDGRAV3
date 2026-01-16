#include "smooth/smooth.h"
#include "pkd.h"

struct bhStepPack {
    blitz::TinyVector<double, 3> position;
    uint8_t uRung;
    uint8_t uNewRung;
    uint8_t iClass;
};

class SmoothHelperBHstep : public SmoothHelperPacked<bhStepPack> {
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    void pack(void* dst, const void* src) final;
    void unpack(void* dst, const void* src, const void* key) final;
    auto shouldSmooth(particleStore::ParticleReference& p) -> bool final { return p.is_bh(); }

  public:
    SmoothHelperBHstep(pkdContext* pkd) : SmoothHelperPacked<bhStepPack>(pkd, false) {}
    auto gasOnly() -> bool override { return true; }
};
