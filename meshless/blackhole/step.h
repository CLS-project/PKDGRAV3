#include "smooth/smooth.h"
#include "pkd.h"

struct bhStepPack {
    blitz::TinyVector<double, 3> position;
    uint8_t uRung;
    uint8_t uNewRung;
    uint8_t iClass;
};

void smBHstep(particleStore::ParticleReference& p, float fBall, int nSmooth, NN* nnList, SMF* smf);
class SmoothHelperBHstep : public SmoothHelperPacked<bhStepPack> {
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;

  public:
    SmoothHelperBHstep(pkdContext* pkd) : SmoothHelperPacked<bhStepPack>(pkd, false) {}
};
