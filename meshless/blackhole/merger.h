#include "smooth/smooth.h"
#include "pkd.h"

class SmoothHelperBHmerger : public SmoothHelper {
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    // TODO: allow for BH mergers across domains
    // virtual void combine(void *dst, const void *src, const void *key) override final;
  public:
    SmoothHelperBHmerger(pkdContext* pkd) : SmoothHelper(pkd, true) {}
};
uint64_t smReSmoothBHNode(SMX smx, SMF* smf, int iSmoothType);
