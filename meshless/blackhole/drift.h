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
    auto range_neighbour(treeStore::NodePointer node) -> std::pair<int, int> override { return SmoothOnlyGas::range_neighbour(node); }
    auto accept_neighbour(particleStore::ParticleReference& p) -> bool override { return SmoothOnlyGas::accept_neighbour(p); }
};

void pkdBHReposition(PKD pkd);
