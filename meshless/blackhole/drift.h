#include "smooth/smooth.h"
#include "pkd.h"

struct bhGasPinPack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    float fMass;
    float fPotential;
    uint8_t iClass;
};

class SmoothHelperBHGasPin : public SmoothHelperPacked<bhGasPinPack> {
    virtual void    pack(void *dst, const void *src) override final;
    virtual void  unpack(void *dst, const void *src, const void *key) override final;
public:
    SmoothHelperBHGasPin(pkdContext *pkd) : SmoothHelperPacked<bhGasPinPack>(pkd,false) {}
};

void smBHGasPin(particleStore::ParticleReference &p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void pkdBHReposition(PKD pkd);

