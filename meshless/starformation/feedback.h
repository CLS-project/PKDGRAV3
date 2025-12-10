#include "pkd.h"
#include "smooth/smooth.h"

static inline void pkdAddFBEnergy(PKD pkd,
                                  particleStore::ParticleReference &p,
                                  decorators::MFM *psph, double dConstGamma) {
    psph->Uint += psph->fAccFBEnergy;
    psph->E += psph->fAccFBEnergy;
#ifdef ENTROPY_SWITCH
    psph->S += psph->fAccFBEnergy*(dConstGamma-1.) *
               pow(p.density(), -dConstGamma+1);
#endif
    psph->fAccFBEnergy = 0.0;
}

struct snFeedbackPack {
    blitz::TinyVector<double,3> position;
    float fMass;
#ifdef ENTROPY_SWITCH
    float fDensity;
#endif
    uint8_t iClass;
};


struct snFeedbackFlush {
    float fAccFBEnergy;
};


void smSNFeedback(particleStore::ParticleReference &p,float fBall,int nSmooth,NN *nnList,SMF *smf) ;
class SmoothHelperSNFeedback : public SmoothHelperPacked<snFeedbackPack,snFeedbackFlush> {
    virtual void    pack(void *dst, const void *src) override final;
    virtual void  unpack(void *dst, const void *src, const void *key) override final;
    virtual void    init(void *dst) override final;
    virtual void   flush(void *dst, const void *src) override final;
    virtual void combine(void *dst, const void *src, const void *key) override final;
    virtual void  init_shared(void *dst) override final {init_shared_gas_decorator(pkd,dst);}
public:
    SmoothHelperSNFeedback(pkdContext *pkd) : SmoothHelperPacked<snFeedbackPack,snFeedbackFlush>(pkd,true) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};


static inline float SNFeedbackEfficiency(float Z, float rho, double dSNFBEfficiency,
        double dSNFBMaxEff, double dSNFBEffIndex,
        double dSNFBEffnH0) {
    if (dSNFBMaxEff > 0.0) {
        const double den = 1.0 + pow(Z/0.00127, dSNFBEffIndex) *
                           pow(rho/dSNFBEffnH0,-dSNFBEffIndex);
        return dSNFBEfficiency + (dSNFBMaxEff - dSNFBEfficiency)/den;
    }
    else {
        return dSNFBEfficiency;
    }
}

