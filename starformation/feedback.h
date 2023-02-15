#include "pkd.h"
#include "smooth/smooth.h"
void smSNFeedback(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) ;
void initSNFeedback(void *vpkd, void *vp);
void combSNFeedback(void *vpkd, void *v1, const void *v2);

static inline void pkdAddFBEnergy(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph, double dConstGamma) {
#ifndef OLD_FB_SCHEME
    psph->Uint += psph->fAccFBEnergy;
    psph->E += psph->fAccFBEnergy;
#ifdef ENTROPY_SWITCH
    psph->S += psph->fAccFBEnergy*(dConstGamma-1.) *
               pow(pkdDensity(pkd,p), -dConstGamma+1);
#endif
    psph->fAccFBEnergy = 0.0;
#endif //OLD_FB_SCHEME
}

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
