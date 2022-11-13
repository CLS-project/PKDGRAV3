#include "pkd.h"
#include "smooth/smooth.h"
#ifdef __cplusplus
extern "C" {
#endif
void smSNFeedback(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) ;
void initSNFeedback(void *vpkd, void *vp);
void combSNFeedback(void *vpkd, void *v1, const void *v2);

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

#ifdef __cplusplus
}
#endif
