#include "pkd.h"
#include "smooth/smooth.h"


struct snFeedbackPack {
    blitz::TinyVector<double,3> position;
    float fMass;
#ifdef ENTROPY_SWITCH
    float fDensity;
#endif
    uint8_t iClass;
};


struct snFeedbackFlush {
#ifdef OLD_FB_SCHEME
    double E;
    double Uint;
#ifdef ENTROPY_SWITCH
    double S;
#endif
#else // OLD_FB_SCHEME
    float fAccFBEnergy;
#endif
};


void smSNFeedback(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) ;
void packSNFeedback(void *vpkd,void *dst,const void *src);
void unpackSNFeedback(void *vpkd,void *dst,const void *src);
void initSNFeedback(void *vpkd,void *dst);
void flushSNFeedback(void *vpkd,void *dst,const void *src);
void combSNFeedback(void *vpkd,void *dst,const void *src);


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

