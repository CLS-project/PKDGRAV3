#include "smooth/smooth.h"
#include "pkd.h"

#define NOT_ACCRETED -1

struct bhEvolvePack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    float c;
    float fMass;
#ifdef ENTROPY_SWITCH
    float fDensity;
#endif
    uint8_t iClass;
};

struct bhEvolveFlush {
#ifdef OLD_FB_SCHEME
    double E;
    double Uint;
#ifdef ENTROPY_SWITCH
    double S;
#endif
#else // OLD_FB_SCHEME
    float fAccFBEnergy;
#endif
    int iAccPid;
    int iAccIndex;
};

void smBHEvolve(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packBHEvolve(void *vpkd,void *dst,const void *src);
void unpackBHEvolve(void *vpkd,void *dst,const void *src);
void initBHEvolve(void *vpkd,void *dst);
void flushBHEvolve(void *vpkd,void *dst,const void *src);
void combBHEvolve(void *vpkd,void *dst,const void *src);

void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference &p, double dTime, double dDelta, double dBHRadiativeEff);

