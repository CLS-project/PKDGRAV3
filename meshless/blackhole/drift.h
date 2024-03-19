#include "smooth/smooth.h"
#include "pkd.h"

struct bhGasPinPack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    float fMass;
    float fPotential;
    uint8_t iClass;
};

void smBHGasPin(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packBHGasPin(void *vpkd,void *dst,const void *src);
void unpackBHGasPin(void *vpkd,void *dst,const void *src);

void pkdBHReposition(PKD pkd);

