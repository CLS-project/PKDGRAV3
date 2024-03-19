#include "smooth/smooth.h"
#include "pkd.h"

struct bhStepPack {
    blitz::TinyVector<double,3> position;
    uint8_t uRung;
    uint8_t uNewRung;
    uint8_t iClass;
};

void smBHstep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packBHstep(void *vpkd,void *dst,const void *src);
void unpackBHstep(void *vpkd,void *dst,const void *src);

