#include "smooth/smooth.h"
#include "pkd.h"

void smBHevolve(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void combBHevolve(void *vpkd, void *p1,const void *p2);
void initBHevolve(void *vpkd,void *vp);

void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference &p, double dTime, double dDelta, double dBHRadiativeEff);
