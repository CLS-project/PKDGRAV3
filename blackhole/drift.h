#include "smooth.h"
#include "pkd.h"
#include "smooth.h"
#include "master.h"

void smBHdrift(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void combBHdrift(void *vpkd, void *p1,void *p2);
void initBHdrift(void *vpkd,void *vp);

