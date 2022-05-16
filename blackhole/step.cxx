#ifdef  BLACKHOLES
#include "blackhole/step.h"
#include "smooth/smooth.h"
#include "master.h"


void MSR::BHStep(double dTime, double dDelta) {
#ifndef DEBUG_BH_ONLY
    Smooth(dTime, dDelta,SMX_BH_STEP,0,param.nSmooth);
#endif
}


#ifdef __cplusplus
extern "C" {
#endif

void smBHstep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

#ifndef DEBUG_BH_ONLY
    PKD pkd = smf->pkd;
    uint8_t uMaxRung = 0;
    for (int i=0; i<nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;
        uMaxRung = (q->uRung > uMaxRung) ? q->uRung : uMaxRung;
    }

    p->uNewRung = uMaxRung;
#endif

}

#ifdef __cplusplus
}
#endif

#endif
