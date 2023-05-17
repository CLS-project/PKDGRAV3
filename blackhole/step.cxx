#ifdef  BLACKHOLES
#include "blackhole/step.h"
#include "smooth/smooth.h"
#include "master.h"


void MSR::BHStep(double dTime, double dDelta) {
#ifndef DEBUG_BH_ONLY
    Smooth(dTime, dDelta,SMX_BH_STEP,0,param.nSmooth);
#endif
}


void smBHstep(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {

#ifndef DEBUG_BH_ONLY
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    uint8_t uMaxRung = 0;
    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        uMaxRung = (q.rung() > uMaxRung) ? q.rung() : uMaxRung;
    }
    p.set_new_rung(uMaxRung);
#endif

}

#endif
