#ifdef  BLACKHOLES
#include "blackhole/step.h"
#include "smooth/smooth.h"
#include "master.h"


void MSR::BHStep(double dTime, double dDelta) {
#ifndef DEBUG_BH_ONLY
    Smooth(dTime, dDelta,SMX_BH_STEP,0,parameters.get_nSmooth());
#endif
}


void packBHstep(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhStepPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->uRung = p2.rung();
    }
}

void unpackBHstep(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhStepPack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.set_rung(p2->uRung);
    }
}


void smBHstep(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
#ifndef DEBUG_BH_ONLY
    PKD pkd = smf->pkd;
    auto ii = std::max_element(nnList, nnList+nSmooth,
    [pkd](const auto &a,const auto &b) {
        auto p = pkd->particles[a.pPart];
        auto q = pkd->particles[b.pPart];
        return p.rung() < q.rung();
    });
    auto p = pkd->particles[pIn];
    p.set_new_rung(pkd->particles[ii->pPart].rung());
#endif
}

#endif
