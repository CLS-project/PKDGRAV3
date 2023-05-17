#ifdef FEEDBACK
#include <numeric>
#include "starformation/feedback.h"
#include "master.h"
#include "imf.h"

void MSR::SetFeedbackParam() {
    param.dSNFBDu = param.dSNFBDT * dTuFacPrimIonised;

    auto IMF = ChooseIMF(param.achIMFType, param.dIMFMinMass, param.dIMFMaxMass);
    const double dCCSNNumPerMass =
        IMF->UnweightedIntegration(param.dCCSNMinMass, param.dCCSNMaxMass);
    param.dCCSNFBSpecEnergy = (param.dCCSNEnergy / MSOLG) * dCCSNNumPerMass /
                              param.units.dErgPerGmUnit;

    param.dSNIaFBSpecEnergy = (param.dSNIaEnergy / MSOLG) * param.dSNIaNumPerMass /
                              param.units.dErgPerGmUnit;

    if (!param.bRestart) {
        param.dCCSNFBDelay *= SECONDSPERYEAR / param.units.dSecUnit;
        param.dSNIaFBDelay *= SECONDSPERYEAR / param.units.dSecUnit;

        const double dnHToRho = MHYDR / param.dInitialH / param.units.dGmPerCcUnit;
        param.dSNFBEffnH0 *= dnHToRho;
    }
}


static inline void snFeedback(PKD pkd, PARTICLE *pIn, int nSmooth, NN *nnList,
                              const float fNgbTotMass, const double dDeltau,
                              const double dStarMass, const double dSpecEnergy,
                              const float fEfficiency, const double dConstGamma);


/* Function that will be called with the information of all the neighbors.
 * The helper snFeedback computes the probability of explosion and adds
 * the energy to the gas particles in nnList
 */
void smSNFeedback(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];

    const auto fNgbTotMass = std::accumulate(nnList, nnList + nSmooth, 0.0f,
    [pkd](auto &a, auto &b) { return a + pkd->particles[b.pPart].mass(); }) - p.mass();

    auto &star = p.star();

#ifdef STELLAR_EVOLUTION
    const double dStarMass = star.fInitialMass;
#else
    const double dStarMass = p.mass();
#endif

    if (!star.bCCSNFBDone && ((smf->dTime - star.fTimer) > smf->dCCSNFBDelay)) {
        snFeedback(pkd, pIn, nSmooth, nnList, fNgbTotMass, smf->dSNFBDu, dStarMass,
                   smf->dCCSNFBSpecEnergy, star.fSNEfficiency, smf->dConstGamma);
        star.bCCSNFBDone = 1;
    }

    if (!star.bSNIaFBDone && ((smf->dTime - star.fTimer) > smf->dSNIaFBDelay)) {
        snFeedback(pkd, pIn, nSmooth, nnList, fNgbTotMass, smf->dSNFBDu, dStarMass,
                   smf->dSNIaFBSpecEnergy, star.fSNEfficiency, smf->dConstGamma);
        star.bSNIaFBDone = 1;
    }
}



void initSNFeedback(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(vp)];

    if (p.is_gas()) {
        auto &sph = p.sph();

#ifdef OLD_FB_SCHEME
        sph.Uint = 0.;
        sph.E = 0.;
#ifdef ENTROPY_SWITCH
        sph.S = 0.;
#endif
#else // OLD_FB_SCHEME
        sph.fAccFBEnergy = 0.;
#endif
    }

}


void combSNFeedback(void *vpkd, void *v1, const void *v2) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(v1)];
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(v2)];

    if (p1.is_gas() && p2.is_gas()) {
        auto &sph1 = p1.sph();
        const auto &sph2 = p2.sph();

#ifdef OLD_FB_SCHEME
        sph1.Uint += sph2.Uint;
        sph1.E += sph2.E;
#ifdef ENTROPY_SWITCH
        sph1.S += sph2.S;
#endif
#else //OLD_FB_SCHEME
        sph1.fAccFBEnergy += sph2.fAccFBEnergy;
#endif

    }

}


static inline void snFeedback(PKD pkd, PARTICLE *pIn, int nSmooth, NN *nnList,
                              const float fNgbTotMass, const double dDeltau,
                              const double dStarMass, const double dSpecEnergy,
                              const float fEfficiency, const double dConstGamma) {

    const double dProb = fEfficiency * dSpecEnergy * dStarMass /
                         (dDeltau * fNgbTotMass);
    const double dCorrFac = dProb > 1.0 ? dProb : 1.0;

    for (auto i = 0; i < nSmooth; ++i) {
        if (nnList[i].pPart == pIn) continue;

        if (rand() < RAND_MAX * dProb) {
            auto q = pkd->particles[nnList[i].pPart];
            auto &qsph = q.sph();
            const double dEnergyInput = dCorrFac * dDeltau * q.mass();

#ifdef OLD_FB_SCHEME
            qsph.Uint += dEnergyInput;
            qsph.E += dEnergyInput;
#ifdef ENTROPY_SWITCH
            qsph.S += dEnergyInput * (dConstGamma - 1.0) *
                      pow(q.density(), 1.0 - dConstGamma);
#endif
#else // OLD_FB_SCHEME
            qsph.fAccFBEnergy += dEnergyInput;
#endif
        }
    }
}

#endif // FEEDBACK
