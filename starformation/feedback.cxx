#ifdef FEEDBACK
#include <numeric>
#include "starformation/feedback.h"
#include "master.h"
#include "imf.h"

void MSR::SetFeedbackParam() {
    calc.dSNFBDu = parameters.get_dSNFBDT() * dTuFacPrimIonised;

    auto IMF = ChooseIMF(parameters.get_achIMFType().data(), parameters.get_dIMFMinMass(), parameters.get_dIMFMaxMass());
    const double dCCSNNumPerMass =
        IMF->UnweightedIntegration(parameters.get_dCCSNMinMass(), parameters.get_dCCSNMaxMass());
    calc.dCCSNFBSpecEnergy = (parameters.get_dCCSNEnergy() / MSOLG) * dCCSNNumPerMass /
                              units.dErgPerGmUnit;

    calc.dSNIaFBSpecEnergy = (parameters.get_dSNIaEnergy() / MSOLG) * parameters.get_dSNIaNumPerMass() /
                              units.dErgPerGmUnit;

    calc.dCCSNFBDelay = parameters.get_dCCSNFBDelay() * SECONDSPERYEAR / units.dSecUnit;
    calc.dSNIaFBDelay = parameters.get_dSNIaFBDelay() * SECONDSPERYEAR / units.dSecUnit;

    const double dnHToRho = MHYDR / parameters.get_dInitialH() / units.dGmPerCcUnit;
    calc.dSNFBEffnH0 = parameters.get_dSNFBEffnH0() * dnHToRho;
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


void packSNFeedback(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<snFeedbackPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->fMass = p2.mass();
#ifdef ENTROPY_SWITCH
        p1->fDensity = p2.density();
#endif
    }
}

void unpackSNFeedback(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const snFeedbackPack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.set_mass(p2->fMass);
#ifdef ENTROPY_SWITCH
        p1.set_density(p2->fDensity);
#endif
    }
}

void initSNFeedback(void *vpkd,void *dst) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(dst)];

    if (p.is_gas()) {
        auto &sph = p.sph();

#ifdef OLD_FB_SCHEME
        sph.E = 0.;
        sph.Uint = 0.;
#ifdef ENTROPY_SWITCH
        sph.S = 0.;
#endif
#else // OLD_FB_SCHEME
        sph.fAccFBEnergy = 0.;
#endif
    }
}

void flushSNFeedback(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<snFeedbackFlush *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    if (p2.is_gas()) {
        const auto &sph = p2.sph();

#ifdef OLD_FB_SCHEME
        p1->E = sph.E;
        p1->Uint = sph.Uint;
#ifdef ENTROPY_SWITCH
        p1->S = sph.S;
#endif
#else // OLD_FB_SCHEME
        p1->fAccFBEnergy = sph.fAccFBEnergy;
#endif
    }
}

void combSNFeedback(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const snFeedbackFlush *>(src);

    if (p1.is_gas()) {
        auto &sph = p1.sph();

#ifdef OLD_FB_SCHEME
        sph.E += p2->E;
        sph.Uint += p2->Uint;
#ifdef ENTROPY_SWITCH
        sph.S += p2->S;
#endif
#else // OLD_FB_SCHEME
        sph.fAccFBEnergy += p2->fAccFBEnergy;
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
