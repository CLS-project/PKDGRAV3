#ifdef FEEDBACK
#include "starformation/feedback.h"
#include "master.h"

void MSR::SetFeedbackParam() {
    param.dSNFBDu = param.dSNFBDT * dTuFac * (param.dMeanMolWeight / 0.5917);
    param.dCCSNFBSpecEnergy = (param.dCCSNEnergy / MSOLG) * param.dCCSNFBNumPerMass /
                              param.units.dErgPerGmUnit;
    param.dSNIaFBSpecEnergy = (param.dSNIaEnergy / MSOLG) * param.dSNIaFBNumPerMass /
                              param.units.dErgPerGmUnit;

    if (!param.bRestart) {
        param.dCCSNFBDelay *= SECONDSPERYEAR / param.units.dSecUnit;
        param.dSNIaFBDelay *= SECONDSPERYEAR / param.units.dSecUnit;

        const double dnHToRho = MHYDR / param.dInitialH / param.units.dGmPerCcUnit;
        param.dSNFBEffnH0 *= dnHToRho;
    }
}

#ifdef __cplusplus
extern "C" {
#endif

static inline void snFeedback(PKD pkd, PARTICLE *p, int nSmooth, NN *nnList,
                              const float fNgbTotMass, const double dDeltau,
                              const double dStarMass, const double dSpecEnergy,
                              const float fEfficiency, const double dConstGamma);


/* Function that will be called with the information of all the neighbors.
 * The helper snFeedback computes the probability of explosion and adds
 * the energy to the gas particles in nnList
 */
void smSNFeedback(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;

    float fNgbTotMass = 0.0f;
    for (int i = 0; i < nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;
        fNgbTotMass += pkdMass(pkd, q);
    }
    fNgbTotMass -= pkdMass(pkd, p);

    STARFIELDS *pStar = pkdStar(pkd, p);

#ifdef STELLAR_EVOLUTION
    const double dStarMass = pStar->fInitialMass;
#else
    const double dStarMass = pkdMass(pkd, p);
#endif

    if (!pStar->bCCSNFBDone && ((smf->dTime - pStar->fTimer) > smf->dCCSNFBDelay)) {
        snFeedback(pkd, p, nSmooth, nnList, fNgbTotMass, smf->dSNFBDu, dStarMass,
                   smf->dCCSNFBSpecEnergy, pStar->fSNEfficiency, smf->dConstGamma);
        pStar->bCCSNFBDone = 1;
    }

    if (!pStar->bSNIaFBDone && ((smf->dTime - pStar->fTimer) > smf->dSNIaFBDelay)) {
        snFeedback(pkd, p, nSmooth, nnList, fNgbTotMass, smf->dSNFBDu, dStarMass,
                   smf->dSNIaFBSpecEnergy, pStar->fSNEfficiency, smf->dConstGamma);
        pStar->bSNIaFBDone = 1;
    }
}



void initSNFeedback(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = (PARTICLE *) vp;

    if (pkdIsGas(pkd,p)) {
        SPHFIELDS *psph = pkdSph(pkd,p);

#ifdef OLD_FB_SCHEME
        psph->Uint = 0.;
        psph->E = 0.;
#ifdef ENTROPY_SWITCH
        psph->S = 0.;
#endif
#else // OLD_FB_SCHEME
        psph->fAccFBEnergy = 0.;
#endif
    }

}


void combSNFeedback(void *vpkd, void *v1, const void *v2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) v1;
    PARTICLE *p2 = (PARTICLE *) v2;

    if (pkdIsGas(pkd,p1) && pkdIsGas(pkd,p2)) {

        SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);

#ifdef OLD_FB_SCHEME
        psph1->Uint += psph2->Uint;
        psph1->E += psph2->E;
#ifdef ENTROPY_SWITCH
        psph1->S += psph2->S;
#endif
#else //OLD_FB_SCHEME
        psph1->fAccFBEnergy += psph2->fAccFBEnergy;
#endif

    }

}


static inline void snFeedback(PKD pkd, PARTICLE *p, int nSmooth, NN *nnList,
                              const float fNgbTotMass, const double dDeltau,
                              const double dStarMass, const double dSpecEnergy,
                              const float fEfficiency, const double dConstGamma) {

    const double dProb = fEfficiency * dSpecEnergy * dStarMass /
                         (dDeltau * fNgbTotMass);
    assert(dProb < 1.0);

    for (int i = 0; i < nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;
        if (q == p) continue;

        if (rand() < RAND_MAX * dProb) {
            SPHFIELDS *qSph = pkdSph(pkd, q);
            const double dEnergyInput = dDeltau * pkdMass(pkd, q);

#ifdef OLD_FB_SCHEME
            qSph->Uint += dEnergyInput;
            qSph->E += dEnergyInput;
#ifdef ENTROPY_SWITCH
            qSph->S += dEnergyInput * (dConstGamma - 1.0) *
                       pow(pkdDensity(pkd, q), 1.0 - dConstGamma);
#endif
#else // OLD_FB_SCHEME
            qSph->fAccFBEnergy += dEnergyInput;
#endif
        }
    }
}

#ifdef __cplusplus
}
#endif
#endif // FEEDBACK
