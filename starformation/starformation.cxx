#ifdef  STAR_FORMATION
#include "starformation/starformation.h"
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#ifdef FEEDBACK
    #include "starformation/feedback.h"
#endif
#include "eEOS/eEOS.h"
#include "master.h"

/*
 * ---------------------
 * HELPER FUNCTIONS
 * ---------------------
 */

/* MSR layer
 */
void MSR::SetStarFormationParam() {
    param.dSFThresholdu = param.dSFThresholdT*dTuFac;

    if (csm->val.bComove) {
        assert(csm->val.dOmegab > 0.);
        // If in PKDGRAV3 units, this should always be unity
        const double rhoCrit0 = 3. * csm->val.dHubble0 * csm->val.dHubble0 /
                                (8. * M_PI);
        param.dSFThresholdOD = rhoCrit0 * csm->val.dOmegab *
                               param.dSFMinOverDensity;
    }
    else {
        param.dSFThresholdOD = 0.0;
    }

    if (!param.bRestart) {
        const double dnHToRho = MHYDR / param.dInitialH / param.units.dGmPerCcUnit;
        param.dSFThresholdDen *= dnHToRho; // Code physical density now

        const double Msolpcm2 = 1. / param.units.dMsolUnit *
                                pow(param.units.dKpcUnit*1e3, 2);
        param.dSFnormalizationKS *= 1. / param.units.dMsolUnit *
                                    param.units.dSecUnit/SECONDSPERYEAR *
                                    pow(param.units.dKpcUnit, 2) *
                                    pow(Msolpcm2,-param.dSFindexKS);
    }
}

void MSR::StarForm(double dTime, double dDelta, int iRung) {
    struct inStarForm in;
    struct outStarForm out;
    double sec,sec1,dsec;

    TimerStart(TIMER_STARFORM);

    const double a = csmTime2Exp(csm,dTime);

    in.dScaleFactor = a;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.dSFindexKS = param.dSFindexKS;
    in.dSFnormalizationKS = param.dSFnormalizationKS;
    in.dConstGamma = param.dConstGamma;
    in.dSFGasFraction = param.dSFGasFraction;
    in.dSFThresholdDen = param.dSFThresholdDen * a * a * a;  // Send in comoving units
    in.dSFThresholdOD = param.dSFThresholdOD;
    in.dSFThresholdu = param.dSFThresholdu;
    in.dSFEfficiency = param.dSFEfficiency;
#ifdef HAVE_METALLICITY
    in.bSFThresholdDenSchaye2004 = param.bSFThresholdDenSchaye2004;
#endif
#ifdef FEEDBACK
    in.bCCSNFeedback = param.bCCSNFeedback;
    in.bSNIaFeedback = param.bSNIaFeedback;
    in.dSNFBEfficiency = param.dSNFBEfficiency;
    in.dSNFBMaxEff = param.dSNFBMaxEff;
    in.dSNFBEffIndex = param.dSNFBEffIndex;
    in.dSNFBEffnH0 = param.dSNFBEffnH0;
#endif
#ifdef EEOS_POLYTROPE
    // In the rare case that not EEOS_POLYTROPE, this will be unused
    in.dEOSPolyFloorIndex = param.dEOSPolyFloorIndex;
    in.dEOSPolyFloorDen = param.dEOSPolyFloorDen;
    in.dEOSPolyFlooru = param.dEOSPolyFlooru;
#endif
#ifdef STELLAR_EVOLUTION
    in.dSNIaMaxMass = param.dSNIaMaxMass;
    in.dCCSNMinMass = param.dCCSNMinMass;
    in.dCCSNMaxMass = param.dCCSNMaxMass;
#endif

    if (param.bVDetails) printf("Star Form (rung: %d) ... ", iRung);

    ActiveRung(iRung,1);
    pstStarForm(pst, &in, sizeof(in), &out, 0);
    if (param.bVDetails)
        printf("%d Stars formed with mass %g, %d gas deleted\n",
               out.nFormed, out.dMassFormed, out.nDeleted);
    massFormed += out.dMassFormed;
    starFormed += out.nFormed;

    nGas -= out.nFormed;
    nStar += out.nFormed;

    TimerStop(TIMER_STARFORM);
    dsec = TimerGet(TIMER_STARFORM);
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);

}

#ifdef __cplusplus
extern "C" {
#endif

static inline double pressure_SFR(const float fMass, const float fDens, const double a_m3,
                                  const double dThreshDen, const double dSFnormalizationKS,
                                  const double dSFindexKS, const double dSFGasFraction,
                                  const double dEOSPolyFloorIndex, const double dEOSPolyFloorDen,
                                  const double dEOSPolyFlooru, const double dConstGamma,
                                  SPHFIELDS *psph);

static inline double density_SFR(const float fMass, const float fDens, const double a_m3,
                                 const double dThreshDen, const double dSFThresholdu,
                                 const double dSFEfficiency, SPHFIELDS *psph);


int pstStarForm(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inStarForm *in = (struct inStarForm *) vin;
    struct outStarForm *out = (struct outStarForm *) vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inStarForm));
    if (pst->nLeaves > 1) {
        struct outStarForm fsStats;

        rID = mdlReqService(pst->mdl,pst->idUpper,PST_STARFORM,in,nIn);
        pstStarForm(pst->pstLower,in,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&fsStats,NULL);
        out->nFormed += fsStats.nFormed;
        out->nDeleted += fsStats.nDeleted;
        out->dMassFormed += fsStats.dMassFormed;
    }
    else {
        pkdStarForm(pst->plcl->pkd, *in,
                    &out->nFormed, &out->dMassFormed, &out->nDeleted);
    }
    return sizeof(struct outStarForm);
}




void pkdStarForm(PKD pkd,
                 struct inStarForm in,
                 int *nFormed, /* number of stars formed */
                 double *dMassFormed,   /* mass of stars formed */
                 int *nDeleted) { /* gas particles deleted */

    PARTICLE *p;
    SPHFIELDS *psph;
    float *pv;
    int i, j;

    assert(pkd->oFieldOffset[oStar]);
    assert(pkd->oFieldOffset[oSph]);
    assert(pkd->oFieldOffset[oMass]);

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;

    const double a_m1 = 1.0/in.dScaleFactor;
    const double a_m2 = a_m1*a_m1;
    const double a_m3 = a_m2*a_m1;

    for (i=0; i<pkdLocal(pkd); ++i) {
        p = pkdParticle(pkd,i);

        if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
            const float fDens = pkdDensity(pkd,p);
            if (fDens < in.dSFThresholdOD) continue;

            const float fMass = pkdMass(pkd, p);
            psph = pkdSph(pkd,p);

#ifdef HAVE_METALLICITY
            double dThreshDen;
            if (in.bSFThresholdDenSchaye2004) {
                const double dMetAbun = (double)psph->fMetalMass / (double)fMass;
                const double dThreshFac = dMetAbun < 1.5e-6 ? 100.0 : pow(2e-3 / dMetAbun, 0.64);
                dThreshDen = in.dSFThresholdDen * dThreshFac;
            }
            else {
                dThreshDen = in.dSFThresholdDen;
            }
#else
            const double dThreshDen = in.dSFThresholdDen;
#endif

            double dmstar;
            if (in.dSFEfficiency > 0.0) {
                dmstar = density_SFR(fMass, fDens, a_m3, dThreshDen, in.dSFThresholdu,
                                     in.dSFEfficiency, psph);
            }
            else {
                dmstar = pressure_SFR(fMass, fDens, a_m3, dThreshDen, in.dSFnormalizationKS,
                                      in.dSFindexKS, in.dSFGasFraction, in.dEOSPolyFloorIndex,
                                      in.dEOSPolyFloorDen, in.dEOSPolyFlooru, in.dConstGamma,
                                      psph);
            }

            psph->SFR = dmstar;
            if (dmstar == 0.0)
                continue;

            const double dt = in.dTime - psph->lastUpdateTime;
            const double prob = 1.0 - exp(-dmstar * dt / fMass);

            // Star formation event?
            if (rand()<RAND_MAX*prob) {
#ifdef STELLAR_EVOLUTION
                float afElemMass[ELEMENT_COUNT];
                float fMetalMass;
                for (j = 0; j < ELEMENT_COUNT; j++)
                    afElemMass[j] = psph->afElemMass[j];
                fMetalMass = psph->fMetalMass;
#endif

                // We just change the class of the particle to stellar one
                pkdSetClass(pkd, fMass, pkdSoft0(pkd,p), FIO_SPECIES_STAR, p);
                STARFIELDS *pStar = pkdStar(pkd, p);

                // When changing the class, we have to take into account that
                // the code velocity has different scale factor dependencies for
                // dm/star particles and gas particles
                pv = pkdVel(pkd,p);
                for (j=0; j<3; j++) {
                    pv[j] *= in.dScaleFactor;
                }

                // We log statistics about the formation time
                pStar->fTimer = in.dTime;

#ifdef STELLAR_EVOLUTION
                const float fMassInv = 1.0f / fMass;
                for (j = 0; j < ELEMENT_COUNT; j++)
                    pStar->afElemAbun[j] = afElemMass[j] * fMassInv;
                pStar->fMetalAbun = fMetalMass * fMassInv;

                pStar->fInitialMass = fMass;
                pStar->fLastEnrichTime = 0.0f;

                stevStarParticleInit(pkd, pStar, in.dSNIaMaxMass, in.dCCSNMinMass,
                                     in.dCCSNMaxMass);
#endif

#ifdef FEEDBACK
                pStar->bCCSNFBDone = in.bCCSNFeedback ? 0 : 1;
                pStar->bSNIaFBDone = in.bSNIaFeedback ? 0 : 1;

                // Compute the feedback efficiency for this particle based
                // on the birth information (eq. 7 Schaye 2015)
#ifdef HAVE_METALLICITY
                const double fMetalAbun = pStar->fMetalAbun;
#else
                const double fMetalAbun = 0.00127; // 0.1 Zsolar
#endif
                pStar->fSNEfficiency = SNFeedbackEfficiency(fMetalAbun,
                                       pkdDensity(pkd,p)*a_m3, in.dSNFBEfficiency,
                                       in.dSNFBMaxEff, in.dSNFBEffIndex, in.dSNFBEffnH0);
#endif

                // Safety check
                assert(pkdIsStar(pkd,p));
                assert(!pkdIsGas(pkd,p));

                (*nFormed)++;
                *dMassFormed += fMass;

                pkd->nGas -= 1;
                pkd->nStar += 1;
            }
        }
    }
}
#endif

static inline double pressure_SFR(const float fMass, const float fDens, const double a_m3,
                                  const double dThreshDen, const double dSFnormalizationKS,
                                  const double dSFindexKS, const double dSFGasFraction,
                                  const double dEOSPolyFloorIndex, const double dEOSPolyFloorDen,
                                  const double dEOSPolyFlooru, const double dConstGamma,
                                  SPHFIELDS *psph) {
    // Two SF thresholds are applied:
    //      a) Minimum density, computed at the master level
    //      b) Maximum temperature of a factor 0.5 dex (i.e., 3.1622)
    //         above the polytropic eEOS
#ifdef EEOS_POLYTROPE
    const double maxUint = 3.16228 * fMass *
                           polytropicEnergyFloor(a_m3, fDens, dEOSPolyFloorIndex,
                                   dEOSPolyFloorDen, dEOSPolyFlooru);
#else
    const double maxUint = INFINITY;
#endif

    if (psph->Uint > maxUint || fDens < dThreshDen) {
        return 0.0;
    }

    const double SFexp = 0.5*(dSFindexKS-1.);
    const double dmstar = dSFnormalizationKS * fMass *
                          pow(dConstGamma * dSFGasFraction * psph->P * a_m3, SFexp);
    return dmstar;
}

static inline double density_SFR(const float fMass, const float fDens, const double a_m3,
                                 const double dThreshDen, const double dSFThresholdu,
                                 const double dSFEfficiency, SPHFIELDS *psph) {
    // Two SF thresholds are applied:
    //      a) Minimum density, computed at the master level
    //      b) Maximum temperature set by dSFThresholdTemp
    if (psph->Uint > dSFThresholdu || fDens < dThreshDen) {
        return 0.0;
    }

    const double tff = sqrt(3.*M_PI/(32.*fDens*a_m3));
    const double dmstar = dSFEfficiency * fMass / tff;
    return dmstar;
}

#ifdef __cplusplus
}
#endif
