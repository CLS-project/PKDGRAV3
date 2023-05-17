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
    param.dSFThresholdu = param.dSFThresholdT*dTuFacPrimNeutral;

    if (csm->val.bComove) {
        // As usual, in code units the critical density is 1.0
        param.dSFThresholdOD = param.dSFMinOverDensity * csm->val.dOmegab;
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

int MSR::ValidateStarFormationParam() {
#if !defined(EEOS_POLYTROPE) && !defined(EEOS_JEANS)
    fprintf(stderr,"WARNING: Star formation is active but no eEOS is selected!\n");
#endif
    return 1;
}

void MSR::StarForm(double dTime, double dDelta, int iRung) {
    struct inStarForm in;
    struct outStarForm out;

    TimerStart(TIMER_STARFORM);

    const double a = csmTime2Exp(csm,dTime);

    in.dHubble = Comove() ? csmTime2Hub(csm,dTime) : 0.0;
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
#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
    eEOSFill(param, &in.eEOS);
#endif
#ifdef STELLAR_EVOLUTION
    in.dSNIaMaxMass = param.dSNIaMaxMass;
    in.dCCSNMinMass = param.dCCSNMinMass;
    in.dCCSNMaxMass = param.dCCSNMaxMass;
#endif

    if (parameters.get_bVDetails()) printf("Star Form (rung: %d) ... ", iRung);

    ActiveRung(iRung,1);
    pstStarForm(pst, &in, sizeof(in), &out, 0);
    if (parameters.get_bVDetails())
        printf("%d Stars formed with mass %g, %d gas deleted\n",
               out.nFormed, out.dMassFormed, out.nDeleted);
    massFormed += out.dMassFormed;
    starFormed += out.nFormed;

    nGas -= out.nFormed;
    nStar += out.nFormed;

    TimerStop(TIMER_STARFORM);
    auto dsec = TimerGet(TIMER_STARFORM);
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);

}

static inline double pressure_SFR(const float fMass, const float fDens,
                                  const double fBall, const double a_m3,
                                  const double dThreshDen, const double dSFnormalizationKS,
                                  const double dSFindexKS, const double dSFGasFraction,
                                  const double dConstGamma, eEOSparam eEOS,
                                  SPHFIELDS &sph);

static inline double density_SFR(const float fMass, const float fDens, const double dHubble,
                                 const double a_m1, const double a_m3, const double dThreshDen,
                                 const double dSFThresholdu, const double dSFEfficiency,
                                 SPHFIELDS &sph);


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

    assert(pkd->particles.present(PKD_FIELD::oStar));
    assert(pkd->particles.present(PKD_FIELD::oSph));
    assert(pkd->particles.present(PKD_FIELD::oMass));

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;

    const double a_m1 = 1.0/in.dScaleFactor;
    const double a_m2 = a_m1*a_m1;
    const double a_m3 = a_m2*a_m1;

    for (auto &p : pkd->particles) {

        if (p.is_active() && p.is_gas()) {
            const auto &mass = p.mass();
            auto &sph = p.sph();

#ifdef HAVE_METALLICITY
            double dThreshDen;
            if (in.bSFThresholdDenSchaye2004) {
                const double dMetAbun = (double)sph.fMetalMass / (double)mass;
                const double dThreshFac = dMetAbun < 1.5e-6 ? 100.0 : pow(2e-3 / dMetAbun, 0.64);
                dThreshDen = in.dSFThresholdDen * dThreshFac;
            }
            else {
                dThreshDen = in.dSFThresholdDen;
            }
#else
            double dThreshDen = in.dSFThresholdDen;
#endif

            dThreshDen = dThreshDen > in.dSFThresholdOD ? dThreshDen : in.dSFThresholdOD;

            double dmstar;
            if (in.dSFEfficiency > 0.0) {
                dmstar = density_SFR(mass, p.density(), in.dHubble, a_m1, a_m3,
                                     dThreshDen, in.dSFThresholdu, in.dSFEfficiency, sph);
            }
            else {
                dmstar = pressure_SFR(mass, p.density(), p.ball(), a_m3,
                                      dThreshDen, in.dSFnormalizationKS, in.dSFindexKS,
                                      in.dSFGasFraction, in.dConstGamma, in.eEOS, sph);
            }

            sph.SFR = dmstar;
            if (dmstar == 0.0)
                continue;

            const double dt = in.dTime - sph.lastUpdateTime;
            const double prob = 1.0 - exp(-dmstar * dt / mass);

            // Star formation event?
            if (rand()<RAND_MAX*prob) {
#ifdef STELLAR_EVOLUTION
                float afElemMass[ELEMENT_COUNT];
                float fMetalMass;
                for (auto i = 0; i < ELEMENT_COUNT; ++i)
                    afElemMass[i] = sph.afElemMass[i];
                fMetalMass = sph.fMetalMass;
#endif

                // We just change the class of the particle to stellar one
                p.set_class(mass, p.soft0(),0, FIO_SPECIES_STAR);
                auto &star = p.star();

                // When changing the class, we have to take into account that
                // the code velocity has different scale factor dependencies for
                // dm/star particles and gas particles
                p.velocity() *= in.dScaleFactor;

                // We log statistics about the formation time
                star.fTimer = in.dTime;

#ifdef STELLAR_EVOLUTION
                const float fMassInv = 1.0f / mass;
                for (auto i = 0; i < ELEMENT_COUNT; ++i)
                    star.afElemAbun[i] = afElemMass[i] * fMassInv;
                star.fMetalAbun = fMetalMass * fMassInv;

                star.fInitialMass = mass;
                star.fLastEnrichTime = 0.0f;

                stevStarParticleInit(pkd, star, in.dSNIaMaxMass, in.dCCSNMinMass,
                                     in.dCCSNMaxMass);
#endif

#ifdef FEEDBACK
                star.bCCSNFBDone = in.bCCSNFeedback ? 0 : 1;
                star.bSNIaFBDone = in.bSNIaFeedback ? 0 : 1;

                // Compute the feedback efficiency for this particle based
                // on the birth information (eq. 7 Schaye 2015)
#ifdef HAVE_METALLICITY
                const double fMetalAbun = star.fMetalAbun;
#else
                const double fMetalAbun = 0.00127; // 0.1 Zsolar
#endif
                star.fSNEfficiency = SNFeedbackEfficiency(fMetalAbun,
                                     p.density()*a_m3, in.dSNFBEfficiency,
                                     in.dSNFBMaxEff, in.dSNFBEffIndex, in.dSNFBEffnH0);
#endif

                // Safety check
                assert(p.is_star());
                assert(!p.is_gas());

                (*nFormed)++;
                *dMassFormed += mass;

                pkd->nGas -= 1;
                pkd->nStar += 1;
            }
        }
    }
}
#endif

static inline double pressure_SFR(const float fMass, const float fDens,
                                  const double fBall, const double a_m3,
                                  const double dThreshDen, const double dSFnormalizationKS,
                                  const double dSFindexKS, const double dSFGasFraction,
                                  const double dConstGamma, eEOSparam eEOS,
                                  SPHFIELDS &sph) {
    // Two SF thresholds are applied:
    //      a) Minimum density, computed at the master level
    //      b) Maximum temperature of a factor 0.5 dex (i.e., 3.1622)
    //         above the polytropic eEOS
#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
    const double maxUint = 3.16228 * fMass *
                           eEOSEnergyFloor(a_m3, fDens, fBall, dConstGamma, eEOS);
#else
    // This configuration is allowed, but can easily produce numerical fragmentation!!!
    const double maxUint = INFINITY;
#endif

#if defined(FEEDBACK) || defined(BLACKHOLES)
#ifdef OLD_FB_SCHEME
    const double dUint = sph.Uint;
#else
    const double dUint = sph.Uint + sph.fAccFBEnergy;
#endif
#else
    const double dUint = sph.Uint;
#endif

    if (dUint > maxUint || fDens < dThreshDen) {
        return 0.0;
    }

    const double SFexp = 0.5*(dSFindexKS-1.);
    const double dmstar = dSFnormalizationKS * fMass *
                          pow(dConstGamma * dSFGasFraction * sph.P * a_m3, SFexp);
    return dmstar;
}

static inline double density_SFR(const float fMass, const float fDens, const double dHubble,
                                 const double a_m1, const double a_m3, const double dThreshDen,
                                 const double dSFThresholdu, const double dSFEfficiency,
                                 SPHFIELDS &sph) {
    // Three SF criteria are enforced:
    //      a) A minimum density, computed at the master level
    //      b) A maximum temperature, set by dSFThresholdTemp
    //      c) Converging flow (negative physical velocity divergence)
    const double dVelDiv = a_m1 * (sph.gradVx[0] + sph.gradVy[1] + sph.gradVz[2]) +
                           3. * dHubble;
    if (fDens < dThreshDen || sph.Uint > dSFThresholdu || dVelDiv > 0.0) {
        return 0.0;
    }

    const double tff = sqrt(3.*M_PI/(32.*fDens*a_m3));
    const double dmstar = dSFEfficiency * fMass / tff;
    return dmstar;
}

