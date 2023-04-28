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
    const double dHydFrac = param.dInitialH;
    const double dnHToRho = MHYDR / dHydFrac / param.units.dGmPerCcUnit;
    param.dSFThresholdu = param.dSFThresholdT*dTuFac;

    const double Msolpcm2 = 1. / param.units.dMsolUnit *
                            pow(param.units.dKpcUnit*1e3, 2);
    if (!param.bRestart) {
        param.dSFThresholdDen *= dnHToRho*dHydFrac; // Code hydrogen density
        param.dSFnormalizationKS *= 1. / param.units.dMsolUnit *
                                    param.units.dSecUnit/SECONDSPERYEAR *
                                    pow(param.units.dKpcUnit, 2) *
                                    pow(Msolpcm2,-param.dSFindexKS);
    }
}

void MSR::StarForm(double dTime, double dDelta, int iRung) {
    struct inStarForm in;
    struct outStarForm out;

    TimerStart(TIMER_STARFORM);

    const double a = csmTime2Exp(csm,dTime);

    in.dScaleFactor = a;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.dSFindexKS = param.dSFindexKS;
    in.dSFnormalizationKS = param.dSFnormalizationKS;
    in.dConstGamma = param.dConstGamma;
    in.dSFGasFraction = param.dSFGasFraction;
    in.dSFThresholdu = param.dSFThresholdu;
    in.dSFEfficiency = param.dSFEfficiency;
#ifdef FEEDBACK
    in.dSNFBEfficiency = param.dSNFBEfficiency;
    in.dSNFBMaxEff = param.dSNFBMaxEff;
    in.dSNFBEffnH0 = param.dSNFBEffnH0;
    in.dSNFBEffIndex = param.dSNFBEffIndex;
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


    // Here we set the minium density a particle must have to be SF
    //  NOTE: We still have to divide by the hydrogen fraction of each particle!
    if (csm->val.bComove) {
        double a3 = a*a*a;
        in.dDenMin = param.dSFThresholdDen*a3;
        assert(csm->val.dOmegab  > 0.);

        // If in PKDGRAV3 units, this should always be unity
        double rhoCrit0 = 3. * csm->val.dHubble0 * csm->val.dHubble0 /
                          (8. * M_PI);
        double denCosmoMin = rhoCrit0 * csm->val.dOmegab *
                             param.dSFMinOverDensity *
                             param.dInitialH;
        in.dDenMin = (in.dDenMin > denCosmoMin) ? in.dDenMin : denCosmoMin;
    }
    else {
        in.dDenMin = param.dSFThresholdDen;
    }


    if (param.get_bVDetails()) printf("Star Form (rung: %d) ... ", iRung);

    ActiveRung(iRung,1);
    pstStarForm(pst, &in, sizeof(in), &out, 0);
    if (param.get_bVDetails())
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

static inline double pressure_SFR(PKD pkd, double a_m3, double dDenMin,
                                  double dSFnormalizationKS, double dSFindexKS, double dSFGasFraction,
                                  double dEOSPolyFloorIndex, double dEOSPolyFloorDen, double dEOSPolyFlooru,
                                  double dConstGamma, particleStore::ParticleReference &p, SPHFIELDS &sph);
static inline double density_SFR(PKD pkd, double a_m3, double dDenMin,
                                 double dSFThresholdu, double dSFEfficiency,
                                 particleStore::ParticleReference &p, SPHFIELDS &sph);


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

    double dt;

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
            auto &psph = p.sph();
            dt = in.dTime - psph.lastUpdateTime;

            double dmstar;
            if (in.dSFEfficiency > 0.0) {
                dmstar = density_SFR(pkd, a_m3, in.dDenMin, in.dSFThresholdu, in.dSFEfficiency, p, psph);
            }
            else {
                dmstar = pressure_SFR(pkd, a_m3, in.dDenMin,
                                      in.dSFnormalizationKS, in.dSFindexKS, in.dSFGasFraction,
                                      in.dEOSPolyFloorIndex, in.dEOSPolyFloorDen, in.dEOSPolyFlooru,
                                      in.dConstGamma, p, psph);
            }

            psph.SFR = dmstar;
            if (dmstar == 0.0)
                continue;

            const double prob = 1.0 - exp(-dmstar*dt/p.mass());

            // Star formation event?
            if (rand()<RAND_MAX*prob) {
                float fMass = p.mass();

                //printf("STARFORM %e %e %e \n", dScaleFactor, rho_H, psph->Uint);

#ifdef STELLAR_EVOLUTION
                float afElemMass[ELEMENT_COUNT];
                float fMetalMass;
                for (auto j = 0; j < ELEMENT_COUNT; j++)
                    afElemMass[j] = psph.afElemMass[j];
                fMetalMass = psph.fMetalMass;
#endif

                // We just change the class of the particle to stellar one
                p.set_class(fMass, p.soft0(),0, FIO_SPECIES_STAR);
                auto &star = p.star();

                // When changing the class, we have to take into account that
                // the code velocity has different scale factor dependencies for
                // dm/star particles and gas particles
                p.velocity() *= in.dScaleFactor;

                // We log statistics about the formation time
                star.fTimer = in.dTime;
                star.hasExploded = 0;

#ifdef STELLAR_EVOLUTION
                const float fMassInv = 1.0f / fMass;
                for (auto j = 0; j < ELEMENT_COUNT; j++)
                    star.afElemAbun[j] = afElemMass[j] * fMassInv;
                star.fMetalAbun = fMetalMass * fMassInv;

                star.fInitialMass = fMass;
                star.fLastEnrichTime = 0.0f;

                stevStarParticleInit(pkd, &star, in.dSNIaMaxMass, in.dCCSNMinMass,
                                     in.dCCSNMaxMass);
#endif

#ifdef FEEDBACK
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
                *dMassFormed += fMass;

                pkd->nGas -= 1;
                pkd->nStar += 1;
            }
        }
    }
}
#endif

static inline double pressure_SFR(PKD pkd, double a_m3, double dDenMin,
                                  double dSFnormalizationKS, double dSFindexKS, double dSFGasFraction,
                                  double dEOSPolyFloorIndex, double dEOSPolyFloorDen,double dEOSPolyFlooru,
                                  double dConstGamma, particleStore::ParticleReference &p, SPHFIELDS &sph) {

    float fMass = p.mass();

    const double hyd_abun = sph.afElemMass[ELEMENT_H] / fMass;
    const double rho_H = p.density() * hyd_abun;


    // Two SF thresholds are applied:
    //      a) Minimum density, computed at the master level
    //      b) Maximum temperature of a
    //            factor 0.5 dex (i.e., 3.1622) above the polytropic eEOS
#ifdef EEOS_POLYTROPE
    const double maxUint = 3.16228 * fMass *
                           polytropicEnergyFloor(a_m3, p.density(),
                                   dEOSPolyFloorIndex, dEOSPolyFloorDen,  dEOSPolyFlooru);
#else
    const double maxUint = INFINITY;
#endif

    if (sph.Uint > maxUint || rho_H < dDenMin) {
        return 0.0;
    }


    const double SFexp = 0.5*(dSFindexKS-1.);

    const double dmstar =  dSFnormalizationKS * p.mass() *
                           pow( dConstGamma*dSFGasFraction*sph.P*a_m3, SFexp);

    return dmstar;
}

static inline double density_SFR(PKD pkd, double a_m3, double dDenMin,
                                 double dSFThresholdu, double dSFEfficiency,
                                 particleStore::ParticleReference &p, SPHFIELDS &sph) {

    float fMass = p.mass();
    float fDens = p.density();

    const double hyd_abun = sph.afElemMass[ELEMENT_H] / fMass;
    const double rho_H = fDens*hyd_abun;


    // Two SF thresholds are applied:
    //      a) Minimum density, computed at the master level
    //      b) Maximum temperature set by dSFThresholdTemp
    // const double dens = fDens*a_m3;
    const double maxUint = dSFThresholdu;

    if (sph.Uint > maxUint || rho_H < dDenMin) {
        return 0.0;
    }

    const double tff = sqrt(3.*M_PI/(32.*fDens*a_m3));

    const double dmstar = dSFEfficiency * fDens / (tff * sph.omega);

    return dmstar;
}
