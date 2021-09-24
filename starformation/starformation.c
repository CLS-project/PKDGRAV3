#ifdef  STAR_FORMATION
#include "master.h"
#include "pkd.h"

#ifdef STELLAR_EVOLUTION
#include "stellarevolution/stellarevolution.h"
#endif

/* IA: MSR layer
 */

void msrStarForm(MSR msr, double dTime, double dDelta, int iRung)
    {
    struct inStarForm in;
    struct outStarForm out;
    double sec,sec1,dsec;

    msrTimerStart(msr, TIMER_STARFORM);

    const double a = csmTime2Exp(msr->csm,dTime);

    in.dScaleFactor = a;
    in.dTime = dTime;
    in.dDelta = dDelta;

    // Here we set the minium density a particle must have to be SF
    //  NOTE: We still have to divide by the hydrogen fraction of each particle!
    if (msr->csm->val.bComove){
       double a3 = a*a*a;
       in.dDenMin = msr->param.dSFThresholdDen*a3;
       assert(msr->csm->val.dOmegab  > 0.);

       // If in PKDGRAV3 units, this should always be unity
       double rhoCrit0 = 3. * msr->csm->val.dHubble0 * msr->csm->val.dHubble0 /
                              (8. * M_PI);
       double denCosmoMin = rhoCrit0 * msr->csm->val.dOmegab *
                                 msr->param.dSFMinOverDensity *
                                 msr->param.dInitialH;
       in.dDenMin = (in.dDenMin > denCosmoMin) ? in.dDenMin : denCosmoMin;
    }else{
       in.dDenMin = msr->param.dSFThresholdDen;
    }


    if (msr->param.bVDetails) printf("Star Form (rung: %d) ... ", iRung);

    msrActiveRung(msr,iRung,1);
    pstStarForm(msr->pst, &in, sizeof(in), &out, 0);
    if (msr->param.bVDetails)
	printf("%d Stars formed with mass %g, %d gas deleted\n",
	       out.nFormed, out.dMassFormed, out.nDeleted);
    msr->massFormed += out.dMassFormed;
    msr->starFormed += out.nFormed;

    msr->nGas -= out.nFormed;
    msr->nStar += out.nFormed;

    msrTimerStop(msr, TIMER_STARFORM);
    dsec = msrTimerGet(msr, TIMER_STARFORM);
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);


    }






void pkdStarForm(PKD pkd,
		 double dTime,
		 double dDelta,
		 double dScaleFactor,
		 double dDenMin, /* co-moving threshold for SF in code units  */
		 int *nFormed, /* number of stars formed */
		 double *dMassFormed,	/* mass of stars formed */
		 int *nDeleted) /* gas particles deleted */ {

    PARTICLE *p;
    SPHFIELDS *psph;
    double dt;
    float* pv;
    int i, j;

    assert(pkd->oStar);
    assert(pkd->oSph);
    assert(pkd->oMass);

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;

    const double a_m1 = 1.0/dScaleFactor;
    const double a_m2 = a_m1*a_m1;
    const double a_m3 = a_m2*a_m1;
    const double SFexp = (pkd->param.dSFindexKS-1.)/2.;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);

	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    psph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung);
          dt = dTime - psph->lastUpdateTime;

          float fMass = pkdMass(pkd, p);
          // If no information, assume primoridal abundance
#ifdef COOLING
          const double hyd_abun = psph->afElemMass[ELEMENT_H] / fMass;
#else
          const double hyd_abun = pkd->param.dInitialH;
#endif

          const double rho_H = pkdDensity(pkd,p) * hyd_abun;


          // Two SF thresholds are applied:
          //      a) minimum density, computed at the master level
          //      b) Maximum temperature of a
          //            factor 0.5 dex (i.e., 3.1622) above the eEOS
          double dens = pkdDensity(pkd,p) * a_m3;
          double maxUint = 3.16228 * fMass * pkd->param.dJeansFlooru *
             pow(dens/pkd->param.dJeansFloorDen, pkd->param.dJeansFloorIndex);

          if (psph->Uint > maxUint || rho_H < dDenMin) {
             psph->SFR = 0.;
             continue;
          }


          const double dmstar =
          pkd->param.dSFnormalizationKS * pkdMass(pkd,p) *
          pow( pkd->param.dConstGamma*pkd->param.dSFGasFraction*psph->P*a_m3,
               SFexp);

          psph->SFR = dmstar;

          const double prob = 1.0 - exp(-dmstar*dt/pkdMass(pkd,p));

          // Star formation event?
          if (rand()<RAND_MAX*prob) {

            //printf("STARFORM %e %e %e \n", dScaleFactor, rho_H, psph->Uint);

#ifdef STELLAR_EVOLUTION
            float afElemMass[ELEMENT_COUNT];
            float fMetalMass;
            for (j = 0; j < ELEMENT_COUNT; j++)
                afElemMass[j] = psph->afElemMass[j];
            fMetalMass = psph->fMetalMass;
#endif

            // We just change the class of the particle to stellar one
            pkdSetClass(pkd, pkdMass(pkd,p), pkdSoft0(pkd,p), FIO_SPECIES_STAR, p);
	      STARFIELDS *pStar = pkdStar(pkd, p);

            // When changing the class, we have to take into account that
            // the code velocity has different scale factor dependencies for
            // dm/star particles and gas particles
            pv = pkdVel(pkd,p);
            for (j=0; j<3; j++){
               pv[j] *= dScaleFactor;
            }

            // We log statistics about the formation time
            pStar->fTimer = dTime;
            pStar->hasExploded = 0;

#ifdef STELLAR_EVOLUTION
            for (j = 0; j < ELEMENT_COUNT; j++)
               pStar->afElemAbun[j] = afElemMass[j] / fMass;
            pStar->fMetalAbun = fMetalMass / fMass;

            pStar->fInitialMass = fMass;
            pStar->fLastEnrichTime = 0.0f;

            if (pkd->param.bChemEnrich)
               stevStarParticleInit(pkd, pStar);
            else
               pStar->fNextEnrichTime = INFINITY;
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
