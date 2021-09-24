#ifdef  STAR_FORMATION
#include "starformation/starformation.h"


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



int pstStarForm(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inStarForm *in = vin;
    struct outStarForm *out = vout;
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
	pkdStarForm(pst->plcl->pkd,
		    in->dTime, in->dDelta, in->dScaleFactor, in->dDenMin,
		     &out->nFormed, &out->dMassFormed, &out->nDeleted);
	}
    return sizeof(struct outStarForm);
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

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);

	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    psph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung);
          dt = dTime - psph->lastUpdateTime;

          double dmstar;
          if (pkd->param.dSFEfficiency > 0.0){
             dmstar = density_SFR(pkd, a_m3, dDenMin, p, psph);
          }else{
             dmstar = pressure_SFR(pkd, a_m3, dDenMin, p, psph);
          }

          psph->SFR = dmstar;

          const double prob = 1.0 - exp(-dmstar*dt/pkdMass(pkd,p));

          // Star formation event?
          if (rand()<RAND_MAX*prob) {
            float fMass = pkdMass(pkd, p);

            //printf("STARFORM %e %e %e \n", dScaleFactor, rho_H, psph->Uint);

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

static inline double pressure_SFR(PKD pkd, double a_m3, double dDenMin,
      PARTICLE *p, SPHFIELDS *psph){

   float fMass = pkdMass(pkd, p);

   // If no information, assume primordial abundance
#ifdef COOLING
   const double hyd_abun = psph->afElemMass[ELEMENT_H] / fMass;
#else
   const double hyd_abun = pkd->param.dInitialH;
#endif

   const double rho_H = pkdDensity(pkd,p) * hyd_abun;


   // Two SF thresholds are applied:
   //      a) Minimum density, computed at the master level
   //      b) Maximum temperature of a
   //            factor 0.5 dex (i.e., 3.1622) above the eEOS
   const double dens = pkdDensity(pkd,p) * a_m3;
   const double maxUint = 3.16228 * fMass * pkd->param.dJeansFlooru *
   pow( dens/pkd->param.dJeansFloorDen , pkd->param.dJeansFloorIndex );

   if (psph->Uint > maxUint || rho_H < dDenMin) {
      return 0.0;
   }


   const double SFexp = 0.5*(pkd->param.dSFindexKS-1.);

   const double dmstar =
   pkd->param.dSFnormalizationKS * pkdMass(pkd,p) *
   pow( pkd->param.dConstGamma*pkd->param.dSFGasFraction*psph->P*a_m3,
      SFexp);

   return dmstar;
}

static inline double density_SFR(PKD pkd, double a_m3, double dDenMin,
      PARTICLE *p, SPHFIELDS *psph){

   float fMass = pkdMass(pkd, p);
   float fDens = pkdDensity(pkd, p);

   // If no information, assume primordial abundance
#ifdef COOLING
   const double hyd_abun = psph->afElemMass[ELEMENT_H] / fMass;
#else
   const double hyd_abun = pkd->param.dInitialH;
#endif

   const double rho_H = fDens*hyd_abun;


   // Two SF thresholds are applied:
   //      a) Minimum density, computed at the master level
   //      b) Maximum temperature of a
   //            factor 0.5 dex (i.e., 3.1622) above the eEOS
   const double dens = fDens*a_m3;
   const double maxUint = pkd->param.dSFThresholdu;

   if (psph->Uint > maxUint || rho_H < dDenMin) {
      return 0.0;
   }

   const double tff = sqrt(3.*M_PI/(32.*fDens*a_m3));

   const double dmstar = pkd->param.dSFEfficiency * fDens /
                              (tff * psph->omega);

   return dmstar;
}
