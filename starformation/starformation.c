
#ifdef  STAR_FORMATION
#include "master.h"
#include "pkd.h"

/* IA: MSR layer
 */
void msrStarFormInit(MSR msr, double dTime){
   struct inStarForm in;
   struct outStarForm out;
   in.dTime = dTime;

   pstStarFormInit(msr->pst, &in, sizeof(in), &out, sizeof(out));

   printf("%d star particles are about to explode in this IC\n", out.nFormed);
}

void msrStarForm(MSR msr, double dTime, double dDelta, int iRung)
    {
    struct inStarForm in;
    struct outStarForm out;
    double sec,sec1,dsec;

    sec = msrTime();

    const double a = csmTime2Exp(msr->csm,dTime);
    const double H = csmTime2Hub(msr->csm,dTime);


    /* Convert input parameters to code units */
    /*
    in.dRateCoeff = msr->param.SFdEfficiency*sqrt(32/(3*M_PI)/pow(a,3)); 
    in.dTMax = msr->param.SFdTMax;
    d1 = msr->param.SFdComovingDenMin;
    d2 = msr->param.SFdPhysDenMin/msr->param.dGmPerCcUnit*pow(a,3);
    in.dDenMin = (d1>d2 ? d1 : d2);
    */
    
    in.dScaleFactor = a;
    in.dTime = dTime;
    in.dDelta = dDelta;
    /*
    in.dInitStarMass = msr->param.SFdInitStarMass;
    in.dESNPerStarMass = msr->param.SFdESNPerStarMass/msr->param.dErgPerGmUnit;
#define SECONDSPERYEAR   31557600.
    in.dtCoolingShutoff = msr->param.SFdtCoolingShutoff*SECONDSPERYEAR/msr->param.dSecUnit;
    in.dtFeedbackDelay = msr->param.SFdtFeedbackDelay*1.0000013254678*SECONDSPERYEAR/msr->param.dSecUnit;
    in.dMassLossPerStarMass = msr->param.SFdMassLossPerStarMass;
    in.dZMassPerStarMass = msr->param.SFdZMassPerStarMass;
    in.dInitStarMass = msr->param.SFdInitStarMass;
    in.dMinGasMass = msr->param.SFdMinGasMass;
    in.bdivv = msr->param.SFbdivv;
    */


    // We convert the threshold density (given in cgs in the parameters file) in code units
    //  NOTE: We still have to divide by the hydrogen fraction of each particle!
    if (msr->csm->val.bComove)
       in.dDenMin = msr->param.dSFThresholdDen*pow(a,3);
    else
       in.dDenMin = msr->param.dSFThresholdDen;

    // Critical density of the Universe, in code units
    in.dDenCrit = 1.*0.048; // TODO: Read from parameters

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
    
    sec1 = msrTime();
    dsec = sec1 - sec;
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);


    }



/* IA: PKD layer
 */
void pkdStarFormInit(PKD pkd, double dTime, int *nFormed){
    *nFormed = 0;
    for (int i=0;i<pkdLocal(pkd);++i) {
       PARTICLE *p = pkdParticle(pkd,i);
      if (pkdIsStar(pkd,p)){
          STARFIELDS* pStar = pkdStar(pkd,p);
          if (pStar->fTimer != 0){
            if ( (dTime-pStar->fTimer) < pkd->param.dFeedbackDelay){
               pStar->hasExploded = 0; // This particle did not explode before the last snapshot
               *nFormed += 1;
            }
          }
      }
    }
}



void pkdStarForm(PKD pkd, 
             double dTime,
             double dDelta,
             double dScaleFactor,
             double dDenMin, /* Threshold for SF in code units  */
             double dDenCrit, /* Multiple of the critical density needed for SF*/
		 int *nFormed, /* number of stars formed */
		 double *dMassFormed,	/* mass of stars formed */
		 int *nDeleted) /* gas particles deleted */ {

    PARTICLE *p;
    SPHFIELDS *psph;
    double T, E, dt, prob;
    float* pv;
    PARTICLE *starp;
    int i;
    
    assert(pkd->oStar);
    assert(pkd->oSph);
    assert(pkd->oMass);

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;

    double a_m2 = 1./(dScaleFactor*dScaleFactor);
    double a_m3 = a_m2/dScaleFactor;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	
	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    psph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung); 
          dt = dTime - psph->lastUpdateTime;

#ifdef COOLING
          const double rho_H = pkdDensity(pkd,p) * psph->chemistry[chemistry_element_H];
#else
          const double rho_H = pkdDensity(pkd,p) * 0.75; // If no information, assume primoridal abundance
#endif

          // Gas too hot is not allowed to form stars (but if it is too collapses, we allow to form stars as it is in the EoS
          if ( (psph->Uint > pkd->param.dSFThresholdu*pkdMass(pkd,p)) && (rho_H*a_m3<10*pkd->param.dCoolingFloorDen) ){
             psph->SFR=0.; 
             continue;
          }
          // Thresholds for star formation
	    if (pkd->csm->val.bComove &&  (pkdDensity(pkd,p) < pkd->param.dSFMinOverDensity*dDenCrit) ){
             psph->SFR=0.; 
             continue;
          }
          if (rho_H < dDenMin){
             psph->SFR=0.; 
             continue;
          }
          const double dmstar = pkd->param.dSFnormalizationKS*  pkdMass(pkd,p) *  pow(a_m2, 1.4) *
                        pow( pkd->param.dConstGamma * pkd->param.dSFGasFraction * psph->P*a_m3, pkd->param.dSFindexKS);
          psph->SFR = dmstar;
	    const double prob = 1.0 - exp(-dmstar*dt/pkdMass(pkd,p)); 
          //printf("%e \n", prob);
	    
	    /* Star formation event? */
	    if (rand()<RAND_MAX*prob) {

             // We just change the class of the particle to stellar one
            pkdSetClass(pkd, pkdMass(pkd,p), 0., FIO_SPECIES_STAR, p);

            // When changing the class, we have to take into account that the code velocity
            // has different scale factor dependencies for dm/star particles and gas particles
            pv = pkdVel(pkd,p);
            for (int j=0; j<3; j++){
               pv[j] *= dScaleFactor;
            }
            // We log statistics about the formation time
            pkdStar(pkd, p)->fTimer = dTime;
            pkdStar(pkd, p)->hasExploded = 0;

            // Safety check, TODO: remove if nothing happens:
            assert(pkdIsStar(pkd,p));
            assert(!pkdIsGas(pkd,p));

		(*nFormed)++;
		*dMassFormed += pkdMass(pkd,p);
		}
	    }
	}

}
#endif 
