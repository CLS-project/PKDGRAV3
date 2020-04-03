
#ifdef  STAR_FORMATION
#include "master.h"
#include "pkd.h"

/* IA: MSR layer
 */
void msrStarForm(MSR msr, double dTime, double dDelta, int iRung)
    {
    struct inStarForm in;
    struct outStarForm out;
    double sec,sec1,dsec;

    sec = msrTime();

    const double a = csmTime2Exp(msr->param.csm,dTime);
    const double H = csmTime2Hub(msr->param.csm,dTime);


    /* Convert input parameters to code units */
    /*
    in.dRateCoeff = msr->param.SFdEfficiency*sqrt(32/(3*M_PI)/pow(a,3)); 
    in.dTMax = msr->param.SFdTMax;
    d1 = msr->param.SFdComovingDenMin;
    d2 = msr->param.SFdPhysDenMin/msr->param.dGmPerCcUnit*pow(a,3);
    in.dDenMin = (d1>d2 ? d1 : d2);
    */
    
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
    if (msr->param.csm->val.bComove)
       in.dDenMin = msr->param.dSFThresholdDen*pow(a,-3);
    else
       in.dDenMin = msr->param.dSFThresholdDen;

    // Critical density of the Universe, in code units
    in.dDenCrit = 3.*H*H*M_1_PI/8.;

    if (msr->param.bVDetails) printf("Star Form ... ");
    
    msrActiveRung(msr,iRung,0); /* important to limit to active gas only */
    pstStarForm(msr->pst, &in, sizeof(in), &out, NULL);
    if (msr->param.bVDetails)
	printf("%d Stars formed with mass %g, %d gas deleted\n",
	       out.nFormed, out.dMassFormed, out.nDeleted);
    msr->massFormed += out.dMassFormed;
    msr->starFormed += out.nFormed;    
    
    sec1 = msrTime();
    dsec = sec1 - sec;
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);

    }



/* IA: PKD layer
 */
void pkdStarForm(PKD pkd, 
             double dTime,
             double dDelta,
             double dDenMin, /* Threshold for SF in code units  */
             double dDenCrit, /* Multiple of the critical density needed for SF*/
		 int *nFormed, /* number of stars formed */
		 double *dMassFormed,	/* mass of stars formed */
		 int *nDeleted) /* gas particles deleted */ {

    PARTICLE *p;
    SPHFIELDS *psph;
    double T, E, dt, prob;
    PARTICLE *starp;
    int i;
    
    assert(pkd->oStar);
    assert(pkd->oSph);
    assert(pkd->oMass);

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;
    //assert(starp != NULL);

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	
	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    psph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung); //* pkd->param.dSecUnit / 3600. / 24. / 365.; 
          //printf("%d %e %e \n", p->uRung, dt, dDelta);
          assert(fabs(dt-dDelta)<1e-6);

#ifdef COOLING
          const double rho_H = pkdDensity(pkd,p) * psph->chemistry[chemistry_element_H];
#else
          const double rho_H = pkdDensity(pkd,p) * 0.75; // If no information, assume primoridal abundance
#endif

          // Gas too hot is not allowed to form stars
          if (psph->Uint > pkd->param.dSFThresholdu*pkdMass(pkd,p)) continue;

          // Thresholds for star formation
	    if (pkd->param.csm->val.bComove &&  (pkdDensity(pkd,p) > pkd->param.dSFMinOverDensity*dDenCrit) ) continue;
          if (rho_H < dDenMin) continue;

          const double mass_Mo = pkdMass(pkd,p) * pkd->param.dMsolUnit;
          const double P_kb = psph->P * pkd->param.dErgUnit * pow(pkd->param.dKpcUnit*KPCCM, -3)/ KBOLTZ; 


/*
	    const double dmstar0 = 5.99e-10*mass_Mo*pow( pkd->param.dSFGasFraction * P_kb * 1e-3  , 0.2);

          const double Msolpcm2 = 1. / pkd->param.dMsolUnit * pow(pkd->param.dKpcUnit*1e3, 2);
          const double A = 2.5e-4 / pkd->param.dMsolUnit * pkd->param.dSecUnit/(3600*24*365) * pow(pkd->param.dKpcUnit, 2);
          const double n = 1.4;
          const double dmstar = A *  pow(Msolpcm2, -n)  *  pkdMass(pkd,p) *
                                     pow( pkd->param.dConstGamma * pkd->param.dSFGasFraction * psph->P, (n-1.)/2.  );
  */           
          const double dmstar = pkd->param.dSFnormalizationKS * pkdMass(pkd,p) * 
                        pow( pkd->param.dConstGamma * pkd->param.dSFGasFraction * psph->P, pkd->param.dSFindexKS);
          //printf("mstar %e mass_Mo %e P_kb %e \t\t prob %e \n", dmstar*dt, mass_Mo, P_kb, dmstar*dt/mass_Mo);
          //printf("dConstGamma %e \t SFGasFraction %e \n", pkd->param.dConstGamma, pkd->param.dSFGasFraction);
	    const double prob = 1.0 - exp(-dmstar*dt/pkdMass(pkd,p)); 
          //printf("%e %e \n", prob, 1. - exp(-dmstar0/mass_Mo * dt* pkd->param.dSecUnit / 3600. / 24. / 365.) ); Dan lo mismo!
	    
	    /* Star formation event? */
	    if (rand()<RAND_MAX*prob) {

             // We just change the class of the particle to stellar one
            pkdSetClass(pkd, pkdMass(pkd,p), pkdSoft(pkd,p), FIO_SPECIES_STAR, p);

            // We log statistics about the formation time
            pkdStar(pkd, p)->fTimer = dTime;
            pkdStar(pkd, p)->hasExploded = 0;

            // Safety check, TODO: remove if nothing happens:
            assert(pkdIsStar(pkd,p));
            assert(!pkdIsGas(pkd,p));
#ifdef FEEDBACK
            // We set the internal energy to that of a exploding particle, which will be used to set the timestep
            //  of the neighbours just before exploding. This way we avoid 'waking up' particles
            //
            //  Note that we can do this because, although the class of p is a star, 
            //  it still holds the SPHFIELD struct (which is suboptimal, but hey, we can take it as an advantage)
            psph->Uint += pkdMass(pkd,p) * pkd->param.dFeedbackDu;
            psph->P = psph->Uint*psph->omega*(pkd->param.dConstGamma -1.);
            psph->c = sqrt(psph->Uint*pkd->param.dConstGamma/pkdMass(pkd,p)/(pkd->param.dConstGamma -1.));
            psph->c *= 4.; // Even more tight constraint
#endif
		(*nFormed)++;
		*dMassFormed += pkdMass(pkd,p);
		}
	    }
	}

}
#endif 
