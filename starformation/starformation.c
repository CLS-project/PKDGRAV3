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

    const double a = csmTime2Exp(msr->csm,dTime);
    const double H = csmTime2Hub(msr->csm,dTime);
    
    in.dScaleFactor = a;
    in.dTime = dTime;
    in.dDelta = dDelta;

    // Here we set the minium density a particle must have to be SF
    //  NOTE: We still have to divide by the hydrogen fraction of each particle!
    if (msr->csm->val.bComove){
       in.dDenMin = msr->param.dSFThresholdDen*pow(a,3);
       assert(msr->csm->val.dOmegab  > 0.);

       // It is assumed that for cosmo runs, rho_crit=1 in code units
       double min_cosmo_den = msr->csm->val.dOmegab * msr->param.dSFMinOverDensity;
       in.dDenMin = (in.dDenMin > min_cosmo_den) ? in.dDenMin : min_cosmo_den;
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
    
    sec1 = msrTime();
    dsec = sec1 - sec;
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);


    }






void pkdStarForm(PKD pkd, 
             double dTime,
             double dDelta,
             double dScaleFactor,
             double dDenMin, /* Threshold for SF in code units  */
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

          // If no information, assume primoridal abundance
#ifdef COOLING
          const double hyd_abun = psph->chemistry[chemistry_element_H];
#else
          const double rho_H = 0.75; 
#endif

          const double rho_H = pkdDensity(pkd,p) * hyd_abun;


          // Two SF thresholds are applied:
          //      a) minimum density, computed at the master level
          //      b) Maximum temperature of a
          //            factor 0.5 dex (i.e., 3.1622) above the eEOS
          double dens = pkdDensity(pkd,p) * a_m3;
          double maxUint = 3.16228 * pkdMass(pkd,p) * pkd->param.dJeansFlooru *
           pow( dens/pkd->param.dJeansFloorDen , pkd->param.dJeansFloorIndex );
            
          if (psph->Uint > maxUint || rho_H < dDenMin) {
             psph->SFR = 0.;
             continue;
          }


          const double dmstar = 
             pkd->param.dSFnormalizationKS * pkdMass(pkd,p) *  pow(a_m2, 1.4) *
             pow( pkd->param.dConstGamma * pkd->param.dSFGasFraction * psph->P*a_m3,
                  pkd->param.dSFindexKS);

          psph->SFR = dmstar;

	    const double prob = 1.0 - exp(-dmstar*dt/pkdMass(pkd,p)); 
          //printf("%e \n", prob);
	    
	    // Star formation event?
	    if (rand()<RAND_MAX*prob) {

            //printf("STARFORM %e %e %e \n", dScaleFactor, rho_H, psph->Uint);

            // We just change the class of the particle to stellar one
            pkdSetClass(pkd, pkdMass(pkd,p), 0., FIO_SPECIES_STAR, p);

            // When changing the class, we have to take into account that 
            // the code velocity has different scale factor dependencies for
            // dm/star particles and gas particles
            pv = pkdVel(pkd,p);
            for (int j=0; j<3; j++){
               pv[j] *= dScaleFactor;
            }
            // We log statistics about the formation time
            pkdStar(pkd, p)->fTimer = dTime;
            pkdStar(pkd, p)->hasExploded = 0;

            // Safety check
            assert(pkdIsStar(pkd,p));
            assert(!pkdIsGas(pkd,p));

		(*nFormed)++;
		*dMassFormed += pkdMass(pkd,p);

            pkd->nGas -= 1;
            pkd->nStar += 1;
		}
	    }
	}

}
#endif 
