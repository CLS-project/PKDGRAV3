
#ifdef  STAR_FORMATION
#include "master.h"
#include "pkd.h"

/* IA: MSR layer
 */
void msrStarForm(MSR msr, double dTime, int iRung)
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
    in.dDenMin = msr->param.SFdPhysDenMin/pow(a,3) * MHYDR / (msr->param.dMsolUnit * MSOLG );

    // Critical density of the Universe, in code units
    in.dDenCrit = 3.*H*H*M_1_PI/8.;

    if (msr->param.bVDetails) printf("Star Form ... ");
    
    msrActiveRung(msr,iRung,1); /* important to limit to active gas only */
    pstStarForm(msr->pst, &in, sizeof(in), &out, NULL);
    if (msr->param.bVDetails)
	printf("%d Stars formed with mass %g, %d gas deleted\n",
	       out.nFormed, out.dMassFormed, out.nDeleted);
    
    
    sec1 = msrTime();
    dsec = sec1 - sec;
    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);

    if (msr->param.bFeedback) {
       //IA: TODO
	msrActiveRung(msr,iRung,1); /* costs nothing -- important to limit to active stars only */
 	//msrSelSrcGas(msr); /* Not really sure what the setting here needs to be */
	//msrSelDstStar(msr,1,dTime); /* Select only stars that have FB to do */ 
/*	msrBuildTree(msr,dTime,msr->param.bEwald);*/
	msrSmooth(msr, dTime, SMX_DIST_SN_ENERGY, 1, msr->param.nSmooth); /* full smooth for stars */
	//msrSelSrcAll(msr);
	//msrSelDstAll(msr);

	dsec = msrTime() - sec1;
	printf("Feedback Calculated, Wallclock: %f secs\n\n",dsec);
	}
    }



/* IA: PKD layer
 */
void pkdStarForm(PKD pkd, 
             double dTime,
             double dDenMin, /* Threshold for SF in code units  */
             double dDenCrit, /* Multiple of the critical density needed for SF*/
		 int *nFormed, /* number of stars formed */
		 double *dMassFormed,	/* mass of stars formed */
		 int *nDeleted) /* gas particles deleted */ {

    PARTICLE *p;
    SPHFIELDS *sph;
    double T, E, dmstar, dt, prob;
    PARTICLE *starp;
    int i;
    
    assert(pkd->oStar);
    assert(pkd->oSph);
    assert(pkd->oMass);

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;
    assert(starp != NULL);

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	
	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    sph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung); 

#ifdef COOLING
          const double rho_H = pkdDensity(pkd,p) * sph->chemistry[chemistry_element_H];
#else
          const double rho_H = pkdDensity(pkd,p) * 0.75; // If no information, assume primoridal abundance
#endif

	    if (rho_H < dDenMin && pkdDensity(pkd,p) > pkd->param.SFdMinOverDensity*dDenCrit) continue;

          const double mass_Mo = pkdMass(pkd,p) * pkd->param.dMsolUnit;
          const double P_kb = sph->P * pkd->param.dErgUnit / KBOLTZ; //TODO: double-check this
	    const double mstar = 5.99e-10*mass_Mo*pow( pkd->param.SFdGasFraction * P_kb * 1e-3  , 0.2);
	    const double prob = 1.0 - exp(-dmstar*dt/mass_Mo); 
	    
	    /* Star formation event? */
	    if (rand()<RAND_MAX*prob) {

             // We just change the class of the particle to stellar one
            pkdSetClass(pkd, pkdMass(pkd,p), pkdSoft(pkd,p), FIO_SPECIES_STAR, p);

            // We log statistics about the formation time
            pkdStar(pkd, p)->fTimer = dTime;

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
