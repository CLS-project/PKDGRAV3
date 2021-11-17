#ifdef FEEDBACK
#include "pkd.h"
#include "smooth.h"
#include "starformation/feedback.h"

void msrSetFeedbackParam(MSR msr){
    const double dHydFrac = msr->param.dInitialH;
    const double dnHToRho = MHYDR / dHydFrac / msr->param.dGmPerCcUnit;
    msr->param.dFeedbackDu *= msr->param.dTuFac;
    msr->param.dFeedbackDelay *=  SECONDSPERYEAR/msr->param.dSecUnit ;
    msr->param.dNumberSNIIperMass *= 8.73e15 / msr->param.dErgPerGmUnit / 1.736e-2;
    msr->param.dFeedbackEffnH0 *= dnHToRho;
}


/* Function that will be called with the information of all the neighbors.
 * Here we compute the probability of explosion, and we add the energy to the
 * surroinding gas particles
 */
void smSNFeedback(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *qsph;
    int i;

    float totMass = 0;
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       totMass += pkdMass(pkd,q);
    }
    totMass -= pkdMass(pkd,p);

    // We could unify this factors in order to avoid more operations
    // (although I think wont make a huge change anyway)
    // There is no need to explicitly convert to solar masses
    // becuse we are doing a ratio
    const double prob = pkdStar(pkd,p)->fSNEfficiency *
                        pkd->param.dNumberSNIIperMass *
                        pkdMass(pkd,p) / (pkd->param.dFeedbackDu*totMass);


    assert(prob<1.0);
    for (i=0; i<nSmooth; ++i){
          if (nnList[i].dx==0 && nnList[i].dy==0 && nnList[i].dz==0) continue;

	    if (rand()<RAND_MAX*prob) { // We have a supernova explosion!
             q = nnList[i].pPart;
             qsph = pkdSph(pkd,q);
             //printf("Uint %e extra %e \n",
             //   qsph->Uint, pkd->param.dFeedbackDu * pkdMass(pkd,q));

             const double feed_energy = pkd->param.dFeedbackDu * pkdMass(pkd,q);
#ifdef OLD_FB_SCHEME
             qsph->Uint += feed_energy;
             qsph->E += feed_energy;
#ifdef ENTROPY_SWITCH
             qsph->S += feed_energy*(pkd->param.dConstGamma-1.) *
                        pow(pkdDensity(pkd,q), -pkd->param.dConstGamma+1);
#endif
#else // OLD_BH_SCHEME
             qsph->fAccFBEnergy += feed_energy;
#endif

             //printf("Adding SN energy! \n");
          }

    }
}



void initSNFeedback(void *vpkd, void *vp){
   PKD pkd = (PKD) vpkd;

   PARTICLE *p = vp;

   if (pkdIsGas(pkd,p)){
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


void combSNFeedback(void *vpkd, void *p1,void *p2){
   PKD pkd = (PKD) vpkd;

   if (pkdIsGas(pkd,p1) && pkdIsGas(pkd,p2)){

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

inline void pkdAddFBEnergy(PKD pkd, PARTICLE* p, SPHFIELDS *psph){
#ifndef OLD_FB_SCHEME
   psph->Uint += psph->fAccFBEnergy;
   psph->E += psph->fAccFBEnergy;
#ifdef ENTROPY_SWITCH
   psph->S += psph->fAccFBEnergy*(pkd->param.dConstGamma-1.) *
            pow(pkdDensity(pkd,p), -pkd->param.dConstGamma+1);
#endif
   psph->fAccFBEnergy = 0.0;
#endif //OLD_FB_SCHEME
}

inline float SNFeedbackEfficiency(PKD pkd, float Z, float rho){
   if (pkd->param.dFeedbackMaxEff > 0.0){
      const double den = 1.0 + pow(Z/0.00127, pkd->param.dFeedbackEffIndex) *
                  pow(rho/pkd->param.dFeedbackEffnH0,-pkd->param.dFeedbackEffIndex);
      return pkd->param.dFeedbackEfficiency +
                  (pkd->param.dFeedbackMaxEff - pkd->param.dFeedbackEfficiency)/
                  den;
   }else{
      return pkd->param.dFeedbackEfficiency;
   }
}

#endif // FEEDBACK
