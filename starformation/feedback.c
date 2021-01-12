#ifdef FEEDBACK
#include "pkd.h"
#include "smooth.h"
#include "starformation/feedback.h"


/* Function that will be called with the information of all the neighbors. 
 * Here we compute the probability of explosion, and we add the energy to the surroinding gas particles
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
    
    // We could unify this factors in order to avoid more operations (although I think wont make a huge change anyway)
    // There is no need to explicitly convert to solar masses as we are doing a ratio
    const double prob = pkd->param.dFeedbackEfficiency * pkd->param.dNumberSNIIperMass * pkdMass(pkd,p) /(pkd->param.dFeedbackDu*totMass);
    //printf("prob %e \t dFeedbackEfficiency %e \t dNumberSNIIperMass %e \t mass/totMass %e \t dFeedbackDu %e \n", prob,
    //        pkd->param.dFeedbackEfficiency,      pkd->param.dNumberSNIIperMass,      pkdMass(pkd,p)/totMass, pkd->param.dFeedbackDu);
    assert(prob<1.0);
    for (i=0; i<nSmooth; ++i){
          if (nnList[i].dx==0 && nnList[i].dy==0 && nnList[i].dz==0) continue;

	    if (rand()<RAND_MAX*prob) { // We have a supernova explosion!
             q = nnList[i].pPart;
             qsph = pkdSph(pkd,q);
             //printf("Uint %e extra %e \n", qsph->Uint, pkd->param.dFeedbackDu * pkdMass(pkd,q));
             
             const double feed_energy = pkd->param.dFeedbackDu * pkdMass(pkd,q);
             qsph->Uint += feed_energy;
             qsph->E += feed_energy;
#ifdef ENTROPY_SWITCH
             qsph->S += feed_energy*(pkd->param.dConstGamma-1.)*pow(pkdDensity(pkd,q), -pkd->param.dConstGamma+1);
#endif
             //
             printf("Adding SN energy! \n");
          }

    }
}



void initSNFeedback(void *vpkd, void *vp){
   PKD pkd = (PKD) vpkd;

   PARTICLE *p = vp;

   SPHFIELDS *psph = pkdSph(pkd,p);

   psph->Uint = 0.;
   psph->E = 0.;
#ifdef ENTROPY_SWITCH
   psph->S = 0.;
#endif

}


void combSNFeedback(void *vpkd, void *p1,void *p2){
   PKD pkd = (PKD) vpkd;

   
   SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);

   psph1->Uint += psph2->Uint;
   psph1->E += psph2->E;
#ifdef ENTROPY_SWITCH
   psph1->S += psph2->S;
#endif


}

#endif
