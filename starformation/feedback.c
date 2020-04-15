#ifdef FEEDBACK
#include "pkd.h"
#include "smooth.h"
#include "starformation/feedback.h"


/* Function that will be called with the information of all the neighbors. 
 * Here we compute the probability of explosion, and we add the energy to the surroinding gas particles
 */
void smFeedback(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
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
             qsph->Uint += pkd->param.dFeedbackDu * pkdMass(pkd,q);
             qsph->E += pkd->param.dFeedbackDu * pkdMass(pkd,q);
             printf("Adding SN energy! \n");
          }

    }
}
#endif
