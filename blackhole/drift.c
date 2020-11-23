#include "blackhole/drift.h"
#include "hydro.h"


void smBHdrift(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf){

    PKD pkd = smf->pkd;
    // We look for the most bounded neighbouring particle
    PARTICLE *pLowPot = NULL;
    float minPot = HUGE_VAL;
    float pMass = pkdMass(pkd,p);
    uint8_t uMaxRung = 0;
    double r_lowPot[3], r_p[3];
    double cs = 0.;
    float inv_a = 1./smf->a;


    for (int i=0; i<nSmooth; ++i){
       PARTICLE* q = nnList[i].pPart;
       assert(pkdIsGas(pkd,q));
       pLowPot = (*pkdPot(pkd,q)<minPot)  ? q : pLowPot;
       minPot = (pLowPot!=NULL) ? *pkdPot(pkd,pLowPot) : minPot; // We could have no potential if gravity is not calculated
       uMaxRung = (q->uRung > uMaxRung) ? q->uRung : uMaxRung;
       cs += pkdSph(pkd,nnList[i].pPart)->c;
    }
    // We do a simple mean, to avoid computing the kernels again
    cs *= 1./nSmooth;


    if (pkd->param.bAccretion || pkd->param.bBHFeedback){
       BHFIELDS* pBH = pkdBH(pkd,p);
       double pDensity;
       double vRel2;


       // First, we gather all the smoothed quantities
       pDensity = 0.0;
       vel_t vx = 0.0;
       vel_t vy = 0.0;
       vel_t vz = 0.0;
       vel_t *pv = pkdVel(pkd,p);
       float kernelSum = 0.0;
       float massSum = 0.0;
       for (int i=0; i<nSmooth; ++i){
          const double rpq = sqrt(nnList[i].fDist2);
          const double kernel = cubicSplineKernel(rpq, fBall);
          kernelSum += kernel;
          massSum += pkdMass(pkd,nnList[i].pPart);

          pDensity += kernel*pkdMass(pkd,nnList[i].pPart);
          vx += kernel*(pkdVel(pkd,nnList[i].pPart)[0]-pv[0]*inv_a);
          vy += kernel*(pkdVel(pkd,nnList[i].pPart)[1]-pv[1]*inv_a);
          vz += kernel*(pkdVel(pkd,nnList[i].pPart)[2]-pv[2]*inv_a);
       }

       kernelSum = 1./kernelSum;
       vRel2 = vx*vx + vy*vy + vz*vz;
       vRel2 *= kernelSum*kernelSum;
       

       // Do we need to convert to physical?
       pDensity *= inv_a*inv_a*inv_a;
       const double dBondiAccretion = pkd->param.dAccretionAlpha * 4.* M_PI * pBH->dInternalMass*pBH->dInternalMass * pDensity / pow(cs*cs + vRel2, 1.5);

       // All the prefactors are computed at the setup phase
       const double dEddingtonAccretion = pkd->param.dEddingtonFactor * pBH->dInternalMass;

       pBH->dEddingtonRatio = dBondiAccretion/dEddingtonAccretion;

       pBH->dAccretionRate = ( dEddingtonAccretion < dBondiAccretion ) ? dEddingtonAccretion : dBondiAccretion;




       //printf("%d cs %e fBall %e \n", nSmooth, cs, fBall);
       //printf("%e %e %e \t %e \n", dBondiAccretion, dEddingtonAccretion, pDensity, pBH->dInternalMass);
       //assert(0);
       // This will be added to the internal mass in pkdUpdatePrims

       if (pkd->param.bAccretion) { // Swallow particles
          double prob_factor = (pBH->dInternalMass - pMass)/pDensity; 
          if (prob_factor > 0.0){
             double prob;
             for (int i=0; i<nSmooth; ++i){


                const double rpq = sqrt(nnList[i].fDist2);
                const double kernel = cubicSplineKernel(rpq, fBall);
                prob = prob_factor * kernel;
                if (rand()<RAND_MAX*prob) {
                   printf("SWALLOW!\n");
                   PARTICLE *q;
                   if (pkd->idSelf != nnList[i].iPid){
                      // TODO: In order to reduce the number of mdlAcquire, we could place this here, and kept the mdlFetch in the smooth operator
                      //q = mdlAcquire(pkd->mdl,CID_PARTICLE,nnList[i].iIndex,nnList[i].iPid);
                   }
                   q = nnList[i].pPart;

                   float qMass = pkdMass(pkd,q);
                   float newMass = pMass + pkdMass(pkd,q);
                   float inv_newMass = 1./newMass;

                   //printf("Mass: internal %e old %e \t new %e \n", pBH->dInternalMass, pMass, newMass);

                   vel_t *qv = pkdVel(pkd,q);
                   //printf("vx: old %e \n", pv[0] );
                   //printf("vy: old %e \n", pv[1] );
                   //printf("vz: old %e \n", pv[2] );
                   for (int j=0; j<3; j++){
                      // To properly conserve momentum, we need to use the hydrodynamic variable, as 
                      //  the pkdVel may not be updated yet
                      //
                      //  In the case of cosmological simulations, they have different scale factors, so we need to correct for that
                      pv[j] = (pMass*pv[j] + smf->a*pkdSph(pkd,q)->mom[j])*inv_newMass;
                   }
                   

                   //printf("vx: new %e \n", pv[0] );
                   //printf("vy: new %e \n", pv[1] );
                   //printf("vz: new %e \n", pv[2] );

                   float *mass_field = (float *)pkdField(p, pkd->oMass);
                   *mass_field = newMass;

                   //mass_field = (float *)pkdField(nnList[i].pPart, pkd->oMass);
                   //*mass_field = 0.0f; // We do this to check that the deleted particles does not enter the gravity computation, as this
                                      //   will activate an assert
                   if (pkd->idSelf != nnList[i].iPid)
                      pkdSetClass(pkd,0.0,0.0,FIO_SPECIES_LAST,q);
                   else
                      pkdDeleteParticle(pkd, q);

                   //if (pkd->idSelf != nnList[i].iPid) mdlRelease(pkd->mdl,CID_PARTICLE,q);

                   // Once we have one event, we stop checking, as our mass will be higher now
                   break;

                }
             }
          }

       }
       if (pkd->param.bBHFeedback){
          pBH->dFeedbackRate = pkd->param.dBHFeedbackEff * pBH->dAccretionRate;

          double meanMass = massSum/nSmooth;
          double Ecrit = pkd->param.dBHFeedbackEcrit * meanMass;
          if (pBH->dAccEnergy > Ecrit) {
             for (int i=0; i<nSmooth; ++i){
                PARTICLE *q;
                q = nnList[i].pPart;

                // We favour adding the feedback to particles with higher masses... does this make sense?
                // Otherwise, just 1./nSmooth would do just fine
                double prob = 1./nSmooth;//pkdMass(pkd,q) / massSum;
                if (rand()<RAND_MAX*prob) {
                   printf("BH eedback event!\n");
                   SPHFIELDS* qsph = pkdSph(pkd,q);
                   //printf("Uint %e extra %e \n", qsph->Uint, pkd->param.dFeedbackDu * pkdMass(pkd,q));
                   double energy = pkd->param.dBHFeedbackEcrit * pkdMass(pkd,q) ;
                   qsph->Uint += energy;
                   qsph->E += energy;

                   pBH->dAccEnergy -= energy;
                   break;
                }
             }
          }


       }
    }


    // This may happen in very special situations, such as test cases where only
    //  BH particles are placed, or when there is no gravity.
    // Disable it at your own risk!
    //printf("%e %d \n", pkdBall(pkd,p), nSmooth);
    //if (pLowPot==NULL)
    //   for (int i=0; i<nSmooth; ++i)
    //      printf("%e \n", *pkdPot(pkd, nnList[i].pPart));
    p->uNewRung = uMaxRung;
    assert(pLowPot!=NULL);
    if (pLowPot==NULL) return;

    // We only follow exactly that particle if the BH does not 
    // have enough mass to dictate the movement of the particles

    //return; //Needed when running TestGrowth
    if (pMass < 10.*pkdMass(pkd,pLowPot)){
       cs *= cs;
       vel_t *lowPotv = pkdVel(pkd, pLowPot);
       vel_t *pv = pkdVel(pkd,p);
       float v2 = 0.0;
       for (int j=0; j<3; j++)
         v2 += (lowPotv[j]-pv[j]*inv_a) * (lowPotv[j]-pv[j]*inv_a);

       // We set a limit of the velocity to avoid being dragged by fast particles (v<0.25cs)
       //  And dont forget that the velocity have different 'a' factors! \propto a for gas, and \propto a^2 for BH/DM/stars
       //
       //  This is overriden when just creating the BH
       BHFIELDS* pBH = pkdBH(pkd,p);
       if (v2 < 0.0625*cs || pBH->newPos[0]==-1){
          pkdGetPos1(pkd,pLowPot,r_lowPot);
          pkdGetPos1(pkd,p,r_p);
          //printf("Old pos \t %e \t %e \t %e\n", pkdPos(pkd,p,0), pkdPos(pkd,p,1), pkdPos(pkd,p,2));
          for (int j=0;j<3;j++){
             //pkdSetPos( pkd, p, j, pkdPos(pkd,pLowPot, j) );
             pBH->newPos[j] = pkdPos(pkd,pLowPot,j);
             pv[j] = smf->a*lowPotv[j];
          }
          //printf("New pos \t %e \t %e \t %e\n", pBH->newPos[0], pBH->newPos[1], pBH->newPos[2]);
       }
    }

}

void combBHdrift(void *vpkd, void *vp1,void *vp2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE* p1 = (PARTICLE*) vp1;
    PARTICLE* p2 = (PARTICLE*) vp2;

    int pSpecies2 = pkdSpecies(pkd, p2);

    if (pSpecies2 == FIO_SPECIES_LAST){
       pkdDeleteParticle(pkd, p1);
    }else if (pSpecies2 == FIO_SPECIES_SPH){
       SPHFIELDS* psph2 = pkdSph(pkd,p2);

       if (psph2->Uint > 0.0){
          SPHFIELDS* psph1 = pkdSph(pkd,p1);

          psph1->Uint += psph2->Uint;
          psph1->E += psph2->E;

       }

    }



    }

// For creating the COcache, we need a init function, even if it is empty
void initBHdrift(void *vpkd,void *vp){
  PARTICLE* p = (PARTICLE*) vp;
  PKD pkd = (PKD) vpkd;

  // This may be too restrictive?
  //assert(pkdIsGas(pkd,p));
  //

  if (pkdIsGas(pkd,p)){
     SPHFIELDS* psph = pkdSph(pkd,p);

     psph->Uint = 0.0;
     psph->E = 0.0;
  }

}
