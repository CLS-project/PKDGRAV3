#include "blackhole/evolve.h"
#include "hydro/hydro.h"


static inline void bhAccretion(PKD pkd, NN *nnList, int nSmooth,
                               PARTICLE *p, BHFIELDS *pBH, float fBall, float pMass,
                               double pDensity, double dScaleFactor) {
    double prob_factor = (pBH->dInternalMass - pMass)/pDensity;
    if (prob_factor > 0.0) {
        double prob;
        for (int i=0; i<nSmooth; ++i) {


            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, fBall);
            prob = prob_factor * kernel;
            if (rand()<RAND_MAX*prob) {
                printf("SWALLOW!\n");
                PARTICLE *q;
                if (pkd->idSelf != nnList[i].iPid) {
                    // TODO: In order to reduce the number of mdlAcquire, could we
                    // place this here, and kept the mdlFetch in the smooth?
                    //q = mdlAcquire(pkd->mdl,CID_PARTICLE,
                    //               nnList[i].iIndex,nnList[i].iPid);
                }
                q = nnList[i].pPart;

                float newMass = pMass + pkdMass(pkd,q);
                float inv_newMass = 1./newMass;

                //printf("Mass: internal %e old %e \t new %e \n",
                //         pBH->dInternalMass, pMass, newMass);

                vel_t *pv = pkdVel(pkd,p);
                for (int j=0; j<3; j++) {
                    // To properly conserve momentum, we need to use the
                    // hydrodynamic variable, as the pkdVel may not be updated yet
                    //
                    // In the case of cosmological simulations, they have
                    // different scale factors, so we need to correct for that
                    pv[j] = (pMass*pv[j] + dScaleFactor*pkdSph(pkd,q)->mom[j]) *
                            inv_newMass;
                }


                float *mass_field = (float *)pkdField(p, pkd->oFieldOffset[oMass]);
                *mass_field = newMass;

                if (pkd->idSelf != nnList[i].iPid)
                    pkdSetClass(pkd,0.0,0.0,0,FIO_SPECIES_LAST,q);
                else
                    pkdDeleteParticle(pkd, q);

                //if (pkd->idSelf != nnList[i].iPid)
                //   mdlRelease(pkd->mdl,CID_PARTICLE,q);

                // Once we have one event, we stop checking, as our mass
                //   will be higher now
                break;

            }
        }
    }
}


static inline void bhFeedback(PKD pkd, NN *nnList, int nSmooth, PARTICLE *p,
                              BHFIELDS *pBH, float massSum, double dConstGamma,
                              double dBHFBEff, double dBHFBEcrit) {

    pBH->dFeedbackRate = dBHFBEff * pBH->dAccretionRate;

    const double meanMass = massSum/nSmooth;
    const double Ecrit = dBHFBEcrit * meanMass;
    if (pBH->dAccEnergy > Ecrit) {
        const double nHeat = pBH->dAccEnergy / Ecrit;
        const double prob = nHeat / nSmooth;
        for (int i=0; i<nSmooth; ++i) {
            if (rand()<RAND_MAX*prob) {
                PARTICLE *q;
                q = nnList[i].pPart;

                printf("BH feedback event!\n");
                SPHFIELDS *qsph = pkdSph(pkd,q);
                const double energy = dBHFBEcrit * pkdMass(pkd,q) ;
                qsph->Uint += energy;
                qsph->E += energy;
#ifdef ENTROPY_SWITCH
                qsph->S += energy*(dConstGamma-1.) *
                           pow(pkdDensity(pkd,q), -dConstGamma+1);
#endif

                pBH->dAccEnergy -= energy;
            }
        }
    }
}

static inline void bhDrift(PKD pkd, PARTICLE *p, float pMass,
                           PARTICLE *pLowPot, uint8_t uMaxRung,
                           double cs, double dScaleFactor) {
    // In situations where only BH are present, or there is no gravity, pLowPot
    // may be null. In that case, the drift operation in this function makes no
    // sense.
    //
    // In normal cosmological simulations, the assert may be activated when
    // something has gone really wrong, such as having no gas particle close to
    // a BH. In those cases, keeping the assert is recommended to detect such
    // pathological cases.
    //if (pLowPot==NULL){
    //   printf("%e %d \n", pkdBall(pkd,p), nSmooth);
    //   for (int i=0; i<nSmooth; ++i)
    //      printf("%e \n", *pkdPot(pkd, nnList[i].pPart));
    //}
#ifndef DEBUG_BH_ONLY
    assert(pLowPot!=NULL);
#endif
    p->uNewRung = uMaxRung;
    if (pLowPot==NULL) return;


#ifdef DEBUG_BH_NODRIFT
    return;
#endif
    // We only follow exactly that particle if the BH does not
    // have enough mass to dictate the movement of the particles
    if (pMass < 10.*pkdMass(pkd,pLowPot)) {
        double inv_a = 1./dScaleFactor;
        cs *= cs;
        vel_t *lowPotv = pkdVel(pkd, pLowPot);
        vel_t *pv = pkdVel(pkd,p);
        float v2 = 0.0;
        for (int j=0; j<3; j++)
            v2 += (lowPotv[j]-pv[j]*inv_a) * (lowPotv[j]-pv[j]*inv_a);

        // We set a limit of the velocity to avoid being dragged
        //  by fast particles (v<0.25cs).
        // And dont forget that the velocity have different 'a' factors!
        //  \propto a for gas, and \propto a^2 for BH/DM/stars.
        //
        //  This is overriden when just creating the BH
        BHFIELDS *pBH = pkdBH(pkd,p);
        if (v2 < 0.0625*cs || pBH->newPos[0]==-1) {
            for (int j=0; j<3; j++) {
                pBH->newPos[j] = pkdPos(pkd,pLowPot,j);
                pv[j] = dScaleFactor*lowPotv[j];
            }
        }
    }

}


void smBHevolve(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    PARTICLE *pLowPot = NULL;
    float minPot = HUGE_VAL;
    float pMass = pkdMass(pkd,p);
    uint8_t uMaxRung = 0;
    double cs = 0.;
    float inv_a = 1./smf->a;


    // We look for the most bounded neighbouring particle
    for (int i=0; i<nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;
        assert(pkdIsGas(pkd,q));
        pLowPot = (*pkdPot(pkd,q)<minPot)  ? q : pLowPot;
        // We could have no potential if gravity is not calculated
        minPot = (pLowPot!=NULL) ? *pkdPot(pkd,pLowPot) : minPot;
        uMaxRung = (q->uRung > uMaxRung) ? q->uRung : uMaxRung;
        cs += pkdSph(pkd,nnList[i].pPart)->c;
    }
    // We do a simple mean, to avoid computing the kernels again
    cs *= 1./nSmooth;


    if (smf->bBHAccretion || smf->bBHFeedback) {
        BHFIELDS *pBH = pkdBH(pkd,p);
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
        for (int i=0; i<nSmooth; ++i) {
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
        const double dBondiAccretion = smf->dBHAccretionAlpha *
                                       4.* M_PI * pBH->dInternalMass*pBH->dInternalMass * pDensity /
                                       pow(cs*cs + vRel2, 1.5);

        // All the prefactors are computed at the setup phase
        const double dEddingtonAccretion = smf->dBHAccretionEddFac *
                                           pBH->dInternalMass;

        pBH->dEddingtonRatio = dBondiAccretion/dEddingtonAccretion;

        pBH->dAccretionRate = ( dEddingtonAccretion < dBondiAccretion ) ?
                              dEddingtonAccretion : dBondiAccretion;




        //printf("%d cs %e fBall %e \n", nSmooth, cs, fBall);
        //printf("%e %e %e \t %e \n",
        //  dBondiAccretion, dEddingtonAccretion, pDensity, pBH->dInternalMass);
        //assert(0);

        if (smf->bBHAccretion) {
            bhAccretion(pkd, nnList, nSmooth, p, pBH,
                        fBall, pMass, pDensity, smf->a);
        }
        if (smf->bBHFeedback) {
            bhFeedback(pkd, nnList, nSmooth, p, pBH, massSum, smf->dConstGamma,
                       smf->dBHFBEff, smf->dBHFBEcrit);
        }
    }

    bhDrift(pkd, p, pMass, pLowPot, uMaxRung, cs, smf->a);


}

void combBHevolve(void *vpkd, void *vp1,const void *vp2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) vp1;
    PARTICLE *p2 = (PARTICLE *) vp2;

    int pSpecies2 = pkdSpecies(pkd, p2);

    if (pSpecies2 == FIO_SPECIES_LAST) {
        pkdDeleteParticle(pkd, p1);
    }
    else if (pSpecies2 == FIO_SPECIES_SPH) {
        SPHFIELDS *psph2 = pkdSph(pkd,p2);

        if (psph2->Uint > 0.0) {
            SPHFIELDS *psph1 = pkdSph(pkd,p1);

            psph1->Uint += psph2->Uint;
            psph1->E += psph2->E;
#ifdef ENTROPY_SWITCH
            psph1->S += psph2->S;
#endif

        }

    }
}


void initBHevolve(void *vpkd,void *vp) {
    PARTICLE *p = (PARTICLE *) vp;
    PKD pkd = (PKD) vpkd;

    if (pkdIsGas(pkd,p)) {
        SPHFIELDS *psph = pkdSph(pkd,p);

        psph->Uint = 0.0;
        psph->E = 0.0;
#ifdef ENTROPY_SWITCH
        psph->S = 0.;
#endif
    }

}


void pkdBHIntegrate(PKD pkd, PARTICLE *p, double dTime, double dDelta, double dBHRadiativeEff) {
    BHFIELDS *pBH = pkdBH(pkd,p);

    double pDelta = 0.0;
    if (dDelta > 0)
        pDelta = dTime - pBH->lastUpdateTime;

    pBH->dInternalMass += pBH->dAccretionRate  * pDelta * (1.-dBHRadiativeEff);
    pBH->dAccEnergy += pBH->dFeedbackRate * pDelta;
    pBH->lastUpdateTime = dTime;

}
