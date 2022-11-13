#include "blackhole/evolve.h"
#include "hydro/hydro.h"
#include "master.h"

void MSR::BHDrift(double dTime, double dDelta) {
    Smooth(dTime,dDelta,SMX_BH_DRIFT,1,param.nSmooth);
    pstRepositionBH(pst, NULL, 0, NULL, 0);

    struct inBHAccretion in;
    in.dScaleFactor = csmTime2Exp(csm,dTime);
    pstBHAccretion(pst, &in, sizeof(in), NULL, 0);
}



#ifdef __cplusplus
extern "C" {
#endif

static inline int bhAccretion(PKD pkd, NN *nnList, int nSmooth,
                              PARTICLE *p, BHFIELDS *pBH, float ph, float pMass,
                              double pDensity, double dScaleFactor) {
    double prob_factor = (pBH->dInternalMass - pMass)/pDensity;
    int naccreted = 0;
    if (prob_factor > 0.0) {
        double prob;
        for (int i=0; i<nSmooth; ++i) {


            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, ph);
            prob = prob_factor * kernel;
            if (rand()<RAND_MAX*prob) {
                naccreted++;
                printf("SWALLOW!\n");
                PARTICLE *q = nnList[i].pPart;
                assert(pkdIsGas(pkd,q));

                pkdSph(pkd,q)->BHAccretor.iPid = pkd->idSelf;
                pkdSph(pkd,q)->BHAccretor.iIndex = pkdParticleIndex(pkd,p);
            }
        }
    }
    return naccreted;
}


static inline void bhFeedback(PKD pkd, NN *nnList, int nSmooth, int naccreted, PARTICLE *p,
                              BHFIELDS *pBH, float massSum, double dConstGamma,
                              double dBHFBEff, double dBHFBEcrit) {

    pBH->dFeedbackRate = dBHFBEff * pBH->dAccretionRate;

    const double meanMass = massSum/nSmooth;
    const double Ecrit = dBHFBEcrit * meanMass;
    if (pBH->dAccEnergy > Ecrit) {
        const double nHeat = pBH->dAccEnergy / Ecrit;
        const double prob = nHeat / (nSmooth-naccreted); // Correct probability for accreted particles
        for (int i=0; i<nSmooth; ++i) {
            if (rand()<RAND_MAX*prob) {
                PARTICLE *q;
                q = nnList[i].pPart;
                assert(pkdIsGas(pkd,q));

                SPHFIELDS *qsph = pkdSph(pkd,q);
                if (qsph->BHAccretor.iPid != NOT_ACCRETED) continue; // Skip accreted particles
                printf("BH feedback event!\n");
                const double dEnergyInput = dBHFBEcrit * pkdMass(pkd,q) ;

#ifdef OLD_FB_SCHEME
                qsph->Uint += dEnergyInput;
                qsph->E += dEnergyInput;
#ifdef ENTROPY_SWITCH
                qsph->S += dEnergyInput*(dConstGamma-1.) *
                           pow(pkdDensity(pkd,q), -dConstGamma+1);
#endif
#else // OLD_FB_SCHEME
                qsph->fAccFBEnergy += dEnergyInput;
#endif

                pBH->dAccEnergy -= dEnergyInput;
            }
        }
    }
}

static inline void bhDrift(PKD pkd, PARTICLE *p, float pMass,
                           PARTICLE *pLowPot, uint8_t uMaxRung,
                           double maxv2, float mv[3], double dScaleFactor) {
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
    //assert(pLowPot!=NULL);
#endif
    if (pLowPot==NULL) return;


#ifdef DEBUG_BH_NODRIFT
    return;
#endif
    // We only follow exactly that particle if the BH does not
    // have enough mass to dictate the movement of the particles
    if (pMass < 10.*pkdMass(pkd,pLowPot)) {
        double inv_a = 1./dScaleFactor;
        vel_t *lowPotv = pkdVel(pkd, pLowPot);
        vel_t *pv = pkdVel(pkd,p);
        float v2 = 0.0;
        for (int j=0; j<3; j++)
            v2 += (lowPotv[j]-pv[j]*inv_a - mv[j]) * (lowPotv[j]-pv[j]*inv_a -mv[j]);

        // We set a limit of the velocity to avoid being dragged
        //  by fast particles (v<0.25cs).
        // And dont forget that the velocity have different 'a' factors!
        //  \propto a for gas, and \propto a^2 for BH/DM/stars.
        //
        //  This is overriden when just creating the BH
        BHFIELDS *pBH = pkdBH(pkd,p);
        //printf("check %e %e \n", v2, cs);
        if (v2 < maxv2 || pBH->doReposition == 2) {
            pBH->doReposition = 1;
            for (int j=0; j<3; j++) {
                pBH->newPos[j] = pkdPos(pkd,pLowPot,j);
                //pv[j] = dScaleFactor*lowPotv[j];
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
    int naccreted = 0;
    float ph = 0.5*fBall;

    // We look for the most bounded neighbouring particle
    vel_t mvx = 0.;
    vel_t mvy = 0.;
    vel_t mvz = 0.;
    for (int i=0; i<nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;
#ifndef DEBUG_BH_ONLY
        assert(pkdIsGas(pkd,q));
#endif
        if (pkdSph(pkd,q)->BHAccretor.iPid != NOT_ACCRETED) naccreted++;
        pLowPot = (*pkdPot(pkd,q)<minPot)  ? q : pLowPot;
        // We could have no potential if gravity is not calculated
        minPot = (pLowPot!=NULL) ? *pkdPot(pkd,pLowPot) : minPot;
        uMaxRung = (q->uRung > uMaxRung) ? q->uRung : uMaxRung;
        cs += pkdSph(pkd,nnList[i].pPart)->c;

        vel_t *pv = pkdVel(pkd,p);
        mvx += pkdVel(pkd,nnList[i].pPart)[0]-pv[0]*inv_a;
        mvy += pkdVel(pkd,nnList[i].pPart)[1]-pv[1]*inv_a;
        mvz += pkdVel(pkd,nnList[i].pPart)[2]-pv[2]*inv_a;
    }
    // We do a simple mean, to avoid computing the kernels again
    const float inv_nSmooth = 1./nSmooth;
    cs *= inv_nSmooth;

    mvx *= inv_nSmooth;
    mvy *= inv_nSmooth;
    mvz *= inv_nSmooth;


    float stdvx = 0.0;
    float stdvy = 0.0;
    float stdvz = 0.0;
    for (int i=0; i<nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;

        vel_t *pv = pkdVel(pkd,p);
        const float vx = pkdVel(pkd,nnList[i].pPart)[0]-pv[0]*inv_a;
        const float vy = pkdVel(pkd,nnList[i].pPart)[1]-pv[1]*inv_a;
        const float vz = pkdVel(pkd,nnList[i].pPart)[2]-pv[2]*inv_a;
        stdvx += (vx-mvx)*(vx-mvx);
        stdvy += (vy-mvy)*(vy-mvy);
        stdvz += (vz-mvz)*(vz-mvz);
    }
    stdvx *= inv_nSmooth;
    stdvy *= inv_nSmooth;
    stdvz *= inv_nSmooth;
    float stdv2 = stdvx + stdvy + stdvz;

    //printf("vx\t %e %e\n", mvx, sqrt(stdvx));
    //printf("vy\t %e %e\n", mvy, sqrt(stdvy));
    //printf("vz\t %e %e\n", mvz, sqrt(stdvz));

    if (smf->bBHAccretion || smf->bBHFeedback) {
        BHFIELDS *pBH = pkdBH(pkd,p);
        double pDensity;
        double vRel2;


        // First, we gather all the smoothed quantities
        pDensity = 0.0;
        vel_t vx = 0.0;
        vel_t vy = 0.0;
        vel_t vz = 0.0;
        double vcircx = 0.0;
        double vcircy = 0.0;
        double vcircz = 0.0;
        vel_t *pv = pkdVel(pkd,p);
        double massSum = 0.0;
        double kernelSum = 0.0;
        for (int i=0; i<nSmooth; ++i) {
#ifndef DEBUG_BH_ONLY
            assert(pkdIsGas(pkd,nnList[i].pPart));
#endif
            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, ph);
            const double  qmass = pkdMass(pkd,nnList[i].pPart);
            kernelSum += kernel;
            massSum += qmass;

            const double weight = qmass * kernel;
            pDensity += weight;

            const double dvx = weight*(pkdVel(pkd,nnList[i].pPart)[0]-pv[0]*inv_a);
            const double dvy = weight*(pkdVel(pkd,nnList[i].pPart)[1]-pv[1]*inv_a);
            const double dvz = weight*(pkdVel(pkd,nnList[i].pPart)[2]-pv[2]*inv_a);

            vx += dvx;
            vy += dvy;
            vz += dvz;

            const double dvcircx = (nnList[i].dy*dvz - nnList[i].dz*dvy);
            const double dvcircy = (nnList[i].dz*dvx - nnList[i].dx*dvz);
            const double dvcircz = (nnList[i].dx*dvy - nnList[i].dy*dvx);
            vcircx += dvcircx;
            vcircy += dvcircy;
            vcircz += dvcircz;
        }
        const double norm = 1./pDensity;
        vRel2 = vx*vx + vy*vy + vz*vz;
        vRel2 *= norm*norm;

        const double vphi = sqrt(vcircx*vcircx + vcircy*vcircy + vcircz*vcircz)*norm/ph;

        // We allow the accretion viscosity to act over a boosted Bondi accretion
        double dBondiPrefactor = smf->dBHAccretionAlpha;
        if (smf->dBHAccretionCvisc) {
            const double fac = cs/vphi;
            dBondiPrefactor *= MIN(1.0, fac*fac*fac/smf->dBHAccretionCvisc);
        }

        // Do we need to convert to physical?
        pDensity *= inv_a*inv_a*inv_a;
        const double dBondiAccretion = dBondiPrefactor *
                                       4.* M_PI * pBH->dInternalMass*pBH->dInternalMass * pDensity /
                                       pow(cs*cs + vRel2, 1.5);

        // All the prefactors are computed at the setup phase
        const double dEddingtonAccretion = smf->dBHAccretionEddFac *
                                           pBH->dInternalMass;

        pBH->dEddingtonRatio = dBondiAccretion/dEddingtonAccretion;

        pBH->dAccretionRate = ( dEddingtonAccretion < dBondiAccretion ) ?
                              dEddingtonAccretion : dBondiAccretion;



        if (smf->bBHAccretion) {
            naccreted += bhAccretion(pkd, nnList, nSmooth, p, pBH,
                                     ph, pMass, pDensity, smf->a);
        }
        if (smf->bBHFeedback) {
            bhFeedback(pkd, nnList, nSmooth, naccreted, p, pBH, massSum, smf->dConstGamma,
                       smf->dBHFBEff, smf->dBHFBEcrit);
        }
    }

    /*
    pLowPot = NULL;
    minPot  = HUGE_VAL;
    float cs2 = cs*cs;
    for (int i=0; i<nSmooth; ++i) {
        PARTICLE *q = nnList[i].pPart;
        vel_t *lowPotv = pkdVel(pkd, q);
        vel_t *pv = pkdVel(pkd,p);
        float v2 = 0.0;
        for (int j=0; j<3; j++)
            v2 += (lowPotv[j]-pv[j]*inv_a) * (lowPotv[j]-pv[j]*inv_a);
        pLowPot = (*pkdPot(pkd,q)<minPot && v2 < 0.0625*cs2)  ? q : pLowPot;
        minPot = (pLowPot!=NULL) ? *pkdPot(pkd,pLowPot) : minPot;
    }
    */
    float mv[3] = {mvx, mvy, mvz};
    bhDrift(pkd, p, pMass, pLowPot, uMaxRung, stdv2, mv, smf->a);


}

void combBHevolve(void *vpkd, void *vp1,const void *vp2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) vp1;
    PARTICLE *p2 = (PARTICLE *) vp2;

    int pSpecies2 = pkdSpecies(pkd, p2);

    if (pSpecies2 == FIO_SPECIES_SPH) {
        SPHFIELDS *psph2 = pkdSph(pkd,p2);

        if (psph2->Uint > 0.0) {
            SPHFIELDS *psph1 = pkdSph(pkd,p1);

            psph1->Uint += psph2->Uint;
            psph1->E += psph2->E;
#ifdef ENTROPY_SWITCH
            psph1->S += psph2->S;
#endif

        }

        if (psph2->BHAccretor.iPid != NOT_ACCRETED) {
            SPHFIELDS *psph1 = pkdSph(pkd,p1);
            if (psph1->BHAccretor.iPid != NOT_ACCRETED) {
                // First try to accrete this particle
                psph1->BHAccretor.iPid   = psph2->BHAccretor.iPid;
                psph1->BHAccretor.iIndex = psph2->BHAccretor.iIndex;
            }// Otherwise just keep the previous attempt
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
        psph->BHAccretor.iPid = NOT_ACCRETED;
    }

}

void initBHAccretion(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = (PARTICLE *) vp;
    float *mass = (float *)pkdField(p, pkd->oFieldOffset[oMass]);
    *mass = 0.0;
    vel_t *v = pkdVel(pkd,p);
    for (int i=0; i<3; i++)
        v[i] = 0.;
}

void combBHAccretion(void *vpkd, void *vp1, const void *vp2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) vp1;
    PARTICLE *p2 = (PARTICLE *) vp2;

    if (pkdIsBH(pkd,p1) && pkdIsBH(pkd,p2)) {
        assert(*pkdParticleID(pkd,p1) == *pkdParticleID(pkd,p2));
        float old_mass = pkdMass(pkd,p1);
        float new_mass = old_mass + pkdMass(pkd, p2);
        float inv_mass = 1./new_mass;

        vel_t *v1 = pkdVel(pkd, p1);
        vel_t *v2 = pkdVel(pkd, p2); // **Momentum** added by the accretion
        for (int i=0; i<3; i++)
            v1[i] = (old_mass*v1[i] + v2[i])*inv_mass;

        float *mass1 = (float *)pkdField(p1, pkd->oFieldOffset[oMass]);
        *mass1 = new_mass;
    }


}


void pkdBHAccretion(PKD pkd, double dScaleFactor) {

    mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),
               pkdLocal(pkd),pkd,initBHAccretion,combBHAccretion);

    for (int i=0; i<pkd->nLocal; i++) {
        PARTICLE *p = pkdParticle(pkd,i);
        if (pkdIsGas(pkd,p)) {
            SPHFIELDS *psph = pkdSph(pkd, p);
            if (psph->BHAccretor.iPid != NOT_ACCRETED) {
                PARTICLE *bh;
                // this particle was accreted!
                //printf("%d,%d accreted by %d,%d\n", i, pkd->idSelf, psph->BHAccretor.iIndex, psph->BHAccretor.iPid);

                if (psph->BHAccretor.iPid != pkd->idSelf) {
                    bh = (PARTICLE * )mdlAcquire(pkd->mdl,CID_PARTICLE,psph->BHAccretor.iIndex,psph->BHAccretor.iPid);
                }
                else {
                    bh = pkdParticle(pkd, psph->BHAccretor.iIndex);
                }
                assert(pkdIsBH(pkd,bh));
                float bhMass = pkdMass(pkd,bh);
                float newMass = bhMass + pkdMass(pkd,p);
                float inv_newMass = 1./newMass;

                //printf("Mass: internal %e old %e \t new %e \n",
                //         pBH->dInternalMass, pMass, newMass);

                vel_t *pv = pkdVel(pkd,bh);

                // To properly conserve momentum, we need to use the
                // hydrodynamic variable, as the pkdVel may not be updated yet
                //
                // We have to consider remote and local particles differently,
                // as for the remotes the momentum is accumulated here but then
                // added in the combine function
                if (psph->BHAccretor.iPid != pkd->idSelf)
                    for (int j=0; j<3; j++)
                        pv[j] = dScaleFactor*psph->mom[j];
                else
                    for (int j=0; j<3; j++)
                        pv[j] = (bhMass*pv[j] + dScaleFactor*psph->mom[j]) *
                                inv_newMass;


                float *mass_field = (float *)pkdField(bh, pkd->oFieldOffset[oMass]);
                *mass_field = newMass;

                pkdDeleteParticle(pkd,p);

                if (psph->BHAccretor.iPid != pkd->idSelf)
                    mdlRelease(pkd->mdl, CID_PARTICLE, bh);

            }

        }

    }
    mdlFinishCache(pkd->mdl,CID_PARTICLE);

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
#ifdef __cplusplus
}
#endif
