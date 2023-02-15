#include "blackhole/evolve.h"
#include "hydro/hydro.h"


static inline void bhAccretion(PKD pkd, NN *nnList, int nSmooth,
                               particleStore::ParticleReference &p, BHFIELDS &BH, float fBall, float pMass,
                               double pDensity, double dScaleFactor) {
    double prob_factor = (BH.dInternalMass - pMass)/pDensity;
    if (prob_factor > 0.0) {
        double prob;
        for (int i=0; i<nSmooth; ++i) {


            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, fBall);
            prob = prob_factor * kernel;
            if (rand()<RAND_MAX*prob) {
                printf("SWALLOW!\n");
                if (pkd->Self() != nnList[i].iPid) {
                    // TODO: In order to reduce the number of mdlAcquire, could we
                    // place this here, and kept the mdlFetch in the smooth?
                    //q = mdlAcquire(pkd->mdl,CID_PARTICLE,
                    //               nnList[i].iIndex,nnList[i].iPid);
                }
                auto q = pkd->particles[nnList[i].pPart];

                float newMass = pMass + q.mass();
                float inv_newMass = 1./newMass;

                //printf("Mass: internal %e old %e \t new %e \n",
                //         BH.dInternalMass, pMass, newMass);

                auto pv = p.velocity();
                for (int j=0; j<3; j++) {
                    // To properly conserve momentum, we need to use the
                    // hydrodynamic variable, as the pkdVel may not be updated yet
                    //
                    // In the case of cosmological simulations, they have
                    // different scale factors, so we need to correct for that
                    pv[j] = (pMass*pv[j] + dScaleFactor*q.sph().mom[j]) *
                            inv_newMass;
                }

                p.set_mass(newMass);

                if (pkd->Self() != nnList[i].iPid)
                    q.set_class(0.0,0.0,0,FIO_SPECIES_UNKNOWN);
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


static inline void bhFeedback(PKD pkd, NN *nnList, int nSmooth, particleStore::ParticleReference &p,
                              BHFIELDS &BH, float massSum, double dConstGamma,
                              double dBHFBEff, double dBHFBEcrit) {
    // auto p = pkd->particles[pIn];
    BH.dFeedbackRate = dBHFBEff * BH.dAccretionRate;

    const double meanMass = massSum/nSmooth;
    const double Ecrit = dBHFBEcrit * meanMass;
    if (BH.dAccEnergy > Ecrit) {
        const double nHeat = BH.dAccEnergy / Ecrit;
        const double prob = nHeat / nSmooth;
        for (int i=0; i<nSmooth; ++i) {
            if (rand()<RAND_MAX*prob) {
                auto q = pkd->particles[nnList[i].pPart];

                printf("BH feedback event!\n");
                auto &qsph = q.sph();
                const double energy = dBHFBEcrit * q.mass();
                qsph.Uint += energy;
                qsph.E += energy;
#ifdef ENTROPY_SWITCH
                qsph.S += energy*(dConstGamma-1.) *
                          pow(q->density(), -dConstGamma+1);
#endif

                BH.dAccEnergy -= energy;
            }
        }
    }
}

static inline void bhDrift(PKD pkd, particleStore::ParticleReference &p, float pMass,
                           PARTICLE *pLowPotIn, uint8_t uMaxRung,
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
    //   printf("%e %d \n", p->ball(), nSmooth);
    //   for (int i=0; i<nSmooth; ++i)
    //      printf("%e \n", *pkdPot(pkd, nnList[i].pPart));
    //}
#ifndef DEBUG_BH_ONLY
    assert(pLowPotIn!=NULL);
#endif
    p.set_new_rung(uMaxRung);
    if (pLowPotIn==NULL) return;
    auto LowPot = pkd->particles[pLowPotIn];

#ifdef DEBUG_BH_NODRIFT
    return;
#endif
    // We only follow exactly that particle if the BH does not
    // have enough mass to dictate the movement of the particles
    if (pMass < 10.*LowPot.mass()) {
        double inv_a = 1./dScaleFactor;
        cs *= cs;
        auto lowPotv = LowPot.velocity();
        auto &pv = p.velocity();
        float v2 = 0.0;
        for (int j=0; j<3; j++)
            v2 += (lowPotv[j]-pv[j]*inv_a) * (lowPotv[j]-pv[j]*inv_a);

        // We set a limit of the velocity to avoid being dragged
        //  by fast particles (v<0.25cs).
        // And dont forget that the velocity have different 'a' factors!
        //  \propto a for gas, and \propto a^2 for BH/DM/stars.
        //
        //  This is overriden when just creating the BH
        auto &BH = p.BH();
        if (v2 < 0.0625*cs || BH.newPos[0]==-1) {
            BH.newPos = LowPot.position();
            pv = dScaleFactor*lowPotv;
        }
    }

}


void smBHevolve(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    particleStore::ParticlePointer pLowPot(pkd->particles);
    float minPot = HUGE_VAL;
    float pMass = p.mass();
    uint8_t uMaxRung = 0;
    double cs = 0.;
    float inv_a = 1./smf->a;


    // We look for the most bounded neighbouring particle
    for (int i=0; i<nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        assert(q.is_gas());
        pLowPot = (q.potential()<minPot)  ? &q : pLowPot;
        // We could have no potential if gravity is not calculated
        minPot = (pLowPot!=nullptr) ? pLowPot->potential() : minPot;
        uMaxRung = (q.rung() > uMaxRung) ? q.rung() : uMaxRung;
        cs += q.sph().c;
    }
    // We do a simple mean, to avoid computing the kernels again
    cs *= 1./nSmooth;


    if (smf->bBHAccretion || smf->bBHFeedback) {
        auto &BH = p.BH();
        double pDensity;
        double vRel2;


        // First, we gather all the smoothed quantities
        pDensity = 0.0;
        blitz::TinyVector<vel_t,3> v(0.0);
        auto &pv = p.velocity();
        float kernelSum = 0.0;
        float massSum = 0.0;
        for (int i=0; i<nSmooth; ++i) {
            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, fBall);
            kernelSum += kernel;
            massSum += pkdMass(pkd,nnList[i].pPart);

            pDensity += kernel*pkdMass(pkd,nnList[i].pPart);
            v += kernel*(pkd->particles[nnList[i].pPart].velocity()-pv*inv_a);
        }

        kernelSum = 1./kernelSum;
        vRel2 = blitz::dot(v,v);
        vRel2 *= kernelSum*kernelSum;


        // Do we need to convert to physical?
        pDensity *= inv_a*inv_a*inv_a;
        const double dBondiAccretion = smf->dBHAccretionAlpha *
                                       4.* M_PI * BH.dInternalMass*BH.dInternalMass * pDensity /
                                       pow(cs*cs + vRel2, 1.5);

        // All the prefactors are computed at the setup phase
        const double dEddingtonAccretion = smf->dBHAccretionEddFac *
                                           BH.dInternalMass;

        BH.dEddingtonRatio = dBondiAccretion/dEddingtonAccretion;

        BH.dAccretionRate = ( dEddingtonAccretion < dBondiAccretion ) ?
                            dEddingtonAccretion : dBondiAccretion;




        //printf("%d cs %e fBall %e \n", nSmooth, cs, fBall);
        //printf("%e %e %e \t %e \n",
        //  dBondiAccretion, dEddingtonAccretion, pDensity, BH.dInternalMass);
        //assert(0);

        if (smf->bBHAccretion) {
            bhAccretion(pkd, nnList, nSmooth, p, BH,
                        fBall, pMass, pDensity, smf->a);
        }
        if (smf->bBHFeedback) {
            bhFeedback(pkd, nnList, nSmooth, p, BH, massSum, smf->dConstGamma,
                       smf->dBHFBEff, smf->dBHFBEcrit);
        }
    }

    bhDrift(pkd, p, pMass, pLowPot, uMaxRung, cs, smf->a);


}

void combBHevolve(void *vpkd, void *vp1,const void *vp2) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(vp1)];
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(vp2)];

    int pSpecies2 = p2.species();

    if (pSpecies2 == FIO_SPECIES_UNKNOWN) {
        pkdDeleteParticle(pkd, p1);
    }
    else if (pSpecies2 == FIO_SPECIES_SPH) {
        const auto &sph2 = p2.sph();

        if (sph2.Uint > 0.0) {
            auto &sph1 = p1.sph();

            sph1.Uint += sph2.Uint;
            sph1.E += sph2.E;
#ifdef ENTROPY_SWITCH
            sph1.S += sph2.S;
#endif

        }

    }
}


void initBHevolve(void *vpkd,void *vp) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(vp)];

    if (p.is_gas()) {
        auto &sph = p.sph();

        sph.Uint = 0.0;
        sph.E = 0.0;
#ifdef ENTROPY_SWITCH
        sph.S = 0.;
#endif
    }

}


void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference &p, double dTime, double dDelta, double dBHRadiativeEff) {
    auto &BH = p.BH();

    double pDelta = 0.0;
    if (dDelta > 0)
        pDelta = dTime - BH.lastUpdateTime;

    BH.dInternalMass += BH.dAccretionRate  * pDelta * (1.-dBHRadiativeEff);
    BH.dAccEnergy += BH.dFeedbackRate * pDelta;
    BH.lastUpdateTime = dTime;

}
