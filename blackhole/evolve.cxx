#include "blackhole/evolve.h"
#include "hydro/hydro.h"
#include "master.h"
using blitz::TinyVector;
using blitz::dot;

void MSR::BHDrift(double dTime, double dDelta) {
    Smooth(dTime,dDelta,SMX_BH_DRIFT,1,param.nSmooth);
    pstRepositionBH(pst, NULL, 0, NULL, 0);

    struct inBHAccretion in;
    in.dScaleFactor = csmTime2Exp(csm,dTime);
    pstBHAccretion(pst, &in, sizeof(in), NULL, 0);
}


static inline int bhAccretion(PKD pkd, NN *nnList, int nSmooth,
                              particleStore::ParticleReference &p, BHFIELDS &bh,
                              float ph, float pMass, double pDensity,
                              double dScaleFactor) {
    const double prob_factor = (bh.dInternalMass - pMass)/pDensity;
    int nAccreted = 0;
    if (prob_factor > 0.0) {
        for (auto i = 0; i < nSmooth; ++i) {
            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, ph);
            const double prob = prob_factor * kernel;
            if (rand()<RAND_MAX*prob) {
                ++nAccreted;
                printf("SWALLOW!\n");
                auto q = pkd->particles[nnList[i].pPart];
                assert(q.is_gas());
                auto &sph = q.sph();

                sph.BHAccretor.iPid = pkd->Self();
                sph.BHAccretor.iIndex = &p - pkd->particles.begin();
            }
        }
    }
    return nAccreted;
}


static inline void bhFeedback(PKD pkd, NN *nnList, int nSmooth, int nAccreted,
                              particleStore::ParticleReference &p, BHFIELDS &bh,
                              float massSum, double dConstGamma, double dBHFBEff,
                              double dBHFBEcrit) {

    bh.dFeedbackRate = dBHFBEff * bh.dAccretionRate;

    const double meanMass = massSum/nSmooth;
    const double Ecrit = dBHFBEcrit * meanMass;
    if (bh.dAccEnergy > Ecrit) {
        const double nHeat = bh.dAccEnergy / Ecrit;
        const double prob = nHeat / (nSmooth-nAccreted); // Correct probability for accreted particles
        for (auto i = 0; i < nSmooth; ++i) {
            if (rand()<RAND_MAX*prob) {
                auto q = pkd->particles[nnList[i].pPart];
                assert(q.is_gas());

                auto &qsph = q.sph();
                if (qsph.BHAccretor.iPid != NOT_ACCRETED) continue; // Skip accreted particles
                printf("BH feedback event!\n");
                const double dEnergyInput = dBHFBEcrit * q.mass();

#ifdef OLD_FB_SCHEME
                qsph.Uint += dEnergyInput;
                qsph.E += dEnergyInput;
#ifdef ENTROPY_SWITCH
                qsph.S += dEnergyInput*(dConstGamma-1.) *
                          pow(q.density(), -dConstGamma+1);
#endif
#else // OLD_FB_SCHEME
                qsph.fAccFBEnergy += dEnergyInput;
#endif

                bh.dAccEnergy -= dEnergyInput;
            }
        }
    }
}

static inline void bhDrift(PKD pkd, particleStore::ParticleReference &p,
                           float pMass, PARTICLE *pLowPotIn, uint8_t uMaxRung,
                           double maxv2, const TinyVector<vel_t,3> &mv,
                           double dScaleFactor) {
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
    //assert(pLowPotIn!=NULL);
#endif
    if (pLowPotIn==NULL) return;

#ifdef DEBUG_BH_NODRIFT
    return;
#endif

    // We only follow exactly that particle if the BH does not
    // have enough mass to dictate the movement of the particles
    auto LowPot = pkd->particles[pLowPotIn];
    if (pMass < 10.*LowPot.mass()) {
        const double inv_a = 1./dScaleFactor;
        const auto &lowPotv = LowPot.velocity();
        const auto &pv = p.velocity();
        const TinyVector<vel_t,3> v {lowPotv - pv *inv_a - mv};

        // We set a limit of the velocity to avoid being dragged
        //  by fast particles (v<0.25cs).
        // And dont forget that the velocity have different 'a' factors!
        //  \propto a for gas, and \propto a^2 for BH/DM/stars.
        //
        //  This is overriden when just creating the BH
        auto &bh = p.BH();
        if (dot(v,v) < maxv2 || bh.doReposition == 2) {
            bh.doReposition = 1;
            bh.newPos = LowPot.position();
        }
    }
}


void smBHevolve(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    const auto &pv = p.velocity();
    const auto &pMass = p.mass();
    particleStore::ParticlePointer pLowPot(pkd->particles);
    float minPot = HUGE_VAL;
    uint8_t uMaxRung = 0;
    double cs = 0.;
    int nAccreted = 0;
    const float inv_a = 1./smf->a;
    const float ph = 0.5*fBall;

    // We look for the most bounded neighbouring particle
    TinyVector<vel_t,3> mv{0.0};
    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
#ifndef DEBUG_BH_ONLY
        assert(q.is_gas());
#endif
        if (q.sph().BHAccretor.iPid != NOT_ACCRETED) ++nAccreted;
        pLowPot = (q.potential()<minPot)  ? &q : pLowPot;
        // We could have no potential if gravity is not calculated
        minPot = (pLowPot!=NULL) ? pLowPot->potential() : minPot;
        uMaxRung = (q.rung() > uMaxRung) ? q.rung() : uMaxRung;

        cs += q.sph().c;
        mv += q.velocity() - pv * inv_a;
    }
    // We do a simple mean, to avoid computing the kernels again
    const float inv_nSmooth = 1./nSmooth;
    cs *= inv_nSmooth;
    mv *= inv_nSmooth;

    float stdv2 = 0.0;
    for (auto i = 0; i < nSmooth; ++i) {
        const auto &qv = pkd->particles[nnList[i].pPart].velocity();
        const TinyVector<vel_t,3> dv {qv - pv *inv_a - mv};
        stdv2 += dot(dv,dv);
    }
    stdv2 *= inv_nSmooth;

    if (smf->bBHAccretion || smf->bBHFeedback) {
        auto &bh = p.BH();

        // First, we gather all the smoothed quantities
        double pDensity = 0.0;
        double massSum = 0.0;
        TinyVector<vel_t,3> v{0.0}, vcirc{0.0};
        for (auto i = 0; i < nSmooth; ++i) {
            auto q = pkd->particles[nnList[i].pPart];
#ifndef DEBUG_BH_ONLY
            assert(q.is_gas());
#endif
            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, ph);
            const auto &qMass = q.mass();
            massSum += qMass;

            const double weight = qMass * kernel;
            pDensity += weight;

            const auto &qv = q.velocity();
            const TinyVector<vel_t,3> dv {weight *(qv - pv*inv_a)};
            v += dv;
            vcirc += blitz::cross(nnList[i].dr, dv);
        }
        const double norm = 1./pDensity;
        const double vRel2 = dot(v,v) * norm * norm;

        p.set_density(pDensity);

        // We allow the accretion viscosity to act over a boosted Bondi accretion
        double dBondiPrefactor = smf->dBHAccretionAlpha;
        if (smf->dBHAccretionCvisc) {
            const double vphi = sqrt(dot(vcirc,vcirc)) * norm / ph;
            const double fac = cs/vphi;
            dBondiPrefactor *= std::min(1.0, fac*fac*fac/smf->dBHAccretionCvisc);
        }

        // Do we need to convert to physical?
        pDensity *= inv_a*inv_a*inv_a;
        const double dBondiAccretion = dBondiPrefactor *
                                       4.* M_PI * bh.dInternalMass*bh.dInternalMass * pDensity /
                                       pow(cs*cs + vRel2, 1.5);

        // All the prefactors are computed at the setup phase
        const double dEddingtonAccretion = smf->dBHAccretionEddFac *
                                           bh.dInternalMass;

        bh.dEddingtonRatio = dBondiAccretion/dEddingtonAccretion;

        bh.dAccretionRate = ( dEddingtonAccretion < dBondiAccretion ) ?
                            dEddingtonAccretion : dBondiAccretion;



        if (smf->bBHAccretion) {
            nAccreted += bhAccretion(pkd, nnList, nSmooth, p, bh,
                                     ph, pMass, pDensity, smf->a);
        }
        if (smf->bBHFeedback) {
            bhFeedback(pkd, nnList, nSmooth, nAccreted, p, bh, massSum, smf->dConstGamma,
                       smf->dBHFBEff, smf->dBHFBEcrit);
        }
    }

    bhDrift(pkd, p, pMass, pLowPot, uMaxRung, stdv2, mv, smf->a);

}

void combBHevolve(void *vpkd, void *vp1,const void *vp2) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(vp1)];
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(vp2)];

    if (p1.is_gas() && p2.is_gas()) {
        auto &sph1 = p1.sph();
        const auto &sph2 = p2.sph();

#ifdef OLD_FB_SCHEME
        sph1.Uint += sph2.Uint;
        sph1.E += sph2.E;
#ifdef ENTROPY_SWITCH
        sph1.S += sph2.S;
#endif
#else //OLD_FB_SCHEME
        sph1.fAccFBEnergy += sph2.fAccFBEnergy;
#endif

        // Is this the first accretion attempt for this particle?
        if (sph1.BHAccretor.iPid == NOT_ACCRETED &&
                sph2.BHAccretor.iPid != NOT_ACCRETED) {
            sph1.BHAccretor.iPid = sph2.BHAccretor.iPid;
            sph1.BHAccretor.iIndex = sph2.BHAccretor.iIndex;
        }
        // If not, just keep the previous attempt
    }
}


void initBHevolve(void *vpkd,void *vp) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(vp)];

    if (p.is_gas()) {
        auto &sph = p.sph();

#ifdef OLD_FB_SCHEME
        sph.Uint = 0.0;
        sph.E = 0.0;
#ifdef ENTROPY_SWITCH
        sph.S = 0.0;
#endif
#else // OLD_FB_SCHEME
        sph.fAccFBEnergy = 0.0;
#endif

        sph.BHAccretor.iPid = NOT_ACCRETED;
    }
}

void initBHAccretion(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(vp)];
    p.set_mass(0.0);
    p.velocity() = 0.0;
}

void combBHAccretion(void *vpkd, void *vp1, const void *vp2) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(vp1)];
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(vp2)];

    if (p1.is_bh() && p2.is_bh()) {
        assert(p1.order() == p2.order());
        float old_mass = p1.mass();
        float new_mass = old_mass + p2.mass();
        float inv_mass = 1./new_mass;

        auto &v1 = p1.velocity();
        const auto &v2 = p2.velocity(); // **Momentum** added by the accretion
        v1 = (old_mass*v1 + v2)*inv_mass;

        p1.set_mass(new_mass);
    }
}


void pkdBHAccretion(PKD pkd, double dScaleFactor) {

    mdlCOcache(pkd->mdl, CID_PARTICLE, NULL, pkd->particles, pkd->particles.ParticleSize(),
               pkd->Local(), pkd, initBHAccretion, combBHAccretion);

    for (auto &p : pkd->particles) {
        if (p.is_gas()) {
            auto &sph = p.sph();
            if (sph.BHAccretor.iPid != NOT_ACCRETED) {
                particleStore::ParticlePointer bh(pkd->particles);
                // this particle was accreted!
                //printf("%d,%d accreted by %d,%d\n", i, pkd->Self(), sph.BHAccretor.iIndex, sph.BHAccretor.iPid);

                if (sph.BHAccretor.iPid != pkd->Self()) {
                    bh = &pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,CID_PARTICLE,sph.BHAccretor.iIndex,sph.BHAccretor.iPid))];
                }
                else {
                    bh = &pkd->particles[sph.BHAccretor.iIndex];
                }
                assert(bh->is_bh());

                const float bhMass = bh->mass();
                bh->set_mass(bhMass + p.mass());

                auto &bhv = bh->velocity();

                // To properly conserve momentum, we need to use the
                // hydrodynamic variable, as the pkdVel may not be updated yet
                //
                // We have to consider remote and local particles differently,
                // as for the remotes the momentum is accumulated here but then
                // added in the combine function
                if (sph.BHAccretor.iPid != pkd->Self()) {
                    bhv += dScaleFactor * sph.mom;
                }
                else {
                    const float inv_newMass = 1. / bh->mass();
                    bhv = (bhMass*bhv + dScaleFactor*sph.mom) * inv_newMass;
                }

                pkdDeleteParticle(pkd,p);

                if (sph.BHAccretor.iPid != pkd->Self())
                    mdlRelease(pkd->mdl, CID_PARTICLE, bh);
            }
        }
    }

    mdlFinishCache(pkd->mdl,CID_PARTICLE);

}

void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference &p, double dTime,
                    double dDelta, double dBHRadiativeEff) {
    auto &bh = p.BH();
    const double pDelta = dDelta > 0. ? dTime - bh.lastUpdateTime : 0.0;

    bh.dInternalMass += bh.dAccretionRate  * pDelta * (1.-dBHRadiativeEff);
    bh.dAccEnergy += bh.dFeedbackRate * pDelta;
    bh.lastUpdateTime = dTime;
}

