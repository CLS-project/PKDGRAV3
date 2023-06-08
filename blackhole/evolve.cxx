#include <numeric>
#include "blackhole/evolve.h"
#include "hydro/hydro.h"
#include "master.h"
using blitz::TinyVector;
using blitz::dot;

void MSR::BHDrift(double dTime, double dDelta) {
    Smooth(dTime,dDelta,SMX_BH_DRIFT,1,param.nSmooth);
    pstRepositionBH(pst, NULL, 0, NULL, 0);

    if (param.bBHAccretion) {
        struct inBHAccretion in;
        in.dScaleFactor = csmTime2Exp(csm,dTime);
        pstBHAccretion(pst, &in, sizeof(in), NULL, 0);
    }
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
                nnList[i].bMarked = 1;
                ++nAccreted;
                printf("SWALLOW!\n");

                auto q = pkd->particles[nnList[i].pPart];
                assert(q.is_gas());

                auto accretor = BHAccretor(pkd,CID_GROUP,nnList[i].iIndex,nnList[i].iPid,true);
                accretor.set_pid(pkd->Self());
                accretor.set_index(&p - pkd->particles.begin());
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
            if (nAccreted > 0 && nnList[i].bMarked) continue; // Skip accreted particles

            if (rand()<RAND_MAX*prob) {
                auto q = pkd->particles[nnList[i].pPart];
                assert(q.is_gas());

                auto &qsph = q.sph();
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
    const float inv_a = 1./smf->a;
    const float ph = 0.5*fBall;

    // We look for the most bounded neighbouring particle
    TinyVector<vel_t,3> mv{0.0};
    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
#ifndef DEBUG_BH_ONLY
        assert(q.is_gas());
#endif
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



        int nAccreted = 0;
        if (smf->bBHAccretion) {
            // We use bMarked to flag accreted neighbours in the BH feedback loop
            nAccreted += std::accumulate(nnList, nnList+nSmooth, 0, [pkd](auto &sum, auto &ii) {
                auto accretor = BHAccretor(pkd,CID_GROUP,ii.iIndex,ii.iPid,false);
                ii.bMarked = accretor.has_accreted();
                return sum + ii.bMarked;
            });
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


void packBHevolve(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhEvolvePack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->velocity = p2.velocity();
        p1->c = p2.sph().c;
        p1->fMass = p2.mass();
#ifdef ENTROPY_SWITCH
        p1->fDensity = p2.density();
#endif
        p1->fPotential = p2.potential();
        p1->uRung = p2.rung();
    }
}

void unpackBHevolve(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhEvolvePack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.velocity() = p2->velocity;
        p1.sph().c = p2->c;
        p1.set_mass(p2->fMass);
#ifdef ENTROPY_SWITCH
        p1.set_density(p2->fDensity);
#endif
        p1.potential() = p2->fPotential;
        p1.set_rung(p2->uRung);
    }
}

void initBHevolve(void *vpkd,void *dst) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(dst)];

    if (p.is_gas()) {
        auto &sph = p.sph();

#ifdef OLD_FB_SCHEME
        sph.E = 0.0;
        sph.Uint = 0.0;
#ifdef ENTROPY_SWITCH
        sph.S = 0.0;
#endif
#else // OLD_FB_SCHEME
        sph.fAccFBEnergy = 0.0;
#endif
    }
}

void flushBHevolve(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhEvolveFlush *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    if (p2.is_gas()) {
        const auto &sph = p2.sph();

#ifdef OLD_FB_SCHEME
        p1->E = sph.E;
        p1->Uint = sph.Uint;
#ifdef ENTROPY_SWITCH
        p1->S = sph.S;
#endif
#else //OLD_FB_SCHEME
        p1->fAccFBEnergy = sph.fAccFBEnergy;
#endif
    }
}

void combBHevolve(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhEvolveFlush *>(src);

    if (p1.is_gas()) {
        auto &sph = p1.sph();

#ifdef OLD_FB_SCHEME
        sph.E += p2->E;
        sph.Uint += p2->Uint;
#ifdef ENTROPY_SWITCH
        sph.S += p2->S;
#endif
#else //OLD_FB_SCHEME
        sph.fAccFBEnergy += p2->fAccFBEnergy;
#endif
    }
}


struct bhAccretionPack {
    uint64_t iOrder;
    uint8_t iClass;
};

struct bhAccretionFlush {
    blitz::TinyVector<double,3> mom;
    uint64_t iOrder;
    float fMass;
};

void packBHAccretion(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhAccretionPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iOrder = p2.order();
    p1->iClass = p2.get_class();
}

void unpackBHAccretion(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhAccretionPack *>(src);

    p1.set_order(p2->iOrder);
    p1.set_class(p2->iClass);
}

void initBHAccretion(void *vpkd,void *dst) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(dst)];

    if (p.is_bh()) {
        p.velocity() = 0.0;
        p.set_mass(0.0);
    }
}

void flushBHAccretion(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhAccretionFlush *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    if (p2.is_bh()) {
        p1->mom = p2.velocity(); // **Momentum** added by the accretion
        p1->iOrder = p2.order();
        p1->fMass = p2.mass();
    }
}

void combBHAccretion(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhAccretionFlush *>(src);

    if (p1.is_bh()) {
        assert(p1.order() == p2->iOrder);
        float old_mass = p1.mass();
        float new_mass = old_mass + p2->fMass;
        float inv_mass = 1./new_mass;

        auto &v1 = p1.velocity();
        v1 = (old_mass*v1 + p2->mom)*inv_mass;

        p1.set_mass(new_mass);
    }
}


void pkdBHAccretion(PKD pkd, double dScaleFactor) {

    mdlPackedCacheCO(pkd->mdl, CID_PARTICLE, NULL, pkd->particles, pkd->Local(),
                     pkd->particles.ParticleSize(), pkd, sizeof(bhAccretionPack),
                     packBHAccretion, unpackBHAccretion, sizeof(bhAccretionFlush),
                     initBHAccretion, flushBHAccretion, combBHAccretion);

    for (auto i = 0; i < pkd->Local(); ++i) {
        auto p = pkd->particles[i];
        if (p.is_gas()) {
            auto accretor = BHAccretor(pkd,i);
            if (accretor.has_accreted()) {
                particleStore::ParticlePointer bh(pkd->particles);

                if (accretor.is_remote()) {
                    bh = &pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,CID_PARTICLE,accretor.get_index(),accretor.get_pid()))];
                }
                else {
                    bh = &pkd->particles[accretor.get_index()];
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
                auto &sph = p.sph();
                if (accretor.is_remote()) {
                    bhv += dScaleFactor * sph.mom;
                }
                else {
                    const float inv_newMass = 1. / bh->mass();
                    bhv = (bhMass*bhv + dScaleFactor*sph.mom) * inv_newMass;
                }

                pkdDeleteParticle(pkd,p);

                if (accretor.is_remote()) mdlRelease(pkd->mdl, CID_PARTICLE, bh);
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

