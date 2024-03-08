#include "blackhole/evolve.h"
#include "hydro/hydro.h"
#include "master.h"
using blitz::TinyVector;
using blitz::dot;

void MSR::BHEvolve(double dTime, double dDelta) {
    Smooth(dTime,dDelta,SMX_BH_EVOLVE,1,parameters.get_nSmooth());
}


static inline int bhAccretion(PKD pkd, NN *nnList, int nSmooth,
                              particleStore::ParticleReference &p, BHFIELDS &bh,
                              float pH, float pMass, double pDensity,
                              double dScaleFactor) {
    const double prob_factor = (bh.dInternalMass - pMass)/pDensity;
    int nAccreted = 0;
    if (prob_factor > 0.0) {
        for (auto i = 0; i < nSmooth; ++i) {
            auto q = pkd->particles[nnList[i].pPart];
            assert(q.is_gas());
            auto &sph = q.sph();
            if (sph.BHAccretor.iPid != NOT_ACCRETED) continue; // Skip accreted particles
            const double rpq = sqrt(nnList[i].fDist2);
            const double kernel = cubicSplineKernel(rpq, static_cast<double>(pH));
            const double prob = prob_factor * kernel;
            if (rand()<RAND_MAX*prob) {
                ++nAccreted;
                printf("SWALLOW!\n");

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


void smBHEvolve(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    const auto &pv = p.velocity();
    const auto &pMass = p.mass();
    const float inv_a = 1./smf->a;
    // CDV - pH = H
    const float ph = 0.5*fBall;
    const float pH = fBall;

    auto &bh = p.BH();

    // First, we gather all the smoothed quantities
    int nAccreted = 0;
    double pDensity = 0.0;
    double massSum = 0.0;
    double cs = 0.0;
    TinyVector<vel_t,3> v{0.0}, vcirc{0.0};
    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        assert(q.is_gas());
        if (q.sph().BHAccretor.iPid != NOT_ACCRETED) ++nAccreted;
        const double rpq = sqrt(nnList[i].fDist2);
        const double kernel = cubicSplineKernel(rpq, static_cast<double>(pH));
        const auto &qMass = q.mass();
        massSum += qMass;

        const double weight = qMass * kernel;
        pDensity += weight;
        cs += weight * q.sph().c;

        const auto &qv = q.velocity();
        const TinyVector<vel_t,3> dv {weight *(qv - pv*inv_a)};
        v += dv;
        vcirc += blitz::cross(nnList[i].dr, dv);
    }
    const double norm = 1./pDensity;
    const double vRel2 = dot(v,v) * norm * norm;
    cs *= norm;

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
                                 pH, pMass, pDensity, smf->a);
    }
    if (smf->bBHFeedback) {
        bhFeedback(pkd, nnList, nSmooth, nAccreted, p, bh, massSum, smf->dConstGamma,
                   smf->dBHFBEff, smf->dBHFBEcrit);
    }
}


void packBHEvolve(void *vpkd,void *dst,const void *src) {
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
    }
}

void unpackBHEvolve(void *vpkd,void *dst,const void *src) {
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
    }
}

void initBHEvolve(void *vpkd,void *dst) {
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

        sph.BHAccretor.iPid = NOT_ACCRETED;
    }
}

void flushBHEvolve(void *vpkd,void *dst,const void *src) {
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

        p1->iAccPid = sph.BHAccretor.iPid;
        p1->iAccIndex = sph.BHAccretor.iIndex;
    }
}

void combBHEvolve(void *vpkd,void *dst,const void *src) {
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

        // Is this the first accretion attempt for this particle?
        if (sph.BHAccretor.iPid == NOT_ACCRETED &&
                p2->iAccPid != NOT_ACCRETED) {
            sph.BHAccretor.iPid = p2->iAccPid;
            sph.BHAccretor.iIndex = p2->iAccIndex;
        }
        // If not, just keep the previous attempt
    }
}


void pkdBHIntegrate(PKD pkd, particleStore::ParticleReference &p, double dTime,
                    double dDelta, double dBHRadiativeEff) {
    auto &bh = p.BH();
    const double pDelta = dDelta > 0. ? dTime - bh.lastUpdateTime : 0.0;

    bh.dInternalMass += bh.dAccretionRate  * pDelta * (1.-dBHRadiativeEff);
    bh.dAccEnergy += bh.dFeedbackRate * pDelta;
    bh.lastUpdateTime = dTime;
}

