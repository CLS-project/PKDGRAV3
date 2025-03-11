#include "master.h"
#include "blackhole/evolve.h"

void MSR::BHAccretion(double dTime) {
    struct inBHAccretion in;
    in.dScaleFactor = csmTime2Exp(csm,dTime);
    pstBHAccretion(pst, &in, sizeof(in), NULL, 0);
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
    for (auto &p : pkd->particles) {
        if (p.is_gas()) {
            auto &sph = p.sph();
            auto &iPid = sph.BHAccretor.iPid;
            auto &iIndex = sph.BHAccretor.iIndex;
            if (iPid != NOT_ACCRETED) {
                particleStore::ParticlePointer bh(pkd->particles);
                if (iPid != pkd->Self()) {
                    bh = &pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,CID_PARTICLE,iIndex,iPid))];
                }
                else {
                    bh = &pkd->particles[iIndex];
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
                if (iPid != pkd->Self()) {
                    bhv += dScaleFactor * sph.mom;
                }
                else {
                    const float inv_newMass = 1. / bh->mass();
                    bhv = (bhMass*bhv + dScaleFactor*sph.mom) * inv_newMass;
                }

                pkdDeleteParticle(pkd,p);

                if (iPid != pkd->Self())
                    mdlRelease(pkd->mdl, CID_PARTICLE, bh);
            }
        }
    }
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
}

