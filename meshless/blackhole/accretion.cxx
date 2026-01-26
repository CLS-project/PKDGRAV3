#include "master.h"
#include "blackhole/evolve.h"
#include "smooth/smooth.h"

void MSR::BHAccretion(double dTime) {
    struct inBHAccretion in;
    in.dScaleFactor = csmTime2Exp(csm, dTime);
    pstBHAccretion(pst, &in, sizeof(in), NULL, 0);
}

struct bhAccretionPack {
    uint64_t iOrder;
    uint8_t iClass;
};

struct bhAccretionFlush {
    blitz::TinyVector<double, 3> mom;
    uint64_t iOrder;
    float fMass;
};

class SmoothHelperAccretionFlush : public SmoothHelperPacked<bhAccretionPack, bhAccretionFlush> {
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final {}
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;
    virtual void init(void* dst) override final;
    virtual void flush(void* dst, const void* src) override final;
    virtual void combine(void* dst, const void* src, const void* key) override final;

  public:
    SmoothHelperAccretionFlush(pkdContext* pkd) : SmoothHelperPacked<bhAccretionPack, bhAccretionFlush>(pkd, true) {}
};

void SmoothHelperAccretionFlush::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->iOrder = p2.order();
    p1->iClass = p2.get_class();
}

void SmoothHelperAccretionFlush::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_order(p2->iOrder);
    p1.set_class(p2->iClass);
}

void SmoothHelperAccretionFlush::init(void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];

    if (p.is_bh()) {
        p.velocity() = 0.0;
        p.set_mass(0.0);
    }
}

void SmoothHelperAccretionFlush::flush(void* dst, const void* src) {
    auto p1 = static_cast<flush_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    if (p2.is_bh()) {
        p1->mom = p2.velocity(); // **Momentum** added by the accretion
        p1->iOrder = p2.order();
        p1->fMass = p2.mass();
    }
}

void SmoothHelperAccretionFlush::combine(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const flush_type*>(src);

    if (p1.is_bh()) {
        assert(p1.order() == p2->iOrder);
        float old_mass = p1.mass();
        float new_mass = old_mass + p2->fMass;
        float inv_mass = 1. / new_mass;

        auto& v1 = p1.velocity();
        v1 = (old_mass * v1 + p2->mom) * inv_mass;

        p1.set_mass(new_mass);
    }
}

void pkdBHAccretion(PKD pkd, double dScaleFactor) {
    pkd->mdl->CacheInitialize(CID_PARTICLE, std::make_shared<SmoothHelperAccretionFlush>(pkd));
    for (auto& p : pkd->particles) {
        if (p.is_gas()) {
            auto& sph = p.sph();
            auto& iPid = sph.BHAccretor.iPid;
            auto& iIndex = sph.BHAccretor.iIndex;
            if (iPid != NOT_ACCRETED) {
                particleStore::ParticlePointer bh(pkd->particles);
                if (iPid != pkd->Self()) {
                    bh = &pkd->particles[pkd->mdl->Acquire<PARTICLE>(CID_PARTICLE, iIndex, iPid)];
                }
                else {
                    bh = &pkd->particles[iIndex];
                }
                assert(bh->is_bh());

                const float bhMass = bh->mass();
                bh->set_mass(bhMass + p.mass());

                auto& bhv = bh->velocity();

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
                    bhv = (bhMass * bhv + dScaleFactor * sph.mom) * inv_newMass;
                }

                pkdDeleteParticle(pkd, p);

                if (iPid != pkd->Self()) pkd->mdl->Release(CID_PARTICLE, bh);
            }
        }
    }
    pkd->mdl->FinishCache(CID_PARTICLE);
}
