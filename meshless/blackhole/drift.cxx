#include "blackhole/drift.h"
#include "hydro/hydro.h"

#define NOT_PINNED -1

using blitz::dot;
using blitz::TinyVector;

auto ServiceBHGasPin::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    auto helper = std::make_shared<SmoothHelperBHGasPin>(pkd, in->a);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
    *out = smx.Smooth<smooth::policy::OnlyBH>();
    return sizeof(output);
}

void SmoothHelperBHGasPin::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->velocity = p2.velocity();
        p1->fMass = p2.mass();
        p1->fPotential = p2.potential();
    }
}

void SmoothHelperBHGasPin::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.velocity() = p2->velocity;
        p1.set_mass(p2->fMass);
        p1.potential() = p2->fPotential;
    }
}

void SmoothHelperBHGasPin::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    const auto pv = p.velocity() / a_;
    const NN* nnLowPot;

    TinyVector<vel_t, 3> meanv{0.0};
    vel_t meanv2 = 0.0;
    float minPot = HUGE_VAL;
    for (auto& nn : nnList) {
        auto q = pkd->particles[nn.pPart];
        assert(q.is_gas());
        if (q.potential() < minPot) {
            minPot = q.potential();
            nnLowPot = &nn;
        }
        const TinyVector<vel_t, 3> v{q.velocity() - pv};
        meanv += v;
        meanv2 += dot(v, v);
    }
    assert(minPot < HUGE_VAL);

    auto pLowPot = pkd->particles[nnLowPot->pPart];
    auto& bh = p.BH();
    if (p.mass() > 10. * pLowPot.mass()) {
        bh.GasPin.iPid = NOT_PINNED;
    }
    else {
        bh.GasPin.iPid = nnLowPot->iPid;
        bh.GasPin.iIndex = nnLowPot->iIndex;

        if (bh.bForceReposition) {
            bh.bForceReposition = false;
        }
        else {
            const double inv_nSmooth = 1. / nnList.size();
            meanv *= inv_nSmooth;
            meanv2 *= inv_nSmooth;
            const auto stdv2 = meanv2 - dot(meanv, meanv);

            const TinyVector<vel_t, 3> v{pLowPot.velocity() - pv - meanv};
            if (dot(v, v) > stdv2) bh.GasPin.iPid = NOT_PINNED;
        }
    }
}

struct bhRepositionPack {
    TinyVector<double, 3> position;
    uint8_t iClass;
};

class SmoothHelperBHReposition : public SmoothHelperPacked<bhRepositionPack> {
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final {}
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;

  public:
    SmoothHelperBHReposition(pkdContext* pkd, bool bSymmetric = false) : SmoothHelperPacked<bhRepositionPack>(pkd, bSymmetric) {}
};

void SmoothHelperBHReposition::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->position = p2.position();
    p1->iClass = p2.get_class();
}

void SmoothHelperBHReposition::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_position(p2->position);
    p1.set_class(p2->iClass);
}

void pkdBHReposition(PKD pkd) {
    pkd->mdl->CacheInitialize(CID_PARTICLE, std::make_shared<SmoothHelperBHReposition>(pkd));
    for (auto& p : pkd->particles) {
        if (p.is_bh()) {
            auto& bh = p.BH();
            auto& iPid = bh.GasPin.iPid;
            auto& iIndex = bh.GasPin.iIndex;
            if (iPid != NOT_PINNED) {
                particleStore::ParticlePointer pin(pkd->particles);
                if (iPid != pkd->Self()) {
                    pin = &pkd->particles[pkd->mdl->Fetch<PARTICLE>(CID_PARTICLE, iIndex, iPid)];
                }
                else {
                    pin = &pkd->particles[iIndex];
                }
                assert(pin->is_gas());
                p.set_position(pin->position());
            }
        }
    }
    pkd->mdl->FinishCache(CID_PARTICLE);
}
