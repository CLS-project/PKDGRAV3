#include "blackhole/drift.h"
#include "master.h"

#define NOT_PINNED -1

using blitz::TinyVector;
using blitz::dot;

void MSR::BHGasPin(double dTime, double dDelta) {
    Smooth(dTime,dDelta,SMX_BH_GASPIN,0,parameters.get_nSmooth());
}

void packBHGasPin(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhGasPinPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->velocity = p2.velocity();
        p1->fMass = p2.mass();
        p1->fPotential = p2.potential();
    }
}

void unpackBHGasPin(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhGasPinPack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.velocity() = p2->velocity;
        p1.set_mass(p2->fMass);
        p1.potential() = p2->fPotential;
    }
}

void smBHGasPin(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    const auto pv = p.velocity() / smf->a;
    const NN *nnLowPot;

    TinyVector<vel_t,3> meanv{0.0};
    vel_t meanv2 = 0.0;
    float minPot = HUGE_VAL;
    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        assert(q.is_gas());
        if (q.potential() < minPot) {
            minPot = q.potential();
            nnLowPot = &nnList[i];
        }
        const TinyVector<vel_t,3> v{q.velocity() - pv};
        meanv += v;
        meanv2 += dot(v,v);
    }
    assert(minPot < HUGE_VAL);

    auto pLowPot = pkd->particles[nnLowPot->pPart];
    auto &bh = p.BH();
    if (p.mass() > 10.*pLowPot.mass()) {
        bh.GasPin.iPid = NOT_PINNED;
    }
    else {
        bh.GasPin.iPid = nnLowPot->iPid;
        bh.GasPin.iIndex = nnLowPot->iIndex;

        if (bh.bForceReposition) {
            bh.bForceReposition = false;
        }
        else {
            const double inv_nSmooth = 1./nSmooth;
            meanv *= inv_nSmooth;
            meanv2 *= inv_nSmooth;
            const auto stdv2 = meanv2 - dot(meanv,meanv);

            const TinyVector<vel_t,3> v{pLowPot.velocity() - pv - meanv};
            if (dot(v,v) > stdv2) bh.GasPin.iPid = NOT_PINNED;
        }
    }
}

void MSR::BHReposition() {
    pstBHReposition(pst, NULL, 0, NULL, 0);
}

struct bhRepositionPack {
    TinyVector<double,3> position;
    uint8_t iClass;
};

void packBHReposition(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<bhRepositionPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->position = p2.position();
    p1->iClass = p2.get_class();
}

void unpackBHReposition(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const bhRepositionPack *>(src);

    p1.set_position(p2->position);
    p1.set_class(p2->iClass);
}

void pkdBHReposition(PKD pkd) {
    mdlPackedCacheRO(pkd->mdl, CID_PARTICLE, NULL, pkd->particles,pkd->Local(),
                     pkd->particles.ParticleSize(), pkd, sizeof(bhRepositionPack),
                     packBHReposition, unpackBHReposition);
    for (auto &p : pkd->particles) {
        if (p.is_bh()) {
            auto &bh = p.BH();
            auto &iPid = bh.GasPin.iPid;
            auto &iIndex = bh.GasPin.iIndex;
            if (iPid != NOT_PINNED) {
                particleStore::ParticlePointer pin(pkd->particles);
                if (iPid != pkd->Self()) {
                    pin = &pkd->particles[static_cast<PARTICLE *>(mdlFetch(pkd->mdl,CID_PARTICLE,iIndex,iPid))];
                }
                else {
                    pin = &pkd->particles[iIndex];
                }
                assert(pin->is_gas());
                p.set_position(pin->position());
            }
        }
    }
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
}

