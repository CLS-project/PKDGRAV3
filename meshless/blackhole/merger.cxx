#include <algorithm>
#include "blackhole/merger.h"
#include "hydro/hydro.h"
#include "smooth/smooth.h"

using blitz::abs;
using blitz::all;
using blitz::any;
using blitz::dot;
using blitz::floor;
using blitz::max;
using blitz::TinyVector;

// This should probably not be a template specialisation; combine with the generic version
template <>
auto smContext::ReSmoothNode<SmoothOnlyBH>(SmoothOnlyBH) -> uint64_t {
    uint64_t nSmoothed = 0;
    std::vector<NN> nnList_p;
    nnList_p.reserve(NNLIST_INCREMENT);
    nnList.clear();

    // Here we store the pointers to the particles whose interaction need
    // to be computed
    std::vector<PARTICLE*> sinks;
    sinks.reserve(64); // At most, the size of the bucket

    for (int i = NRESERVED_NODES; i < pkd->Nodes() - 1; ++i) {
        auto node = pkd->tree[i];
        bool bBHinBucket;
#ifdef OPTIM_REORDER_IN_NODES
        bBHinBucket = node->is_bucket() && node->Nbh() > 0;
#else
        bBHinBucket = node->is_bucket() && std::any_of(node->begin(), node->end(), [](auto& p) { return p.is_bh(); });
#endif
        if (bBHinBucket) { // We are in a bucket containing a BH
            auto bnd_node = getNodeSinks<SmoothOnlyBH>(node, sinks);
            buildInteractionList<SmoothOnlyBH>(bnd_node);

            // IA: Now we should have inside nnList all the particles in the bucket (sinks) and those of which can
            //  interact with them from other buckets (smx->nnList)
            //

            for (auto& P : sinks) {
                auto partj = pkd->particles[P];
                // We need to double check this, as the particle may have
                // been marked for deletion just before this
                if (partj.is_marked()) {
                    const double fBall2_p = 4. * partj.soft() * partj.soft();
                    TinyVector<double, 3> dr_node{bnd_node.center() - partj.position()};
                    nnList_p.clear();
                    for (auto& pk : nnList) {
                        if (P == pk.pPart) continue;
                        TinyVector<double, 3> dr{pk.dr - dr_node};
                        const auto fDist2 = dot(dr, dr);
                        if (fDist2 < fBall2_p) {
                            nnList_p.push_back({.pPart = pk.pPart, .fDist2 = fDist2, .dr = dr, .iIndex = pk.iIndex, .iPid = pk.iPid});
                        }
                    }

                    function->Smooth(partj, partj.soft(), nnList_p);
                }
            }

            nSmoothed += sinks.size();

            for (auto& pk : nnList) {
                if (pk.iPid != pkd->Self()) {
                    pkd->mdl->Release(CID_PARTICLE, pk.pPart);
                }
            }
        }
    }
    return nSmoothed;
}

auto ServiceHydroMergerBH::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    auto helper = std::make_shared<SmoothHelperBHmerger>(pkd);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
    *out = smx.ReSmoothNode<SmoothOnlyBH>();
    return sizeof(output);
}

void SmoothHelperBHmerger::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    for (auto& nn : nnList) {
        auto q = pkd->particles[nn.pPart];
        /* IA: Merging across domain boundaries is a tricky thing,
         * This would require careful atomic operations in order to maintain
         * conservation *and* a deterministic behaviour independent of
         * the number of cores.
         *
         * The easiest solution is to relax the deterministic and invariant
         * behaviour and just merge BH that are in the same domain.
         * This should be conservative, and in standard cosmological simulations
         * this situation will be very rare, and eventually both BH will lie
         * in the same computational domain, so I think this won't be a
         * critical problem.
         *
         * Anyway, the most important thing is *conservation*, so that is what
         * I should check for in the test cases
         */
        if (pkd->Self() != nn.iPid) continue;

        if (&p != &q && p.is_bh() && q.is_bh() && p.marked() && q.marked()) {
            // At this point, we are sure that both particles are black holes and
            //  the distance between them is, at most, p.soft()
            //
            // So, the following should not give any problem
            // (remove for performance)
            assert(nn.fDist2 < 4. * p.soft() * p.soft());

            if (p.position(0) >= q.position(0)) {
                const auto& pmass = p.mass();
                auto& pv = p.velocity();
                const auto& qv = q.velocity();
                const auto dv2 = dot(pv - qv, pv - qv);

                const auto& dist2 = nn.fDist2;
                const float inv_dist = dist2 > 0.0 ? sqrt(1. / dist2) : FLT_MAX;
                if (dv2 < pmass * inv_dist) {
                    // We have a merger!!
                    const auto& qmass = q.mass();
                    printf("Merger!!! %e %e \n", pmass, qmass);

                    const float newmass = pmass + qmass;
                    const float inv_newmass = 1. / newmass;

                    auto& pbh = p.BH();
                    p.set_position(p.position() - qmass * inv_newmass * nn.dr);
                    pbh.dInternalMass += q.BH().dInternalMass;
                    pbh.dAccEnergy += q.BH().dAccEnergy;
                    pv = (pmass * pv + qmass * qv) * inv_newmass;
                    p.set_mass(newmass);
                    q.set_mass(0.);

                    // We do this to check that the deleted particles are not
                    // entering the gravity computation, as this will activate
                    // an assert in the gravity loop

                    pkdDeleteParticle(pkd, q);

                    // We only allow for one merger per particle per timestep
                    p.set_marked(false);
                    return;
                }
            }
        }
    }
}
