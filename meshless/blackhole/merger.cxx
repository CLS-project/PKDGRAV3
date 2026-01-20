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

struct SmoothOnlyBH : SmoothEverything {
#ifdef OPTIM_REORDER_IN_NODES
    static inline auto range(treeStore::NodePointer node) {
        auto start = node->lower();
        start += node->Ngas();
        if (node->have_Nstar()) start += node->Nstar();
        auto end = start + node->Nbh();
        return std::pair{start, end};
    }
    static inline auto range_neighbour(treeStore::NodePointer node) { return range(node); }
    static inline auto may_have_active(treeStore::NodePointer& node) -> bool {
    #ifdef OPTIM_REORDER_IN_NODES
        return node->Nbh() > 0;
    #else
        return std::any_of(node->begin(), node->end(), [](auto& p) { return p.is_bh(); });
    #endif
    }
    // The raange() method filters particles, so is_active() will always check only BH particles
    static inline auto is_bh(particleStore::ParticleReference& p) -> bool { return true; }
#else
    static inline auto is_bh(particleStore::ParticleReference& p) -> bool { return p.is_bh(); }
#endif
    // DP: It looks like blackholes have to be marked and get unmarked during smoothing
    //     for example when there is a merger and one of them is deleted.
    static inline auto is_active(particleStore::ParticleReference& p) -> bool { return is_bh(p) && p.marked(); }
    static inline auto accept_neighbour(particleStore::ParticleReference& p) -> bool { return p.is_bh(); }
    static inline auto accept_sink(particleStore::ParticleReference& p) -> bool { return p.marked(); }
    // Use twice the softening length for BHs
    static constexpr auto ball2_sink(particleStore::ParticleReference& p) -> float { return 4 * p.soft() * p.soft(); }
};

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
