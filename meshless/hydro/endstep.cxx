#include "endstep.h"
#include "hydro.h"
#include "meshless/cooling/cooling.h"
#include "meshless/starformation/feedback.h"
#include "meshless/blackhole/evolve.h"
#include "eEOS/eEOS.h"
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif

int ServiceEndTimestepIntegration::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    double pDelta, dScaleFactor, dHubble;
    blitz::TinyVector<double, 3> pa;

    int bComove = pkd->csm->val.bComove;

    assert(pkd->particles.present(PKD_FIELD::oVelocity));
    assert(pkd->particles.present(PKD_FIELD::oMass) || !pkd->particles.present(PKD_FIELD::oMFV));

    if (bComove) {
        dScaleFactor = csmTime2Exp(pkd->csm, in->dTime);
        dHubble = csmTime2Hub(pkd->csm, in->dTime);
    }
    else {
        dScaleFactor = 1.0;
        dHubble = 0.0;
    }
#ifdef GRACKLE
    pkdGrackleUpdate(pkd, dScaleFactor, in->achCoolingTable, in->units);
#endif
    for (auto& p : pkd->particles) {
        if (p.is_gas() && p.is_active()) {
            auto& sph = p.sph();

            // ##### Add ejecta from stellar evolution
#ifdef STELLAR_EVOLUTION
            if (in->bChemEnrich && sph.fReceivedMass > 0.0f) {
                pkdAddStellarEjecta(pkd, p, sph, in->dConstGamma);
            }
#endif

            if (in->dDelta > 0) {
                pDelta = in->dTime - sph.lastUpdateTime;
            }
            else {
                pDelta = 0.0;
            }

            // ##### Gravity
            hydroSourceGravity(pkd, p, &sph, pDelta, pa, dScaleFactor, bComove);

            // ##### Expansion effects
            hydroSourceExpansion(pkd, p, &sph, pDelta, dScaleFactor, dHubble, bComove, in->dConstGamma);

            // ##### Synchronize Uint, Etot (and possibly S)
            hydroSyncEnergies(pkd, p, &sph, pa, in->dConstGamma);

            // ##### Cooling
            if (pkd->cooling) {
                double dRedshift = 1. / dScaleFactor - 1.;
                const float delta_redshift = -pDelta * dHubble * (dRedshift + 1.);
                cooling_cool_part(pkd, pkd->cooling, p, &sph, pDelta, in->dTime, delta_redshift, dRedshift);
            }
#ifdef GRACKLE
            pkdGrackleCooling(pkd, p, pDelta, in->dTuFac);
#endif

            // ##### Effective Equation Of State
            if (in->eEOS.bDoEOS) {
                const double a_inv3 = 1. / (dScaleFactor * dScaleFactor * dScaleFactor);
                const double dFlooru = eEOSEnergyFloor<vec<double, double>, mmask<bool>>(a_inv3, p.density(), p.ball(), in->dConstGamma, in->eEOS);
                if (dFlooru != NOT_IN_EEOS) {
                    const double dEOSUint = p.mass() * dFlooru;
                    if (sph.Uint < dEOSUint) {
                        sph.E = sph.E - sph.Uint;
                        sph.Uint = dEOSUint;
                        sph.E = sph.E + sph.Uint;
                    }
                }
            }

            // ##### Apply feedback
            pkdAddFBEnergy(pkd, p, &sph, in->dConstGamma);

            // Actually set the primitive variables
            hydroSetPrimitives(pkd, p, &sph, in->dTuFac, in->dConstGamma);

            // Set 'last*' variables for next timestep
            hydroSetLastVars(pkd, p, &sph, pa, dScaleFactor, in->dTime, in->dDelta, in->dConstGamma);

            hydroResetFluxes(&sph);
        }
        else if (p.is_bh() && p.is_active()) {
            pkdBHIntegrate(pkd, p, in->dTime, in->dDelta, in->dBHRadiativeEff);
        }
    }
    return 0;
}
