#include "blackhole/init.h"
#include "master.h"

void MSR::SetBlackholeParam() {
    if (!param.bRestart) {
        param.dBHAccretionEddFac *= (1e3/MSOLG/param.units.dMsolUnit )  /
                                    pow( 100. / KPCCM / param.units.dKpcUnit, 3) /
                                    param.units.dSecUnit / param.dBHRadiativeEff ;

        // We precompute the factor such that we only need to multiply
        // AccretionRate by this amount to get E_feed
        param.dBHFBEff *= param.dBHRadiativeEff *
                          pow( LIGHTSPEED * 1e-5 / param.units.dKmPerSecUnit, 2);
    }

    // This, in principle, will not be a parameter
    double n_heat = 1.0;

    // We convert from Delta T to energy per mass.
    // This needs to be multiplied by the mass of the gas particle
    param.dBHFBEcrit = param.dBHFBDT * dTuFacPrimIonised * n_heat;

}

int MSR::ValidateBlackholeParam() {
    if (param.bBHAccretion) {
        if (param.dBHAccretionAlpha <= 0) {
            fprintf(stderr,"ERROR: dBHAccretionAlpha should be positive."
                    "If you want to avoid boosting the Bondi accretion rate,"
                    "just set dBHAccretionAlpha=1.0\n");
            return 0;
        }
    }
    return 1;
}

void MSR::BlackholeInit(uint8_t uRungMax) {
    // We reuse this struct for simplicity
    struct inPlaceBHSeed in;

    in.uRungMax = uRungMax;

    pstBHInit(pst, &in, sizeof(in), NULL, 0);
}

void pkdBHInit(PKD pkd, uint8_t uRungMax) {
    for (auto &p : pkd->particles) {
        if (p.is_bh()) {
            p.set_rung(uRungMax);
            p.set_new_rung(uRungMax);
        }
    }
}
