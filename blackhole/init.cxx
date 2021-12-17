#include "blackhole/init.h"
#include "master.h"

void MSR::SetBlackholeParam(){
    param.dBHAccretionEddFac *= (1e3/MSOLG/param.units.dMsolUnit )  /
                        pow( 100. / KPCCM / param.units.dKpcUnit, 3) /
                        param.units.dSecUnit / param.dBHRadiativeEff ;

    // We precompute the factor such that we only need to multiply
    // AccretionRate by this amount to get E_feed
    param.dBHFBEff = param.dBHFBEff *
        param.dBHRadiativeEff * (1. - param.dBHRadiativeEff) *
        pow( LIGHTSPEED * 1e-5 /param.units.dKmPerSecUnit ,2);

    // This, in principle, will not be a parameter
    double n_heat = 1.0;

    // We convert from Delta T to energy per mass.
    // This needs to be multiplied by the mass of the gas particle
    param.dBHFBEcrit = param.dBHFBDT * param.units.dGasConst/
        (param.dConstGamma - 1.)/0.58 * n_heat;

}

void MSR::BlackholeInit(uint8_t uRungMax) {
   // We reuse this struct for simplicity
   struct inPlaceBHSeed in;

   in.uRungMax = uRungMax;

   pstBHInit(pst, &in, sizeof(in), NULL, 0);
}

#ifdef __cplusplus
extern "C" {
#endif


void pkdBHInit(PKD pkd, uint8_t uRungMax) {
   for (int i=0; i<pkdLocal(pkd);i++){
      PARTICLE *p = pkdParticle(pkd,i);
      if (pkdIsBH(pkd,p)){
         p->uRung = uRungMax;
         p->uNewRung = uRungMax;
      }
   }
}
#ifdef __cplusplus
}
#endif
