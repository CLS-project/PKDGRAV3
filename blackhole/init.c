#include "blackhole/init.h"

void msrSetBlackholeParam(MSR msr){
    msr->param.dEddingtonFactor *= (1e3/MSOLG/msr->param.dMsolUnit )  /
                        pow( 100. / KPCCM / msr->param.dKpcUnit, 3) /
                        msr->param.dSecUnit / msr->param.dBHRadiativeEff ;

    // We precompute the factor such that we only need to multiply
    // AccretionRate by this amount to get E_feed
    msr->param.dBHFeedbackEff = msr->param.dBHFeedbackEff *
        msr->param.dBHRadiativeEff * (1. - msr->param.dBHRadiativeEff) *
        pow( LIGHTSPEED * 1e-5 /msr->param.dKmPerSecUnit ,2);

    // This, in principle, will not be a parameter
    double n_heat = 1.0;

    // We convert from Delta T to energy per mass.
    // This needs to be multiplied by the mass of the gas particle
    msr->param.dBHFeedbackEcrit *= msr->param.dGasConst/
        (msr->param.dConstGamma - 1.)/0.58 * n_heat;

}

void msrBHInit(MSR msr, uint8_t uRungMax) {
   // We reuse this struct for simplicity
   struct inPlaceBHSeed in;

   in.uRungMax = uRungMax;

   pstBHInit(msr->pst, &in, sizeof(in), NULL, 0);

}



void pkdBHInit(PKD pkd, uint8_t uRungMax) {
   for (int i=0; i<pkdLocal(pkd);i++){
      PARTICLE *p = pkdParticle(pkd,i);
      if (pkdIsBH(pkd,p)){
         p->uRung = uRungMax;
         p->uNewRung = uRungMax;
      }
   }
}
