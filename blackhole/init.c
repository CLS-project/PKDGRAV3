#include "blackhole/init.h"

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
