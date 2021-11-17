#include "pkd.h"
#include "master.h"

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
void msrSetEOSParam(MSR msr);
#endif


#ifdef EEOS_POLYTROPE
/*
 * General polytropic eEOS floor
 */
inline static double polytropicEnergyFloor(PKD pkd, float a_inv3, float fDens){
   return pkd->param.dEOSPolyFlooru *
          pow( fDens*a_inv3/pkd->param.dEOSPolyFloorDen,
               pkd->param.dEOSPolyFloorIndex );
}

inline static double polytropicPressureFloor(PKD pkd, float a_inv3, float fDens){
   return polytropicEnergyFloor(pkd, a_inv3, fDens) *
            fDens * (pkd->param.dConstGamma -1.);
}
#endif



#ifdef EEOS_JEANS
/*
 * Pressure floor based on the Jeans criterion, where the only parameter is
 *  the minimum "number" of resolution elements per Jeans length
 *
 * It is similar to polytropicPressureFloor with gamma=4/3, but it takes into
 *  account explicitly the smoothing length
 */
inline static double jeansPressureFloor(PKD pkd, float fDens, float fBall){
   return 1.2/pkd->param.dConstGamma * pkd->param.dEOSNJeans *
                   fDens * fDens * fBall * fBall;

}

inline static double jeansEnergyFloor(PKD pkd, float fDens, float fBall){
   return jeansPressureFloor(pkd, fDens, fBall) /
          (fDens * (pkd->param.dConstGamma - 1.));

}
#endif

inline static void internalEnergyFloor(PKD pkd, PARTICLE* p,
                                       SPHFIELDS* psph, const double a){

   const double a_inv = 1./a;
   const double a_inv3 = a_inv*a_inv*a_inv;
   const float mass = pkdMass(pkd,p);
   const float dens = pkdDensity(pkd,p);
   const float dens_phys = pkdDensity(pkd,p)*a_inv3;
   const float fBall = 2.*pkdBall(pkd,p);
   psph->E -= psph->Uint;

   // First, the cooling floor,
   // which is only applied if the gas is overdense enough */
#ifdef COOLING
   double denMin = pkd->param.dCoolingFloorDen;
#ifdef STAR_FORMATION
   double minOverDens = pkd->param.dSFMinOverDensity;
#else
   double minOverDens = 57.7;
#endif
   if (pkd->csm->val.bComove){
       double rhoCrit0 = 3. * pkd->csm->val.dHubble0 * pkd->csm->val.dHubble0 /
                              (8. * M_PI);
       double denCosmoMin = rhoCrit0 * pkd->csm->val.dOmegab *
                                 minOverDens *
                                 a_inv3; // We do this in proper density

       denMin = ( denCosmoMin > denMin) ? denCosmoMin : denMin;
   }

   if ( (dens_phys > denMin) &&
        (psph->Uint < pkd->param.dCoolingFlooru*mass ) ){
      psph->Uint = pkd->param.dCoolingFlooru*mass;
   }
#endif


#ifdef EEOS_POLYTROPE
   /* Second, the polytropic EoS */
   if (dens_phys > pkd->param.dEOSPolyFloorDen){

       const double minUint =  polytropicEnergyFloor(pkd, a_inv3, dens)*mass;

       if (psph->Uint < minUint) psph->Uint = minUint;
   }
#endif

#ifdef  EEOS_JEANS
   const double Ujeans = jeansEnergyFloor(pkd, dens, fBall)*mass;

   if (psph->Uint < Ujeans)
      psph->Uint = Ujeans;
#endif



   psph->E += psph->Uint;

}
