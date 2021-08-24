
inline static void internalEnergyFloor(PKD pkd, PARTICLE* p,
                                       SPHFIELDS* psph, const float a_m3){

   const float mass = pkdMass(pkd,p);
   const float dens = pkdDensity(pkd,p)*a_m3;
   psph->E -= psph->Uint;

   // First, the cooling floor,
   // which is only applied if the gas is overdense enough */
#ifdef COOLING
   double denMin = pkd->param.dCoolingFloorDen;
   if (pkd->csm->val.bComove){
       double rhoCrit0 = 3. * pkd->csm->val.dHubble0 * pkd->csm->val.dHubble0 /
                              (8. * M_PI);
       double denCosmoMin = rhoCrit0 * pkd->csm->val.dOmegab *
                                 pkd->param.dSFMinOverDensity*
                                 a_m3; // We do this in proper density

       denMin = ( denCosmoMin > denMin) ? denCosmoMin : denMin;
   }

   if ( (dens > denMin) &&
        (psph->Uint < pkd->param.dCoolingFlooru*mass ) ){
      psph->Uint = pkd->param.dCoolingFlooru*mass;
   }
#endif


   /* Second, the polytropic EoS */
   if (dens > pkd->param.dJeansFloorDen){

       const double minUint = pkd->param.dJeansFlooru*
          pow( dens/pkd->param.dJeansFloorDen , pkd->param.dJeansFloorIndex )*mass;

       if (psph->Uint < minUint) psph->Uint = minUint;
   }



   psph->E += psph->Uint;

}
