
inline static void internalEnergyFloor(PKD pkd, PARTICLE* p, SPHFIELDS* psph, const float a_m3){

    /* First, the cooling floor */
   const float mass = pkdMass(pkd,p);
   const float dens = pkdDensity(pkd,p)*a_m3;
   psph->E -= psph->Uint;
   
    if ( (dens > pkd->param.dCoolingFloorDen) &&  
         (psph->Uint < pkd->param.dCoolingFlooru*mass ) ){
       //printf("Increasing Uint from %e to %e \n", psph->Uint, pkd->param.dCoolingFlooru*mass);
       psph->Uint = pkd->param.dCoolingFlooru*mass;
    }

    
    /* Second, the polytropic EoS */
    if (dens > pkd->param.dJeansFloorDen){

        const double minUint = pkd->param.dJeansFlooru*
           pow( dens/pkd->param.dJeansFloorDen , pkd->param.dJeansFloorIndex )*mass;

       //printf("Increasing Uint from %e to %e \n", psph->Uint, minUint);
        if (psph->Uint < minUint) psph->Uint = minUint;
    }



   psph->E += psph->Uint;

}
