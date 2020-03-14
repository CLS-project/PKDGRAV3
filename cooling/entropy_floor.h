
inline static void internalEnergyFloor(PKD pkd, PARTICLE* p, SPHFIELDS* psph){

    /* First, the cooling floor */
   const float mass = pkdMass(pkd,p);
   const float dens = pkdDensity(pkd,p);
    if ( (dens < pkd->param.dCoolingFloorDen) &&  
         (psph->Uint < pkd->param.dCoolingFlooru*mass ) )
       psph->Uint = pkd->param.dCoolingFlooru*mass;

    /* Second, the polytropic EoS */
    if (dens > pkd->param.dJeansFloorDen){
        const double minUint = pkd->param.dJeansFlooru*pow( dens/pkd->param.dJeansFloorDen , pkd->param.dJeansFloorIndex )*mass;
        if (psph->Uint < minUint) psph->Uint = minUint;
    }

}
