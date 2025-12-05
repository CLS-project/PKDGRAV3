#include "cooling.h"


int ServiceCoolingInit::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<input *>(vin);

    pkd_cooling_init_backend(pst->plcl->pkd, in->in_cooling_data,
                             in->Redshifts,
                             in->nH,
                             in->Temp,
                             in->HeFrac,
                             in->Therm,
                             in->SolarAbundances,
                             in->SolarAbundances_inv);

    return 0;
}


int ServiceCoolingUpdate::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<input *>(vin);
    pkd_cooling_update(pst->plcl->pkd, in);
    return 0;
}

int ServiceCoolingUpdateZ::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<float *>(vin);
    pst->plcl->pkd->cooling->dz = *in;
    return 0;
}

int ServiceCoolingHydReionization::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    cooling_Hydrogen_reionization(pst->plcl->pkd);
    return 0;
}
