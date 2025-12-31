#include "cooling.h"

int ServiceCoolingInit::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());

    pkd_cooling_init_backend(pst->lcl()->pkd, in->in_cooling_data, in->Redshifts, in->nH, in->Temp, in->HeFrac, in->Therm, in->SolarAbundances,
                             in->SolarAbundances_inv);

    return 0;
}

int ServiceCoolingUpdate::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    pkd_cooling_update(pst->lcl()->pkd, in);
    return 0;
}

int ServiceCoolingUpdateZ::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const float*>(vin.data());
    pst->lcl()->pkd->cooling->dz = *in;
    return 0;
}

int ServiceCoolingHydReionization::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    cooling_Hydrogen_reionization(pst->lcl()->pkd);
    return 0;
}
