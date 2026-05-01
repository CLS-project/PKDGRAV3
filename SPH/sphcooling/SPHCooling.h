//
// Created by tobias on 11.03.24.
//
#ifndef PKDGRAV3_SPHCOOLING_H
#define PKDGRAV3_SPHCOOLING_H
#include "pkd.h"
#include "pst.h"
#include "../SPHOptions.h"
struct coolingVariables {
    struct UNITS units;
    char *achCoolingTable;
};
//namespace SPHCooling {
/// @brief Initializes SPH cooling
/// @param pkd the pkd context
/// @param cool_vars struct with all necessary variables (including units and cooling table)
/// Initial setup of units and chemistry objects. This should be done at simulation start.
void pkdSPHCoolingInit(PKD pkd, struct coolingVariables *cool_vars);
/// @brief Updates the scale factor of SPH cooling when running cosmology simulation
/// @param pkd the pkd context
/// @param dScaleFactor new cosmological scale factor
void pkdSPHCoolingUpdate(PKD pkd, double dScaleFactor);
/// @brief Applies cooling to particle
/// @param pkd the pkd context
/// @param p the particle to apply cooling to
/// @param T the temperature of the particle
/// @return uDot from cooling
double pkdSPHCooling(PKD pkd, particleStore::ParticleReference &p, float T, SPHOptions *SPHoptions);
//} // namespace SPHCooling
#endif //PKDGRAV3_SPHCOOLING_H
