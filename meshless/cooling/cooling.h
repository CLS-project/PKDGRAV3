/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_COOLING_EAGLE_H
#define SWIFT_COOLING_EAGLE_H

#include "cooling_struct.h"
#include "services.h"
#include "pkd.h"

namespace cooling {
/* Maximum number of iterations for bisection scheme */
static const int bisection_max_iterations = 150;

/* Tolerances for termination criteria. */
static const float explicit_tolerance = 0.05;
static const float bisection_tolerance = 1.0e-6;
static const double bracket_factor = 1.5;

void get_redshift_index(const float z, int* z_index, float* dz, cooling_function_data* restrict cooling);
} // namespace cooling

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function declarations
 */

#ifdef __cplusplus
extern "C" {
#endif
void cooling_cool_part(PKD pkd, const struct cooling_function_data* cooling, particleStore::ParticleReference& p, decorators::MFM* psph, const float dt,
                       const double time, const float delta_redshift, const double redshift);
float cooling_get_temperature(PKD pkd, const float redshift, const struct cooling_function_data* restrict cooling, particleStore::Particle& p,
                              decorators::MFM* psph);
#ifdef __cplusplus
}
#endif

void pkd_cooling_update(PKD pkd, const ServiceCoolingUpdate::input* in);
void pkd_cooling_init_backend(PKD pkd, struct cooling_function_data in_cooling_data, const float Redshifts[eagle_cooling_N_redshifts],
                              const float nH[eagle_cooling_N_density], const float Temp[eagle_cooling_N_temperature],
                              const float HeFrac[eagle_cooling_N_He_frac], const float Therm[eagle_cooling_N_temperature],
                              const float SolarAbundances[eagle_cooling_N_abundances], const float SolarAbundances_inv[eagle_cooling_N_abundances]);
void cooling_Hydrogen_reionization(PKD pkd);
void cooling_clean(struct cooling_function_data* data);

#endif /* SWIFT_COOLING_EAGLE_H */
