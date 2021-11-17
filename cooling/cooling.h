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

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function declarations
 */

/* Local includes. */
#include "cooling_struct.h"
/* IA: PKDGRAV3 includes */
#include "pkd.h"
#include "master.h"


void msrSetCoolingParam(MSR msr);
void cooling_update(MSR msr, const float redshift, int sync);
void pkd_cooling_update(PKD pkd, struct inCoolUpdate *in);

void cooling_cool_part(PKD pkd,
                       const struct cooling_function_data *cooling,
                       PARTICLE* p, SPHFIELDS* psph, 
                       const float dt, const double time,
                       const float delta_redshift, const double redshift);



float cooling_get_temperature(PKD pkd, const float redshift,
    const struct cooling_function_data *restrict cooling,
    PARTICLE* p, SPHFIELDS* psph);


void cooling_Hydrogen_reionization(PKD pkd);

void cooling_init_backend(MSR msr);
void pkd_cooling_init_backend(PKD pkd, struct cooling_function_data in_cooling_data,
  float Redshifts[eagle_cooling_N_redshifts],
  float nH[eagle_cooling_N_density],
  float Temp[eagle_cooling_N_temperature],
  float HeFrac[eagle_cooling_N_He_frac],
  float Therm[eagle_cooling_N_temperature],
  float SolarAbundances[eagle_cooling_N_abundances],
  float SolarAbundances_inv[eagle_cooling_N_abundances]); 

void cooling_clean(struct cooling_function_data *data);

#endif /* SWIFT_COOLING_EAGLE_H */
