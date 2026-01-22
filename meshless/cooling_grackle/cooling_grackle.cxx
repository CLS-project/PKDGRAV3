#include "cooling_grackle/cooling_grackle.h"

void pkdContext::GrackleInit(int bComove, double dScaleFactor, char* achCoolingTable, UNITS units) {

    grackle_data = (chemistry_data*)malloc(sizeof(chemistry_data));

    chemistry_data temp_data = _set_default_chemistry_parameters();

    memcpy(grackle_data, &temp_data, sizeof(chemistry_data));

    grackle_data->use_grackle = 1;            // chemistry on
    grackle_data->with_radiative_cooling = 1; // cooling on
    grackle_data->primordial_chemistry = 0;   // molecular network with H, He, D
    grackle_data->metal_cooling = 1;          // metal cooling on
    grackle_data->UVbackground = 1;           // UV background on
    grackle_data->cmb_temperature_floor = 100;
    grackle_data->grackle_data_file = achCoolingTable; // data file

    grackle_units = (code_units*)malloc(sizeof(code_units));
    grackle_rates = (chemistry_data_storage*)malloc(sizeof(chemistry_data_storage));
    grackle_field = (grackle_field_data*)malloc(sizeof(grackle_field_data));

    grackle_field->grid_rank = 1;
    grackle_field->grid_dimension = (int*)malloc(sizeof(int));
    grackle_field->grid_start = (int*)malloc(sizeof(int));
    grackle_field->grid_end = (int*)malloc(sizeof(int));

    grackle_field->grid_dimension[0] = 1;
    grackle_field->grid_dimension[1] = 0;
    grackle_field->grid_dimension[2] = 0;
    grackle_field->grid_start[0] = 0;
    grackle_field->grid_end[0] = 0;

    grackle_field->grid_start[1] = 0;
    grackle_field->grid_end[1] = 0;

    grackle_field->grid_start[2] = 0;
    grackle_field->grid_end[2] = 0;

    // Set field arrays.
    grackle_field->density = (gr_float*)malloc(sizeof(gr_float));
    grackle_field->internal_energy = (gr_float*)malloc(sizeof(gr_float));
    grackle_field->x_velocity = (gr_float*)malloc(sizeof(gr_float));
    grackle_field->y_velocity = (gr_float*)malloc(sizeof(gr_float));
    grackle_field->z_velocity = (gr_float*)malloc(sizeof(gr_float));
    // for metal_cooling = 1
    grackle_field->metal_density = (gr_float*)malloc(sizeof(gr_float));
    grackle_units->comoving_coordinates = 1.; // bComove;

    // I do not understand why, but I need to pass the unit conversion directly
    //  to proper, although from the documentation it seems that it should be
    //  comoving...
    // Otherwise, altough the density is correctly converted to proper at
    //   solve_rate_cool.F:342, the hydrogen is converted back to comoving because
    //   dom \propto a^3 (see solve_rate_cool.F:342 and cool1d_cloudy_f.F:124)
    // Furthermore, I need to pass an extra dScaleFactor to length_units to
    //  compensate for the definition of comoving internal energy inside grackle,
    //  which has a a^2.
    // In the current version, length_units is only used for this factor, so it
    //  is safe, *for now* (changeset 59db82b, 14 Jun 2021)
    grackle_units->density_units = units.dGmPerCcUnit * pow(dScaleFactor, -3);
    grackle_units->length_units = units.dKpcUnit * KPCCM * dScaleFactor;
    grackle_units->time_units = units.dSecUnit;
    grackle_units->velocity_units = units.dKmPerSecUnit * 1e5;
    grackle_units->a_units = 1.;
    grackle_units->a_value = dScaleFactor;

    _initialize_chemistry_data(grackle_data, grackle_rates, grackle_units);
}

/* Each timestep, the Grackle data is inititialized.
 * This is suboptimal. But this is the way that I found that can easily
 * take into account the varying scale-factor.
 */
void pkdContext::GrackleUpdate(double dScaleFactor, char* achCoolingTable, UNITS units) {
    _free_chemistry_data(grackle_data, grackle_rates);
    free(grackle_data);
    free(grackle_field->density);
    free(grackle_field->internal_energy);
    free(grackle_field->metal_density);
    free(grackle_field->x_velocity);
    free(grackle_field->y_velocity);
    free(grackle_field->z_velocity);
    free(grackle_field->grid_dimension);
    free(grackle_field->grid_start);
    free(grackle_field->grid_end);
    free(grackle_field);

    GrackleInit(1, dScaleFactor, achCoolingTable, units);

    /* This should work and be way faster... Instead of doing the init all the time*/
    // grackle_units->a_value = dScaleFactor;
}

void pkdContext::GrackleCooling(particleStore::ParticleReference& p, double pDelta, double dTuFac) {
    auto& sph = p.sph();
    gr_float fDensity = p.density();
    gr_float fMetalDensity = sph.fMetalMass * sph.omega;
    gr_float minUint = 100. * dTuFac;
    gr_float fSpecificUint = sph.Uint / p.mass();

    // Set field arrays.
    grackle_field->density[0] = fDensity;
    grackle_field->internal_energy[0] = (fSpecificUint > minUint) ? fSpecificUint : minUint;
    grackle_field->x_velocity[0] = 1.; // Velocity input is not used
    grackle_field->y_velocity[0] = 1.;
    grackle_field->z_velocity[0] = 1.;
    // for metal_cooling = 1
    grackle_field->metal_density[0] = fMetalDensity;

    int err;
    gr_float temperature;
    err = local_calculate_temperature(grackle_data, grackle_rates, grackle_units, grackle_field, &temperature);
    if (err == 0) fprintf(stderr, "Error in calculate_temperature.\n");

    if (pDelta > 0) {
        err = local_solve_chemistry(grackle_data, grackle_rates, grackle_units, grackle_field, pDelta);
        if (err == 0) fprintf(stderr, "Error in calculate_cooling_time.\n");
        sph.E -= sph.Uint;
        sph.Uint = grackle_field->internal_energy[0] * p.mass();
        sph.E += sph.Uint;
    }
    double diff = (fSpecificUint - grackle_field->internal_energy[0]) / fSpecificUint;
}
