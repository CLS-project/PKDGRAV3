#include "master.h"
#include "cooling_grackle/cooling_grackle.h"

void MSR::GrackleInit(int bComove, double dScaleFactor) {
    struct inGrackleInit in;
    in.bComove = bComove;
    in.dScaleFactor = dScaleFactor;
    strcpy(in.achCoolingTable,param.achCoolingTables);
    in.units = units;

    pstGrackleInit(pst,&in,sizeof(struct inGrackleInit),NULL,0);
}

#ifdef __cplusplus
extern "C" {
#endif

void pkdGrackleInit(PKD pkd, int bComove, double dScaleFactor, char *achCoolingTable,
                    UNITS units) {

    pkd->grackle_data = (chemistry_data *) malloc(sizeof(chemistry_data));

    chemistry_data temp_data = _set_default_chemistry_parameters();

    memcpy(pkd->grackle_data, &temp_data, sizeof(chemistry_data));

    pkd->grackle_data->use_grackle = 1;            // chemistry on
    pkd->grackle_data->with_radiative_cooling = 1; // cooling on
    pkd->grackle_data->primordial_chemistry = 0;   // molecular network with H, He, D
    pkd->grackle_data->metal_cooling = 1;          // metal cooling on
    pkd->grackle_data->UVbackground = 1;           // UV background on
    pkd->grackle_data->cmb_temperature_floor = 100;
    pkd->grackle_data->grackle_data_file = achCoolingTable; // data file

    pkd->grackle_units = (code_units *)malloc(sizeof(code_units));
    pkd->grackle_rates = (chemistry_data_storage *)malloc(sizeof(chemistry_data_storage));
    pkd->grackle_field = (grackle_field_data *)malloc(sizeof(grackle_field_data));


    pkd->grackle_field->grid_rank = 1;
    pkd->grackle_field->grid_dimension = (int *)malloc(sizeof(int));
    pkd->grackle_field->grid_start = (int *)malloc(sizeof(int));
    pkd->grackle_field->grid_end = (int *)malloc(sizeof(int));

    pkd->grackle_field->grid_dimension[0] = 1;
    pkd->grackle_field->grid_dimension[1] = 0;
    pkd->grackle_field->grid_dimension[2] = 0;

    pkd->grackle_field->grid_start[0] = 0;
    pkd->grackle_field->grid_end[0] = 0;

    pkd->grackle_field->grid_start[1] = 0;
    pkd->grackle_field->grid_end[1] = 0;

    pkd->grackle_field->grid_start[2] = 0;
    pkd->grackle_field->grid_end[2] = 0;

    // Set field arrays.
    pkd->grackle_field->density         = (gr_float *)malloc(sizeof(gr_float));
    pkd->grackle_field->internal_energy = (gr_float *)malloc(sizeof(gr_float));
    pkd->grackle_field->x_velocity      = (gr_float *)malloc(sizeof(gr_float));
    pkd->grackle_field->y_velocity      = (gr_float *)malloc(sizeof(gr_float));
    pkd->grackle_field->z_velocity      = (gr_float *)malloc(sizeof(gr_float));
    // for metal_cooling = 1
    pkd->grackle_field->metal_density   = (gr_float *)malloc(sizeof(gr_float));


    pkd->grackle_units->comoving_coordinates = 1.;// bComove;

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
    pkd->grackle_units->density_units = units.dGmPerCcUnit*pow(dScaleFactor,-3);
    pkd->grackle_units->length_units = units.dKpcUnit * KPCCM * dScaleFactor;
    pkd->grackle_units->time_units = units.dSecUnit;
    pkd->grackle_units->velocity_units = units.dKmPerSecUnit*1e5;
    pkd->grackle_units->a_units = 1.;
    pkd->grackle_units->a_value = dScaleFactor;

    _initialize_chemistry_data(pkd->grackle_data, pkd->grackle_rates, pkd->grackle_units);
}


/* Each timestep, the Grackle data is inititialized.
 * This is suboptimal. But this is the way that I found that can easily
 * take into account the varying scale-factor.
 */
void pkdGrackleUpdate(PKD pkd, double dScaleFactor, char *achCoolingTable, UNITS units) {
    _free_chemistry_data(pkd->grackle_data, pkd->grackle_rates);
    free(pkd->grackle_data);
    free(pkd->grackle_field->density);
    free(pkd->grackle_field->internal_energy);
    free(pkd->grackle_field->metal_density);
    free(pkd->grackle_field->x_velocity);
    free(pkd->grackle_field->y_velocity);
    free(pkd->grackle_field->z_velocity);
    free(pkd->grackle_field->grid_dimension);
    free(pkd->grackle_field->grid_start);
    free(pkd->grackle_field->grid_end);
    free(pkd->grackle_field);

    pkdGrackleInit(pkd, 1, dScaleFactor, achCoolingTable, units);

    /* This should work and be way faster... Instead of doing the init all the time*/
    //pkd->grackle_units->a_value = dScaleFactor;

}



void pkdGrackleCooling(PKD pkd, particleStore::ParticleReference &p, double pDelta,
                       double dTuFac) {
    auto &sph = p.sph();
    gr_float fDensity = p.density();
    gr_float fMetalDensity = sph.fMetalMass*sph.omega;
    gr_float minUint = 100. * dTuFac;
    gr_float fSpecificUint = sph.Uint/p.mass();


    // Set field arrays.
    pkd->grackle_field->density[0]         = fDensity;
    pkd->grackle_field->internal_energy[0] = (fSpecificUint > minUint) ? fSpecificUint : minUint;
    pkd->grackle_field->x_velocity[0]      = 1.; // Velocity input is not used
    pkd->grackle_field->y_velocity[0]      = 1.;
    pkd->grackle_field->z_velocity[0]      = 1.;
    // for metal_cooling = 1
    pkd->grackle_field->metal_density[0]   = fMetalDensity;

    int err;
    gr_float temperature;
    err = local_calculate_temperature(pkd->grackle_data, pkd->grackle_rates, pkd->grackle_units, pkd->grackle_field,
                                      &temperature);
    if (err == 0) fprintf(stderr, "Error in calculate_temperature.\n");



    if (pDelta > 0) {
        err = local_solve_chemistry( pkd->grackle_data, pkd->grackle_rates, pkd->grackle_units,
                                     pkd->grackle_field, pDelta );
        if (err == 0) fprintf(stderr, "Error in calculate_cooling_time.\n");
        sph.E -= sph.Uint;
        sph.Uint = pkd->grackle_field->internal_energy[0]*p.mass();
        sph.E += sph.Uint;
    }
    double diff = (fSpecificUint-pkd->grackle_field->internal_energy[0])/fSpecificUint;


}

#ifdef __cplusplus
}
#endif
