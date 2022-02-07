/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_STRUCT_EAGLE_H
#define SWIFT_COOLING_STRUCT_EAGLE_H

#include "units.h"

#define eagle_table_path_name_length 256




/*! Number of different bins along the redhsift axis of the tables */
#define eagle_cooling_N_redshifts 49

/*! Number of redshift bins loaded at any given point int time */
#define eagle_cooling_N_loaded_redshifts 2

/*! Number of different bins along the temperature axis of the tables */
#define eagle_cooling_N_temperature 176

/*! Number of different bins along the density axis of the tables */
#define eagle_cooling_N_density 41

/*! Number of different bins along the metal axis of the tables */
#define eagle_cooling_N_metal 9

/*! Number of different bins along the metal axis of the tables */
#define eagle_cooling_N_He_frac 7

/*! Number of different bins along the abundances axis of the tables */
#define eagle_cooling_N_abundances 11


/**
 * @brief Names of the elements in the order they are stored in the files
 */
static const char *eagle_tables_element_names[eagle_cooling_N_metal] = {
    "Carbon",  "Nitrogen", "Oxygen",  "Neon", "Magnesium",
    "Silicon", "Sulphur",  "Calcium", "Iron"
};

/*! Number of elements in a z-slice of the H+He cooling rate tables */
#define num_elements_cooling_rate ( eagle_cooling_N_temperature * eagle_cooling_N_density)

/*! Number of elements in a z-slice of the metal cooling rate tables */
#define num_elements_metal_heating (eagle_cooling_N_metal * eagle_cooling_N_temperature * eagle_cooling_N_density)

/*! Number of elements in a z-slice of the metal electron abundance tables */
#define num_elements_electron_abundance (eagle_cooling_N_temperature * eagle_cooling_N_density)

/*! Number of elements in a z-slice of the temperature tables */
#define num_elements_temperature (eagle_cooling_N_He_frac * eagle_cooling_N_temperature * eagle_cooling_N_density)

/*! Number of elements in a z-slice of the H+He cooling rate tables */
#define num_elements_HpHe_heating (eagle_cooling_N_He_frac * eagle_cooling_N_temperature * eagle_cooling_N_density)

/*! Number of elements in a z-slice of the H+He electron abundance tables */
#define num_elements_HpHe_electron_abundance (eagle_cooling_N_He_frac * eagle_cooling_N_temperature * eagle_cooling_N_density)








/**
 * @brief struct containing cooling tables
 */
struct cooling_tables {

    /* array of heating rates due to metals */
    float *metal_heating;

    /* array of heating rates due to hydrogen and helium */
    float *H_plus_He_heating;

    /* array of electron abundances due to hydrogen and helium */
    float *H_plus_He_electron_abundance;

    /* array of temperatures */
    float *temperature;

    /* array of electron abundances due to metals */
    float *electron_abundance;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {
    UNITS units;

    /*! Cooling tables */
    struct cooling_tables table;

    /*! Redshift bins */
    float *Redshifts;

    /*! Hydrogen number density bins */
    float *nH;

    /*! Temperature bins */
    float *Temp;

    /*! Helium fraction bins */
    float *HeFrac;

    /*! Internal energy bins */
    float *Therm;

    /*! Mass fractions of elements for solar abundances (from the tables) */
    float *SolarAbundances;

    /*! Inverse of the solar mass fractions */
    float *SolarAbundances_inv;

    /*! Filepath to the directory containing the HDF5 cooling tables */
    char cooling_table_path[eagle_table_path_name_length];

    /*! Redshift of H reionization */
    float H_reion_z;

    /*! H reionization energy in CGS units */
    float H_reion_heat_cgs;

    /*! Have we already done H reioisation? */
    int H_reion_done;

    /*! Ca over Si abundance divided by the solar ratio for these elements */
    float Ca_over_Si_ratio_in_solar;

    /*! S over Si abundance divided by the solar ratio for these elements */
    float S_over_Si_ratio_in_solar;

    /*! Redshift of He reionization */
    float He_reion_z_centre;

    /*! Spread of the He reionization */
    float He_reion_z_sigma;

    /*! He reionization energy in CGS units */
    float He_reion_heat_cgs;

    double dConstGamma;

    /*! Internal energy conversion from internal units to CGS (for quick access)
     */
    double internal_energy_to_cgs;

    /*! Internal energy conversion from CGS to internal units (for quick access)
     */
    double internal_energy_from_cgs;

    /*! Number density conversion from internal units to CGS (for quick access) */
    double number_density_to_cgs;

    /*! Inverse of proton mass in cgs (for quick access) */
    double inv_proton_mass_cgs;

    /*! Temperature of the CMB at present day (for quick access) */
    double T_CMB_0;

    /*! Compton rate in cgs units */
    double compton_rate_cgs;

    /*! Index of the current redshift along the redshift index of the tables */
    int z_index;

    /*! Distance between the current redshift and table[z_index] */
    float dz;

    /*! Index of the previous tables along the redshift index of the tables */
    int previous_z_index;
};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {

    /*! Cumulative energy radiated by the particle */
    float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_EAGLE_H */
