//
// Created by tobias on 16.03.24.
//
#ifndef PKDGRAV3_COOLINGDATA_H
#define PKDGRAV3_COOLINGDATA_H
extern "C" {
//#include <grackle.h>
}
typedef struct {
    // structs needed for grackle cooling
    //chemistry_data *grackle_data;
    //chemistry_data_storage *grackle_rates;
   // grackle_field_data *grackle_field;
    //code_units *grackle_units;
    // conversion factors needed for KWH cooling
    double cooling_rate_to_cgs;
    double cgs_to_cooling_rate;
    double density_to_cgs;
    double cgs_to_density;
} CoolingData;
#endif //PKDGRAV3_COOLINGDATA_H
