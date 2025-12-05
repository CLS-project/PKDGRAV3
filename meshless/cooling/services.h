#ifndef COOLING_SERVICES
#define COOLING_SERVICES
#include "TraversePST.h"
#include "cooling_struct.h"

class ServiceCoolingInit : public TraversePST {
public:
    static constexpr std::string_view name = "CoolingInit";
    struct input {
        struct cooling_function_data in_cooling_data;
        float Redshifts[eagle_cooling_N_redshifts];
        float nH[eagle_cooling_N_density];
        float Temp[eagle_cooling_N_temperature];
        float HeFrac[eagle_cooling_N_He_frac];
        float Therm[eagle_cooling_N_temperature];
        float SolarAbundances[eagle_cooling_N_temperature];
        float SolarAbundances_inv[eagle_cooling_N_temperature];
    };
    explicit ServiceCoolingInit(PST pst)
        : TraversePST(pst,name) {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceCoolingUpdate : public TraversePST {
public:
    static constexpr std::string_view name = "CoolingUpdate";
    struct input {
        int z_index;
        int previous_z_index;
        float dz;
        float metal_heating[eagle_cooling_N_loaded_redshifts * num_elements_metal_heating ];
        float H_plus_He_heating[eagle_cooling_N_loaded_redshifts * num_elements_HpHe_heating];
        float H_plus_He_electron_abundance[eagle_cooling_N_loaded_redshifts * num_elements_HpHe_electron_abundance];
        float temperature[eagle_cooling_N_loaded_redshifts * num_elements_temperature];
        float electron_abundance[eagle_cooling_N_loaded_redshifts * num_elements_electron_abundance];
    };
    explicit ServiceCoolingUpdate(PST pst)
        : TraversePST(pst,name) {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceCoolingUpdateZ : public TraversePST {
public:
    static constexpr std::string_view name = "CoolingUpdateZ";
    struct input {
        double dTime;
    };
    explicit ServiceCoolingUpdateZ(PST pst)
        : TraversePST(pst,name) {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceCoolingHydReionization : public TraversePST {
public:
    static constexpr std::string_view name = "CoolingHydReionization";
    struct input {
        double dTime;
    };
    explicit ServiceCoolingHydReionization(PST pst)
        : TraversePST(pst,name) {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};
#endif
