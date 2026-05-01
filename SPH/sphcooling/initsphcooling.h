//
// Created by tobias on 14.03.24.
//
#ifndef PKDGRAV3_INITSPHCOOLING_H
#define PKDGRAV3_INITSPHCOOLING_H
#include "TraversePST.h"
#include "SPH/sphcooling/SPHCooling.h"
class ServiceInitSPHCooling : public TraversePST {
public:
    static constexpr std::string_view name = "InitSPHCooling";
    typedef struct coolingVariables input;
    typedef void output;
    explicit ServiceInitSPHCooling(PST pst)
        : TraversePST(pst, name, sizeof(input)) {}
protected:
    int Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) override;
};
#endif //PKDGRAV3_INITSPHCOOLING_H
