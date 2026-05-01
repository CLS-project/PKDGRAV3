//
// Created by tobias on 14.03.24.
//
#include "initsphcooling.h"
// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceInitSPHCooling::input> || std::is_trivial_v<ServiceInitSPHCooling::input>);
static_assert(std::is_void_v<ServiceInitSPHCooling::output> || std::is_trivial_v<ServiceInitSPHCooling::output>);
int ServiceInitSPHCooling::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<input*>(vin.data());
    static_assert(std::is_void_v<output>);
    assert(vin.size() == sizeof(input));
    (void)vout;
    pkdSPHCoolingInit(pst->pkd(), in);
    return 0;
}
