// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServicePlaceBHSeed : public TraverseCount<ServicePlaceBHSeed> {
  public:
    static constexpr std::string_view name = "PlaceBHSeed";
    struct input {
        double dTime;
        double dScaleFactor;
        double dDenMin;
        double dBHMhaloMin;
        double dTau;
        double dBHSeedMass;
        uint8_t uRungMax;
    };
    using output = uint64_t;
    explicit ServicePlaceBHSeed(PST pst) : TraverseCount(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<const input*>(vin.data());
        auto out = reinterpret_cast<output*>(vout.data());
        auto pkd = pst->lcl()->pkd;
        *out = pkd->PlaceBHSeed(in->dTime, in->dScaleFactor, in->uRungMax, in->dDenMin, in->dBHMhaloMin, in->dTau, in->dBHSeedMass);
        return 0;
    }
};