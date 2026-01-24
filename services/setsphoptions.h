// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
#include "SPH/SPHOptions.h"

class ServiceSetSPHoptions : public TraversePST {
  public:
    static constexpr std::string_view name = "SetSPHoptions";
    using input = SPHOptions;
    using output = void;
    explicit ServiceSetSPHoptions(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<input*>(vin.data());
        auto pkd = pst->lcl()->pkd;
        pkd->SPHoptions = *in;
        return 0;
    }
};
