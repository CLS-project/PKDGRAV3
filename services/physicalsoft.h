// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServicePhysicalSoft : public TraversePST {
  public:
    static constexpr std::string_view name = "PhysicalSoft";
    struct input {
        double dSoftMax;
        double dFac;
        bool bSoftMaxMul;
    };
    using output = void;
    explicit ServicePhysicalSoft(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<input*>(vin.data());
        auto pkd = pst->lcl()->pkd;
        pkd->PhysicalSoft(in->dSoftMax, in->dFac, in->bSoftMaxMul);
        return 0;
    }
};