// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceBHInit : public TraversePST {
  public:
    static constexpr std::string_view name = "BHInit";
    struct input {
        uint8_t uRungMax;
    };
    using output = void;
    explicit ServiceBHInit(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<const input*>(vin.data());
        auto pkd = pst->lcl()->pkd;
        pkd->BHInit(in->uRungMax);
        return 0;
    }
};