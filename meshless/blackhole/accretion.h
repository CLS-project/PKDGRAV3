// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceBHAccretion : public TraversePST {
  public:
    static constexpr std::string_view name = "BHAccretion";
    struct input {
        double dScaleFactor;
    };
    using output = void;
    explicit ServiceBHAccretion(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<const input*>(vin.data());
        auto pkd = pst->pkd();
        pkd->BHAccretion(in->dScaleFactor);
        return 0;
    }
};