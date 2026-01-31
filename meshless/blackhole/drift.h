// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceBHReposition : public TraversePST {
  public:
    static constexpr std::string_view name = "BHReposition";
    using input = void;
    using output = void;
    explicit ServiceBHReposition(PST pst) : TraversePST(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<const input*>(vin.data());
        auto pkd = pst->pkd();
        pkd->BHReposition();
        return 0;
    }
};
