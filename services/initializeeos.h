// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceInitializeEOS : public TraversePST {
  public:
    static constexpr std::string_view name = "InitializeEOS";
    using input = void;
    using output = void;
    explicit ServiceInitializeEOS(PST pst) : TraversePST(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        pst->lcl()->pkd->InitializeEOS();
        return 0;
    }
};