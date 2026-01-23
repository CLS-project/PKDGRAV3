// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceSmooth : public TraverseCount<ServiceSmooth> {
  public:
    static constexpr std::string_view name = "Smooth";
    struct input {
        int nSmooth;
        int bPeriodic;
        int bSymmetric;
        int iSmoothType;
    };
    using output = uint64_t;
    explicit ServiceSmooth(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceReSmooth : public TraverseCount<ServiceReSmooth> {
  public:
    static constexpr std::string_view name = "ReSmooth";
    using input = ServiceSmooth::input;
    using output = ServiceSmooth::output;
    explicit ServiceReSmooth(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};
