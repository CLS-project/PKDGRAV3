// SPDX-License-Identifier: GPL-3.0-only
// Copyright (C) Douglas Potter and Joachim Stadel
#include "io/service.h"
class ServiceCheckpoint : public ServiceIO {
  public:
    static constexpr std::string_view name = "Checkpoint";
    using output = void;
    explicit ServiceCheckpoint(PST pst) : ServiceIO(pst, name, sizeof(input)) {}

  protected:
    void IO(PST pst, std::span<std::byte> vin, int iGroup, int iSegment, int nSegment) final;
};
