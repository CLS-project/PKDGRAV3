// SPDX-License-Identifier: GPL-3.0-only
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
#ifdef GRACKLE
class ServiceGrackleInit : public TraversePST {
  public:
    static constexpr std::string_view name = "GrackleInit";
    struct input {
        char achCoolingTable[PATH_MAX];
        UNITS units;
        int bComove;
        double dScaleFactor;
    };
    using output = void;
    explicit ServiceGrackleInit(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override {
        pst->lcl()->pkd->GrackleInit(in->bComove, in->dScaleFactor, in->achCoolingTable, in->units);
        return 0;
    }
};
#endif