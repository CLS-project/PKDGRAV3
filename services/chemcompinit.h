// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
#ifdef OPTIM_REORDER_IN_NODES
class ServiceChemCompInit : public TraversePST {
  public:
    static constexpr std::string_view name = "ChemCompInit";
    struct input {
        double H;
        double He;
        double C;
        double N;
        double O;
        double Ne;
        double Mg;
        double Si;
        double Fe;
        double Metallicity;
    };
    using output = void;
    explicit ServiceChemCompInit(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto* in = reinterpret_cast<const input*>(vin.data());
        pst->lcl()->pkd->ChemCompInit(in->H, in->He, in->C, in->N, in->O, in->Ne, in->Mg, in->Si, in->Fe, in->Metallicity);
        return 0;
    }
};
#endif