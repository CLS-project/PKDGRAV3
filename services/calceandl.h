// SPDX-License-Identifier: GPL-3.0-only
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceCalcEandL : public TraverseCount<ServiceCalcEandL> {
  public:
    static constexpr std::string_view name = "CalcEandL";
    using input = void;
    struct output {
        double T;
        double U;
        double Eth;
        blitz::TinyVector<double, 3> L;
        blitz::TinyVector<double, 3> F;
        double W;

        auto operator+=(const output& rhs) -> output& {
            T += rhs.T;
            U += rhs.U;
            Eth += rhs.Eth;
            L += rhs.L;
            F += rhs.F;
            W += rhs.W;
            return *this;
        }
    };
    explicit ServiceCalcEandL(PST pst) : TraverseCount(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto* out = reinterpret_cast<output*>(vout.data());
        pst->lcl()->pkd->CalcEandL(out->T, out->U, out->Eth, out->L, out->F, out->W);
        return sizeof(output);
    }
};