// SPDX-License-Identifier: GPL-3.0-only
// Copyright (C) Douglas Potter and Joachim Stadel
#include "pst.h"

class ServiceSendArray : public mdl::BasicService {
    PST node_pst;

  public:
    static constexpr std::string_view name = "SendArray";
    struct input {
        double dvFac;
        PKD_FIELD field;
        int iTo;
        int iUnitSize;
        int bMarked;
    };
    using output = void;
    explicit ServiceSendArray(PST node_pst) : BasicService(name, sizeof(input)), node_pst(node_pst) {}

  protected:
    auto operator()(std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto in = reinterpret_cast<input*>(vin.data());
        node_pst->lcl()->pkd->SendArray(in->iTo, in->field, in->iUnitSize, in->dvFac, in->bMarked);
        return 0;
    }
};
