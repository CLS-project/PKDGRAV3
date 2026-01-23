// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
#ifdef OPTIM_REORDER_IN_NODES
class ServiceReorderInNodes : public TraversePST {
  public:
    static constexpr std::string_view name = "ReorderInNodes";
    using input = void;
    using output = void;
    explicit ServiceReorderInNodes(PST pst) : TraversePST(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override {
        pst->lcl()->pkd->ReorderInNodes();
        return 0;
    }
};
#endif