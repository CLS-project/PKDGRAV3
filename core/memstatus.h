// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"

class ServiceMemStatus : public TraverseGatherPST<ServiceMemStatus> {
  public:
    static constexpr std::string_view name = "MemStatus";
    using input = void;
    struct output {
#ifdef __linux__
        uint64_t rss;
        uint64_t freeMemory;
#endif
        uint64_t nBytesTree;
        uint64_t nBytesCl;
        uint64_t nBytesIlp;
        uint64_t nBytesIlc;
    };
    explicit ServiceMemStatus(PST pst) : TraverseGatherPST(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int final;
};
