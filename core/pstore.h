// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#pragma once

#include "TraversePST.h"

class ServiceInitializePStore : public TraversePST {
  public:
    static constexpr std::string_view name = "InitializePStore";
    struct input {
        uint64_t mMemoryModel;
        uint64_t nStore;
        uint64_t nDecoratorStore;
        fioSpeciesList nSpecies;
        blitz::TinyVector<double, 3> fPeriod;
        uint64_t nMinEphemeral;
        uint64_t nMinTotalStore;
        uint32_t nIntegerFactor;
        uint32_t nBucket;
        int nEphemeralBytes;
        int nTreeBitsLo;
        int nTreeBitsHi;
        int iCacheSize;
        int iCacheMaxInflight;
        int iWorkQueueSize;
        mdl::MPI_RECV_MODE mpi_recv_mode;
        int mpi_recv_count;
        bool mpi_RDMA;
    };
    struct output {
        size_t nSizeParticle;
        size_t nSizeDecorator;
        size_t nSizeNode;
    };
    explicit ServiceInitializePStore(PST pst, std::string_view service_name = name, int nInBytes = 0)
        : TraversePST(pst, service_name, nInBytes + sizeof(input), sizeof(output)) {}

  private:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int final;
};