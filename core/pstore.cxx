// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "pstore.h"

auto ServiceInitializePStore::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto mdl = pst->mdl();

    auto pkd = pst->set_pkd(std::make_unique<pkdContext>(mdl, in->nStore, in->nDecoratorStore, in->nBucket, in->nMinTotalStore, in->nMinEphemeral,
                                                         in->nEphemeralBytes, in->nTreeBitsLo, in->nTreeBitsHi, in->iCacheSize, in->iCacheMaxInflight,
                                                         in->iWorkQueueSize, in->mpi_recv_mode, in->mpi_recv_count, in->mpi_RDMA, in->fPeriod,
                                                         in->nSpecies[FIO_SPECIES_DARK], in->nSpecies[FIO_SPECIES_SPH], in->nSpecies[FIO_SPECIES_STAR],
                                                         in->nSpecies[FIO_SPECIES_BH], in->mMemoryModel, in->nIntegerFactor));
    out->nSizeParticle = pkd->particles.ElementSize();
    out->nSizeDecorator = pkd->decorators.ElementSize();
    out->nSizeNode = pkd->tree.ElementSize();
    return sizeof(output);
}