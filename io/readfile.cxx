// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include <vector>

#include "readfile.h"
#include "core/swapall.h"

auto ServiceOneNodeReadInit::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());
    auto pout = span_cast<output>(vout);
    output nFileStart, nFileEnd, nFileTotal, nFileSplit;

    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        nFileStart[s] = in->nSpeciesStart[s];
        nFileEnd[s] = in->nSpeciesEnd[s];
        nFileTotal[s] = nFileEnd[s] - nFileStart[s] + 1;
    }
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        nFileSplit[s] = nFileStart[s] + pst->n_lower() * (nFileTotal[s] / pst->n_total());
        in->nSpeciesStart[s] = nFileSplit[s];
    }
    auto rID = ReqService(pst, vin);
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        in->nSpeciesStart[s] = nFileStart[s];
        in->nSpeciesEnd[s] = nFileSplit[s] - 1;
    }
    Traverse(pst->lower(), vin, pout.subspan(0, pst->n_lower()));
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        in->nSpeciesEnd[s] = nFileEnd[s];
    }
    GetReply(rID, pout.subspan(pst->n_lower(), pst->n_upper()));

    return sizeof(output) * pst->n_total();
}

auto ServiceOneNodeReadInit::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());
    auto& out = *reinterpret_cast<output*>(vout.data());
    output nFileStart, nFileEnd, nFileTotal, nFileSplit;

    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        nFileStart[s] = in->nSpeciesStart[s];
        nFileEnd[s] = in->nSpeciesEnd[s];
        nFileTotal[s] = nFileEnd[s] - nFileStart[s] + 1;
    }

    // Determine the size of the local particle store.
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        out[s] = nFileTotal[s]; /* Truncated: okay */
    }

    return sizeof(output);
}

auto ServiceReadFile::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());

    // If we only write with one processor at this point, then stop recursing
    if (in->nProcessors <= 1) return Service(pst, vin, vout);

    SpeciesList nNodeSpeciesStart, nNodeSpeciesEnd, nNodeSpeciesTotal, nNodeSpeciesSplit;
    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        nNodeSpeciesStart[s] = in->nSpeciesStart[s];                          // First element to read
        nNodeSpeciesEnd[s] = in->nSpeciesEnd[s];                              // Last element to read
        nNodeSpeciesTotal[s] = nNodeSpeciesEnd[s] - nNodeSpeciesStart[s] + 1; // Total count of read elements
    }

    int nProcessors = in->nProcessors;
    int nProcUpper = pst->n_upper() * nProcessors / pst->n_total();
    int nProcLower = nProcessors - nProcUpper;

    // Index of the first element in the second half of the split
    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1))
        nNodeSpeciesSplit[s] = nNodeSpeciesStart[s] + pst->n_lower() * (nNodeSpeciesTotal[s] / pst->n_total());

    // Upper split
    in->nProcessors = nProcUpper;
    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        in->nSpeciesStart[s] = nNodeSpeciesSplit[s];
        in->nSpeciesEnd[s] = nNodeSpeciesEnd[s];
    }
    auto rID = ReqService(pst, vin);

    // Lower split
    in->nProcessors = nProcLower;
    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        in->nSpeciesStart[s] = nNodeSpeciesStart[s];
        in->nSpeciesEnd[s] = nNodeSpeciesSplit[s] - 1;
    }
    Traverse(pst->lower(), vin);
    GetReply(rID);

    return 0;
}

auto ServiceReadFile::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());
    SpeciesList nNodeSpeciesStart, nNodeSpeciesEnd, nNodeSpeciesTotal, nNodeSpeciesSplit;
    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        nNodeSpeciesStart[s] = in->nSpeciesStart[s];                          // First element to read
        nNodeSpeciesEnd[s] = in->nSpeciesEnd[s];                              // Last element to read
        nNodeSpeciesTotal[s] = nNodeSpeciesEnd[s] - nNodeSpeciesStart[s] + 1; // Total count of read elements
    }

    // Number of particles to be read per core per species
    std::vector<SpeciesList> nSpeciesParts(pst->n_total());

    PKD pkd = pst->lcl()->pkd;
    auto mdl = pkd->mdl;

    Traverse<ServiceOneNodeReadInit>(pst, *in, nSpeciesParts);

    FIO fio = fioLoad(in + 1, in->dOmega0, in->dOmegab);
    assert(fio != nullptr);

    SpeciesList nRead, nStart, nBasePart;
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        // We take care of the first nSpeciesParts later
        nStart[s] = in->nSpeciesStart[s] + nSpeciesParts[0][s];
    }
    for (auto i = 1; i < pst->n_total(); ++i) {
        int id = mdl->Self() + i;
        int inswap;
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            if (in->nSpeciesEnd[s] == 0) continue;
            if (nSpeciesParts[i][s] == 0) continue;

            // Read particles into the local storage.
            assert(pkd->FreeStore() >= nSpeciesParts[i][s]);
            pkdReadFIO(pkd, fio, nStart[s], nSpeciesParts[i][s], in->dvFac, in->dTuFac, FIO_SPECIES(s));
            nStart[s] += nSpeciesParts[i][s];
        }

        // Now shove them over to the remote processor.
        pkd->swap_classes(id);
        inswap = mdl->Self();
        auto rID = pst->mdl()->ReqService<ServiceSwapAll>(id, inswap);
        pkdSwapAll(pkd, id);
        pst->mdl()->GetReply(rID);
    }
    // Finish with our own particles
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        if (nSpeciesParts[0][s] == 0) continue;
        pkdReadFIO(pkd, fio, in->nSpeciesStart[s], nSpeciesParts[0][s], in->dvFac, in->dTuFac, FIO_SPECIES(s));
    }
    fioClose(fio);

    return 0;
}
