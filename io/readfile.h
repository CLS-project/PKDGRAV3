// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#pragma once
#include "TraversePST.h"
#include "io/fio.h"
#include <cstddef>

class ServiceReadFile : public TraversePST {
  public:
    static constexpr std::string_view name = "ReadFile";
    using SpeciesList = std::array<uint64_t, FIO_SPECIES_LAST + 1>;
    struct input {
        SpeciesList nSpeciesStart;
        SpeciesList nSpeciesEnd;
        double dvFac;
        double dTuFac;
        double dOmega0;
        double dOmegab;
        int nProcessors;
        /*char achFilename[PATH_MAX];*/
    };
    static constexpr size_t max_files = 16384;
    static constexpr size_t max_input_size_bytes = sizeof(input) + max_files * (sizeof(fioSpeciesList) + PATH_MAX);
    explicit ServiceReadFile(PST pst) : TraversePST(pst, name, max_input_size_bytes) {}
    using output = void;

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override;
    auto Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceOneNodeReadInit : public TraversePST {
  public:
    static constexpr std::string_view name = "OneNodeReadInit";
    using input = ServiceReadFile::input;
    using output = ServiceReadFile::SpeciesList;
    explicit ServiceOneNodeReadInit(PST pst) : TraversePST(pst, name, sizeof(input), pst->mdl()->Threads() * sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override;
    auto Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};