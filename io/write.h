// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
class ServiceWrite : public TraversePST {
  public:
    static constexpr std::string_view name = "Write";
    struct input {
        Bound bnd;
        UNITS units;
        double dTime;
        double dExp;
        double dEcosmo;
        double dTimeOld;
        double dUOld;
        double dvFac;
        double dTuFac;
        double dBoxSize;
        double Omega0;
        double OmegaLambda;
        double HubbleParam;
        uint64_t nGas;
        uint64_t nDark;
        uint64_t nStar;
        uint64_t nBH;
        int bStandard;
        int iIndex;
        int nProcessors;
        int nProcessorsTotal;
        int iLower, iUpper;
        int bHDF5;
        int mFlags;
        char achOutFile[PATH_MAX];
    };
    using output = void;
    explicit ServiceWrite(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override;
    auto Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};