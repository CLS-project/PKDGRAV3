/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "TraversePST.h"
#include <cstddef>

class ServiceGridCreateFFT : public TraversePST {
  public:
    static constexpr std::string_view name = "GridCreateFFT";
    struct input {
        int nGrid;
    };
    using output = void;
    explicit ServiceGridCreateFFT(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceGridDeleteFFT : public TraversePST {
  public:
    static constexpr std::string_view name = "GridDeleteFFT";
    using input = void;
    using output = void;
    explicit ServiceGridDeleteFFT(PST pst) : TraversePST(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceLinearKick : public TraversePST {
  public:
    static constexpr std::string_view name = "LinearKick";
    struct input {
        vel_t dtOpen, dtClose;
    };
    using output = void;
    explicit ServiceLinearKick(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceMeasureLinearPk : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "MeasureLinearPk";
    static constexpr std::size_t max_k_bins = 2500;
    struct input {
        double dA;
        double dBoxSize;
        float fPhase;
        int iSeed;
        int nGrid;
        int nBins;
        bool bFixed;
    };
    struct output {
        double fK[max_k_bins];
        double fPower[max_k_bins];
        uint64_t nPower[max_k_bins];
    };
    explicit ServiceMeasureLinearPk(PST pst) : TraverseCombinePST(pst, name, sizeof(input), sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) -> int override;
};

class ServiceSetLinearGrid : public TraversePST {
  public:
    static constexpr std::string_view name = "SetLinearGrid";
    struct input {
        double a0;
        double a;
        double a1;
        double dBSize;
        int nGrid;
        /* Noise generation */
        int iSeed;
        int bFixed;
        float fPhase;
    };
    using output = void;
    explicit ServiceSetLinearGrid(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};
