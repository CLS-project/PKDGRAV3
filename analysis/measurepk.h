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
#include <array>

class ServiceAddLinearSignal : public TraversePST {
  public:
    static constexpr std::string_view name = "AddLinearSignal";
    struct input {
        int iGrid;
        int iSeed;
        int bFixed;
        float fPhase;
        double Lbox;
        double a;
    };
    using output = void;
    explicit ServiceAddLinearSignal(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceDensityContrast : public TraversePST {
  public:
    static constexpr std::string_view name = "DensityContrast";
    struct input {
        double dTotalMass;
        int iGrid;
        int k;
    };
    using output = void;
    explicit ServiceDensityContrast(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceInterlace : public TraversePST {
  public:
    static constexpr std::string_view name = "Interlace";
    struct input {
        int iGridTarget;
        int iGridSource;
    };
    using output = void;
    explicit ServiceInterlace(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceGridBinK : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "GridBinK";
    struct input {
        int nBins;
        int iGrid;
    };
    static constexpr size_t max_k_bins = 2500;
    struct output {
        std::array<double, max_k_bins> fK;
        std::array<double, max_k_bins> fPower;
        std::array<uint64_t, max_k_bins> nPower;
    };
    explicit ServiceGridBinK(PST pst) : TraverseCombinePST(pst, name, sizeof(input), sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte> in, PST) -> int final;
};
