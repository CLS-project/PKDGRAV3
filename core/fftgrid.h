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

class ServiceFftSizes : public TraverseCount<ServiceFftSizes> {
  public:
    static constexpr std::string_view name = "FftSizes";
    struct input {
        int nx, ny, nz;
    };
    struct output {
        uint64_t nMaxLocal;
        int nMaxZ;
        int nMaxY;
        // Combine by taking maximums
        auto operator+=(const output& o) {
            nMaxLocal = std::max(nMaxLocal, o.nMaxLocal);
            nMaxZ = std::max(nMaxZ, o.nMaxZ);
            nMaxY = std::max(nMaxY, o.nMaxY);
            return *this;
        }
    };
    explicit ServiceFftSizes(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) override;
};

class ServiceGridCreateFFT : public TraverseCount<ServiceGridCreateFFT> {
  public:
    static constexpr std::string_view name = "GridCreateFFT";
    struct input {
        int nGrid;
    };
    struct output {
        size_t nMaxLocal;
        int nMaxZ;
        int nMaxY;
        auto operator+=(const output& o) { // Combine by taking maximums
            nMaxLocal = std::max(nMaxLocal, o.nMaxLocal);
            nMaxZ = std::max(nMaxZ, o.nMaxZ);
            nMaxY = std::max(nMaxY, o.nMaxY);
            return *this;
        }
    };
    explicit ServiceGridCreateFFT(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

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
