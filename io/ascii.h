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
#include <cstdint>

class ServiceCompressASCII : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "CompressASCII";
    struct input {
        uint64_t nTotal;
        int iFile;
        int iType;
        int iDim;
    };
    using output = uint64_t;
    explicit ServiceCompressASCII(PST pst) : TraverseCombinePST(pst, name, sizeof(input), sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) -> int override;
};

class ServiceWriteASCII : public TraversePST {
  public:
    static constexpr std::string_view name = "WriteASCII";
    struct input {
        uint64_t nFileOffset;
        char achOutFile[PATH_MAX];
    };
    using output = void;
    explicit ServiceWriteASCII(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};
