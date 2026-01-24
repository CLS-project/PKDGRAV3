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
#ifndef IO_CLASSES_H
#define IO_CLASSES_H
#include "TraversePST.h"
#include "core/partclass.h"

class ServiceSwapClasses : public TraversePST {
  public:
    static constexpr std::string_view name = "SwapClasses";
    using input = PARTCLASS;  // vector of PARTCLASS
    using output = PARTCLASS; // vector of PARTCLASS
    explicit ServiceSwapClasses(PST pst) : TraversePST(pst, name, PKD_MAX_CLASSES * sizeof(input), PKD_MAX_CLASSES * sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceGetClasses : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "GetClasses";
    using input = void;
    using output = PARTCLASS; // vector of PARTCLASS
    explicit ServiceGetClasses(PST pst) : TraverseCombinePST(pst, name, 0, PKD_MAX_CLASSES * sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout, std::span<const std::byte> vin, PST) -> int override;
};

class ServiceSetClasses : public TraversePST {
  public:
    static constexpr std::string_view name = "SetClasses";
    using input = PARTCLASS; // vector of PARTCLASS
    using output = void;
    explicit ServiceSetClasses(PST pst) : TraversePST(pst, name, PKD_MAX_CLASSES * sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif // IO_CLASSES_H