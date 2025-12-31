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

class ServiceMaxRung : public TraverseCount<ServiceMaxRung> {
  public:
    static constexpr std::string_view name = "MaxRung";
    using input = void;
    struct output {
        uint8_t uMaxRung;
        auto operator+=(const output& rhs) -> output& {
            uMaxRung = std::max(uMaxRung, rhs.uMaxRung);
            return *this;
        }
    };
    explicit ServiceMaxRung(PST pst) : TraverseCount(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceSetGlobalRung : public TraversePST {
  public:
    static constexpr std::string_view name = "SetGlobalRung";
    struct input {
        uint8_t uMaxRung;
    };
    using output = void;
    explicit ServiceSetGlobalRung(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};
