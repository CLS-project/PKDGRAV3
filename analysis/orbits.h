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

#ifndef ORBITS_H
#define ORBITS_H

#include "TraversePST.h"
#include <compare>

class ServiceGetOrbits : public TraversePST {
  public:
    static constexpr std::string_view name = "GetOrbits";
    using input = uint64_t;
    struct output { // Array of these
        output() = default;
        uint64_t id;
        float mass, phi;
        blitz::TinyVector<float, 3> r, v;
    };
    static constexpr size_t max_orbits = 1000;

    explicit ServiceGetOrbits(PST pst) : TraversePST(pst, name, sizeof(input) * max_orbits, sizeof(output) * max_orbits) {}
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif // ORBITS_H