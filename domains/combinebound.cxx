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
#include "combinebound.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceCombineBound::input> || std::is_standard_layout_v<ServiceCombineBound::input>);
static_assert(std::is_void_v<ServiceCombineBound::output> || std::is_standard_layout_v<ServiceCombineBound::output>);

int ServiceCombineBound::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    static_assert(std::is_void_v<input>);
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == 0);
    assert(vout.size() == sizeof(output));
    *out = pst->pkd()->bnd;
    return sizeof(output);
}

int ServiceCombineBound::Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) {
    auto out1 = reinterpret_cast<output*>(vout1.data());
    auto out2 = reinterpret_cast<const output*>(vout2.data());
    *out1 = out1->combine(*out2);
    return sizeof(output);
}
