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
#include "swapall.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceSwapAll::input> || std::is_trivial_v<ServiceSwapAll::input>);
static_assert(std::is_void_v<ServiceSwapAll::output> || std::is_trivial_v<ServiceSwapAll::output>);

// Routine to swap all particles.  Note that this does not walk the pst
// but simply works with one other processor.
int ServiceSwapAll::operator()(std::span<std::byte> In, std::span<std::byte> Out) {
    auto in = reinterpret_cast<input*>(In.data());
    static_assert(std::is_void_v<output>);
    auto pkd = node_pst->plcl->pkd;
    assert(In.size() == sizeof(input));
    pkdSwapAll(pkd, in->idSwap);
    return 0;
}
