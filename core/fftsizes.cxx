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
#include "fftsizes.h"

// Ensure the communication structures are "standard" so that they can be moved around with "memcpy" (required by MDL)
static_assert(std::is_void_v<ServiceFftSizes::input> || std::is_standard_layout_v<ServiceFftSizes::input>);
static_assert(std::is_void_v<ServiceFftSizes::output> || std::is_standard_layout_v<ServiceFftSizes::output>);

int ServiceFftSizes::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == sizeof(output));
    out->nMaxLocal = pst->mdl->FFTlocalCount(in->nx, in->ny, in->nz, &out->nMaxZ, 0, &out->nMaxY, 0);
    return sizeof(output);
}

int ServiceFftSizes::Combine(std::span<std::byte> vout, std::span<const std::byte> vin) {
    auto out = reinterpret_cast<output*>(vout.data());
    auto in = reinterpret_cast<const output*>(vin.data());
    out->nMaxLocal = std::max(out->nMaxLocal, in->nMaxLocal);
    out->nMaxZ = std::max(out->nMaxZ, in->nMaxZ);
    out->nMaxY = std::max(out->nMaxY, in->nMaxY);
    return sizeof(output);
}
