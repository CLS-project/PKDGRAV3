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
#include "cachestats.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceCacheStats::input> || std::is_standard_layout_v<ServiceCacheStats::input>);
static_assert(std::is_void_v<ServiceCacheStats::output> || std::is_standard_layout_v<ServiceCacheStats::output>);

int ServiceCacheStats::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    static_assert(std::is_void_v<input>);
    auto& out = *reinterpret_cast<output*>(vout.data());
    auto mpi = pst->mdl()->mpi;
    std::tie(out.mpi_loop, out.mpi_test, out.gpu_test, out.msg_test) = mpi->get_cache_stats();
    return sizeof(output);
}