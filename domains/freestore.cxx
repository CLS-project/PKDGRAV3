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
#include "freestore.h"

static_assert(std::is_void_v<ServiceFreeStore::input> || std::is_trivial_v<ServiceFreeStore::input>);
static_assert(std::is_void_v<ServiceFreeStore::output> || std::is_trivial_v<ServiceFreeStore::output>);
int ServiceFreeStore::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    static_assert(std::is_void_v<input>);
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    assert(vout.size() == sizeof(output));
    *out = pkd->FreeStore();
    return sizeof(output);
}
