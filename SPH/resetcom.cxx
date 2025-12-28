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
#include "resetcom.h"
#include "pkd.h"

auto ServiceResetCOM::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == 0);
    pst->lcl()->pkd->ResetCOM(in->r_com, in->v_com);
    return 0;
}

void pkdContext::ResetCOM(blitz::TinyVector<double, 3> r_com, blitz::TinyVector<double, 3> v_com) {
    for (auto& p : particles) {
        auto r = p.position();
        p.set_position(p.position() - r_com);
        auto& v = p.velocity();
        v -= v_com;
    }
}
