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
#include "updaterung.h"
#include "core/particle.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceUpdateRung::input> || std::is_standard_layout_v<ServiceUpdateRung::input>);
static_assert(std::is_void_v<ServiceUpdateRung::output> || std::is_standard_layout_v<ServiceUpdateRung::output>);

int ServiceUpdateRung::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto const& in = *reinterpret_cast<const input*>(vin.data());
    auto& out = *reinterpret_cast<output*>(vout.data());
    auto pkd = pst->pkd();

    assert(pkd->particles.present(PKD_FIELD::oNewRung) || !pkd->bNoParticleOrder);
    out = 0;
    for (auto& p : pkd->particles) {
        if (p.is_active()) {
            if (p.new_rung() > in.uMaxRung) p.set_new_rung(in.uMaxRung);
            if (p.new_rung() >= in.uMinRung) p.set_rung(p.new_rung());
            else if (p.rung() > in.uMinRung) p.set_rung(in.uMinRung);
        }
        out[p.rung()] += 1;
    }
    pkd->nRung = out;
    return sizeof(output);
}
