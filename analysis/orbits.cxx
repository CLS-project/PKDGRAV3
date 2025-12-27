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
#include "orbits.h"
#include <algorithm>
#include <ranges>

int ServiceGetOrbits::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = span_cast<const input>(vin);
    auto out = span_cast<output>(vout, in.size());
    auto pkd = pst->lcl()->pkd;
    assert(!pkd->bNoParticleOrder);     // We need particle IDs
    assert(std::ranges::is_sorted(in)); // input IDs must be sorted
    size_t nOut = 0;
    for (auto& p : pkd->particles) {
        if (std::ranges::binary_search(in, p.order()) == false) continue;
        assert(nOut < out.size());
        out[nOut].id = p.order();
        out[nOut].mass = p.mass();
        out[nOut].phi = p.have_potential() ? p.potential() : 0.0;
        out[nOut].r = p.position();
        out[nOut].v = p.velocity();
        ++nOut;
    }
    return nOut * sizeof(output);
}