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
#include "pkd_config.h"
#include "timestep.h"
#include <algorithm>
#include <numeric>

auto ServiceMaxRung::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto& out = *reinterpret_cast<output*>(vout.data());
    out.uMaxRung = pkd->MaxRung();
    return sizeof(output);
}

auto pkdContext::MaxRung() -> uint8_t {
    if (particles.empty()) return 0;
    // return std::ranges::max(particles | std::views::transform([](const auto& p) { return p.new_rung(); }));
    return std::accumulate(particles.begin(), particles.end(), uint8_t{0}, [](uint8_t a, auto& p) { return std::max(a, p.new_rung()); });
}

auto ServiceSetGlobalRung::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    pkd->SetGlobalRung(in->uMaxRung);
    return 0;
}

void pkdContext::SetGlobalRung(uint8_t minDt) {
    for (auto& p : particles) p.set_new_rung(minDt);
}
