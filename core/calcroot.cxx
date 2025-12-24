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
#include "calcroot.h"

/*
** Hopefully we can bypass this step once we figure out how to do the
** Multipole Ewald with reduced multipoles.
*/
static void pkdCalcRoot(PKD pkd, uint32_t uRoot, blitz::TinyVector<double, 3> com, MOMC& mom) {
    MOMC mc;
    auto kdn = pkd->tree[uRoot];
    momClearMomc(&mom);
    for (auto& p : *kdn) {
        blitz::TinyVector<double, 3> r = p.position() - com;
        auto m = p.mass();
        momMakeMomc(&mc, m, r[0], r[1], r[2]);
        momAddMomc(&mom, &mc);
    }
}

// Ensure the communication structures are "standard" so that they can be moved around with "memcpy" (required by MDL)
static_assert(std::is_void_v<ServiceCalcRoot::input> || std::is_standard_layout_v<ServiceCalcRoot::input>);
static_assert(std::is_void_v<ServiceCalcRoot::output> || std::is_standard_layout_v<ServiceCalcRoot::output>);

int ServiceCalcRoot::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == sizeof(output));
    pkdCalcRoot(pst->plcl->pkd, in->uRoot, in->com, out->momc);
    return sizeof(output);
}

int ServiceCalcRoot::Combine(std::span<std::byte> vout, std::span<std::byte> vout2) {
    auto out = reinterpret_cast<output*>(vout.data());
    auto out2 = reinterpret_cast<output*>(vout2.data());
    momAddMomc(&out->momc, &out2->momc);
    return sizeof(output);
}
