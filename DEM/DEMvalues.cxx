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
#include "DEMvalues.h"
#include "SPH/SPHpredict.h"
#include "pkd.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServicePredictDEMValues::input> || std::is_trivial_v<ServicePredictDEMValues::input>);
static_assert(std::is_void_v<ServicePredictDEMValues::output> || std::is_trivial_v<ServicePredictDEMValues::output>);

int ServicePredictDEMValues::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast</*const*/ input*>(vin.data());
    static_assert(std::is_void_v<output>);
    auto pkd = pst->pkd();
    pkd->PredictDEMValues(&in->kick, &in->SPHoptions);
    return 0;
}

void pkdContext::PredictDEMValues(struct pkdKickParameters* kick, SPHOptions* SPHoptions) {
    for (auto& p : particles) {
        DEMpredict(this, p, kick, SPHoptions->nPredictRung, SPHoptions);
    }
}