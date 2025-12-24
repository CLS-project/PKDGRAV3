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
#include "hostname.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceHostname::input> || std::is_trivial_v<ServiceHostname::input>);
static_assert(std::is_void_v<ServiceHostname::output> || std::is_trivial_v<ServiceHostname::output>);

int ServiceHostname::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    static_assert(std::is_void_v<input>);
    auto out = reinterpret_cast<output*>(vout.data());
    auto iUpper = pst->idUpper - pst->idSelf;
    auto outUp = out + iUpper;
    auto rID = pst->mdl->ReqService(pst->idUpper, getServiceID(), vin);
    auto nLower = Traverse(pst->pstLower, vin, vout);
    pst->mdl->GetReply(rID, vout.size() - nLower, outUp);
    return pst->nLeaves * sizeof(output);
}

int ServiceHostname::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    static_assert(std::is_void_v<input>);
    auto out = reinterpret_cast<output*>(vout.data());
    out->iMpiID = pst->mdl->Self();
    strncpy(out->szHostname, mdlName(pst->mdl), sizeof(out->szHostname));
    out->szHostname[sizeof(out->szHostname) - 1] = 0;
    auto p = strchr(out->szHostname, '.');
    if (p) *p = 0;
    return sizeof(output);
}
