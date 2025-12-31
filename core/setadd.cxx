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
#include "setadd.h"
#include "pst.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceSetAdd::input> || std::is_trivial_v<ServiceSetAdd::input>);
static_assert(std::is_void_v<ServiceSetAdd::output> || std::is_trivial_v<ServiceSetAdd::output>);

// This is a weird service. Normally we use TraversePST, but here
// we are actually setting up the PST structure.
int ServiceSetAdd::operator()(std::span<std::byte> in, std::span<std::byte> out) {
    assert(in.size() == sizeof(input));
    SetAdd(node_pst, *reinterpret_cast<input*>(in.data()));
    return 0;
}

void ServiceSetAdd::SetAdd(PST pst, input in) {
    auto mdl = pst->mdl();
    assert(in.idLower == mdl->Self());
    auto n = in.idUpper - in.idLower;
    if (n > 1) {
        auto idMiddle = (in.idUpper + in.idLower) / 2;
        /* Make sure that the pst lands on core zero */
        auto iProcLower = mdl->ThreadToProc(in.idLower);
        auto iProcUpper = mdl->ThreadToProc(in.idUpper - 1);
        if (iProcLower != iProcUpper) {
            idMiddle = mdl->ProcToThread(mdl->ThreadToProc(idMiddle));
        }
        input lopart(in.idLower, idMiddle);
        input hipart(idMiddle, in.idUpper);
        auto low = pst->create_lower(idMiddle, in.idUpper);
        auto rID = mdl->ReqService(pst->upper(), getServiceID(), hipart);
        SetAdd(low, lopart);
        mdl->GetReply(rID);
    }
}