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
#include "TraversePST.h"

int TraversePST::operator()(std::span<std::byte> In, std::span<std::byte> Out) {
    // This is not really true; the output buffer may be less than the max size
    // assert(Out.size() >= size_t(getMaxBytesOut()) && "Output buffer too small for TraversePST service");
    return Traverse(node_pst, In, Out.subspan(0, getMaxBytesOut()));
}

int TraversePST::Traverse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    if (pst->AmCore()) return Service(pst, vin, vout);
    else if (pst->OffNode()) return OffNode(pst, vin, vout);
    else if (pst->AmNode()) return AtNode(pst, vin, vout);
    else return Recurse(pst, vin, vout);
}

// This is a static version of traverse that we use when we want to traverse a subtree,
// but for a difference service. We look it up in the MDL and call the real traverse.
// This relies on the service being a "TraversePST" and if PST was converted to C++ then
// there would be a more elegant way of doing this.
int TraversePST::Traverse(unsigned sid, PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto service = dynamic_cast<TraversePST*>(mdl->GetService(sid));
    assert(service); // It would be a serious error if this was not a TraversePST service.
    return service->Traverse(pst, vin, vout);
}
int TraversePST::Traverse(std::string_view service_name, PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto service = dynamic_cast<TraversePST*>(mdl->GetService(service_name));
    assert(service); // It would be a serious error if this was not a TraversePST service.
    return service->Traverse(pst, vin, vout);
}

// This is the default function for when we are still traversing the PST
// (OffNode,AtNode,Recurse) and will just naively recurse. Once we reach
// the leaf (AmCore) then Service() is called and must be provided.
int TraversePST::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto rID = ReqService(pst, vin);
    Traverse(pst->lower(), vin, vout);
    return mdl->GetReply(rID, vout);
}

// This class is used for services that pass the same input (or none),
// and need to combine output results of the same size and type.
int TraverseCombinePST::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto rID = ReqService(pst, vin);
    auto nOut1 = Traverse(pst->lower(), vin, vout);
    std::span<std::byte> vout2{reinterpret_cast<std::byte*>(alloca(vout.size())), vout.size()};
    auto nOut2 = mdl->GetReply(rID, vout2);
    return Combine(vout.subspan(0, nOut1), vout2.subspan(0, nOut2), pst);
}