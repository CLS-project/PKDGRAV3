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
#include "rsreorder.h"
#include <algorithm>
#include <tuple>

int ServiceRsReorderIds::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto pOrds = reinterpret_cast<input *>(vin);
    pst->plcl->pkd->RsReorder(pOrds);
    return 0;
}
using dd_offset_type = mdl::mdlClass::dd_offset_type;
using ordinal_t = uint64_t;
void pkdContext::RsReorder(uint64_t *globalOrds) {
    auto pOrders_first = static_cast<ordinal_t *>(pLite);
    auto pOrders_last = pOrders_first + nRsElements;
    auto N = EphemeralBytes() * FreeStore() / sizeof(ordinal_t);

    // Order particles (this is in domain order)
    if (!std::is_sorted(pOrders_first,pOrders_last)) std::sort(pOrders_first,pOrders_last);
    pOrders_last = std::unique(pOrders_first,pOrders_last);
    nRsElements = pOrders_last - pOrders_first;

    std::vector<dd_offset_type> counts;
    counts.reserve(std::max(mdl->Procs(),mdl->Cores())+1);

    // First do a global partition
    if (mdl->Procs()>1) {
        std::vector<ordinal_t> procOrds;
        procOrds.reserve(mdl->Procs()+1);
        for (auto i=0; i<=mdl->Procs(); ++i) procOrds.push_back(globalOrds[mdl->ProcToThread(i)]);

        int iProc = 0;
        counts.clear();
        dd_offset_type count = 0;
        for (auto pOrder = pOrders_first; pOrder!=pOrders_last; ++pOrder) {
            while (*pOrder >= procOrds[iProc+1]) {
                counts.push_back(count);
                count = 0;
                ++iProc;
                assert(iProc<mdl->Procs());
            }
            ++count;
        }
        counts.push_back(count); // might be zero, usually not
        while (counts.size()<mdl->Procs()) counts.push_back(0);

        dd_offset_type xi = 0;
        for (auto i=0; i<mdl->Proc(); ++i) xi += counts[i];
        auto p1 = pOrders_first + xi;
        auto p2 = p1 + counts[mdl->Proc()];
        for (auto pOrder = p1; pOrder!=p2; ++pOrder) {
            assert(*pOrder >= procOrds[mdl->Proc()] && *pOrder < procOrds[mdl->Proc()+1]);
        }

        auto n = mdl->swapglobal(pOrders_first,N,sizeof(ordinal_t),counts.data());

        pOrders_last = pOrders_first + n;

        for (auto pOrder = pOrders_first; pOrder!=pOrders_last; ++pOrder) {
            assert(*pOrder >= procOrds[mdl->Proc()] && *pOrder < procOrds[mdl->Proc()+1]);
        }

        std::sort(pOrders_first,pOrders_last);
        pOrders_last = std::unique(pOrders_first,pOrders_last);
        nRsElements = pOrders_last - pOrders_first;
    }

    // Should only be elements for this rank, so start:
    int id = mdl->ProcToThread(mdl->Proc());
    if (nRsElements) { assert(*pOrders_first >= globalOrds[id]); }
    counts.clear();
    dd_offset_type count = 0;
    for (auto pOrder = pOrders_first; pOrder!=pOrders_last; ++pOrder) {
        auto iOrder = *pOrder;
        while (iOrder >= globalOrds[id+1]) {
            counts.push_back(count);
            count = 0;
            ++id;
        }
        count++;
    }
    counts.push_back(count); // might be zero, usually not
    while (counts.size()<mdl->Cores()) counts.push_back(0);
    assert(counts.size()==mdl->Cores());

    nRsElements = mdl->swaplocal(pOrders_first,N,sizeof(ordinal_t),counts.data());
    pOrders_last = pOrders_first + nRsElements;
    std::sort(pOrders_first,pOrders_last);
    pOrders_last = std::unique(pOrders_first,pOrders_last);
    nRsElements = pOrders_last - pOrders_first;

}