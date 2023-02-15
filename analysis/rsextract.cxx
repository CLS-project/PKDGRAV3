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
#include "rsextract.h"
#include <algorithm>
#include "blitz/array.h"
using blitz::TinyVector;
#include <fstream>
using std::ofstream;
using std::ios_base;

struct output_particle {
    uint64_t iOrder;
    TinyVector<float,3> r, v;
};

void ServiceRsExtract::Write(PST pst,void *vin,int nIn,int iGroup,const std::string &filename,int iSegment,int nSegment) {
    auto hdr = static_cast<header *>(vin);
    auto pOrds = reinterpret_cast<input *>(hdr+1);
    pst->plcl->pkd->RsExtract(filename, iSegment, pOrds);
}

void pkdContext::RsExtract(const std::string &filename, int iSegment, uint64_t *pOrds)  {
    auto pOrders_first = static_cast<uint64_t *>(pLite);
    auto pOrders_last = pOrders_first + nRsElements;
    if (!std::is_sorted(pOrders_first,pOrders_last)) std::sort(pOrders_first,pOrders_last);

    ofstream file(filename,ios_base::binary | (iSegment?ios_base::app : ios_base::openmode(0)));
    assert(file.good());

    int id = 0;
    for (auto pOrder = pOrders_first; pOrder!=pOrders_last; ++pOrder) {
        auto iOrder = *pOrder;
        while (iOrder >= pOrds[id+1]) ++id;
        assert(id<mdl->Threads());
        int iIndex = iOrder - pOrds[id];
        assert(id==mdl->Self());
        auto p = particles[iIndex];
        assert(iOrder == p.order());
        output_particle out;
        out.iOrder = p.order();
        out.r = p.position();
        out.v = p.velocity();
        file.write(reinterpret_cast<ofstream::char_type *>(&out),sizeof(out));
    }
}
