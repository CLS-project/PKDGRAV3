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
#include "rshaloloadids.h"
#include <fstream>
using std::ifstream;
using std::ios_base;
namespace rs {
#include "rockstar/io/io_internal.h"
#include "rockstar/halo.h"
}

void ServiceRsHaloLoadIds::Read(PST pst,uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd) {
    pst->plcl->pkd->RsHaloIdRead(iElement,filename,iBeg,iEnd);
}
void ServiceRsHaloLoadIds::start(PST pst,uint64_t nElements,void *vin,int nIn) {
    auto in = static_cast<input *>(vin);
    pst->plcl->pkd->RsHaloIdStart(nElements,in->bAppend);
}
void ServiceRsHaloLoadIds::finish(PST pst,uint64_t nElements,void *vin,int nIn) {
    pst->plcl->pkd->RsHaloIdFinish(nElements);
}
void pkdContext::RsHaloIdStart(uint64_t nElements,bool bAppend)  {
    if (!bAppend) nRsElements = 0;
}
void pkdContext::RsHaloIdFinish(uint64_t nElements) { }
void pkdContext::RsHaloIdRead(uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd) {
    auto ids = static_cast<uint64_t *>(pLite);
    ifstream fhalo(filename,ios_base::in | ios_base::binary);                   assert(fhalo.good());
    ifstream fpart(filename,ios_base::in | ios_base::binary);                   assert(fpart.good());

    // Read the header, then skip to the first halo that we are supposed to read
    rs::binary_output_header hdr;
    fhalo.read(reinterpret_cast<ifstream::char_type *>(&hdr),sizeof(hdr));      assert(fhalo.good());
    assert(iEnd <= hdr.num_halos);
    fhalo.seekg(iBeg * sizeof(rs::halo),ios_base::cur);

    auto base = sizeof(rs::binary_output_header) + hdr.num_halos * sizeof(rs::halo);

    for (auto i=iBeg; i<iEnd; ++i) {
        rs::halo h;
        // Read the halo, then seek to the correct particle and read it
        fhalo.read(reinterpret_cast<ifstream::char_type *>(&h),sizeof(h));      assert(fhalo.good());
        fpart.seekg(base + h.p_start*sizeof(int64_t));
        fpart.read(reinterpret_cast<ifstream::char_type *>(&ids[nRsElements]),sizeof(ids[nRsElements]));
        ++nRsElements;
        assert(fhalo.good());
    }
}
