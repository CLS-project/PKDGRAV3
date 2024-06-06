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
#include "fftsizes.h"

// Ensure the communication structures are "standard" so that they can be moved around with "memcpy" (required by MDL)
static_assert(std::is_void<ServiceFftSizes::input>()  || std::is_standard_layout<ServiceFftSizes::input>());
static_assert(std::is_void<ServiceFftSizes::output>() || std::is_standard_layout<ServiceFftSizes::output>());

int ServiceFftSizes::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    out->nMaxLocal = mdlFFTlocalCount(pst->mdl,in->nx,in->ny,in->nz,&out->nMaxZ,0,&out->nMaxY,0);
    return sizeof(output);
}

int ServiceFftSizes::Combine(void *vout,void *vout2) {
    auto out  = static_cast<output *>(vout);
    auto out2 = static_cast<output *>(vout2);
    out->nMaxLocal = std::max(out->nMaxLocal,out2->nMaxLocal);
    out->nMaxZ = std::max(out->nMaxZ,out2->nMaxZ);
    out->nMaxY = std::max(out->nMaxY,out2->nMaxY);
    return sizeof(output);
}
