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
#include "combinebound.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void<ServiceCombineBound::input>()  || std::is_trivial<ServiceCombineBound::input>());
static_assert(std::is_void<ServiceCombineBound::output>() || std::is_trivial<ServiceCombineBound::output>());

int ServiceCombineBound::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    static_assert(std::is_void<input>());
    auto out = static_cast<output *>(vout);
    assert(nIn==0);
    assert(nOut==sizeof(output));
    *out = pst->plcl->pkd->bnd;
    return sizeof(output);
}

int ServiceCombineBound::Combine(void *vout,void *vout2) {
    auto out  = static_cast<output *>(vout);
    auto out2 = static_cast<output *>(vout2);
    *out = out->combine(*out2);
    return sizeof(output);
}
