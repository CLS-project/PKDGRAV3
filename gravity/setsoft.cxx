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
#include "setsoft.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void<ServiceSetSoft::input>()  || std::is_trivial<ServiceSetSoft::input>());
static_assert(std::is_void<ServiceSetSoft::output>() || std::is_trivial<ServiceSetSoft::output>());

int ServiceSetSoft::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    pkdSetSoft(pkd,in->dSoft);
    return 0;
    }
