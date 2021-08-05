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
#include "calcroot.h"

// Make sure that the communication structure is a simple data structure
#include <type_traits>
static_assert(std::is_standard_layout<ServiceCalcRoot::input>(),"not POD");
static_assert(std::is_standard_layout<ServiceCalcRoot::output>(),"not POD");

int ServiceCalcRoot::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    pkdCalcRoot(pst->plcl->pkd,in->uRoot,in->com,&out->momc);
    return sizeof(output);
    }

int ServiceCalcRoot::Combine(void *vout,void *vout2) {
    auto out  = static_cast<output*>(vout);
    auto out2 = static_cast<output*>(vout2);
    momAddMomc(&out->momc,&out2->momc);
    return sizeof(output);
    }
