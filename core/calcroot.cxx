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

/*
** Hopefully we can bypass this step once we figure out how to do the
** Multipole Ewald with reduced multipoles.
*/
static void pkdCalcRoot(PKD pkd,uint32_t uRoot,blitz::TinyVector<double,3> com,MOMC &mom) {
    MOMC mc;
    auto kdn = pkd->tree[uRoot];
    momClearMomc(&mom);
    for (auto &p : *kdn) {
        blitz::TinyVector<double,3> r = p.position() - com;
        auto m = p.mass();
        momMakeMomc(&mc,m,r[0],r[1],r[2]);
        momAddMomc(&mom,&mc);
    }
}

// Ensure the communication structures are "standard" so that they can be moved around with "memcpy" (required by MDL)
static_assert(std::is_void<ServiceCalcRoot::input>()  || std::is_standard_layout<ServiceCalcRoot::input>());
static_assert(std::is_void<ServiceCalcRoot::output>() || std::is_standard_layout<ServiceCalcRoot::output>());

int ServiceCalcRoot::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    pkdCalcRoot(pst->plcl->pkd,in->uRoot,in->com,out->momc);
    return sizeof(output);
}

int ServiceCalcRoot::Combine(void *vout,void *vout2) {
    auto out  = static_cast<output *>(vout);
    auto out2 = static_cast<output *>(vout2);
    momAddMomc(&out->momc,&out2->momc);
    return sizeof(output);
}
