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
#include "select.h"

static_assert(std::is_void<ServiceCountSelected::input>()  || std::is_standard_layout<ServiceCountSelected::input>());
static_assert(std::is_void<ServiceCountSelected::output>() || std::is_standard_layout<ServiceCountSelected::output>());
int ServiceCountSelected::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    static_assert(std::is_void<input>());
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nOut==sizeof(output));
    *out = pkd->particles.CountSelected();
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelActives::input>()  || std::is_standard_layout<ServiceSelActives::input>());
static_assert(std::is_void<ServiceSelActives::output>() || std::is_standard_layout<ServiceSelActives::output>());
int ServiceSelActives::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelActive(in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelBlackholes::input>()  || std::is_standard_layout<ServiceSelBlackholes::input>());
static_assert(std::is_void<ServiceSelBlackholes::output>() || std::is_standard_layout<ServiceSelBlackholes::output>());
int ServiceSelBlackholes::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelBlackholes(in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelSpecies::input>()  || std::is_standard_layout<ServiceSelSpecies::input>());
static_assert(std::is_void<ServiceSelSpecies::output>() || std::is_standard_layout<ServiceSelSpecies::output>());
int ServiceSelSpecies::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelSpecies(in->mSpecies,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelGroup::input>()  || std::is_standard_layout<ServiceSelGroup::input>());
static_assert(std::is_void<ServiceSelGroup::output>() || std::is_standard_layout<ServiceSelGroup::output>());
int ServiceSelGroup::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelGroup(in->iGroup,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelMass::input>()  || std::is_standard_layout<ServiceSelMass::input>());
static_assert(std::is_void<ServiceSelMass::output>() || std::is_standard_layout<ServiceSelMass::output>());
int ServiceSelMass::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelMass(in->dMinMass,in->dMaxMass,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelPhaseDensity::input>()  || std::is_standard_layout<ServiceSelPhaseDensity::input>());
static_assert(std::is_void<ServiceSelPhaseDensity::output>() || std::is_standard_layout<ServiceSelPhaseDensity::output>());
int ServiceSelPhaseDensity::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelPhaseDensity(in->dMinDensity,in->dMaxDensity,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelById::input>()  || std::is_standard_layout<ServiceSelById::input>());
static_assert(std::is_void<ServiceSelById::output>() || std::is_standard_layout<ServiceSelById::output>());
int ServiceSelById::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelById(in->idStart,in->idEnd,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelBox::input>()  || std::is_standard_layout<ServiceSelBox::input>());
static_assert(std::is_void<ServiceSelBox::output>() || std::is_standard_layout<ServiceSelBox::output>());
int ServiceSelBox::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelBox(in->dCenter,in->dSize,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelSphere::input>()  || std::is_standard_layout<ServiceSelSphere::input>());
static_assert(std::is_void<ServiceSelSphere::output>() || std::is_standard_layout<ServiceSelSphere::output>());
int ServiceSelSphere::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelSphere(in->r,in->dRadius,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void<ServiceSelCylinder::input>()  || std::is_standard_layout<ServiceSelCylinder::input>());
static_assert(std::is_void<ServiceSelCylinder::output>() || std::is_standard_layout<ServiceSelCylinder::output>());
int ServiceSelCylinder::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input *>(vin);
    auto out  = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkd->particles.SelCylinder(in->dP1,in->dP2,in->dRadius,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
}
