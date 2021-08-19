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

static_assert(std::is_void<ServiceCountSelected::input>()  || std::is_trivial<ServiceCountSelected::input>());
static_assert(std::is_void<ServiceCountSelected::output>() || std::is_trivial<ServiceCountSelected::output>());
int ServiceCountSelected::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nOut==sizeof(output));
    *out = pkdCountSelected(pkd);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelBlackholes::input>()  || std::is_trivial<ServiceSelBlackholes::input>());
static_assert(std::is_void<ServiceSelBlackholes::output>() || std::is_trivial<ServiceSelBlackholes::output>());
int ServiceSelBlackholes::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelBlackholes(pkd,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelSpecies::input>()  || std::is_trivial<ServiceSelSpecies::input>());
static_assert(std::is_void<ServiceSelSpecies::output>() || std::is_trivial<ServiceSelSpecies::output>());
int ServiceSelSpecies::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelSpecies(pkd,in->mSpecies,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelGroup::input>()  || std::is_trivial<ServiceSelGroup::input>());
static_assert(std::is_void<ServiceSelGroup::output>() || std::is_trivial<ServiceSelGroup::output>());
int ServiceSelGroup::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelGroup(pkd,in->iGroup,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelMass::input>()  || std::is_trivial<ServiceSelMass::input>());
static_assert(std::is_void<ServiceSelMass::output>() || std::is_trivial<ServiceSelMass::output>());
int ServiceSelMass::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelMass(pkd,in->dMinMass,in->dMaxMass,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelPhaseDensity::input>()  || std::is_trivial<ServiceSelPhaseDensity::input>());
static_assert(std::is_void<ServiceSelPhaseDensity::output>() || std::is_trivial<ServiceSelPhaseDensity::output>());
int ServiceSelPhaseDensity::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelPhaseDensity(pkd,in->dMinDensity,in->dMaxDensity,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelById::input>()  || std::is_trivial<ServiceSelById::input>());
static_assert(std::is_void<ServiceSelById::output>() || std::is_trivial<ServiceSelById::output>());
int ServiceSelById::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelById(pkd,in->idStart,in->idEnd,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelBox::input>()  || std::is_trivial<ServiceSelBox::input>());
static_assert(std::is_void<ServiceSelBox::output>() || std::is_trivial<ServiceSelBox::output>());
int ServiceSelBox::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelBox(pkd,in->dCenter,in->dSize,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelSphere::input>()  || std::is_trivial<ServiceSelSphere::input>());
static_assert(std::is_void<ServiceSelSphere::output>() || std::is_trivial<ServiceSelSphere::output>());
int ServiceSelSphere::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelSphere(pkd,in->r,in->dRadius,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }

static_assert(std::is_void<ServiceSelCylinder::input>()  || std::is_trivial<ServiceSelCylinder::input>());
static_assert(std::is_void<ServiceSelCylinder::output>() || std::is_trivial<ServiceSelCylinder::output>());
int ServiceSelCylinder::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn==sizeof(input));
    assert(nOut==sizeof(output));
    *out = pkdSelCylinder(pkd,in->dP1,in->dP2,in->dRadius,in->setIfTrue,in->clearIfFalse);
    return sizeof(output);
    }
