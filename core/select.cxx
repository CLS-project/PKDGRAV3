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

static_assert(std::is_void_v<ServiceCountSelected::input> || std::is_standard_layout_v<ServiceCountSelected::input>);
static_assert(std::is_void_v<ServiceCountSelected::output> || std::is_standard_layout_v<ServiceCountSelected::output>);
int ServiceCountSelected::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    static_assert(std::is_void_v<input>);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nOut == sizeof(output));
    *out = pkd->particles.CountSelected();
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelActives::input> || std::is_standard_layout_v<ServiceSelActives::input>);
static_assert(std::is_void_v<ServiceSelActives::output> || std::is_standard_layout_v<ServiceSelActives::output>);
int ServiceSelActives::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelActive(in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelBlackholes::input> || std::is_standard_layout_v<ServiceSelBlackholes::input>);
static_assert(std::is_void_v<ServiceSelBlackholes::output> || std::is_standard_layout_v<ServiceSelBlackholes::output>);
int ServiceSelBlackholes::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelBlackholes(in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelSpecies::input> || std::is_standard_layout_v<ServiceSelSpecies::input>);
static_assert(std::is_void_v<ServiceSelSpecies::output> || std::is_standard_layout_v<ServiceSelSpecies::output>);
int ServiceSelSpecies::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelSpecies(in->mSpecies, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelGroup::input> || std::is_standard_layout_v<ServiceSelGroup::input>);
static_assert(std::is_void_v<ServiceSelGroup::output> || std::is_standard_layout_v<ServiceSelGroup::output>);
int ServiceSelGroup::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelGroup(in->iGroup, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelMass::input> || std::is_standard_layout_v<ServiceSelMass::input>);
static_assert(std::is_void_v<ServiceSelMass::output> || std::is_standard_layout_v<ServiceSelMass::output>);
int ServiceSelMass::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelMass(in->dMinMass, in->dMaxMass, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelPhaseDensity::input> || std::is_standard_layout_v<ServiceSelPhaseDensity::input>);
static_assert(std::is_void_v<ServiceSelPhaseDensity::output> || std::is_standard_layout_v<ServiceSelPhaseDensity::output>);
int ServiceSelPhaseDensity::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelPhaseDensity(in->dMinDensity, in->dMaxDensity, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelById::input> || std::is_standard_layout_v<ServiceSelById::input>);
static_assert(std::is_void_v<ServiceSelById::output> || std::is_standard_layout_v<ServiceSelById::output>);
int ServiceSelById::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelById(in->idStart, in->idEnd, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelBox::input> || std::is_standard_layout_v<ServiceSelBox::input>);
static_assert(std::is_void_v<ServiceSelBox::output> || std::is_standard_layout_v<ServiceSelBox::output>);
int ServiceSelBox::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelBox(in->dCenter, in->dSize, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelSphere::input> || std::is_standard_layout_v<ServiceSelSphere::input>);
static_assert(std::is_void_v<ServiceSelSphere::output> || std::is_standard_layout_v<ServiceSelSphere::output>);
int ServiceSelSphere::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelSphere(in->r, in->dRadius, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}

static_assert(std::is_void_v<ServiceSelCylinder::input> || std::is_standard_layout_v<ServiceSelCylinder::input>);
static_assert(std::is_void_v<ServiceSelCylinder::output> || std::is_standard_layout_v<ServiceSelCylinder::output>);
int ServiceSelCylinder::Service(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));
    assert(nOut == sizeof(output));
    *out = pkd->particles.SelCylinder(in->dP1, in->dP2, in->dRadius, in->setIfTrue, in->clearIfFalse);
    return sizeof(output);
}
