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
#ifndef CORE_PARTCLASS_H
#define CORE_PARTCLASS_H
#include "io/fio.h"
#include <compare>
#include <type_traits>

#define PKD_MAX_CLASSES 256

struct PARTCLASS {
    float fMass;          /* Particle mass */
    float fSoft;          /* Current softening */
    int iMat;             /* newSPH material id */
    float accFac;         /* Acceleration factor to facilitate fixed particles*/
    FIO_SPECIES eSpecies; /* Species: dark, star, etc. */
    auto operator<=>(PARTCLASS const& rhs) const = default;
};
static_assert(std::is_trivially_copyable_v<PARTCLASS>);

#endif // CORE_PARTCLASS_H