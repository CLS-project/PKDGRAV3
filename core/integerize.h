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

#ifndef CORE_INTEGERIZE_H
#define CORE_INTEGERIZE_H
#include <cstdint>
#include "blitz/array.h"

class Integerize {
    static constexpr std::uint32_t factor = 0x80000000u;
public:
    double  convert(std::int32_t i) const {return i * 1.0/factor; }
    template<int n>
    blitz::TinyVector<double,n>  convert(blitz::TinyVector<std::int32_t,n> i) const {return i * 1.0/factor; }
    std::int32_t convert(double d) const {return d * factor; }
    template<int n>
    blitz::TinyVector<std::int32_t,n> convert(blitz::TinyVector<double,n> d) const {return d * factor; }
};

#endif
